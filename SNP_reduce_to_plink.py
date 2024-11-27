import argparse
from timeit import default_timer as timer
from datetime import timedelta
from pathlib import Path, PurePath
import subprocess
import os
import shutil
from lib.filter import snp_map_filter, dataset_ID_filter, pheno_ID_filter
from lib.select import select_reduce_by_id,  selected_chr_SNP_map, select_top_list
from lib.merge import merge_temp_file, merge_selected_SNP_by_range
from lib.transform import make_map_file, make_ped_file, extract_raw_to_csv

def main():
    parser = argparse.ArgumentParser(description="Convert reduced file to plink")
    parser.add_argument("--SNP_folder", required=True, type=Path, help="Folder that contains reduced SNP files.")
    parser.add_argument("--SNP_map", required=True, type=Path, help="CSV file that contains list of SNP, with at least chr, SNP_id, BP.")
    parser.add_argument("--patient_file", required=True, type=Path, help="File that contains patients information (including phenotype).")
    parser.add_argument("--id_name", required=True, type=str, help="Column name of patients' identity in both patient file and SNP file.")
    parser.add_argument("--phenotype", required=True, type=str, help="Phenotype's name in the file")
    parser.add_argument("--type", required=True, type=str, help="Phenotype's type in the file, eg category or continuous")
    parser.add_argument("--chr", type=str, help="Select specific chromosomes to analyze (without filter). Could use single or comma-delimited string (eg 1,7,11)")
    parser.add_argument("--merge_size", type=int, default=800, help="Merge temp file by how many cases (default: 800)")
    parser.add_argument("--chunk_size", type=int, default=50000, help="Generate ped file with how many snp per chunk (default: 50000)")
    parser.add_argument("--keep_temp", type=bool, default=False, help="Keep all temp files (need more disk space).")
    parser.add_argument("--make_bed", type=bool, default=True, help="Use plink to make bed file.")
    parser.add_argument("--prune", type=bool, default=True, help="Use plink to QC and prune bed file.")
    parser.add_argument("--make_plot", type=bool, default=True, help="Use R to make manhattan plot and qq plot.")
    parser.add_argument("--rscript", type=str, default="C:\\Program Files\\R\\R-4.3.1\\bin\\Rscript", help="Path to R script file.")
    parser.add_argument("--rfile", type=str, default="SNP_manhattan_top.R", help="Path to R plot script file.")
    parser.add_argument("--snp", type=str, help="Select specific snps to extract. Could use single or comma-delimited string (eg rs1883832,rs11569323)")
    
    args = parser.parse_args()

    # 1-22, X, Y, XY (pseudo-autosomal region of X), and MT
    SNP_map_file = args.SNP_map # chr,SNP_id,BP
    reduce_folder = PurePath(args.SNP_folder)
    patient_data_file = args.patient_file
    patient_data_file_prefix = Path(patient_data_file).stem
    patient_data_id = args.id_name
    phenotype = [args.phenotype, args.type] # category, continuous
    
    if args.type == "category":
        assoc_analysis = "logistic"
    else:
        assoc_analysis = "linear"
    
    if args.chr:
        selected_chr = args.chr.split(",")
    else:
        selected_chr = ["1","2","3","4","5","6","7","8","9","10",
                        "11","12","13","14","15","16","17","18","19","20",
                        "21","22","MT","X","XY","Y"]
    
    # Part 0
    # Preparing folder and sorted SNP files
    dataset_folder = PurePath(patient_data_file_prefix)
    snp_folder = PurePath(dataset_folder, "SNP_files")
    temp_folder = PurePath(snp_folder, "temp")
    merged_folder = PurePath(snp_folder, "merged")
    final_folder = PurePath(dataset_folder, phenotype[0])
    
    # Part 1
    def make_dataset_folder():
        # Generate working folder
        if not os.path.exists(dataset_folder):
            print("Part 1: Generate working folders")
            Path(dataset_folder).mkdir(parents=True, exist_ok=True)    
    
    # Part 2
    dataseta_ID_file = PurePath(dataset_folder, patient_data_file_prefix + r"_ID_list.csv")
    ID_phenotype_file = PurePath(final_folder, patient_data_file_prefix + r"_" + phenotype[0] + r"_list.csv")
    def make_id_files():
        # Use all ID in dataset to generate dataset-level ID filter 
        # Labels: ID S
        if not os.path.exists(dataseta_ID_file):
            print("Part 2: Generate dataset ID filter.")
            Path(Path(dataseta_ID_file).parent.absolute()).mkdir(parents=True, exist_ok=True)
            dataset_ID_filter(patient_data_file, patient_data_id, dataseta_ID_file)
        
        # Generate phenotype-based ID filter (phenotype=-9 if missing)
        # Labels: FID ID P
        if not os.path.exists(ID_phenotype_file):
            print("Part 2: Generate phenotype-based patient ID file.")
            Path(Path(ID_phenotype_file).parent.absolute()).mkdir(parents=True, exist_ok=True)
            pheno_ID_filter(patient_data_file, patient_data_id, phenotype, ID_phenotype_file)
    
    # Part 3
    # Merge dataset-level SNP files
    new_SNP_map_file = PurePath(snp_folder, r"SNP_map_sorted.csv")
    new_SNP_map_range_file = PurePath(snp_folder, r"SNP_map_sorted_range.csv")
    def merge_reduce_snp_files():
        # Generate sorted SNP file
        if not os.path.exists(new_SNP_map_range_file):
            print("Part 3: Generate sorted SNP file and file for split merged file.")
            Path(Path(new_SNP_map_range_file).parent.absolute()).mkdir(parents=True, exist_ok=True)
            snp_map_filter(SNP_map_file, new_SNP_map_file, new_SNP_map_range_file)

        # Generate merged snp files
        if not os.path.exists(merged_folder): 
            if not os.path.exists(temp_folder):
                print("Part 3: Select SNP files by patient ID filter.")
                Path(temp_folder).mkdir(parents=True, exist_ok=True)
                select_reduce_by_id(reduce_folder, dataseta_ID_file, temp_folder)
            
            print("Part 3: Merge selected SNP files.")
            Path(merged_folder).mkdir(parents=True, exist_ok=True)
            merge_temp_file(temp_folder, merged_folder, new_SNP_map_file, args.merge_size)

        # Remove temp folder
        if not args.keep_temp and os.path.exists(temp_folder):
            print("Part 3: Remove temp folder...")
            shutil.rmtree(temp_folder)

    # Part 4
    # Split dataset-level SNP files to chr-level and merge them
    def chr_split_merge():
        for chr in selected_chr:
            # Make selected chr filter 
            chr_folder = PurePath(snp_folder, "chr", str(chr))
            chr_map_file = PurePath(chr_folder, patient_data_file_prefix + r"_chr" + str(chr) + r"_SNP_map.csv")
            if not os.path.exists(chr_map_file):
                print("Part 4: Select chr: " + str(chr))
                Path(Path(chr_map_file).parent.absolute()).mkdir(parents=True, exist_ok=True)
                selected_chr_SNP_map(chr, new_SNP_map_file, chr_map_file)
            
            # Split dataset-level SNP files and merge to chr-level file
            # and merge chr map and patient genotype together
            chr_map_patient_file = PurePath(chr_folder, patient_data_file_prefix + r"_chr" + str(chr) + r"_SNP_map_patient.csv")
            if not os.path.exists(chr_map_patient_file):
                print("Part 4: Generate chr map patient file of chr: "+ str(chr))
                Path(Path(chr_map_patient_file).parent.absolute()).mkdir(parents=True, exist_ok=True)
                merge_selected_SNP_by_range(merged_folder, chr, new_SNP_map_range_file, chr_map_file, chr_map_patient_file)

    # Part 5
    # Make dataset-level Map and Ped file
    def chr_map_ped():
        for chr in selected_chr:
            chr_folder = PurePath(snp_folder, "chr", str(chr))
            chr_map_patient_file = PurePath(chr_folder, patient_data_file_prefix + r"_chr" + str(chr) + r"_SNP_map_patient.csv")
            
            # Generate map file
            chr_map_file = PurePath(snp_folder, "map", patient_data_file_prefix + r"_chr" + str(chr) + r".map")
            if not os.path.exists(chr_map_file):
                print("Part 5: Generate map file of chr: "+ str(chr))
                Path(Path(chr_map_file).parent.absolute()).mkdir(parents=True, exist_ok=True)
                make_map_file(chr_map_patient_file, chr_map_file)

            # Generate dataset-based ped file (without phenotype)
            chr_ped_file = PurePath(snp_folder, "ped", patient_data_file_prefix + r"_chr" + str(chr) + r".ped")
            if not os.path.exists(chr_ped_file):
                print("Part 5: Generate ped file of chr: "+ str(chr))
                Path(Path(chr_ped_file).parent.absolute()).mkdir(parents=True, exist_ok=True)
                make_ped_file(chr_map_patient_file, dataseta_ID_file, args.chunk_size, chr_ped_file)

    # Part 6
    # Generate binary files with plink and write to mergefiles.txt
    merge_files_txt = PurePath(final_folder, patient_data_file_prefix + r"_" + phenotype[0] + r"_mergefiles.txt")
    def make_bed():
        for chr in selected_chr:
            chr_map_file = PurePath(snp_folder, "map", patient_data_file_prefix + r"_chr" + str(chr) + r".map")
            chr_ped_file = PurePath(snp_folder, "ped", patient_data_file_prefix + r"_chr" + str(chr) + r".ped")
            chr_bed = PurePath(snp_folder, "bed", phenotype[0], patient_data_file_prefix + r"_" + phenotype[0] + r"_chr" + str(chr))
            
            if os.path.exists(chr_ped_file):
                if not os.path.exists(PurePath(str(chr_bed) + r".bed")):
                    print("Part 6: Generate bed files, combine phenotyp file, now on chr: "+ str(chr))
                    Path(Path(chr_bed).parent.absolute()).mkdir(parents=True, exist_ok=True)
                    plink_makebed = r".\\plink --silent --map .\\" + str(chr_map_file) + r" --ped .\\" + str(chr_ped_file) + r" --pheno .\\" + str(ID_phenotype_file) + r" --make-bed --out .\\" + str(chr_bed)
                    os.system(plink_makebed)

                    print("Part 6: Write bed path to mergefiles.txt")
                    with open(merge_files_txt, 'a+') as f:
                        f.write(f'{chr_bed}\n')
            else:
                raise FileNotFoundError(chr_ped_file)

    # Part 7
    # Merge all plink binary files in mergefiles.txt and generate merged plink file
    merged_final_name = PurePath(final_folder, patient_data_file_prefix + r"_" + phenotype[0] + r"_merged")
    merged_final_file = PurePath(str(merged_final_name) + r".bed")
    def merge_binary_files():
        if os.path.exists(merge_files_txt):
            # merge binary files + freqx + missing
            if not os.path.exists(merged_final_file):
                print("Part 7: Merge binary files, calculate frequency and missing of phenotype and SNPs.")
                plink_merge = r".\\plink --silent --merge-list .\\" + str(merge_files_txt) + r" --make-bed --freqx --missing --out .\\" + str(merged_final_name)
                os.system(plink_merge)
        else:
            raise FileNotFoundError(merge_files_txt)
        
    # Part 8
    # Prune merged binary file
    merged_final_QC_name = PurePath(final_folder, patient_data_file_prefix + r"_" + phenotype[0] + r"_merged_QC_prune")
    merged_final_QC_file = PurePath(str(merged_final_QC_name) + r".assoc." + assoc_analysis)
    def prune_binary_files():
        if os.path.exists(merged_final_file):
           # QC + prune
            if not os.path.exists(merged_final_QC_file):
                print("Part 8: Prune merged binary file.")
                plink_prune = r".\\plink --silent --bfile .\\" + str(merged_final_name) + r" --not-chr 0 x y xy --maf 0.05 --hwe 0.000001 --geno 0.05 --mind 0.05 --indep-pairwise 50 5 0.2 --" + assoc_analysis + r" --ci 0.95 --out .\\" + str(merged_final_QC_name)
                os.system(plink_prune)
        else:
            raise FileNotFoundError(merged_final_file)

    # Part 9
    # Use R to make plots
    # change rscript path if different version or path to folder
    top_file_name = PurePath(final_folder, patient_data_file_prefix + r"_" + phenotype[0] + r"_Top_SNPs.csv")
    def R_plot():
        if os.path.exists(merged_final_QC_file):
            if not os.path.exists(top_file_name):
                print("Part 9: Use R to make plots.")
                r_bat = [args.rscript, r"--vanilla", args.rfile, r"--dataset=" + patient_data_file_prefix, r"--phenotype=" + phenotype[0], r"--type=" + phenotype[1]]
                subprocess.call(r_bat, shell=True)
        else:
            raise FileNotFoundError(merged_final_QC_file)
        
    # Part 10
    # Additional functions
    # Select specific snps and extract genotype
    extract_snp_name = PurePath(final_folder, patient_data_file_prefix + r"_" + phenotype[0] + r"_extract")
    def extract_snp():
        if os.path.exists(merged_final_file):
            if args.snp:
                extract_snp_list = args.snp
            else:
                extract_snp_list = select_top_list(top_file_name)
            
            if len(extract_snp_list) != 0:
                print("Part 10: Extract genotype of specific snps.")
                plink_extract = r".\\plink --silent --bfile .\\" + str(merged_final_name) + r" --recodeA include-alt --d ? --snps " + extract_snp_list + r" --out .\\" + str(extract_snp_name)
                os.system(plink_extract)

                print("Part 10: Convert raw file to csv file.")
                extract_raw_to_csv(PurePath(str(extract_snp_name) + r".raw"), patient_data_id, phenotype[0], PurePath(str(extract_snp_name) + r".csv"))
            else:
                print('Part 10: No designated snp or significant snp found in the top snps file!')
        else:
            raise FileNotFoundError(merged_final_file)
    
    # Main workflow
    if args.chr:
        print("Progress selected chr: " + str(args.chr))
    else:
        print("Loop over all chr.")
    
    make_dataset_folder()
    make_id_files()
    merge_reduce_snp_files()
    chr_split_merge()
    chr_map_ped()

    if args.make_bed:
        make_bed()
        merge_binary_files()
    if args.prune:
        prune_binary_files()
    if args.make_plot:
        R_plot()
        extract_snp()

if __name__ == "__main__":
    try:
        start = timer()
        main()
        end = timer()
        print("All tasks were done! This round spent " + str(timedelta(seconds=end-start)).split(".")[0] + ".")
    except Exception as e:
        print(e)