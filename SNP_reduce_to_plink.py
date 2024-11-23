import argparse
from timeit import default_timer as timer
from datetime import timedelta
from pathlib import Path, PurePath
import os
import shutil
from lib.filter import generate_ID_filter, remove_out_of_range_SNP, generate_divers_SNP_filter
from lib.select import select_reduce_by_id,  selected_chr_SNP, intersection_of_SNP, out_range_out_intersection_SNP
from lib.merge import merge_temp_file, filted_SNP_merge, divers_SNP_merge
from lib.transform import make_map_file, make_ped_file

def main():
    parser = argparse.ArgumentParser(description="Convert reduced file to plink")
    parser.add_argument("--folder", required=True, type=str, help="Folder that contains reduced SNP files.")
    parser.add_argument("--list", required=True, type=str, help="CSV file that contains list of SNP, with at least chr, SNP_id, BP.")
    parser.add_argument("--file", required=True, type=str, help="File that contains patients information (including phenotype).")
    parser.add_argument("--id", required=True, type=str, help="Column name of patients' identity in both patient file and SNP file.")
    parser.add_argument("--pheno", required=True, type=str, help="Phenotype's name in the file")
    parser.add_argument("--type", required=True, type=str, help="Phenotype's type in the file, eg category or continuous")
    parser.add_argument("--chr", type=str, help="Select specific chromosomes to analyze (without filter). Could use single or comma-delimited string (eg 1,7,11)")
    parser.add_argument("--use_filter", type=bool, default=True, help="Use filter 1~16 to remove indels")
    parser.add_argument("--check_indel", type=bool, help="Also check indel SNPs (it might take 2x~3x times)")
    parser.add_argument("--only_check_indel", type=bool, default=False, help="ONLY check indel SNPs, need merged files (it might take 2x~3x times)")
    parser.add_argument("--keep_temp", type=bool, default=False, help="Keep all temp files (need more disk space).")
    parser.add_argument("--make_bed", type=bool, default=True, help="Use plink to make bed file.")
    parser.add_argument("--QC", type=bool, default=True, help="Use plink to QC bed file.")
    parser.add_argument("--prune", type=bool, default=True, help="Use plink to QC and prune bed file.")
    args = parser.parse_args()

    start = timer()

    # 1-22, X, Y, XY (pseudo-autosomal region of X), and MT
    SNP_list = args.list # chr,SNP_id,BP
    reduce_folder = PurePath(args.folder).as_posix()
    patient_data_file = args.file
    patient_data_file_prefix = Path(patient_data_file).stem
    patient_data_id = args.id
    phenotype = [args.pheno, args.type] # category, continuous

    if args.chr:
        selected_chr = args.chr.split(",")
    else:
        selected_chr = ["1","2","3","4","5","6","7","8","9","10",
                        "11","12","13","14","15","16","17","18","19","20",
                        "21","22","MT","X","XY","Y"]
    if args.only_check_indel:
        args.check_indel = True

    temp_folder = PurePath("Selected", phenotype[0], "temp").as_posix()
    merged_folder = PurePath("Selected", phenotype[0], "merged").as_posix()
    filted_folder = PurePath("Selected", phenotype[0], "filted").as_posix()
    chr_folder = PurePath("Selected", phenotype[0], "chr").as_posix()
    
    ID_filter_file_raw = PurePath("Selected", phenotype[0], patient_data_file_prefix + r"_list_" + phenotype[0] + r"_raw.csv").as_posix()
    
    common_SNP_file = PurePath("Selected", phenotype[0], patient_data_file_prefix + r"_common_SNP_" + phenotype[0] + r".csv").as_posix()
    all_SNP_file = PurePath("Selected", phenotype[0], patient_data_file_prefix + r"_all_SNP_" + phenotype[0] + r".csv").as_posix()
    out_range_SNP_file = PurePath("Selected", phenotype[0], patient_data_file_prefix + r"_out_range_uncommon_SNP_" + phenotype[0] + r".csv").as_posix()
    patient_common_SNP_file = PurePath("Selected", phenotype[0], patient_data_file_prefix + r"_patient_common_SNP_" + phenotype[0] + r".csv").as_posix()
    patient_outrange_SNP_file = PurePath("Selected", phenotype[0], patient_data_file_prefix + r"_patient_outrange_SNP_" + phenotype[0] + r".csv").as_posix()
    
    if not os.path.exists(ID_filter_file_raw):
        print("Step: generate patient ID filter.")
        Path(Path(ID_filter_file_raw).parent.absolute()).mkdir(parents=True, exist_ok=True)
        generate_ID_filter(patient_data_file, patient_data_id, phenotype, ID_filter_file_raw)

    if not os.path.exists(merged_folder):
        print("Step: select SNP files by patient ID filter.")
        Path(temp_folder).mkdir(parents=True, exist_ok=True)
        select_reduce_by_id(reduce_folder, ID_filter_file_raw, temp_folder)

        print("Step: merge selected SNP files.")
        Path(merged_folder).mkdir(parents=True, exist_ok=True)
        merge_temp_file(temp_folder, merged_folder, SNP_list, 800)

        if not args.keep_temp and os.path.exists(temp_folder):
            print("Removing temp folder...")
            shutil.rmtree(temp_folder)
    
    def chr_of_interest(selected_chr):
        for chr in selected_chr:
            selected_chr_file = PurePath("Selected", phenotype[0], "chr", patient_data_file_prefix + r"_chr" + str(chr) + r"_SNP_" + phenotype[0] + r".csv").as_posix()
            patient_select_chr_file = PurePath("Selected", phenotype[0], "chr", patient_data_file_prefix + r"_patient_chr" + str(chr) + r"_SNP_" + phenotype[0] + r".csv").as_posix()
            patient_chr_snp_map = PurePath("Selected", phenotype[0], "chr", patient_data_file_prefix + r"_chr" + str(chr) + r"_SNP_map_" + phenotype[0] + r".csv").as_posix()
            patient_chr_snp_file = PurePath("Selected", phenotype[0], "chr", patient_data_file_prefix + r"_patient_chr" + str(chr) + r"_SNP_map_" + phenotype[0] + r".csv").as_posix()
            Path(Path(selected_chr_file).parent.absolute()).mkdir(parents=True, exist_ok=True)
            
            chr_map_file = PurePath(phenotype[0], "chr", patient_data_file_prefix + r"_chr" + str(chr) + r"_" + phenotype[0] + r".map").as_posix()
            chr_ped_file = PurePath(phenotype[0], "chr", patient_data_file_prefix + r"_chr" + str(chr) + r"_" + phenotype[0] + r".ped").as_posix()
            chr_bed_file = PurePath(phenotype[0], "chr", patient_data_file_prefix + r"_chr" + str(chr) + r"_" + phenotype[0] + r".bed").as_posix()
            Path(Path(chr_ped_file).parent.absolute()).mkdir(parents=True, exist_ok=True)
            
            if not os.path.exists(selected_chr_file):
                print("select chr: " + str(chr))
                selected_chr_SNP(SNP_list, chr, selected_chr_file)

            if not os.path.exists(patient_select_chr_file):
                print("generate patient chr file...")
                filted_SNP_merge(merged_folder, selected_chr_file, patient_select_chr_file)

            if not os.path.exists(patient_chr_snp_map):
                print("generate patient SNP map...")
                generate_divers_SNP_filter(patient_select_chr_file, SNP_list,  patient_chr_snp_map)

            if not os.path.exists(patient_chr_snp_file):
                print("generate patient SNP file...")
                divers_SNP_merge(patient_select_chr_file, patient_chr_snp_map, patient_chr_snp_file)

            if not os.path.exists(chr_map_file):
                print("generate map file...")
                make_map_file(patient_chr_snp_file, chr_map_file)

            if not os.path.exists(chr_ped_file):
                print("generate ped file...")
                make_ped_file(patient_chr_snp_file, ID_filter_file_raw, chr_ped_file)

            if args.make_bed:
                if not os.path.exists(chr_bed_file):
                    print("generate bed file")
                    map_ped_file_prefix = Path(chr_ped_file).stem
                    map_ped_file_relative = PurePath(phenotype[0], "chr", Path(chr_ped_file).stem).as_posix()

                    plink_makebed = ".\plink --file .\\" + phenotype[0] + "\\chr\\" + map_ped_file_prefix + " --make-bed --out .\\" + phenotype[0] + "\\chr\\" + map_ped_file_prefix
                    os.system(plink_makebed)

                    print("write to mergefiles.txt")
                    merge_files_txt = PurePath(phenotype[0], "chr", "mergefiles.txt")
                    with open(merge_files_txt, 'a+') as f:
                        f.write(f'{map_ped_file_relative}\n')

        # merge binary files
        print("merging binary files")
        plink_merge = ".\plink --merge-list .\\" + phenotype[0] + "\\chr\\mergefiles.txt --make-bed --out .\\" + phenotype[0] + "\\merged"
        os.system(plink_merge)

        # freqx
        print("calculate frequency of each SNP")
        plink_freq = ".\plink --bfile .\\" + phenotype[0] + "\\merged --freqx --out .\\" + phenotype[0] + "\\merged"
        os.system(plink_freq)

        if args.type == "category":
            assoc_analysis = "--logistic"
        else:
            assoc_analysis = "--linear"
        
        # no QC
        plink_whole = ".\plink --bfile .\\" + phenotype[0] + "\\merged --not-chr 0 x y xy " + assoc_analysis + " --ci 0.95 --out .\\" + phenotype[0] + "\\merged_wo_QC"
        os.system(plink_whole)

        # QC
        if args.QC:
            plink_QC = ".\plink --bfile .\\" + phenotype[0] + "\\merged --not-chr 0 x y xy --maf 0.05 --hwe 0.000001 --geno 0.05 --mind 0.05 " + assoc_analysis + " --ci 0.95 --out .\\" + phenotype[0] + "\\merged_QC"
            os.system(plink_QC)

        # QC + prune
        if args.prune:
            plink_prune = ".\plink --bfile .\\" + phenotype[0] + "\\merged --not-chr 0 x y xy --maf 0.05 --hwe 0.000001 --geno 0.05 --mind 0.05 --indep-pairwise 50 5 0.2 " + assoc_analysis + " --ci 0.95 --out .\\" + phenotype[0] + "\\merged_QC_prune"
            os.system(plink_prune)

    if args.chr:
        print("Step: progressing chr of interest.")
        chr_of_interest(selected_chr)
    else:
        print("Step: loop over all chr.")
        chr_of_interest(selected_chr)

        # if args.use_filter and not os.path.exists(common_SNP_file):
        #     print("Step: remove SNPs not in range.")
        #     Path(filted_folder).mkdir(parents=True, exist_ok=True)
        #     remove_out_of_range_SNP(merged_folder, filted_folder)

        #     print("Step: make common SNPs list.")
        #     intersection_of_SNP(filted_folder, common_SNP_file)
        # elif not args.use_filter and not os.path.exists(all_SNP_file):
        #     print("Step: make common SNPs without filter.")
        #     intersection_of_SNP(merged_folder, all_SNP_file)

        # if args.check_indel:
        #     print("Step: make outrange SNPs list.")
        #     out_range_out_intersection_SNP(SNP_list, common_SNP_file, out_range_SNP_file)

        # if not args.only_check_indel and not os.path.exists(patient_common_SNP_file):
        #     print("Step 6: merge filted SNP files with common SNPs list.")
        #     filted_SNP_merge(filted_folder, common_SNP_file, patient_common_SNP_file)
        
        # if args.check_indel and not os.path.exists(patient_outrange_SNP_file):
        #     print("Step 6: merge filted SNP files with outrange SNPs list.")
        #     filted_SNP_merge(merged_folder, out_range_SNP_file, patient_outrange_SNP_file)

        # if not args.keep_temp and os.path.exists(filted_folder):
        #     print("Removing filted folder...")
        #     shutil.rmtree(filted_folder)

        # divers_snp_map = PurePath("Selected", phenotype[0], patient_data_file_prefix + r"_divers_SNP_map_" + phenotype[0] + r".csv").as_posix()
        # outrange_divers_snp_map = PurePath("Selected", phenotype[0], patient_data_file_prefix + r"_divers_outrange_SNP_map_" + phenotype[0] + r".csv").as_posix()
        # if not args.only_check_indel and not os.path.exists(divers_snp_map):
        #     generate_divers_SNP_filter(patient_common_SNP_file, SNP_list, divers_snp_map)
        # if args.check_indel and not os.path.exists(outrange_divers_snp_map):
        #     generate_divers_SNP_filter(patient_outrange_SNP_file, SNP_list, outrange_divers_snp_map)

        # patient_divers_SNP_file = PurePath("Selected", phenotype[0], patient_data_file_prefix + r"_divers_SNP_" + phenotype[0] + r".csv").as_posix()
        # patient_divers_outrange_SNP_file = PurePath("Selected", phenotype[0], patient_data_file_prefix + r"_divers_outrange_SNP_" + phenotype[0] + r".csv").as_posix()
        # if not args.only_check_indel and not os.path.exists(patient_divers_SNP_file):
        #     divers_SNP_merge(patient_common_SNP_file, divers_snp_map, patient_divers_SNP_file)
        # if args.check_indel and not os.path.exists(patient_divers_outrange_SNP_file):
        #     divers_SNP_merge(patient_outrange_SNP_file, outrange_divers_snp_map, patient_divers_outrange_SNP_file)

        # map_file = PurePath(phenotype[0], patient_data_file_prefix + r"_" + phenotype[0] + r".map").as_posix()
        # outrange_map_file = PurePath(phenotype[0], patient_data_file_prefix + r"_" + phenotype[0] + r"_outrange.map").as_posix()
        # Path(Path(map_file).parent.absolute()).mkdir(parents=True, exist_ok=True)
        # if not args.only_check_indel and not os.path.exists(map_file):
        #     make_map_file(patient_divers_SNP_file, map_file)
        # if args.check_indel and not os.path.exists(outrange_map_file):
        #     make_map_file(patient_divers_outrange_SNP_file, outrange_map_file)

        # ped_file= PurePath(phenotype[0], patient_data_file_prefix + r"_" + phenotype[0] + r".ped").as_posix()
        # outrange_ped_file= PurePath(phenotype[0], patient_data_file_prefix + r"_" + phenotype[0] + r"_outrange.ped").as_posix()
        # Path(Path(ped_file).parent.absolute()).mkdir(parents=True, exist_ok=True)
        # if not args.only_check_indel and not os.path.exists(ped_file):
        #     make_ped_file(patient_divers_SNP_file, ID_filter_file_raw, ped_file)
        # if args.check_indel and not os.path.exists(outrange_ped_file):
        #     make_ped_file(patient_divers_outrange_SNP_file, ID_filter_file_raw, outrange_ped_file)

        # map_ped_file_prefix = Path(ped_file).stem
        # # make bed. bim. fam
        # if args.make_bed:
        #     if not args.only_check_indel:
        #         plink_string1 = ".\plink --file .\\" + phenotype[0] + "\\" + map_ped_file_prefix + " --make-bed --out .\\" + phenotype[0] + "\\" + map_ped_file_prefix
        #         os.system(plink_string1)
        #     if args.check_indel:
        #         plink_string1 = ".\plink --file .\\" + phenotype[0] + "\\" + map_ped_file_prefix + "_outrange --make-bed --out .\\" + phenotype[0] + "\\" + map_ped_file_prefix + "_outrange"
        #         os.system(plink_string1)
        #     # freqx
        #     if not args.only_check_indel:
        #         plink_freq = ".\plink --bfile .\\" + phenotype[0] + "\\" + map_ped_file_prefix + " --freqx --out .\\" + phenotype[0] + "\\" + map_ped_file_prefix
        #         os.system(plink_freq)
        #     if args.check_indel:
        #         plink_freq = ".\plink --file .\\" + phenotype[0] + "\\" + map_ped_file_prefix + "_outrange --freqx --out .\\" + phenotype[0] + "\\" + map_ped_file_prefix + "_outrange"
        #         os.system(plink_freq)  
                
        # # QC
        # if args.type == "category":
        #     assoc_analysis = "--logistic"
        # else:
        #     assoc_analysis = "--linear"
        
        # if args.QC:
        #     if not args.only_check_indel:
        #         plink_string2 = ".\plink --bfile .\\" + phenotype[0] + "\\" + map_ped_file_prefix + " --not-chr 0 x y xy --maf 0.05 --hwe 0.000001 --geno 0.05 --mind 0.05 " + assoc_analysis + " --ci 0.95 --out .\\" + phenotype[0] + "\\" + map_ped_file_prefix + "_QC"
        #     if args.check_indel:
        #         plink_string2 = ".\plink --bfile .\\" + phenotype[0] + "\\" + map_ped_file_prefix + "_outrange --not-chr 0 x y xy --maf 0.05 --hwe 0.000001 --geno 0.05 --mind 0.05 " + assoc_analysis + " --ci 0.95 --out .\\" + phenotype[0] + "\\" + map_ped_file_prefix + "_outrange_QC"
        #     os.system(plink_string2)

        # # QC + prune
        # if args.prune:
        #     if not args.only_check_indel:
        #         plink_string3 = ".\plink --bfile .\\" + phenotype[0] + "\\" + map_ped_file_prefix + " --not-chr 0 x y xy --maf 0.05 --hwe 0.000001 --geno 0.05 --mind 0.05 --indep-pairwise 50 5 0.2 " + assoc_analysis + " --ci 0.95 --out .\\" + phenotype[0] + "\\" + map_ped_file_prefix + "_QC_prune"
        #         os.system(plink_string3)
        #     if args.check_indel:
        #         plink_string3 = ".\plink --bfile .\\" + phenotype[0] + "\\" + map_ped_file_prefix + "_outrange --not-chr 0 x y xy --maf 0.05 --hwe 0.000001 --geno 0.05 --mind 0.05 --indep-pairwise 50 5 0.2 " + assoc_analysis + " --ci 0.95 --out .\\" + phenotype[0] + "\\" + map_ped_file_prefix + "_outrange_QC_prune"
        #         os.system(plink_string3)

    end = timer()
    print("All tasks were done! This round spent " + str(timedelta(seconds=end-start)).split(".")[0] + ".")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(e)