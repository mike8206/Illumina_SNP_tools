import argparse
from timeit import default_timer as timer
from datetime import timedelta
from pathlib import Path, PurePath
import os
import shutil
from lib.filter import generate_ID_filter, remove_out_of_range_SNP, generate_divers_SNP_filter
from lib.select import select_reduce_by_id, intersection_of_SNP, out_range_out_intersection_SNP
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
    parser.add_argument("--check_indel", type=bool, help="Also check indel SNPs (it might take 2x~3x times)")
    parser.add_argument("--only_check_indel", type=bool, help="Also check indel SNPs (it might take 2x~3x times)")
    parser.add_argument("--keep_temp", default=False, type=bool, help="Keep all temp files (need more disk space).")
    args = parser.parse_args()

    start = timer()

    # 1-22, X, Y, XY (pseudo-autosomal region of X), and MT
    SNP_list = args.list # chr,SNP_id,BP
    reduce_folder = PurePath(args.folder).as_posix()
    patient_data_file = args.file
    patient_data_file_prefix = Path(patient_data_file).stem
    patient_data_id = args.id
    phenotype = [args.pheno, args.type] # category, continuous

    function_list = [[True,False],
                    True,True,True,True,True,
                    True,True,True,True,True,
                    True,True,True]
    
    if args.check_indel:
        function_list[0][1] == True

    if args.only_check_indel:
            function_list[0][0] == False

    ID_filter_file_raw = PurePath("Selected", phenotype[0], patient_data_file_prefix + r"_list_" + phenotype[0] + r"_raw.csv").as_posix()
    if function_list[1]:
        print("Step 1: generate patient ID filter.")
        Path(Path(ID_filter_file_raw).parent.absolute()).mkdir(parents=True, exist_ok=True)
        generate_ID_filter(patient_data_file, patient_data_id, phenotype, ID_filter_file_raw)

    temp_folder = PurePath("Selected", phenotype[0], "temp").as_posix()
    if function_list[2]:
        print("Step 2: select SNP files by patient ID filter.")
        Path(temp_folder).mkdir(parents=True, exist_ok=True)
        select_reduce_by_id(reduce_folder, ID_filter_file_raw, temp_folder)

    merged_folder = PurePath("Selected", phenotype[0], "merged").as_posix()
    if function_list[3]:
        print("Step 3: merge selected SNP files.")
        Path(merged_folder).mkdir(parents=True, exist_ok=True)
        merge_temp_file(temp_folder, merged_folder, SNP_list, 800)
        print("Removing temp folder...")
        if not args.keep_temp:
            shutil.rmtree(temp_folder)

    filted_folder = PurePath("Selected", phenotype[0], "filted").as_posix()
    if function_list[4]:
        print("Step 4: remove SNPs not in range.")
        Path(filted_folder).mkdir(parents=True, exist_ok=True)
        remove_out_of_range_SNP(merged_folder, filted_folder)

    common_SNP_file = PurePath("Selected", phenotype[0], patient_data_file_prefix + r"_common_SNP_" + phenotype[0] + r".csv").as_posix()
    out_range_SNP_file = PurePath("Selected", phenotype[0], patient_data_file_prefix + r"_out_range_uncommon_SNP_" + phenotype[0] + r".csv").as_posix()

    if function_list[5]:
        print("Step 5: make common SNPs list.")
        if function_list[0][0]:
            intersection_of_SNP(filted_folder, common_SNP_file)
        if function_list[0][1]:
            out_range_out_intersection_SNP(SNP_list, common_SNP_file, out_range_SNP_file)

    patient_common_SNP_file = PurePath("Selected", phenotype[0], patient_data_file_prefix + r"_common_SNP_" + phenotype[0] + r".csv").as_posix()
    patient_outrange_SNP_file = PurePath("Selected", phenotype[0], patient_data_file_prefix + r"_outrange_SNP_" + phenotype[0] + r".csv").as_posix()
    if function_list[6]:
        print("Step 6: merge filted SNP files with common SNPs list.")
        if function_list[0][0]:
            filted_SNP_merge(filted_folder, common_SNP_file, patient_common_SNP_file)
            if not args.keep_temp:
                print("Removing filted folder...")
                shutil.rmtree(filted_folder)
        if function_list[0][1]:
            filted_SNP_merge(merged_folder, out_range_SNP_file, patient_outrange_SNP_file)

    divers_snp_map = PurePath("Selected", phenotype[0], patient_data_file_prefix + r"_divers_SNP_map_" + phenotype[0] + r".csv").as_posix()
    outrange_divers_snp_map = PurePath("Selected", phenotype[0], patient_data_file_prefix + r"_divers_outrange_SNP_map_" + phenotype[0] + r".csv").as_posix()
    if function_list[7]:
        if function_list[0][0]:
            generate_divers_SNP_filter(patient_common_SNP_file, SNP_list, divers_snp_map)
        if function_list[0][1]:
            generate_divers_SNP_filter(patient_outrange_SNP_file, SNP_list, outrange_divers_snp_map)

    patient_divers_SNP_file = PurePath("Selected", phenotype[0], patient_data_file_prefix + r"_divers_SNP_" + phenotype[0] + r".csv").as_posix()
    patient_divers_outrange_SNP_file = PurePath("Selected", phenotype[0], patient_data_file_prefix + r"_divers_outrange_SNP_" + phenotype[0] + r".csv").as_posix()
    if function_list[8]:
        if function_list[0][0]:
            divers_SNP_merge(patient_common_SNP_file, divers_snp_map, patient_divers_SNP_file)
        if function_list[0][1]:
            divers_SNP_merge(patient_outrange_SNP_file, outrange_divers_snp_map, patient_divers_outrange_SNP_file)

    map_file = PurePath(phenotype[0], patient_data_file_prefix + r"_" + phenotype[0] + r".map").as_posix()
    outrange_map_file = PurePath(phenotype[0], patient_data_file_prefix + r"_" + phenotype[0] + r"_outrange.map").as_posix()
    if function_list[9]:
        Path(Path(map_file).parent.absolute()).mkdir(parents=True, exist_ok=True)
        if function_list[0][0]:
            make_map_file(patient_divers_SNP_file, map_file)
        if function_list[0][1]:
            make_map_file(patient_divers_outrange_SNP_file, outrange_map_file)

    ped_file= PurePath(phenotype[0], patient_data_file_prefix + r"_" + phenotype[0] + r".ped").as_posix()
    outrange_ped_file= PurePath(phenotype[0], patient_data_file_prefix + r"_" + phenotype[0] + r"_outrange.ped").as_posix()
    if function_list[10]:
        Path(Path(ped_file).parent.absolute()).mkdir(parents=True, exist_ok=True)
        if function_list[0][0]:
            make_ped_file(patient_divers_SNP_file, ID_filter_file_raw, ped_file)
        if function_list[0][1]:
            make_ped_file(patient_divers_outrange_SNP_file, ID_filter_file_raw, outrange_ped_file)

    map_ped_file_prefix = Path(ped_file).stem
    # make bed. bim. fam
    if function_list[11]:
        if function_list[0][0]:
            plink_string1 = ".\plink --file .\\" + phenotype[0] + "\\" + map_ped_file_prefix + " --make-bed --out .\\" + phenotype[0] + "\\" + map_ped_file_prefix
            os.system(plink_string1)
        if function_list[0][1]:
            plink_string1 = ".\plink --file .\\" + phenotype[0] + "\\" + map_ped_file_prefix + "_outrange --make-bed --out .\\" + phenotype[0] + "\\" + map_ped_file_prefix + "_outrange"
            os.system(plink_string1)
        # freqx
        plink_freq = ".\plink --bfile .\\" + phenotype[0] + "\\" + map_ped_file_prefix + " --freqx --out .\\" + phenotype[0] + "\\" + map_ped_file_prefix
        os.system(plink_freq)
            
    # QC
    if args.type == "category":
        assoc_analysis = "--logistic"
    else:
        assoc_analysis = "--linear"
    
    if function_list[12]:
        if function_list[0][0]:
            plink_string2 = ".\plink --bfile .\\" + phenotype[0] + "\\" + map_ped_file_prefix + " --not-chr 0 x y xy --maf 0.05 --hwe 0.000001 --geno 0.05 --mind 0.05 " + assoc_analysis + " --ci 0.95 --out .\\" + phenotype[0] + "\\" + map_ped_file_prefix + "_QC"
        if function_list[0][1]:
            plink_string2 = ".\plink --bfile .\\" + phenotype[0] + "\\" + map_ped_file_prefix + "_outrange --not-chr 0 x y xy --maf 0.05 --hwe 0.000001 --geno 0.05 --mind 0.05 " + assoc_analysis + " --ci 0.95 --out .\\" + phenotype[0] + "\\" + map_ped_file_prefix + "_outrange_QC"
        os.system(plink_string2)

    # QC + prune
    if function_list[13]:
        if function_list[0][0]:
            plink_string3 = ".\plink --bfile .\\" + phenotype[0] + "\\" + map_ped_file_prefix + " --not-chr 0 x y xy --maf 0.05 --hwe 0.000001 --geno 0.05 --mind 0.05 --indep-pairwise 50 5 0.2 " + assoc_analysis + " --ci 0.95 --out .\\" + phenotype[0] + "\\" + map_ped_file_prefix + "_QC_prune"
            os.system(plink_string3)
        if function_list[0][1]:
            plink_string3 = ".\plink --bfile .\\" + phenotype[0] + "\\" + map_ped_file_prefix + "_outrange --not-chr 0 x y xy --maf 0.05 --hwe 0.000001 --geno 0.05 --mind 0.05 --indep-pairwise 50 5 0.2 " + assoc_analysis + " --ci 0.95 --out .\\" + phenotype[0] + "\\" + map_ped_file_prefix + "_outrange_QC_prune"
            os.system(plink_string3)

    end = timer()
    print("All tasks were done! This round spent " + str(timedelta(seconds=end-start)).split(".")[0] + ".")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(e)