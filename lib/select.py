from pathlib import Path
import pandas as pd
import gc

def select_reduce_by_id(reduce_folder: str, filter_file: str, temp_folder: str):
    filter_df = pd.read_csv(filter_file)
    id_list = filter_df['ID'].tolist()
    print("There are " + str(len(id_list)) + " cases in the ID filter.")

    files = [str(f) for f in Path(reduce_folder).iterdir() if f.match("*.parquet")]
    files.sort()
    print("There are " + str(len(files)) + " files in the folder.")

    selected_id_df = pd.DataFrame()
    zfill_num = len(str(len(files)))

    for i in range(len(files)):
        df = pd.read_parquet(files[i])
        try:
            matched_id_df = df[df.columns.intersection(id_list)]
            if len(matched_id_df.columns) > 0:
                print("Found "+ str(len(matched_id_df.columns)) + " cases in file: " + Path(files[i]).stem +". Saving temp file " + str(i).zfill(zfill_num) + " ...")
                matched_id_df.to_parquet(temp_folder + r"/temp_" + str(i).zfill(zfill_num) + r".parquet", index=False)
                # remove matched ID from id_list
                ids_to_remove = matched_id_df.columns.to_list()
                id_list = [item for item in id_list if item not in ids_to_remove]
                print("There are " + str(len(id_list)) + " cases not found yet.")
        except KeyError:
            continue
        
    del filter_df, id_list, df, matched_id_df, selected_id_df
    gc.collect()
    print("Done!")

def intersection_of_SNP(filted_folder: str, outputfile: str):
    file_list = [str(f) for f in Path(filted_folder).iterdir() if f.match("*.csv")]
    file_list.sort()
    print("There are " + str(len(file_list)) + " filted files in the folder.")

    for i in range(len(file_list)):
        print('Loading file: '+ Path(file_list[i]).stem)
        df = pd.read_csv(file_list[i], usecols=['SNP_id'])
        SNPs_temp = df['SNP_id'].to_list()
        set_temp = set(SNPs_temp)
        print('Make intersection of file: '+file_list[i])
        if i < 1:
            SNPs_list = set_temp
        else:
            SNPs_list = SNPs_list.intersection(set_temp)

    common_SNPs_df = pd.DataFrame(SNPs_list)
    common_SNPs_df = common_SNPs_df.rename(columns={0: 'SNP_id'})

    print('Exporting common SNPs list file: '+ outputfile)
    common_SNPs_df.to_csv(outputfile, index=False)
    del common_SNPs_df
    gc.collect()
    print('Done!')

def out_range_out_intersection_SNP(SNP_list_file: str, common_SNP_file: str, out_range_incommon_SNP_file: str):
    SNP_df = pd.read_csv(SNP_list_file, usecols=['SNP_id'])
    set1 = set(SNP_df['SNP_id'].to_list())
    common_SNP_df = pd.read_csv(common_SNP_file, usecols=['SNP_id'])
    set2 = set(common_SNP_df['SNP_id'].to_list())
    del SNP_df, common_SNP_df
    out_range_df = pd.DataFrame(list(set1 - set2), columns=['SNP_id'])
    del set1, set2
    out_range_df.to_csv(out_range_incommon_SNP_file, index=False)
    del out_range_df
    gc.collect()
    