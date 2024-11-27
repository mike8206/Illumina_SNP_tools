from pathlib import Path
import pandas as pd
import gc

def select_reduce_by_id(reduce_folder: Path, filter_file: Path, temp_folder: Path):
    filter_df = pd.read_csv(filter_file)
    id_list = filter_df['ID'].tolist()
    print("There are " + str(len(id_list)) + " cases in the ID filter.")

    files = [str(f) for f in Path(reduce_folder).iterdir() if f.match("*.parquet")]
    files.sort()
    print("There are " + str(len(files)) + " files in the folder.")

    selected_id_df = pd.DataFrame()
    num = 1
    zfill_num = len(str(len(files)))

    for i in range(len(files)):
        df = pd.read_parquet(files[i])
        try:
            matched_id_df = df[df.columns.intersection(id_list)]
            if len(matched_id_df.columns) > 0:
                # Add back column SNP_id
                matched_id_df.insert(0, "SNP_id", df['SNPs'])
            
                # remove matched ID from id_list
                ids_to_remove = matched_id_df.columns.to_list()
                id_list = [item for item in id_list if item not in ids_to_remove]
                print("Found "+ str(len(matched_id_df.columns)) + " cases in file: " + Path(files[i]).stem +". Saving temp file " + str(num).zfill(zfill_num) + " ... There are " + str(len(id_list)) + " cases not found yet.")
                matched_id_df.to_parquet(str(temp_folder) + r"/temp_" + str(num).zfill(zfill_num) + r".parquet", index=False)
                num += 1
        except KeyError:
            continue
        
    del filter_df, id_list, df, matched_id_df, selected_id_df
    gc.collect()

def selected_chr_SNP_map(selected_chr: str, snp_map_file: Path, selected_chr_file: Path):
    SNP_df = pd.read_csv(snp_map_file, low_memory=False)
    selected_chr_df = SNP_df[SNP_df['chr'] == str(selected_chr)]
    selected_chr_df.to_csv(selected_chr_file, index=False)
    del SNP_df, selected_chr_df
    gc.collect()

def select_top_list(top_list_file: Path):
    top_list_df = pd.read_csv(top_list_file)
    extract_snp_list = top_list_df['SNP'].to_list()
    result = ",".join(map(str, extract_snp_list))
    return result