from pathlib import Path
import pandas as pd
import gc

def merge_temp_file(temp_folder: Path, merged_folder: Path, SNP_list_file: Path, case_num: int):
    temp_files = [str(f) for f in Path(temp_folder).iterdir() if f.match("*.parquet")]
    temp_files.sort()
    num = 1
    zfill_num = len(str(len(temp_files)))
    print("There are " + str(len(temp_files)) + " temp files in the folder.")

    # Add SNP map (first column)
    result_df = pd.read_csv(SNP_list_file, usecols=['SNP_id'])

    # Add columns with loop
    for i in range(len(temp_files)):
        print("Merging file: "+ temp_files[i])
        temp_df = pd.read_parquet(temp_files[i])
        result_df = pd.merge(result_df, temp_df, on="SNP_id")
        if i < len(temp_files)-1:
            if len(result_df.columns) > case_num:
                print("Collect cases number: " + str(len(result_df.columns)) + ", generating merged file "+ str(num).zfill(zfill_num) +"...")
                Path(merged_folder).mkdir(parents=True, exist_ok=True)
                result_df.to_csv(str(merged_folder) + r'/merged_' + str(num).zfill(zfill_num) + r'_' + str(len(result_df.columns)) + r'.csv', index=False)
                num += 1
                result_df = pd.read_csv(SNP_list_file, usecols=['SNP_id'])
        else:
            print("Finishing... Collect cases number: " + str(len(result_df.columns)) + ", generating merged file "+ str(num).zfill(zfill_num) +"...")
            result_df.to_csv(str(merged_folder) + r'/merged_' + str(num).zfill(zfill_num) + r'_' + str(len(result_df.columns)) + r'.csv', index=False)
            del temp_df, result_df
            gc.collect()

def merge_selected_SNP_by_range(filted_folder: Path, chr:str, input_snp_range_file: Path, chr_map_file: Path, outputfilename: Path):
    file_list = [str(f) for f in Path(filted_folder).iterdir() if f.match("*.csv")]
    file_list.sort()
    print("There are " + str(len(file_list)) + " filted files in the folder.")

    snp_range_df = pd.read_csv(input_snp_range_file)
    selected_chr = snp_range_df[snp_range_df["chr"] == chr]
    start, end = selected_chr["start"].values[0], selected_chr["end"].values[0]
    del snp_range_df, selected_chr
    gc.collect()

    for i in range(len(file_list)):
        print('Loading file: '+ file_list[i])
        temp_df = pd.read_csv(file_list[i], skiprows=range(1, start + 1), nrows = end - start + 1)
        row_count, col_count = temp_df.shape
        print('The size of df is ' + str(row_count) + ' by ' + str(col_count) + ', merging file: ' + file_list[i])
        if i == 0:
            merged_df = temp_df
        else:
            merged_df = pd.merge(merged_df, temp_df, on="SNP_id")
            del temp_df
            gc.collect()

    # Insert chr and BP from SNP map
    print("Inserting SNP map data...")
    snp_map_df = pd.read_csv(chr_map_file, usecols=['chr', 'BP'])
    merged_df.insert(1, "chr", snp_map_df['chr'])
    merged_df.insert(2, "BP", snp_map_df['BP'])

    row_count, col_count = merged_df.shape
    print('The size of merged df is ' + str(row_count) + ' by ' + str(col_count) + ', exporting file: ' + str(outputfilename))
    merged_df.to_csv(outputfilename, index=False)
    del snp_map_df, merged_df
    gc.collect()