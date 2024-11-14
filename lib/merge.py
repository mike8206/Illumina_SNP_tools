from pathlib import Path
import pandas as pd
import gc

def merge_temp_file(temp_folder: str, merged_folder: str, SNP_list_file: str, case_num: int):
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
        result_df = pd.concat([result_df, temp_df], axis=1)
        if i < len(temp_files)-1:
            if len(result_df.columns) > case_num:
                print("Collect cases number: " + str(len(result_df.columns)) + ", generating merged file "+ str(num).zfill(zfill_num) +"...")
                Path(merged_folder).mkdir(parents=True, exist_ok=True)
                result_df.to_csv(merged_folder + r'/merged_' + str(num).zfill(zfill_num) +r'_' + str(len(result_df.columns)) + r'.csv', index=False)
                num += 1
                result_df = pd.read_csv(SNP_list_file, usecols=['SNP_id'])
        else:
            print("Finishing... Collect cases number: " + str(len(result_df.columns)) + ", generating merged file "+ str(num).zfill(zfill_num) +"...")
            result_df.to_csv(merged_folder + r'/merged_' + str(num).zfill(zfill_num) +r'_' + str(len(result_df.columns)) + r'.csv', index=False)
            del temp_df, result_df
            gc.collect()
    print("Done!")

def filted_SNP_merge(filted_folder: str, input_SNP_file: str, outputfilename: str):
    file_list = [str(f) for f in Path(filted_folder).iterdir() if f.match("*.csv")]
    file_list.sort()
    print("There are " + str(len(file_list)) + " filted files in the folder.")

    merged_df = pd.read_csv(input_SNP_file)
    filter_list = merged_df['SNP_id'].to_list()
    print('There are ' + str(len(filter_list)) + " SNPs in the list.")

    for i in range(len(file_list)):
        print('Loading file: '+file_list[i])
        df = pd.read_csv(file_list[i], usecols=['SNP_id'])
        filtered_df = df[df['SNP_id'].isin(filter_list)]
        index_list = list(filtered_df.index)
        del df, filtered_df
        gc.collect()
        print('There are ' + str(len(index_list)) + ' found in this file.')

        # Read and process the CSV file in chunks
        print('Processing file in chunks...')
        filtered_chunks = []
        for chunk in pd.read_csv(file_list[i], chunksize=30000):
            filtered_chunk = chunk[chunk.index.isin(index_list)]
            filtered_chunks.append(filtered_chunk)
            del chunk
        filtered_df = pd.concat(filtered_chunks, ignore_index=True)
        del filtered_chunks
        gc.collect()
        filtered_df.set_index('SNP_id', inplace=True)
        filtered_df = filtered_df.reindex(filter_list)
        filtered_df.reset_index(inplace=True)
        row_count, col_count = filtered_df.shape
        print('The size of df is ' + str(row_count) + ' by ' + str(col_count) + ', merging file: ' + file_list[i])

        merged_df = pd.concat([merged_df, filtered_df.iloc[:, 1:]], axis=1)
        del filtered_df
        gc.collect()

    print('Exporting file: '+outputfilename)
    merged_df.to_csv(outputfilename, index=False)
    del merged_df
    gc.collect()
    print('Done!')

def divers_SNP_merge(patient_common_SNP_file: str, divers_snp_map: str, patient_divers_SNP_file: str):
    # Merge with original file
    print("Merge with original file...")
    divers_snp_map_df = pd.read_csv(divers_snp_map, usecols=['SNP_id'])
    divers_filter_list = divers_snp_map_df['SNP_id'].to_list()
    print('There are ' + str(len(divers_filter_list)) + " SNPs in the list.")
    del divers_snp_map_df
    gc.collect()
    
    print("Loading the merged file...")
    filtered_chunks = []
    for chunk in pd.read_csv(patient_common_SNP_file, chunksize=30000):
        print('Processing file in chunks...')
        filtered_chunk = chunk[chunk['SNP_id'].isin(divers_filter_list)]
        filtered_chunks.append(filtered_chunk)
        del chunk
    final_df = pd.concat(filtered_chunks, ignore_index=True)
    del filtered_chunks
    gc.collect()
    final_df.reset_index(drop=True, inplace=True)
    row_count, col_count = final_df.shape
    print('The size of the df is ' + str(row_count) + ' by ' + str(col_count) + '.')

    # Insert chr and BP from SNP map
    print("Inserting SNP map data...")
    divers_snp_map_df = pd.read_csv(divers_snp_map, usecols=['chr', 'BP'])
    final_df.insert(1, "chr", divers_snp_map_df['chr'])
    final_df.insert(2, "BP", divers_snp_map_df['BP'])
    del divers_snp_map_df
    gc.collect()

    print("Sorting dataframe by Chr and BP")
    final_df.sort_values(by=['chr','BP'], inplace=True)
    final_df.reset_index(drop=True, inplace=True)

    print('Exporting file: '+ patient_divers_SNP_file)
    final_df.to_csv(patient_divers_SNP_file, index=False)
    del final_df
    gc.collect()
    print('Done!')