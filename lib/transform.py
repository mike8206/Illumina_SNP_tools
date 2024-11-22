import pandas as pd
import gc

def make_map_file(patient_divers_SNP_file: str, map_file: str):
    print("Loading file...")
    map_df = pd.read_csv(patient_divers_SNP_file, usecols=['SNP_id', 'chr', 'BP'], dtype={'SNP_id': str, 'chr': str, 'BP': str})

    # remove chr or BP == '0'
    print("Removing chr or BP == 0")
    map_df = map_df.query('chr != "0" & BP != "0"')

    print("Found " + str(len(map_df)) + " SNPs in the map!")

    # *Labels: Chr SNP SNP_Position Base-Pair Coordinate
    map_df.columns = ['SNP', 'Chr', 'Base-Pair Coordinate']
    map_df.loc[:, 'SNP_Position'] = ""
    map_df = map_df[['Chr', 'SNP', 'SNP_Position', 'Base-Pair Coordinate']]

    print('Exporting file: ' + map_file)
    map_df.to_csv(map_file, sep="\t", index=False, header=False)
    del map_df
    gc.collect()
    print('Finished exporting!')

# the mapping function
def transform_row(row):
    # Mapping dictionary
    mapping_dict = {
        '0': '0 0',
        '1': 'A A', '2': 'A T', '3': 'A C', '4': 'A G',
        '5': 'T A', '6': 'T T', '7': 'T C', '8': 'T G',
        '9': 'C A', '10': 'C T', '11': 'C C', '12': 'C G',
        '13': 'G A', '14': 'G T', '15': 'G C', '16': 'G G',
        '17': 'I I', '18': 'I D', '19': 'D I', '20': 'D D',
        '21': '0 0'
        }
    # Map each element in the row using the dictionary, make sure to convert to string if necessary
    return ' '.join(row.apply(lambda x: mapping_dict.get(str(int(x)))))

def make_ped_file(patient_divers_SNP_file, ID_filter_file_raw, ped_file):
    print("Loading file...")
    post_T_chunks = []
    chunk_num = 1
    for chunk in pd.read_csv(patient_divers_SNP_file, chunksize = 50000, dtype={'SNP_id': str, 'chr': str, 'BP': str}):
        print("Working on chunk: " + str(chunk_num) + ". Removing chr or BP == 0 ...")
        chunk = chunk.query('chr != "0" & BP != "0"')
        print("Transposing columns and rows...")
        chunk_T = chunk.iloc[:, 3:].transpose()
        del chunk
        print("Applying transformer...")
        chunk_T.reset_index(inplace=True)
        chunk_T['Combined'] = chunk_T.iloc[:, 1:].apply(transform_row, axis=1)
        chunk_T = chunk_T[['index', 'Combined']] # Remain necessary columns
        chunk_T.columns = ['ID', 'Genotype']  # Rename for clarity
        post_T_chunks.append(chunk_T['Genotype'])
        chunk_num = chunk_num + 1
        gc.collect()

    # Concatenating all the transposed chunks together
    print("concate chunks...")
    df_transposed = pd.concat(post_T_chunks, axis=1)
    del post_T_chunks
    # Reset the index of the transposed DataFrame
    df_transposed.reset_index(inplace=True, drop=True)

    # insert ID from last chunk_T
    df_transposed.insert(0, "ID", chunk_T['ID'])
    del chunk_T
    gc.collect()

    # Combine transposed data into 'Combined'
    df_transposed['Combined'] = df_transposed.iloc[:, 1:].apply(lambda row: ' '.join(row.astype(str)), axis=1)
    df_transposed = df_transposed[['ID', 'Combined']] # Remain necessary columns
    df_transposed.columns = ['ID', 'Genotype'] # Rename for clarity

    # *Labels: FID ID F M S P Genotype
    print("Adding columns...")
    # df_transposed['S'] = df_transposed['ID'].apply(lambda x: extract_first_numeric(x)) # Apply the function to create the new column 'S'
    df_transposed.loc[:, 'FID'] = 0
    df_transposed.loc[:, 'F'] = 0
    df_transposed.loc[:, 'M'] = 0
        
    print("Re-arranging...")
    df_transposed = df_transposed[['FID', 'ID', 'F', 'M', 'Genotype']]

    pt_df_raw = pd.read_csv(ID_filter_file_raw) # ID, S, P
    final_df = df_transposed.merge(pt_df_raw, how='left', on='ID')
    del pt_df_raw, df_transposed
    gc.collect()

    print('Exporting file: ' + ped_file)
    final_df = final_df[['FID', 'ID', 'F', 'M', 'S', 'P', 'Genotype']]
    final_df.columns = ['FID', 'ID', 'F', 'M', 'S', 'P', 'Genotype'] 
    final_df.to_csv(ped_file, sep="\t", index=False, header=False)
    del final_df
    gc.collect()
    print('Finished exporting!')