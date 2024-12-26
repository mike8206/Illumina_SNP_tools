from pathlib import Path
import pandas as pd
import gc

def split_transform_row(row):
    # Mapping dictionary
    mapping_dict = {
        '--': 0,
        'AA': 1, 'AT': 2, 'AC': 3, 'AG': 4,
        'TA': 5, 'TT': 6, 'TC': 7, 'TG': 8,
        'CA': 9, 'CT': 10, 'CC': 11, 'CG': 12,
        'GA': 13, 'GT': 14, 'GC': 15, 'GG': 16,
        'II': 17, 'ID': 18, 'DI': 19, 'DD': 20,
        'UU': 21}
    
    # Map each element in the row using the dictionary, make sure to convert to string if necessary
    return row.apply(lambda x: mapping_dict.get(str(x).split('|')[0]))

def make_reduce_file(matrix_folder: Path, skip_rows: int, chunk_size: int, reduce_folder: Path):
    temp_files = [str(f) for f in Path(matrix_folder).iterdir() if f.match("*.txt")]
    temp_files.sort()
    print("There are " + str(len(temp_files)) + " temp files in the matrix folder.")

    for i in range(len(temp_files)):
        post_filename = str(Path(temp_files[i]).stem) + r".parquet"
        post_map_chunks = []
        chunk_num = 1
        for chunk in pd.read_csv(temp_files[i], skiprows=skip_rows, sep="\t", chunksize=chunk_size):
            print("Working on chunk: " + str(chunk_num) + ". Applying transformer ...")
            chunk.iloc[:, 1:] = chunk.iloc[:, 1:].apply(split_transform_row, axis=1)
            post_map_chunks.append(chunk)
            chunk_num += 1
            del chunk
            gc.collect()

        df = pd.concat(post_map_chunks, axis=0)
        df.rename(columns={'Unnamed: 0': 'SNPs'}, inplace=True)
        print("Exporting reduced file: " + str(post_filename))
        df.to_parquet(str(reduce_folder) + r'/' + post_filename, index=False)
        del post_map_chunks, df
        gc.collect()

def make_map_file(patient_snp_map: Path, map_file: Path):
    print("Loading file...")
    map_df = pd.read_csv(patient_snp_map, usecols=['SNP_id', 'chr', 'BP'], dtype={'SNP_id': str, 'chr': str, 'BP': str})

    # remove chr or BP == '0'
    print("Removing chr or BP == 0")
    map_df = map_df.query('chr != "0" & BP != "0"')

    print("Found " + str(len(map_df)) + " SNPs in the map!")

    # *Labels: Chr SNP SNP_Position Base-Pair Coordinate
    map_df.columns = ['SNP', 'Chr', 'Base-Pair Coordinate']
    map_df.loc[:, 'SNP_Position'] = ""
    map_df = map_df[['Chr', 'SNP', 'SNP_Position', 'Base-Pair Coordinate']]

    print('Exporting file: ' + str(map_file))
    map_df.to_csv(map_file, sep="\t", index=False, header=False)
    del map_df
    gc.collect()

# the mapping function
def reverse_transform_row(row):
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

def make_ped_file(patient_snp_map: Path, dataseta_ID_file: Path, chunksize: int, ped_file: Path):
    print("Loading file...")
    post_T_chunks = []
    chunk_num = 1
    for chunk in pd.read_csv(patient_snp_map, chunksize = chunksize, dtype={'SNP_id': str, 'chr': str, 'BP': str}):
        print("Working on chunk: " + str(chunk_num) + ". Removing chr or BP == 0 ...")
        chunk = chunk.query('chr != "0" & BP != "0"')
        print("Transposing columns and rows...")
        chunk_T = chunk.iloc[:, 3:].transpose()
        del chunk
        print("Applying transformer...")
        chunk_T.reset_index(inplace=True)
        chunk_T['Combined'] = chunk_T.iloc[:, 1:].apply(reverse_transform_row, axis=1)
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

    print("Adding S column...")
    pt_df_raw = pd.read_csv(dataseta_ID_file) # *Labels: ID S
    df_transposed = df_transposed.merge(pt_df_raw, how='left', on='ID')
    del pt_df_raw
    gc.collect()

    # *Labels: FID ID F M S P Genotype
    print("Adding FID, F, M, P columns...")
    df_transposed.loc[:, 'FID'] = 0
    df_transposed.loc[:, 'F'] = 0
    df_transposed.loc[:, 'M'] = 0
    df_transposed.loc[:, 'P'] = -9 # default = missing (-9)
        
    print("Re-arranging...")
    df_transposed = df_transposed[['FID', 'ID', 'F', 'M', 'S', 'P', 'Genotype']]

    print('Exporting file: ' + str(ped_file))
    df_transposed.to_csv(ped_file, sep="\t", index=False, header=False)
    del df_transposed
    gc.collect()

def extract_raw_to_csv(raw_file: Path, id_name: str, pheno_name:str, output_file: Path):
    raw_df = pd.read_csv(raw_file, sep=" ")
    raw_df.drop(columns=['FID', 'SEX', 'PAT', 'MAT'], inplace=True)
    raw_df.rename(columns={'IID': id_name, 'PHENOTYPE': pheno_name}, inplace=True)
    raw_df = raw_df[raw_df[pheno_name] != -9]
    # Convert columns from the 3rd column to the end (zero-indexed, so 2 means the 3rd column)
    # raw_df.iloc[:, 2:] = raw_df.iloc[:, 2:].astype("Int64")
    raw_df.to_csv(output_file, index=False)
    del raw_df
    gc.collect()
