from pathlib import Path
import pandas as pd
import gc
import os

def generate_ID_filter(patient_file: str, patient_id: str, phenotype: list, outputfilename: str):
    
    gender_verb_list = ['GENDER', 'SEX', 'S']
    male_verb_list = ['MALE', 'M']
    case_verb_list = ['CASE', 'YES', 'Y', '1']

    # function to return the file extension
    file_extension = Path(patient_file).suffix

    print("Opening patient's file...")
    if file_extension == '.xlsx':
        patient_df = pd.read_excel(patient_file)
    elif file_extension == '.sav':
        patient_df = pd.read_spss(patient_file, convert_categoricals=False)
    elif file_extension == '.csv':
        patient_df = pd.read_csv(patient_file)

    # Create a mapping from uppercase column names to the original column names
    uppercase_to_original_col = {col.upper(): col for col in patient_df.columns}

    new_patient_df = pd.DataFrame()
    new_patient_df['ID'] = patient_df[patient_id]

    # Find matching columns using the mapping
    matching_gender_column = [uppercase_to_original_col[uc_col] for uc_col in gender_verb_list if uc_col in uppercase_to_original_col][0]

    if matching_gender_column:
        # Male = 1, Female = 2
        max_gender = patient_df[matching_gender_column].max()
        try: 
            # if gender in numeric
            if max_gender >= 2: # if female is already 2
                new_patient_df['S'] = patient_df[matching_gender_column]
                print("The gender is in mixed form, keep as original form.")
            else:
                new_patient_df['Gender_int'] = patient_df[matching_gender_column].astype(int)
                new_patient_df['S'] = new_patient_df['Gender_int'].apply(lambda x: 2 if x == 0 else 1)
                new_patient_df.drop('Gender_int', axis=1, inplace=True)
                print("The gender is in 0/1 form, changing to Male = 1, Female = 2.")
        except:
            # if gender in verb
            new_patient_df['Gender_upper'] = patient_df[matching_gender_column].astype(str).str.upper()
            new_patient_df['S'] = new_patient_df['Gender_upper'].isin(male_verb_list).replace({True: 1, False: 2})
            new_patient_df.drop('Gender_upper', axis=1, inplace=True)
            print("The gender is string-like, changing to Male = 1, Female = 2.")

        if phenotype[0] in patient_df.columns:
            new_patient_df[phenotype[0]] = patient_df[phenotype[0]]
            new_patient_df.dropna(inplace=True)
            if phenotype[1] == "category":
                new_patient_df[phenotype[0]] = new_patient_df[phenotype[0]].astype(int)
                min_pheno = new_patient_df[phenotype[0]].min()
                max_pheno = new_patient_df[phenotype[0]].max()
                # control = 1, case = 2
                try: # if phenotype in numeric
                    if max_pheno > 2:
                        # if phenotype is more than 2
                        print("The control/case group is not binary (eg 0/1), stop program.")
                        os._exit(1)
                    elif min_pheno == 0: 
                        new_patient_df['P_int'] = new_patient_df[phenotype[0]]
                        new_patient_df['P'] = new_patient_df['P_int'].apply(lambda x: 2 if x == 1 else 1)
                        new_patient_df.drop('P_int', axis=1, inplace=True)
                        print("Change the phenotype 0/1 to control = 1, case = 2.")
                    else:
                        new_patient_df['P'] = new_patient_df[phenotype[0]]
                        print("The control/case group might be 1/2, keep original serial.")
                except: # if phenotype in verb
                    new_patient_df['Pheno_upper'] = new_patient_df[phenotype[0]].str.upper()
                    new_patient_df['P'] = new_patient_df['Pheno_upper'].isin(case_verb_list).replace({True: 2, False: 1})
                    new_patient_df.drop('Pheno_upper', axis=1, inplace=True)
                    print("The phenotype is string-like, changing the phenotype to control = 1, case = 2.")
                
                new_patient_df = new_patient_df.astype({"ID": str, "S": int, "P": int})
            else: 
                # if the phenotype is continuous
                new_patient_df['P'] = patient_df[phenotype[0]]
                new_patient_df = new_patient_df.astype({"ID": str, "S": int, "P": float})
            del patient_df
            
            new_patient_df.drop(phenotype[0], axis=1, inplace=True)    
            new_patient_df.drop_duplicates(subset=["ID"], inplace=True)
            new_patient_df.to_csv(outputfilename, index=False)
            print("Done!")
            del new_patient_df
            gc.collect()
        else:
            print('Cant find the phenotype in column! Current phenotype: ' + phenotype[0])
            os._exit(1)
    else:
        print('Cant find the gender column!')
        os._exit(1)

# Filter that judge from col_2 to col_end and remove rows not meet the criteria
def filter_row(row):
    return row.iloc[1:].apply(lambda x: 1 <= x <= 16).all()

def remove_out_of_range_SNP(mergedfilefolder: str, filtedfilefolder: str):
    file_list = [str(f) for f in Path(mergedfilefolder).iterdir() if f.match("*.csv")]
    file_list.sort()
    num = 1
    zfill_num = len(str(len(file_list)))
    print("There are " + str(len(file_list)) + " merged files in the folder.")

    for i in range(len(file_list)):
        print("Open merged file: " + file_list[i])
        result_df = pd.DataFrame()
        chunk_num = 1

        for chunk in pd.read_csv(file_list[i], chunksize=100000):
            print("Working on chunk: " + str(chunk_num))
            print("Applying filter...")
            result_df_temp = chunk[chunk.apply(filter_row, axis=1)]
            result_df = pd.concat([result_df, result_df_temp], axis=0, ignore_index=True)
            del chunk, result_df_temp
            gc.collect()
            chunk_num = chunk_num+1
        print("Finish this file and start output...")
        result_df.to_csv(filtedfilefolder + r"/filted_" + str(num).zfill(zfill_num) + r".csv", index=False)
        del result_df
        gc.collect()
        num += 1
    print("Done!")

def generate_divers_SNP_filter(patient_common_SNP_file: str, SNP_list: str, divers_snp_map: str):
    # Count unique genotype
    print("Counting SNPs's unique genotype...")
    mutant_df = pd.DataFrame()
    chunk_num = 1
    mutant_chunks = []
    for chunk in pd.read_csv(patient_common_SNP_file, chunksize = 50000):
        print("Working on chunk: " + str(chunk_num))
        df_mutant_temp = chunk.iloc[:,1:].nunique(axis=1)
        mutant_chunks.append(df_mutant_temp)
        chunk_num = chunk_num + 1
        del chunk, df_mutant_temp
        gc.collect()
    mutant_df = pd.concat(mutant_chunks, ignore_index=True)
    merged_snp_df = pd.read_csv(patient_common_SNP_file, usecols=['SNP_id'])
    unique_df = pd.concat([merged_snp_df, mutant_df], axis=1)
    unique_df.columns=['SNP_id','Unique']
    del merged_snp_df, mutant_df
    gc.collect()

    print("Select SNPs not only 1 genotype...")
    # Select SNPs not only 1 genotype
    unique_df = unique_df[unique_df.Unique > 1]

    # Merge with SNP map
    print("Merge with SNP map...")
    SNP_map_df = pd.read_csv(SNP_list, dtype={'chr': str, 'SNP_id': str, 'BP': str})
    map_df = pd.merge(unique_df, SNP_map_df, on='SNP_id', how='left')
    map_df = map_df.drop(columns=['Unique'])
    map_df.to_csv(divers_snp_map, index=False)
    del unique_df, map_df, SNP_map_df
    gc.collect()

    print("Done!")