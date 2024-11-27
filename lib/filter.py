from pathlib import Path
import pandas as pd
import gc
import os

def snp_map_filter(SNP_map_file: Path, new_SNP_map_file: Path, new_SNP_map_range_file: Path):
    print("Opening snp map file...")
    snp_df = pd.read_csv(SNP_map_file)
    snp_df['chr'] = snp_df['chr'].astype(object)
    s = pd.to_numeric(snp_df['chr'],errors='coerce')
    snp_df.loc[s.notna(),'chr'] = s.dropna().astype(int)
    snp_df.sort_values(by=["chr","BP"], inplace=True, ascending=True)
    snp_df.reset_index(drop=True, inplace=True)
    print("Saving sorted snp map file...")
    snp_df.to_csv(new_SNP_map_file, index=False)
    result = (
        snp_df.groupby("chr")
        .apply(lambda group: pd.Series({"start": group.index[0], "end": group.index[-1]}))
        .reset_index()
    )
    result.to_csv(new_SNP_map_range_file, index=False)

def dataset_ID_filter(patient_file: Path, patient_id: str, outputfilename: Path):
    # function to return the file extension
    print("Opening patient's file...")
    file_extension = Path(patient_file).suffix
    if file_extension == '.xlsx':
        patient_df = pd.read_excel(patient_file)
    elif file_extension == '.sav':
        patient_df = pd.read_spss(patient_file, convert_categoricals=False)
    elif file_extension == '.csv':
        patient_df = pd.read_csv(patient_file)
    
    new_patient_df = pd.DataFrame()
    new_patient_df['ID'] = patient_df[patient_id]
    
    # Create a mapping from uppercase column names to the original column names
    uppercase_to_original_col = {col.upper(): col for col in patient_df.columns}
    
    # Find matching columns using the mapping
    gender_verb_list = ['GENDER', 'SEX', 'S']
    male_verb_list = ['MALE', 'M']
    matching_gender_column = [uppercase_to_original_col[uc_col] for uc_col in gender_verb_list if uc_col in uppercase_to_original_col][0]

    if not matching_gender_column:
        print('Cant find the gender column!')
        os._exit(1)
    else:
        # Male = 1, Female = 2
        list_gender = patient_df[matching_gender_column].dropna().unique()
        # if gender in numeric
        if list_gender.max() == 1: # gender is in 0/1 form
            print("The gender is in 0/1 form, changing to Male = 1, Female = 2. Fill null = -9.")
            new_patient_df['S'] = patient_df[matching_gender_column].apply(lambda x: 2 if x == 0 else (1 if x == 1 else -9))
        elif list_gender.max() == 2: # female is already 2
            print("The gender is in mixed form, keep as original form. Fill null = -9.")
            new_patient_df['S'] = patient_df[matching_gender_column]
            new_patient_df['S'] = new_patient_df['S'].fillna(-9)
        else: # gender in verb
            print("The gender is string-like, changing to Male = 1, Female = 2.")
            new_patient_df['Gender_upper'] = patient_df[matching_gender_column].astype(str).str.upper()
            new_patient_df['S'] = new_patient_df['Gender_upper'].apply(lambda x: 2 if x in male_verb_list else (-9 if pd.isna(x) else 1))
            new_patient_df.drop('Gender_upper', axis=1, inplace=True)

    # *Labels: ID S
    new_patient_df = new_patient_df.astype({"ID": str, "S": int})
    new_patient_df.drop_duplicates(subset=["ID"], inplace=True)
    print("Exporting file...")
    new_patient_df.to_csv(outputfilename, index=False)
    del new_patient_df
    gc.collect()

def pheno_ID_filter(patient_file: Path, patient_id: str, phenotype: list, outputfilename: Path):
    # function to return the file extension
    print("Opening patient's file...")
    file_extension = Path(patient_file).suffix
    if file_extension == '.xlsx':
        patient_df = pd.read_excel(patient_file)
    elif file_extension == '.sav':
        patient_df = pd.read_spss(patient_file, convert_categoricals=False)
    elif file_extension == '.csv':
        patient_df = pd.read_csv(patient_file)
    
    new_patient_df = pd.DataFrame()
    new_patient_df['ID'] = patient_df[patient_id]

    case_verb_list = ['CASE', 'YES', 'Y', '1']
    if phenotype[0] not in patient_df.columns:
        print('Cant find the phenotype in column! Current phenotype: ' + phenotype[0])
        os._exit(1)
    else:
        new_patient_df[phenotype[0]] = patient_df[phenotype[0]]
        if phenotype[1] == "category":
            new_patient_df[phenotype[0]] = new_patient_df[phenotype[0]].astype(float)
            list_pheno = new_patient_df[phenotype[0]].dropna().unique()
            if len(list_pheno) > 2: # if phenotype is more than 2
                print("The control/case group is not binary (eg 0/1), stop program.")
                os._exit(1)
            else:
                # control = 1, case = 2
                if list_pheno.max() == 1: 
                    # if phenotype in numeric and 0/1 form
                    print("Change the phenotype 0/1 to control = 1, case = 2, missing = -9.")
                    new_patient_df['P'] = new_patient_df[phenotype[0]].apply(lambda x: 2 if x == 1 else (1 if x == 0 else -9))
                elif list_pheno.max() == 2:
                    print("The control/case group might be 1/2, keep original serial.")
                    new_patient_df['P'] = new_patient_df[phenotype[0]]
                else: 
                    # if phenotype in verb
                    print("The phenotype is string-like, changing the phenotype to control = 1, case = 2.")
                    new_patient_df['Pheno_upper'] = new_patient_df[phenotype[0]].str.upper()
                    new_patient_df['P'] = new_patient_df['Pheno_upper'].isin(case_verb_list).replace({True: 2, False: 1})
                    new_patient_df.drop('Pheno_upper', axis=1, inplace=True)
        elif phenotype[1] == "continuous":  
            # if the phenotype is continuous
            new_patient_df['P'] = new_patient_df[phenotype[0]]
        else:
            print("Wrong argument of type of phenotype!! Must be category or continuous!!")
            os._exit(1)
    
    # fill phenotype = -9 if missing
    new_patient_df['P'] = new_patient_df['P'].fillna(-9)
    # *Labels: FID ID P
    print("Adding columns...")
    new_patient_df.insert(0, "FID", 0)
    new_patient_df.drop(phenotype[0], axis=1, inplace=True)    
    new_patient_df.drop_duplicates(subset=["ID"], inplace=True)
    print("Exporting file...")
    new_patient_df.to_csv(outputfilename, sep="\t", index=False, header=False)
    del patient_df, new_patient_df
    gc.collect()