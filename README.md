# 使用前請先讀過此文件!!
## 執行環境及必要程式：
* 執行環境：Windows、記憶體至少12G、足夠之硬碟空間
* Python > 3.9 (須安裝openpyxl, pyreadstat, pandas, pandas pyarrow)
```
pip install openpyxl, pyreadstat, pandas, pandas pyarrow
```
* R > 4.3.1 (須安裝optparse, data.table, qqman)
```
install.packages(c("optparse", "data.table", "qqman"))
```
* Plink 1.9：https://www.cog-genomics.org/plink/
## 檔案結構：
* 同資料夾內須有lib資料夾、SNP_matrix_to_reduce.py、SNP_reduce_to_plink.py、SNP_manhattan_top.R、plink.exe
# From Illumina Matrix files to Reduced SNP files
## 目的：
處理Illumina的matrix檔，擷取出Genotype並縮減檔案大小，產生parquet檔方便保存及後續處理。
## 範例：
```
python SNP_matrix_to_reduce.py --matrix_folder Matrix --SNP_folder Reduced
```
## 程式執行步驟：
1. 製作SNP_folder資料夾
2. 逐一轉換matrix_folder內文字檔(*.txt)成簡化過的parquet檔
## 基本參數：
* --matrix_folder：Illumina的Matrix檔的資料夾位置
* --SNP_folder：簡化過的parquet檔的資料夾位置
## 可選參數：
* --skip_rows：Matrix檔須跳過的行數 (預設為9行)
* --chunk_size：影響第2步批次處理SNP的單位量 (預設為100000個)
# From Reduced SNP files to PLINK Binary files and Making Plots
## 目的：
處理簡化過的parquet檔，並根據輸入的dataset來進行篩選、合併、分析，最終產生一系列可進一步處理及分析的檔案。
## 範例：
```
python SNP_reduce_to_plink.py --SNP_folder Reduced --SNP_map SNP_map_final.csv --patient_file patient.xlsx --id_name ID --phenotype BMI --type continuous
```
## 程式執行步驟：
1. 製作dataset資料夾
2. 使用dataset製作篩選ID的檔案、具有表現型(phenotype)的檔案
3. 篩選簡化過的parquet檔，並將parquet檔以單位個案量合併處理
4. 將合併之檔案以chr為單位分割，並合併成以chr為單位之子檔案
5. 將子檔案與dataset內必要的資料合併，製作以dataset為主、chr為單位之MAP檔與PED檔(此時不具有表現型)
6. 使用PLINK 1.9將以chr為單位之MAP檔、PED檔與具有表現型的檔案合併為具有表現型的binary檔(BED)
7. 使用PLINK 1.9將以chr為單位、具有表現型之binary檔(BED)合併
8. 使用PLINK 1.9將binary檔清理 (參數為--not-chr 0 x y xy --maf 0.05 --hwe 0.000001 --geno 0.05 --mind 0.05 --indep-pairwise 50 5 0.2 --ci 0.95)
9. 使用R來製作Manhattan plot、Q-Q plot、篩選出p<10^-5的SNPs (SOI, snp of interest)
10. 使用PLINK 1.9將SOI的Genotype篩選出來，並製作具有ID、表現型及SOI的檔案
## 基本參數：
* --SNP_folder：簡化過的parquet檔的資料夾位置
* --SNP_map：具備相同SNP的檔案，須包含chr, SNP_id, BP
* --patient_file：dataset的檔案名稱(須包含副檔名)
* --id_name：dataset中個案的欄位名稱
* --phenotype：dataset中表現型的欄位名稱
* --type：表現型的類型(category或continuous)
## 可選參數：
* --chr：選定特定之chr進行分析，未給定則產生全部chr (影響第4步之後之檔案)
* --merge_size：影響第3步的單位個案合併量 (預設為超過800個)
* --chunk_size：影響第5步製作PED檔的批次量 (預設為50000個SNP)
* --keep_temp：預設為False，更改為True可不刪除第3步之暫存檔
* --make_bed：預設為True，更改為False則不執行第6步、第7步
* --prune：預設為True，更改為False則不執行第8步
* --make_plot：預設為True，更改為False則不執行第9步
* --rscript：變更R執行檔的位置 (Win 4.3.1版預設為C:\\Program Files\\R\\R-4.3.1\\bin\\Rscript)
* --rfile：變更SNP_manhattan_top.R檔位置 (預設為SNP_manhattan_top.R)
* --snp：可自訂欲篩選Genotype之SNP (預設自動根據SOI的檔案篩選Genotype)
