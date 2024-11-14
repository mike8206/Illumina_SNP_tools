library(data.table)
library(tidyverse)
library(qqman)

getwd()

setwd("D:/Nextcloud/SNP/SNP_data/MetSf") ## You will need to change this in your own scripts

# Read in file
example_basic_logistic_GWAS_results<-read.table(file="SNP_manucript20221122_MetSf_QC.assoc.linear", 
                                                header=TRUE, na="NA")
glimpse(example_basic_logistic_GWAS_results)

### Manhattan plot (qqman version)
manhattan(example_basic_logistic_GWAS_results %>% na.omit(), 
          main = "Manhattan Plot of Metabolic factors", cex.axis = 0.7,
          chr = "CHR", bp = "BP", p = "P", chrlabs = c(1:22),
          col = c("blue4", "orange3"), annotatePval = 0.00001)

## Q-Q plot (qqman version)
qq(example_basic_logistic_GWAS_results$P, 
   main = "Q-Q plot of SNP GWAS p-values", xlim = c(0, 7), ylim = c(0,12), 
   pch = 18, col = "blue4", cex = 1.5, las = 1)

top_results_linear<-example_basic_logistic_GWAS_results[order(example_basic_logistic_GWAS_results$P),]
fwrite(top_results_linear[1:25,], "Top_25_SNPs.csv", sep = ",")

# .\plink --file .\Lean_HL\SNP_5000_HL --r2 --ld-snp rs12592963 --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0