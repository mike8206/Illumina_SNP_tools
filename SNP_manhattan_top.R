#!/usr/bin/env Rscript
library("optparse")
library("data.table")
library("qqman")

option_list = list(
  make_option(c("-d", "--dataset"), type="character", default=NULL, 
              help="dataset name", metavar="character"),
  make_option(c("-p", "--phenotype"), type="character", default=NULL, 
              help="name of phenotype", metavar="character"),
  make_option(c("-t", "--type"), type="character", default=NULL, 
              help="type of phenotype", metavar="character"),
  make_option(c("-m", "--manhattan"), type="character", default=NULL, 
              help="phenotype name of manhattan plot", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (length(opt) == 1){
  print_help(opt_parser)
  stop("At least one argument must be supplied.n", call.=FALSE)
}

# setup general variables
dataset <- opt$dataset
phenotype <- opt$phenotype
type <- opt$type
if (type == "category") {
  assoc_analysis <- "logistic"
} else {
  assoc_analysis <- "linear"
}
if (is.null(opt$manhattan)) {
  manhattan_name <- phenotype
}

dataset_folder <- paste(dataset, phenotype, sep = "/")
filename <- paste(paste(dataset, phenotype, "merged_QC_prune", sep = "_"), "assoc", assoc_analysis , sep = ".")
top_filename <- paste(paste(dataset, phenotype, "Top_SNPs", sep = "_"), "csv", sep = ".")
manhattan_title <- gsub("_", " ", paste("Manhattan Plot of", manhattan_name, sep = " "))
manhattan_filename <- paste(paste(dataset, phenotype, "plot", "manhattan", sep = "_"), "tiff", sep = ".")
QQ_title <- "Q-Q plot of SNP GWAS p-values"
QQ_filename <- paste(paste(dataset, phenotype, "plot", "qq", sep = "_"), "tiff", sep = ".")

getwd()
setwd(dataset_folder) ## You will need to change this in your own scripts

# Read in file
example_basic_logistic_GWAS_results <- read.table(file=filename, header=TRUE, na="NA")
# glimpse(example_basic_logistic_GWAS_results)

### Manhattan plot (qqman version)
tiff(filename = manhattan_filename, width = 1080, height = 608, units = "px", pointsize = 18)
manhattan(example_basic_logistic_GWAS_results, 
          main = manhattan_title, cex.axis = 0.7,
          chr = "CHR", bp = "BP", p = "P", chrlabs = c(1:22, "MT"),
          col = c("blue4", "orange3"), annotatePval = 0.00001)
dev.off()

## Q-Q plot (qqman version)
tiff(filename = QQ_filename, width = 1080, height = 608, units = "px", pointsize = 24)
qq(example_basic_logistic_GWAS_results$P, 
   main = QQ_title, xlim = c(0, 7), ylim = c(0,12), 
   pch = 18, col = "blue4", cex = 1.5, las = 1)
dev.off()

top_results_linear <- example_basic_logistic_GWAS_results[order(example_basic_logistic_GWAS_results$P),]
file_df <- top_results_linear[top_results_linear$P <= 0.00001, ] 
fwrite(file_df, top_filename, sep = ",")