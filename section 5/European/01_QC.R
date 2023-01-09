#### Output the full gwas result of hypertension on European
library(stringr)
library(data.table)
library(tidyverse)
library(bigsnpr)
library(bigstatsr)
library(progress)
library(foreach)
library(parallel)
library(qqman)
library(HardyWeinberg)
library(SPCompute)


#####################################
#####################################
##### Step 0: Read-in the full data
#####################################
#####################################
bd <- read.table("/home/ziang/Documents/UKB/ukb47570.tab", header=TRUE, sep="\t", nrows = 3)
bd_names <- names(bd)
names_want <- grep("f.20002.0.",bd_names,value = T)
bd_want <- fread("/home/ziang/Documents/UKB/ukb47570.tab", header=TRUE, sep="\t", select = c("f.eid",names_want))
check_asthma <- function(row){
  if(any(na.omit(row) %in% 1111)) 1
  else 0
}
check_hypertension <- function(row){
  if(any(na.omit(row) %in% 1065)) 1
  else 0
}
asthma <- apply(bd_want[,-1], 1, check_asthma)
ncol_bd <- ncol(bd_want)
hypertension <- apply(bd_want[,-1], 1, check_hypertension)
pheno_data <- data.frame(ID = bd_want$f.eid, asthma = asthma, hypertension = hypertension)
names_want2 <- grep("f.21001.0.",bd_names,value = T)
BMI_data <- fread("/home/ziang/Documents/UKB/ukb47570.tab", header=TRUE, sep="\t", select = c("f.eid",names_want2))
names(BMI_data) <- c("ID", "BMI")
pheno_data <- merge(pheno_data, BMI_data, by = "ID")
names_want3 <- grep("f.21022.",bd_names,value = T)
age_data <- fread("/home/ziang/Documents/UKB/ukb47570.tab", header=TRUE, sep="\t", select = c("f.eid",names_want3))
names(age_data) <- c("ID", "age")
pheno_data <- merge(pheno_data, age_data, by = "ID")
names_want_sex <- grep("f.22001.",bd_names,value = T)
sex_data <- fread("/home/ziang/Documents/UKB/ukb47570.tab", header=TRUE, sep="\t", select = c("f.eid",names_want_sex))
names(sex_data) <- c("ID", "sex")
pheno_data <- merge(pheno_data, sex_data, by = "ID")
### automated reading of Diastolic blood pressure, initial visit
names_want4 <- grep("f.4079.0",bd_names,value = T)
bp_data <- fread("/home/ziang/Documents/UKB/ukb47570.tab", header=TRUE, sep="\t", select = c("f.eid",names_want4))
### There are two measures of bp taken two moments apart, we take their average value:
bp_data$bp <- 0.5*(bp_data$f.4079.0.0 + bp_data$f.4079.0.1)
bp_data <- bp_data %>% select("f.eid", "bp")
names(bp_data)[1] <- "ID"
pheno_data <- merge(pheno_data, bp_data, by = "ID")
save(pheno_data, file = "/home/ziang/Documents/Stats-Gene/UKB/European_GWAS_hypertension_bp/data/pheno_data.rda")


#################################################
#################################################
##### Step 1: Subset out self-reported European
#################################################
#################################################
names_want_pop <- grep("f.21000.0",bd_names,value = T)
names_want_pop_genetic_confirmed <- grep("f.22006",bd_names,value = T)
names_want_pop_kin <- grep("f.22021.0.0", bd_names,value = T)
pop_want <- fread("/home/ziang/Documents/UKB/ukb47570.tab", header=TRUE, sep="\t", select = c("f.eid",names_want_pop, names_want_pop_genetic_confirmed, names_want_pop_kin))
pop_want_European <- pop_want %>% filter(f.21000.0.0 == 1001, !is.na(f.22006.0.0))
pop_want_European <- pop_want_European %>% filter(f.22021.0.0 == 0)
European_ID <- pop_want_European$f.eid
pheno_data <- pheno_data %>% filter(ID %in% European_ID)
### Number of self-reported European individuals to start with:
nrow(pheno_data)




#####################################################################
#####################################################################
##### Step 2: Merge the plink file for all self_report_European
#####################################################################
#####################################################################
geno_path = "/home/ziang/Documents/UKB/"
plink = "/home/ziang/Documents/UKB/plink2"
plink19 = "/home/ziang/Documents/UKB/plink19"
matrix_to_write <- cbind(European_ID,European_ID)
write.table(x = matrix_to_write, file = "/home/ziang/Documents/UKB/self_report_European.txt", row.names = F, col.names = F)
result <- NULL
setwd("/home/ziang/Documents/UKB")
for (i in 1:22) {
  i_chromo <- paste0("c",i)
  i_selected_name <- paste0(paste0("selected_",i_chromo),"_European")
  system(paste0(plink,
                paste(" --bfile", i_chromo),
                " --keep self_report_European.txt",
                " --make-bed",
                paste(" --out", i_selected_name)))
}
system(paste0(plink19,
              " --bfile ", "selected_c1_European",
              " --merge-list ", "merge_list_European.txt",
              " --make-bed",
              " --out ", "self_European_merged_data"))

#####################################################################
#####################################################################
##### Step 3: Filter Indivds by MIND
#####################################################################
#####################################################################
system(paste0(plink,
              paste(" --bfile", "self_European_merged_data"),
              " --mind 0.2",
              " --make-bed ",
              paste(" --out", "self_European_filtered_data")
))
# 276682 samples (147522 females, 129160 males; 276682 founders) remaining

#####################################################################
#####################################################################
##### Step 4: Filter SNPs and record the number
#####################################################################
#####################################################################
system(paste0(plink,
              paste(" --bfile", "self_European_filtered_data"),
              " --hwe 1e-10",
              " --maf 0.01",
              " --geno 0.2",
              " --make-bed ",
              paste(" --out", "self_European_filtered_data_afterSNP_QC")
))
# 784256 variants loaded from self_European_filtered_data.bim.
# 123610 variants removed due to allele frequency threshold(s)
# 17934 variants removed due to missing genotype data
# --hwe: 41037 variants removed due to Hardy-Weinberg exact test (founders only).
# 601675 variants remaining after main filters.











