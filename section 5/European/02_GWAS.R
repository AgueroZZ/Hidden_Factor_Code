#### Output the full gwas result of hypertension on african
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


#####################################################################
#####################################################################
##### Step 7: GWAS for the binary Hypertension traits
#####################################################################
#####################################################################
geno_path = "/home/ziang/Documents/UKB/"
data_path = "/home/ziang/Documents/Stats-Gene/UKB/African_GWAS_hypertension_bp/data/"
plink = "/home/ziang/Documents/UKB/plink2"
plink19 = "/home/ziang/Documents/UKB/plink19"
result_path = "/home/ziang/Documents/Stats-Gene/UKB/African_GWAS_hypertension_bp/results/"
figure_path = "/home/ziang/Documents/Stats-Gene/UKB/African_GWAS_hypertension_bp/figures/"
load("/home/ziang/Documents/Stats-Gene/UKB/African_GWAS_hypertension_bp/data/final_data.rda")
colnames(final_data)[1] <- "IID"
final_data <- cbind(FID = final_data$IID ,final_data)
final_data$hypertension <- ifelse(final_data$hypertension == 1, 2,1)
final_data$sex <- ifelse(final_data$sex == 1, 2,1)
write.table(x = final_data, file = "/home/ziang/Documents/Stats-Gene/UKB/African_GWAS_hypertension_bp/data/pheno.txt", row.names = F, 
            col.names = T, quote = F, sep = "\t")
system(paste0(plink,
              " --bfile ", geno_path,"self_africa_filtered_data_afterSNP_QC",
              " --allow-no-sex ",
              " --freq ",
              " --pheno ", data_path ,"pheno.txt",
              " --pheno-name hypertension",
              " --covar ", data_path,"pheno.txt",
              " --covar-name age,sex",
              " --glm",
              " --out ", result_path, "self_africa_result_hyper"
))




#####################################################################
#####################################################################
##### Step 8: GWAS for the continuous blood pressure traits
#####################################################################
#####################################################################
system(paste0(plink,
              " --bfile ", geno_path,"self_africa_filtered_data_afterSNP_QC",
              " --allow-no-sex ",
              " --freq ",
              " --pheno ", data_path ,"pheno.txt",
              " --pheno-name bp",
              " --covar ", data_path,"pheno.txt",
              " --covar-name age,sex",
              " --glm",
              " --out ", result_path, "self_africa_result_bp"
))





#####################################################################
#####################################################################
##### Step 9: GWAS results summary
#####################################################################
#####################################################################
result_bp <- fread(file = paste0(result_path,"self_africa_result_bp.bp.glm.linear"))
result_hyper <- fread(file = paste0(result_path,"self_africa_result_hyper.hypertension.glm.logistic"))
names(result_bp)[1:3] <- c("CHR", "BP", "SNP")
names(result_hyper)[1:3] <- c("CHR", "BP", "SNP")
result_bp_G <- result_bp %>% filter(TEST == "ADD")
result_hyper_G <- result_hyper %>% filter(TEST == "ADD")
png(paste0(figure_path,"gwas_hyper.png")) 
manhattan(x = result_hyper_G, annotatePval = -log10(5e-08), col = 1:22, main = "Hypertension", suggestiveline = F)
dev.off()
png(paste0(figure_path,"gwas_bp.png")) 
manhattan(x = result_bp_G, annotatePval = -log10(5e-08), col = 1:22, main = "Diastolic Blood Pressure", suggestiveline = F)
dev.off()






#####################################################################
#####################################################################
##### Step 10: Extract the Top 5 SNPs in each trait
#####################################################################
#####################################################################
top_result_hyper <- result_hyper_G %>% arrange(by = P) %>% head(n = 5)
top_result_bp <- result_bp_G %>% arrange(by = P) %>% head(n = 5)
top_result_hyper_summary <- top_result_hyper[,c(1,2,3)]
top_result_hyper_summary$BetaG <- log(top_result_hyper$OR)
top_result_bp_summary <- top_result_bp[,c(1,2,3)]
top_result_bp_summary$BetaG <- top_result_bp$BETA
result_hyper_age <- result_hyper %>% filter(TEST == "age", SNP %in% top_result_hyper_summary$SNP) %>% mutate(BetaG = log(OR)) %>% select("SNP", "BetaG")
result_bp_age <- result_bp %>% filter(TEST == "age", SNP %in% top_result_bp_summary$SNP) %>% select("SNP", "BETA")
result_hyper_sex <- result_hyper %>% filter(TEST == "sex", SNP %in% top_result_hyper_summary$SNP) %>% mutate(BetaG = log(OR)) %>% select("SNP", "BetaG")
result_bp_sex <- result_bp %>% filter(TEST == "sex", SNP %in% top_result_bp_summary$SNP) %>% select("SNP", "BETA")
names(result_bp_sex)[2] <- "Beta_Sex"
names(result_bp_age)[2] <- "Beta_Age"
names(result_hyper_sex)[2] <- "Beta_Sex"
names(result_hyper_age)[2] <- "Beta_Age"
top_result_bp_summary <- merge(merge(top_result_bp_summary, result_bp_sex, by = "SNP"), result_bp_age, by = "SNP")
top_result_hyper_summary <- merge(merge(top_result_hyper_summary, result_hyper_sex, by = "SNP"), result_hyper_age, by = "SNP")
save(top_result_bp_summary, file = paste0(data_path, "top_result_bp_summary.rda"))
save(top_result_hyper_summary, file = paste0(data_path, "top_result_hyper_summary.rda"))




#####################################################################
#####################################################################
##### Step 11: Extract the results of Gamma_age, Gamma_sex
#####################################################################
#####################################################################
system(paste0(plink,
              " --bfile ", geno_path,"self_africa_filtered_data_afterSNP_QC",
              " --snps ", "rs472771, rs514400, rs596875, rs7158982, rs9313506, rs163913, rs1648707, rs17069257, rs765528, rs9540255",
              " --allow-no-sex ",
              " --pheno ", data_path ,"pheno.txt",
              " --pheno-name age",
              " --glm",
              " --out ", result_path, "self_africa_gamma_age_result"
))
system(paste0(plink,
              " --bfile ", geno_path,"self_africa_filtered_data_afterSNP_QC",
              " --snps ", "rs472771, rs514400, rs596875, rs7158982, rs9313506, rs163913, rs1648707, rs17069257, rs765528, rs9540255",
              " --allow-no-sex ",
              " --pheno ", data_path ,"pheno.txt",
              " --pheno-name sex",
              " --glm",
              " --out ", result_path, "self_africa_gamma_sex_result"
))
result_gamma_age <- fread(file = paste0(result_path,"self_africa_gamma_age_result.age.glm.linear"))
result_gamma_age$Gamma_Age <- result_gamma_age$BETA
result_gamma_sex <- fread(file = paste0(result_path,"self_africa_gamma_sex_result.sex.glm.logistic"))
result_gamma_sex$Gamma_Sex <- log(result_gamma_sex$OR)
save(result_gamma_sex, file = paste0(data_path, "result_gamma_sex.rda"))
save(result_gamma_age, file = paste0(data_path, "result_gamma_age.rda"))






