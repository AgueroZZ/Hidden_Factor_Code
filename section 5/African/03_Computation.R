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

geno_path = "/home/ziang/Documents/UKB/"
data_path = "/home/ziang/Documents/Stats-Gene/UKB/African_GWAS_hypertension_bp/data/"
plink = "/home/ziang/Documents/UKB/plink2"
plink19 = "/home/ziang/Documents/UKB/plink19"
result_path = "/home/ziang/Documents/Stats-Gene/UKB/African_GWAS_hypertension_bp/results/"
figure_path = "/home/ziang/Documents/Stats-Gene/UKB/African_GWAS_hypertension_bp/figures/"


# geno_path = "/Users/ziangzhang/Documents/GitHub/Stats-Gene/UKB"
# data_path = "/Users/ziangzhang/Documents/GitHub/Stats-Gene/UKB/African_GWAS_hypertension_bp/data/"
# result_path = "/Users/ziangzhang/Documents/GitHub/Stats-Gene/UKB/African_GWAS_hypertension_bp/results/"
# figure_path = "/Users/ziangzhang/Documents/GitHub/Stats-Gene/UKB/African_GWAS_hypertension_bp/figures/"



# geno_path = "D:/Stats-Gene/UKB/"
# data_path = "D:/Stats-Gene/UKB/African_GWAS_hypertension_bp/data/"
# result_path = "D:/Stats-Gene/UKB/African_GWAS_hypertension_bp/results/"
# figure_path = "D:/Stats-Gene/UKB/African_GWAS_hypertension_bp/figures/"

load(paste0(data_path, "final_data.rda"))
colnames(final_data)[1] <- "IID"
final_data <- cbind(FID = final_data$IID ,final_data)
final_data$hypertension <- ifelse(final_data$hypertension == 1, 2,1)
final_data$sex <- ifelse(final_data$sex == 1, 2,1)

load(paste0(data_path, "result_gamma_sex.rda"))
load(paste0(data_path, "result_gamma_age.rda"))
load(paste0(data_path, "top_result_bp_summary.rda"))
load(paste0(data_path, "top_result_hyper_summary.rda"))
MAF_bp <- fread(file = paste0(result_path, "self_africa_result_bp.afreq"))
MAF_bp <- MAF_bp %>% filter(ID %in% top_result_bp_summary$SNP) %>% mutate(SNP = ID, MAF = ALT_FREQS) %>% select(SNP, MAF)
MAF_hyper <- fread(file = paste0(result_path, "self_africa_result_hyper.afreq"))
MAF_hyper <- MAF_hyper %>% filter(ID %in% top_result_hyper_summary$SNP) %>% mutate(SNP = ID, MAF = ALT_FREQS) %>% select(SNP, MAF)
Gamma_sex_bp <- result_gamma_sex %>% filter(ID %in% top_result_bp_summary$SNP) %>% mutate(SNP = ID) %>% select(SNP, Gamma_Sex)
Gamma_sex_hyper <- result_gamma_sex %>% filter(ID %in% top_result_hyper_summary$SNP) %>% mutate(SNP = ID) %>% select(SNP, Gamma_Sex)
Gamma_age_bp <- result_gamma_age %>% filter(ID %in% top_result_bp_summary$SNP) %>% mutate(SNP = ID) %>% select(SNP, Gamma_Age)
Gamma_age_hyper <- result_gamma_age %>% filter(ID %in% top_result_hyper_summary$SNP) %>% mutate(SNP = ID) %>% select(SNP, Gamma_Age)

top_result_bp_summary <- merge(top_result_bp_summary, MAF_bp, by = "SNP")
top_result_hyper_summary <- merge(top_result_hyper_summary, MAF_hyper, by = "SNP")
top_result_bp_summary <- merge(top_result_bp_summary, Gamma_sex_bp, by = "SNP")
top_result_hyper_summary <- merge(top_result_hyper_summary, Gamma_sex_hyper, by = "SNP")
top_result_bp_summary <- merge(top_result_bp_summary, Gamma_age_bp, by = "SNP")
top_result_hyper_summary <- merge(top_result_hyper_summary, Gamma_age_hyper, by = "SNP")

#### Obtain the top SNP data:
## SNP names:
SNP_names <- top_result_bp_summary$SNP
## SNP chromosomes:
result <- c()
SNP_chro <- top_result_bp_summary$CHR
for (i in 1:length(SNP_chro)) {
  system(paste0(plink,
                " --bfile ", geno_path,"self_africa_filtered_data_afterSNP_QC",
                " --snp ", SNP_names[i],
                " --make-bed ",
                " --out ", result_path, "top_SNP_bp"))
  path <- paste(result_path, paste0("top_SNP_bp",".bed"), sep = '')
  tempfile <- tempfile()
  snp_readBed(path, backingfile = tempfile)
  SNP_data_full <- snp_attach(paste0(tempfile,".rds",sep=''))
  SNP_data <- data.frame(ID = SNP_data_full$fam$family.ID, sex = SNP_data_full$fam$sex)
  G <- SNP_data_full$genotypes
  result <- cbind(result, as.matrix(G[,]))
}
result <- cbind(SNP_data$ID, result)
SNP_data <- as.data.frame(result)
names(SNP_data) <- c("IID",top_result_bp_summary$SNP)
save(file = paste0(result_path, "top_SNP_data_bp.rda"), SNP_data)
ResidualSD <- c()
for (i in 1:length(top_result_bp_summary$SNP)) {
  full_data <- na.omit(left_join(x = SNP_data[,c(1,(i+1))], final_data[,c("IID", "age", "sex", "bp")], by = "IID"))
  names(full_data)[2] <- "SNP"
  reg1 <- glm(bp ~ SNP + age + sex, family = gaussian(), data = full_data)
  ResidualSD[i] <- sd(reg1$residual)
}







#####################################################################
#####################################################################
##### Step 12: Power Computation for the finished GWAS
#####################################################################
#####################################################################
Power_bp_with_E <- c()
Power_bp_without_E <- c()
for (i in 1:5) {
  parameters <- list(pG = top_result_bp_summary$MAF[i], pE = mean(final_data$sex-1), 
                     muE = mean(final_data$age), sigmaE = sd(final_data$age), 
                     betaG = top_result_bp_summary$BetaG[i], TraitMean = mean(final_data$bp), ResidualSD = ResidualSD[i],
                     betaE = c(top_result_bp_summary$Beta_Sex[i], top_result_bp_summary$Beta_Age[i]),
                     gammaG = c(top_result_bp_summary$Gamma_Sex[i], top_result_bp_summary$Gamma_Sex[i]))
  Power_bp_with_E[i] <- Compute_Power_multi(parameters = parameters, n = nrow(final_data), response = "continuous",
                                     covariate = c("binary", "continuous"), alpha = 5e-8, B = 100000)
  Power_bp_without_E[i] <- Compute_Power(parameters = parameters, n = nrow(final_data), response = "continuous",
                                            covariate = "none", alpha = 5e-8, B = 100000)
}

bp_result <- data.frame(SNP = top_result_bp_summary$SNP, Power_withE = Power_bp_with_E, Power_without_E = Power_bp_without_E)
bp_result %>% pivot_longer(cols = c(Power_without_E, Power_withE), values_to = "Power", names_to = "condition") %>%
  ggplot() + 
  geom_bar(aes(x = reorder(SNP, -Power), y = Power,fill = condition), position = "dodge", stat="identity") + 
  ylab("Computed Power") +
  xlab("SNP") +
  scale_fill_manual(labels = c("with E", "without E"), values = c("2", "4")) +
  ggtitle("Blood Pressure") +
  coord_cartesian(ylim=c(0.1,0.6)) +
  theme_classic()

ggsave(filename = paste0(figure_path, "bp_barplot.pdf"), device = "pdf", width = 5,
       height = 5)


Power_hyper_with_E <- c()
Power_hyper_without_E <- c()
for (i in 1:5) {
  parameters <- list(pG = top_result_hyper_summary$MAF[i], pE = mean(final_data$sex-1), 
                     muE = mean(final_data$age), sigmaE = sd(final_data$age), 
                     betaG = top_result_hyper_summary$BetaG[i], preva = mean(final_data$hypertension - 1),
                     betaE = c(top_result_hyper_summary$Beta_Sex[i], top_result_hyper_summary$Beta_Age[i]),
                     gammaG = c(top_result_hyper_summary$Gamma_Sex[i], top_result_hyper_summary$Gamma_Sex[i]))
  Power_hyper_with_E[i] <- Compute_Power_multi(parameters = parameters, n = nrow(final_data), response = "binary",
                                            covariate = c("binary", "continuous"), alpha = 5e-8, B = 100000)
  Power_hyper_without_E[i] <- Compute_Power(parameters = parameters, n = nrow(final_data), response = "binary",
                                         covariate = "none", alpha = 5e-8, B = 100000)
}

hypertension_result <- data.frame(SNP = top_result_hyper_summary$SNP, Power_withE = Power_hyper_with_E, Power_without_E = Power_hyper_without_E)

hypertension_result %>% pivot_longer(cols = c(Power_without_E, Power_withE), values_to = "Power", names_to = "condition") %>%
  ggplot() + 
  geom_bar(aes(x = reorder(SNP, -Power), y = Power,fill = condition), position = "dodge", stat="identity") + 
  ylab("Computed Power") +
  xlab("SNP") +
  scale_fill_manual(labels = c("with E", "without E"), values = c("2", "4")) +
  ggtitle("Hypertension") +
  coord_cartesian(ylim=c(0.1,0.6)) +
  theme_classic()

ggsave(filename = paste0(figure_path, "hypertension_barplot.pdf"), device = "pdf", width = 5,
       height = 5)







#####################################################################
#####################################################################
##### Step 13: Power Computation for the replication GWAS
#####################################################################
#####################################################################
Power_bp_with_E <- c()
Power_bp_without_E <- c()
for (i in 1:5) {
  parameters <- list(pG = top_result_bp_summary$MAF[i], pE = mean(final_data$sex-1), 
                     muE = mean(final_data$age), sigmaE = sd(final_data$age), 
                     betaG = top_result_bp_summary$BetaG[i], TraitMean = mean(final_data$bp), ResidualSD = ResidualSD[i],
                     betaE = c(top_result_bp_summary$Beta_Sex[i], top_result_bp_summary$Beta_Age[i]),
                     gammaG = c(top_result_bp_summary$Gamma_Sex[i], top_result_bp_summary$Gamma_Sex[i]))
  Power_bp_with_E[i] <- Compute_Power_multi(parameters = parameters, n = nrow(final_data), response = "continuous",
                                            covariate = c("binary", "continuous"), alpha = 5e-3, B = 100000)
  Power_bp_without_E[i] <- Compute_Power(parameters = parameters, n = nrow(final_data), response = "continuous",
                                         covariate = "none", alpha = 5e-3, B = 100000)
}

bp_result <- data.frame(SNP = top_result_bp_summary$SNP, Power_withE = Power_bp_with_E, Power_without_E = Power_bp_without_E)
bp_result %>% pivot_longer(cols = c(Power_without_E, Power_withE), values_to = "Power", names_to = "condition") %>%
  ggplot() + 
  geom_bar(aes(x = reorder(SNP, -Power), y = Power,fill = condition), position = "dodge", stat="identity") + 
  ylab("Computed Power") +
  xlab("SNP") +
  scale_fill_manual(labels = c("with E", "without E"), values = c("2", "4")) +
  ggtitle("Blood Pressure") +
  coord_cartesian(ylim=c(0.1,1)) +
  theme_classic()

ggsave(filename = paste0(figure_path, "bp_barplot_replicate.pdf"), device = "pdf", width = 5,
       height = 5)


Power_hyper_with_E <- c()
Power_hyper_without_E <- c()
for (i in 1:5) {
  parameters <- list(pG = top_result_hyper_summary$MAF[i], pE = mean(final_data$sex-1), 
                     muE = mean(final_data$age), sigmaE = sd(final_data$age), 
                     betaG = top_result_hyper_summary$BetaG[i], preva = mean(final_data$hypertension - 1),
                     betaE = c(top_result_hyper_summary$Beta_Sex[i], top_result_hyper_summary$Beta_Age[i]),
                     gammaG = c(top_result_hyper_summary$Gamma_Sex[i], top_result_hyper_summary$Gamma_Sex[i]))
  Power_hyper_with_E[i] <- Compute_Power_multi(parameters = parameters, n = nrow(final_data), response = "binary",
                                               covariate = c("binary", "continuous"), alpha = 5e-3, B = 100000)
  Power_hyper_without_E[i] <- Compute_Power(parameters = parameters, n = nrow(final_data), response = "binary",
                                            covariate = "none", alpha = 5e-3, B = 100000)
}

hypertension_result <- data.frame(SNP = top_result_hyper_summary$SNP, Power_withE = Power_hyper_with_E, Power_without_E = Power_hyper_without_E)

hypertension_result %>% pivot_longer(cols = c(Power_without_E, Power_withE), values_to = "Power", names_to = "condition") %>%
  ggplot() + 
  geom_bar(aes(x = reorder(SNP, -Power), y = Power,fill = condition), position = "dodge", stat="identity") + 
  ylab("Computed Power") +
  xlab("SNP") +
  scale_fill_manual(labels = c("with E", "without E"), values = c("2", "4")) +
  ggtitle("Hypertension") +
  coord_cartesian(ylim=c(0.1,1)) +
  theme_classic()

ggsave(filename = paste0(figure_path, "hypertension_barplot_replicate.pdf"), device = "pdf", width = 5,
       height = 5)




#####################################################################
#####################################################################
##### Step 14: Sample Size Computation for the replication GWAS
#####################################################################
#####################################################################
Size_bp_with_E <- c()
Size_bp_without_E <- c()
for (i in 1:5) {
  parameters <- list(pG = top_result_bp_summary$MAF[i], pE = mean(final_data$sex-1), 
                     muE = mean(final_data$age), sigmaE = sd(final_data$age), 
                     betaG = top_result_bp_summary$BetaG[i], TraitMean = mean(final_data$bp), ResidualSD = ResidualSD[i],
                     betaE = c(top_result_bp_summary$Beta_Sex[i], top_result_bp_summary$Beta_Age[i]),
                     gammaG = c(top_result_bp_summary$Gamma_Sex[i], top_result_bp_summary$Gamma_Sex[i]))
  Size_bp_with_E[i] <- Compute_Size_multi(parameters = parameters, PowerAim = 0.8, response = "continuous",
                                            covariate = c("binary", "continuous"), alpha = 5e-3, B = 100000)
  Size_bp_without_E[i] <- round(Compute_Size(parameters = parameters, PowerAim = 0.8, response = "continuous",
                                         covariate = "none", alpha = 5e-3, B = 100000),0)
}

bp_result <- data.frame(SNP = top_result_bp_summary$SNP, Size_withE = Size_bp_with_E, Size_without_E = Size_bp_without_E)
bp_result %>% pivot_longer(cols = c(Size_without_E, Size_withE), values_to = "Sample_Size", names_to = "condition") %>%
  ggplot() + 
  geom_bar(aes(x = reorder(SNP, -Sample_Size), y = Sample_Size,fill = condition), position = "dodge", stat="identity") + 
  ylab("Replication Sample Size") +
  xlab("SNP") +
  scale_fill_manual(labels = c("with E", "without E"), values = c("2", "4")) +
  ggtitle("Blood Pressure") +
  coord_cartesian(ylim=c(1000,2000)) +
  theme_classic()

ggsave(filename = paste0(figure_path, "bp_barplot_size.pdf"), device = "pdf", width = 5,
       height = 5)


Size_hyper_with_E <- c()
Size_hyper_without_E <- c()
for (i in 1:5) {
  parameters <- list(pG = top_result_hyper_summary$MAF[i], pE = mean(final_data$sex-1), 
                     muE = mean(final_data$age), sigmaE = sd(final_data$age), 
                     betaG = top_result_hyper_summary$BetaG[i], preva = mean(final_data$hypertension - 1),
                     betaE = c(top_result_hyper_summary$Beta_Sex[i], top_result_hyper_summary$Beta_Age[i]),
                     gammaG = c(top_result_hyper_summary$Gamma_Sex[i], top_result_hyper_summary$Gamma_Sex[i]))
  Size_hyper_with_E[i] <- Compute_Size_multi(parameters = parameters, PowerAim = 0.8, response = "binary",
                                               covariate = c("binary", "continuous"), alpha = 5e-3, B = 100000)
  Size_hyper_without_E[i] <- round(Compute_Size(parameters = parameters, PowerAim = 0.8, response = "binary",
                                            covariate = "none", alpha = 5e-3, B = 100000), 0)
}

hypertension_result <- data.frame(SNP = top_result_hyper_summary$SNP, Size_withE = Size_hyper_with_E, Size_without_E = Size_hyper_without_E)

hypertension_result %>% pivot_longer(cols = c(Size_without_E, Size_withE), values_to = "size", names_to = "condition") %>%
  ggplot() + 
  geom_bar(aes(x = reorder(SNP, -size), y = size,fill = condition), position = "dodge", stat="identity") + 
  ylab("Replication Sample Size") +
  xlab("SNP") +
  scale_fill_manual(labels = c("with E", "without E"), values = c("2", "4")) +
  ggtitle("Hypertension") +
  coord_cartesian(ylim=c(1000, 2000)) +
  theme_classic()

ggsave(filename = paste0(figure_path, "hypertension_barplot_size.pdf"), device = "pdf", width = 5,
       height = 5)







