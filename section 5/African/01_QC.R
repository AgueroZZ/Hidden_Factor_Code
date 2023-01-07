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
save(pheno_data, file = "/home/ziang/Documents/Stats-Gene/UKB/African_GWAS_hypertension_bp/data/pheno_data.rda")


#################################################
#################################################
##### Step 1: Subset out self-reported African
#################################################
#################################################
names_want_pop <- grep("f.21000.0",bd_names,value = T)
pop_want <- fread("/home/ziang/Documents/UKB/ukb47570.tab", header=TRUE, sep="\t", select = c("f.eid",names_want_pop))
pop_want_african <- pop_want %>% filter(f.21000.0.0 == 4002)
african_ID <- pop_want_african$f.eid
pheno_data <- pheno_data %>% filter(ID %in% african_ID)
### Number of self-reported African individuals to start with:
nrow(pheno_data)




#####################################################################
#####################################################################
##### Step 2: Merge the plink file for all self_report_african
#####################################################################
#####################################################################
geno_path = "/home/ziang/Documents/UKB/"
plink = "/home/ziang/Documents/UKB/plink2"
plink19 = "/home/ziang/Documents/UKB/plink19"
matrix_to_write <- cbind(african_ID,african_ID)
write.table(x = matrix_to_write, file = "/home/ziang/Documents/UKB/self_report_african.txt", row.names = F, col.names = F)
result <- NULL
setwd("/home/ziang/Documents/UKB")
for (i in 1:22) {
  i_chromo <- paste0("c",i)
  i_selected_name <- paste0(paste0("selected_",i_chromo),"_african")
  system(paste0(plink,
                paste(" --bfile", i_chromo),
                " --keep self_report_african.txt",
                " --make-bed",
                paste(" --out", i_selected_name)))
}
system(paste0(plink19,
              " --bfile ", "selected_c1_african",
              " --merge-list ", "merge_list.txt",
              " --make-bed",
              " --out ", "self_africa_merged_data"))



#####################################################################
#####################################################################
##### Step 3: Filter Indivds by KIN and MIND
#####################################################################
#####################################################################
system(paste0(plink,
              paste(" --bfile", "self_africa_merged_data"),
              " --king-cutoff 0.25",
              " --mind 0.2",
              " --make-bed ",
              paste(" --out", "self_africa_filtered_data")
))
# 0 samples removed due to missing genotype data (--mind).
# self_africa_filtered_data.king.cutoff.out.id , and 3182 remaining sample IDs





#####################################################################
#####################################################################
##### Step 4: Filter Indivds by overall PC with Kmeans
#####################################################################
#####################################################################
names_want_PC <- grep("f.22009", bd_names,value = T)
names_want_PCi <- names_want_PC[c(1:4,15)]
PCi_want <- fread("/home/ziang/Documents/UKB/ukb47570.tab", header=TRUE, sep="\t", select = c("f.eid",names_want_PCi))
names(PCi_want) <- c("ID", "PC1", "PC2", "PC3", "PC4", "PC15")
PC_data_full <- left_join(x = PCi_want, y = pheno_data, by = "ID")
PC_data_full$Ancestry <- ifelse(is.na(PC_data_full$asthma), "Other", "African (Self-Reported)")
PC_data_full %>% ggplot() + geom_point(aes(x = PC1, y = PC2, color = Ancestry, shape = Ancestry), alpha = 0.3, size = 2) + theme_classic() +
  xlim(0,500) + theme(legend.position = c(0.8,0.3))
ggsave(filename = "/home/ziang/Documents/Stats-Gene/UKB/African_GWAS_hypertension_bp/figures/All_PC_Figure12.pdf", device = "pdf")
PC_data_full %>% ggplot() + geom_point(aes(x = PC3, y = PC4, color = Ancestry, shape = Ancestry), alpha = 0.3, size = 2) + theme_classic()
ggsave(filename = "/home/ziang/Documents/Stats-Gene/UKB/African_GWAS_hypertension_bp/figures/All_PC_Figure34.pdf", device = "pdf")
pheno_data_merged <- merge(pheno_data, PCi_want, by = "ID")
select_ID <- read.table("self_africa_filtered_data.king.cutoff.in.id")
select_ID <- select_ID$V1
### Only use unrelated individuals for clustering
pheno_data_selected <- pheno_data_merged %>% filter(ID %in% select_ID)
Restrict_Data <- pheno_data_selected %>% na.omit()
Restrict_Data12 <- Restrict_Data[, c("PC1", "PC2")]
set.seed(123)
kmeans_mod12 <- kmeans(Restrict_Data12, centers = 4, nstart = 100)
Restrict_Data$Clusters <- as.character(kmeans_mod12$cluster)
Restrict_Data %>% ggplot() + geom_point(aes(x = PC1, y = PC2, color = Clusters), alpha = 0.3, size = 2) + theme_classic()
ggsave(filename = "/home/ziang/Documents/Stats-Gene/UKB/African_GWAS_hypertension_bp/figures/African_PC_cluster_Figure12.pdf", device = "pdf")
table(Restrict_Data$Clusters)/sum(table(Restrict_Data$Clusters))
largest_Cluster_ID <- which.is.max(table(Restrict_Data$Clusters))
Selected_ID <- Restrict_Data$ID[Restrict_Data$Clusters == largest_Cluster_ID]
Cluster_Restrict_Data <- Restrict_Data[Restrict_Data$Clusters == largest_Cluster_ID,]
nrow(Cluster_Restrict_Data) #2601
save(Cluster_Restrict_Data, file = "/home/ziang/Documents/Stats-Gene/UKB/African_GWAS_hypertension_bp/data/Cluster_Restrict_Data.rda")

Restrict_Data %>% ggplot() + geom_point(aes(x = PC1, y = PC2, color = Clusters), alpha = 0.3, size = 2) + theme_classic() + xlim(c(0,500)) + theme(legend.position = c(0.9,0.3))
ggsave(filename = "/home/ziang/Documents/Stats-Gene/UKB/African_GWAS_hypertension_bp/figures/African_PC_Figure12.pdf", device = "pdf")
Restrict_Data %>% ggplot() + geom_point(aes(x = PC1, y = PC3, color = Clusters), alpha = 0.3, size = 2) + theme_classic() + xlim(c(0,500)) + theme(legend.position = c(0.9,0.3))
ggsave(filename = "/home/ziang/Documents/Stats-Gene/UKB/African_GWAS_hypertension_bp/figures/African_PC_Figure13.pdf", device = "pdf")
Restrict_Data %>% ggplot() + geom_point(aes(x = PC2, y = PC3, color = Clusters), alpha = 0.3, size = 2) + theme_classic() + xlim(c(-150,100)) + theme(legend.position = c(0.9,0.8))
ggsave(filename = "/home/ziang/Documents/Stats-Gene/UKB/African_GWAS_hypertension_bp/figures/African_PC_Figure23.pdf", device = "pdf")
Restrict_Data %>% ggplot() + geom_point(aes(x = PC3, y = PC4, color = Clusters), alpha = 0.3, size = 2) 
ggsave(filename = "/home/ziang/Documents/Stats-Gene/UKB/African_GWAS_hypertension_bp/figures/African_PC_Figure34.pdf", device = "pdf")
Restrict_Data %>% ggplot() + geom_point(aes(x = PC1, y = PC15, color = Clusters), alpha = 0.3, size = 2) 
ggsave(filename = "/home/ziang/Documents/Stats-Gene/UKB/African_GWAS_hypertension_bp/figures/African_PC_Figure115.pdf", device = "pdf")






#####################################################################
#####################################################################
##### Step 5: Compute new PCs and remove outliers
#####################################################################
#####################################################################
system(paste0(plink19,
              " --bfile ", "self_africa_filtered_data",
              " --pca 20 header ",
              " --make-bed",
              " --out ", "self_africa_PCA"))
newPCs <- read.table(file = "/home/ziang/Documents/Stats-Gene/UKB/African_GWAS_hypertension_bp/data/self_africa_PCA.eigenvec", header = T)
newPCs$ID <- newPCs$IID
Cluster_Data <- data.frame(ID = (pheno_data_selected %>% na.omit())$ID, Cluster = as.character(kmeans_mod12$cluster))
newPCs_clustered <- inner_join(newPCs, Cluster_Data, by = "ID")
newPCS_cluster1 <- newPCs %>% filter(ID %in% Selected_ID)
newPCS_cluster1 <- newPCS_cluster1[,-c(1:2)]
newPCs_clustered %>% ggplot() + geom_point(aes(x = PC1, y = PC2, color = Cluster), alpha = 0.3, size = 2) + theme_classic() + theme(legend.position = c(0.1,0.3))
ggsave(filename = "/home/ziang/Documents/Stats-Gene/UKB/African_GWAS_hypertension_bp/figures/African_newPC_Figure12.pdf", device = "pdf")
newPCs_clustered %>% ggplot() + geom_point(aes(x = PC1, y = PC3, color = Cluster), alpha = 0.3, size = 2) + theme_classic() + theme(legend.position = c(0.1,0.3))
ggsave(filename = "/home/ziang/Documents/Stats-Gene/UKB/African_GWAS_hypertension_bp/figures/African_newPC_Figure13.pdf", device = "pdf")
newPCs_clustered %>% ggplot() + geom_point(aes(x = PC2, y = PC3, color = Cluster), alpha = 0.3, size = 2) + theme_classic() + theme(legend.position = c(0.8,0.3))
ggsave(filename = "/home/ziang/Documents/Stats-Gene/UKB/African_GWAS_hypertension_bp/figures/African_newPC_Figure23.pdf", device = "pdf")
newPCs_clustered %>% ggplot() + geom_point(aes(x = PC3, y = PC4, color = Cluster), alpha = 0.3, size = 2) + theme_classic() + theme(legend.position = c(0.1,0.3))
ggsave(filename = "/home/ziang/Documents/Stats-Gene/UKB/African_GWAS_hypertension_bp/figures/African_newPC_Figure34.pdf", device = "pdf")
Identify_outliers <- function(all_data, index){
  PCi <- all_data[,index]
  ID <- all_data$ID
  PC_mu <- mean(PCi)
  PC_sd <- sd(PCi)
  removed_ID <- ID[PCi>=(PC_mu + 4*PC_sd) | PCi<=(PC_mu - 4*PC_sd)]
  removed_ID
}
ID_to_removes <- c()
for (i in 1:20) {
  ID_to_removes <- unique(c(ID_to_removes,Identify_outliers(newPCS_cluster1, index = i)))
}
length(ID_to_removes) #91
final_ID <- newPCS_cluster1$ID[!newPCS_cluster1$ID %in% ID_to_removes]
length(final_ID) #2510
final_data <- Cluster_Restrict_Data[Cluster_Restrict_Data$ID %in% final_ID,]
save(final_data, file = "/home/ziang/Documents/Stats-Gene/UKB/African_GWAS_hypertension_bp/data/final_data.rda")

### Elbow Plot using new PCs:
newPC_eigenvec <- read.table(file = paste0(data_path, "self_africa_PCA.eigenvec"), header = T)
newPC_eigen <- read.table(file = paste0(data_path, "self_africa_PCA.eigenval"), header = F)
newPC_eigen <- as.data.frame(newPC_eigen)
newPC_eigen$PC <- 1:nrow(newPC_eigen)
names(newPC_eigen) <- c("Variance", "PC")
newPC_eigen$variance_explained <- newPC_eigen$Variance/sum(newPC_eigen$Variance)
newPC_eigen$cumu_var <- cumsum(newPC_eigen$Variance)
newPC_eigen %>% ggplot(aes(x = PC, y = variance_explained)) + geom_point() + 
  geom_line() + 
  ylab("Variance Percentage") +
  ggtitle("Percentage of variance explained by each PC") +
  theme_classic()
ggsave(filename = paste0(figure_path,"elbow_plot.pdf"), device = "pdf", width = 5,
       height = 5)
ggsave(filename = paste0(figure_path,"elbow_plot.png"), device = "png", width = 5,
       height = 5)






#####################################################################
#####################################################################
##### Step 6: Filter SNPs and record the number
#####################################################################
#####################################################################
matrix_to_write <- cbind(final_ID,final_ID)
write.table(x = matrix_to_write, file = "/home/ziang/Documents/UKB/final_ID.txt", row.names = F, col.names = F)
system(paste0(plink,
              paste(" --bfile", "self_africa_filtered_data"),
              " --keep final_ID.txt",
              " --hwe 1e-10",
              " --maf 0.01",
              " --geno 0.2",
              " --make-bed ",
              paste(" --out", "self_africa_filtered_data_afterSNP_QC")
))
# 784256 variants loaded from self_africa_filtered_data.bim.
# 384645 variants removed due to allele frequency threshold(s)
# --hwe: 2603 variants removed due to Hardy-Weinberg exact test (founders only).
# --geno: 18005 variants removed due to missing genotype data.
# 379003 variants remaining after main filters.











