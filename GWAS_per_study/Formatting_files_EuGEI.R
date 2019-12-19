## This script takes in the genotype and methylation data to estimate cell types (phenotype files)
## And the covariates (age,sex,pcs) for PLINK GWAS
setwd('/mnt/data1/GWAS_bloodcelltypes/Scripts/')
library(methylumi)
library(wateRmelon)
require(gdata)
library(minfi)
library(ggplot2)
require(gridExtra)
library(plyr)
require(IlluminaHumanMethylationEPICmanifest)
library(dplyr)
library(tidyr)
library(snpStats)
library(factoextra)
load("/mnt/data1/EuGEI/QC/GeorginasQC/All_Plates_Blood_WithRepeats/JustEuGEIresults/EuGEIBloodSamples_Normalised.rdat") #934 samples

#############       Phenotype txt file for PLINK          ##################
dim(pheno) #934 samples
#Remove cases
pheno <- pheno[,-which(names(pheno) %in% c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"))]
# pheno <- pheno[which(pheno$Phenotype.EuGEI == "Control"),] #521 samples
# betas <- betas[,which(pheno$Basename %in% colnames(betas))]
dim(pheno) #521 samples

#Add Cell counts
source("/mnt/data1/reference_files/BloodCellPropCalc/cellpropfunctions.r") #loads functions
load('/mnt/data1/reference_files/BloodCellPropCalc/Bloodcoefs_withNeu.rdat') #loads the bloodcoefs to calculate cell proportion
library(genefilter)
library(quadprog)
library(matrixStats)

coefs2 <- coefs[rownames(coefs) %in% rownames(betas),]
betas2 <- betas[rownames(coefs2),]
counts <- projectCellType(Y = betas2, coefCellType = coefs2)
rownames(counts) <- colnames(betas2)
save(counts, file = '/mnt/data1/GWAS_bloodcelltypes/EuGEI/cellcounts.rdat')
boxplot(counts, las = 2)


#############       Covariates txt file for PLINK         ###################
#Imputed
pathM <- paste("/mnt/data1/EuGEI/Genotypes/Imputed/EUGEI_GROUP_STrios_TrulyFinal_XVII_BestGuess", c(".bed", ".bim", ".fam"), sep = "")
##This contains European samples only?
SNP_M <- read.plink(pathM[1], pathM[2], pathM[3])
fam_imputed <- SNP_M$fam #Take the fam file to add the IID and FID to covaraites table and make sure the methylation and genotype files match
dim(fam_imputed) #899 samples genotype info

#Not imputed
pathM <- paste("/mnt/data1/EuGEI/Genotypes/EUGEI_Methylation_Samples_PostQC_OUT", c(".bed", ".bim", ".fam"), sep = "")
##This contains European samples only?
SNP_M <- read.plink(pathM[1], pathM[2], pathM[3])
fam_unimputed <- SNP_M$fam #Take the fam file to add the IID and FID to covaraites table and make sure the methylation and genotype files match
dim(fam_unimputed) #899 samples genotype info

#check which samples overlap and keep these from bed files
length(intersect(fam_imputed$pedigree, fam_unimputed$pedigree)) #1203 overlap
keep <- intersect(fam_imputed$pedigree, fam_unimputed$pedigree)
keep <- as.data.frame(cbind(keep, keep))
colnames(keep) <- c("FID", "IID")
keep$FID <- as.character(keep$FID)
keep$IID <- as.character(keep$IID)
write.table(keep, '/mnt/data1/GWAS_bloodcelltypes/EuGEI/keep.txt',sep = "\t", quote = FALSE, row.names = F)


##PCA analysis here!
#read in the PCA eigenvectors and visualise if samples include other ethnicities
pca <- read.table("/mnt/data1/GWAS_bloodcelltypes/EuGEI/EUGEI_GROUP_STrios_TrulyFinal_XVII_BestGuess_filtered.pca.eigenvec", header = T)
pcs <-pca[,3:22]
pdf("/mnt/data1/GWAS_bloodcelltypes/EuGEI/pca.pdf")
plot(pcs[,1], pcs[,2], xlab = "PC1", ylab = "PC2")
dev.off()

#The pca showed that there were a samples that are african/asian
#remove those samples not in keep file from the pheno file
length(intersect(pheno$Geno.CHIP.Location, keep$FID)) #922 samples that have methy and geno
pheno_fam <- pheno[which(pheno$Geno.CHIP.Location %in% keep$FID),]

##There is a spreadsheet which includes ethnic info - do they overlap with the 1203 samples
ethnichart <- read.csv("/mnt/data1/EuGEI/Phenotypes/Marta_AdditonalPhenoInfo.csv", header = T)
#need to rename the Eilis sample name 
pheno_fam$Sample_Name <- sub('\\.', '', pheno_fam$Eilis.Sample_Name)
pheno_fam$Sample_Name <- gsub("\\..*","",pheno_fam$Sample_Name) #remove dots after it
pheno_fam$Sample_Name <- sub('\\-', '', pheno_fam$Sample_Name) #remove dashes

#cbind pheno with clinical info from ethnichart
pheno_ethnic <- cbind(pheno_fam, ethnichart[match(pheno_fam$Sample_Name, ethnichart$st_subjid), 
                                             c("mrc1_socde03", "mrc1_socde04_amnd")])
pheno_white <- pheno_ethnic[which(pheno_ethnic$mrc1_socde03 == "White"),] #730 samples are white samples
pheno_white$mrc1_socde04_amnd <- as.character(pheno_white$mrc1_socde04_amnd)
pheno_white$mrc1_socde04_amnd <- as.factor(pheno_white$mrc1_socde04_amnd)
levels(pheno_white$mrc1_socde04_amnd)

#Now put these ID as keep_european then filter do a PCA to see any outliers from white samples
keep_europe <- as.data.frame(cbind(pheno_white$Geno.CHIP.Location, pheno_white$Geno.CHIP.Location))
colnames(keep_europe) <- c("FID", "IID")
keep_europe$FID <- as.character(keep_europe$FID)
keep_europe$IID <- as.character(keep_europe$IID)
write.table(keep_europe, '/mnt/data1/GWAS_bloodcelltypes/EuGEI/keep_european.txt',sep = "\t", quote = FALSE, row.names = F)


#match the methylation and genotype information
conver <- read.csv("/mnt/data1/EuGEI/Genotypes/EUGEI_SNP_genotypes_matchedIDs.csv", row.names = 1)
pheno$Geno.Plate.ID <- rep(NA, nrow(pheno))
pheno$Geno.Plate.ID <- conver$Geno.Plate.ID[match(pheno$Eilis.Sample_Name, conver$Sample.Name)]
fam2 <- fam[which(pheno$Geno.Plate.ID %in% fam$pedigree),] #this drops the number of samples to 1111
keep <- as.data.frame(cbind(pheno$Geno.Plate.ID, pheno$Geno.Plate.ID))

write.table(keep, '/mnt/data1/GWAS_bloodcelltypes/EuGEI/keep.txt',sep = "\t", quote = FALSE, row.names = F)

## On plink update the samples and use this code: 
## plink --bfile data --keep keep.txt --make-bed --out data_keep

pheno2 <- pheno[which(pheno$Geno.CHIP.Location %in% fam2$pedigree),]
rownames(pheno) <- pheno$barcode
betas <- betas[,which(colnames(betas) %in% rownames(pheno))]
rm(fam2)

##Check PCA of european samples
pca <- read.table("/mnt/data1/GWAS_bloodcelltypes/EuGEI/EUGEI_filtered_white.pca.eigenvec", header = T)
pcs <-pca[,3:22]
pdf("/mnt/data1/GWAS_bloodcelltypes/EuGEI/pca_whiteeuropean.pdf")
plot(pcs[,1], pcs[,2], xlab = "PC1", ylab = "PC2")
dev.off()
pcs2 <- pcs[which(pcs$X.0.00476049 <= 0.1),] #there is some trailing samples, filter these
pdf("/mnt/data1/GWAS_bloodcelltypes/EuGEI/pca_whiteeuropeanfiltered.pdf")
plot(pcs2[,1], pcs2[,2], xlab = "PC1", ylab = "PC2")
dev.off()

pcs3 <- pca[which(rownames(pca) %in% rownames(pcs2)),]
#filter these names as true european whites
pcs3$X3999870009_R01C01 <- as.character(pcs3$X3999870009_R01C01)
pcs3$X3999870009_R01C01.1 <- as.character(pcs3$X3999870009_R01C01.1)
keep_white <- as.data.frame(cbind(pcs3$X3999870009_R01C01, pcs3$X3999870009_R01C01))
colnames(keep_white) <- c("FID", "IID")
keep_white$FID <- as.character(keep_white$FID)
keep_white$IID <- as.character(keep_white$IID)
write.table(keep_white, '/mnt/data1/GWAS_bloodcelltypes/EuGEI/keep_whiteeuropean.txt',sep = "\t", quote = FALSE, row.names = F)

#update pheno with white europeans
pheno_whiteeuro <- pheno_white[which(pcs3$X3999870009_R01C01 %in% pheno_white$Geno.CHIP.Location),]
pheno_whiteeuro$mrc1_socde04_amnd <- as.character(pheno_whiteeuro$mrc1_socde04_amnd)
pheno_whiteeuro$mrc1_socde04_amnd <- as.factor(pheno_whiteeuro$mrc1_socde04_amnd)
levels(pheno_whiteeuro$mrc1_socde04_amnd)

##Filter out Non european samples 
non.european <- c("Asia", "Argentia", "America", "Chile", "Colombia", "Ecuador", 
                  "Non-Western immigrant (Africa)",
                  "Non-Western immigrant (Latin&South America ", 
                  "Peru", "Western immigrant", "Other white", "Other")
pheno_whiteeuro$mrc1_socde04_amnd <- as.character(pheno_whiteeuro$mrc1_socde04_amnd)
pheno_whiteeuro2 <- pheno_whiteeuro[which(!pheno_whiteeuro$mrc1_socde04_amnd %in% 
                                            non.european),]
pheno_whiteeuro2$mrc1_socde04_amnd <- as.character(pheno_whiteeuro2$mrc1_socde04_amnd)
pheno_whiteeuro2$mrc1_socde04_amnd <- as.factor(pheno_whiteeuro2$mrc1_socde04_amnd)
levels(pheno_whiteeuro2$mrc1_socde04_amnd)
#Now we have 671 samples
#


#cbind the pheno,cellcounts,pca
#rename the sample IDs in the first two columns of pca file
#match the number of samples from pheno and cellcoutns
counts <- counts[which(rownames(counts) %in% pheno_whiteeuro2$Basename),]
cov <- cbind(pheno_whiteeuro2, counts[match(rownames(counts), pheno_whiteeuro2$Basename),])
cov <- cbind(cov, pcs3[match(cov$Geno.CHIP.Location, pcs3$X3999870009_R01C01), 3:12])
#remove columns which are not needed
remove_cols <- c("Basename","Eilis.Sample_Name", "Cohort","Country","MethPlate","Meth.CHIP.ID","Meth.CHIP.Location","Batch",             
"SampleInfo","Phenotype.EuGEI","Psychotic.22q11",
"Smoking","PredictedAge","SmokingScore", "Sample_Name","mrc1_socde03",
"mrc1_socde04_amnd")
cov <- cov[ , -which(names(cov) %in% remove_cols)]
#MAKE the IID and FID columns match the fam style with _
cov$IID <- cov$Geno.CHIP.Location
cov$FID <- cov$Geno.CHIP.Location
cov <- cov[,-1]
cov <- cov[,c("FID","IID","Age","Sex","CD8T","CD4T","NK" ,          
              "Bcell","Mono","Gran","Eos","Neu" ,  "X.0.00476049", "X0.0330666",  
              "X0.00232756","X.0.00612395","X.0.0088076","X0.00624999","X0.00308993",
              "X.0.00246849","X0.0135685","X.0.00700942"
              
)]
colnames(cov)<-c("FID","IID", "Age", "Sex", "CD8T","CD4T", "NK",
                 "Bcell","Mono","Gran",
                 "Eos", "Neu", paste("PC", 1:10, sep = ""))
cov2 <- cov[,c("FID","IID","Age","Sex", paste("PC", 1:10, sep = ""))]
#covariates file
cov2$Age <- as.numeric(as.character(cov2$Age))
cov2$Sex[cov2$Sex == 'Female'] <- 2
cov2$Sex[cov2$Sex == 'Male'] <- 1
write.table(cov2, '/mnt/data1/GWAS_bloodcelltypes/EuGEI/cov_agesexpcs.txt',sep = "\t", quote = FALSE, row.names = F)

#pheno file for cells
cellphen <- cov[,c("FID","IID", "CD8T","CD4T", "NK",
                   "Bcell","Mono","Gran",
                   "Eos", "Neu")]
write.table(cellphen, '/mnt/data1/GWAS_bloodcelltypes/EuGEI/cell_phenotype.txt',sep = "\t", quote = FALSE, row.names = F)

#update white samples to keep in bed files
keep_white <- as.data.frame(cbind(cov$FID, cov$IID))
colnames(keep_white) <- c("FID", "IID")
keep_white$FID <- as.character(keep_white$FID)
keep_white$IID <- as.character(keep_white$IID)
write.table(keep_white, '/mnt/data1/GWAS_bloodcelltypes/EuGEI/keep_whiteeuropean.txt',sep = "\t", quote = FALSE, row.names = F)

#sex update
sexupd <- as.data.frame(cbind(cov2$FID, cov2$IID, cov2$Sex))
colnames(sexupd) <- c("FID", "IID", "Sex")
sexupd$FID <- as.character(sexupd$FID)
sexupd$IID <- as.character(sexupd$IID)
sexupd$Sex <- as.character(sexupd$Sex)
write.table(sexupd, '/mnt/data1/GWAS_bloodcelltypes/EuGEI/EuGEI_sex_update.txt',sep = "\t", quote = FALSE, row.names = F)

#case vs control?
pheno_whiteeuro3 <- pheno_whiteeuro[which(pheno_whiteeuro$Phenotype.EuGEI == "Control"),]
