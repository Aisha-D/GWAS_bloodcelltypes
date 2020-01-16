## Read Meta CD4T output
## Plot Minimanhattan plot - use Bex's funtion
## Run LD on the significant on the same chr3
library(tidyr)
library(qqman)

setwd("/mnt/data1/GWAS_bloodcelltypes/MetaAnalysis/")
cd4t = read.table("/mnt/data1/GWAS_bloodcelltypes/MetaAnalysis/METAANALYSIS_CD4T_1.txt",
                stringsAsFactors = F, header = T)
cd4t_n <- cd4t
cd4t_n$ChrBp <- cd4t$MarkerName
cd4t_n <- separate(cd4t_n, ChrBp, into= c("Chr", "Bp"), sep = ":")
cd4t_n$Chr <- as.numeric(as.character(cd4t_n$Chr))
cd4t_n$Bp <- as.numeric(as.character(cd4t_n$Bp))
cd4t_n <- cd4t_n[order(cd4t_n$P.value),]

cd4t_n2 = cd4t_n

colnames(cd4t_n2)[colnames(cd4t_n2) == 'Chr'] <- 'CHR'
colnames(cd4t_n2)[colnames(cd4t_n2) == 'Bp'] <- 'BP'
colnames(cd4t_n2)[colnames(cd4t_n2) == 'MarkerName'] <- 'SNP'
colnames(cd4t_n2)[colnames(cd4t_n2) == 'P.value'] <- 'P'
rownames(cd4t_n2) = cd4t_n2$MarkerName

snpsOfInterest = cd4t_n2[1:10,1]
tensnps = cd4t_n2[1:10,]
chr3snps = tensnps[which(tensnps$CHR == 3),]
ldwindow = chr3snps[1,9] - chr3snps[5,9]#19397

#Calculate upper and lower points of snp of interests
snp_loc = cd4t_n2[1, 'BP']
padding = 150000
snp_low = snp_loc + padding
snp_high = snp_loc - padding

png("MetaAnalysis_cd4t_Chr3.png")
manhattan(subset(cd4t_n2, CHR == 3), highlight = snpsOfInterest, 
          xlim = c(snp_low, snp_high), main = "CD4T Chr3")

dev.off()

#Now check what the ld r score if for the other snps around the most significant snp
#Check this in extend first
#Run this part in Plink2. 

library(getmstatistic) 
getmstatistic_results <- getmstatistic(heartgenes214$beta_flipped, 
                                       heartgenes214$gcse, 
                                       heartgenes214$variants, 
                                       heartgenes214$studies)
getmstatistic_results


## PCA of blood cell types???
load("/mnt/data1/EXTEND/Methylation/QC/EXTEND_batch1_2_merged/EUR_unlrelated_QCd/EXTEND_batch1_2_genoQCd_Normalised.rdat")
load('/mnt/data1/reference_files/BloodCellPropCalc/Bloodcoefs_withNeu.rdat') #Probes used to calculate the cell types

celltypes_betas = betas[rownames(betas) %in% rownames(coefs),]
cellprcomp <- prcomp(celltypes_betas, center = TRUE, scale. = TRUE)                     
summary(cellprcomp)
library(ggfortify)
autoplot(cellprcomp, data = celltypes_betas) #colour by cell type?

