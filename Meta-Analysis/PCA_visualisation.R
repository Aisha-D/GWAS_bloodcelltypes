##### PCA for ethnicity overview
setwd('/mnt/data1/GWAS_bloodcelltypes/MetaAnalysis_IVW/Ethnicity_Check/')
#use the .eigenvec file for the pcs
pcs = read.table("AllCohort_Merged1KG.pca.eigenvec", stringsAsFactors = F)
library(RColorBrewer)
library(snpStats)

## Label our populations
# 1000Genomes
pop <- read.csv('/mnt/data1/reference_files/1000G/phase3_sample_info.csv', 
                stringsAsFactors = F, header = T)
# EXTEND
extend_fam <- read.table("/mnt/data1/GWAS_bloodcelltypes/MetaAnalysis_IVW/Ethnicity_Check/Extend_1kgIDs.fam", 
                         stringsAsFactors = FALSE)
colnames(extend_fam) <- c("pedigree","member","father","mother", "sex","affected")
extend_fam$Sample = rep('Extend', nrow(extend_fam))

# EuGEI
eugei_fam  <- read.table("/mnt/data1/GWAS_bloodcelltypes/MetaAnalysis_IVW/Ethnicity_Check/EuGEI_1kgIDs.fam",
                         stringsAsFactors = FALSE)
colnames(eugei_fam) <- c("pedigree","member","father","mother", "sex","affected")
eugei_fam$Sample = rep('EuGEI', nrow(eugei_fam))

# UnderSoc
undersoc_fam <- read.table("/mnt/data1/GWAS_bloodcelltypes/MetaAnalysis_IVW/Ethnicity_Check/UnderSoc_1kgIDs.fam",
                           stringsAsFactors = FALSE)
colnames(undersoc_fam) <- c("pedigree","member","father","mother", "sex","affected")
undersoc_fam$Sample = rep('UnderSoc', nrow(undersoc_fam))


# Combine all four datasets to label correctly
allcohort = rbind(extend_fam, eugei_fam,undersoc_fam)
allcohort2 = allcohort[,which(names(allcohort) %in% c("pedigree","Sample"))]
allcohort2 = allcohort2[,c(2,1)]
colnames(allcohort2)[colnames(allcohort2) == 'Sample'] <- 'Population'
colnames(allcohort2)[colnames(allcohort2) == 'pedigree'] <- 'Sample'
allcohort_1Kg = rbind(pop,allcohort2) #merge with 1KG dataset


# Merge the popualtion names with pcs files of pca eigenvalues
allcohort_1Kg$Sample = trimws(allcohort_1Kg$Sample)
pcs$V2 = trimws(pcs$V2)
pcs1 = cbind(pcs, allcohort_1Kg[match(pcs$V2, allcohort_1Kg$Sample), 'Population'])
colnames(pcs1)[colnames(pcs1)=="allcohort_1Kg[match(pcs$V2, allcohort_1Kg$Sample), c(\"Population\")]"] <- "Population"


##Rainbow runs out of color
#Remove NAs
library(ggplot2)
library(RColorBrewer)
myColors1 <- brewer.pal(8,"Dark2")
myColors2 <- brewer.pal(8,"Set2")
myColors3 <- brewer.pal(12,"Set3")
myColors4 <- brewer.pal(12,"Paired")
myColors = c(myColors1, myColors2, myColors3, myColors4)
names(myColors) <- levels(pcs1$Population)
colScale <- scale_colour_manual(name = "Population",values = myColors)


pdf("Ethinic check of All Cohorts using 1000 Genomes [ggplot].pdf")
ggplot(pcs1, aes(pcs1[,3], pcs1[,4], colour = pcs1$Population)) +
  geom_point()+
  labs(x = 'PC1', y = 'PC2') +
  colScale
dev.off()


#default R plot
plot(pcs1[,3], pcs1[,4], xlab = "PC1", ylab = "PC2",cex = .9,
     col = rainbow(nlevels(pcs1$Population))[pcs1$Population]) #third and fourth values are PC1 and Pc2 in matrix
legend("topright", levels(pcs1$Population), col = rainbow(nlevels(pcs1$Population)), pch = 10, cex=0.5)


