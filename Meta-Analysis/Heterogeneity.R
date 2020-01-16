## Plot forest plots for META-Analysis
library(dplyr)
library(qqman)
setwd("/mnt/data1/GWAS_bloodcelltypes/MetaAnalysis_IVW/")

## We will quickly look at CD8T
#Read in CD8T MetaAnalysis
cd8t = read.table("/mnt/data1/GWAS_bloodcelltypes/MetaAnalysis_IVW/METAANALYSIS_IVWCD8T_1.txt",
                  stringsAsFactors = F, header = T)

#Read in cd8t for each study
ext_cd8t = read.table("/mnt/data1/GWAS_bloodcelltypes/EXTEND/SEinc/plink.CD8T.assoc.linear.2", 
                    stringsAsFactors = F, header = T)

eugei_cd8t = read.table("/mnt/data1/GWAS_bloodcelltypes/EuGEI/SEinc/plink.CD8T.assoc.linear.2",
                   stringsAsFact= F, header = T)

undersoc_cd8t = read.table("/mnt/data1/GWAS_bloodcelltypes/Understanding_Society/SEinc/plink.CD8T.assoc.linear.2",
                      stringsAsFactors = F, header = T)

##Check heterohgenity using I2 

#First need to do Cochran Q test of extend?
library(nonpar)

tmp = full_join(extend_sig, eugei_sig,by = 'MarkerName')
tmp = full_joint(tmp, undersoc, by = 'MarkerName')

cochrans.q(c(extend$BETA, eugei$BETA))

## Follow tutorial:https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/fixed.html
#Look at the variance of cd8t zscore
pdf('CD8T_MetaAnalysis_variance.pdf')
hist(cd8t_het$Zscore, main = "Zscore Distribution CD8T", xlab = "Zscore")
dev.off()

#Look at the I scores of the top P.values
cd8t_het = cd8t_het[order(cd8t_het$P.value),]
head(cd8t_het)

#take ones where the symbols of direction all match
pos = '+++'
neg = '---'

cd8t_het_dir = cd8t_het[which(cd8t_het$Direction %in% c(pos, neg)),]
head(cd8t_het_dir) 
nrow(cd8t_het) -  nrow(cd8t_het_dir) # 7,596,855 snps didn't have same direction
plot(density(cd8t_het_dir$P.value))

# Select 10 most significant Pvalues
cd8t_het_dir = cd8t_het_dir[order(cd8t_het_dir$P.value),]
head(cd8t_het_dir)
cd8t_sig = cd8t_het_dir[1:10,1]

# Filter these snps from extend
extend_hetsig = extend_sig2[which(extend_sig2$MarkerName %in% cd8t_sig),]

# Plot STAT distribution
library(ggplot2)
ggplot(extend_hetsig, aes(extend_hetsig$MarkerName,extend_hetsig$BETA))+
  geom_jitter()

plot(extend_hetsig$BETA)
hist(extend_hetsig$BETA)


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


##########################################################################
### Comparing Neu and CD8T heterogeneity as qqplot for Neu was inflated ##
##########################################################################
## Read in PLINK data with SE
ex_cd8t = read.table("/mnt/data1/GWAS_bloodcelltypes/EXTEND/SEinc/plink.CD8T.assoc.linear",
                     stringsAsFactors = F, header = T)

#Look at the meta-analysis data more
cd8t_het = read.table("/mnt/data1/GWAS_bloodcelltypes/MetaAnalysis_Heterogeneity/METAANALYSIS_CD8T_1.txt",
                      stringsAsFactors = F, header = T)
neu_het = read.table("/mnt/data1/GWAS_bloodcelltypes/MetaAnalysis_Heterogeneity/METAANALYSIS_Neu_Het_1.txt",
                 stringsAsFactors = F, header = T)
neu_het = neu_het[order(neu_het$Zscore),]


## First Let us check ex_cd8t M value
library(getmstatistic)  # for calculating M statistics
library(gridExtra)      # for generating tables
library(metafor) # for conducting meta-analysis
library(dplyr)   # for sorting data.frames

ex_cd8t_top10 = as.data.frame(ex_cd8t[which(ex_cd8t$TEST == "ADD"),])
ex_cd8t_top10 = ex_cd8t_top10[order(ex_cd8t_top10$P),]
ex_cd8t_top10 = ex_cd8t_top10[1:10,]

getmstatistic_results <- getmstatistic(beta_in = ex_cd8t_top10$BETA,
                                       lambda_se_in = ex_cd8t_top10$SE,
                                       variant_names_in = ex_cd8t_top10$SNP,
                                       study_names_in = 'EXTEND')
head(getmstatistic_results)

#visualise results
dframe <- getmstatistic_results$M_dataset
head(dframe)


# Retrieve dataset of stronger than average studies (significant at 5% level)
getmstatistic_results$influential_studies_0_05

# Retrieve dataset of weaker than average studies (significant at 5% level)
getmstatistic_results$weaker_studies_0_05

# Retrieve number of studies and variants
getmstatistic_results$number_studies
getmstatistic_results$number_variants

# Retrieve expected mean, sd and critical M value at 5% significance level
getmstatistic_results$M_expected_mean
getmstatistic_results$M_expected_sd
getmstatistic_results$M_crit_alpha_0_05

# Sort getmstatistic_results dataframe by M statistics
getm_res_srtd <- dplyr::arrange(getmstatistic_results$M_dataset, M)
head(getm_res_srtd)


#### 13/01/2020 plot forest plot
cd8t_ex = read.table("/mnt/data1/GWAS_bloodcelltypes/EXTEND/SEinc/plink.CD8T.assoc.linear",
                     stringsAsFactors = F, header = T)
cd8t_eu = read.table("/mnt/data1/GWAS_bloodcelltypes/EuGEI/SEinc/plink.CD8T.assoc.linear",
                     stringsAsFactors = F, header = T)
cd8t_us = read.table("/mnt/data1/GWAS_bloodcelltypes/Understanding_Society/SEinc/plink.CD8T.assoc.linear",
                     stringsAsFactors = F, header = T)

extend_sig = cd8t_us
extend_sig$MarkerName = paste(extend_sig$CHR, extend_sig$BP, sep = ":")
extend_sig2 = as.data.frame(extend_sig[which(extend_sig$TEST == "ADD"),])
rownames(extend_sig2) <- extend_sig2$MarkerName
extend_sig3 =  extend_sig2[order(extend_sig2$P),]
cd8t_us <- extend_sig3

#Filter for most significant from meta-analysis
meta_sig = rbind(cd8t_ex[which(cd8t_ex$MarkerName == '3:186666337'),],
                 cd8t_eu[which(cd8t_eu$MarkerName == '3:186666337'),],
                 cd8t_us[which(cd8t_us$MarkerName == '3:186666337'),])
meta_sig$Study = c("Extend", "EuGEI", "UnderSoc")

##### ForestES function #######
#https://git.exeter.ac.uk/ejh243/ExeterEWASPipeline/blob/master/R/forestES.r
forestES<- function (x, data = NULL, xlim = NULL, ylim = NULL, xlab = NULL, 
                     ylab = NULL, pch = NULL, main = NULL, col = NULL, multiply = NULL, cex = NULL, 
                     names = NULL, poly = NULL, polycol = NULL, sepline = NULL, seplinecol = NULL,
                     septitle = NULL, septitlepos = NULL, septitleposx = NULL, namesx = NULL){
  
  # 'data' must be list of 1) dataframe of effect sizes 2) dataframe of 
  # standard errors with columns as analyses/cohorts and rows of cpgs or locations. 
  # Columns should be in same order. Will put first column on bottom and last 
  # column on top
  # x is cpg site (must be rownames)
  # multiply is multiplication factor added to ES and SE
  # names are what lines of forest should be labelled
  # poly are lines you want to be represented by diamonds rather than points 
  # and lines
  # polycol what colours do you want the poly to be (from bottom to top)
  # sepline is where you want to put separation lines (0.5 between points)
  # septitle is section labels for separation
  # septitlepos where do you want section labels
  # namesx position of names on x axis (is autimatically calculated but may 
  # need adjustment)
  
  Est<-as.numeric(data[[1]][x,])
  
  SE<-as.numeric(data[[2]][x,])
  
  if(is.null(ylim)) {
    ylim=c(0.6,length(Est))
  }
  
  if(is.null(xlab)) {
    xlab="Effect Size"
  }
  
  if(is.null(ylab)) {
    ylab=""
  }	
  
  if(is.null(pch)) {
    pch = 19 
  }
  
  if(is.null(main)) {
    main = x
  }
  
  if(is.null(col)) {
    col = "black"
  }
  
  if(is.null(multiply)) {
    multiply = 1
  }
  
  if(is.null(seplinecol)) {
    seplinecol = "gray85"
  }
  
  if(is.null(names)) {
    names=colnames(data[[1]])
  }
  
  if(is.null(poly)) {
    poly=NULL
  }
  
  if(is.null(polycol)) {
    polycol=replicate(length(poly),"black")
  }
  
  Est*multiply->Est2
  SE*multiply->SE2
  
  if(is.null(xlim)) {
    xlim <- c(min(Est2-SE2)-0.05,max(Est2+SE2)+0.02)
  }
  
  if(is.null(septitleposx)) {
    septitleposx = min(xlim)
  }
  
  if(is.null(namesx)) {
    namesx<-min(xlim) + 0.02
  }
  
  par(xpd=TRUE)
  
  plot(1:length(Est2)~Est2, xlim=xlim, axes=F, ylim=ylim, pch=pch, main=main, 
       col=col, cex=cex, ylab=ylab, xlab=xlab)
  segments((Est2 - SE2), 1:length(Est2), (Est2 + SE2), 1:length(Est2), col=col)
  abline(v=0, xpd=FALSE)
  axis(1)
  
  for(i in 1:length(poly)){
    polygon(c(Est2[poly[i]] - SE2[poly[i]], Est2[poly[i]],Est2[poly[i]] + 
                SE2[poly[i]],Est2[poly[i]]), c(poly[i], poly[i]-0.2, poly[i], poly[i]+0.2), 
            col=polycol[i], border=polycol[i])
  }
  
  for (i in 1:length(names)){
    text(namesx, i, names[i])
  }
  
  for (i in 1:length(sepline)){
    abline(h=sepline[i], xpd=FALSE, col=seplinecol)
  }
  
  for (i in 1:length(septitle)){
    text(septitleposx, septitlepos[i], septitle[i], srt=90, cex=1)
  }
  
}
###### Plot forest ######
ES <- meta_sig[1:3,"BETA"]
ES <- as.data.frame(t(ES))
colnames(ES) <-c("EXTEND", "EuGEI", "Understanding_Society")
rownames(ES) <- 'cpg1'

SE <- meta_sig[1:3,"SE"]
SE <- as.data.frame(t(SE))
colnames(SE) <-c("EXTEND", "EuGEI", "Understanding_Society")
rownames(SE) <- 'cpg1'

dat <- list(ES = ES, SE = SE)
pdf("forest plot cd8t.pdf")
forestES(x = 'cpg1', data = dat, xlim = c(-0.02,0.02), namesx = -0.02,
         main = '3:186666337 (Most Significant CD8T SNP)',
         names = c("Extend(n=983)", "EuGEI(n=671)", "UnderSoc(n=1112)"))
dev.off()


