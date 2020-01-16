## Plot forest plots for META-Analysis
library(dplyr)
#Read in cd8t MetaAnalysis
cd8t = read.table("/mnt/data1/GWAS_bloodcelltypes/MetaAnalysis/METAANALYSIS_CD8T_1.txt",
                    stringsAsFactors = F, header = T)

#Read in cd8t for each study
extend = read.table("/mnt/data1/GWAS_bloodcelltypes/EXTEND/plink.CD8T.assoc.linear", 
                    stringsAsFactors = F, header = T)

eugei = read.table("/mnt/data1/GWAS_bloodcelltypes/EuGEI/plink.CD8T.assoc.linear",
                   stringsAsFact= F, header = T)

undersoc = read.table("/mnt/data1/GWAS_bloodcelltypes/Understanding_Society/plink.CD8T.assoc.linear",
                      stringsAsFactors = F, header = T)

## Filter for top significant p.values from MetaAnalysis
cd8t = cd8t[order(cd8t$P.value),]
sig_cd8t <- cd8t[1:5,]
lssig_cd8t <- sig_cd8t[,1]

##Take the SE from each of the cogs from sig-cd8t list
extend_sig = extend
extend_sig$MarkerName = paste(extend_sig$CHR, extend_sig$BP, sep = ":")
extend_sig2 = as.data.frame(extend_sig[which(extend_sig$TEST == "ADD"),])
extend_sig3 = as.data.frame(extend_sig2[extend_sig2$MarkerName %in% lssig_cd8t, c(7,8,10)])
rownames(extend_sig3) <- extend_sig3$MarkerName
colnames(extend_sig3)[colnames(extend_sig3) == 'STAT'] <- 'EXTEND'

eugei_sig = eugei
eugei_sig$MarkerName = eugei$SNP
eugei_sig2 = as.data.frame(eugei_sig[which(eugei_sig$TEST == "ADD"),])
eugei_sig3 = as.data.frame(eugei_sig2[eugei_sig2$MarkerName %in% lssig_cd8t, c(7,8,10)])
rownames(eugei_sig3) <- eugei_sig3$MarkerName
colnames(eugei_sig3)[colnames(eugei_sig3) == 'STAT'] <- 'EuGEI'

undersoc_sig = undersoc
undersoc_sig$MarkerName = paste(undersoc_sig$CHR, undersoc_sig$BP, sep = ":")
undersoc_sig2 = as.data.frame(undersoc_sig[which(undersoc_sig$TEST == "ADD"),])
undersoc_sig3 = as.data.frame(undersoc_sig2[undersoc_sig2$MarkerName %in% lssig_cd8t, c(7,8,10)])
rownames(undersoc_sig3) <- undersoc_sig3$MarkerName
colnames(undersoc_sig3)[colnames(undersoc_sig3) == 'STAT'] <- 'Understanding_Society'

#combine the three studies z-scores in one
allstudies = full_join(extend_sig3, eugei_sig3,by = 'MarkerName')
allstudies = merge(allstudies, undersoc_sig3, by = 'MarkerName')

allstudies = merge(allstudies, sig_cd8t[,5], by = 'MarkerName')
rownames(allstudies) <- allstudies$MarkerName
allstudies = allstudies[,-1]
allstudies = allstudies[1,]

allstudies


ES <- allstudies[,c("BETA.x", "BETA.y", "BETA")]
colnames(ES) <-c("EXTEND", "EuGEI", "Understanding_Society")

SE <- allstudies[,c("EXTEND", "EuGEI", "Understanding_Society")]

dat <- list(ES = ES, SE = SE)


library(forestplot)
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
pdf("forest plot1.pdf")
forestES(x = rownames(allstudies[1,]), data = dat)
dev.off()

##Check heterohgenity using I2 

#First need to do Cochran Q test of extend?

















