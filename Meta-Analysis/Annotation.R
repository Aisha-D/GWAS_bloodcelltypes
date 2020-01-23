## Check the VEP score of significant same directional snps
## Read in txt file of snps per cell type
library(tidyr)

cd8t = read.table("/mnt/data1/GWAS_bloodcelltypes/MetaAnalysis_IVW/METAANALYSIS_IVWBcell_1.txt",
                  stringsAsFactors = F, header = T)

cd8t2 = separate(data = cd8t, col = MarkerName, 
                 into = c("Chr", "Bp"), sep=":")
## Use BiomaRt
library(biomaRt)

snp_mart = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", 
                   path="/biomart/martservice", dataset="hsapiens_snp")
results <- c()
# Initialise storage vector

cd8t2$refsnp_id = rep(NA, nrow(cd8t2))
cd8t2$allele = rep(NA, nrow(cd8t2))
cd8t2$chrom_start = rep(NA, nrow(cd8t2))
cd8t2$chrom_strand = rep(NA, nrow(cd8t2))

for (snp in 1:5) { 
  #df is the object name of the SNP file

  temp <- getBM(attributes = c('refsnp_id', 'allele', 'chrom_start', 'chrom_strand'), 
                filters = c('chr_name', 'start', 'end'), 
                values = list(cd8t2[snp,1], cd8t2[snp, 2] , cd8t2[snp, 2]), 
                mart = snp_mart)
  
  temp[, 5] <- snp 
  # Store SNP file index row next to each answer
  
  results <- rbind(results, temp)
  
}


cd8t3 = cd8t2[1:10,1:2]
mart=useMart(biomart="ENSEMBL_MART_SNP", 
             host="grch37.ensembl.org", 
             path="/biomart/martservice", 
             dataset="hsapiens_snp")
system.time(do.call(rbind,(apply(cd8t3, 1, function (x) getBM(attributes = c('refsnp_id','chrom_start','chrom_end', 'chrom_strand','allele'), 
        filters = c('chr_name','start'), values = as.list(x), mart = mart)))))

cd8t3$end <- cd8t3$Bp
cd8t3$start <- cd8t3$Bp
cd8t3 = cd8t3[,-2]
cd8t3$ChrBp <- paste(cd8t3$Chr, cd8t3$start, sep = ":")
cd8t3$ChrBpstr <- paste(cd8t3$ChrBp, cd8t3$end, sep = ":")
cords <- list(cd8t3$ChrBpstr)

coords <- apply(cd8t3, 1, paste, collapse = ":")
coords3 <- apply(coords, 1, paste, collapse = ":")
coords2 <- as.vector(coords)
snp_mart = useMart(biomart="ENSEMBL_MART_SNP", 
                   path="/biomart/martservice", dataset="hsapiens_snp",
                   host="useast.ensembl.org")

rsIDs_1 <- getBM(attributes = c('refsnp_id','chrom_start','chrom_strand','allele'),
              filters = c("chr_name","start","end"), 
              values = list(as.numeric(cd8t3$Chr), 
                            as.numeric(cd8t3$start), 
                            as.numeric(cd8t3$end)), 
              mart = snp_mart)

rsIDs$ChrBp <- paste(rsIDs$chrom_strand, rsIDs$chrom_start, sep = ":")
cd8t3$ChrBp <- paste(cd8t3$Chr, cd8t3$start, sep = ":")


## loop to get every snp on the choromsome
coords4 = as.character(cd8t3$ChrBpstr)
names(coords4) <- 1:nrow(cd8t3)
coords5 <- coords4[1:3]
coords6 <- as.numeric(as.character(list(coords5)))

rsIDs_2 <- getBM(attributes = c('refsnp_id','allele','chrom_start','chrom_strand'), 
                   filters = c('chromosomal_region'), 
                   values = list(coords4),
                   mart = snp_mart)


#### 16/01/2020
listEnsembl()
variation = useEnsembl(biomart="snp")
listDatasets(variation)

variation = useEnsembl(biomart="snp", dataset="hsapiens_snp")
listFilters(variation)
listAttributes(variation)

rsIDs_2 <- getBM(attributes = c('refsnp_id','allele','chrom_start', 'chrom_end', 'chrom_strand'), 
                 filters = c('chromosomal_region'), 
                 values = list(cd8t4$Chr, cd8t4$start, cd8t4$end),
                 mart = variation)




ensembl <- useMart("ensembl")
listMarts(ensembl)

##Another way:https://www.biostars.org/p/323319/
cd8t4 <- cd8t3[,1:3]
cd8t4$Chr <- as.numeric(cd8t4$Chr)
cd8t4$end <- as.numeric(cd8t4$end)
cd8t4$start <- as.numeric(cd8t4$start)
mart=useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_snp")
tmp = do.call(rbind,(apply(cd8t4,1, function (x) getBM(attributes = c('refsnp_id','chrom_start','chrom_end', 'chrom_strand','allele'), 
                                                  filters = c('chr_name','start','end'), 
                                                  values = as.list(x), mart = mart))))
### Ensembl didn't work so lets us the UCSC outputs instead
cd8t2$alle <- paste(cd8t2$Allele1, cd8t2$Allele2, sep = "/")
cd8tv <- cd8t2[,which(names(cd8t2) %in% c("Chr", "Bp", "Bp", "alle"))]
cd8tv$Bp2 = cd8tv$Bp
cd8tv <- cd8tv[,c("Chr", "Bp", "Bp2", "alle")]
cd8tv2 <- cd8tv[1:25,]
write.table(cd8tv2, "CD8T_topsig.txt",row.names = F, sep = "\t", quote = FALSE)
