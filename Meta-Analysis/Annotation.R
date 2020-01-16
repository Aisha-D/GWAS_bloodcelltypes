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
                   host="grch37.ensembl.org", 
                   path="/biomart/martservice", dataset="hsapiens_snp")

rsIDs_1 <- getBM(attributes = c('refsnp_id','chrom_start','chrom_strand','allele'),
              filters = c("chr_name","start","end"), values = list(cd8t3$Chr, cd8t3$start, cd8t3$end), mart = snp_mart)

rsIDs$ChrBp <- paste(rsIDs$chrom_strand, rsIDs$chrom_start, sep = ":")
cd8t3$ChrBp <- paste(cd8t3$Chr, cd8t3$start, sep = ":")


## loop to get every snp on the choromsome
coords4 = cd8t3$ChrBpstr
names(coords4) <- 1:nrow(cd8t3)
rsIDs_2 <- getBM(attributes = c('refsnp_id','allele','chrom_start','chrom_strand'), 
                   filters = c('chromosomal_region'), 
                   values = coords4,
                   mart = snp_mart)


#### 16/01/2020
listEnsembl()
variation = useEnsembl(biomart="snp")
listDatasets(variation)

variation = useEnsembl(biomart="snp", dataset="hsapiens_snp")
listFilters(variation)
listAttributes(variation)
rsIDs_2 <- getBM(attributes = c('refsnp_id','allele','chrom_start','chrom_strand'), 
                 filters = c('chromosomal_region'), 
                 values = coords4,
                 mart = variation)




ensembl <- useMart("ensembl")
listMarts(ensembl)



### Ensembl didn't work so lets us the UCSC outputs instead
track = read.table("/mnt/data1/reference_files/hg38/All_SNPS_151/All_SNPs_hg19.txt",
                   stringsAsFactors = F, header = T)
