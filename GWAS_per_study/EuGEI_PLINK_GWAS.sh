
#### keep samples  ####
plink2 --bfile /mnt/data1/GWAS_bloodcelltypes/EuGEI/EUGEI_filtered_white --keep keep_whiteeuropean.txt --make-bed --out /mnt/data1/GWAS_bloodcelltypes/EuGEI/EUGEI_filtered_whitesamples

##Update sex info
plink2 --bfile /mnt/data1/GWAS_bloodcelltypes/EuGEI/EUGEI_filtered_white --update-sex EuGEI_sex_update.txt --make-bed --out /mnt/data1/GWAS_bloodcelltypes/EuGEI/EUGEI_filtered_white_sexupd


##GWAS with covariates
plink2 --bfile /mnt/data1/GWAS_bloodcelltypes/EuGEI/EUGEI_filtered_white_sexupd --pheno cell_phenotype.txt --all-pheno --covar cov_agesexpcs.txt --linear
