##Update sex info
plink2 --bfile /mnt/data1/EXTEND/Genotypes/Imputed/EXTEND_Unrelated_EUR_QCd --update-sex EXTEND_sex_update.txt --make-bed --out /mnt/data1/GWAS_bloodcelltypes/EXTEND/EXTEND_Unrelated_EUR_QCd_sex_upd

##GWAS with covariates

plink --bfile /mnt/data1/GWAS_bloodcelltypes/EXTEND/EXTEND_Unrelated_EUR_QCd_sex_upd --pheno cell_phenotype.txt --all-pheno 
--covar cov_agesexpcs.txt --linear 
