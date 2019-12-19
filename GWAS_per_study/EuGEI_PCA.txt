LD prune??
plink2 --bfile EUGEI_GROUP_STrios_TrulyFinal_XVII_BestGuess_filtered --indep 50 5 1.5 --out EUGEI_GROUP_STrios_TrulyFinal_XVII_BestGuess_filtered.ld
plink2 --bfile EUGEI_GROUP_STrios_TrulyFinal_XVII_BestGuess_filtered --extract EUGEI_GROUP_STrios_TrulyFinal_XVII_BestGuess_filtered.ld.prune.in --make-bed --out EUGEI_GROUP_STrios_TrulyFinal_XVII_BestGuess_filtered.ld.prune

# use GCTA to calc PCs
/mnt/data1/programs/gcta_v24_installation/gcta64 --bfile EUGEI_GROUP_STrios_TrulyFinal_XVII_BestGuess_filtered.ld.prune --make-grm-bin --autosome --out EUGEI_GROUP_STrios_TrulyFinal_XVII_BestGuess_filtered
/mnt/data1/programs/gcta_v24_installation/gcta64 --grm EUGEI_GROUP_STrios_TrulyFinal_XVII_BestGuess_filtered --pca --out EUGEI_GROUP_STrios_TrulyFinal_XVII_BestGuess_filtered.pca



####KEEEP some imputed samples that were in the pheno file
plink2 --bfile /mnt/data1/EuGEI/Genotypes/Imputed/EUGEI_GROUP_STrios_TrulyFinal_XVII_BestGuess --keep keep.txt --make-bed --out EUGEI_GROUP_STrios_TrulyFinal_XVII_BestGuess_filtered

####Keep european samples then visualise it as a pca
plink2 --bfile EUGEI_GROUP_STrios_TrulyFinal_XVII_BestGuess_filtered --keep keep_european.txt --make-bed --out EUGEI_filtered_white