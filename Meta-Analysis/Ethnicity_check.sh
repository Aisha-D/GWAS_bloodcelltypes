## This uses the 1KG to compare ethnicity samples from all three cohorts(EXTEND, EuGEI, Understanding Society)
cd /mnt/data1/GWAS_bloodcelltypes/MetaAnalysis_IVW/Ethnicity_Check/

##############   Cchange variant ids to chr:bp for each cohort    ##################
awk '{if ($1 != 0) print $2,"chr"$1":"$4}' /mnt/data1/GWAS_bloodcelltypes/EXTEND/QC/EXTEND_Unrelated_EUR_QCd_sex_upd_6.bim > updateTo1KGFormatEx.txt

awk '{if ($1 != 0) print $2,"chr"$1":"$4}' /mnt/data1/GWAS_bloodcelltypes/EuGEI/QC/EUGEI_filtered_white_sexupd_6.bim > updateTo1KGFormatEu.txt

awk '{if ($1 != 0) print $2,"chr"$1":"$4}' /mnt/data1/GWAS_bloodcelltypes/Understanding_Society/QC/data_filtered_1_sex_upd_5.bim > updateTo1KGFormatUS.txt

#############    Change names to correct format so there it is easier to match    #################
plink2 --bfile /mnt/data1/GWAS_bloodcelltypes/EXTEND/QC/EXTEND_Unrelated_EUR_QCd_sex_upd_6 --update-name updateTo1KGFormatEx.txt --make-bed --out Extend_1kgIDs

plink2 --bfile /mnt/data1/GWAS_bloodcelltypes/EuGEI/QC/EUGEI_filtered_white_sexupd_6 --update-name updateTo1KGFormatEu.txt --make-bed --out EuGEI_1kgIDs

plink2 --bfile /mnt/data1/GWAS_bloodcelltypes/Understanding_Society/QC/data_filtered_1_sex_upd_5 --update-name updateTo1KGFormatUS.txt --make-bed --out UnderSoc_1kgIDs


############    Merge Bed files (all cohorts + 1000Genomes)    ##########################
## To compare all three ethnic backgrounds we need to merge the bed files

#EXTEND and EuGEI will be merged first
plink2 --bfile Extend_1kgIDs --bmerge EuGEI_1kgIDs --maf 0.1 --geno 0.1 --make-bed --out EXTEND_EuGEI

#EXTEND_EuGEI will mege with UnderSoc
plink2 --bfile EXTEND_EuGEI --bmerge UnderSoc_1kgIDs --maf 0.1 --geno 0.1 --make-bed --out EXTEND_EuGEI_UnderSoc


#UnderSoc bed files have variants with 3+ alleles present. Remove these snps from UnderSoc Bed files.
plink2 --bfile EXTEND_EuGEI --exclude EXTEND_EuGEI_UnderSoc-merge.missnp  --make-bed --out EXTEND_EUGEI_exludetriallelic

plink2 --bfile UnderSoc_1kgIDs --exclude EXTEND_EuGEI_UnderSoc-merge.missnp  --make-bed --out UnderSoc_1kgIDs_exludetriallelic

plink2 --bfile EXTEND_EUGEI_exludetriallelic --bmerge UnderSoc_1kgIDs_exludetriallelic -make-bed  --maf 0.1 --geno 0.1 --out EXTEND_EuGEI_UnderSoc


#1000 genomes merge
plink2 --bfile EXTEND_EuGEI_UnderSoc --bmerge /mnt/data1/reference_files/1000G/Phase3/1000G_gr38 --maf 0.2 --geno 0.05 --make-bed --out AllCohort_Merged1KG

plink2 --bfile /mnt/data1/reference_files/1000G/Phase3/1000G_gr38 --exclude AllCohort_Merged1KG-merge.missnp --make-bed --out 1000G_gr38

plink2 --bfile EXTEND_EuGEI_UnderSoc --bmerge 1000G_gr38 --maf 0.2 --geno 0.05 --make-bed --out AllCohort_Merged1KG


# LD prune
plink2 --bfile AllCohort_Merged1KG --indep 50 5 1.5 --out AllCohort_Merged1KG.ld
plink2 --bfile AllCohort_Merged1KG --extract AllCohort_Merged1KG.ld.prune.in --make-bed --out AllCohort_Merged1KG.ld.prune

# use GCTA to calc PCs
/mnt/data1/programs/gcta_v24_installation/gcta64 --bfile AllCohort_Merged1KG.ld.prune --make-grm-bin --autosome --out AllCohort_Merged1KG
/mnt/data1/programs/gcta_v24_installation/gcta64 --grm AllCohort_Merged1KG --pca --out AllCohort_Merged1KG.pca


#PCA for 3 cohorts only
# LD prune
plink2 --bfile EXTEND_EuGEI_UnderSoc --indep 50 5 1.5 --out EXTEND_EuGEI_UnderSoc.ld
plink2 --bfile EXTEND_EuGEI_UnderSoc --extract EXTEND_EuGEI_UnderSoc.ld.prune.in --make-bed --out EXTEND_EuGEI_UnderSoc.ld.prune

# use GCTA to calc PCs
/mnt/data1/programs/gcta_v24_installation/gcta64 --bfile EXTEND_EuGEI_UnderSoc.ld.prune --make-grm-bin --autosome --out EXTEND_EuGEI_UnderSoc
/mnt/data1/programs/gcta_v24_installation/gcta64 --grm EXTEND_EuGEI_UnderSoc --pca --out EXTEND_EuGEI_UnderSoc.pca
