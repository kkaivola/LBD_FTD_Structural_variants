#!/bin/bash

#Make a new directory
mkdir /data/ALS_50k/karri/LBD_FTD_SV_project/beta_beta_plot

##Download Belenguez harmonized summary statistics

wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/harmonised/35379992-GCST90027158-MONDO_0004975.h.tsv.gz
gunzip  35379992-GCST90027158-MONDO_0004975.h.tsv.gz


#Extract analyzed SNPs/indels from the TPCN1 LD block in LBD data
cd /data/ALS_50k/karri/LBD_FTD_SV_project/beta_beta_plot

##TPCN1 haplotype region, results from all variants, hwe filtered MAF controls > 1%
awk '($1=="#CHROM" || ($1==12 && $2>113150000 && $2<113350000 && $3 !~/LBD/))' /data/ALS_50k/DementiaSeq_SV.LBD.FTD/merged.SV.SNV_1perc_fdr/Analysis.GLM.hg38_minGQ300/LBD/LBD.controls.UNRELATED.SNVindels.SV.glm.hwe1e-6.maf001.ALLchr.txt > TPCN1_analyzed_snps_indels_glm.txt

##Get IDs
awk '(NR>1) {print $3}' TPCN1_analyzed_snps_indels_glm.txt > TPCN1_analyzed_snps_indels.txt

##Format to match Belenguez data, i.e. change : separator to _ and remove chr prefix
sed -i 's/chr//g' TPCN1_analyzed_snps_indels.txt
sed -i 's/:/_/g' TPCN1_analyzed_snps_indels.txt

#Grep TPCN1 haplotype SNPs indentified in LBD data from harmonized Belenguez summary stats
cd /data/ALS_50k/karri/LBD_FTD_SV_project/beta_beta_plot

head -n1 35379992-GCST90027158-MONDO_0004975.h.tsv > 35379992-GCST90027158-MONDO_0004975.h._TPCN1_SNPs.tsv
grep -Ff TPCN1_analyzed_snps_indels.txt 35379992-GCST90027158-MONDO_0004975.h.tsv >> 35379992-GCST90027158-MONDO_0004975.h._TPCN1_SNPs.tsv

##Check how many were found (410 out of 426)
echo "number of searched SNPs:"
wc -l TPCN1_analyzed_snps_indels.txt
echo "Number of found SNPs"
wc -l 35379992-GCST90027158-MONDO_0004975.h._TPCN1_SNPs.tsv
