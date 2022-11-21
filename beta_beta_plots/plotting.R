#Set wd
setwd("/data/ALS_50k/karri/LBD_FTD_SV_project/beta_beta_plot")

#Load libraries
library(tidyverse)
library(ggplot2)
library(scales)

#Import data
our <- read.delim("TPCN1_analyzed_snps_indels_glm.txt")
belenguez <- read.delim("35379992-GCST90027158-MONDO_0004975.h._TPCN1_SNPs.tsv")

#Harmonize SNP IDs, i.e. change : separator to _ and remove chr prefix in our data
our$ID <- gsub("chr","",our$ID)
our$ID <- gsub(":","_",our$ID)

#Merge dfs by SNP ID
df <- merge(our, belenguez, by.x="ID", by.y="hm_variant_id", how="union")

#Harmonize allele for which effect was tested
print(table(df$A1 == df$effect_allele))
df$tmp_allele <- df$A1
df$tmp_allele[df$A1 != df$effect_allele] <- df$REF[df$A1 != df$effect_allele]
print(table(df$tmp_allele == df$effect_allele))
###All match now, so no strand problems

#For variants where A1 is not effect allele, make belenguez beta *-1
df$beta[df$A1 != df$effect_allele] <- -1*(df$beta[df$A1 != df$effect_allele])

#Calculate beta for our LBD data ORs
df$beta_our <- log(df$OR)

#Create variable that shows only suggestive p value variants
df$suggestive <- ifelse((df$P < 0.00001 | df$p_value < 0.00001), TRUE, FALSE)


##Correlation tests
cor1 <- cor.test(df$beta_our, df$beta, 
                 method = "pearson")
print(cor1)

#Plot
print(ggplot(df, aes(x=beta, y=beta_our)) +
        geom_point(aes(color=suggestive)) + 
        scale_color_manual(values = c('black', 'red')) +
        theme_bw() + 
        theme(legend.position="none") +
        geom_smooth(se = T, method = lm) + 
        geom_vline(xintercept = 0, linetype = 2) + 
        geom_hline(yintercept = 0, linetype = 2) +
        xlab("TPCN1 locus beta in AD") +
        ylab("TPCN1 locus beta in LBD") +
        ylim(-0.4,0.4))
