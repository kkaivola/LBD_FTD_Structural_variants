library("Gviz")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("tidyverse")

setwd("~/Documents/LBD_FTD_project_2021")

#Import files and do formatting
LBD <- read.delim("LBD_all_info_unfiltered_with_CADDSV_genehancer.txt")
FTD <- read.delim("FTD_all_info_unfiltered_with_CADDSV_genehancer.txt")
genes <- read.csv("NDD_ensembl_104_canonical_gene_regions.txt", sep="", stringsAsFactors=FALSE)

LBD$CHROM <- paste0("chr",LBD$CHROM)
FTD$CHROM <- paste0("chr", FTD$CHROM)

Dementia_seq_LBD_case_IDs <- read.table("Dementia_seq_LBD_case_IDs", quote="\"", comment.char="", stringsAsFactors=FALSE)
Dementia_seq_LBD_case_IDs <- Dementia_seq_LBD_case_IDs$V1
Dementiaseq_LBD_controls_IDs <- read.table("Dementiaseq_LBD_controls_IDs", quote="\"", comment.char="", stringsAsFactors=FALSE)
Dementiaseq_LBD_controls_IDs <- Dementiaseq_LBD_controls_IDs$V1

Dementia_seq_FTD_case_IDs <- read.table("Dementiaseq_FTD_cases_IDs", quote="\"", comment.char="", stringsAsFactors=FALSE)
Dementia_seq_FTD_case_IDs <- Dementia_seq_FTD_case_IDs$V1
Dementiaseq_FTD_controls_IDs <- read.table("Dementiaseq_FTD_control_IDs", quote="\"", comment.char="", stringsAsFactors=FALSE)
Dementiaseq_FTD_controls_IDs <- Dementiaseq_FTD_controls_IDs$V1


#Convert carrier information from factor to chracter, then remove brackets, apostrophes and genotypes
LBD$minor_allele_samples <- gsub(pattern = "[][]", replacement = "", LBD$minor_allele_samples)
LBD$minor_allele_samples <- gsub(pattern = "'", replacement = "", LBD$minor_allele_samples)
LBD$minor_allele_samples <- gsub(pattern = ":0/0", replacement = "", LBD$minor_allele_samples)
LBD$minor_allele_samples <- gsub(pattern = ":0/1", replacement = "", LBD$minor_allele_samples)
LBD$minor_allele_samples <- gsub(pattern = ":1/1", replacement = "", LBD$minor_allele_samples)
LBD$minor_allele_samples <- strsplit(LBD$minor_allele_samples, ",")

FTD$minor_allele_samples <- gsub(pattern = "[][]", replacement = "", FTD$minor_allele_samples)
FTD$minor_allele_samples <- gsub(pattern = "'", replacement = "", FTD$minor_allele_samples)
FTD$minor_allele_samples <- gsub(pattern = ":0/0", replacement = "", FTD$minor_allele_samples)
FTD$minor_allele_samples <- gsub(pattern = ":0/1", replacement = "", FTD$minor_allele_samples)
FTD$minor_allele_samples <- gsub(pattern = ":1/1", replacement = "", FTD$minor_allele_samples)
FTD$minor_allele_samples <- strsplit(FTD$minor_allele_samples, ",")

LBD$ALT <- gsub(pattern = ":", replacement = "_", LBD$ALT)
FTD$ALT <- gsub(pattern = ":", replacement = "_", FTD$ALT)

#Compute number of case and control carriers per SV

for (i in 1:nrow(LBD)) {
  IDs <- unlist(LBD$minor_allele_samples[i])
  IDs <- gsub(" ", "", IDs)
  LBD$n_case_carriers[i] <- length(intersect(IDs, Dementia_seq_LBD_case_IDs))
  LBD$n_control_carriers[i] <- length(intersect(IDs, Dementiaseq_LBD_controls_IDs))
}

for (i in 1:nrow(FTD)) {
  IDs <- unlist(FTD$minor_allele_samples[i])
  IDs <- gsub(" ", "", IDs)
  FTD$n_case_carriers[i] <- length(intersect(IDs, Dementia_seq_FTD_case_IDs))
  FTD$n_control_carriers[i] <- length(intersect(IDs, Dementiaseq_FTD_controls_IDs))
}


#make all SV END >= POS, since it doesn't matter in visualisation and otherwise need extra complication in
#visualisation function where would test if END > POS (gviz errors if from > to)
print("LBD")
for (i in 1:nrow(LBD)) {
  if (LBD$END[i] >= LBD$POS[i]) {
    next
  } else {
    print(paste("Initial pos:", LBD$POS[i], "and initial end:", LBD$END[i]))
    tmp <- LBD$END[i]
    LBD$END[i] <- LBD$POS[i]
    LBD$POS[i] <- tmp
    print(paste("Updated pos:", LBD$POS[i], "and updated end:", LBD$END[i]))
  }
}


print("FTD")
for (i in 1:nrow(FTD)) {
  if (FTD$END[i] >= FTD$POS[i]) {
    next
  } else {
    print(paste("Initial pos:", FTD$POS[i], "and initial end:", FTD$END[i]))
    tmp <- FTD$END[i]
    FTD$END[i] <- FTD$POS[i]
    FTD$POS[i] <- tmp
    print(paste("Updated pos:", FTD$POS[i], "and updated end:", FTD$END[i]))
  }
}


#Make info column that has count of carriers and in parentheses QUAL value
LBD$info_case <- paste0(LBD$n_case_carriers, " (", LBD$QUAL, ")")
LBD$info_control <- paste0(LBD$n_control_carriers, " (", LBD$QUAL, ")")

FTD$info_case <- paste0(FTD$n_case_carriers, " (", FTD$QUAL, ")")
FTD$info_control <- paste0(FTD$n_control_carriers, " (", FTD$QUAL, ")")

#Run Gviz
##LBD candidate gene exonic overlaps
for (i in 1:nrow(genes)) {
  gene = genes$Gene[i]
  start = genes$txStart[i] -1000
  stop = genes$txEnd[i] + 1000
  
  LBD_cases <- LBD %>%
    filter(GENE == gene) %>%
    filter(LBD_ALT_FREQS > 0) %>%
    filter(candidate_gene_exon_overlap == "Exonic") 
  LBD_controls <- LBD %>%
    filter(GENE == gene) %>%
    filter(controls_ALT_FREQS > 0) %>%
    filter(candidate_gene_exon_overlap == "Exonic") 
  
  LBD_cases$POS[LBD_cases$POS < start] <- start - (0.2*(stop-start))
  LBD_cases$END[LBD_cases$END > stop] <- stop + (0.2*(stop-start))
  LBD_controls$POS[LBD_controls$POS < start] <- start - (0.2*(stop-start))
  LBD_controls$END[LBD_controls$END > stop] <- stop + (0.2*(stop-start))
  
  png(paste0("NIH_LBD_FTD_SV_shinyapp/NIH_SV_images/LBD/",gene,"_LBD_exonic.png"), height = 1180, width = 1180, units = "px", res=100)
  
  if ((nrow(LBD_cases) > 0) & (nrow(LBD_controls) > 0)) {
    gtrack <- GenomeAxisTrack()
    itrack <- IdeogramTrack(genome = "hg38", chromosome = LBD_cases$CHROM[1])
    grtrack <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg38.knownGene, genome = "hg38", chromosome = LBD_cases$CHROM[1], name = gene, transcriptAnnotation = "symbol", size=0.4)
    atrack_LBD <- AnnotationTrack(start=LBD_cases$POS, end=LBD_cases$END, chromosome = LBD_cases$CHROM[1], genome="hg38", name="LBD", feature = LBD_cases$ALT, id=LBD_cases$info_case, showFeatureId = TRUE, fontcolor.item="black")
    atrack_controls <- AnnotationTrack(start=LBD_controls$POS, end=LBD_controls$END, chromosome = LBD_controls$CHROM[1], genome="hg38", name="controls", feature = LBD_controls$ALT, id=LBD_controls$info_control, showFeatureId = TRUE, fontcolor.item="black")
    plotTracks(list(itrack,gtrack, atrack_LBD, atrack_controls, grtrack), from = start, to = stop, fontcolor.title="#808080", extend.right = 0.2, extend.left = 0.2, fontsize.group=10, fontsize.title=10, BND="#377246", CPX="#86E496", DEL="#D6402D", DUP="#5282AF", INS="#D579E1", INS_ME="#F6AAFF", INS_ME_ALU="#F6AAFF", INS_ME_LINE1="#F6AAFF", INS_ME_SVA="#F6AAFF", INV="#F69A57")
  } else if ((nrow(LBD_cases) > 0) & (nrow(LBD_controls) == 0)) {
    gtrack <- GenomeAxisTrack()
    itrack <- IdeogramTrack(genome = "hg38", chromosome = LBD_cases$CHROM[1])
    grtrack <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg38.knownGene, genome = "hg38", chromosome = LBD_cases$CHROM[1], name = gene, transcriptAnnotation = "symbol", size=0.4)
    atrack_LBD <- AnnotationTrack(start=LBD_cases$POS, end=LBD_cases$END, chromosome = LBD_cases$CHROM, genome="hg38", name="LBD", feature = LBD_cases$ALT, id=LBD_cases$info_case, showFeatureId = TRUE, fontcolor.item="black")
    plotTracks(list(itrack,gtrack, atrack_LBD, grtrack), from = start, to = stop, fontcolor.title="#808080", extend.right = 0.2, extend.left = 0.2, fontsize.group=10, fontsize.title=10, BND="#377246", CPX="#86E496", DEL="#D6402D", DUP="#5282AF", INS="#D579E1", INS_ME="#F6AAFF", INS_ME_ALU="#F6AAFF", INS_ME_LINE1="#F6AAFF", INS_ME_SVA="#F6AAFF", INV="#F69A57")
  } else if ((nrow(LBD_cases) == 0) & (nrow(LBD_controls) == 0)) {
    gtrack <- GenomeAxisTrack(showTitle=TRUE, name="No SVs", fontcolor.title="#808080", rotation.title=0)
    plotTracks(gtrack, from = 100, to = 110, title.width = 5)
  } else {
    gtrack <- GenomeAxisTrack()
    itrack <- IdeogramTrack(genome = "hg38", chromosome = LBD_controls$CHROM[1])
    grtrack <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg38.knownGene, genome = "hg38", chromosome = LBD_controls$CHROM[1], name = gene, transcriptAnnotation = "symbol", size=0.4)
    atrack_controls <- AnnotationTrack(start=LBD_controls$POS, end=LBD_controls$END, chromosome = LBD_controls$CHROM, genome="hg38", name="controls", feature = LBD_controls$ALT, id=LBD_controls$info_control, showFeatureId = TRUE, fontcolor.item="black")
    plotTracks(list(itrack,gtrack, atrack_controls, grtrack), from = start, to = stop, fontcolor.title="#808080",extend.right = 0.2, extend.left = 0.2, fontsize.group=10, fontsize.title=10, BND="#377246", CPX="#86E496", DEL="#D6402D", DUP="#5282AF", INS="#D579E1", INS_ME="#F6AAFF", INS_ME_ALU="#F6AAFF", INS_ME_LINE1="#F6AAFF", INS_ME_SVA="#F6AAFF", INV="#F69A57")
  }
  dev.off()
}

##FTD exonic 
for (i in 1:nrow(genes)) {
  gene = genes$Gene[i]
  start = genes$txStart[i] -1000
  stop = genes$txEnd[i] + 1000
  
  FTD_cases <- FTD %>%
    filter(GENE == gene) %>%
    filter(FTD_ALT_FREQS > 0) %>%
    filter(Candidate_gene_exon_overlap == "Exonic")
  
  FTD_controls <- FTD %>%
    filter(GENE == gene) %>%
    filter(controls_ALT_FREQS > 0) %>%
    filter(Candidate_gene_exon_overlap == "Exonic")
  
  FTD_cases$POS[FTD_cases$POS < start] <- start - (0.2*(stop-start))
  FTD_cases$END[FTD_cases$END > stop] <- stop + (0.2*(stop-start))
  FTD_controls$POS[FTD_controls$POS < start] <- start - (0.2*(stop-start))
  FTD_controls$END[FTD_controls$END > stop] <- stop + (0.2*(stop-start))
  
  png(paste0("NIH_LBD_FTD_SV_shinyapp/NIH_SV_images/FTD/",gene,"_FTD_exonic.png"), height = 1180, width = 1180, units = "px", res=100)
  
  if ((nrow(FTD_cases) > 0) & (nrow(FTD_controls) > 0)) {
    gtrack <- GenomeAxisTrack()
    itrack <- IdeogramTrack(genome = "hg38", chromosome = FTD_cases$CHROM[1])
    grtrack <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg38.knownGene, genome = "hg38", chromosome = FTD_cases$CHROM[1], name = gene, transcriptAnnotation = "symbol", size=0.4)
    atrack_FTD <- AnnotationTrack(start=FTD_cases$POS, end=FTD_cases$END, chromosome = FTD_cases$CHROM[1], genome="hg38", name="FTD", feature = FTD_cases$ALT, id=FTD_cases$info_case, showFeatureId = TRUE, fontcolor.item="black")
    atrack_controls <- AnnotationTrack(start=FTD_controls$POS, end=FTD_controls$END, chromosome = FTD_controls$CHROM[1], genome="hg38", name="controls", feature = FTD_controls$ALT, id=FTD_controls$info_control, showFeatureId = TRUE, fontcolor.item="black")
    plotTracks(list(itrack,gtrack, atrack_FTD, atrack_controls, grtrack), from = start, to = stop, fontcolor.title="#808080", extend.right = 0.2, extend.left = 0.2, fontsize.group=10, fontsize.title=10, BND="#377246", CPX="#86E496", DEL="#D6402D", DUP="#5282AF", INS="#D579E1", INS_ME="#F6AAFF", INS_ME_ALU="#F6AAFF", INS_ME_LINE1="#F6AAFF", INS_ME_SVA="#F6AAFF", INV="#F69A57")
  } else if ((nrow(FTD_cases) > 0) & (nrow(FTD_controls) == 0)) {
    gtrack <- GenomeAxisTrack()
    itrack <- IdeogramTrack(genome = "hg38", chromosome = FTD_cases$CHROM[1])
    grtrack <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg38.knownGene, genome = "hg38", chromosome = FTD_cases$CHROM[1], name = gene, transcriptAnnotation = "symbol", size=0.4)
    atrack_FTD <- AnnotationTrack(start=FTD_cases$POS, end=FTD_cases$END, chromosome = FTD_cases$CHROM, genome="hg38", name="FTD", feature = FTD_cases$ALT, id=FTD_cases$info_case, showFeatureId = TRUE, fontcolor.item="black")
    plotTracks(list(itrack,gtrack, atrack_FTD, grtrack), from = start, to = stop, fontcolor.title="#808080", extend.right = 0.2, extend.left = 0.2, fontsize.group=10, fontsize.title=10, BND="#377246", CPX="#86E496", DEL="#D6402D", DUP="#5282AF", INS="#D579E1", INS_ME="#F6AAFF", INS_ME_ALU="#F6AAFF", INS_ME_LINE1="#F6AAFF", INS_ME_SVA="#F6AAFF", INV="#F69A57")
  } else if ((nrow(FTD_cases) == 0) & (nrow(FTD_controls) == 0)) {
    gtrack <- GenomeAxisTrack(showTitle=TRUE, name="No SVs", fontcolor.title="#808080", rotation.title=0)
    plotTracks(gtrack, from = 100, to = 110, title.width = 5)
  } else {
    gtrack <- GenomeAxisTrack()
    itrack <- IdeogramTrack(genome = "hg38", chromosome = FTD_controls$CHROM[1])
    grtrack <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg38.knownGene, genome = "hg38", chromosome = FTD_controls$CHROM[1], name = gene, transcriptAnnotation = "symbol", size=0.4)
    atrack_controls <- AnnotationTrack(start=FTD_controls$POS, end=FTD_controls$END, chromosome = FTD_controls$CHROM, genome="hg38", name="controls", feature = FTD_controls$ALT, id=FTD_controls$info_control, showFeatureId = TRUE, fontcolor.item="black")
    plotTracks(list(itrack,gtrack, atrack_controls, grtrack), from = start, to = stop, fontcolor.title="#808080",extend.right = 0.2, extend.left = 0.2, fontsize.group=10, fontsize.title=10, BND="#377246", CPX="#86E496", DEL="#D6402D", DUP="#5282AF", INS="#D579E1", INS_ME="#F6AAFF", INS_ME_ALU="#F6AAFF", INS_ME_LINE1="#F6AAFF", INS_ME_SVA="#F6AAFF", INV="#F69A57")
  }
  dev.off()
}

##LBD gene region overlaps
for (i in 1:nrow(genes)) {
  gene = genes$Gene[i]
  start = genes$txStart[i] -1000
  stop = genes$txEnd[i] + 1000
  
  LBD_cases <- LBD %>%
    filter(GENE == gene) %>%
    filter(LBD_ALT_FREQS > 0) %>%
    filter((Candidate_gene_overlap == "Candidate_gene_overlap") | (Candidate_gene_promoter_overlap == "Promoter_overlap"))
  
  LBD_controls <- LBD %>%
    filter(GENE == gene) %>%
    filter(controls_ALT_FREQS > 0) %>%
    filter((Candidate_gene_overlap == "Candidate_gene_overlap") | (Candidate_gene_promoter_overlap == "Promoter_overlap"))
  
  LBD_cases$POS[LBD_cases$POS < start] <- start - (0.2*(stop-start))
  LBD_cases$END[LBD_cases$END > stop] <- stop + (0.2*(stop-start))
  LBD_controls$POS[LBD_controls$POS < start] <- start - (0.2*(stop-start))
  LBD_controls$END[LBD_controls$END > stop] <- stop + (0.2*(stop-start))
  
  png(paste0("NIH_LBD_FTD_SV_shinyapp/NIH_SV_images/LBD/",gene,"_LBD_candidate_gene_overlap.png"), height = 1180, width = 1180, units = "px", res=100)
  
  if ((nrow(LBD_cases) > 0) & (nrow(LBD_controls) > 0)) {
    gtrack <- GenomeAxisTrack()
    itrack <- IdeogramTrack(genome = "hg38", chromosome = LBD_cases$CHROM[1])
    grtrack <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg38.knownGene, genome = "hg38", chromosome = LBD_cases$CHROM[1], name = gene, transcriptAnnotation = "symbol", size=0.4)
    atrack_LBD <- AnnotationTrack(start=LBD_cases$POS, end=LBD_cases$END, chromosome = LBD_cases$CHROM[1], genome="hg38", name="LBD", feature = LBD_cases$ALT, id=LBD_cases$info_case, showFeatureId = TRUE, fontcolor.item="black")
    atrack_controls <- AnnotationTrack(start=LBD_controls$POS, end=LBD_controls$END, chromosome = LBD_controls$CHROM[1], genome="hg38", name="controls", feature = LBD_controls$ALT, id=LBD_controls$info_control, showFeatureId = TRUE, fontcolor.item="black")
    plotTracks(list(itrack,gtrack, atrack_LBD, atrack_controls, grtrack), from = start, to = stop, fontcolor.title="#808080", extend.right = 0.2, extend.left = 0.2, fontsize.group=10, fontsize.title=10, BND="#377246", CPX="#86E496", DEL="#D6402D", DUP="#5282AF", INS="#D579E1", INS_ME="#F6AAFF", INS_ME_ALU="#F6AAFF", INS_ME_LINE1="#F6AAFF", INS_ME_SVA="#F6AAFF", INV="#F69A57")
  } else if ((nrow(LBD_cases) > 0) & (nrow(LBD_controls) == 0)) {
    gtrack <- GenomeAxisTrack()
    itrack <- IdeogramTrack(genome = "hg38", chromosome = LBD_cases$CHROM[1])
    grtrack <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg38.knownGene, genome = "hg38", chromosome = LBD_cases$CHROM[1], name = gene, transcriptAnnotation = "symbol", size=0.4)
    atrack_LBD <- AnnotationTrack(start=LBD_cases$POS, end=LBD_cases$END, chromosome = LBD_cases$CHROM, genome="hg38", name="LBD", feature = LBD_cases$ALT, id=LBD_cases$info_case, showFeatureId = TRUE, fontcolor.item="black")
    plotTracks(list(itrack,gtrack, atrack_LBD, grtrack), from = start, to = stop, fontcolor.title="#808080", extend.right = 0.2, extend.left = 0.2, fontsize.group=10, fontsize.title=10, BND="#377246", CPX="#86E496", DEL="#D6402D", DUP="#5282AF", INS="#D579E1", INS_ME="#F6AAFF", INS_ME_ALU="#F6AAFF", INS_ME_LINE1="#F6AAFF", INS_ME_SVA="#F6AAFF", INV="#F69A57")
  } else if ((nrow(LBD_cases) == 0) & (nrow(LBD_controls) == 0)) {
    gtrack <- GenomeAxisTrack(showTitle=TRUE, name="No SVs", fontcolor.title="#808080", rotation.title=0)
    plotTracks(gtrack, from = 100, to = 110, title.width = 5)
  } else {
    gtrack <- GenomeAxisTrack()
    itrack <- IdeogramTrack(genome = "hg38", chromosome = LBD_controls$CHROM[1])
    grtrack <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg38.knownGene, genome = "hg38", chromosome = LBD_controls$CHROM[1], name = gene, transcriptAnnotation = "symbol", size=0.4)
    atrack_controls <- AnnotationTrack(start=LBD_controls$POS, end=LBD_controls$END, chromosome = LBD_controls$CHROM, genome="hg38", name="controls", feature = LBD_controls$ALT, id=LBD_controls$info_control, showFeatureId = TRUE, fontcolor.item="black")
    plotTracks(list(itrack,gtrack, atrack_controls, grtrack), from = start, to = stop, fontcolor.title="#808080",extend.right = 0.2, extend.left = 0.2, fontsize.group=10, fontsize.title=10, BND="#377246", CPX="#86E496", DEL="#D6402D", DUP="#5282AF", INS="#D579E1", INS_ME="#F6AAFF", INS_ME_ALU="#F6AAFF", INS_ME_LINE1="#F6AAFF", INS_ME_SVA="#F6AAFF", INV="#F69A57")
  }
  dev.off()
}

##FTD gene region overlaps
for (i in 1:nrow(genes)) {
  gene = genes$Gene[i]
  start = genes$txStart[i] -1000
  stop = genes$txEnd[i] + 1000
  
  FTD_cases <- FTD %>%
    filter(GENE == gene) %>%
    filter(FTD_ALT_FREQS > 0) %>%
    filter((Candidate_gene_overlap == "Candidate_gene_overlap") | (Candidate_gene_promoter_overlap == "Promoter_overlap"))
  
  FTD_controls <- FTD %>%
    filter(GENE == gene) %>%
    filter(controls_ALT_FREQS > 0) %>%
    filter((Candidate_gene_overlap == "Candidate_gene_overlap") | (Candidate_gene_promoter_overlap == "Promoter_overlap"))
  
  FTD_cases$POS[FTD_cases$POS < start] <- start - (0.2*(stop-start))
  FTD_cases$END[FTD_cases$END > stop] <- stop + (0.2*(stop-start))
  FTD_controls$POS[FTD_controls$POS < start] <- start - (0.2*(stop-start))
  FTD_controls$END[FTD_controls$END > stop] <- stop + (0.2*(stop-start))
  
  png(paste0("NIH_LBD_FTD_SV_shinyapp/NIH_SV_images/FTD/",gene,"_FTD_candidate_gene_overlap.png"), height = 1180, width = 1180, units = "px", res=100)
  
  if ((nrow(FTD_cases) > 0) & (nrow(FTD_controls) > 0)) {
    gtrack <- GenomeAxisTrack()
    itrack <- IdeogramTrack(genome = "hg38", chromosome = FTD_cases$CHROM[1])
    grtrack <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg38.knownGene, genome = "hg38", chromosome = FTD_cases$CHROM[1], name = gene, transcriptAnnotation = "symbol", size=0.4)
    atrack_FTD <- AnnotationTrack(start=FTD_cases$POS, end=FTD_cases$END, chromosome = FTD_cases$CHROM[1], genome="hg38", name="FTD", feature = FTD_cases$ALT, id=FTD_cases$info_case, showFeatureId = TRUE, fontcolor.item="black")
    atrack_controls <- AnnotationTrack(start=FTD_controls$POS, end=FTD_controls$END, chromosome = FTD_controls$CHROM[1], genome="hg38", name="controls", feature = FTD_controls$ALT, id=FTD_controls$info_control, showFeatureId = TRUE, fontcolor.item="black")
    plotTracks(list(itrack,gtrack, atrack_FTD, atrack_controls, grtrack), from = start, to = stop, fontcolor.title="#808080", extend.right = 0.2, extend.left = 0.2, fontsize.group=10, fontsize.title=10, BND="#377246", CPX="#86E496", DEL="#D6402D", DUP="#5282AF", INS="#D579E1", INS_ME="#F6AAFF", INS_ME_ALU="#F6AAFF", INS_ME_LINE1="#F6AAFF", INS_ME_SVA="#F6AAFF", INV="#F69A57")
  } else if ((nrow(FTD_cases) > 0) & (nrow(FTD_controls) == 0)) {
    gtrack <- GenomeAxisTrack()
    itrack <- IdeogramTrack(genome = "hg38", chromosome = FTD_cases$CHROM[1])
    grtrack <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg38.knownGene, genome = "hg38", chromosome = FTD_cases$CHROM[1], name = gene, transcriptAnnotation = "symbol", size=0.4)
    atrack_FTD <- AnnotationTrack(start=FTD_cases$POS, end=FTD_cases$END, chromosome = FTD_cases$CHROM, genome="hg38", name="FTD", feature = FTD_cases$ALT, id=FTD_cases$info_case, showFeatureId = TRUE, fontcolor.item="black")
    plotTracks(list(itrack,gtrack, atrack_FTD, grtrack), from = start, to = stop, fontcolor.title="#808080", extend.right = 0.2, extend.left = 0.2, fontsize.group=10, fontsize.title=10, BND="#377246", CPX="#86E496", DEL="#D6402D", DUP="#5282AF", INS="#D579E1", INS_ME="#F6AAFF", INS_ME_ALU="#F6AAFF", INS_ME_LINE1="#F6AAFF", INS_ME_SVA="#F6AAFF", INV="#F69A57")
  } else if ((nrow(FTD_cases) == 0) & (nrow(FTD_controls) == 0)) {
    gtrack <- GenomeAxisTrack(showTitle=TRUE, name="No SVs", fontcolor.title="#808080", rotation.title=0)
    plotTracks(gtrack, from = 100, to = 110, title.width = 5)
  } else {
    gtrack <- GenomeAxisTrack()
    itrack <- IdeogramTrack(genome = "hg38", chromosome = FTD_controls$CHROM[1])
    grtrack <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg38.knownGene, genome = "hg38", chromosome = FTD_controls$CHROM[1], name = gene, transcriptAnnotation = "symbol", size=0.4)
    atrack_controls <- AnnotationTrack(start=FTD_controls$POS, end=FTD_controls$END, chromosome = FTD_controls$CHROM, genome="hg38", name="controls", feature = FTD_controls$ALT, id=FTD_controls$info_control, showFeatureId = TRUE, fontcolor.item="black")
    plotTracks(list(itrack,gtrack, atrack_controls, grtrack), from = start, to = stop, fontcolor.title="#808080",extend.right = 0.2, extend.left = 0.2, fontsize.group=10, fontsize.title=10, BND="#377246", CPX="#86E496", DEL="#D6402D", DUP="#5282AF", INS="#D579E1", INS_ME="#F6AAFF", INS_ME_ALU="#F6AAFF", INS_ME_LINE1="#F6AAFF", INS_ME_SVA="#F6AAFF", INV="#F69A57")
  }
  dev.off()
}


