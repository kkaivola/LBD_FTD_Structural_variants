#!/bin/bash

#Extract each LBD sample of interest separately from 1% GATK-SV (all variants), include only SVs with an alternative allele (=exclude 0/0s) 

SR_SAMPLES="RES01827 RES01897 RES02284 UMARY-5077 RES08227 UMARY-1571 UMARY-634 RES01859 RES01877 UMARY-1544 RES00032 RES00064 RES00275 UMARY-4542 RES00071 RES00048"

for SAMPLE in $SR_SAMPLES; do
    bcftools view \
    --threads $SLURM_CPUS_PER_TASK \
    -s $SAMPLE \
    -Ou \
    /data/ALS_50k/DementiaSeq_SV.LBD.FTD/CLEAN.SV.fromJinHui/LBD_gatksv_1perc_fdr_cleaned_filters_qual_recalibrated/LBD_1perc_fdr.cleaned_filters_qual_recalibrated.LNGsampleID.vcf.gz \
    | \
    bcftools view -i 'GT[*]="alt"' \
    --threads $SLURM_CPUS_PER_TASK \
    -Oz \
    -o ${SAMPLE}_altgtonly_LBD_1perc_fdr.cleaned_filters_qual_recalibrated.LNGsampleID.vcf.gz 
done

#Index vcfs
parallel -j1 tabix -p vcf {} ::: *_altgtonly_LBD_1perc_fdr.cleaned_filters_qual_recalibrated.LNGsampleID.vcf.gz

#Split by SVTYPE into separate vcfs
module load samtools

cd /data/ALS_50k/karri/LBD_FTD_SV_project/nanopore_vs_sr_validation/SR_vcfs
mkdir split_SR_vcfs

SR_SAMPLES="RES01827 RES01897 RES02284 UMARY-5077 RES08227 UMARY-1571 UMARY-634 RES01859 RES01877 UMARY-1544 RES00032 RES00064 RES00275 UMARY-4542 RES00071 RES00048"

for SAMPLE in $SR_SAMPLES; do
    for svtype in "DEL" "INS" "DUP" "INV" "CPX"; do 
        filt="SVTYPE=\"$svtype\""
        bcftools view -i $filt \
        --threads $SLURM_CPUS_PER_TASK \
        -O z \
        -o split_SR_vcfs/${SAMPLE}_${svtype}_1percfdr_pass_multialleic_alt_SVs.vcf.gz \
        ${SAMPLE}_altgtonly_LBD_1perc_fdr.cleaned_filters_qual_recalibrated.LNGsampleID.vcf.gz 
    done
done


#Rename DUP SVTYPE to INS since that is what sniffles2 for Nanopore (mostly) uses. But keep DUP file separate from INS to be able to differ DUP from INS

SR_SAMPLES="RES01827 RES01897 RES02284 UMARY-5077 RES08227 UMARY-1571 UMARY-634 RES01859 RES01877 UMARY-1544 RES00032 RES00064 RES00275 UMARY-4542 RES00071 RES00048"

for SAMPLE in $SR_SAMPLES; do
    zcat ${SAMPLE}_DUP_1percfdr_pass_multialleic_alt_SVs.vcf.gz | \
    sed 's/SVTYPE=DUP/SVTYPE=INS/g' > ${SAMPLE}_DUP_as_INS_1percfdr_pass_multialleic_alt_SVs.vcf
    bgzip ${SAMPLE}_DUP_as_INS_1percfdr_pass_multialleic_alt_SVs.vcf
done



#Index vcfs 
parallel -j $SLURM_CPUS_PER_TASK tabix -f -p vcf {} ::: *vcf.gz


#To match Nanopore format, for INS short read vcfs, make %END = %POS
module load bcftools
cd /data/ALS_50k/karri/LBD_FTD_SV_project/nanopore_vs_sr_validation/SR_vcfs/split_SR_vcfs

##Extract needed info
SR_SAMPLES="RES01827 RES01897 RES02284 UMARY-5077 RES08227 UMARY-1571 UMARY-634 RES01859 RES01877 UMARY-1544 RES00032 RES00064 RES00275 UMARY-4542 RES00071 RES00048"

for SAMPLE in $SR_SAMPLES; do
    bcftools query -f '%CHROM\t%POS\t%END\t%POS\n' ${SAMPLE}_INS_1percfdr_pass_multialleic_alt_SVs.vcf.gz >  ${SAMPLE}_end_annotation.txt
    bgzip -f ${SAMPLE}_end_annotation.txt
    tabix -f -s1 -b2 -e3 ${SAMPLE}_end_annotation.txt.gz
done

##Update END flag: but make it STOP named since END cannot be directly altered.
echo "##INFO=<ID=STOP,Number=1,Type=Integer,Description=""End position of the structural variant"">"  > tmp_stop_header.txt

##Add STOP flag and remove END flag
SR_SAMPLES="RES01827 RES01897 RES02284 UMARY-5077 RES08227 UMARY-1571 UMARY-634 RES01859 RES01877 UMARY-1544 RES00032 RES00064 RES00275 UMARY-4542 RES00071 RES00048"

for SAMPLE in $SR_SAMPLES; do
    bcftools annotate -a ${SAMPLE}_end_annotation.txt.gz -c CHROM,FROM,TO,STOP -x INFO/END -h tmp_stop_header.txt -Oz -o ${SAMPLE}_tmp.vcf.gz ${SAMPLE}_INS_1percfdr_pass_multialleic_alt_SVs.vcf.gz
    zcat ${SAMPLE}_tmp.vcf.gz | sed 's/STOP/END/g' > ${SAMPLE}_INS_end_is_pos_1percfdr_pass_multialleic_alt_SVs.vcf
    bgzip ${SAMPLE}_INS_end_is_pos_1percfdr_pass_multialleic_alt_SVs.vcf
    tabix -p vcf ${SAMPLE}_INS_end_is_pos_1percfdr_pass_multialleic_alt_SVs.vcf.gz
done



#Split LBD Nanopore vcfs by SV type

SAMPLES_NANOPORE="UMARY634 UMARY1544 UMARY4542 UMARY5077 UMARY1571 UMIC1205 UMIC1527 UMIC1601 UMIC959 BARCS1192 BARCS1429 BARCS1720 BARCS1777 INDU2004038 OREG2142 OREG2595"

for SAMPLE in $SAMPLES_NANOPORE; do
    for svtype in "DEL" "INS" "DUP" "INV"; do 
        filt="SVTYPE=\"$svtype\""
        bcftools view -i $filt \
        --threads $SLURM_CPUS_PER_TASK \
        -O z \
        -o split_vcfs/${SAMPLE}_${svtype}_fastq_pass.wf_sv.vcf.gz \
        ${SAMPLE}_fastq_pass.wf_sv.vcf.gz 
    done
done


#Index split Nanopore vcfs by SV type
parallel -j $SLURM_CPUS_PER_TASK tabix -p vcf {} ::: *vcf.gz


#Run truvari bench per SV type per sample
##Validation criterion are: refdist:200bp, pctsize=0.70, pctsim=0,typeignore=FALSE
##Nanopore sample data in same order than short read sample data
SR_SAMPLES="RES01827 RES01897 RES02284 UMARY-5077 RES08227 UMARY-1571 UMARY-634 RES01859 RES01877 UMARY-1544 RES00032 RES00064 RES00275 UMARY-4542 RES00071 RES00048"
SAMPLES_NANOPORE="BARCS1720 BARCS1192 INDU2004038 UMARY5077 UMIC959 UMARY1571 UMARY634 BARCS1777 BARCS1429 UMARY1544 OREG2142 UMIC1205 UMIC1527 UMARY4542 UMIC1601 OREG2595"


##DEL,INV
i=0
for SAMPLE in $SR_SAMPLES; do
    array=(BARCS1720 BARCS1192 INDU2004038 UMARY5077 UMIC959 UMARY1571 UMARY634 BARCS1777 BARCS1429 UMARY1544 OREG2142 UMIC1205 UMIC1527 UMARY4542 UMIC1601 OREG2595)
    for svtype in DEL INV; do
        truvari bench \
        --refdist=200 \
        --pctsize=0.70 \
        --pctsim=0 \
        --sizemax=999999999 \
        -b vcfs_sniffles2/split_vcfs/${array[${i}]}_${svtype}_fastq_pass.wf_sv.vcf.gz \
        -c SR_vcfs/split_SR_vcfs/${SAMPLE}_${svtype}_1percfdr_pass_multialleic_alt_SVs.vcf.gz \
        -f /data/kaivolakk/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta \
        -o truvari_bench/truvari_bench_${SAMPLE}_${svtype} 
    done
    i=${i}+1
done

##INS, separately since SR path name is different
i=0
for SAMPLE in $SR_SAMPLES; do
    array=($SAMPLES_NANOPORE)
    for svtype in INS; do
        truvari bench \
        --refdist=200 \
        --pctsize=0.70 \
        --pctsim=0 \
        --sizemax=999999999 \
        -b vcfs_sniffles2/split_vcfs/${array[${i}]}_${svtype}_fastq_pass.wf_sv.vcf.gz \
        -c SR_vcfs/split_SR_vcfs/${SAMPLE}_INS_end_is_pos_1percfdr_pass_multialleic_alt_SVs.vcf.gz \
        -f /data/kaivolakk/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta \
        -o truvari_bench/truvari_bench_${SAMPLE}_${svtype} 
    done
    i=${i}+1
done


##Run SR DUPs as INS
i=0
for SAMPLE in $SR_SAMPLES; do
    array=($SAMPLES_NANOPORE)
    for svtype in DUP; do
        truvari bench \
        --refdist=200 \
        --pctsize=0.70 \
        --pctsim=0 \
        --sizemax=999999999 \
        -b vcfs_sniffles2/split_vcfs/${array[${i}]}_INS_fastq_pass.wf_sv.vcf.gz \
        -c SR_vcfs/split_SR_vcfs/${SAMPLE}_DUP_as_INS_1percfdr_pass_multialleic_alt_SVs.vcf.gz \
        -f /data/kaivolakk/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta \
        -o truvari_bench/truvari_bench_${SAMPLE}_${svtype} 
    done
    i=${i}+1
done

conda deactivate
