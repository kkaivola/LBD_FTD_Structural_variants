# Background
BigBed files were created according to instructions by UCSC at https://genome.ucsc.edu/goldenPath/help/bigBed.html .  
Color codes follow UCSC color coding https://genome.ucsc.edu/goldenPath/help/hgCnvColoring.html .

_Input files_: gzipped and tabix/indexed vcfs with clean samples and high-quality subset of structural variants, split by case-control status
_Output files_: One bigBed file per phenotype-status
```
##Fetch chrom sizes (UCSC binary)
fetchChromSizes hg38 > hg38.chrom.sizes

##Create bigbeds
for PHENO in LBD FTD; do
  for STATUS in cases controls; do
        bcftools query -f '%CHROM\t%POS\t%END\t"INFO"\t%QUAL\t"."\t%POS\t%END\t"RGB"\t%AF\t%SVLEN\t%SVTYPE\t%ID\n' /data/ALS_50k/karri/LBD_FTD_SV_project/filtered_SV_GATKSV_style_barplots/${PHENO}_${STATUS}_filtered_SV_Samples_GQ300_to_missing.vcf.gz > tmp_${PHENO}_${STATUS}_analyzed.bed #Extract info per variant
        sed -i 's/"//g' tmp_${PHENO}_${STATUS}_analyzed.bed #Remove quote marks
        awk '($10 > 0)' tmp_${PHENO}_${STATUS}_analyzed.bed > tmp.txt #Filter out variants with AF = 0
        awk '{if($12=="INS") {$3=$2 ; $8=$2; print} else{print $0}}' tmp.txt > tmp2.txt #If SVTYPE is INS, make END=POS
        awk '{if($12=="INS") {$11=1; print} else{print $0}}' tmp2.txt > tmp3.txt #For INS, set SVLEN = 1
        awk '{if($11<0) {$11=1; print} else{print $0}}' tmp3.txt > tmp4.txt #For those SVs (BNDs, CTX) where SVLEN is under -1 (undetermined), set SVLEN to 1
        awk '{if($12 == "DUP") {$9="0,0,255"; print} else{print $0}}' tmp4.txt > tmp5.txt #Add color code
        awk '{if($12 == "DEL") {$9="255,0,0"; print} else{print $0}}' tmp5.txt > tmp6.txt #Add color code
        awk '{if($12 == "INS") {$9="255,153,51"; print} else{print $0}}' tmp6.txt > tmp7.txt #Add color code
        awk '{if($12 == "INV") {$9="102,0,204"; print} else{print $0}}' tmp7.txt > tmp8.txt #Add color code
        awk '{if($12 == "INV") {$9="102,0,204"; print} else{print $0}}' tmp8.txt > tmp9.txt #Add color code
        awk '{if($9 == "RGB") {$9="128,128,128"; print} else{print $0}}' tmp9.txt > tmp10.txt #Add color code
        awk '{printf "%.2f%\n",$10*100}' tmp10.txt > tmp11.txt #Make INFO field as "SVTYPE_AF(rounded % two decimals)_QUAL"
        paste tmp10.txt tmp11.txt > tmp12.txt
        awk '$4=$12"_AF="$14"_QUAL="$5' tmp12.txt > tmp13.txt
        awk '$4=$12"_AF="$14"_QUAL="$5' tmp12.txt > tmp13.txt
        awk '$2=$2-1' tmp13.txt > tmp14.txt #Make 0-based bed by substracting 1 from POS
        cut -d " " -f1-13 tmp14.txt > ${PHENO}_${STATUS}_analyzed.bed
        sort -k1,1 -k2,2n ${PHENO}_${STATUS}_analyzed.bed -o ${PHENO}_${STATUS}_analyzed.bed #Sort bed file
        bedToBigBed  -as=fieldnames.as -type=bed9+4 ${PHENO}_${STATUS}_analyzed.bed hg38.chrom.sizes ${PHENO}_${STATUS}_analyzed.bb #Make bed into bigbed file
    done
done
    
```

Autosql file fieldnames.as

```
table X
"Clean samples with high_quality_SVs"
(
string  chrom;		"Reference sequence chromosome or scaffold"
uint    chromStart;	"Start position of feature on chromosome"
uint    chromEnd;	"End position of feature on chromosome"
string  name;		"SVTYPE, AF in % and QUAL"
uint    score;		"QUAL score"
char[1] strand;		"no strand info"
uint    thickStart;	"Coding region start"
uint    thickEnd;	"Coding region end"
uint  	reserved;	"Green on + strand, Red on - strand"
float  AF;  "Allele frequency"
uint    SVLEN;  "Structural variant lenght"
string  SVTYPE;		"Structural variant type"
string	ID;	"Structural variant ID"
)
```
