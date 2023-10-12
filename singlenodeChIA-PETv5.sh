echo `bash --version`
module load macs/041014
module load bedtools/2.30.0
module load samtools/1.15.1
module load bowtie/1.2.2
module load R/4.2.2
module load conda3/202210
conda activate cutadapt_global
## parameters Queue
## cutadapt (v 4.3)
## scripts plotpositions.R interaction_pvalue.R


if [ "$1" == "-h" ]; then
  echo "usage: script.sh -f file_in.fastq -r file2.fastq -b bowtieindexname"
  exit 0
fi

ARGS=$(getopt -a --options f:r:b: --long "fastq1:,fastq2:,bowtiepath:" -- "$@")

eval set -- "$ARGS"

while true; do
  case "$1" in
    -f|--fastq1)
      fastq1="$2"
      shift 2;;
    -r|--fastq2)
      fastq2="$2"
      shift 2;;
    -b|--bowtiepath)
      bowtiepath="$2"
      shift 2;;
    --)
      break;;
  esac
done
echo "$fastq1"
echo "$fastq2"
echo "$bowtiepath"

#===========#
#format ChIAPET fastq
# 1.3 /home/avelucha/.local/bin/cutadapt
#4.3 /home/avelucha/miniconda3/envs/cutadapt-env/bin/cutadapt
#===========#

#===========#
#find the distribution of linker in reads
#===========#
cutadapt -j 40 --info-file positions.tsv --no-indels --pair-adapters --pair-filter any -A "CTGCTGTCATN{13}X" -a "CTGCTGTCATN{13}X"  -m 27:27 -M 36:36 -n 1 -O 10 -e 0 -o R1trimmed_1.fastq -p R1trimmed_2.fastq --untrimmed-output R1untrimmed.fastq --untrimmed-paired-output R1untrimmed.paired.fastq --too-short-output R1short_1.fastq --too-short-paired-output R1short_2.fastq --too-long-output R1long_1.fastq --too-long-paired-output R1long_2.fastq $fastq1 $fastq2

cutadapt -j 40 --no-indels --pair-adapters --pair-filter any -A "CTGCTGTCCGN{13}X" -a "CTGCTGTCCGN{13}X"  -m 27:27 -M 36:36 -n 1 -O 10 -e 0 -o R2trimmed_1.fastq -p R2trimmed_2.fastq --untrimmed-output R2untrimmed.fastq --untrimmed-paired-output R2untrimmed.paired.fastq --too-short-output R2short_1.fastq --too-short-paired-output R2short_2.fastq --too-long-output R2long_1.fastq --too-long-paired-output R2long_2.fastq $fastq1 $fastq2

cat R1trimmed_1.fastq R2trimmed_1.fastq > R1_formatted_1.fastq
cat R1trimmed_2.fastq R2trimmed_2.fastq > R2_formatted_2.fastq

#===========#
# plot linker positions in in reads
#===========#

cut -f3 readposition.tsv | grep -v -P "A|T|G|C" > positions.tsv
R CMD BATCH plotpositions.R 

#===========#
#read alignment for non-chimeric reads
# mine: bowtie index: /research_jude/rgs01_jude/groups/abrahgrp/projects/3D_GENOME_CONSORTIUM/abrahgrp/Raj/INS_mESC/mm10Bowtie2/mm10
#===========#

bowtie -e 70 -k 1 -m 1 -v 1 -p 4 --best --strata -S $bowtiepath R1_formatted_1.fastq > R1_formatted.SAM
bowtie -e 70 -k 1 -m 1 -v 1 -p 4 --best --strata -S $bowtiepath R2_formatted_2.fastq > R2_formatted.SAM


##===========#
## Alignment formatting
##===========#

samtools view -F 4 -bS 'R1_formatted.SAM' > R1_formatted.bam
samtools view -F 4 -bS 'R2_formatted.SAM' > R2_formatted.bam

##===========#
## Merging Alignment and fixing the mate
##===========#

samtools merge -@ 16 R12merged.bam R1_formatted.bam R2_formatted.bam
samtools sort -@ 16 -n R12merged.bam > R12merged.sorted.bam
samtools fixmate R12merged.sorted.bam R12merged.fixedmate
bedtools bamtobed -i R12merged.fixedmate -bedpe > R12merged.fixedmate.bedpe

##===========#
## Prepare fixed mate/PET for peak finding
##===========#
awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$9,$10}' R12merged.fixedmate.bedpe | sort | uniq  > R12_unique.bedpe
awk 'BEGIN{FS=OFS="\t"}{if(($1==$4) && ($5>=$3)){print $0}}' R12_unique.bedpe > R12_unique_preanchor.bedpe
awk 'BEGIN{FS=OFS="\t"}{if($1!=$4){print $0}}' R12_unique.bedpe >> R12_unique_preanchor.bedpe
awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$1$2$3$4$5$6,"255",$7"\n"$4,$5,$6,$1$2$3$6$5$6,"255",$8}' R12_unique_preanchor.bedpe > R12_unique_anchor.bed

##===========#
## peak finding with macs14
##===========#
macs14 --nomodel --nolambda -g mm -t R12_unique_anchor.bed -f BED -n R12_unique_anchor_1e-9_nolambda -p 1e-9 -m 10,30 -w -S --space=50 --keep-dup=2
awk 'BEGIN{FS="	";OFS="	"} {print $1,$2,$3,$4,$5,$1":"$2"-"$3}' R12_unique_anchor_1e-9_nolambda_peaks.bed >| R12_peaks.bed 

##===========#
# Filter 5kb and Pooling intra and inter
##===========#

awk 'BEGIN{FS=OFS="\t"}{if(($1==$4) && ($5 - $3 > 4999)){print $1,$2,$3,$4,$5,$6,$7$8}}' R12_unique.bedpe > R12_unique.intra5000andinter
awk 'BEGIN{FS=OFS="\t"}{if($1!=$4){print $1,$2,$3,$4,$5,$6,$7$8}}' R12_unique.bedpe >> R12_unique.intra5000andinter

##===========#
## 1. Overlap with peaks: Peaks PET succinct done
##===========#
pairToBed -a R12_unique.intra5000andinter -b R12_peaks.bed -type both -f 0.75  | awk -F"\t" '{print $8"\t"$9"\t"$10}' | sort | uniq -c | sed 's/^\s*//' | sed 's/\s/\t/' | awk -F"\t" '{print $2"\t"$3"\t"$4"\t"$2":"$3"-"$4"\t"$1}' > R12_peaks_PETsummary_succinct.txt

##===========#
##2. Creating interaction succinct
##===========#

pairToBed -a R12_unique.intra5000andinter -b R12_peaks.bed -type both -f 0.75  > R12_peaks_formatted.bed
awk -F"\t" '{print $1"-"$2"-"$3"-"$4"-"$5"-"$6"\t"$NF}' R12_peaks_formatted.bed > R12_peaks_formatted_rearr.bed
awk -F"\t" '$1==last {printf "\t%s",$2; next} NR>1 {print "";} {last=$1; printf "%s",$0;} END{print "";}' R12_peaks_formatted_rearr.bed | awk -F"\t" '{print $2"\t"$3}' | sed 's/\:/\t/g' | sed 's/\-/\t/g' | sort | uniq > R12_peaks_formatted_rearr.sort.bed

pairToBed -a R12_unique.intra5000andinter -b R12_peaks.bed -type both -f 0.75 | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7}' | sort | uniq > R12_unique.intra5000andinter_f0.75.bedpe

pairToPair -a R12_peaks_formatted_rearr.sort.bed -b R12_unique.intra5000andinter_f0.75.bedpe -type both | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6}' | sort | uniq -c | sed 's/^\s*//' | sed 's/\s/\t/' | awk -F"\t" '{print $2":"$3"-"$4"\t"$5":"$6"-"$7"\t"$2":"$3"-"$4"=="$5":"$6"-"$7"\t"$1"\t"$6-$4}' > R12_PET_interactionSummary_succinct.txt


##===========#
## Find statistically significant interaction
##===========#
R --no-save R12_PET_interactionSummary_succinct.txt R12_peaks_PETsummary_succinct.txt 0.01 2 R12_filtered_PET < interaction_pvalue.R

#### End of execution ####


##===========#
## calling Insulated neighborhood
##===========#

#pairToBed -a R12_filtered_PET_interactionSummary_n2FDR0.01.bedpe -b CTCF_ChIP_peaks.bed -type both | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7}' | sort | uniq > R12_SignificantInteractionVsCTCF.bedpe
#pairToBed -a R12_SignificantInteractionVsCTCF.bedpe -b Annotation_genes.bed -type ospan | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7}' | sort | uniq > R12_INS.bedpe
