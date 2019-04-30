# ULLL: NEW simplified parsing homoz blocks
# script to call snps for a given region

# !/bin/bash -x
# M Perez-Enciso & J. Leno
# miguel.perez@uab.es
#
# OPTIONS IN BATCH ARE
# -sra2fq : downloads and converts into fastaq, check quality
# -bam2fq
# -bwa    : aligns clean and realigns with GATK
# -vcf    : call snps, filter by max and minimum depth, extracts correct homozygous blocks
# -fasta  : converts into fasta, which can be merged into multiple vcf for later applications
# -qual   : do a few graphs and statistics
# -annot  : uses vep for annotation
# -merge  : merges two realigned bam files
# -clean  : clean files related to sample
#
# STANDARD SOFTWARE REQUIRED
# bwa, samtools (both 1.2.1 and 0.1.19), gatk, picard, bedtools, fastqc, R
#
# CITATION
# If you find useful the pipeline, please cite
#   Bianco et al (2015)

date

DIRDATA=/home/miguel/Projects/Chacuba/data
DIRDATAINRA=/home/miguel/Projects/Chacuba/datainra
DIRBIN=/home/miguel/Projects/Chacuba/bin
DIRBAM=/home/miguel/Projects/Chacuba/bamfiles
DIRVCF=/home/miguel/Projects/Chacuba/varfiles
DIRASSEMBLY=/home/miguel/Projects/Chacuba/assembly
ASSEMBLY=Bos_taurus.UMD3.1.1

bwa=bwa
samtools=$DIRBIN/samtools
GATK=$DIRBIN/GenomeAnalysisTK.jar
picard=$DIRBIN/picard.jar
bedtools=bedtools
bcftools=$DIRBIN/bcftools
fastqdump=fasterq-dump
ASPERA=~/.aspera

##############################################################################################
#                              SAMPLE NAME AND MAIN PARAMETERS                               #
##############################################################################################
OUT=$1
CHR=$2
echo 'sample chr' $OUT $CHR

# global variables
MINCOV=5        # min coverage
MAXCOV=30       # maximum coverage computed later as 2xmean_depth+1
SNPQ=10         # min snp quality
MAPQ=20         # min map quality
BASEQ=20        # min base quality
NP=12           # no. of threads

   cd $DIRVCF
   Q=`awk '$2>0 {a+=$1*$2;n+=$1} END {print a/n}' "$DIRBAM/$OUT/$OUT.realigned.depth" `
   MAXCOV=`echo "tmp=2*($Q+0.5); tmp /= 1; tmp" | bc`
   echo 'depth maxcov' $Q $MAXCOV

   # new in Chacuba
   $bcftools mpileup -r $CHR -Ov -g $MINCOV -q $MAPQ -Q $BASEQ -d $MAXCOV -f $DIRASSEMBLY/$ASSEMBLY.fa \
      $DIRBAM/$OUT/$OUT.realigned.bam | \
      $bcftools call -g$MINCOV -mO z -o $OUT.$CHR.gvcf.gz
   $bcftools filter -O v -g3 -s LOWQUAL -e"%QUAL<$SNPQ || %MAX(INFO/DP)<$MINCOV || %MAX(INFO/DP)>$MAXCOV" \
      $OUT.$CHR.gvcf.gz | \
      awk 'FS="\t" {if ($8 ~ /END/)  $7="PASS"} {OFS="\t"; print $0}' | \
      $bcftools view -f PASS -O z > $OUT.$CHR.flt.gvcf.gz

   # Homozygous blocks are parsed to remove regions with high depth higher than MAXCOV
   # This is done in two steps, first the regions fulfilling the criteria are written into *.block.bed
   # Second, the original gvcf and block.bed are intersected

   # blocks are merged from consecutive positions, a small perl script is used
   $samtools depth -r $CHR -q $BASEQ -Q $MAPQ $DIRBAM/$OUT/$OUT.realigned.bam | \
       	awk -v FS="\t" -v OFS="\t" "(\$3>=$MINCOV && \$3<=$MAXCOV) {print \$1,\$2} " | \
       	perl $DIRBIN/joinConsecutivePos.pl > $OUT.$CHR.blocks.bed

   # gvcf is converted into bed-like, parsed with intersect and back to gvcf,
   # header is saved (no need to parse the whole file) and attached at the end
   zcat $OUT.$CHR.flt.gvcf.gz | head -n500 | grep "#"  > $OUT.header

   zcat $OUT.$CHR.flt.gvcf.gz | grep -v "#" |  \
   awk -v FS="\t" -v OFS="\t" '
      { if($10 ~ "0/0" && $8~"END") {
         split($8,pos,";") 
         split(pos[1],bp,"=")
         print $1,$2-1,bp[2],$3,$4,$5,$6,$7,$8,$9,$10 ;
      } else {
         print $1,$2-1,$2,$3,$4,$5,$6,$7,$8,$9,$10
      } }' |

   $bedtools intersect -a stdin -b $OUT.$CHR.blocks.bed | \

   awk -v FS="\t" -v OFS="\t" '
      {if($11 ~ "0/0" )
         print $1,$2+1,$4,$5,$6,$7,$8,"END="$3,$10,$11
      else
         print $1,$2+1,$4,$5,$6,$7,$8,$9,$10,$11
      }' > $OUT.$CHR.noheader


   # sort option in reverse was to avoid SNP being positioned before the gblock and missed in fasta calling
   #sort -k2,2n -k10,10r $OUT.$CHR.noheader | cat $OUT.header - | bgzip -c > $OUT.$CHR.final.gvcf.gz
   sort -k2,2n -k10,10 $OUT.$CHR.noheader | bgzip -c > $OUT.$CHR.noheader.gvcf.gz

   rm $OUT.$CHR.blocks.bed $OUT.$CHR.noheader $OUT.$CHR.gvcf.gz $OUT.$CHR.flt.gvcf.gz
exit 0

Q=`awk '$2>0 {a+=$1*$2;n+=$1} END {print a/n}' "$DIRBAM/$OUT/$OUT.realigned.depth" `
MAXCOV=`echo "tmp=2*($Q+0.5); tmp /= 1; tmp" | bc`
samtools depth -r $CHR -q $MAPQ -Q $BASEQ $DIRBAM/$OUT/$OUT.realigned.bam | awk '$3 > 5 && $3 < 13' | wc
cat $OUT.$CHR.blocks.bed | awk '{n00+=($3-$2)} END {print n00}'
echo 'raw'
zcat $OUT.$CHR.gvcf.gz | grep '0/0' | grep -v END| wc
zcat $OUT.$CHR.gvcf.gz | grep '0/0' | grep  END | sed s/END=// | sed s/";"/" "/  | awk '{n00+=($8-$2+1)} END {print n00}'
echo 'flt'
zcat $OUT.$CHR.gvcf.gz | grep '0/0' | grep -v END| wc
zcat $OUT.$CHR.flt.gvcf.gz | grep '0/0' | grep  END | sed s/END=// | sed s/";"/" "/  | awk '{n00+=($8-$2+1)} END {print n00}'

zcat $OUT.final.gvcf.gz | grep '0/0' | sed s/"END="// | sed s/";"/" "/  | awk '{n00+=($8-$2+1)} END {print n00}'
### test


zcat $OUT.$CHR.flt.gvcf.gz | grep '0/0' |  \
   awk -v FS="\t" -v OFS="\t" '
      { if($10 ~ "0/0" && $8~"END") {
         split($8,pos,";") 
         split(pos[1],bp,"=")
         print $1,$2-1,bp[2],$3,$4,$5,$6,$7,$8,$9,$10 ;
      } else {
         print $1,$2-1,$2,$3,$4,$5,$6,$7,$8,$9,$10
      } }' | awk '{n00+=($3-$2)} END {print n00}'
 


zcat $OUT.$CHR.flt.gvcf.gz | grep '0/0' |  \
   awk -v FS="\t" -v OFS="\t" '
      { if($10 ~ "0/0" && $8~"END") {
         split($8,pos,";") 
         split(pos[1],bp,"=")
         print $1,$2-1,bp[2],$3,$4,$5,$6,$7,$8,$9,$10 ;
      } else {
         print $1,$2-1,$2,$3,$4,$5,$6,$7,$8,$9,$10
      } }' |

   $bedtools intersect -a stdin -b $OUT.$CHR.blocks.bed |

   awk -v FS="\t" -v OFS="\t" '
      {if($11 ~ "0/0" && $9 !~ "END")
         print $1,$2+1,$4,$5,$6,$7,$8,"END="$3";"$9,$10,$11
      else
         print $1,$2+1,$4,$5,$6,$7,$8,$9,$10,$11
      }' | sed s/"END="// | sed s/";"/" "/  | awk '{n00+=($8-$2+1)} END {print n00}'

