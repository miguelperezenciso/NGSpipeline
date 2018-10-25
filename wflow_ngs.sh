# !/bin/bash -x
# NGS generic pipeline to 
# M Perez-Enciso with J. Leno-Colorado
# miguel.perez@uab.es
#
# OPTIONS IN BATCH ARE
# -sra2fq  : downloads and converts into fastq, check quality
# -bam2fq  : extract reads from bam
# -bwa     : aligns clean and realigns with GATK and picardtools
# -vcf     : call snps, filter by max and minimum depth, extracts correct homozygous blocks
# -fasta   : converts into fasta, which can be merged into multiple vcf for later applications
# -qual    : do a few statistics
# -pdf     : do a plot with depths and variability
# -clean   : clean files related to sample
# -fas2vcf : combines several fasta files into a single vcf file
#
# STANDARD SOFTWARE REQUIRED
# fasterq-dump, bwa, samtools, bcftools (1.8), bgzip, gatk, picard, bedtools, R
#
# NON-STANDARD SOFTWARE REQUIRED (depending on option) they should be in $DIRBIN
# joinConsecutivePos.pl
# covXwin.pl
# ngs_theta.f90 & fact-m.f90; to compile f95 -O4 fact-m.f90 ngs_theta.f90 -o ngs_theta
# plotNGS.R
# fas2vcf.f90; to compile f95 -O4 fas2vcf.f90 -o fas2vcf
#
# CITATION
# If you find useful the pipeline, please cite
#   Bianco et al (2015)
# If you use the estimations of variability (-qual), please cite
#   Esteve-Codina et al. (2011) Partial short-read sequencing of a highly inbred Iberian pig and genomics inference thereof.
#   Heredity, 107:256-264. doi: 10.1038/hdy.2011.13


# to transform fasta into vcf
# $DIRBIN/fas2vcf -i falist.txt -f $DIRASSEMBLY/$ASSEMBLY.fa | sed s/' '//g | $bgzip -c > file.vcf.gz

date

##############################################################################################
#            FOLDER STRUCTURE    modify to suit your taste                                   #
##############################################################################################
WDIR=/home/miguel/NGS    # working directory
DIRDATA=$WDIR/data       # contains reads
DIRBIN=$WDIR/bin         # contains executables
DIRBAM=$WDIR/bamfiles    # contains bamfiles and depth files, one folder per sample
DIRVCF=$WDIR/varfiles    # contains vcf and other output files
DIRFAS=$WDIR/fasfiles    # contains fasta files
DIRASSEMBLY=$WDIR/assembly # contains asssembly and indices

# actual assembly name should be $ASSEMBLY.fa
ASSEMBLY=Bos_taurus.UMD3.1.1

##### links to software required
bwa=bwa
samtools=$DIRBIN/samtools
GATK=$DIRBIN/GenomeAnalysisTK.jar
picard=$DIRBIN/picard.jar
bedtools=bedtools
bcftools=$DIRBIN/bcftools
bgzip=bgzip
fastqdump=fasterq-dump
ASPERA=~/.aspera

##############################################################################################
#                              SAMPLE NAME AND MAIN PARAMETERS                               #
##############################################################################################
OPTION=$1
SRA=$2        # sample analyzed
OUT=$SRA
echo 'samplex' $OUT

# reads from sample, they are in $DIRDATA
# WARNING: assumed that are uncompressed, as output from fasterq-dump
READS_PE1=${SRA}.sra_1.fastq
READS_PE2=${SRA}.sra_2.fastq

# global variables
MINCOV=5        # min coverage
MAXCOV=         # maximum coverage computed later as 2xmean_depth+1
SNPQ=10         # min snp quality
MAPQ=20         # min map quality
BASEQ=20        # min base quality
NP=10           # no. of threads
WINSIZE=100000  # window size used to comple plots (-pdf option)


################################################################################
# 0. INDEXING THE GENOME Build index sequence archive for reference sequence
################################################################################
if [ $OPTION = "-index" ]; then
   cd $DIRASSEMBLY
   time $bwa index -a bwtsw $DIRASSEMBLY/$ASSEMBLY.fa

   # creating index with samtools
   time samtools faidx $ASSEMBLY.fa

   # dictionary with picard tools' CreateSequenceDictionary (same name -> dict=reference)
   time java -jar $picard CreateSequenceDictionary R=$DIRASSEMBLY/$ASSEMBLY.fa O=$DIRASSEMBLY/$ASSEMBLY.dict
fi


################################################################################
# 1a. CONVERT SRA TO FASTQ (allows any number of lanes / sample)
################################################################################
if [ $OPTION = "-sra2fq" ]; then
    cd $DIRDATA
    SRA1=$OUT
    for SRA in "$@"; do
        if [ $SRA != "-sra2fq" ]; then
           echo $SRA
           DIRSRA=/sra/sra-instant/reads/ByRun/sra/${SRA:0:3}/${SRA:0:6}/$SRA
           $ASPERA/connect/bin/ascp -i  $ASPERA/connect/etc/asperaweb_id_dsa.openssh \
                                    -k1 -Tr -l100m anonftp@ftp-private.ncbi.nlm.nih.gov:$DIRSRA $DIRDATA
           cd $SRA
           $fastqdump -e $NP --split-files $SRA.sra -O $DIRDATA/$SRA1
           rm $SRA.sra
           cd ..
        fi 
    done
    cd $DIRDATA
    cat $SRA1/SRR*.sra_1.fastq > $SRA1.sra_1.fastq
    cat $SRA1/SRR*.sra_2.fastq > $SRA1.sra_2.fastq
    echo 'No. reads ' $SRA1
    wc -l $SRA1.sra_2.fastq
    rm $SRA1/*
fi

################################################################################
# 1b. Extract paired reads from bam file
################################################################################
if [ $OPTION = "-bam2fq" ]; then
    cd $DIRDATA
    time java -jar $picard SamToFastq I=$OUT.bam FASTQ=$OUT.sra_1.fastq SECOND_END_FASTQ=$OUT.sra_2.fastq
    mv $OUT.sra_?.fastq $DIRDATA/.
fi

################################################################################
# 2. BWA ALIGNMENT AND GATK REALIGNMENT
################################################################################
if [ $OPTION = "-bwa" ]; then
   mkdir $DIRBAM/$OUT
   cd $DIRBAM/$OUT
   TAG="@RG\tID:$OUT\tSM:$OUT"

   # align and sort
   time $bwa mem -t $NP -R $TAG $DIRASSEMBLY/$ASSEMBLY.fa $DIRDATA/$READS_PE1 $DIRDATA/$READS_PE2 | \
        $samtools view -b -  > $OUT.bam
   #time $bwa mem -t $NP -R $TAG $DIRASSEMBLY/$ASSEMBLY.fa $DIRDATA/$READS_PE1 $DIRDATA/$READS_PE2 $OUT.raw.sam
   #$samtools fixmate -O bam $OUT.raw.sam $OUT.bam
   time $samtools sort -m 4G  -@ $NP -O bam -T tmp $OUT.bam > $OUT.sort.bam

   # index
   $samtools index $OUT.sort.bam
   echo $OUT 'bam index finished'

   # BAM REFINEMENT WITH GATK
   # realign reads around indels with GATK
   time java -Xmx4G -jar $GATK -T RealignerTargetCreator -R $DIRASSEMBLY/$ASSEMBLY.fa \
                  -I $OUT.sort.bam -o $OUT.intervals -nt $NP --allow_potentially_misencoded_quality_scores

   time java -Xmx4G -jar $GATK \
	          -T IndelRealigner \
	          -R $DIRASSEMBLY/$ASSEMBLY.fa \
	          -I $OUT.sort.bam \
	          -targetIntervals $OUT.intervals \
	          -o $OUT.rl.bam --allow_potentially_misencoded_quality_scores
   echo $OUT 'bam realigned finished'

   # rm duplicates with picard
   time java -Xmx4g -jar $picard MarkDuplicates \
                         REMOVE_DUPLICATES=true \
                         INPUT=$OUT.rl.bam \
                         OUTPUT=$OUT.realigned.bam \
                         METRICS_FILE=metrics.out
   echo $OUT 'bam rmdup finished'

   # depth file
   time $samtools depth -q $BASEQ -Q $MAPQ $OUT.realigned.bam | awk '{print $3}'  | \
                 sort | uniq -c | sort -n -k2 > $OUT.realigned.depth

   $samtools index $OUT.realigned.bam
   echo $OUT 'bam index finished'
fi


################################################################################
# 4a. SNP CALLING with gVCF blocks for each chr separately                     #
################################################################################
# requires joinConsecutivePos.pl
if [ $OPTION = "-cvcf" ]; then
   cd $DIRVCF
   $samtools view -H $DIRBAM/$OUT/$OUT.realigned.bam | grep "@SQ"| awk '{print $2}' | sed  s/SN:// > chr.names
   for chr in $(cat chr.names); do
       time sh ../wflow_ngs_vcf_chr.sh $OUT $chr &
       # consider editing this if too many chrs simultaneously
       sleep 1s
   done
done

if [ $OPTION = "-cmerge" ]; then
   # header must be generated by -cvcf option
   bgzip $OUT.header 
   mv $OUT.header.gz $OUT.final.gvcf.gz
   for chr in $(cat chr.names); do
       cat $OUT.$chr.noheader.gvcf.gz >> $OUT.final.gvcf.gz
       rm $OUT.$chr.noheader.gvcf.gz
   done
fi

################################################################################
# 4b. SNP CALLING with gVCF blocks all chrs                                    #
################################################################################
if [ $OPTION = "-vcf" ]; then
   cd $DIRVCF
   Q=`awk '$2>0 {a+=$1*$2;n+=$1} END {print a/n}' "$DIRBAM/$OUT/$OUT.realigned.depth" `
   MAXCOV=`echo "tmp=2*($Q+0.5); tmp /= 1; tmp" | bc`
   echo 'sample meanDepth maxDepth ' $OUT $Q $MAXCOV

   # samtools 1.8
   $bcftools mpileup -Ov -g $MINCOV -q $MAPQ -Q $BASEQ -d $MAXCOV -f $DIRASSEMBLY/$ASSEMBLY.fa \
      $DIRBAM/$OUT/$OUT.realigned.bam | \
      $bcftools call -g$MINCOV -mOz -o $OUT.gvcf.gz
   $bcftools filter -O v -g3 -s LOWQUAL -e"%QUAL<$SNPQ || %MAX(INFO/DP)<$MINCOV || %MAX(INFO/DP)>$MAXCOV" \
      $OUT.gvcf.gz | \
      awk 'FS="\t" {if ($8 ~ /END/)  $7="PASS"} {OFS="\t"; print $0}' | \
      $bcftools view -f PASS -O z > $OUT.flt.gvcf.gz

   # Homozygous blocks are parsed to remove regions with high depth higher than MAXCOV
   # This is done in two steps, first the regions fulfilling the criteria are written into *.block.bed
   # Second, the original gvcf and block.bed are intersected

   # blocks are merged from consecutive positions, a small perl script is used
   $samtools depth -q $BASEQ -Q $MAPQ $DIRBAM/$OUT/$OUT.realigned.bam | \
       	awk -v FS="\t" -v OFS="\t" "(\$3>=$MINCOV && \$3<=$MAXCOV) {print \$1,\$2} " | \
       	perl $DIRBIN/joinConsecutivePos.pl > $OUT.blocks.bed

   # gvcf is converted into bed-like, parsed with intersect and back to gvcf,
   # header is saved (no need to parse the whole file) and attached at the end
   zcat $OUT.flt.gvcf.gz | head -n500 | grep "#"  > $OUT.header

   zcat $OUT.flt.gvcf.gz | grep -v "#" |  \
   awk -v FS="\t" -v OFS="\t" '
      { if($10 ~ "0/0" && $8~"END") {
         split($8,pos,";") 
         split(pos[1],bp,"=")
         print $1,$2-1,bp[2],$3,$4,$5,$6,$7,$8,$9,$10 ;
      } else {
         print $1,$2-1,$2,$3,$4,$5,$6,$7,$8,$9,$10
      } }' |

   $bedtools intersect -a stdin -b $OUT.blocks.bed | \

   awk -v FS="\t" -v OFS="\t" '
      {if($11 ~ "0/0" )
         print $1,$2+1,$4,$5,$6,$7,$8,"END="$3,$10,$11
      else
         print $1,$2+1,$4,$5,$6,$7,$8,$9,$10,$11
      }' > $OUT.noheader


   # sort option in reverse was to avoid SNP being positioned before the gblock and missed in fasta calling
   $samtools view -H $DIRBAM/$OUT/$OUT.realigned.bam | grep "@SQ"| awk '{print $2}' | sed  s/SN:// > chr.names
   $bgzip $OUT.header
   mv $OUT.header.gz $OUT.final.gvcf.gz
   for chrom in $(cat chr.names); do
       awk -v chrom=$chrom '$1==chrom' $OUT.noheader | sort -k2,2n -k10,10 | $bgzip -c >> $OUT.final.gvcf.gz
   done

   rm $OUT.blocks.bed $OUT.noheader $OUT.gvcf.gz $OUT.flt.gvcf.gz chr.names
fi



################################################################################
# 5. SUMMARY STATISTICS
################################################################################
if [ $OPTION = "-qual" ]; then
   cd $DIRVCF
   Q=`awk '$2>0 {a+=$1*$2;n+=$1} END {print a/n}' "$DIRBAM/$OUT/$OUT.realigned.depth" `
   MAXCOV=`echo "tmp=2*($Q+0.5); tmp /= 1; tmp" | bc`
   echo $Q $MAXCOV

   # total no. of positions aligned, 0/1 and 1/1 SNPs
   N00=` zcat $OUT.final.gvcf.gz | grep '0/0' | sed s/"END="// | sed s/";"/" "/  | awk '{n00+=($8-$2+1)} END {print n00}'`
   N01=` zcat $OUT.final.gvcf.gz | grep '0/1' | grep -v 'INDEL' | wc -l`
   N11=` zcat $OUT.final.gvcf.gz | grep '1/1'  | grep -v 'INDEL' | wc -l`
   ND=`  zcat $OUT.final.gvcf.gz | grep 'INDEL' | wc -l`
   echo 'N Homozygous ref: ' $N00
   L=`echo "($N00+$N01+$N11)"|bc`
   q1=`echo "scale=6; $N01/($N00+$N01+$N11)"|bc`
   echo 'N Heterozygous:   ' $N01 $q1
   q2=`echo "scale=6; $N11/($N00+$N01+$N11)"|bc`
   q3=`echo "scale=6; $ND/($N00+$N01+$N11)"|bc`
   echo 'N Homozygous alt: ' $N11 $q2

   #-->  computes some SNP statistics by chr: total snps, heteroz snps, fixed snp, indels
   zcat $OUT.final.gvcf.gz | grep -v "INDEL" | grep -v '#'| grep "0/0" | awk '{print $1}' | sort | uniq -c | sort -k2 > 1.$OUT
   zcat $OUT.final.gvcf.gz | grep -v "INDEL" | grep -v '#'| grep "0/1" | awk '{print $1}' | sort | uniq -c | sort -k2 > 2.$OUT
   zcat $OUT.final.gvcf.gz | grep -v "INDEL" | grep -v '#'| grep "1/1" | awk '{print $1}' | sort | uniq -c | sort -k2 > 3.$OUT
   zcat $OUT.final.gvcf.gz | grep    "INDEL" | grep -v '#'|              awk '{print $1}' | sort | uniq -c | sort -k2 > 4.$OUT
   $samtools depth -q $BASEQ -Q $MAPQ $DIRBAM/$OUT/$OUT.realigned.bam | \
        awk '($3>='$MINCOV' && $3<='$MAXCOV') {print $1}' | sort | uniq -c | sort -k2 > 5.$OUT

   # merging
   join -1 2 -2 2  1.$OUT 2.$OUT > q1.$OUT
   join -1 2 -2 2  3.$OUT 4.$OUT > q2.$OUT
   join q1.$OUT q2.$OUT > q.$OUT
   join -1 2 5.$OUT q.$OUT | awk '{print $0, $3/$2,$4/$2,$5/$2,$6/$2}' | sort -n > $OUT.stat
   echo 'All' $L $N01 $N11 $ND $q1 $q2 $q3 >> $OUT.stat
   rm q?.$OUT ?.$OUT
fi


if [ $OPTION = "-pdf" ]; then
   cd $DIRVCF
   Q=`awk '$2>0 {a+=$1*$2;n+=$1} END {print a/n}' "$DIRBAM/$OUT/$OUT.realigned.depth" `
   MAXCOV=`echo "tmp=2*($Q+0.5); tmp /= 1; tmp" | bc`
   echo '-qual' $Q $MAXCOV

   zcat $OUT.final.gvcf | grep -v '0/0' | grep -v '\./\.' > $OUT.tmp.gvcf

   #--> computes diversity by windows (v3)
   time $samtools mpileup -Bq $MAPQ -Q $BASEQ $DIRBAM/$OUT/$OUT.realigned.bam  | \
        awk -v OFS="\t" "(\$4>=$MINCOV && \$4<=$MAXCOV) {print \$1,\$2,\$3,\$4,\$5} " | \
        perl $DIRBIN/covXwin-v3.1a.pl -v $OUT.tmp.gvcf -w $WINSIZE -d $MINCOV -m $MAXCOV -b $DIRBAM/$OUT/$OUT.realigned.bam | \
        $DIRBIN/ngs_theta -d $MINCOV -m $MAXCOV > $OUT.wintheta

   # plots
   cp $OUT.wintheta sample.wintheta
   R CMD BATCH $DIRBIN/plotNGS.R
   mv plot.pdf $OUT.wintheta.pdf
   rm sample.wintheta $OUT.tmp.gvcf
fi

################################################################################
# 6. CONVERTS TO FASTA (RM INDELS)
################################################################################
# SNPs next to indels are removed
if [ $OPTION = "-fasta" ]; then
   cd $DIRFAS
   time zcat $DIRVCF/$OUT.final.gvcf.gz | grep -v "#" | grep -v INDEL | grep -v "\./\." | sed s/"END="// | \
        awk -v FS="\t" -v OFS="\t" '{print $1,$2,$4,$5,$10,$8}' | cut -d";" -f1 | \
        awk -v FS="\t" -v OFS="\t" '
        {if( $5 ~ "0/0" )
           print $1,$2-1,$6,$3,$3,$5
        else
       	   print $1,$2-1,$2,$3,$4,$5
        }'  | sed s/':'/' '/ | sed s/"\/"/' '/ | sed s/"|"/' '/ | sed s/,//g | \
        awk '($4!=$5 && ($6+$7)>0) || ($6==0 && $7==0)' | \
        $DIRBIN/gvcf2fas -f $DIRASSEMBLY/$ASSEMBLY.fa -nocheck | $bgzip -c > $OUT.fa.gz
        $samtools faidx $OUT.fa.gz
fi


################################################################################
# 7. MERGES FASTA FROM SEVERAL INDIVIDUALS INTO VCF 
################################################################################
# falist is a file containing, for each row, the name of the fasta file and the name of the sample
if [ $OPTION = "-fas2vcf" ]; then
   cd $DIRFAS
   falist=$2
   outfile=$3
   $DIRBIN/fas2vcf -i $falist -f $DIRASSEMBLY/$ASSEMBLY.fa | sed s/' '//g | $bgzip -c > $outfile.vcf.gz
fi

################################################################################
# 8. CLEAN (excl bam file)
################################################################################
if [ $OPTION = "-clean" ]; then
   rm $DIRDATA/$OUT*
   cd $DIRBAM/$OUT
   rm $OUT.ready.* $OUT.intervals $OUT.bam $OUT.rl.* $OUT.rmdup.bam $OUT.sort.* metrics.out
fi

exit 0
