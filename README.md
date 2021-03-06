##  NGSpipeline: NGS pipeline to infer variability and merge gvcf files
#### Scripts to analyze NGS data and multiple individual coherent SNP calling
#### Miguel Pérez-Enciso (miguel.perez@uab.es) and Jordi Leno-Colorado (jordi.leno@cragenomica.es)

A problem with simultaneous SNP calling across samples is to distinguish between missing positions and bases equal to the reference. Furthermore, if depth is very different, it is desirable to filter all positions by the same restrictions, say in depth. This NGS data pipeline was developed to merge vcf files from separate individuals by keeping control on homozygous regions as well. This is done by filtering the raw gvcf file such that the same restrictions on SNP, base, map qualities, minimum and maximum depth are applied to both SNPs and regions where the sample is equal to the reference genome. A second utility converts this gvcf file into a fasta file, with N's where the sample was not sequenced with the restrictions set by the user. Finally, a merged vcf for several samples can be generated back. In this process, individual SNP quality is lost although we plan to develop a new program that prints SNP depth and likelihood in the merged vcf file. INDELs are not considered. A few plots with depth and other statistics can also be generated that help in visualizing possible weird results.

### Warnings and main usage
Some program options used in this pipeline seem to be deprecated, ie, in GATK. The main use of this pipeline is to show
how to merge individual gvcf files such that the same filters are applied to all samples. This is a bit tricky specially
in the homozygous blocks. The main options of the pipeline are therefore (SRR1 is the sample prefix):

* `sh wflow_ngs -vcf SRR1`  (obtains filtered gvcf file)
* `sh wflow_ngs -fasta SRR1` (obtains fasta file from gvcf)
* steps above are repeated for every sample
*  `sh wflow_ngs -fas2vcf FASTA_LIST OUTFILE` 
(obtains a merged vcf file where missing genotypes ./. are distinguished from 0/0)

See below for details.

### Standard software required
 - ascp
 - fasterq-dump
 - bwa
 - samtools (1.8)
 - bcftools (1.8)
 - bgzip
 - gatk 
 - picardtools
 - bedtools
 - R

### Non standard software provided
 - joinConsecutivePos.pl
 - covXwin.pl
 - gcvf2fas.f90
 - ngs_theta.f90 & fact-m.f90
 - plotNGS.R
 - fas2vcf.f90 

### Installing
Within the working directory, you should have an ASSEMBLY folder with reference genome in fasta format, BIN folder with executables, a DATA folder with fq reads, BAMFILES folder, a VARFILES and a FASTAFILES folders. This structure can be changed in the script.
To compile f90 programs:

```
f95 -O4 gvcf2fas.f90 -o gvcf2fas
f95 -O4 fact-m.f90 ngs_theta.f90 -o ngs_theta
f95 -O4 fas2vcf.f90 -o fas2vcf
```

Paths to all programs should be specified in the shell script, eg,

```
bwa=bwa
samtools=$DIRBIN/samtools
GATK=$DIRBIN/GenomeAnalysisTK.jar
picard=$DIRBIN/picard.jar
bedtools=bedtools
bcftools=$DIRBIN/bcftools
bgzip=bgzip
fastqdump=fasterq-dump
ASPERA=~/.aspera
```

### Parameters
```
MINCOV=5        # minimum depth required
MAXCOV=         # maximum depth, by default is (2 x mean_depth + 1)
SNPQ=10         # min snp quality
MAPQ=20         # min map quality
BASEQ=20        # min base quality
NP=10           # no. of threads
WINSIZE=100000  # window size used to comple plots (only for -pdf option)
```


### Running the pipeline

#### Indexing

   `sh wflow_ngs -index`

#### To download reads

   `sh wflow_ngs -sra2fq SRR_1 SRR_2 ... SRR_n`

SRR_i are the SRR ids for a given sample, which are all merged in a single fq paired end file with names **SRR1\_sra\_1.fastq** and **SRR1\_sra\_2.fastq**. Reads are stored in DATA folder. 

#### To align with bwa
In the following, SRR1 represents the sample thta is being analyzed. To align with bwa and refine the alignment do

   `sh wflow_ngs -bwa SRR1`

Aligns with bwa, realigns around indels with GATK and remove duplicates with picard. It also computes a file with number of bases sequenced at a given depth. Produces files **SRR1.realigned.bam**, **SRR1.realigned.bam.bai** and **SRR1.realigned.depth** in directory BAMFILES/SRR1

#### To obtain and filter gvcf file (SNP calling)

For the whole genome jointly:

   `sh wflow_ngs -vcf SRR1`

For each chromosome separately:

```
sh wflow_ngs -cvcf SRR1
# Once finished
sh wflow_ngs -cmerge SRR1
```

This option requires the additional script **wflow\_ngs\_vcf\_chr.sh**

Produces file **SRR1.final.gvcf.gz** in folder VARFILES

#### To get a fasta file from gvcf

   `sh wflow_ngs -fasta SRR1`

Produces **SRR1.fa.gz** file in FASTAFILES folder.

#### To get a vcf file from multiple individual gvcf files
First you need to produce an **uncompressed** fasta file for each sample using the command above.

   `sh wflow_ngs -fas2vcf FASTA_LIST OUTFILE`

Produces **OUTFILE.vcf.gz** file in FASTAFILES folder. File FASTA_LIST contains the name of the fasta file and sample name for each sample, eg,

```
sample1.fa  sample1
sample2.fa  sample2
...
```

#### To get some stats and quality metrics

   `sh wflow_ngs -qual SRR1`

Produces a file **SRR1.stats**, stored in VARFILES folder, with total lengths sequenced, fixed positions (0/0), heteroz (0/1), fixed differences (1/1), and indels per chromosome.

   `sh wflow_ngs -pdf SRR1`

Produces a text file **SRR1.wintheta** and plots (**SRR1.pdf**) with depths per chr per window of size WINSIZE togther with several rough estimates of variability as in [Esteve-Codina et al. (2013)](https://www.ncbi.nlm.nih.gov/pubmed/23497037). These plots can be mainly be used to inspect whether depth is uniform or detect some weird patterns in terms of variability, etc.

### Citations
If you find these scripts useful, please cite:

[Leno-Colorado J, Hudson NJ, Reverter A, Pérez-Enciso M. 2017, A Pathway-Centered Analysis of Pig Domestication and Breeding in Eurasia. G3 7(7):2171-2184](http://www.g3journal.org/content/7/7/2171.long)

[E Bianco, B Nevado, SE Ramos-Onsins, M Pérez-Enciso, 2015. A Deep Catalog of Autosomal Single Nucleotide Variation in the Pig, PlosOne](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0118867)

If you happen to employ the statistics produced in option -pdf:

[Esteve-Codina A, Paudel Y, Ferretti L, Raineri E, Megens HJ, Silió L, Rodríguez MC, Groenen MA, Ramos-Onsins SE, Pérez-Enciso M. 2013. Dissecting structural and nucleotide genome-wide variation in inbred Iberian pigs, BMC Genomics](https://www.ncbi.nlm.nih.gov/pubmed/23497037)
