my $USAGE = <<END_USAGE;

 CovXwin v-3

 M Perez-Enciso (miguel.perez@uab.es) & A Esteve
 last mod:
 last mod: Tue 22 Oct 2013 01:47:00 PM CEST
 last  v: Tue 12 Feb 2013 18:37:57 CET
 first v: 20110712

 Computes average coverage, depth and no. of variants per window, as described in
    Partial short-read sequencing of a highly inbred Iberian pig and genomics inference thereof
    Esteve-Codina et al. 2011, Heredity, 107:256-64. doi: 10.1038/hdy.2011.13
    http://www.nature.com/hdy/journal/v107/n3/full/hdy201113a.html
 Please cite the reference if you use this program.
 To be used jointly with ngs_theta

 USAGE

  samtools mpileup -BQ \$MINQ -d 100000 file.bam | perl covXwin [..options] | ngs_theta -d min_depth -m max_depth > file.wintheta

 OPTIONS

    -v=s   vcf file (required)
    -b=s   bamfile, used to get chr names and lengths (required)
    -w:i   window size in bp [100000]
    -d:i   minimum depth [3]
    -m:i   maximum depth [30]
    -q:i   minimum SNP quality

 NOTES

   - vcf file is taken 'as is' ie no further filtering is done,
     filtering by quality, ambiguouos calls, multiallelic snps is highly recommended
   - output is in format
        chr window n_heterozygous_snps n_fixed_snps n_indels {bp_at_depth[i], i=min_deth..max_depth}
   - this output is further processed by ngs_theta and needs not be saved
   - analyzes only one sample at a time
   - non overlapping windows
   - indel statistics ususally not too reliable, printed but not considered in theta estimates
   - care with interpreting HKA, as it takes the reference genome as an outgroup

END_USAGE

#! /usr/local/bin/perl -w
use warnings;
use strict;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

# default values
my $vcffile;
my $bamfile;
my $window_size=100000;
my $min_cov=3;
my $max_cov=30;
my $help;
my $qsnp=0;
GetOptions('v=s'       => \$vcffile,
           'w=i'       => \$window_size,
           'd=i'       => \$min_cov,
           'm=i'       => \$max_cov,
           'bamfile=s' => \$bamfile,
           'q:i'       => \$qsnp,
	       'help'      => \$help) ;

#-- prints usage if no var file is defined or both formats defined

die $USAGE if ( !$vcffile || $help);

# all data are stored in a hash of hashes where primary key has form $chr $win
my %HoH;
my @chrs;

#--> extracts chr names and lengths from bam file
my @header =qx(samtools view -H $bamfile);
foreach (@header) {
    next if ($_ !~ m/\@SQ/);
    $_ =~ /SN:(\w+)\tLN:(\d+)/;
    my $chr=$1;
    push (@chrs,$chr);
    my $nwin=int($2/$window_size)+1;
    #--> initializes each window
    foreach my $iwin (1..$nwin) {
        $HoH{$chr}{$iwin}{SNPh}=0;
        $HoH{$chr}{$iwin}{SNPf}=0;
        $HoH{$chr}{$iwin}{indel}=0;
        foreach ($min_cov..$max_cov) {$HoH{$chr}{$iwin}{$_}=0}
    }
}

#----------> reads pileup file from STID
while (<>) {
       chomp;
       my @info = split;
       my $iwin = int($info[1]/$window_size)+1;
       my $chr = $info[0];
       my $d = $info[3];  #--> depth
       $HoH{$chr}{$iwin}{$d}++ if ($min_cov<=$d && $d<=$max_cov);
}

#-----------> variant file is read (only mpileup)
open (INFILE,'<',$vcffile) || die "cannot open $vcffile for reading: $!";

while (<INFILE>){
       my $line = $_;
       next if ($line =~ /^#/); #--> skips comments
       my @info = split;
       if (looks_like_number($info[5])) {next if ($info[5] < $qsnp)};
       my $iwin = int($info[1]/$window_size)+1;
       my $chr = $info[0];
       if ($line =~ /INDEL/) {
	      $HoH{$chr}{$iwin}{indel}++;
       }
       elsif ($line =~ /0\/1/ ) {
	      $HoH{$chr}{$iwin}{SNPh}++;
       }
       elsif ($line =~ /1\/1/ ) {
	      $HoH{$chr}{$iwin}{SNPf}++;
       }
}


#--------> prints output that can be processed by program ngs_theta
foreach my $chr (@chrs) {
   #--> sort by pos
   foreach my $iwin (sort {$a<=>$b} keys %{$HoH{$chr}} ) {
       print "$chr\t$iwin\t$HoH{$chr}{$iwin}{SNPh}\t$HoH{$chr}{$iwin}{SNPf}\t$HoH{$chr}{$iwin}{indel}\t" ;
       for my $i ($min_cov..$max_cov) {
           print "$HoH{$chr}{$iwin}{$i}\t";
       }
       print "\n";
   }
}
