#!/usr/bin/perl
# joins in bed-format blocks (chr begin-1 end), tab separated, a list of consecutive positions
# MPE: miguel.perez@uab.es
use warnings;
use strict;

my ($begin, $chr, $last_chr);
my $last_pos=0;
my $pos=1;
while (<>) {
   chomp;
   ($chr,$pos) = split("\t");
   if (!$begin) {
      $begin = $pos-1;
      $last_chr = $chr;
   } elsif ($pos != $last_pos+1) {
      print "$last_chr\t$begin\t$last_pos\n";
      $begin=$pos-1;
   }
   $last_pos = $pos;
   $last_chr = $chr;
}
print "$last_chr\t$begin\t$last_pos\n";
