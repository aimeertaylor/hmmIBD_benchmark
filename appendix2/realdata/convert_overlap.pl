#!/usr/bin/perl -w
use strict;

# Change overlap output to have start and stop on one line

my $infile = shift;
my $outfile = shift;
if (!defined $infile || !defined $outfile) {
  die "usage: convert_overlap.pl <input overlap file> <output file>\n";
}

open(IN, "$infile") || die "could not open $infile\n";
open(OUT, ">$outfile") || die "could not open $outfile\n";

<IN>;
print OUT "chrom\tstart\tstop\tN_ident_pair\tlocal_start\tlocal_stop\n";
while (my $line1 = <IN>) {
  my $line2 = <IN>;
  chomp($line1);
  chomp($line2);
  my @pieces1 = split /\t/, $line1;
  my @pieces2 = split /\t/, $line2;
  if ($pieces1[3] != $pieces2[3]) {die "$line1:$line2\n";}
  print OUT "$pieces1[0]\t$pieces1[2]\t$pieces2[2]\t$pieces1[3]\t$pieces1[1]\t$pieces2[1]\n";
}
