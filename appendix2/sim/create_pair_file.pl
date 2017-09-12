#!/usr/bin/perl
use strict;
use warnings;

my $file = shift;
if (!defined $file) {die "gimme a file name\n";}
open(IN, $file) || die "could not open $file\n";
open(OUT, ">scratch/good_sim_pair.txt") || die "could not open...\n";

my $head = <IN>;
chomp $head;
my @pieces = split /\t/, $head;
for (my $ip = 2; $ip < @pieces; $ip += 2) {
  print OUT "$pieces[$ip]\t$pieces[$ip+1]\n";
}
