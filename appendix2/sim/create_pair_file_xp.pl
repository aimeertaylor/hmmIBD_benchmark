#!/usr/bin/perl
use strict;
use warnings;

my $file1 = shift;
my $file2 = shift;
if (!defined $file1 || !defined $file2) {die "gimme 2 file names\n";}
open(IN1, $file1) || die "could not open $file1\n";
open(IN2, $file2) || die "could not open $file2\n";
open(OUT, ">scratch/good_sim_pair_xp.txt") || die "could not open...\n";

my $head1 = <IN1>;
chomp $head1;
my $head2 = <IN2>;
chomp $head2;
my @pieces1 = split /\t/, $head1;
my @pieces2 = split /\t/, $head2;
for (my $ip = 2; $ip < @pieces1; $ip++) {
  print OUT "$pieces1[$ip]\t$pieces2[$ip]\n";
}
