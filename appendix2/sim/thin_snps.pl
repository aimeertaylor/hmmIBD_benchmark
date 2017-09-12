#!/usr/bin/perl
use strict;
use warnings;

my $thinfact = shift;
my $infile = shift;
my $outfile = shift;
if (!defined $thinfact || !defined $infile || !defined $outfile) {
  die "usage: thin_snps.pl <N, will keep every 1/N snps> <input freq file> <output freq file>\n";
}

open(STATE, "sim_state.txt") || die "could not open state file\n";
open(OSTATE, ">scratch/sim_state_thin.txt") || die "could not open output state file\n";
open(SIM, "sim_model.txt") || die "could not open sim file\n";
open(OSIM, ">scratch/sim_model_thin.txt") || die "could not open output sim file\n";
open(FREQ, $infile) || die "could not open freq file\n";
open(OFREQ, ">$outfile") || die "could not open output freq file\n";

my $head = <SIM>;
print OSIM $head;
$head = <STATE>;
print OSTATE $head;

my $iline = 0;
while (my $simline = <SIM>) {
  my $stateline = <STATE>;
  my $freqline = <FREQ>;
  if ($iline % $thinfact == 0) {
    print OSIM $simline;
    print OSTATE $stateline;
    print OFREQ $freqline;
  }
  $iline++;
}
