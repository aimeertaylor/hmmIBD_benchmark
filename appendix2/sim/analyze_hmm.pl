#!/usr/bin/perl -w
use strict;

my $sumfile = shift;
my $hmmfile = shift;
my $vitfile = shift;
my $outfile = shift;
if (!defined $sumfile || !defined $outfile || !defined $hmmfile || !defined $vitfile) {
  die "usage: analyze_hmm.pl <sim sum file> <hmm outfile> <viterbi comp file< <outfile (will append results)>\n";
}
open(OUT, ">>$outfile") || die "could not open outfile $outfile\n";

# "Truth" (i.e. simulated) data
my $file = $sumfile;
open(IN, $file) || die "could not open sim $file\n";
my $head = <IN>;
my %index;
my $ind = 0;
my (@k, @eps, @pi, @ntrans);
while (my $line = <IN>) {
  chomp($line);
  my @pieces = split /\t/, $line;
  my $samp1 = $pieces[0];
  my $samp2 = $pieces[1];
  $index{$samp1}{$samp2} = $ind;
  push @k, $pieces[3];
  push @eps, $pieces[4];
  push @pi, $pieces[5];
  push @ntrans, $pieces[6];
  $ind++;
}
print "n pairs generated: ", scalar(@pi), "\n";

# Viterbi comparison
my %nVitRight;
my %nVitWrong;
open(IN, $vitfile) || die "could not open Viterbi file $vitfile\n";
$head = <IN>;
while (my $line = <IN>) {
  chomp $line;
  my @pieces = split /\t/, $line;
  my $samp1 = $pieces[0];
  my $samp2 = $pieces [1];
  $nVitRight{$samp1}->{$samp2} = $pieces[2];
  $nVitWrong{$samp1}->{$samp2} = $pieces[3];
}

# HMM results
$file = $hmmfile;
open(IN, "$file") || die "could not open fract $file\n";
$head = <IN>;
while (my $line = <IN>) {
  chomp($line);
  my @pieces = split /\t/, $line;
  my @samps = @pieces[0..1];
  my @subsamp = split /_/, $samps[1];
#  if ($samps[0] eq $subsamp[0]) {
    my $piH = $pieces[9];
    my $kH = $pieces[6];
    my $ntransH = $pieces[7];
#    my $delpi = $pieces[8];
#    my $del_k = $pieces[9];
#    my $del_prob = $pieces[10];
    my $nfit = $pieces[5];
    my $index = $index{$samps[0]}{$samps[1]};
    if (!defined $index) {print "samps[0]: $samps[0], subsamp: $subsamp[1]\n";}
    print OUT "$pieces[0]\t$eps[$index]\t$k[$index]\t$pi[$index]\t$ntrans[$index]\t$kH\t$piH\t$ntransH\t$nfit";
    my $nvr = $nVitRight{$samps[0]}->{$samps[1]};
    my $nvw = $nVitWrong{$samps[0]}->{$samps[1]};
    printf OUT "\t%.4f\n",  $nvr / ($nvr + $nvw); 
#  }
}
