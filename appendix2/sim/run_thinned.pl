#!/usr/bin/perl -w
use strict;

my @kval = (10);
my @pival = (0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.95);
#my @pival = (0.5);
#my @thinval = (100);
my @thinval = (1, 2, 4, 10, 20, 100, 500);
my $eps = 0.001;
my $freqfile = "data/freqs_Ghana_biall.txt";

open(STDERR, ">results_thinned/timing.txt") || die "could not open timing file\n";
open(STDOUT, ">results_thinned/log.txt") || die "could not open log file\n";
my $filebase = "results_thinned/model";
my @fileh;
for (my $ifile = 0; $ifile < @thinval; $ifile++) {
  #localize the file glob, so FILE is unique to
    #    the inner loop.
  local *FILE;
  my $filename = "$filebase$ifile.results.txt";
  open(FILE, ">$filename") || die "could not open $filename\n";
  push @fileh, *FILE;
  print FILE "sample\terror_gen\tk_gen\tpi_gen\tntrans_gen\tk_obs";
  print FILE "\tpi_obs\tntrans_obs\tn_fit_iter\tfract_right_state\n";
}

for (my $ik = 0; $ik < @kval; $ik++) {
  for (my $ip = 0; $ip < @pival; $ip++) {
    my $comm = "./sim_model -p $pival[$ip] -k $kval[$ik] -e $eps -f $freqfile";
    my $resp = `$comm`;
    print $resp;
    $comm = "./create_pair_file.pl sim_model.txt";
    $resp = `$comm`;
    print $resp;

    for (my $ifile = 0; $ifile < @thinval; $ifile++) {
      $comm = "thin_snps.pl $thinval[$ifile] $freqfile scratch/freqs_Ghana_thin.txt";
      $resp = `$comm`;

      $comm = "time ../hmmIBD -i scratch/sim_model_thin.txt -o scratch/sim$ifile -m 40 -f scratch/freqs_Ghana_thin.txt";
      $comm .= " -g scratch/good_sim_pair.txt";
      $resp = `$comm`;
      print $resp;
      $comm = "./eval_viterbi.pl scratch/sim_state_thin.txt scratch/sim$ifile.hmm.txt scratch/vit_match.txt";
      $resp = `$comm`;
      print $resp;
      $comm = "analyze_hmm.pl simsum_model.txt scratch/sim$ifile.hmm_fract.txt scratch/vit_match.txt $filebase$ifile.results.txt";
      $resp = `$comm`;
      print $resp;
    }
  }
  `rm -f sim_model.txt`;
  `rm -f simsum_model.txt`;
  `rm -f sim_state.txt`;
  `rm -f scratch/freqs_Ghana_thin.txt`;
  `rm -f scratch/sim_model_thin.txt`;
  `rm -f scratch/sim_state_thin.txt`;
}
