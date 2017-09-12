#!/usr/bin/perl -w
use strict;

my @pival = (0.5);
my @kval = (10);
my @smudgeval = (1000, 400, 200, 100, 50, 25, 10);

my $eps = 0.001;

my $filebase = "results_smudged/model";
my @fileh;
for (my $ifile = 0; $ifile < @smudgeval; $ifile++) {
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
    my $comm = "./sim_model -p $pival[$ip] -k $kval[$ik] -e $eps -f data/freqs_Ghana_biall.txt";
    my $resp = `$comm`;
    print $resp;
    $comm = "./create_pair_file.pl sim_model.txt";
    $resp = `$comm`;
    print $resp;

    for (my $ifile = 0; $ifile < @smudgeval; $ifile++) {
      $comm = "../hmmIBD -i sim_model.txt -o scratch/sim$ifile -m 40 -f data/freqs_Ghana_$smudgeval[$ifile].txt";
      $comm .= " -g scratch/good_sim_pair.txt";
      $resp = `$comm`;
      print $resp;
      $comm = "./eval_viterbi.pl sim_state.txt scratch/sim$ifile.hmm.txt scratch/vit_match.txt";
      $resp = `$comm`;
      print $resp;
      $comm = "analyze_hmm.pl simsum_model.txt scratch/sim$ifile.hmm_fract.txt scratch/vit_match.txt $filebase$ifile.results.txt";
      $resp = `$comm`;
      print $resp;
    }
  }
  `rm -f scratch/good_sim_pair.txt`;
  `rm -f sim_model.txt`;
  `rm -f simsum_model.txt`;
  `rm -f sim_state.txt`;
}
