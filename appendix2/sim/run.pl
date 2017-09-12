#!/usr/bin/perl -w
use strict;

my @kval = (1,3,5,10,20,50);
my @pival = (0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.95);
my @iterval = (1, 2, 3, 4, 5, 10, 20, 40);

my $eps = 0.001;

my $freqfile = "data/freqs_Ghana_biall.txt";
my $filebase = "results/model";
my @fileh;
for (my $ifile = 0; $ifile < @iterval; $ifile++) {
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

    for (my $ifile = 0; $ifile < @iterval; $ifile++) {
      $comm = "../hmmIBD -i sim_model.txt -o scratch/sim$ifile -m $iterval[$ifile] -f $freqfile";
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
