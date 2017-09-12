#!/usr/bin/perl -w
use strict;

my @kval = (10, 50);
my @pival = (0.3, 0.5, 0.7);
my @iterval = (1, 2, 3, 4, 5, 10, 20, 40);

my $eps = 0.001;

my $filebase = "results/xp";
my @fileh;
for (my $ifile = 0; $ifile < @iterval; $ifile++) {
  #localize the file glob, so FILE is unique to
    #    the inner loop.
  local *FILE;
  my $filename = "$filebase$ifile.results.txt";
  open(FILE, ">$filename") || die "could not open $filename\n";
  push @fileh, *FILE;
  print FILE "sample\terror_gen\tk_gen\tpi_gen\tntrans_gen\tk_obs";
  print FILE "\tpi_obs\tntrans_obs\tn_fit_iter\tfract_right_state\n";
}

my $pop1 = "Senegal";
my $pop2 = "Thailand";
my $file1 = "simXP_$pop1.txt";
my $file2 = "simXP_$pop2.txt";
my $freq1 = "data/freqs_pf3k_".$pop1."_biall.txt";
my $freq2 = "data/freqs_pf3k_".$pop2."_biall.txt";

for (my $ik = 0; $ik < @kval; $ik++) {
  for (my $ip = 0; $ip < @pival; $ip++) {
    my $comm = "./sim_modelXP -p $pival[$ip] -k $kval[$ik] -e $eps -f $freq1 -F $freq2 -o $file1 -O $file2";
    my $resp = `$comm`;
    print $resp;
    $comm = "./create_pair_file_xp.pl $file1 $file2";
    $resp = `$comm`;
    print $resp;

    for (my $ifile = 0; $ifile < @iterval; $ifile++) {
      $comm = "../hmmIBD -i $file1 -I $file2 -o scratch/simXP$ifile -m $iterval[$ifile] -f $freq1 -F $freq2";
      $comm .= " -g scratch/good_sim_pair_xp.txt";
      $resp = `$comm`;
      print $resp;
      $comm = "./eval_viterbi.pl sim_state.txt scratch/simXP$ifile.hmm.txt scratch/vit_match.txt";
      $resp = `$comm`;
      print $resp;
      $comm = "analyze_hmm.pl simsumXP.txt scratch/simXP$ifile.hmm_fract.txt scratch/vit_match.txt $filebase$ifile.results.txt";
      $resp = `$comm`;
      print $resp;
    }
  }
  `rm -f scratch/good_sim_pair_xp.txt`;
  `rm -f $file1`;
  `rm -f $file2`;
  `rm -f simsumXP.txt`;
  `rm -f sim_state.txt`;
}
