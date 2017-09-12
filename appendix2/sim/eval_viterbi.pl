#!/usr/bin/perl
use strict;
use warnings;

my $statefile = shift;
my $hmmfile = shift;
my $outfile = shift;
if (!defined $statefile || !defined $outfile || !defined $hmmfile) {
  die "usage: eval_viterbi.pl <sim state file> <hmm outfile> <outfile (will append results)>\n";
}
open(HMM, $hmmfile) || die "could not open HMM file $hmmfile\n";
my $head = <HMM>;

# samp2 is unique for each samp1
my %seg_start; # $seg{$samp1}->[$chr][$index]
my %seg_end;
my %seg_state;   # 0 => ibd
while (my $line = <HMM>) {
  chomp $line;
  my @pieces = split /\t/, $line;
  my $samp1 = $pieces[0];
  my $samp2 = $pieces[1];
  my $chr = $pieces[2];
  my $start = $pieces[3];
  my $end = $pieces[4];
  my $state = $pieces[5];
  push @{$seg_start{$samp1}->[$chr]}, $start;
  push @{$seg_end{$samp1}->[$chr]}, $end;
  push @{$seg_state{$samp1}->[$chr]}, $state;
}
close HMM;

my @samp1;
my @samp2;
open(IN, $statefile) || die "could not open state file $statefile\n";
$head = <IN>;
chomp $head;
my @pieces = split /\t/, $head;
for (my $ind = 0; $ind < @pieces - 2; $ind++) {
  $pieces[$ind+2] =~ m/(\w+)\:(\w+)/;
  $samp1[$ind] = $1;
  $samp2[$ind] = $2;
#  print "$samp1[$ind], $samp2[$ind]\n";
}

my @nright = (0) x scalar(@samp1);
my @nwrong = (0) x scalar(@samp1);
my $oldchr = -1;
my @segind ;
while (my $line = <IN>) {
  chomp $line;
  my @pieces = split /\t/, $line;
  my $chr = $pieces[0];
  my $pos = $pieces[1];
  for (my $ipair = 0; $ipair < @samp1; $ipair++) {
    my $s1 = $samp1[$ipair];
    if ($oldchr != $chr) {$segind[$ipair] = 0;}
    if (defined $seg_end{$s1}->[$chr][$segind[$ipair]] && $pos > $seg_end{$s1}->[$chr][$segind[$ipair]]) {
      # have walked past end of previous segment
      $segind[$ipair]++;
    }
    # There can be snps outside any segment: skipped (for spacing) snps at end of chrom or between segs
    if (defined $seg_state{$s1}->[$chr][$segind[$ipair]] && $pos >= $seg_start{$s1}->[$chr][$segind[$ipair]]) {
      if ($pieces[2+$ipair] == $seg_state{$s1}->[$chr][$segind[$ipair]]) {$nright[$ipair]++;}
      else {$nwrong[$ipair]++;}
    }
  }
  $oldchr = $chr;
}

open(OUT, ">$outfile") || die "could not open outfile $outfile\n";
print OUT "samp1\tsamp2\tNright\tNwrong\n";
for (my $ipair = 0; $ipair < @samp1; $ipair++) {
  print OUT "$samp1[$ipair]\t$samp2[$ipair]\t$nright[$ipair]\t$nwrong[$ipair]\n";
}


