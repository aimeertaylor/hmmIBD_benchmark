#!/usr/bin/perl
use strict;
use warnings;

my $comm = "../hmmIBD -i data/Cambodia_Ghana_merged_seq.txt -f data/freqs_Cambodia_Ghana.txt -o onepop -m 40 -n 100 -g data/Cambodia_Ghana_sample_pairs.txt";
my $resp = `$comm`;
print $resp;

$resp = `./sum_overlaps onepop.hmm.txt`;
print $resp;

$resp = `./convert_overlap.pl onepop.hmm.txt.summed onepop.overlaps`;
print $resp;
