#!/usr/bin/perl
use strict;
use warnings;

# Note: should give the same results without providing the frequency files
my $comm = "../hmmIBD -i data/pf3k_Cambodia_seq_7.txt -f data/freqs_Cambodia_seq_7.txt -I data/pf3k_Ghana_seq_7.txt -F data/freqs_Ghana_seq_7.txt -o twopop -m 40 -n 100 -g data/Cambodia_Ghana_sample_pairs.txt";
my $resp = `$comm`;
print $resp;

$resp = `./sum_overlaps twopop.hmm.txt`;
print $resp;

$resp = `./convert_overlap.pl twopop.hmm.txt.summed twopop.overlaps`;
print $resp;

