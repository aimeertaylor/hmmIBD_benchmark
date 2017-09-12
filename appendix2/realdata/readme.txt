Steps needed to recreate IBD analysis of P. falciparum chromosome 7. Working executable of hmmIBD 
assumed to be in the parent directory.

From command line:

- Compile sum_overlaps.c. A sample command to do so is in script 'complink'.
- Execute run_2pop.pl. Runs hmmIBD on the two populations separately and tabulates results in twopop.overlaps.
- Execute run_1pop.pl. Runs hmmIBD on the merged file, with highly divergent sites filtered out, and using 
  averaged allele frequencies. Tabulates results in onepop.overlaps.

Display results in R (if you want to). The R functions are hacky and contain hard-wired values. Commands in R:

- setwd("[path]/realdata")
- source("plot_chrom7.R")
- plot_both(to_file). The argument to_file determines whether results are sent to pdf file (to_file=1) or plotted 
on screen (anything else).

Note: hmmIBD in run_Xpop.pl is invoked with a limit of 100 generations set via the -n option, which was not 
done in the version used in the paper. It makes very little difference.
