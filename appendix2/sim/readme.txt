Steps needed to recreate IBD analysis of data simulated under the HMM. Working executable of hmmIBD 
assumed to be in the parent directory.

From command line:

- Compile simulation code by executing script 'complink'.
- Generate and analyze single-population simulation data by executing script 'run.pl'. 
Output should appear in results/.
- Generate and analyze cross-population simulation data by executing 'runXP.pl'. Results also 
in results/.
- Same for testing results of reducing number of SNPs, by executing 'run_thinned.pl'. Results
appear in results_thinned/. Timing info should be appear in timing.txt and stdout output (including
number of variants used) in log.txt.
- Test effect of error in estimated allele frequency: 'run_smudged.pl'. The re-estimated allele frequencies 
are already in files in data/.

Most of the plots in the paper can be recreated using the various functions in plotSims.R. Commands in R:
- setwd("[path]/sim")
- source("plotSims.R")
- execute any of the plotting functions, e.g.
- plot_piErr(to_file). The argument to_file determines whether results are sent to pdf file (to_file=1) or plotted 
on screen (anything else).
