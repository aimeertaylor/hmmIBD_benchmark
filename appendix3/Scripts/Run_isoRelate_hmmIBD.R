############################################################################
# Script to run IsoRelate and hmmIBD on chimeric data and samples for timing
# Since, it is not possible to specify comparisons to analyse under isoRelate, 
# all possible pairwise comparisons of the children and their parents were analysed. 
# It is possible to filter by comparisons under hmmIBD, via b_site.txt file.
# In total takes around 9 hours without timing experiments (plus ~ 5 to do timiming experiments)
############################################################################

# -------------------- Set up -------------------- 
rm(list = ls())
library(isoRelate) # load isoRelate library (reinstalled on Aug 9th 2017)
library(data.table) # for fread - to load subsample 
require(tictoc) # to time
tic()
load('../pf3k_chimeric_data/sites.RData')
erroneous <- c(FALSE, TRUE) 

# Magic numbers
Magic_numbers <- list(n_children_to_analyse = 50, # Selected at random 
                      n_timing_samples_to_run = 50, # First n_timing_samples_to_run in list
                      n_timing_repetition = 3, 
                      num_cores = 1,
                      min_snps = 0, 
                      minimum_length_bp = 0, 
                      error_prob = 0.001,  
                      seed_run_isoRelate = 1,
                      isolate_max_missing = 1, # Don't filter isolates
                      MAF_my_genotypes = 0,# Don't filter (already done)
                      snp_max_missing = 1) # Don't filter snps
attach(Magic_numbers, warn.conflicts = FALSE)
save(Magic_numbers, file = '../pf3k_chimeric_data/Run_isoRelate_hmmIBD_magic_numbers.RData')

# -------------------- isoRelate: chimeric data analysis -------------------- 
tic()
for(err in erroneous){
  
  for(site in sites){
    
    suffix <- ifelse(err, paste(site, 'erroneous', sep = '_'), site)
    
    # Load data
    load(sprintf('../pf3k_chimeric_data/parents_children_%s.RData', suffix))
    parents_children <- pedmap
    load(sprintf('../pf3k_chimeric_data/reference_%s.RData', site))
    reference <- pedmap
    
    # Select a sample of children at random 
    set.seed(seed_run_isoRelate) 
    children_to_analyse_ind <- sample(grep(':', parents_children[[1]]$iid), size = n_children_to_analyse)
    children_to_analyse_names <- as.character(parents_children[[1]]$iid[children_to_analyse_ind])
    parents_to_analyse_names <- unique(as.vector(do.call(rbind, strsplit(as.character(children_to_analyse_names), split = ":"))))
    ind <- as.character(parents_children[[1]]$iid) %in% c(children_to_analyse_names, parents_to_analyse_names)
    parents_children_subset <- list(parents_children[[1]][ind,], parents_children[[2]]) 
    
    # Save a list of samples to exclude from hmmIBD analysis
    to_exclude_hmmIBD <- parents_children[[1]]$iid[(!parents_children[[1]]$iid %in% parents_children_subset[[1]]$iid)]
    write.table(to_exclude_hmmIBD, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE,
                file = sprintf('../pf3k_chimeric_data/b_%s.txt', suffix))
    
    # Reformat 
    my_genotypes <- getGenotypes(ped.map = parents_children_subset,
                                 reference.ped.map = reference,
                                 maf = MAF_my_genotypes,
                                 isolate.max.missing = isolate_max_missing, 
                                 snp.max.missing = snp_max_missing, 
                                 chromosomes = NULL,
                                 input.map.distance = "cM",
                                 reference.map.distance = "cM")
    
    # Aside: check reference used to calculate frequencies
    X <- read.delim(sprintf('../pf3k_chimeric_data/pf3k_data_biallelic_frequencies_%s.txt', site),  
                    header = FALSE)
    plot(y = my_genotypes$genotypes$freq, x = X[,'V4'], 
         bty = 'n', pch = 20, col = 'gray')
    abline(a = 0, b = 1, col = 'red')
    
    # Estimate parameters
    my_parameters <- getIBDparameters(ped.genotypes = my_genotypes, number.cores = 1)

    # Detect IBD segments and save
    time_my_ibd <- system.time(
      my_ibd <- getIBDsegments(ped.genotypes = my_genotypes,
                               parameters = my_parameters,
                               number.cores = num_cores,
                               minimum.snps = min_snps,
                               minimum.length.bp = minimum_length_bp,
                               error = error_prob))

    results <- list(my_parameters =  my_parameters, my_ibd = my_ibd, time_my_ibd = time_my_ibd)
    save(results, file = sprintf('../pf3k_chimeric_data/IsoRelate_output_%s.RData', suffix))
  } 
}
toc()


# -------------------- isoRelate: Standalone timing experiment --------------------
tic()
for(site in sites){
  for(repetition in 1:n_timing_repetition){

    # Load data
    load(sprintf('../pf3k_chimeric_data/parents_children_%s.RData', site))
    parents_children <- pedmap
    load(sprintf('../pf3k_chimeric_data/reference_%s.RData', site))
    reference <- pedmap

    # Select a sample of the first n_timing_samples_to_run individuals
    samples_to_analyse_names <- parents_children[[1]]$iid[1:n_timing_samples_to_run]
    ind <- parents_children[[1]]$iid %in% samples_to_analyse_names
    parents_subset <- list(parents_children[[1]][ind,], parents_children[[2]])

    # Save a data set of samples for hmmIBD
    X <- fread(sprintf('../pf3k_chimeric_data/parents_children_%s.txt', site),
               header = TRUE, sep = '\t', select = c('chrom', 'pos', as.character(samples_to_analyse_names)))
    write.table(X, sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE,
                file = sprintf('../pf3k_chimeric_data/hmmIBD_timing_input_%s_%s.txt', site, repetition))

    # Reformat
    my_genotypes <- getGenotypes(ped.map = parents_subset,
                                 reference.ped.map = reference,
                                 maf = MAF_my_genotypes,
                                 isolate.max.missing = isolate_max_missing,
                                 snp.max.missing = snp_max_missing,
                                 chromosomes = NULL,
                                 input.map.distance = "cM",
                                 reference.map.distance = "cM")

    # Estimate parameters
    my_parameters <- getIBDparameters(ped.genotypes = my_genotypes, number.cores = 1)

    # Detect IBD segments and save
    time_my_ibd <- system.time(
      my_ibd <- getIBDsegments(ped.genotypes = my_genotypes,
                               parameters = my_parameters,
                               number.cores = num_cores,
                               minimum.snps = min_snps,
                               minimum.length.bp = minimum_length_bp,
                               error = error_prob))

    results <- list(my_parameters =  my_parameters, my_ibd = my_ibd, time_my_ibd = time_my_ibd)
    save(results, file = sprintf('../pf3k_chimeric_data/IsoRelate_timing_output_%s_%s.RData', site, repetition))
  }
}
toc()
 

# -------------------- hmmIBD: Standalone timing experiment --------------------
tic()
time_my_ibd <- vector('list', length(sites) * Magic_numbers$n_timing_repetition)
names(time_my_ibd) <- paste(rep(sites, each = Magic_numbers$n_timing_repetition),
                            rep(1:Magic_numbers$n_timing_repetition, length(sites)))
for(site in sites){
  for(i in 1:Magic_numbers$n_timing_repetition){
    time_my_ibd[[paste(site, i)]] <- system.time(
      system(sprintf('../HMM/hmmIBD -i ../pf3k_chimeric_data/hmmIBD_timing_input_%s_%s.txt -f ../pf3k_chimeric_data/pf3k_data_biallelic_frequencies_%s.txt -o hmmIBD_timing_output_%s_%s && mv ./hmmIBD_timing_output_* ../pf3k_chimeric_data/',
                     site, i, site, site, i, site, i)))
  }
}
save(time_my_ibd, file = '../pf3k_chimeric_data/hmmIBD_timing_output.RData')
toc()


# -------------------- hmmIBD: chimeric data analysis ----------------------
# README: Check magic numbers in hmm.c before running, change and recompile if neccessary
rm(list = ls())
tic()
load('../pf3k_chimeric_data/sites.RData')
erroneous <- c(FALSE, TRUE) 
for(err in erroneous){ # Run hmmIBD on chimeric samples 
  for(site in sites){
    suffix <- ifelse(err, paste(site, 'erroneous', sep = '_'), site)
    system(sprintf('../HMM/hmmIBD_default_rho -i ../pf3k_chimeric_data/parents_children_%s.txt -b ../pf3k_chimeric_data/b_%s.txt -f ../pf3k_chimeric_data/pf3k_data_biallelic_frequencies_%s.txt -o ../pf3k_chimeric_data/hmmIBD_%s', 
                   suffix, suffix, site, suffix))
  }  
}
toc()

# -------------------- hmmIBD: chimeric piecewise rho data analysis ---------------
# README: use the hmm with default rho 
rm(list = ls())
tic()
load('../pf3k_chimeric_data/sites.RData')
for(site in sites){
  system(sprintf('../HMM/hmmIBD_isoRelate_rho -i ../pf3k_chimeric_data/parents_children_nonuniform_%s.txt -b ../pf3k_chimeric_data/b_%s.txt -f ../pf3k_chimeric_data/pf3k_data_biallelic_frequencies_%s.txt -o ../pf3k_chimeric_data/hmmIBD_nonuniform_%s', 
                 site, site, site, site))
}  
toc()
toc()









