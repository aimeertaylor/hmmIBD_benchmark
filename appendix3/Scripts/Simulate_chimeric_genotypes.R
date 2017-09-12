################################################################################
# Script to 
# 1) extract biallelic Pf3k data per site for sites with more than 100 isolates; 
# 2) calculate, plot and save pairwise IBS per site; 
# 3) extract unrelated sample pair names per site and save; 
# 4) make chimeric genotypes out of unrelated sample pairs and plot crossovers;
################################################################################

# Clear workspace and load libraries/functions
rm(list = ls()) 

# Packages and ancillary files 
inst_pkg_list <- rownames(installed.packages())
if("data.table" %in% inst_pkg_list){
  library(data.table) # for fread (allows one to read a subset of column names)  
} else {
  install.packages("data.table")
  library(data.table)
}
source(file = './functions.R') # functions

# Magic numbers
Magic_numbers <- list(min_num = 100, # minimum number of isolates per site
                      MAF = 0.01, # minimum minor allele frequency  
                      cutoff = 0.01, # the identity-by-state percentile cut off
                      nchrom = 14, # number of Plasmodium falciparum chromosomes 
                      rho = 5.833965e-07, # this is based on a comparison of positions in isoRelates png_pedmap
                      seed_simulate = 1) # recombination rate per base pair 

attach(Magic_numbers)
save(Magic_numbers, file = '../pf3k_chimeric_data/Simulate_chimeric_genotypes_Magic_numbers.RData')
set.seed(seed_simulate)

# Load data and save sites 
sample_meta_data <- read.delim('../pf3k_data/monogenomic_samples_touse.txt', check.names=FALSE)
sites <- names(which(table(sample_meta_data$site) > min_num)) # select sites with more than min_num isolates
save(sites, file = '../pf3k_chimeric_data/sites.RData')

#  --------------------- Extract isolates per site --------------------------
for(site in sites){
  
  # Extract isolate names
  site_ind <- sample_meta_data$site == site
  sample_names_site <- as.character(sample_meta_data$isolate[site_ind])
  
  # Load data per site
  site_data <- NULL # Set to NUll since not sure of array dim ahead of data extraction
  for(chrom in 1:nchrom){
    site_data <- rbind(site_data, 
                       fread(sprintf('../pf3k_data/pf3k_%s_seq.txt', chrom), 
                             select = c('chrom', 'pos', sample_names_site), sep = '\t', header = TRUE))
  }
  
  # Remove multiallelic SNPs (not supported by isoRelate) 
  biallelic <- apply(site_data[,-(1:2)], 1, FUN = function(x){max(x) <= 1})
  site_data_biallelic <- site_data[biallelic, ]
  
  # Remove SNPs with frequency <= MAF
  temp <-  site_data_biallelic # Replace missing with NA
  temp[temp == -1] <- NA
  allele_frequencies <- apply(temp[,-(1:2)], 1, mean, na.rm = TRUE)
  MAFplus <- allele_frequencies > MAF
  site_data_biallelic_MAFplus <- as.matrix(site_data_biallelic[MAFplus, ])
  site_data_biallelic_MAFplus_frequencies <- cbind(site_data_biallelic_MAFplus[,c('chrom', 'pos')], 
                                                   1-allele_frequencies[MAFplus], 
                                                   allele_frequencies[MAFplus])
  print(sprintf('%s: %s SNPs', site, nrow(site_data_biallelic_MAFplus_frequencies)))
  
  # Reformat to avoid scientific notation
  site_data_biallelic_MAFplus_frequencies[,'pos'] <- format(site_data_biallelic_MAFplus_frequencies[,'pos'], scientific = FALSE)
  
  # Save site specific biallelic data and allele frequencies
  save(site_data_biallelic_MAFplus,
       file = sprintf('../pf3k_chimeric_data/pf3k_data_biallelic_%s.RData', site))
  
  # Save site specific allele frequencies for hmmIBD
  write.table(site_data_biallelic_MAFplus_frequencies, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t', # to prevent 8e+05
              file = sprintf('../pf3k_chimeric_data/pf3k_data_biallelic_frequencies_%s.txt', site))
  
  # Save site reference pedmap for isoRelate
  pedmap <- reformat_isorelate(site_data_biallelic_MAFplus, rho)
  save(pedmap, file = sprintf('../pf3k_chimeric_data/reference_%s.RData',site))
}




# -------------------- Calculate IBS dist per site --------------------
for(site in sites){
  
  # Reload data
  load(sprintf('../pf3k_chimeric_data/pf3k_data_biallelic_%s.RData', site))
  
  # Replace missing with NA
  site_data_biallelic_MAFplus[site_data_biallelic_MAFplus == -1] <- NA
  
  # Calculate IBS for each sample pair per site using manhatten. Note that, 
  # "If some columns are excluded in calculating... Manhattan, ...,
  # the sum is scaled up proportionally to the number of columns used."
  n_snps <- nrow(site_data_biallelic_MAFplus)
  pairwise_IBS <- (n_snps - dist(t(site_data_biallelic_MAFplus[,-(1:2)]), method = 'manhattan'))/n_snps
  
  # Save IBS results
  save(pairwise_IBS, file = sprintf('../pf3k_chimeric_data/pf3k_data_IBS_%s.RData', site))
}  



# --------------------  Plot IBS dist per site -------------------- 
par(mfrow = c(length(sites),1))
for(site in sites){
  
  # Reload IBS to plot
  load(sprintf('../pf3k_chimeric_data/pf3k_data_IBS_%s.RData', site))
  
  # Plot IBS distribution per site with cutoff
  hist(pairwise_IBS, col = adjustcolor('blue', alpha.f = 0.5), freq = FALSE, main = site, 
       xlim = c(0.5,1), xaxt = 'n',  breaks = 50)
  axis(side = 1, at = c(0.5, 0.75, 1), cex.axis = 0.7)
  abline(v = quantile(pairwise_IBS, probs = cutoff), col = 'red')
}


# -------------------- Extract unrelated sample pair names -------------------- 
for(site in sites){
  
  # Load IBS results
  load(sprintf('../pf3k_chimeric_data/pf3k_data_IBS_%s.RData', site))
  
  # Reformat dist
  sample_names <- attr(pairwise_IBS, which = 'Labels')
  pairwise_IBS <- as.vector(pairwise_IBS) # c(column 1, column 2, column 3... )
  
  # Name by sample comparison
  n_samples <- length(sample_names)
  sample_j <- rep(sample_names, (n_samples:1)-1)
  sample_i <- c()
  for(i in 2:n_samples){
    sample_i <- c(sample_i, sample_names[i:n_samples])
  }
  names(pairwise_IBS) <- paste(sample_i, sample_j, sep = "_")
  
  # Extract unrelated
  unrelated <- names(pairwise_IBS[pairwise_IBS < quantile(pairwise_IBS, probs = cutoff)])
  unrelated <- do.call(rbind, strsplit(unrelated, split = "_"))
  
  print(sprintf('%s: %s unrelated sample comparisons, %s unique parents', 
                site, nrow(unrelated), length(unique(as.vector(unrelated)))))

  # Save
  save(unrelated, file = sprintf('../pf3k_chimeric_data/pf3k_data_unrelated_%s.RData', site)) 
}



# -------------- Make chimeric genotypes out of unrelated sample pairs  --------------------
# Load chromosome lengths in bp
chrom_lengths <- read.csv(file = '../pf3k_data/Pf_v3_chrom_length.csv', 
                          skip = 2, header = FALSE)

for(site in sites){
  
  # Load data and unrelated names
  load(sprintf('../pf3k_chimeric_data/pf3k_data_biallelic_%s.RData', site))
  load(sprintf('../pf3k_chimeric_data/pf3k_data_unrelated_%s.RData', site))
  parents <- unique(as.vector(unrelated))
  children <- apply(unrelated, 1, paste, collapse = ":") 
  n_parents <- length(parents)
  n_children <- length(children)
  n_snps <- nrow(site_data_biallelic_MAFplus)
  
  # Children
  children_store <- array(dim = c(n_snps, n_children), 
                          dimnames = list(NULL, children))
  
  # Store to save assignment of chimeric child
  parent_assignment <- data.frame(chrom = site_data_biallelic_MAFplus[,'chrom'], 
                                  pos = site_data_biallelic_MAFplus[,'pos'], 
                                  array(dim = c(n_snps, n_children), 
                                        dimnames = list(NULL, children)), 
                                  check.names = FALSE) # otherwise ":" being converted to "."
  
  # For each unrelated sample pair generate a chimeric child
  for(i in 1:n_children){ 
    
    trio <- data.frame(assignment = NA, 
                       chrom = site_data_biallelic_MAFplus[,'chrom'], 
                       pos = site_data_biallelic_MAFplus[,'pos'], 
                       parent1 = site_data_biallelic_MAFplus[, unrelated[i,1]], 
                       parent2 = site_data_biallelic_MAFplus[, unrelated[i,2]], 
                       child = NA)
    
    # create chimeric
    result <- recombine(trio, rho, nchrom)
    
    # Record parent assignment
    parent_assignment[, children[i]] <- result[,1]
    children_store[, children[i]] <- result[,'child']
  }
  
  # Collate genomes
  parents_children <- cbind(site_data_biallelic_MAFplus[,c('chrom', 'pos', parents)], children_store)
  
  # Save parent assignment
  save(parent_assignment, file = sprintf('../pf3k_chimeric_data/parent_assignment_%s.RData',site))
  
  # Save data set for hmmIBD
  write.table(parents_children, sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE, 
              file = sprintf('../pf3k_chimeric_data/parents_children_%s.txt', 
                             site))
  
  # Save all parent children in pedmap for isoRelate
  pedmap <- reformat_isorelate(parents_children, rho)
  save(pedmap, file = sprintf('../pf3k_chimeric_data/parents_children_%s.RData',
                              site))
}




# -------------- Plot average number of cross-overs per chromosome --------------
par(mfrow = c(length(sites),1))
for(site in sites){
  
  # load chimeric data 
  load(sprintf('../pf3k_chimeric_data/parent_assignment_%s.RData', site))
  unrelated_pairs <- colnames(parent_assignment[,-(1:2)])
  crossovers <- array(dim = c(length(unrelated_pairs), 14), dimnames = list(unrelated_pairs, as.character(1:14)))
  
  # Extract crossovers from parent assignment
  for(i in unrelated_pairs){
    X <- table(unique(parent_assignment[,c('chrom',i)])[,1])
    crossovers[i, names(X)] <- X
  }
  
  # Plot
  barplot(colMeans(crossovers), las = 2, ylim = c(0,2), main = site, 
          xlab = 'Chromosome', ylab = 'Average no. of crossovers')
}