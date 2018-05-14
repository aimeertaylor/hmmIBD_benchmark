############################################################################
# Script to extract and post-process both isoRelate and hmmIBD results
# Takes approx. 400 sec to run 
############################################################################
rm(list = ls())
require(tictoc)
tic()
load('../pf3k_chimeric_data/sites.RData')
erroneous <- c(FALSE, TRUE) 

for(err in erroneous){
  for(site in sites){
    
    suffix <- ifelse(err, paste(site, 'erroneous', sep = '_'), site)
    
    # Load assignments for comparisons of interest
    load(sprintf('../pf3k_chimeric_data/parent_assignment_%s.RData', site))
    n_snps <- nrow(parent_assignment)
    
    # Load results from isoRelate and hmmIBD 
    load(sprintf('../pf3k_chimeric_data/IsoRelate_output_%s.RData', suffix))
    hmmIBD <- read.delim(sprintf('../pf3k_chimeric_data/hmmIBD_%s.hmm.txt', suffix)) 
    
    # List comparisons of interest 
    sample_names <- unique(c(results$my_ibd$iid1, results$my_ibd$iid2))
    children_names <- sample_names[grep(':', sample_names)]
    comparisons_of_interest <- paste(child = rep(children_names,2), as.vector(do.call(rbind,(strsplit(children_names, split = ":")))), sep = ' ')
    n_comparisons <- length(comparisons_of_interest)
    
    # Create stores in which to extract results for comparison with parent assignment
    isoRelate_assignment <- cbind(parent_assignment[, c('chrom', 'pos')], 
                                  array(dim = c(n_snps, n_comparisons), dimnames = list(NULL, comparisons_of_interest)))
    isoRelate_assignment[,-(1:2)] <- NA
    hmmIBD_assignment <- isoRelate_assignment
    isoRelate_assignment[is.na(isoRelate_assignment)] <- 'DBD'
    
    # Create comparison names for isoRelate and hmmIBD for look up with comparisons of interest
    isoRelate_comp_12 <- apply(results$my_ibd[, c('iid1', 'iid2')], 1, paste, collapse = ' ')
    isoRelate_comp_21 <- apply(results$my_ibd[, c('iid2', 'iid1')], 1, paste, collapse = ' ')
    hmmIBD_comp_12 <- apply(hmmIBD[, c('sample1', 'sample2')], 1, paste, collapse = ' ')
    hmmIBD_comp_21 <- apply(hmmIBD[, c('sample2', 'sample1')], 1, paste, collapse = ' ')
    
    # For each comparison of interest, extract inferred segments 
    Results <- vector('list', length = n_comparisons)
    names(Results) <- comparisons_of_interest
    for(i in 1:n_comparisons){ 
      X <- results$my_ibd[isoRelate_comp_21 %in% comparisons_of_interest[i] | isoRelate_comp_12 %in% comparisons_of_interest[i],]
      Y <- hmmIBD[hmmIBD_comp_12 %in% comparisons_of_interest[i] | hmmIBD_comp_21 %in% comparisons_of_interest[i], ]
      if(nrow(X)== 0){print(sprintf('no results for isorelate %s %s', i, site))} 
      if(nrow(Y)== 0){print(sprintf('no results for hmmIBD %s %s', i, site))} 
      Results[[i]] <- list(isoRelate = X, hmmIBD = Y)
    }
    
    # Make memory space
    rm(list = c('hmmIBD', 'hmmIBD_comp_12', 'hmmIBD_comp_21', 'isoRelate_comp_12', 'isoRelate_comp_21'))
    
    # For each comparison of interest, reformat inferred segments
    for(i in 1:n_comparisons){ 
      
      X <- Results[[i]]$isoRelate
      Y <- Results[[i]]$hmmIBD
      
      samples <- strsplit(comparisons_of_interest[i], split = ' ')[[1]]
      ind <- grepl(':', samples)
      child <- samples[ind]
      parent <- sprintf('parent%s', which(strsplit(as.character(child), split = ":")[[1]] == samples[!ind]))
      
      if(nrow(X) > 0){
        for(j in 1:nrow(X)){ # Populate isoRelate_assignment$chrom
          ind_chr <- isoRelate_assignment$chrom == X[j,'chr']
          ind_pos <- (isoRelate_assignment$pos >= X[j,'start_position_bp']) & 
            (isoRelate_assignment$pos <= X[j,'end_position_bp']) 
          isoRelate_assignment[ind_chr & ind_pos, comparisons_of_interest[i]] <- parent
        }}
      
      if(nrow(Y) > 0){
        for(j in 1:nrow(Y)){ # Populate hmmIBD_assignment$chrom
          
          # If IBD
          if(Y[j,'different'] == 0){ 
            ind_chr <- hmmIBD_assignment$chrom == Y[j,'chr']
            ind_pos <- (hmmIBD_assignment$pos >= Y[j,'start']) & (hmmIBD_assignment$pos <= Y[j,'end']) 
            hmmIBD_assignment[ind_chr & ind_pos, comparisons_of_interest[i]] <- parent  
          } 
          
          # If not IBD
          if(Y[j,'different'] == 1){ 
            ind_chr <- hmmIBD_assignment$chrom == Y[j,'chr']
            ind_pos <- (hmmIBD_assignment$pos >= Y[j,'start']) & (hmmIBD_assignment$pos <= Y[j,'end']) 
            hmmIBD_assignment[ind_chr & ind_pos, comparisons_of_interest[i]] <- 'DBD'  
          } 
        }
      }
      
    }
    
    # Not strictly necessary to save different comparisons_of_interest for non-erroneous and erroneous
    # but do so anyway in case of future discrepancies that may arise upon modifications of the code
    save(Results, file = sprintf('../pf3k_chimeric_data/Results_comparisons_of_interest_%s.RData', suffix))
    save(comparisons_of_interest, file = sprintf('../pf3k_chimeric_data/comparisons_of_interest_%s.RData', suffix))
    save(isoRelate_assignment, file = sprintf('../pf3k_chimeric_data/isoRelate_assignment_%s.RData', suffix))
    save(hmmIBD_assignment, file = sprintf('../pf3k_chimeric_data/hmmIBD_assignment_%s.RData', suffix))
  }
}
toc()


############################################################################
# Addition to postprocess results generated under hmmIBD given nonuniform  #
# recombination rates.                                                      
############################################################################
rm(list = ls())
require(tictoc)
tic()
load('../pf3k_chimeric_data/sites.RData')

for(site in sites){
  
  # Load assignments for comparisons of interest
  load(sprintf('../pf3k_chimeric_data/parent_assignment_nonuniform_%s.RData', site))
  n_snps <- nrow(parent_assignment_nonuniform)
  
  # Load results from isoRelate and hmmIBD 
  hmmIBD <- read.delim(sprintf('../pf3k_chimeric_data/hmmIBD_nonuniform_%s.hmm.txt', site)) 
  
  # List comparisons of interest 
  sample_names <- unique(c(as.character(hmmIBD$sample1), as.character(hmmIBD$sample2)))
  children_names <- sample_names[grep(':', sample_names)]
  comparisons_of_interest <- paste(child = rep(children_names,2), as.vector(do.call(rbind,(strsplit(children_names, split = ":")))), sep = ' ')
  n_comparisons <- length(comparisons_of_interest)
  
  # Create stores in which to extract results for comparison with parent assignment
  hmmIBD_assignment <- cbind(parent_assignment_nonuniform[, c('chrom', 'pos')], 
                             array(dim = c(n_snps, n_comparisons), dimnames = list(NULL, comparisons_of_interest)))
  hmmIBD_assignment[,-(1:2)] <- NA
  
  # Create comparison names for hmmIBD for look up with comparisons of interest
  hmmIBD_comp_12 <- apply(hmmIBD[, c('sample1', 'sample2')], 1, paste, collapse = ' ')
  hmmIBD_comp_21 <- apply(hmmIBD[, c('sample2', 'sample1')], 1, paste, collapse = ' ')
  
  # For each comparison of interest, extract inferred segments 
  Results <- vector('list', length = n_comparisons)
  names(Results) <- comparisons_of_interest
  for(i in 1:n_comparisons){ 
    Y <- hmmIBD[hmmIBD_comp_12 %in% comparisons_of_interest[i] | hmmIBD_comp_21 %in% comparisons_of_interest[i], ]
    if(nrow(Y) == 0){print(sprintf('no results for hmmIBD %s %s', i, site))} 
    Results[[i]] <- list(hmmIBD = Y)
  }
  
  # Make memory space
  rm(list = c('hmmIBD', 'hmmIBD_comp_12', 'hmmIBD_comp_21'))
  
  # For each comparison of interest, reformat inferred segments
  for(i in 1:n_comparisons){ 
    
    Y <- Results[[i]]$hmmIBD
    
    samples <- strsplit(comparisons_of_interest[i], split = ' ')[[1]]
    ind <- grepl(':', samples)
    child <- samples[ind]
    parent <- sprintf('parent%s', which(strsplit(as.character(child), split = ":")[[1]] == samples[!ind]))
    
    if(nrow(Y) > 0){
      for(j in 1:nrow(Y)){ # Populate hmmIBD_assignment$chrom
        
        # If IBD
        if(Y[j,'different'] == 0){ 
          ind_chr <- hmmIBD_assignment$chrom == Y[j,'chr']
          ind_pos <- (hmmIBD_assignment$pos >= Y[j,'start']) & (hmmIBD_assignment$pos <= Y[j,'end']) 
          hmmIBD_assignment[ind_chr & ind_pos, comparisons_of_interest[i]] <- parent  
        } 
        
        # If DBD
        if(Y[j,'different'] == 1){ 
          ind_chr <- hmmIBD_assignment$chrom == Y[j,'chr']
          ind_pos <- (hmmIBD_assignment$pos >= Y[j,'start']) & (hmmIBD_assignment$pos <= Y[j,'end']) 
          hmmIBD_assignment[ind_chr & ind_pos, comparisons_of_interest[i]] <- 'DBD'  
        } 
      }
    }
  }
  
  # Not strictly necessary to save different comparisons_of_interest for non-erroneous and erroneous
  # but do so anyway in case of future discrepancies that may arise upon modifications of the code
  save(Results, file = sprintf('../pf3k_chimeric_data/Results_comparisons_of_interest_nonuniform_%s.RData', site))
  save(comparisons_of_interest, file = sprintf('../pf3k_chimeric_data/comparisons_of_interest_nonuniform_%s.RData', site))
  save(hmmIBD_assignment, file = sprintf('../pf3k_chimeric_data/hmmIBD_assignment_nonuniform_%s.RData', site))
}





