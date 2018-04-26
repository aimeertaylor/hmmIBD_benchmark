#######################################################################
# Function to recombine two parent genomes into a chimeric child      #
# using a piecewise constant function based on Miles A, et al. (2016) #                             #
#######################################################################
recombine_piecewise <- function(trio, nchrom, rho_av, rho_low, rho_high){
  
  parents <- colnames(trio)[grep('parent', colnames(trio))]
  
  for(chrom in 1:nchrom){
    
    # Extract length of chromosome and centromere start and end points
    chrom_length <- chrom_lengths[chrom,'V3']
    centromere <- centromeres[chrom,c('V2','V3')]
    
    # Change points for piecewise constant function 
    c1 <- min(centromere)-120000
    c2 <- min(centromere)-80000
    c3 <- min(centromere)-30000
    c4 <- min(centromere)
    c5 <- max(centromere)
    c6 <- max(centromere)+30000
    c7 <- max(centromere)+80000
    c8 <- max(centromere)+120000
    
    # Create vector of rhos for entire chrom (inc. non polymorphic sites)
    rhos <- rep(rho_av, chrom_length) 
    rhos[c(c1:c2, c7:c8)] <- rho_high
    rhos[c(c3:c4, c5:c6)] <- rho_low
    
    # Sample crossover events with prob according to rhos
    cross_overs <- runif(chrom_length) <= rhos
    
    # Choose a parent to start with probability 0.5
    parent_first <- sample(parents, size = 1) 
    parent_second <- parents[which(parents != parent_first)] 
    ind_parent <- as.logical(cumsum(cross_overs)%%2)  # TRUE if odd, FALSE if even
    
    # Populate chrom 
    ind_chrom <- trio$chrom == chrom  # Extract chrom ind
    chrom_pos <- trio[ind_chrom, 'pos'] # Extract pos of polymorphic on chrom
    ind_parent_pos <- ind_parent[chrom_pos] # Extract parent ind for polymorphic only
    
    trio[ind_chrom, 'child'][ind_parent_pos] <- trio[ind_chrom, parent_first][ind_parent_pos]
    trio[ind_chrom, 'child'][!ind_parent_pos] <- trio[ind_chrom, parent_second][!ind_parent_pos]
    trio[ind_chrom, 'assignment'][ind_parent_pos] <- parent_first
    trio[ind_chrom, 'assignment'][!ind_parent_pos] <- parent_second
  }
  return(trio)
}



#######################################################################
# Function to recombine two parent genomes into a chimeric child     #                                  
#######################################################################
recombine_uniform <- function(trio, nchrom, rho){
  
  parents <- colnames(trio)[grep('parent', colnames(trio))]
  
  for(chrom in 1:nchrom){
    
    parent <- sample(parents, size = 1) # choose a parent with probability 0.5
    chrom_length <- chrom_lengths[chrom,'V3']
    time <- rexp(1, rho)
    
    if(time >= chrom_length){ # no recombination 
      
      ind <- (trio$chrom == chrom) & (trio$pos <= time)
      trio[ind, 'child'] <- trio[ind, parent]
      trio[ind, 'assignment'] <- rep(parent, sum(ind))
      
    } else { # recombine
      
      while(time < chrom_length){
        
        # populate upto crossover time 
        ind <- (trio$chrom == chrom) & (trio$pos <= time)
        trio[ind, 'child'] <- trio[ind, parent]
        trio[ind, 'assignment'] <- rep(parent, sum(ind))
        
        # switch parent and redraw a new crossover time
        parent <- parents[which(parents != parent)] 
        time <- time + rexp(1, rho)
        
      } 
      
      # fill in remaining
      ind <- (trio$chrom == chrom) & is.na(trio$child)
      trio[ind, 'child'] <- trio[ind, parent]
      trio[ind, 'assignment'] <- rep(parent, sum(ind))
    }
  }
  return(trio)
}



#######################################################################
# Function to format data for isolate                                 #
#######################################################################
reformat_isorelate <- function(Data, rho){
  
  # Translate genotype set from {-1,0,1} to {0,1,2}    
  Y <- as.matrix(Data[,-(1:2)]) + 1  
  
  # Duplicate each row (since isoRelate expects two adjacent columns per SNP)
  SNPData <- array(dim = c(ncol(Y), nrow(Y)*2))
  inds <- cbind(seq(1, ncol(SNPData), by = 2), seq(2, ncol(SNPData), by = 2)) 
  for(i in 1:nrow(Y)){
    SNPData[,inds[i,1]] <- Y[i,]
    SNPData[,inds[i,2]] <- Y[i,]
  }
  
  X <- data.frame(fid = site, 
                  iid = colnames(Data)[-(1:2)], 
                  pid = 0, mid = 0, 
                  moi = 1, aff = 2,
                  SNPData)
  
  Z <- data.frame(chr = Data[,'chrom'], 
                  snp_id = paste(Data[,'chrom'], ":", Data[,'pos'], sep = ''), 
                  pos_cM = Data[,'pos'] * rho * 100, # divide by 100 since centi Morgans, not Morgans 
                  pos_bp = Data[,'pos'])  
  
  return(list(X, Z))
}

