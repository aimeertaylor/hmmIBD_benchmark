#######################################################################
# Functions to recombine two parent genomes into a chimeric child     #                                  
#######################################################################
recombine <- function(trio, rho, nchrom){
  
  parents <- colnames(trio)[grep('parent', colnames(trio))]
  parent <- sample(parents, size = 1) # choose a parent with probability 0.5
  for(chrom in 1:nchrom){
    
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
# June 2017                                                           #
# Aimee Taylor                                                        #
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
                  pos_cM = Data[,'pos'] * rho*100, # divide by 100 since centi Morgans, not Morgans 
                  pos_bp = Data[,'pos'])  
  
  return(list(X, Z))
}

