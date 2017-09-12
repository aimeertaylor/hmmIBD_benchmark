plot_piErr <- function(to_file) {
  iterval = c(1, 2, 3, 4, 5, 10, 20, 40)
  niter = length(iterval)

  if (to_file == 1) {
    outfile = "piErr.pdf"
    pdf(outfile)
  }

  ngen = 10
  rms_10 = numeric(niter)
  for (i in 0:(niter-1)) {
    infile = paste("results/model", i, sep="")
    infile = paste(infile, ".results.txt", sep="")
    d = read.delim(infile)
    pio = d$pi_obs[d$k_gen == ngen & d$pi_gen >= 0.3 & d$pi_gen <= 0.7]
    pig = d$pi_gen[d$k_gen == ngen & d$pi_gen >= 0.3 & d$pi_gen <= 0.7]
    dpi = pig - pio
    rms_10[i+1] = sqrt(mean(dpi^2))
  }
  print(rms_10)

  ngen = 50
  rms_50 = numeric(niter)
  for (i in 0:(niter-1)) {
    infile = paste("results/model", i, sep="")
    infile = paste(infile, ".results.txt", sep="")
    d = read.delim(infile)
    pio = d$pi_obs[d$pi_gen >= 0.3 & d$pi_gen <= 0.7 & d$k_gen == ngen]
    pig = d$pi_gen[d$pi_gen >= 0.3 & d$pi_gen <= 0.7 & d$k_gen == ngen]
    dpi = pig - pio
    rms_50[i+1] = sqrt(mean(dpi^2))
  }

  ymax = max(rms_10, rms_50)
  ymax = 1.3*ymax
#  xmax = max(iterval)
  xmax = 25

  plot(iterval,rms_10,type="b",col="blue",pch=16,ylab=expression(paste("RMS error in IBD fraction")),
     xlab="Maximum fit iterations", ylim=c(0,ymax), 
     main=paste("RMS error in IBD fraction, ", ngen, " generations", sep=""),
     xlim=c(0,xmax))
  points(iterval, rms_50, type="b",col="red",pch=16)

  legend("topright", c("N gen = 10","N gen = 50"), 
		      col=c("blue","red"), pch=rep(16,4))
  if(to_file == 1) {dev.off()}
}

plot_genErr <- function(to_file) {
  iterval = c(1, 2, 3, 4, 5, 10, 20, 40)
  niter = length(iterval)

  if (to_file == 1) {
    outfile = "genErr.pdf"
    pdf(outfile)
  }

  mingen = 10
  maxgen = 10
  rms_10 = numeric(niter)
  err_norm = numeric(niter)
  for (i in 0:(niter-1)) {
    infile = paste("results/model", i, sep="")
    infile = paste(infile, ".results.txt", sep="")
    d = read.delim(infile)
    ko = d$k_obs[d$pi_gen >= 0.3 & d$pi_gen <= 0.7 & d$k_gen<=maxgen & d$k_gen >= mingen]
    kg = d$k_gen[d$pi_gen >= 0.3 & d$pi_gen <= 0.7 & d$k_gen<=maxgen & d$k_gen >= mingen]
    dk = kg - ko
    rms_10[i+1] = sqrt(mean(dk^2))
    err_norm[i+1] = mean(dk) / mean(kg)
  }
  print(rms_10)

  mingen = 50
  maxgen = 50
  rms_50 = numeric(niter)
  for (i in 0:(niter-1)) {
    infile = paste("results/model", i, sep="")
    infile = paste(infile, ".results.txt", sep="")
    d = read.delim(infile)
    ko = d$k_obs[d$pi_gen >= 0.3 & d$pi_gen <= 0.7 & d$k_gen<=maxgen & d$k_gen >= mingen]
    kg = d$k_gen[d$pi_gen >= 0.3 & d$pi_gen <= 0.7 & d$k_gen<=maxgen & d$k_gen >= mingen]
    dk = kg - ko
    rms_50[i+1] = sqrt(mean(dk^2))
  }

  ymax = max(rms_10, rms_50)
  xmax = max(iterval)
  ymin = 0
  xmax = 25

  plot(iterval,rms_10,type="b",col="blue",pch=16,ylab=expression(paste("RMS  error on k")),
     xlab="Maximum fit iterations", ylim=c(ymin,ymax), 
     main=expression(paste("Normalized error in number of generations")),
     xlim=c(0,xmax))
  points(iterval, rms_50, type="b",col="red",pch=16)
#  points(iterval, err_norm, type="b",col="green",pch=16)

  legend("topright", c("10 gen","50 gen"
#     ,"(gen - obs)/gen"
     ),
     col=c("blue", "red"
 #    ,"green"
     ), pch=rep(16,2))
  if(to_file == 1) {dev.off()}
}

plot_pi <- function(fileno, ngen, to_file) {
# fileno controls which output file is read, and therefore what maximum fit iterations were used by HMM
# ngen specifies number of generations (k) in simulating the data

  iterval = c(1, 2, 3, 4, 5, 10, 20, 40)
  kval = c(1,3,5,10,20,50)
  if (!ngen %in% kval) {
    print("ngen must be one of")
    print(kval)
    return()
  }

  infile = paste("results/model", fileno, sep="")
  infile = paste(infile, ".results.txt", sep="")
  d = read.delim(infile)
  pio = d$pi_obs[d$k_gen==ngen]
  pig = d$pi_gen[d$k_gen==ngen]
  if (to_file == 1) {
    outfile = "pi"
    outfile = paste(outfile,"_",sep="")
    outfile = paste(outfile, ngen, sep="")
    outfile = paste(outfile, "gen.pdf",sep="")
    pdf(outfile)
  }
  plot(pig,pio,
  main=paste("Fraction IBD, ", ngen," gen. Max fit iter: ", iterval[fileno+1], sep=""),
  xlab="Generated IBD fraction",ylab="Reconstructed IBD fraction",pch=16, 
  cex=.5,col="darkblue")
  points(c(0,1),c(0,1),col="black",type="l")
  if(to_file == 1) {dev.off()}
}

plot_piErrXP <- function(to_file) {
  iterval = c(1, 2, 3, 4, 5, 10, 20, 40)
  niter = length(iterval)
  if (to_file == 1) {
    outfile = "piErrXP.pdf"
    pdf(outfile)
  }

  ngen = 10
  rms_10 = numeric(niter)
  for (i in 0:(niter-1)) {
    infile = paste("results/xp", i, sep="")
    infile = paste(infile, ".results.txt", sep="")
    d = read.delim(infile)
    pio = d$pi_obs[d$k_gen == ngen]
    pig = d$pi_gen[d$k_gen == ngen]
    dpi = pig - pio
    rms_10[i+1] = sqrt(mean(dpi^2))
  }

  ngen = 50
  rms_50 = numeric(niter)
  for (i in 0:(niter-1)) {
    infile = paste("results/xp", i, sep="")
    infile = paste(infile, ".results.txt", sep="")
    d = read.delim(infile)
    pio = d$pi_obs[d$k_gen == ngen]
    pig = d$pi_gen[d$k_gen == ngen]
    dpi = pig - pio
    rms_50[i+1] = sqrt(mean(dpi^2))
  }

  ymax = max(rms_10, rms_50)
  ymax = 1.3*ymax
  xmax = max(iterval)
  xmax = 25

  plot(iterval,rms_10,type="b",col="blue",pch=16,ylab=expression(paste("RMS error in IBD fraction")),
     xlab="Maximum fit iterations", ylim=c(0,ymax), 
     main=paste("RMS error in IBD fraction, ", ngen, " generations", sep=""),
     xlim=c(0,xmax))
  points(iterval, rms_50, type="b",col="red",pch=16)

  legend("topright", c("N gen = 10","N gen = 50"), 
		      col=c("blue","red"), pch=rep(16,4))
  if(to_file == 1) {dev.off()}
}

plot_genErrXP <- function(to_file) {
  iterval = c(1, 2, 3, 4, 5, 10, 20, 40)
  niter = length(iterval)

  if (to_file == 1) {
    outfile = "genErrXP.pdf"
    pdf(outfile)
  }

  mingen = 10
  maxgen = 10
  rms_10 = numeric(niter)
  err_norm = numeric(niter)
  for (i in 0:(niter-1)) {
    infile = paste("results/xp", i, sep="")
    infile = paste(infile, ".results.txt", sep="")
    d = read.delim(infile)
    ko = d$k_obs[d$pi_gen >= 0.3 & d$pi_gen <= 0.7 & d$k_gen<=maxgen & d$k_gen >= mingen]
    kg = d$k_gen[d$pi_gen >= 0.3 & d$pi_gen <= 0.7 & d$k_gen<=maxgen & d$k_gen >= mingen]
    dk = kg - ko
    rms_10[i+1] = sqrt(mean(dk^2))
    err_norm[i+1] = mean(dk) / mean(kg)
  }
  print(rms_10)
#  print(err_norm)

  mingen = 50
  maxgen = 50
  rms_50 = numeric(niter)
  for (i in 0:(niter-1)) {
    infile = paste("results/xp", i, sep="")
    infile = paste(infile, ".results.txt", sep="")
    d = read.delim(infile)
    ko = d$k_obs[d$pi_gen >= 0.3 & d$pi_gen <= 0.7 & d$k_gen<=maxgen & d$k_gen >= mingen]
    kg = d$k_gen[d$pi_gen >= 0.3 & d$pi_gen <= 0.7 & d$k_gen<=maxgen & d$k_gen >= mingen]
    dk = kg - ko
    rms_50[i+1] = sqrt(mean(dk^2))
  }

  ymax = max(rms_10, rms_50)
  xmax = max(iterval)
  ymin = 0
  xmax = 25
  plot(iterval,rms_10,type="b",col="blue",pch=16,ylab=expression(paste("RMS  error on k")),
     xlab="Maximum fit iterations", ylim=c(ymin,ymax), 
     main=expression(paste("Normalized error in number of generations")),
     xlim=c(0,xmax))
  points(iterval, rms_50, type="b",col="red",pch=16)
#  points(iterval, err_norm, type="b",col="green",pch=16)

  legend("topright", c("10 gen","50 gen"), 
#     "(gen - obs)/gen"), 
		      col=c("blue", "red"), pch=rep(16,2))
  if(to_file == 1) {dev.off()}
}

plot_piErrT <- function(to_file) {
  thinval = c(1, 2, 4, 10, 20, 100, 500)
  nmarker = c(132610, 73231, 37190, 14880, 7440, 1488, 298)

  if (to_file == 1) {
    outfile = "piErr_thin.pdf"
    pdf(outfile)
  }

  ngen = 10
  nthin = length(thinval)
  rms = numeric(nthin)
  for (i in 0:(nthin-1)) {
    infile = paste("results_thinned/model", i, sep="")
    infile = paste(infile, ".results.txt", sep="")
    d = read.delim(infile)
    pio = d$pi_obs[d$k_gen == ngen & d$pi_gen >= 0.3 & d$pi_gen <= 0.7]
    pig = d$pi_gen[d$k_gen == ngen & d$pi_gen >= 0.3 & d$pi_gen <= 0.7]
    dpi = pig - pio
    rms[i+1] = sqrt(mean(dpi^2))
  }

  ymax = max(rms)
  ymax = .2
  xmax = max(nmarker)

  plot(nmarker,rms,type="b",col="blue",pch=16,ylab=expression(paste("RMS error in IBD fraction")),
     xlab="Number of variant sites", ylim=c(0,ymax),
     main="RMS error in IBD fraction",
     xlim=c(0,xmax))

  if(to_file == 1) {dev.off()}
}

plot_piT <- function(fileno, to_file) {
  thinval = c(1, 2, 4, 10, 20, 100, 500)
  nmarker = c(132610, 73231, 37190, 14880, 7440, 1488, 298)

  ngen = 10
  infile = paste("results_thinned/model", fileno, sep="")
  infile = paste(infile, ".results.txt", sep="")
  print(infile)
  d = read.delim(infile)
  pio = d$pi_obs[d$k_gen==ngen]
  pig = d$pi_gen[d$k_gen==ngen]
  if (to_file == 1) {
    outfile = "pi_thin"
    outfile = paste(outfile, thinval[fileno+1], sep="")
    outfile = paste(outfile, ".pdf",sep="")
    pdf(outfile)
  }
  plot(pig,pio,
   main=paste("Fraction IBD, N variants: ", nmarker[fileno+1], sep=""),
  xlab=expression(Generated~IBD~fraction),ylab=expression(Reconstructed~IBD~fraction),pch=16, 
  cex=.5,col="darkblue")
  points(c(0,1),c(0,1),col="black",type="l")
  if(to_file == 1) {dev.off()}
}

plot_piErrS <- function(to_file) {
  smudgeval = c(1000, 400, 200, 100, 50, 25, 10)

  if (to_file == 1) {
    outfile = "piErr_smudged.pdf"
    pdf(outfile)
  }

  nsmudge = length(smudgeval)
  rms = numeric(nsmudge)
  for (i in 0:(nsmudge-1)) {
    infile = paste("results_smudged/model", i, sep="")
    infile = paste(infile, ".results.txt", sep="")
    d = read.delim(infile)
    pio = d$pi_obs[d$pi_gen >= 0.3 & d$pi_gen <= 0.7]
    pig = d$pi_gen[d$pi_gen >= 0.3 & d$pi_gen <= 0.7]
    dpi = pig - pio
    rms[i+1] = sqrt(mean(dpi^2))
  }

  ymax = max(rms)
  ymin = min(rms)
#  ymax = .2
  xmax = max(smudgeval)
  xmin = .8 * min(smudgeval)


  plot(smudgeval,rms,type="b",col="blue",pch=16,ylab=expression(paste("RMS error in IBD fraction")),
     xlab="Sample size", ylim=c(ymin,ymax), log="xy",
     main="Effect of error allele frequency on RMS error in IBD fraction",
     xlim=c(xmin,xmax))
  if(to_file == 1) {dev.off()}
}

