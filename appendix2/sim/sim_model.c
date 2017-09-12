#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "random.h"

// Create simulated emission data purely from the model. Frequencies and spacing are 
//  drawn from real data. States and transitions are assigned from model parameters, and outputs 
//  generated from emission probabilities. 

int main(int argc, char **argv) {
  // Input parameters
  double target_pi[2];       // [0] = fraction IBD
  double k_rec;        // recombination scaling
  double eps;          // genotype error rate
  
  int max_snp = 30000;
  int linesize = 4000;
  const int nchrom = 14;            // 14 for falciparum
  //const double rec_rate = 7.4e-7; // 7.4e-5 cM/bp or 13.5 kb/cM Miles et al, Genome Res 26:1288-1299 (2016)
  const double rec_rate = 5.8e-7;   // 5.8e-5 cM/bp, or 17kb/cM
  const int niter = 100;

  FILE *outf=NULL, *sumf=NULL, *ff=NULL, *statef=NULL;
  int pflag, kflag, iflag, c, itoken, **geno, chr, *pos, nsnp, **sampstate=NULL;
  int *start_chr=NULL, *end_chr=NULL, state, ksamp;
  int isnp, ntrans, dist, nstate[2], eflag, fflag, iter;
  char outfile[256], *newLine=NULL, *token, *running;
  char sumfile[256], statefile[256];
  double x, ptrans, pi, *ffreq=NULL, pr, p11, p12, p22;
  
  newLine = malloc((linesize+1) * sizeof(char));
  start_chr = malloc((nchrom+1) * sizeof(int));
  end_chr = calloc(nchrom+1, sizeof(int));
  for (chr = 1; chr <= nchrom; chr++) {start_chr[chr] = 10000000;}
  pflag = kflag = iflag = fflag = 0;
  seed_rng();   
  //  set_rng_seed(111111);
  while ( (c = getopt(argc, argv, ":p:k:i:e:f:")) != -1) {
    switch(c) {
    case 'p':
      target_pi[0] = strtod(optarg, NULL);
      target_pi[1] = 1 - target_pi[0];
      fprintf(stdout, "pi set to %f\n", target_pi[0]);
      if (target_pi[0] > 1 || target_pi[0] < 0) {
	fprintf(stderr, "Illegal value for target_pi: %.4f\n", target_pi[0]);
	abort();
      }
      pflag = 1;
      break;
    case 'k':
      k_rec = strtod(optarg, NULL);
      fprintf(stdout, "krec set to %f\n", k_rec);
      kflag = 1;
      break;
    case 'e':
      eps = strtod(optarg, NULL);
      eflag = 1;
      break;
    case 'f':
      ff = fopen(optarg, "r");
      if (ff == NULL) {fprintf(stderr, "could not open frequency file %s\n", optarg); exit(0);}
      fflag = 1;
      break;
    case ':':
      fprintf(stdout, "argument %c requires an argument\n", optopt);
      abort();
    }
  }
  if (optind != argc || kflag == 0 || pflag == 0 || eflag == 0  || fflag == 0) {
    fprintf(stderr, "Usage: sim_model -k <recomb coeff> -p <fraction IBD> ");
    fprintf(stderr, "-e <genotype err rate> -f <allele freq file>\n");
    exit(0);
  }
  //  sprintf(outfile, "sims_%.1f_%.0f.txt", k_rec, 100*target_pi[0]);
  sprintf(outfile, "sim_model.txt");
  sprintf(sumfile, "simsum_model.txt");
  sprintf(statefile, "sim_state.txt");
  pr = 1 - eps;
    
  // Read frequencies from file
  nsnp = 0;
  pos = malloc(max_snp * sizeof(int));
  ffreq = malloc(max_snp * sizeof(double));
  geno = malloc(niter * 2 * sizeof(int*));
  sampstate = malloc(niter * sizeof(int*));
  assert(ffreq != NULL);
  assert(pos != NULL);
  for (ksamp = 0; ksamp < niter; ksamp++) {
    sampstate[ksamp] = malloc(max_snp * sizeof(int));
    assert(sampstate[ksamp] != NULL);
  }

  for (ksamp = 0; ksamp < 2*niter; ksamp++) {
    geno[ksamp] = malloc(max_snp * sizeof(int));
    assert(geno[ksamp] != NULL);
  }

  while (fgets(newLine, linesize, ff) != NULL) {
    for (running = newLine, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
      if (itoken == 0) {
	chr = strtol(token, NULL, 10);
      }
      else if (itoken == 1) {
	pos[nsnp] = strtol(token, NULL, 10);
      }
      else if (itoken == 2) {
	ffreq[nsnp] = strtod(token, NULL);
	break;
      }
    }
    if (chr > nchrom) {continue;}
    if (nsnp < start_chr[chr]) {
      start_chr[chr] = nsnp;
    }
    if (nsnp > end_chr[chr]) {
      end_chr[chr] = nsnp;
    }
    nsnp++;
    if (nsnp == max_snp) {
      pos = realloc(pos, 2*max_snp*sizeof(int));
      ffreq = realloc(ffreq, 2*max_snp*sizeof(double));
      assert(pos != NULL);
      assert(ffreq != NULL);
      for (ksamp = 0; ksamp < niter; ksamp++) {
	sampstate[ksamp] = realloc(sampstate[ksamp], 2*max_snp*sizeof(int));
	assert(sampstate[ksamp] != NULL);
      }
      for (ksamp = 0; ksamp < 2*niter; ksamp++) {
	geno[ksamp] = realloc(geno[ksamp], 2*max_snp*sizeof(int));
	assert(geno[ksamp] != NULL);
      }
      max_snp *= 2;
    }
  }
  
  outf = fopen(outfile, "w");
  assert(outf != NULL);
  fprintf(outf, "chrom\tpos");
  sumf = fopen(sumfile, "w");
  assert(sumf != NULL);
  fprintf(sumf, "samp1\tsamp2\ttarget_fract_ibd\trec_coeff\tgeno_error\tpi\tN_trans\n");
  statef = fopen(statefile, "w");
  assert(statef != NULL);
  fprintf(statef, "chr\tpos");
  for (iter = 0; iter < niter; iter++) {
    fprintf(outf, "\t%ds\t%ds_%ds", 2*iter, 2*iter, 2*iter+1);
    fprintf(statef, "\t%ds:%ds_%ds", 2*iter, 2*iter, 2*iter+1);
  }
  fprintf(outf, "\n");
  fprintf(statef, "\n");
  for (iter = 0; iter < niter; iter++) {
    ntrans = nstate[0] = nstate[1] = 0;
    for (chr = 1; chr <= 14; chr++) {
      // Pick start state (0 = ibd)
      state = 1;   // ibd or not
      x = random_double();
      if (x < target_pi[0]) {state = 0;}
      for (isnp = start_chr[chr]; isnp <= end_chr[chr]; isnp++) {
	nstate[state]++;
	// IBD (state == 0)
	p11 = pr * pr * ffreq[isnp] + eps * eps * (1. - ffreq[isnp]);
	p22 = pr * pr * (1. - ffreq[isnp]) + eps * eps * ffreq[isnp];
	p12 = 2 * pr * eps;
	if (state == 1) {  // DBD
	  p11 = pr * pr * ffreq[isnp] * ffreq[isnp] 
	    + 2 * eps * pr * ffreq[isnp] * (1 - ffreq[isnp]) + 
	    eps * eps * (1 - ffreq[isnp]) * (1 - ffreq[isnp]);
	  p22 = pr * pr * (1 - ffreq[isnp]) * (1 - ffreq[isnp]) 
	    + 2 * eps * pr * (1 - ffreq[isnp]) * ffreq[isnp] +
	    eps * eps * ffreq[isnp] * ffreq[isnp];
	  p12 = 2 * pr * pr * ffreq[isnp] * (1 - ffreq[isnp]) 
	    + 2 * eps * pr * ffreq[isnp] * ffreq[isnp] 
	    + 2 * eps * pr * (1 - ffreq[isnp]) * (1 - ffreq[isnp])
	    + 2 * eps * eps * ffreq[isnp] * (1 - ffreq[isnp]);
	}
	x = random_double();
	sampstate[iter][isnp] = state;
	if (x < p11) {
	  geno[2*iter][isnp] = geno[2*iter+1][isnp] = 0;
	}
	else if (x < p11 + p22) {
	  geno[2*iter][isnp] = geno[2*iter+1][isnp] = 1;
	}
	else {
	  geno[2*iter][isnp] = 1;
	  geno[2*iter+1][isnp] = 0;
	  x = random_double();
	  if (x < 0.5) {
	    geno[2*iter][isnp] = 0;
	    geno[2*iter+1][isnp] = 1;
	  }
	}
	// add genotyping error
	//   correction, should already be included in emission probs
	//	x = random_double();
	//	if (x < eps) {
	//	  geno[2*iter][isnp] = geno[2*iter][isnp] % 2;
	//	}
	//	x = random_double();
	//	if (x < eps) {
	//	  geno[2*iter+1][isnp] = geno[2*iter+1][isnp] % 2;
	//	}

	// Change state?
	if (isnp < end_chr[chr]) {	
	  dist = pos[isnp+1] - pos[isnp];
	  ptrans = 1 - target_pi[state] - (1-target_pi[state]) * exp(-rec_rate*k_rec*dist);
	  //	    fprintf(stderr, "state: %d  pi: %.2e dist: %d ptrans: %.5e\n", which, target_pi[which], 
	  //	    		    dist, ptrans);
	  x = random_double();
	  if (x < ptrans) {
	    state = (state+1)%2;
	    ntrans++;
	  }
	}
      }  // end loop over snps
    }  // end loop over chroms
    pi = (double) nstate[0] / (nstate[0] + nstate[1]);
    fprintf(sumf, "%ds\t%ds_%ds\t%.5f\t%.3f\t%.5f\t%.5f\t%d\n", 
	    2*iter, 2*iter, 2*iter+1, target_pi[0], k_rec, eps, pi, ntrans);
  }   // end iter loop
  
  for (chr = 1; chr <= nchrom; chr++) {
    for (isnp = start_chr[chr]; isnp <= end_chr[chr]; isnp++) {
      fprintf(outf, "%d\t%d", chr, pos[isnp]);
      fprintf(statef, "%d\t%d", chr, pos[isnp]);
      for (iter = 0; iter < niter; iter++) {
	fprintf(outf, "\t%d\t%d", geno[2*iter][isnp], geno[2*iter+1][isnp]);
	fprintf(statef, "\t%d", sampstate[iter][isnp]);
      }
      fprintf(outf, "\n");
      fprintf(statef, "\n");
    }    
  }

  exit(0);
}
