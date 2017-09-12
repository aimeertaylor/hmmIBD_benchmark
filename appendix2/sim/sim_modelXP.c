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
  const double rec_rate = 5.8e-7;   // 5.8e-5 cM/bp, or 17kb/cM
  const int niter = 100;

  FILE *outf1=NULL, *outf2=NULL, *sumf=NULL, *ff1=NULL, *ff2=NULL, *statef=NULL;
  int pflag, kflag, iflag, c, itoken, **geno, chr, *pos, nsnp, **sampstate=NULL;
  int *start_chr=NULL, *end_chr=NULL, state, ksamp;
  int isnp, ntrans, dist, nstate[2], eflag, freq_flag1, freq_flag2, iter;
  int oflag1, oflag2;
  char *newLine=NULL, *token, *running;
  char sumfile[256], statefile[256];
  double x, ptrans, pi, *freq1=NULL, *freq2=NULL, pr, p11, p12, p22, p21, fIBD;
  
  newLine = malloc((linesize+1) * sizeof(char));
  start_chr = malloc((nchrom+1) * sizeof(int));
  end_chr = calloc(nchrom+1, sizeof(int));
  for (chr = 1; chr <= nchrom; chr++) {start_chr[chr] = 10000000;}
  pflag = kflag = iflag = freq_flag1 = freq_flag2 = oflag1 = oflag2 = 0;
  seed_rng();   
  //  set_rng_seed(111111);
  while ( (c = getopt(argc, argv, ":p:k:F:e:f:o:O:")) != -1) {
    switch(c) {
    case 'p':
      target_pi[0] = strtod(optarg, NULL);
      target_pi[1] = 1 - target_pi[0];
      fprintf(stderr, "pi set to %f\n", target_pi[0]);
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
      ff1 = fopen(optarg, "r");
      if (ff1 == NULL) {fprintf(stderr, "could not open frequency file %s\n", optarg); exit(0);}
      freq_flag1 = 1;
      break;
    case 'F':
      ff2 = fopen(optarg, "r");
      if (ff2 == NULL) {fprintf(stderr, "could not open frequency file %s\n", optarg); exit(0);}
      freq_flag2 = 1;
      break;
    case 'o':
      outf1 = fopen(optarg, "w");
      if (outf1 == NULL) {fprintf(stderr, "could not open output file %s\n", optarg); exit(0);}
      oflag1 = 1;
      break;
    case 'O':
      outf2 = fopen(optarg, "w");
      if (outf2 == NULL) {fprintf(stderr, "could not open output file %s\n", optarg); exit(0);}
      oflag2 = 1;
      break;
    case ':':
      fprintf(stdout, "argument %c requires an argument\n", optopt);
      abort();
    }
  }
  if (optind != argc || kflag == 0 || pflag == 0 || eflag == 0  || freq_flag1 == 0 || freq_flag2 == 0 
      || oflag1 == 0 || oflag2 == 0) {
    //    fprintf(stderr, "%c/%c: k: %d p: %d e: %d f: %d F: %d\n", optind, argc, kflag, pflag, eflag, freq_flag1, freq_flag2);
    fprintf(stderr, "Usage: sim_model -k <recomb coeff> -p <fraction IBD> -o <outfile 1> -O <outfile 2> ");
    fprintf(stderr, "-e <genotype err rate> -f <allele freq file 1> -F <freq file 2>\n");
    exit(0);
  }
  sprintf(sumfile, "simsumXP.txt");
  sprintf(statefile, "sim_state.txt");
  pr = 1 - eps;
    
  // Read frequencies from files
  nsnp = 0;
  pos = malloc(max_snp * sizeof(int));
  freq1 = malloc(max_snp * sizeof(double));
  freq2 = malloc(max_snp * sizeof(double));
  geno = malloc(niter * 2 * sizeof(int*));
  sampstate = malloc(niter * sizeof(int*));
  assert(freq1 != NULL);
  assert(freq2 != NULL);
  assert(pos != NULL);
  for (ksamp = 0; ksamp < niter; ksamp++) {
    sampstate[ksamp] = malloc(max_snp * sizeof(int));
    assert(sampstate[ksamp] != NULL);
  }
  for (ksamp = 0; ksamp < 2*niter; ksamp++) {
    geno[ksamp] = malloc(max_snp * sizeof(int));
    assert(geno[ksamp] != NULL);
  }

  while (fgets(newLine, linesize, ff1) != NULL) {
    for (running = newLine, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
      if (itoken == 0) {
	chr = strtol(token, NULL, 10);
      }
      else if (itoken == 1) {
	pos[nsnp] = strtol(token, NULL, 10);
      }
      else if (itoken == 2) {
	freq1[nsnp] = strtod(token, NULL);
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
      freq1 = realloc(freq1, 2*max_snp*sizeof(double));
      freq2 = realloc(freq2, 2*max_snp*sizeof(double));
      assert(pos != NULL);
      assert(freq1 != NULL);
      assert(freq2 != NULL);
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
  isnp = 0;
  while (fgets(newLine, linesize, ff2) != NULL) {
    for (running = newLine, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
      if (itoken == 0) {
      }
      else if (itoken == 1) {
	assert(pos[isnp] == strtol(token, NULL, 10));
      }
      else if (itoken == 2) {
	freq2[isnp] = strtod(token, NULL);
	break;
      }
    }
    isnp++;
  }
  
  fprintf(outf1, "chrom\tpos");
  fprintf(outf2, "chrom\tpos");
  sumf = fopen(sumfile, "w ");
  assert(sumf != NULL);
  fprintf(sumf, "samp1\tsamp2\ttarget_fract_ibd\trec_coeff\tgeno_error\tpi\tN_trans\n");
  statef = fopen(statefile, "w");
  assert(statef != NULL);
  fprintf(statef, "chr\tpos");
  for (iter = 0; iter < niter; iter++) {
    fprintf(outf1, "\t%ds1", iter);
    fprintf(outf2, "\t%ds2", iter);
    fprintf(statef, "\t%ds1:%ds2", iter, iter);
  }
  fprintf(outf1, "\n");
  fprintf(outf2, "\n");
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
	fIBD = 0.5 * (freq1[isnp] + freq2[isnp]);
	// IBD (state == 0)
	p11 = pr * pr * fIBD + eps * eps * (1. - fIBD);
	p22 = pr * pr * (1. - fIBD) + eps * eps * fIBD;
	p12 = p21 = pr * eps;   
	if (state == 1) {  // DBD
	  p11 = pr * pr * freq1[isnp] * freq2[isnp] 
	    + eps * pr * freq1[isnp] * (1 - freq2[isnp]) 
	    + eps * pr * freq2[isnp] * (1 - freq1[isnp]) + 
	    eps * eps * (1 - freq1[isnp]) * (1 - freq2[isnp]);
	  p22 = pr * pr * (1 - freq1[isnp]) * (1 - freq2[isnp]) 
	    + eps * pr * (1 - freq1[isnp]) * freq2[isnp]
	    + eps * pr * (1 - freq2[isnp]) * freq1[isnp] +
	    eps * eps * freq1[isnp] * freq2[isnp];
	  p12 = pr * pr * freq1[isnp] * (1 - freq2[isnp]) 
	    + pr * eps * (freq1[isnp] * freq2[isnp] + (1 - freq1[isnp]) * (1 - freq2[isnp]))
	    + eps * eps * (1 - freq1[isnp]) * freq2[isnp];
	  p21 = pr * pr * (1 - freq1[isnp]) * freq2[isnp] 
	    + pr * eps * (freq1[isnp] * freq2[isnp] + (1 - freq1[isnp]) * (1 - freq2[isnp]))
	    + eps * eps * freq1[isnp] * (1 - freq2[isnp]);
	}
	//	fprintf(stderr, "%d/%d state: %d p11: %.6f p22 %.6f p12: %.6f p21: %.6f", 
	//		chr, pos[isnp], state, p11, p22, p12, p21);
	x = random_double();
	sampstate[iter][isnp] = state;
	if (x < p11) {
	  geno[2*iter][isnp] = geno[2*iter+1][isnp] = 0;
	  //	  fprintf(stderr, " p11 freq1: %.5f freq2: %.5f", freq1[isnp], freq2[isnp]);
	}
	else if (x < p11 + p22) {
	  geno[2*iter][isnp] = geno[2*iter+1][isnp] = 1;
	  //	  fprintf(stderr, " p22 freq1: %.5f freq2: %.5f", freq1[isnp], freq2[isnp] );
	}
	else if (x < p11 + p22 + p12) {
	  geno[2*iter][isnp] = 0;
	  geno[2*iter+1][isnp] = 1;
	  //	  fprintf(stderr, " p12 freq1: %.5f freq2: %.5f", freq1[isnp], freq2[isnp] );
	}
	else {
	  geno[2*iter][isnp] = 1;
	  geno[2*iter+1][isnp] = 0;
	  //	  fprintf(stderr, " p21 freq1: %.5f freq2: %.5f", freq1[isnp], freq2[isnp] );
	}
	//	fprintf(stderr, " geno: %d/%d\n", geno[2*iter][isnp], geno[2*iter+1][isnp]);
	//	if (fabs(1 - (p11+p22+p12) - p21) > .002) {fprintf(stderr, "sum: %.6f p21: %.6f\n", p11 + p22 + p12, p21);}
	// add genotyping error
	//  genotyping error should already be included in the emission probs
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
	  //	    dist, ptrans);
	  x = random_double();
	  if (x < ptrans) {
	    state = (state+1)%2;
	    ntrans++;
	  }
	}
      }  // end loop over snps
    }  // end loop over chroms
    pi = (double) nstate[0] / (nstate[0] + nstate[1]);
    fprintf(sumf, "%ds1\t%ds2\t%.5f\t%.3f\t%.5f\t%.5f\t%d\n", 
	    iter, iter, target_pi[0], k_rec, eps, pi, ntrans);
  }   // end iter loop
  
  for (chr = 1; chr <= nchrom; chr++) {
    for (isnp = start_chr[chr]; isnp <= end_chr[chr]; isnp++) {
      fprintf(outf1, "%d\t%d", chr, pos[isnp]);
      fprintf(outf2, "%d\t%d", chr, pos[isnp]);
      fprintf(statef, "%d\t%d", chr, pos[isnp]);
      for (iter = 0; iter < niter; iter++) {
	fprintf(outf1, "\t%d", geno[2*iter][isnp]);
	fprintf(outf2, "\t%d", geno[2*iter+1][isnp]);
	fprintf(statef, "\t%d", sampstate[iter][isnp]);
      }
      fprintf(outf1, "\n");
      fprintf(outf2, "\n");
      fprintf(statef, "\n");
    }    
  }

  exit(0);
}
