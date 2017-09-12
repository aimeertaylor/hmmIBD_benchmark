#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

// special version to handle isorelate output

int main(int argc, char **argv) {
  const int min_size = 2000;
  const int min_snp = 100;
  const int max_chrom = 40, linesize = 500;
  FILE *inf=NULL, *outf=NULL;
  int itoken, pos, chrom, start, end, nchrom=0, *chrom_start, *chrom_end, ichr;
  int **ibd_sum, count, *running_start, runtot, state, old_val, gpos, nsnp;
  char *newLine=NULL, *token, *running, outfile[128];

  running_start = malloc(max_chrom * sizeof(int));
  chrom_start = malloc(max_chrom * sizeof(int));
  chrom_end = malloc(max_chrom * sizeof(int));
  newLine = malloc((linesize+1) * sizeof(char));
  assert(newLine != NULL);
  if (argc != 2) {fprintf(stderr, "usage: sum_overlaps <hmm output file>\n"); exit(0);}
  inf = fopen(argv[1], "r");
  assert(inf != NULL);
  sprintf(outfile, "%s.summed", argv[1]);
  outf = fopen(outfile, "w");
  assert(outf != NULL);
  fprintf(outf, "chrom\tpos\tgenome_pos\tN_ident_pair\n");

  for (ichr = 0; ichr < max_chrom; ichr++) {
    chrom_end[ichr] = 0;
    chrom_start[ichr] = -1;
  }
  fgets(newLine, linesize, inf);   // header
  while (fgets(newLine, linesize, inf) != NULL) {
    for (running = newLine, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
      if (itoken == 2) {
	chrom = strtol(token, NULL, 10);
	if (chrom > max_chrom) {fprintf(stderr, "bad: %s\n", newLine);}
      }
      else if (itoken == 3) {start = strtol(token, NULL, 10);}
      else if (itoken == 4) {
	end = strtol(token, NULL, 10);
	if (chrom_start[chrom] == -1) {chrom_start[chrom] = start;}
	if (chrom_end[chrom] < end) {chrom_end[chrom] = end;}
	if (chrom > nchrom) {nchrom = chrom;}
	break;
      }
    }
  }
  ibd_sum = malloc((nchrom+1) * sizeof(int*));
  for (ichr = 1; ichr <= nchrom; ichr++) {
    chrom_start[ichr] = 1;
    ibd_sum[ichr] = calloc((chrom_end[ichr]+1), sizeof(int));
    assert(ibd_sum[ichr] != NULL);
  }
  fseek(inf, 0, 0);
  fgets(newLine, linesize, inf);   // header
  count = 0;
  while (fgets(newLine, linesize, inf) != NULL) {
    for (running = newLine, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
      if (itoken == 2) {chrom = strtol(token, NULL, 10);}
      else if (itoken == 3) {start = strtol(token, NULL, 10);}
      else if (itoken == 4) {end = strtol(token, NULL, 10);}
      else if (itoken == 5) {
	if (end - start + 1 >= min_size) {
	  state = (1+strtol(token, NULL, 10)) % 2;
	}
      }
      else if (itoken == 6) {
	nsnp = strtol(token, NULL, 10);
	if (nsnp > min_snp) {
	  for (pos = start; pos <= end; pos++) {
	    ibd_sum[chrom][pos] += state;
	  }
	}
	break;
      }
    }
    //    if (count % 100 == 0) {fprintf(stderr, "%d\n", count);}
    count++;
  }
  runtot = 0;
  for (chrom = 1; chrom <= nchrom; chrom++) {
    running_start[chrom] = runtot;
    runtot += chrom_end[chrom] - chrom_start[chrom] + 1;
  }

  for (chrom = 1; chrom <= nchrom; chrom++) {
    old_val = 0;
    for (pos = chrom_start[chrom]; pos <= chrom_end[chrom]; pos++) {
      if (pos == chrom_start[chrom]) {
	gpos = pos - chrom_start[chrom] + running_start[chrom];
	fprintf(outf, "%d\t%d\t%d\t%d\n", chrom, pos, gpos, ibd_sum[chrom][pos]);
      }
      else if (pos == chrom_end[chrom]) {
	gpos = pos - chrom_start[chrom] + running_start[chrom];	
	if (ibd_sum[chrom][pos] == ibd_sum[chrom][pos-1]) {
	  fprintf(outf, "%d\t%d\t%d\t%d\n", chrom, pos, gpos, ibd_sum[chrom][pos]);
	}
	else {
	  fprintf(outf, "%d\t%d\t%d\t%d\n", chrom, pos-1, gpos-1, ibd_sum[chrom][pos-1]);
	  fprintf(outf, "%d\t%d\t%d\t%d\n", chrom, pos, gpos, ibd_sum[chrom][pos]);
	  fprintf(outf, "%d\t%d\t%d\t%d\n", chrom, pos, gpos, ibd_sum[chrom][pos]);
	}
      }
      else if (ibd_sum[chrom][pos] != old_val) {
	gpos = pos - chrom_start[chrom] + running_start[chrom];
	fprintf(outf, "%d\t%d\t%d\t%d\n", chrom, pos-1, gpos-1, ibd_sum[chrom][pos-1]);	
	fprintf(outf, "%d\t%d\t%d\t%d\n", chrom, pos, gpos, ibd_sum[chrom][pos]);
      }
      old_val = ibd_sum[chrom][pos];
    }
  }  

  return 0;
}
