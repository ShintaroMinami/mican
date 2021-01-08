#include "mican.h"

/** function for comparison in quick sort **/
PRIVATE int compare(const void *a, const void *b);

typedef struct {
  int naa_align;
  int ich;
  int jch;
} CHALN_DAT;

float C3score(int nch1, int nch2, int naa_align, ALIGN **align_chain) {
  int ich, jch;
  int ichaln, nchaln;
  float C3score = 0.0;
  CHALN_DAT *chaln;
  int *check_ch1;
  int *check_ch2;
  int naa_C3score = 0;
  
  chaln = (CHALN_DAT *)calloc((size_t)(nch1*nch2), sizeof(CHALN_DAT));
  check_ch1 = (int *)calloc((size_t)(nch1+1), sizeof(int));
  check_ch2 = (int *)calloc((size_t)(nch2+1), sizeof(int));

  /**  Copy chain-align info  **/
  nchaln = -1;
  for(ich=1; ich<=nch1; ich++){
    for(jch=1; jch<=nch2; jch++){
      nchaln ++;
      chaln[nchaln].ich = ich;
      chaln[nchaln].jch = jch;
      chaln[nchaln].naa_align = align_chain[ich][jch].naa_align;
    }
  }
  /**  Quick Sort  **/
  qsort(chaln, (size_t)(nchaln+1), sizeof(CHALN_DAT), compare);

  /**  count aligned residues for C3 score  **/
  for(ich=1; ich<=nch1; ich++){ check_ch1[ich] = OFF; }
  for(ich=1; ich<=nch2; ich++){ check_ch2[ich] = OFF; }
  naa_C3score = 0;

  for(ichaln=0; ichaln<=nchaln; ichaln++){
    if(check_ch1[chaln[ichaln].ich] == OFF && check_ch2[chaln[ichaln].jch] == OFF){
      naa_C3score += chaln[ichaln].naa_align;
      check_ch1[chaln[ichaln].ich] = ON;
      check_ch2[chaln[ichaln].jch] = ON;
    }
  }
  /**  C3 score  **/
  C3score = (float)(naa_C3score) / (float)(naa_align);

  /** Free **/
  free(chaln);
  free(check_ch1);
  free(check_ch2);
  
  return(C3score);
}


/****** function for comparison ******/
PRIVATE int compare(const void *a, const void *b) {
  CHALN_DAT *x, *y;
  float delta;
  x = (CHALN_DAT *)a;
  y = (CHALN_DAT *)b;
  delta = (*y).naa_align - (*x).naa_align;

  if (delta < 0) {
    return -1;
  }
  if (delta > 0) {
    return 1;
  }
  return 0;
}

