#include "mican.h"
#include "main.h"

typedef struct {
  int iaa;
  int jaa;
  int segmode;
  float dist2;
} TERMINI;

PRIVATE int compare_termini(const void *a, const void *b);

void extend_termini(int naa_q, int naa_t, ALIGN align, float **co_q, float **co_t) {

  int iaa, jaa;
  int itermini, ntermini;
  int mode;
  int *align_t;
  float dist2;
  TERMINI *termini;

  align_t = (int *)calloc((size_t)(naa_t + 1), sizeof(int));
  termini = (TERMINI *)calloc((size_t)(naa_q+1), sizeof(TERMINI));

  /**  initialize  **/
  for(iaa = 1; iaa <= naa_q; iaa ++) align.termini[iaa] = OFF;

  /**  copy alignment  **/
  for(iaa = 1; iaa <= naa_q; iaa ++) align_t[align.align_q[iaa]] = iaa;

  /**  search termini  **/
  ntermini = -1;
  for(iaa = 2; iaa < naa_q; iaa ++) {
    if(align.align_q[iaa] != OFF) {
      continue;
    }
    else {
      if(align.align_q[iaa - 1] != OFF) {
	/** termini 1 **/
	jaa = OFF;
	if(align.segmode[iaa - 1] == FORWARD){
	  jaa = align.align_q[iaa - 1] + 1;
	  mode = FORWARD;
	}
	else if(align.segmode[iaa - 1] == REVERSE) {
	  jaa = align.align_q[iaa - 1] - 1;
	  mode = REVERSE;
	}
	/** check **/
	if(1 <= jaa && jaa <= naa_t) {
	  dist2 = distance_2(co_q[iaa], co_t[jaa]);
	  if(dist2 < (input.TM_d0*input.TM_d0)) {
	    ntermini ++;
	    termini[ntermini].iaa = iaa;
	    termini[ntermini].jaa = jaa;
	    termini[ntermini].dist2 = dist2;
	    termini[ntermini].segmode = mode;
	  }
	}
      }
      if(align.align_q[iaa + 1] != OFF) {
	/** termini 2 **/
	jaa = OFF;
	if(align.segmode[iaa + 1] == FORWARD){
	  jaa = align.align_q[iaa + 1] - 1;
	  mode = FORWARD;
	}
	else if(align.segmode[iaa + 1] == REVERSE) {
	  jaa = align.align_q[iaa + 1] + 1;
	  mode = REVERSE;
	}
	/** check **/
	if(1 <= jaa && jaa <= naa_t) {
	  dist2 = distance_2(co_q[iaa], co_t[jaa]);
	  if(dist2 < (input.TM_d0*input.TM_d0)) {
	    ntermini ++;
	    termini[ntermini].iaa = iaa;
	    termini[ntermini].jaa = jaa;
	    termini[ntermini].dist2 = dist2;
	    termini[ntermini].segmode = mode;
	  }
	}
      }
    }
  }

  /****  qsort()  ****/
  qsort(termini, (size_t)(ntermini+1), sizeof(TERMINI), compare_termini);

  /****  extend termini  ****/
  for(itermini = 0; itermini <= ntermini; itermini ++) {
    iaa = termini[itermini].iaa;
    jaa = termini[itermini].jaa;
    if(align.align_q[iaa] == OFF && align_t[jaa] == OFF) {
      align.align_q[iaa] = jaa;
      align.segmode[iaa] = termini[itermini].segmode;
      align.termini[iaa] = ON;
      align_t[jaa] = iaa;
    }
  }

  free(align_t);
  free(termini);

  return;
}




/**** function for qsort() ****/
PRIVATE int compare_termini(const void *a, const void *b) {
  TERMINI *x, *y;
  float delta;
  x = (TERMINI *)a;
  y = (TERMINI *)b;
  delta = (*y).dist2 - (*x).dist2;

  if(delta < 0){
    return 1;
  }
  if(delta > 0){
    return -1;
  }
  return 0;
}
