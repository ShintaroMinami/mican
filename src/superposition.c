#include "mican.h"

int superposition(int naa_q, RESDAT *resdat_q, RESDAT *resdat_t, ALIGN align, int sw) {

  int iaa, jaa, naa;
  int ico;
  float *RMSD_weight;
  RESDAT *resdat1;
  RESDAT *resdat2;

  resdat1 = (RESDAT *)calloc((size_t)(naa_q + 1), sizeof(RESDAT));
  resdat2 = (RESDAT *)calloc((size_t)(naa_q + 1), sizeof(RESDAT));

  RMSD_weight = (float *)calloc(naa_q+1, sizeof(float));
  /**  get aligned coordinate  **/
  naa = 0;
  for (iaa = 1; iaa <= naa_q; iaa++) {
    if (align.align_q[iaa] > 0) {
      naa++;
      jaa = align.align_q[iaa];
      for (ico = 1; ico <= Ncoord; ico++) {
	resdat1[naa].cacoord[ico] = resdat_q[iaa].cacoord[ico];
	resdat2[naa].cacoord[ico] = resdat_t[jaa].cacoord[ico];
      }
      RMSD_weight[naa] = align.local_score[iaa];
    }
  }
  /**  calc best fit rotation matrix  **/
  if (naa > 2) {
    if (sw == ON) {
      wbstft(naa, resdat1, resdat2, align.vec, align.rot, RMSD_weight);
    } else {
      bstft(naa, resdat1, resdat2, align.vec, align.rot);
    }
  }

  /**  free  **/
  free(resdat1);
  free(resdat2);
  free(RMSD_weight);

  if (naa > 2) {
    return (TRUE);
  } else {
    return (FALSE);
  }
}
