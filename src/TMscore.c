#include "mican.h"
#include "main.h"

float TMscore(int naa0, int naa1, RESDAT *resdat1, RESDAT *resdat2, ALIGN align) {
  int iaa, jaa;
  int ialn, naln;
  int ico;
  int i;
  float **co1;
  float **co2;
  float TMscore;
  float d, d0;

  naln = align.naa_align;

  /*******************/
  /**     calloc    **/
  /*******************/
  co1 = calloc_coords(naln);
  co2 = calloc_coords(naln);

  /***********************/
  /**   copy coordinate  **/
  /***********************/
  ialn = 0;
  for (iaa = 1; iaa <= naa1; iaa++) {
    jaa = align.align_q[iaa];
    if(jaa < 1) continue;
    ialn ++;
    for (ico = 1; ico <= Ncoord; ico++) {
      co1[ialn][ico] = resdat1[iaa].cacoord[ico];
    }
    for (i = 1; i <= 3; i++) {
      co2[ialn][i] = 0;
      for (ico = 1; ico <= Ncoord; ico++) {
        co2[ialn][i] += align.rot[i][ico] * resdat2[jaa].cacoord[ico];
      }
      co2[ialn][i] += align.vec[i];
    }
  }

  /**************************/
  /**       TM score       **/
  /**************************/
  d0 = input.TM_d0;
  d0 = d0 * d0;

  /****  calc  ****/
  TMscore = 0;
  for (ialn = 1; ialn <= naln; ialn++) {
    d = 0;
    for (i = 1; i <= 3; i++) {
      d += (co1[ialn][i] - co2[ialn][i]) * (co1[ialn][i] - co2[ialn][i]);
    }

    TMscore += 1 / (1 + d / d0);
  }
  TMscore /= (float)naa0;

  /**************/
  /**   free   **/
  /**************/
  free_coords(co1);
  free_coords(co2);

  return (TMscore);
}
