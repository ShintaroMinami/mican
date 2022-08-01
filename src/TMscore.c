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
  co1 = (float **)calloc((size_t)(naln + 1), sizeof(float *));
  co2 = (float **)calloc((size_t)(naln + 1), sizeof(float *));
  for (ialn = 1; ialn <= naln; ialn++) {
    co1[ialn] = (float *)calloc(3 + 1, sizeof(float));
    co2[ialn] = (float *)calloc(3 + 1, sizeof(float));
  }

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
  for (ialn = 1; ialn <= naln; ialn++) {
    free(co1[ialn]);
    free(co2[ialn]);
  }
  free(co1);
  free(co2);

  return (TMscore);
}
