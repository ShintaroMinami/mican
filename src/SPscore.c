#include "mican.h"

float SPscore(int naa1, int naa2, RESDAT *resdat1, RESDAT *resdat2,
              ALIGN align) {
  int i;
  int iaa, jaa;
  int irem, nrem1, nrem2;
  int icore, ncore;
  int ico;
  int check;
  float N_e;
  float SP_core;
  float SPscore;
  float d, d02, dcore2, dsurr2;
  int *ALcheck1, *ALcheck2;
  int *COREcheck1, *COREcheck2;
  int *iaa1_core, *iaa2_core;
  int *iaa1_rem, *iaa2_rem;
  float **co1;
  float **co2;

  d02 = SP_d0 * SP_d0;
  dcore2 = SP_dcore * SP_dcore;
  dsurr2 = SP_dsurr * SP_dsurr;

  /*******************/
  /**     calloc    **/
  /*******************/
  /**  coordinates  **/
  co1 = (float **)calloc((size_t)(naa1 + 1), sizeof(float *));
  co2 = (float **)calloc((size_t)(naa2 + 1), sizeof(float *));
  for (iaa = 1; iaa <= naa1; iaa++) {
    co1[iaa] = (float *)calloc(3 + 1, sizeof(float));
  }
  for (iaa = 1; iaa <= naa2; iaa++) {
    co2[iaa] = (float *)calloc(3 + 1, sizeof(float));
  }
  /**  core residues  **/
  iaa1_core = (int *)calloc((size_t)(naa1 + 1), sizeof(int));
  iaa2_core = (int *)calloc((size_t)(naa2 + 1), sizeof(int));
  /**  surrounding residues  **/
  iaa1_rem = (int *)calloc((size_t)(naa1 + 1), sizeof(int));
  iaa2_rem = (int *)calloc((size_t)(naa2 + 1), sizeof(int));
  /**  alignment checker  **/
  ALcheck1 = (int *)calloc((size_t)(naa1 + 1), sizeof(int));
  ALcheck2 = (int *)calloc((size_t)(naa2 + 1), sizeof(int));
  COREcheck1 = (int *)calloc((size_t)(naa1 + 1), sizeof(int));
  COREcheck2 = (int *)calloc((size_t)(naa2 + 1), sizeof(int));

  /***********************/
  /**  copy coordinate  **/
  /***********************/
  for (iaa = 1; iaa <= naa1; iaa++) {
    for (ico = 1; ico <= Ncoord; ico++) {
      co1[iaa][ico] = resdat1[iaa].cacoord[ico];
    }
  }
  for (iaa = 1; iaa <= naa2; iaa++) {
    for (i = 1; i <= Ncoord; i++) {
      co2[iaa][i] = 0;
      for (ico = 1; ico <= Ncoord; ico++) {
        co2[iaa][i] += align.rot[i][ico] * resdat2[iaa].cacoord[ico];
      }
      co2[iaa][i] += align.vec[i];
    }
  }

  /**************************/
  /**      core score      **/
  /**************************/
  SP_core = 0;
  ncore = 0;
  for (iaa = 1; iaa <= naa1; iaa++) {
    jaa = align.align_q[iaa];
    if (jaa != 0) {
      /**  check  **/
      ALcheck1[iaa] = 1;
      ALcheck2[jaa] = 1;

      /**  calc score  **/
      d = 0;
      for (i = 1; i <= 3; i++) {
        d += (co1[iaa][i] - co2[jaa][i]) * (co1[iaa][i] - co2[jaa][i]);
      }
      if (d <= dcore2) {
        SP_core += 1 / (1 + d / d02) - SP_const;
        ncore++;
        iaa1_core[ncore] = iaa;
        iaa2_core[ncore] = jaa;
        COREcheck1[iaa] = 1;
        COREcheck2[jaa] = 1;
      }
    }
  }

  /**************************/
  /**  remaining residues  **/
  /**************************/
  nrem1 = 0;
  for (iaa = 1; iaa <= naa1; iaa++) {
    if (COREcheck1[iaa] != 1) {
      nrem1++;
      iaa1_rem[nrem1] = iaa;
    }
  }
  nrem2 = 0;
  for (iaa = 1; iaa <= naa2; iaa++) {
    if (COREcheck2[iaa] != 1) {
      nrem2++;
      iaa2_rem[nrem2] = iaa;
    }
  }

  /****************************/
  /**    effective length    **/
  /****************************/
  /**  for protein1  **/
  N_e = 0;
  for (irem = 1; irem <= nrem1; irem++) {
    check = 0;
    for (icore = 1; icore <= ncore; icore++) {
      d = 0;
      for (i = 1; i <= 3; i++) {
        d += (co1[iaa1_rem[irem]][i] - co1[iaa1_core[icore]][i]) *
             (co1[iaa1_rem[irem]][i] - co1[iaa1_core[icore]][i]);
      }
      if (d <= dsurr2) {
        check = 1;
        break;
      }
      d = 0;
      for (i = 1; i <= 3; i++) {
        d += (co1[iaa1_rem[irem]][i] - co2[iaa2_core[icore]][i]) *
             (co1[iaa1_rem[irem]][i] - co2[iaa2_core[icore]][i]);
      }
      if (d <= dsurr2) {
        check = 1;
        break;
      }
    }
    N_e += (float)check;
  }

  /**  for protein2  **/
  for (irem = 1; irem <= nrem2; irem++) {
    check = 0;
    for (icore = 1; icore <= ncore; icore++) {
      d = 0;
      for (i = 1; i <= 3; i++) {
        d += (co2[iaa2_rem[irem]][i] - co1[iaa1_core[icore]][i]) *
             (co2[iaa2_rem[irem]][i] - co1[iaa1_core[icore]][i]);
      }
      if (d <= dsurr2) {
        check = 1;
        break;
      }
      d = 0;
      for (i = 1; i <= 3; i++) {
        d += (co2[iaa2_rem[irem]][i] - co2[iaa2_core[icore]][i]) *
             (co2[iaa2_rem[irem]][i] - co2[iaa2_core[icore]][i]);
      }
      if (d <= dsurr2) {
        check = 1;
        break;
      }
    }
    N_e += (float)check;
  }

  N_e = 3 * (float)pow((float)align.naa_align + N_e / 2, SP_power);

  /*************************/
  /**    calc SP-score    **/
  /*************************/
  SPscore = SP_core / N_e;

  /**************/
  /**   free   **/
  /**************/
  for (iaa = 1; iaa <= naa1; iaa++) {
    free(co1[iaa]);
  }
  for (iaa = 1; iaa <= naa2; iaa++) {
    free(co2[iaa]);
  }
  free(co1);
  free(co2);

  free(ALcheck1);
  free(ALcheck2);
  free(COREcheck1);
  free(COREcheck2);
  free(iaa1_core);
  free(iaa2_core);
  free(iaa1_rem);
  free(iaa2_rem);

  return (SPscore);
}
