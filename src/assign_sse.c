#include "mican.h"

/************************************************/
/**     FUNCTION  assignment sse               **/
/************************************************/
int assign_sse(int naa, SSEDAT *ssedat, RESDAT *resdat, int coil_as_strand) {
  int i, j;
  int iaa, jaa;
  int icoord;
  int nsse, isse;
  int nsse_ok, isse_ok;
  int sw_H, sw_E;
  int iaa_ini[MX_SSE];
  int iaa_end[MX_SSE];
  int SSEsw;
  int len_seg;
  int init, end;
  int dl;
  int *sse_iaa;
  int sse_coil, sse_helix, sse_strand;
  char sse_coil_char, sse_helix_char, sse_strand_char;
  float d;
  float param_H1[4 + 1], param_H2;
  float param_E1[4 + 1], param_E2;
  float **dist_ca;

  /********************************/
  /**          calloc            **/
  /********************************/
  sse_iaa = (int *)calloc((size_t)(naa + 1), sizeof(int));
  dist_ca = (float **)calloc((size_t)(naa + 1), sizeof(float *));
  for (iaa = 0; iaa <= naa; iaa++) {
    dist_ca[iaa] = (float *)calloc((size_t)(4 + 1), sizeof(float));
  }

  /************************************************************************
   * parameters : 
   *
   *  Zhang Y. and Skolnick J., TM-align: a protein structure alignment
   *  algorithm based on the TM-score, Nucleic Acid Research, 2005, Vol.33,
   *  No.7
   ************************************************************************/
  param_H1[2] = 5.45f;
  param_H1[3] = 5.18f;
  param_H1[4] = 6.37f;
  param_H2    = 2.10f;
  param_E1[2] = 6.10f;
  param_E1[3] = 10.40f;
  param_E1[4] = 13.00f;
  param_E2    = 1.60f; // original=1.42 ( Zhang.et.al )

  /** let all coils be strand in the case there is too few SSSs **/
  if (coil_as_strand == TRUE) {
    sse_coil = SSE_STRAND;
    sse_coil_char = 'C';
  }
  /** default vehavior **/
  else{
    sse_coil = SSE_COIL;
    sse_coil_char = 'C';
  }
  sse_helix = SSE_HELIX_;
  sse_helix_char = 'H';
  sse_strand = SSE_STRAND;
  sse_strand_char = 'E';
  /***************************************/
  /**       initial calcuration         **/
  /***************************************/
  /**  initialize  **/
  for (iaa = 1; iaa <= naa; iaa++) {
    resdat[iaa].ssechar = sse_coil_char;
    resdat[iaa].ssetype = sse_coil;
    sse_iaa[iaa] = sse_coil;
  }

  /**  calc distance  **/
  for (iaa = 1; iaa <= naa - 4; iaa++) {
    for (i = 2; i <= 4; i++) {
      for (icoord = 1; icoord <= Ncoord; icoord++) {
        dist_ca[iaa][i] +=
            (resdat[iaa].cacoord[icoord] - resdat[iaa + i].cacoord[icoord]) *
            (resdat[iaa].cacoord[icoord] - resdat[iaa + i].cacoord[icoord]);
      }
      dist_ca[iaa][i] = (float)sqrt(dist_ca[iaa][i]);
    }
  }

  /**************************************/
  /**            assignment            **/
  /**************************************/
  for (iaa = 3; iaa <= naa - 3; iaa++) {
    /** initialize **/
    sw_H = ON;
    sw_E = ON;
    /** assign HELIX **/
    for (jaa = iaa - 2; jaa <= iaa - 1; jaa++) {
      for (i = 2; i <= 4; i++) {
        d = (float)fabs(dist_ca[jaa][i] - param_H1[i]);
        if (d > param_H2) {
          sw_H = OFF;
        }
      }
    }
    /** assign STRAND **/
    for (jaa = iaa - 2; jaa <= iaa - 1; jaa++) {
      for (i = 2; i <= 4; i++) {
        d = (float)fabs(dist_ca[jaa][i] - param_E1[i]);
        if (d > param_E2) {
          sw_E = OFF;
        }
      }
    }

    if (sw_H == ON) {
      /****  assign HELIX  ****/
      sse_iaa[iaa - 1] = sse_helix;
      sse_iaa[iaa] = sse_helix;
      sse_iaa[iaa + 1] = sse_helix;
      resdat[iaa - 1].ssetype = sse_helix;
      resdat[iaa].ssetype = sse_helix;
      resdat[iaa + 1].ssetype = sse_helix;
      resdat[iaa - 1].ssechar = sse_helix_char;
      resdat[iaa].ssechar = sse_helix_char;
      resdat[iaa + 1].ssechar = sse_helix_char;
    }
    if (sw_E == ON) {
      /****  assign STRAND  ****/
      sse_iaa[iaa - 1] = sse_strand;
      sse_iaa[iaa] = sse_strand;
      sse_iaa[iaa + 1] = sse_strand;
      resdat[iaa - 1].ssetype = sse_strand;
      resdat[iaa].ssetype = sse_strand;
      resdat[iaa + 1].ssetype = sse_strand;
      resdat[iaa - 1].ssechar = sse_strand_char;
      resdat[iaa].ssechar = sse_strand_char;
      resdat[iaa + 1].ssechar = sse_strand_char;
    }
  }

  /***************************/
  /**  define SSE segments  **/
  /***************************/
  /****  find SSE start & end  ****/
  isse = 0;
  SSEsw = OFF;
  for (iaa = 1; iaa <= naa; iaa++) {
    if (sse_iaa[iaa] == sse_strand || sse_iaa[iaa] == sse_helix) {
      if (SSEsw == OFF) {
        isse++;
        iaa_ini[isse] = iaa;
      } else {
        if (sse_iaa[iaa] != sse_iaa[iaa - 1]) {
          iaa_end[isse] = iaa - 1;
          isse++;
          iaa_ini[isse] = iaa;
        } else {
          iaa_end[isse] = iaa - 1;
        }
      }
      SSEsw = ON;
    } else {
      if (SSEsw == ON) {
        iaa_end[isse] = iaa - 1;
      }
      SSEsw = OFF;
    }
  }
  /*  for last aa  */
  if (sse_iaa[naa] == sse_strand || sse_iaa[naa] == sse_helix) {
    iaa_end[isse] = naa;
  }
  nsse = isse;

  /****  check SSE length  ****/
  isse_ok = 0;
  for (isse = 1; isse <= nsse; isse++) {
    init = iaa_ini[isse];
    end = iaa_end[isse];
    len_seg = 0;
    if (sse_iaa[init] == sse_strand) {
      len_seg = seg_STRAND;
    }
    if (sse_iaa[init] == sse_helix) {
      len_seg = seg_HELIX_;
    }
    dl = end - init - len_seg + 2;

    if (1 <= dl) {
      for (i = 1; i <= dl; i++) {
        isse_ok++;
        ssedat[isse_ok].ssetype = sse_iaa[init];
        if (ssedat[isse_ok].ssetype == sse_strand) {
          ssedat[isse_ok].ssechar = sse_helix_char;
        } else if (ssedat[isse_ok].ssetype == sse_helix) {
          ssedat[isse_ok].ssechar = sse_strand_char;
        } else {
          ssedat[isse_ok].ssetype = sse_coil;
        }
        ssedat[isse_ok].naa = len_seg;
        for (j = 1; j <= len_seg; j++) {
          ssedat[isse_ok].iaa[j] = init + i + j - 2;
        }
      }
    }
  }
  nsse_ok = isse_ok;

  /****************************/
  /**          free          **/
  /****************************/
  free(sse_iaa);
  for (iaa = 0; iaa <= naa; iaa++) {
    free(dist_ca[iaa]);
  }
  free(dist_ca);

  return (nsse_ok);
}
