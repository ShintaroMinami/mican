#include "mican.h"
#include "main.h"

#define USED -1000.0
#define ALIGNED -2000.0

PRIVATE void check_close(int naa1, int naa2, int *ncons, int **id_iaa,
                         float **mat, int mode);
PRIVATE void takeover_mat(float **mat1, float **mat2, int ncons, int **id_iaa);
PRIVATE void loc_ali_iterate(int naa1, int naa2, float **mat, float **score,
                             float **score_r, float **length, float **length_r,
			     int *align, int *segmode, int ncons, int **fw_iaa, int **rv_iaa,
			     RESDAT *resdat_q, RESDAT *resdat_t);

void global_alignment(int naa1, int naa2, float ***mat, int *align, int *segmode,
                      float **score, float **score_r, RESDAT *resdat_q,
                      RESDAT *resdat_t) {
  int i;
  int iaa, jaa;
  int icons, ncons;
  int **fw_iaa, **rv_iaa;
  float **length, **length_r;

  fw_iaa = (int **)calloc(2, sizeof(int *));
  rv_iaa = (int **)calloc(2, sizeof(int *));
  for (i = 0; i <= 1; i++) {
    fw_iaa[i] = (int *)calloc((size_t)(naa1 * naa2 + 1), sizeof(int));
    rv_iaa[i] = (int *)calloc((size_t)(naa1 * naa2 + 1), sizeof(int));
  }
  length = (float **)calloc((size_t)(naa1 + 2), sizeof(float *));
  length_r = (float **)calloc((size_t)(naa1 + 2), sizeof(float *));
  for(iaa = 0; iaa <= naa1 + 1; iaa++){
    length[iaa] = (float *)calloc((size_t)(naa2 + 2), sizeof(float));
    length_r[iaa] = (float *)calloc((size_t)(naa2 + 2), sizeof(float));
  }

  /**  check close atom  **/
  check_close(naa1, naa2, &ncons, fw_iaa, mat[Nite_loc], FORWARD);
  check_close(naa1, naa2, &ncons, rv_iaa, mat[Nite_loc], REVERSE);

  /**********************************/
  /**         initialize           **/
  /**********************************/
  for (icons = 1; icons <= ncons; icons++) {
    iaa = fw_iaa[0][icons];
    jaa = fw_iaa[1][icons];
    score[iaa][jaa] = 0;
    score_r[iaa][jaa] = 0;
    length[iaa][jaa] = 0;
    length_r[iaa][jaa] = 0;
  }

  /**************************/
  /**      alignment       **/
  /**************************/
  for (i = 1; i <= Nite_loc; i++) {
    /**  recalc mat  **/
    takeover_mat(mat[0], mat[i], ncons, fw_iaa);
    /**  get segments  **/
    loc_ali_iterate(naa1, naa2, mat[0], score, score_r, length, length_r,
		    align, segmode, ncons, fw_iaa, rv_iaa, resdat_q, resdat_t);
  }

  /**********************/
  /**       free       **/
  /**********************/
  for (i = 0; i <= 1; i++) {
    free(fw_iaa[i]);
    free(rv_iaa[i]);
  }
  free(fw_iaa);
  free(rv_iaa);

  for(iaa = 0; iaa <= naa1 + 1; iaa++){
    free(length[iaa]);
    free(length_r[iaa]);
  }
  free(length);
  free(length_r);

  return;
} /* end of function */

/**  function : takeover matrix  **/
PRIVATE void takeover_mat(float **mat1, float **mat2, int ncons, int **id_iaa) {
  int icons, iaa, jaa;
  for (icons = 1; icons <= ncons; icons++) {
    iaa = id_iaa[0][icons];
    jaa = id_iaa[1][icons];
    if (mat1[iaa][jaa] == USED) {
      continue;
    }
    mat1[iaa][jaa] = mat2[iaa][jaa];
  }
  return;
}

/**  check close pair a.a.  **/
PRIVATE void check_close(int naa1, int naa2, int *ncons, int **id_iaa,
                         float **mat, int mode) {
  int iaa, jaa;

  *ncons = 0;
  if (mode == FORWARD) {
    for (iaa = 1; iaa <= naa1; iaa++) {
      for (jaa = 1; jaa <= naa2; jaa++) {
        if (mat[iaa][jaa] > 0) {
          *ncons = *ncons + 1;
          id_iaa[0][*ncons] = iaa;
          id_iaa[1][*ncons] = jaa;
        }
      }
    }
  } else {
    for (iaa = 1; iaa <= naa1; iaa++) {
      for (jaa = naa2; jaa >= 1; jaa--) {
        if (mat[iaa][jaa] > 0) {
          *ncons = *ncons + 1;
          id_iaa[0][*ncons] = iaa;
          id_iaa[1][*ncons] = jaa;
        }
      }
    }
  }

  return;
}

PRIVATE void get_alignment(int iaa, int jaa, float **mat, int *align, int *segmode, int mode);
PRIVATE void recalc_mat(int naa1, int naa2, float **mat, int iaa_max,
                        int jaa_max, int mode);
PRIVATE void recalc_mat_seq(int naa1, int naa2, float **mat, int iaa_max,
                            int jaa_max, RESDAT *resdat_q, RESDAT *resdat_t,
                            int mode);
PRIVATE void calc_score(int ncons, float **score, float **length, float **mat,
			int **id_iaa, int *iaa_max, int *jaa_max, int mode);

PRIVATE void loc_ali_iterate(int naa1, int naa2, float **mat, float **score,
                             float **score_r, float **length, float **length_r,
			     int *align, int *segmode, int ncons, int **fw_iaa, int **rv_iaa,
			     RESDAT *resdat_q, RESDAT *resdat_t) {
  int iaa_max, jaa_max;
  int iaa_max_r, jaa_max_r;
  int mode;
  float score_max, score_max_r;
  float length_max, length_max_r;
  float cutoff;

  cutoff = input.smin;

  iaa_max = 0;
  jaa_max = 0;
  iaa_max_r = 0;
  jaa_max_r = 0;

  /*********************************/
  /**    find local alignment     **/
  /*********************************/
  /**  calc score  **/
  if (input.mode == FWandRV || input.mode == SEQUENT || input.mode == FORWARD) {
    calc_score(ncons, score, length, mat, fw_iaa, &iaa_max, &jaa_max, FORWARD);
  }
  if (input.mode == FWandRV || input.mode == REVERSE) {
    calc_score(ncons, score_r, length_r, mat, rv_iaa, &iaa_max_r, &jaa_max_r, REVERSE);
  }

  if (input.mode == FORWARD) {
    /**  forward only mode  **/
    score_max = score[iaa_max][jaa_max];
    length_max = length[iaa_max][jaa_max];
    while (score_max >= cutoff && length_max >= input.lmin) {
      get_alignment(iaa_max, jaa_max, mat, align, segmode, FORWARD);
      recalc_mat(naa1, naa2, mat, iaa_max, jaa_max, FORWARD);
      calc_score(ncons, score, length, mat, fw_iaa, &iaa_max, &jaa_max, FORWARD);
      score_max = score[iaa_max][jaa_max];
      length_max = length[iaa_max][jaa_max];
    }
  } else if (input.mode == SEQUENT) {
    score_max = score[iaa_max][jaa_max];
    length_max = length[iaa_max][jaa_max];
    while (score_max >= cutoff && length_max >= input.lmin) {
      get_alignment(iaa_max, jaa_max, mat, align, segmode, FORWARD);
      recalc_mat_seq(naa1, naa2, mat, iaa_max, jaa_max, resdat_q, resdat_t,
                     FORWARD);
      calc_score(ncons, score, length, mat, fw_iaa, &iaa_max, &jaa_max, FORWARD);
      score_max = score[iaa_max][jaa_max];
      length_max = length[iaa_max][jaa_max];
    }
  } else if (input.mode == REVERSE) {
    /**  reverse only mode  **/
    score_max = score_r[iaa_max_r][jaa_max_r];
    length_max = length_r[iaa_max_r][jaa_max_r];
    while (score_max >= cutoff && length_max >= input.lmin) {
      get_alignment(iaa_max_r, jaa_max_r, mat, align, segmode, REVERSE);
      recalc_mat(naa1, naa2, mat, iaa_max_r, jaa_max_r, REVERSE);
      calc_score(ncons, score_r, length_r, mat, rv_iaa, &iaa_max_r, &jaa_max_r, REVERSE);
      score_max = score_r[iaa_max_r][jaa_max_r];
      length_max = length_r[iaa_max_r][jaa_max_r];
    }

  } else {
    /**  forward & reverse mode  **/
    score_max = score[iaa_max][jaa_max];
    score_max_r = score_r[iaa_max_r][jaa_max_r];
    length_max = length[iaa_max][jaa_max];
    length_max_r = length_r[iaa_max_r][jaa_max_r];
    mode = FORWARD;
    if (score_max < score_max_r) {
      mode = REVERSE;
      score_max = score_max_r;
      length_max = length_max_r;
    }
    while (score_max >= cutoff && length_max >= input.lmin) {
      if (mode == FORWARD) {
        get_alignment(iaa_max, jaa_max, mat, align, segmode, FORWARD);
        recalc_mat(naa1, naa2, mat, iaa_max, jaa_max, FORWARD);
      } else {
        get_alignment(iaa_max_r, jaa_max_r, mat, align, segmode, REVERSE);
        recalc_mat(naa1, naa2, mat, iaa_max_r, jaa_max_r, REVERSE);
      }
      calc_score(ncons, score, length, mat, fw_iaa, &iaa_max, &jaa_max, FORWARD);
      calc_score(ncons, score_r, length_r, mat, rv_iaa, &iaa_max_r, &jaa_max_r, REVERSE);
      score_max = score[iaa_max][jaa_max];
      score_max_r = score_r[iaa_max_r][jaa_max_r];
      length_max = length[iaa_max][jaa_max];
      length_max_r = length_r[iaa_max_r][jaa_max_r];
      mode = FORWARD;
      if (score_max < score_max_r) {
        mode = REVERSE;
        score_max = score_max_r;
	length_max = length_max_r;
      }
    }
  }

  return;
}

/**  function : calc_score  **/
PRIVATE void calc_score(int ncons, float **score, float **length, float **mat,
			int **id_iaa, int *iaa_max, int *jaa_max, int mode) {
  int icons;
  int iaa, jaa;
  float max_score;

  max_score = 0;
  for (icons = 1; icons <= ncons; icons++) {
    iaa = id_iaa[0][icons];
    jaa = id_iaa[1][icons];
    if (mat[iaa][jaa] < 0) {
      score[iaa][jaa] = 0;
      length[iaa][jaa] = 0;
      continue;
    }

    if (mode == FORWARD) {
      score[iaa][jaa] = score[iaa - 1][jaa - 1] + mat[iaa][jaa];
      length[iaa][jaa] = length[iaa - 1][jaa - 1] + 1;
    }
    if (mode == REVERSE) {
      score[iaa][jaa] = score[iaa - 1][jaa + 1] + mat[iaa][jaa];
      length[iaa][jaa] = length[iaa - 1][jaa + 1] + 1;
    }

    if (score[iaa][jaa] > max_score && length[iaa][jaa] >= input.lmin) {
      max_score = score[iaa][jaa];
      *iaa_max = iaa;
      *jaa_max = jaa;
    }
  }

  return;
}

/**  function : get alignment  **/
PRIVATE void get_alignment(int iaa, int jaa, float **mat, int *align, int *segmode,
                           int mode) {
  float sc;

  sc = mat[iaa][jaa];
  while (sc > 0) {
    align[iaa] = jaa;
    if (mode == FORWARD) {
      segmode[iaa] = FORWARD;
      jaa--;
    } else {
      segmode[iaa] = REVERSE;
      jaa++;
    }
    iaa--;
    sc = mat[iaa][jaa];
  }
  return;
}

/**  function : recalc matrix  **/
PRIVATE void recalc_mat(int naa1, int naa2, float **mat, int iaa_max,
                        int jaa_max, int mode) {
  int iaa, jaa;
  int i, j;
  float sc;

  iaa = iaa_max;
  jaa = jaa_max;
  sc = mat[iaa][jaa];
  while (sc > 0) {
    for (i = 1; i <= naa1; i++) {
      mat[i][jaa] = USED;
    }
    for (j = 1; j <= naa2; j++) {
      mat[iaa][j] = USED;
    }
    mat[iaa][jaa] = ALIGNED;
    iaa--;
    if (mode == FORWARD) {
      jaa--;
    } else {
      jaa++;
    }
    sc = mat[iaa][jaa];
  }
  return;
}

PRIVATE void recalc_mat_seq(int naa1, int naa2, float **mat, int iaa_max,
                            int jaa_max, RESDAT *resdat_q, RESDAT *resdat_t,
                            int mode) {
  int iaa, jaa;
  int i, j;
  float sc;
  int iaa_min, jaa_min;
  int chainID_q;
  int chainID_t;
  int tempID_q;
  int tempID_t;

  iaa = iaa_max;
  jaa = jaa_max;
  sc = mat[iaa][jaa];
  while (sc >= 0) {
    for (i = 1; i <= naa1; i++) {
      mat[i][jaa] = USED;
    }
    for (j = 1; j <= naa2; j++) {
      mat[iaa][j] = USED;
    }
    mat[iaa][jaa] = ALIGNED;
    iaa--;
    if (mode == FORWARD) {
      jaa--;
    } else {
      jaa++;
    }
    sc = mat[iaa][jaa];
  }

  iaa_min = iaa + 1;
  jaa_min = jaa + 1;

  jaa = jaa_min;
  tempID_q = 0;
  tempID_t = 0;

  for (iaa = iaa_min; iaa <= iaa_max; iaa++) {

    chainID_q = resdat_q[iaa].chainID;
    chainID_t = resdat_t[jaa].chainID;
    if (tempID_q != chainID_q || tempID_t != chainID_t) {
      for (i = 1; i <= iaa_min; i++) {
        for (j = jaa_max; j <= naa2; j++) {
          if (mat[i][j] != USED) {
            if (chainID_q == resdat_q[i].chainID &&
                chainID_t == resdat_t[j].chainID) {
              mat[i][j] = USED;
            }
          }
        }
      }
      for (i = iaa_max; i <= naa1; i++) {
        for (j = 1; j <= jaa_min; j++) {
          if (mat[i][j] != USED) {
            if (chainID_q == resdat_q[i].chainID &&
                chainID_t == resdat_t[j].chainID) {
              mat[i][j] = USED;
            }
          }
        }
      }
      tempID_q = chainID_q;
      tempID_t = chainID_t;
    }
    jaa++;
  }

  return;
}
