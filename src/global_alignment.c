#include "mican.h"
#include "main.h"

#define USED -1000.0
#define ALIGNED -2000.0
#define NOT_SET -100.0f

/*
 * Sparse implementation of the greedy local-segment alignment.
 *
 * Only "contact" cells (pairs whose loosest-level matrix value mat[Nite_loc]
 * is positive) can ever hold a positive value, so the dense (naa1 x naa2)
 * matrices of the original implementation are replaced by arrays indexed by
 * contact index c.  Contacts are given in forward order (iaa ascending, then
 * jaa ascending).  Every operation of the dense version is reproduced on the
 * contact cells in the same order, so results are bit-identical:
 *   - a non-contact cell always holds a value <= 0 and score/length 0, hence
 *     a missing diagonal predecessor points to a sentinel slot (index ncons)
 *     whose value is NOT_SET and whose score/length stay 0;
 *   - USED sweeps over a row/column only need to touch contact cells.
 *
 * Example: contacts (1,1) (1,2) (2,2) (2,3) give row_start = {_,0,2,4},
 *          fw_pred = {4,4,0,1} (4 = sentinel), rv_order = {1,0,3,2}.
 */

typedef struct {
  int ncons;
  const int *ci, *cj;  /* residue indices of contact c */
  const float *cv;     /* raw matrix value of contact c */
  int *row_start;      /* contacts of row i: [row_start[i], row_start[i+1]) */
  int *col_start;      /* contacts of column j: col_list[col_start[j] ..] */
  int *col_list;
  int *fw_pred;        /* contact index of (i-1, j-1), or ncons (sentinel) */
  int *rv_pred;        /* contact index of (i-1, j+1), or ncons (sentinel) */
  int *rv_order;       /* contacts in reverse order (i asc, j desc) */
  float *m0;           /* working matrix value */
  float *score, *score_r, *length, *length_r;
} CONTACTS;

PRIVATE void loc_ali_iterate(int naa1, int naa2, CONTACTS *cs, int *align,
                             int *segmode, RESDAT *resdat_q, RESDAT *resdat_t);

/* Grow a buffer to hold at least n elements of size sz. */
PRIVATE void *grow(void *p, int *cap, int n, size_t sz) {
  if (n <= *cap) {
    return p;
  }
  *cap = n + n / 2 + 16;
  p = realloc(p, (size_t)(*cap) * sz);
  if (p == NULL) {
    fprintf(stderr, "!! Error   : memory allocation failed\n");
    exit(1);
  }
  return p;
}

void global_alignment(int naa1, int naa2, int ncons, const int *ci,
                      const int *cj, const float *cv, const float *d0,
                      int *align, int *segmode, RESDAT *resdat_q,
                      RESDAT *resdat_t) {
  /* Work buffers are reused across calls. */
  static int *row_start, *col_start, *col_list, *fw_pred, *rv_pred, *rv_order;
  static int *prev_row;
  static float *buf;
  static int cap_row, cap_col, cap_prev, cap_cl, cap_fp, cap_rp, cap_ro, cap_buf;
  CONTACTS cs;
  int c, k, j, iaa, ite;

  row_start = grow(row_start, &cap_row, naa1 + 2, sizeof(int));
  col_start = grow(col_start, &cap_col, naa2 + 3, sizeof(int));
  prev_row = grow(prev_row, &cap_prev, naa2 + 2, sizeof(int));
  col_list = grow(col_list, &cap_cl, ncons + 1, sizeof(int));
  fw_pred = grow(fw_pred, &cap_fp, ncons + 1, sizeof(int));
  rv_pred = grow(rv_pred, &cap_rp, ncons + 1, sizeof(int));
  rv_order = grow(rv_order, &cap_ro, ncons + 1, sizeof(int));
  buf = grow(buf, &cap_buf, 5 * (ncons + 1), sizeof(float));

  cs.ncons = ncons;
  cs.ci = ci;
  cs.cj = cj;
  cs.cv = cv;
  cs.row_start = row_start;
  cs.col_start = col_start;
  cs.col_list = col_list;
  cs.fw_pred = fw_pred;
  cs.rv_pred = rv_pred;
  cs.rv_order = rv_order;
  cs.m0 = buf;
  cs.score = buf + (ncons + 1);
  cs.score_r = buf + 2 * (ncons + 1);
  cs.length = buf + 3 * (ncons + 1);
  cs.length_r = buf + 4 * (ncons + 1);

  /**  row ranges  **/
  c = 0;
  for (iaa = 0; iaa <= naa1 + 1; iaa++) {
    while (c < ncons && ci[c] < iaa) {
      c++;
    }
    row_start[iaa] = c;
  }

  /**  column lists  **/
  for (j = 0; j <= naa2 + 2; j++) {
    col_start[j] = 0;
  }
  for (c = 0; c < ncons; c++) {
    col_start[cj[c] + 1]++;
  }
  for (j = 1; j <= naa2 + 2; j++) {
    col_start[j] += col_start[j - 1];
  }
  for (c = 0; c < ncons; c++) {
    col_list[col_start[cj[c]]++] = c;
  }
  for (j = naa2 + 2; j >= 1; j--) {
    col_start[j] = col_start[j - 1];
  }
  col_start[0] = 0;

  /**  diagonal predecessors & reverse order  **/
  for (j = 0; j <= naa2 + 1; j++) {
    prev_row[j] = ncons;
  }
  k = 0;
  for (iaa = 1; iaa <= naa1; iaa++) {
    for (c = row_start[iaa]; c < row_start[iaa + 1]; c++) {
      fw_pred[c] = prev_row[cj[c] - 1];
      rv_pred[c] = prev_row[cj[c] + 1];
    }
    for (c = row_start[iaa - 1]; c < row_start[iaa]; c++) {
      prev_row[cj[c]] = ncons;
    }
    for (c = row_start[iaa]; c < row_start[iaa + 1]; c++) {
      prev_row[cj[c]] = c;
    }
    for (c = row_start[iaa + 1] - 1; c >= row_start[iaa]; c--) {
      rv_order[k++] = c;
    }
  }
  for (c = row_start[naa1]; c < row_start[naa1 + 1]; c++) {
    prev_row[cj[c]] = ncons;
  }

  /**  initialize (c == ncons is the sentinel)  **/
  for (c = 0; c <= ncons; c++) {
    cs.m0[c] = NOT_SET;
  }
  cs.score[ncons] = cs.score_r[ncons] = 0;
  cs.length[ncons] = cs.length_r[ncons] = 0;

  /**************************/
  /**      alignment       **/
  /**************************/
  for (ite = 1; ite <= Nite_loc; ite++) {
    /**  recalc mat (take over level-ite values)  **/
    for (c = 0; c < ncons; c++) {
      if (cs.m0[c] == USED) {
        continue;
      }
      cs.m0[c] = (cv[c] >= d0[ite]) ? cv[c] : NOT_SET;
    }
    /**  get segments  **/
    loc_ali_iterate(naa1, naa2, &cs, align, segmode, resdat_q, resdat_t);
  }
  return;
} /* end of function */

PRIVATE void get_alignment(CONTACTS *cs, int c, int *align, int *segmode,
                           int mode);
PRIVATE void recalc_mat(CONTACTS *cs, int c, int mode);
PRIVATE void recalc_mat_seq(int naa1, int naa2, CONTACTS *cs, int c,
                            RESDAT *resdat_q, RESDAT *resdat_t);
PRIVATE void calc_score(CONTACTS *cs, float *score, float *length, int *c_max,
                        int mode);

/* score/length of contact c; c == -1 stands for the (0,0) cell (always 0) */
#define SC(arr, c) ((c) >= 0 ? (arr)[c] : 0.0f)

PRIVATE void loc_ali_iterate(int naa1, int naa2, CONTACTS *cs, int *align,
                             int *segmode, RESDAT *resdat_q, RESDAT *resdat_t) {
  int c_max, c_max_r;
  int mode;
  float score_max, score_max_r;
  float length_max, length_max_r;
  float cutoff;

  cutoff = input.smin;

  c_max = -1;
  c_max_r = -1;

  /*********************************/
  /**    find local alignment     **/
  /*********************************/
  /**  calc score  **/
  if (input.mode == FWandRV || input.mode == SEQUENT || input.mode == FORWARD) {
    calc_score(cs, cs->score, cs->length, &c_max, FORWARD);
  }
  if (input.mode == FWandRV || input.mode == REVERSE) {
    calc_score(cs, cs->score_r, cs->length_r, &c_max_r, REVERSE);
  }

  if (input.mode == FORWARD) {
    /**  forward only mode  **/
    score_max = SC(cs->score, c_max);
    length_max = SC(cs->length, c_max);
    while (score_max >= cutoff && length_max >= input.lmin) {
      get_alignment(cs, c_max, align, segmode, FORWARD);
      recalc_mat(cs, c_max, FORWARD);
      calc_score(cs, cs->score, cs->length, &c_max, FORWARD);
      score_max = SC(cs->score, c_max);
      length_max = SC(cs->length, c_max);
    }
  } else if (input.mode == SEQUENT) {
    score_max = SC(cs->score, c_max);
    length_max = SC(cs->length, c_max);
    while (score_max >= cutoff && length_max >= input.lmin) {
      get_alignment(cs, c_max, align, segmode, FORWARD);
      recalc_mat_seq(naa1, naa2, cs, c_max, resdat_q, resdat_t);
      calc_score(cs, cs->score, cs->length, &c_max, FORWARD);
      score_max = SC(cs->score, c_max);
      length_max = SC(cs->length, c_max);
    }
  } else if (input.mode == REVERSE) {
    /**  reverse only mode  **/
    score_max = SC(cs->score_r, c_max_r);
    length_max = SC(cs->length_r, c_max_r);
    while (score_max >= cutoff && length_max >= input.lmin) {
      get_alignment(cs, c_max_r, align, segmode, REVERSE);
      recalc_mat(cs, c_max_r, REVERSE);
      calc_score(cs, cs->score_r, cs->length_r, &c_max_r, REVERSE);
      score_max = SC(cs->score_r, c_max_r);
      length_max = SC(cs->length_r, c_max_r);
    }

  } else {
    /**  forward & reverse mode  **/
    score_max = SC(cs->score, c_max);
    score_max_r = SC(cs->score_r, c_max_r);
    length_max = SC(cs->length, c_max);
    length_max_r = SC(cs->length_r, c_max_r);
    mode = FORWARD;
    if (score_max < score_max_r) {
      mode = REVERSE;
      score_max = score_max_r;
      length_max = length_max_r;
    }
    while (score_max >= cutoff && length_max >= input.lmin) {
      if (mode == FORWARD) {
        get_alignment(cs, c_max, align, segmode, FORWARD);
        recalc_mat(cs, c_max, FORWARD);
      } else {
        get_alignment(cs, c_max_r, align, segmode, REVERSE);
        recalc_mat(cs, c_max_r, REVERSE);
      }
      calc_score(cs, cs->score, cs->length, &c_max, FORWARD);
      calc_score(cs, cs->score_r, cs->length_r, &c_max_r, REVERSE);
      score_max = SC(cs->score, c_max);
      score_max_r = SC(cs->score_r, c_max_r);
      length_max = SC(cs->length, c_max);
      length_max_r = SC(cs->length_r, c_max_r);
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
/* Keeps the previous *c_max when no cell exceeds 0 (as the dense version). */
#define CALC_SCORE_CELL(c)                                \
  do {                                                    \
    if (m0[c] < 0) {                                      \
      score[c] = 0;                                       \
      length[c] = 0;                                      \
      break;                                              \
    }                                                     \
    p = pred[c];                                          \
    score[c] = score[p] + m0[c];                          \
    length[c] = length[p] + 1;                            \
    if (score[c] > max_score && length[c] >= lmin) {      \
      max_score = score[c];                               \
      *c_max = c;                                         \
    }                                                     \
  } while (0)

PRIVATE void calc_score(CONTACTS *cs, float *score, float *length, int *c_max,
                        int mode) {
  int k, c, p;
  int ncons = cs->ncons;
  const int *pred = (mode == FORWARD) ? cs->fw_pred : cs->rv_pred;
  const int *order = cs->rv_order;
  const float *m0 = cs->m0;
  float lmin = input.lmin;
  float max_score;

  max_score = 0;
  if (mode == FORWARD) {
    for (c = 0; c < ncons; c++) {
      CALC_SCORE_CELL(c);
    }
  } else {
    for (k = 0; k < ncons; k++) {
      c = order[k];
      CALC_SCORE_CELL(c);
    }
  }

  return;
}

/**  function : get alignment  **/
PRIVATE void get_alignment(CONTACTS *cs, int c, int *align, int *segmode,
                           int mode) {
  const int *pred = (mode == FORWARD) ? cs->fw_pred : cs->rv_pred;

  while (c >= 0 && cs->m0[c] > 0) {
    align[cs->ci[c]] = cs->cj[c];
    segmode[cs->ci[c]] = (mode == FORWARD) ? FORWARD : REVERSE;
    c = pred[c];
  }
  return;
}

/* Mark every contact in row iaa and column jaa as USED. */
PRIVATE void sweep_row_col(CONTACTS *cs, int iaa, int jaa) {
  int k;
  for (k = cs->col_start[jaa]; k < cs->col_start[jaa + 1]; k++) {
    cs->m0[cs->col_list[k]] = USED;
  }
  for (k = cs->row_start[iaa]; k < cs->row_start[iaa + 1]; k++) {
    cs->m0[k] = USED;
  }
}

/**  function : recalc matrix  **/
PRIVATE void recalc_mat(CONTACTS *cs, int c, int mode) {
  const int *pred = (mode == FORWARD) ? cs->fw_pred : cs->rv_pred;

  while (c >= 0 && cs->m0[c] > 0) {
    sweep_row_col(cs, cs->ci[c], cs->cj[c]);
    cs->m0[c] = ALIGNED;
    c = pred[c];
  }
  return;
}

/* Mark contacts in rows [r0, r1] and columns [j0, j1] whose chain pair is
   (chainID_q, chainID_t) as USED. */
PRIVATE void block_used(CONTACTS *cs, int r0, int r1, int j0, int j1,
                        int chainID_q, int chainID_t, RESDAT *resdat_q,
                        RESDAT *resdat_t) {
  int c, j;
  for (c = cs->row_start[r0]; c < cs->row_start[r1 + 1]; c++) {
    j = cs->cj[c];
    if (j < j0 || j > j1 || cs->m0[c] == USED) {
      continue;
    }
    if (chainID_q == resdat_q[cs->ci[c]].chainID &&
        chainID_t == resdat_t[j].chainID) {
      cs->m0[c] = USED;
    }
  }
}

PRIVATE void recalc_mat_seq(int naa1, int naa2, CONTACTS *cs, int c,
                            RESDAT *resdat_q, RESDAT *resdat_t) {
  int iaa, jaa;
  int iaa_max, jaa_max;
  int iaa_min, jaa_min;
  int chainID_q;
  int chainID_t;
  int tempID_q;
  int tempID_t;

  iaa_max = (c >= 0) ? cs->ci[c] : 0;
  jaa_max = (c >= 0) ? cs->cj[c] : 0;
  iaa_min = iaa_max + 1;
  jaa_min = jaa_max + 1;
  while (c >= 0 && cs->m0[c] >= 0) {
    iaa_min = cs->ci[c];
    jaa_min = cs->cj[c];
    sweep_row_col(cs, iaa_min, jaa_min);
    cs->m0[c] = ALIGNED;
    c = cs->fw_pred[c];
  }

  jaa = jaa_min;
  tempID_q = 0;
  tempID_t = 0;

  for (iaa = iaa_min; iaa <= iaa_max; iaa++) {

    chainID_q = resdat_q[iaa].chainID;
    chainID_t = resdat_t[jaa].chainID;
    if (tempID_q != chainID_q || tempID_t != chainID_t) {
      block_used(cs, 1, iaa_min, jaa_max, naa2, chainID_q, chainID_t,
                 resdat_q, resdat_t);
      block_used(cs, iaa_max, naa1, 1, jaa_min, chainID_q, chainID_t,
                 resdat_q, resdat_t);
      tempID_q = chainID_q;
      tempID_t = chainID_t;
    }
    jaa++;
  }

  return;
}
