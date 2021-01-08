#include "mican.h"

typedef struct {
  int iali, jali;
  float score;
} PAIR;

PRIVATE float similarity(int naa_q, int naa_t, ALIGN *align1, ALIGN *align2);
/** functions for qsort() **/
PRIVATE int compare_align(const void *a, const void *b);
PRIVATE int compare_scorecmp(const void *a, const void *b);

void clustering_align(int naa_q, int naa_t, int nalign, ALIGN *align,
                      float cutoff) {
  int ialign, jalign;
  int icount, ncount;

  PAIR *pair;
  pair = (PAIR *)calloc((size_t)(nalign * (nalign + 1) / 2), sizeof(PAIR));

  /******   calc similarity   ******/
  ncount = -1;
  for (ialign = 0; ialign < nalign; ialign++) {
    if (align[ialign].check == FALSE) {
      continue;
    }
    for (jalign = ialign + 1; jalign <= nalign; jalign++) {
      if (align[jalign].check == FALSE) {
        continue;
      }
      ncount++;
      pair[ncount].iali = ialign;
      pair[ncount].jali = jalign;
      pair[ncount].score =
          similarity(naa_q, naa_t, &(align[ialign]), &(align[jalign]));
    }
  }

  /******   quick sort   *********/
  qsort(pair, (size_t)(ncount), sizeof(PAIR), compare_align);

  /******   clustering   ******/
  for (icount = 0; icount <= ncount; icount++) {
    if (pair[icount].score < cutoff) {
      break;
    }
    /****  for similar alignment  ****/
    ialign = pair[icount].iali;
    jalign = pair[icount].jali;

    if (align[ialign].score_cmp >= align[jalign].score_cmp) {
      align[jalign].score_cmp = 0;
      align[jalign].check = FALSE;
    } else {
      align[ialign].score_cmp = 0;
      align[ialign].check = FALSE;
    }
  }

  /******   quick sort by structural similarity score   ******/
  qsort(align, (size_t)(nalign), sizeof(ALIGN), compare_scorecmp);

  /******   free   ******/
  free(pair);

  return;
}

/****** function for qsort() alignment ******/
PRIVATE int compare_align(const void *a, const void *b) {
  PAIR *x, *y;
  float delta;
  x = (PAIR *)a;
  y = (PAIR *)b;
  delta = (*y).score - (*x).score;

  if (delta < 0) {
    return -1;
  }
  if (delta > 0) {
    return 1;
  }
  return 0;
}

/****** function for qsort() TMscore ******/
PRIVATE int compare_scorecmp(const void *a, const void *b) {
  ALIGN *x, *y;
  float delta;
  x = (ALIGN *)a;
  y = (ALIGN *)b;
  delta = (*y).score_cmp - (*x).score_cmp;

  if (delta < 0) {
    return -1;
  }
  if (delta > 0) {
    return 1;
  }
  return 0;
}

/*********************************************************/
PRIVATE float similarity(int naa_q, int naa_t, ALIGN *align1, ALIGN *align2) {
  int iaa;
  float calign, score1, score2;
  int delta;

  calign = 0;
  for (iaa = 1; iaa <= naa_q; iaa++) {
    if (align1->align_q[iaa] == 0) {
      continue;
    }
    delta = (int)abs(align1->align_q[iaa] - align2->align_q[iaa]);
    if (delta <= 4) {
      calign += 1.0f;
    } else if (delta <= 6) {
      calign += 0.8f;
    }
  }
  score1 = calign / (float)align1->naa_align;

  calign = 0;
  for (iaa = 1; iaa <= naa_t; iaa++) {
    if (align1->align_t[iaa] == 0) {
      continue;
    }
    delta = (int)abs(align1->align_t[iaa] - align2->align_t[iaa]);
    if (delta <= 4) {
      calign += 1.0f;
    } else if (delta <= 6) {
      calign += 0.8f;
    }
  }
  score2 = calign / (float)align1->naa_align;

  return ((score1 + score2) / 2);
}
