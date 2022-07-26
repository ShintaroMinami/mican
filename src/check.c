#include "mican.h"
#include "main.h"

/** function for comparison of hash results **/
PRIVATE int compare(const void *a, const void *b);
PRIVATE int check_redundancy(ALIGN align1, ALIGN align2);

void check_GHresult(int maxalign, HASH_RESULT *hres, ALIGN *align,
                    SSEDAT *ssedat_q, SSEDAT *ssedat_t) {
  int ialign;
  int num;

  /**********************************/
  /**          Quick Sort          **/
  /**********************************/
  qsort(hres, (size_t)(maxalign + 1), sizeof(HASH_RESULT), compare);

  /***********************************/
  /**         Select Top Ns         **/
  /***********************************/
  num = -1;
  for (ialign = 0; ialign <= maxalign; ialign++) {
    if (hres[ialign].count < input.min_GHcount) {
      break;
    } else {
      num++;

      align[num].count = hres[ialign].count;
      align[num].isse_q = hres[ialign].isse_q;
      align[num].isse_t = hres[ialign].isse_t;
      align[num].coord_mode = hres[ialign].coord_mode;

      align[num].check = hres2trans(ssedat_q, ssedat_t, align[num]);

      if (num != 0) {
        if (check_redundancy(align[num - 1], align[num]) == FALSE) {
          align[num].check = FALSE;
          num--;
        }
      }

      if (num >= input.GHcutnum) {
        return;
      }
    }
  }

  if (num < 0) {
    printf("No alignment results.\n");
    exit(0);
  }

  return;
}

/****** function for comparison ******/
PRIVATE int compare(const void *a, const void *b) {
  HASH_RESULT *x, *y;
  float delta;
  x = (HASH_RESULT *)a;
  y = (HASH_RESULT *)b;
  delta = (*y).count - (*x).count;

  if (delta < 0) {
    return -1;
  }
  if (delta > 0) {
    return 1;
  }
  return 0;
}

/****** check redundancy ******/
PRIVATE int check_redundancy(ALIGN align1, ALIGN align2) {
  int i, j;
  float diff1, diff2;

  diff1 = 0;
  diff2 = 0;

  for (i = 1; i <= 3; i++) {
    for (j = 1; j <= 3; j++) {
      diff1 += (align1.rot[i][j] - align2.rot[i][j]) *
               (align1.rot[i][j] - align2.rot[i][j]);
    }
    diff2 += (align1.vec[i] - align2.vec[i]) * (align1.vec[i] - align2.vec[i]);
  }

  if (diff1 + diff2 < EPS) {
    return (FALSE);
  }
  return (TRUE);
}

/****** function for check alignment ******/
int check_align(int naa_q, int naa_t, ALIGN *align, RESDAT *resdat_q, RESDAT *resdat_t) {
  int iaa;
  float TM_ang;

  align->naa_align = 0;
  for (iaa = 1; iaa <= naa_q; iaa++) {
    if (align->align_q[iaa] > 0) align->naa_align++;
  }

  TM_ang = TMscore_mod(naa_q, naa_q, naa_t, resdat_q, resdat_t, *align, ON, OFF);

  /******************/
  /**    return    **/
  /******************/
  if (align->naa_align == 0) {
    return (FALSE);
  }
  if (TM_ang > align->TMscore_ang + EPS) {
    align->TMscore_ang = TM_ang;
    return (TRUE);
  }

  return (FALSE);
}
