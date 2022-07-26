#include "mican.h"
#include "main.h"

/** function for comparison of hash results **/
PRIVATE int compare(const void *a, const void *b);

void score_calculation(int naa_q, int naa_t, int nalign, ALIGN *align,
                       RESDAT *resdat_q, RESDAT *resdat_t) {
  int iaa;
  int ialign, nalign_true;

  nalign_true = 0;

  for (ialign = 0; ialign < nalign; ialign++) {
    align[ialign].naa_align = 0;
    align[ialign].TMscore_q = 0.0;
    align[ialign].TMscore_t = 0.0;
    align[ialign].TMscore_mod = 0.0;
    align[ialign].TMscore_ang = 0.0;
    align[ialign].SPscore = 0.0;
    align[ialign].rmsd = 0.0;
    align[ialign].score_cmp = 0.0;

    if (align[ialign].check == FALSE) {
      continue;
    }

    nalign_true++;

    for (iaa = 1; iaa <= naa_q; iaa++) {
      if (align[ialign].align_q[iaa] > 0) {
        align[ialign].naa_align++;
      }
    }

    /** TM-score query **/
    align[ialign].TMscore_q =
      TMscore(naa_q, naa_q, resdat_q, resdat_t, align[ialign]);
    /** TM-score template **/
    align[ialign].TMscore_t =
      TMscore(naa_t, naa_q, resdat_q, resdat_t, align[ialign]);
    /** modified TM-score **/
    align[ialign].TMscore_mod =
      TMscore_mod(naa_q, naa_q, naa_t, resdat_q, resdat_t, align[ialign], OFF, OFF);
    /** modified TM-score local **/
    align[ialign].TMloc_mod =
      TMscore_mod(align[ialign].naa_align, naa_q, naa_t, resdat_q, resdat_t, align[ialign], OFF, ON);
    /** modified TM-score + angle **/
    align[ialign].TMscore_ang =
      TMscore_mod(naa_q, naa_q, naa_t, resdat_q, resdat_t, align[ialign], ON, OFF);
    /** SP-score **/
    align[ialign].SPscore =
        SPscore(naa_q, naa_t, resdat_q, resdat_t, align[ialign]);
    /** RMSD **/
    align[ialign].rmsd = rmsd(naa_q, resdat_q, resdat_t, align[ialign]);

    /** coverage **/
    align[ialign].cover_q =
        100 * (float)(align[ialign].naa_align) / (float)(naa_q);
    align[ialign].cover_t =
        100 * (float)(align[ialign].naa_align) / (float)(naa_t);

    /**  Comparing score  **/
    if (input.type_select == SELE_SPscore) {
      align[ialign].score_cmp = align[ialign].SPscore;
    } else if (input.type_select == SELE_TMscore) {
      align[ialign].score_cmp = align[ialign].TMscore_q;
    } else if (input.type_select == SELE_aTMscore) {
      align[ialign].score_cmp = align[ialign].TMscore_ang;
    } else if (input.type_select == SELE_Daliscore) {
      align[ialign].score_cmp = align[ialign].Daliscore;
    } else {
      align[ialign].score_cmp = align[ialign].TMscore_mod;
    }
  }

  if (nalign_true < 1) {
    printf("No alignment results\n");
    exit(0);
  }

  /**  Quick Sort  **/
  if (nalign_true > 1) {
    qsort(align, (size_t)(nalign + 1), sizeof(ALIGN), compare);
  }

  return;
}

/****** function for comparison ******/
PRIVATE int compare(const void *a, const void *b) {
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



void score_calc_chain(int naa_q, int naa_t, ALIGN align, ALIGN **align_chain,
                      RESDAT *resdat_q, RESDAT *resdat_t, PDBDAT pdbdat_q,
                      PDBDAT pdbdat_t) {
  int iaa, jaa;
  int nchain_q = pdbdat_q.nchain;
  int nchain_t = pdbdat_t.nchain;
  int ichain, jchain;


  for (ichain = 1; ichain <= nchain_q; ichain++) {
    for (jchain = 1; jchain <= nchain_t; jchain++) {
      align_chain[ichain][jchain].naa_align = 0;
      for (iaa = 1; iaa <= naa_q; iaa++) {
        if (align.align_q[iaa] > 0) {
          jaa = align.align_q[iaa];
          if (resdat_q[iaa].chainID == pdbdat_q.chainID[ichain] &&
              resdat_t[jaa].chainID == pdbdat_t.chainID[jchain]) {
            align_chain[ichain][jchain].naa_align++;
          }
        }
      }
      align_chain[ichain][jchain].cover_q =
          100 * (float)(align_chain[ichain][jchain].naa_align) /
          (float)(pdbdat_q.naa[ichain]);
      align_chain[ichain][jchain].cover_t =
          100 * (float)(align_chain[ichain][jchain].naa_align) /
          (float)(pdbdat_t.naa[jchain]);
    }
  }

  jchain = 0;
  for (ichain = 1; ichain <= nchain_q; ichain++) {
    align_chain[ichain][jchain].naa_align = 0;
    for (iaa = 1; iaa <= naa_q; iaa++) {
      if (align.align_q[iaa] > 0) {
        jaa = align.align_q[iaa];
        if (resdat_q[iaa].chainID == pdbdat_q.chainID[ichain]) {
          align_chain[ichain][jchain].naa_align++;
        }
      }
    }
    align_chain[ichain][jchain].cover_q =
        100 * (float)(align_chain[ichain][jchain].naa_align) /
        (float)(pdbdat_q.naa[ichain]);
    align_chain[ichain][jchain].cover_t =
        100 * (float)(align_chain[ichain][jchain].naa_align) /
        (float)(pdbdat_t.naa[jchain]);
  }

  ichain = 0;
  for (jchain = 1; jchain <= nchain_t; jchain++) {
    align_chain[ichain][jchain].naa_align = 0;
    for (iaa = 1; iaa <= naa_q; iaa++) {
      if (align.align_q[iaa] > 0) {
        jaa = align.align_q[iaa];
        if (resdat_t[jaa].chainID == pdbdat_t.chainID[jchain]) {
          align_chain[ichain][jchain].naa_align++;
        }
      }
    }
    align_chain[ichain][jchain].cover_q =
        100 * (float)(align_chain[ichain][jchain].naa_align) /
        (float)(pdbdat_q.naa[ichain]);
    align_chain[ichain][jchain].cover_t =
        100 * (float)(align_chain[ichain][jchain].naa_align) /
        (float)(pdbdat_t.naa[jchain]);
  }

  return;
}




void calc_seq_identity(int naa_q, int naa_t, int nalign, ALIGN *align,
                       RESDAT *resdat_q, RESDAT *resdat_t) {

  int ialign;
  int iaa, jaa;
  int count;

  for (ialign = 0; ialign <= nalign; ialign++) {
    if (align[ialign].check == FALSE) {
      continue;
    }
    count = 0;
    for (iaa = 1; iaa <= naa_q; iaa++) {
      if (align[ialign].align_q[iaa] > 0) {
        jaa = align[ialign].align_q[iaa];
        if (strcmp(resdat_q[iaa].resname, resdat_t[jaa].resname) == 0) {
          count++;
        }
      }
    }
    align[ialign].seqID = 100 * (float)count / (float)align[ialign].naa_align;
  }

  return;
}
