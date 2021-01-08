#include "mican.h"

int hres2trans(SSEDAT *ssedat_q, SSEDAT *ssedat_t, ALIGN align) {
  int i, j;
  int isse_q, isse_t;
  float xq[4], yq[4], zq[4];
  float xt[4], yt[4], zt[4];

  isse_q = align.isse_q;
  isse_t = align.isse_t;

  /****  coordinate difinition  ****/
  coord_definition(xq, yq, zq, ssedat_q, isse_q);
  if (align.coord_mode == FORWARD) {
    coord_definition(xt, yt, zt, ssedat_t, isse_t);
  } else {
    coord_definition_r(xt, yt, zt, ssedat_t, isse_t);
  }

  /****  rotation matrix  ****/
  for (i = 1; i <= 3; i++) {
    for (j = 1; j <= 3; j++) {
      align.rot[i][j] = xq[i] * xt[j] + yq[i] * yt[j] + zq[i] * zt[j];
    }
  }
  /****  translation vector  ****/
  for (i = 1; i <= 3; i++) {
    align.vec[i] = ssedat_q[isse_q].mid[i];
    for (j = 1; j <= 3; j++) {
      align.vec[i] += -align.rot[i][j] * ssedat_t[isse_t].mid[j];
    }
  }

  return (TRUE);
}
