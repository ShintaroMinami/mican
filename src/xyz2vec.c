#include "mican.h"

PRIVATE void calc_perp(float *n, float *v, float *a, float *b);

void xyz2vec(int nsse, SSEDAT *ssedat, RESDAT *resdat) {
  int i;
  int isse;
  int naa;
  int icoord;
  int init, end;
  float a[3 + 1];
  float b[3 + 1];
  float c[3 + 1];
  float sum;

  /**  initialize  **/
  for (i = 0; i <= 3; i++) {
    a[i] = 0;
    b[i] = 0;
    c[i] = 0;
  }

  for (isse = 1; isse <= nsse; isse++) {
    naa = ssedat[isse].naa;
    init = ssedat[isse].iaa[1];
    end = ssedat[isse].iaa[naa];
    if (ssedat[isse].ssetype == SSE_HELIX_) {
      /**  sse start point for HELIX  **/
      for (icoord = 1; icoord <= Ncoord; icoord++) {
        a[icoord] = (0.74f * resdat[init].cacoord[icoord] +
                     resdat[init + 1].cacoord[icoord] +
                     resdat[init + 2].cacoord[icoord] +
                     0.74f * resdat[init + 3].cacoord[icoord]) /
                    3.48f;
      }
      /**  sse end point for HELIX  **/
      for (icoord = 1; icoord <= Ncoord; icoord++) {
        b[icoord] =
            (0.74f * resdat[end - 3].cacoord[icoord] +
             resdat[end - 2].cacoord[icoord] + resdat[end - 1].cacoord[icoord] +
             0.74f * resdat[end].cacoord[icoord]) /
            3.48f;
      }
      /**  sse perp point for HELIX  **/
      for (icoord = 1; icoord <= Ncoord; icoord++) {
        c[icoord] = (resdat[init + 2].cacoord[icoord] +
                     resdat[end - 2].cacoord[icoord]) /
                    2;
      }
    } else if (ssedat[isse].ssetype == SSE_STRAND) {
      /**  sse start point for STRAND  **/
      for (icoord = 1; icoord <= Ncoord; icoord++) {
        a[icoord] =
            (resdat[init].cacoord[icoord] + resdat[init + 1].cacoord[icoord]) /
            2;
      }
      /**  sse end point for STRAND  **/
      for (icoord = 1; icoord <= Ncoord; icoord++) {
        b[icoord] =
            (resdat[end - 1].cacoord[icoord] + resdat[end].cacoord[icoord]) / 2;
      }
      /**  sse perp point for STRAND  **/
      for (icoord = 1; icoord <= Ncoord; icoord++) {
        c[icoord] = resdat[init + 1].cacoord[icoord];
      }
    }

    /**  calc parallel vector  **/
    sum = 0.0;
    for (i = 1; i <= 3; i++) {
      ssedat[isse].para[i] = (b[i] - a[i]);
      sum += ssedat[isse].para[i] * ssedat[isse].para[i];
    }
    sum = (float)sqrt(sum);
    for (i = 1; i <= 3; i++) {
      ssedat[isse].para[i] /= sum;
    }

    /**  calc mid point  **/
    for (i = 1; i <= 3; i++) {
      ssedat[isse].mid[i] = (a[i] + b[i]) / 2;
    }

    /**  calc perpendicular vector  **/
    calc_perp(ssedat[isse].perp, ssedat[isse].para, ssedat[isse].mid, c);
  }

  return;
}

PRIVATE void calc_perp(float *n, float *v, float *a, float *b) {
  int i;
  float b_a[3 + 1];
  float pro;
  float sum;

  pro = 0.0;
  for (i = 1; i <= 3; i++) {
    b_a[i] = b[i] - a[i];
    pro += b_a[i] * v[i];
  }

  sum = 0.0;
  for (i = 1; i <= 3; i++) {
    n[i] = b_a[i] - pro * v[i];
    sum += n[i] * n[i];
  }
  sum = (float)sqrt(sum);

  for (i = 1; i <= 3; i++) {
    n[i] /= sum;
  }

  return;
}
