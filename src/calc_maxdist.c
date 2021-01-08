#include "mican.h"

float calc_maxdist(int nsse, SSEDAT *ssedat) {
  int icoord;
  int isse, jsse;
  float dist;
  float max_dist = 0.0f;

  for (isse = 1; isse <= nsse - 1; isse++) {
    for (jsse = isse + 1; jsse <= nsse; jsse++) {
      dist = 0;
      for (icoord = 1; icoord <= Ncoord; icoord++) {
        dist += (ssedat[isse].mid[icoord] - ssedat[jsse].mid[icoord]) *
                (ssedat[isse].mid[icoord] - ssedat[jsse].mid[icoord]);
      }
      if (dist > max_dist) {
        max_dist = dist;
      }
    }
  }

  return ((float)sqrt(max_dist));
}
