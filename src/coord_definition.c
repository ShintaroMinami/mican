#include "mican.h"

void coord_definition(float *x, float *y, float *z, SSEDAT *ssedat, int isse) {
  int i;

  /**  coordinate x  **/
  for (i = 1; i <= 3; i++) {
    x[i] = ssedat[isse].para[i];
  }
  /**  coordinate y  **/
  for (i = 1; i <= 3; i++) {
    y[i] = ssedat[isse].perp[i];
  }
  /**  coordinate z  **/
  z[1] = x[2] * y[3] - x[3] * y[2];
  z[2] = x[3] * y[1] - x[1] * y[3];
  z[3] = x[1] * y[2] - x[2] * y[1];

  return;
}

void coord_definition_r(float *x, float *y, float *z, SSEDAT *ssedat,
                        int isse) {
  int i;

  /**  coordinate x  **/
  for (i = 1; i <= 3; i++) {
    x[i] = -ssedat[isse].para[i];
  }
  /**  coordinate y  **/
  for (i = 1; i <= 3; i++) {
    y[i] = ssedat[isse].perp[i];
  }
  /**  coordinate z  **/
  z[1] = x[2] * y[3] - x[3] * y[2];
  z[2] = x[3] * y[1] - x[1] * y[3];
  z[3] = x[1] * y[2] - x[2] * y[1];

  return;
}
