#include "mican.h"

/**  function : distance^2  **/
float distance_2(float *a, float *b) {

  float sum = 0.0;

  sum = (a[1] - b[1]) * (a[1] - b[1])
    + (a[2] - b[2]) * (a[2] - b[2])
    + (a[3] - b[3]) * (a[3] - b[3]);

  return (sum);
}

/** function : angle **/
float angle(float *a1, float *a2, float *a3, float *b1, float *b2, float *b3) {
  float v1[3+1];
  float v2[3+1];
  float size_v1;
  float size_v2;
  float ip;
  float ang;
  int i;

  size_v1 = 0.0;
  size_v2 = 0.0;
  ip = 0;
  for (i = 1; i <= 3; i++) {
    v1[i] = (a2[i] - a1[i]) + (a2[i] - a3[i]);
    v2[i] = (b2[i] - b1[i]) + (b2[i] - b3[i]);
    size_v1 += v1[i] * v1[i];
    size_v2 += v2[i] * v2[i];
    ip += v1[i] * v2[i];
  }
  ip /= (float)sqrt(size_v1 * size_v2);

  if(ip > 1.0) ip = 1.0;
  if(ip < -1.0) ip = -1.0;

  ang = (float)acos(ip);

  return (ang);
}
