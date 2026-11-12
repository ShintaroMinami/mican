#include "mican.h"

/**  function : zeroed (n+1) x (3+1) coordinate matrix in one block  **/
/* Usage: float **co = calloc_coords(n); co[i][1..3] ...; free_coords(co); */
float **calloc_coords(int n) {
  int i;
  float **co = (float **)malloc((size_t)(n + 1) * sizeof(float *));
  float *mem = (float *)calloc((size_t)(n + 1) * 4, sizeof(float));
  if (co == NULL || mem == NULL) {
    fprintf(stderr, "!! Error   : memory allocation failed\n");
    exit(1);
  }
  for (i = 0; i <= n; i++) {
    co[i] = mem + 4 * i;
  }
  return (co);
}

void free_coords(float **co) {
  if (co != NULL) {
    free(co[0]);
    free(co);
  }
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
