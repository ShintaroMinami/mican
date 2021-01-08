#include "mican.h"

PRIVATE float determinant(double *a);






void bstft(int naa, RESDAT *resdat1, RESDAT *resdat2, float *vec, float **rot) {
  int i, j;
  int iaa;
  float **co1, **co2;
  float *xc1, *xc2;
  float det; /* determinant of U */
  float *sign;
  int smallest_w_i;

  double *a, *w, *u, *v, *vt;

  /*********************/
  /*      calloc       */
  /*********************/
  co1 = (float **)calloc((size_t)(naa + 1), sizeof(float *));
  co2 = (float **)calloc((size_t)(naa + 1), sizeof(float *));
  for (iaa = 0; iaa <= naa; iaa++) {
    co1[iaa] = (float *)calloc((size_t)(3 + 1), sizeof(float));
    co2[iaa] = (float *)calloc((size_t)(3 + 1), sizeof(float));
  }

  xc1 = (float *)calloc((size_t)(naa + 1), sizeof(float));
  xc2 = (float *)calloc((size_t)(naa + 1), sizeof(float));
  sign = (float *)calloc((size_t)(naa + 1), sizeof(float));

  a = (double *)calloc((size_t)(3 * 3), sizeof(double));
  w = (double *)calloc((size_t)(3 * 3), sizeof(double));
  u = (double *)calloc((size_t)(3 * 3), sizeof(double));
  v = (double *)calloc((size_t)(3 * 3), sizeof(double));
  vt = (double *)calloc((size_t)(3 * 3), sizeof(double));

  /****  copy coordinate  ****/
  for (iaa = 1; iaa <= naa; iaa++) {
    for (i = 1; i <= 3; i++) {
      co1[iaa][i] = resdat1[iaa].cacoord[i];
      co2[iaa][i] = resdat2[iaa].cacoord[i];
    }
  }

  /****  initialize  ****/
  for (i = 1; i <= 3; i++) {
    sign[i] = 1.0;
  }

  for (i = 1; i <= 3; i++) {
    xc1[i] = 0.0;
    xc2[i] = 0.0;
  }

  /****  calculate centroids  *****/
  for (iaa = 1; iaa <= naa; iaa++) {
    for (i = 1; i <= 3; i++) {
      xc1[i] += co1[iaa][i];
      xc2[i] += co2[iaa][i];
    }
  }

  for (i = 1; i <= 3; i++) {
    xc1[i] = xc1[i] / (float)(naa);
    xc2[i] = xc2[i] / (float)(naa);
  }

  /****  translation  ****/
  for (iaa = 1; iaa <= naa; iaa++) {
    for (i = 1; i <= 3; i++) {
      co1[iaa][i] = co1[iaa][i] - xc1[i];
      co2[iaa][i] = co2[iaa][i] - xc2[i];
    }
  }

  /******************************/
  /*     a[i][j]: U matrix      */
  /******************************/
  /****  initialize  ****/
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      a[i * 3 + j] = 0.0;
    }
  }

  /****  calcurate U matrix  ****/
  for (iaa = 1; iaa <= naa; iaa++) {
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        a[i * 3 + j] = a[i * 3 + j] + co1[iaa][i + 1] * co2[iaa][j + 1];
      }
    }
  }

  /****  calculate determinant A  ****/
  det = determinant(a);

  /**********************************/
  /*  singular value decomposition  */
  /**********************************/
  svd(3, 3, 1, 1, EPS, EPS, a, w, u, v, vt);

  /****  sign  ****/
  if (det < 0.0) {
    smallest_w_i = 0;
    if (w[0] > w[1]) {
      smallest_w_i = 1;
      if (w[1] > w[2]) {
        smallest_w_i = 2;
      }
    } else {
      if (w[0] > w[2]) {
        smallest_w_i = 2;
      }
    }
    sign[smallest_w_i] = -1.0;
  }

  /**************************************/
  /*  best-fit rotation matric r[3][3]  */
  /**************************************/
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      rot[i + 1][j + 1] = (float)(sign[1] * u[i * 3 + 0] * v[j * 3 + 0] +
                                  sign[2] * u[i * 3 + 1] * v[j * 3 + 1] +
                                  sign[3] * u[i * 3 + 2] * v[j * 3 + 2]);
    }
  }

  /****************************************/
  /*  best-fit translation vector vec[3]  */
  /****************************************/
  for (i = 1; i <= 3; i++) {
    vec[i] = xc1[i];
    for (j = 1; j <= 3; j++) {
      vec[i] += -rot[i][j] * xc2[j];
    }
  }

  /******************/
  /*      free      */
  /******************/
  for (iaa = 0; iaa <= naa; iaa++) {
    free(co1[iaa]);
    free(co2[iaa]);
  }
  free(co1);
  free(co2);
  free(xc1);
  free(xc2);
  free(sign);
  free(a);
  free(w);
  free(u);
  free(v);
  free(vt);

  return;
} /* end of function */



void wbstft(int naa, RESDAT *resdat1, RESDAT *resdat2, float *vec, float **rot, float *weight)
{
  int     i,j;
  int     iaa;
  float   **co1, **co2;
  float   *xc1, *xc2;
  float   det;     /* determinant of U */
  float   *sign;
  float   weight_sum;
  int     smallest_w_i;

  double  *a, *w, *u, *v,*vt;

  /*********************/
  /*      calloc       */
  /*********************/
  co1 = (float **)calloc((size_t)(naa+1), sizeof(float *));
  co2 = (float **)calloc((size_t)(naa+1), sizeof(float *));
  for(iaa=0; iaa<=naa; iaa++){
    co1[iaa] = (float *)calloc((size_t)(3+1), sizeof(float));
    co2[iaa] = (float *)calloc((size_t)(3+1), sizeof(float));
  }

  xc1  = (float *)calloc((size_t)(naa+1), sizeof(float));
  xc2  = (float *)calloc((size_t)(naa+1), sizeof(float));
  sign = (float *)calloc((size_t)(naa+1), sizeof(float));

  a  = (double *)calloc((size_t)(3*3), sizeof(double));
  w  = (double *)calloc((size_t)(3*3), sizeof(double));
  u  = (double *)calloc((size_t)(3*3), sizeof(double));
  v  = (double *)calloc((size_t)(3*3), sizeof(double));
  vt = (double *)calloc((size_t)(3*3), sizeof(double));


  /****  copy coordinate  ****/
  for(iaa=1; iaa<=naa; iaa++){
    for(i=1; i<=3; i++){
      co1[iaa][i] = resdat1[iaa].cacoord[i];
      co2[iaa][i] = resdat2[iaa].cacoord[i];
    }
  }
 

  /****  initialize  ****/
  for(i=1;i<=3;i++){ sign[i] = 1.0; }

  for(i=1; i<=3; i++){
    xc1[i] = 0.0;
    xc2[i] = 0.0;
  }

  /****  calculate centroids  *****/
  weight_sum = 0;
  for(iaa=1; iaa<=naa; iaa++){
    for(i=1; i<=3; i++){
      xc1[i] += weight[iaa]*co1[iaa][i]; // weighted
      xc2[i] += weight[iaa]*co2[iaa][i]; // weighted
    }
    weight_sum += weight[iaa];
  }
    
  for(i=1;i<=3;i++){
    xc1[i] = xc1[i] / weight_sum; // weighted
    xc2[i] = xc2[i] / weight_sum; // weighted
  }

  /****  translation  ****/
  for(iaa=1;iaa<=naa;iaa++){
    for(i=1;i<=3;i++){
      co1[iaa][i] = co1[iaa][i] - xc1[i];
      co2[iaa][i] = co2[iaa][i] - xc2[i];
    }
  }

  /******************************/    
  /*     a[i][j]: U matrix      */
  /******************************/
  /****  initialize  ****/
  for(i=0; i<3; i++){
    for(j=0; j<3; j++){ a[i*3+j] = 0.0; }
  }

  /****  calcurate U matrix  ****/
  for(iaa=1; iaa<=naa; iaa++){
    for(i=0;i<3;i++){
      for(j=0;j<3;j++){
	a[i*3+j] = a[i*3+j] + weight[iaa]*co1[iaa][i+1]*co2[iaa][j+1]; // weighted
      }
    }
  }
  
  /****  calculate determinant A  ****/
  det = determinant(a);


  /**********************************/
  /*  singular value decomposition  */
  /**********************************/
  svd(3,3,1,1,EPS,EPS,a,w,u,v,vt);


  /****  sign  ****/
  if( det < 0.0){
    smallest_w_i = 0;
    if(w[0] > w[1]){
      smallest_w_i = 1;    
      if(w[1] > w[2]){ smallest_w_i = 2; }
    }
    else{
      if(w[0] > w[2]){ smallest_w_i = 2; }
    }
    sign[smallest_w_i] = -1.0;
  }

  /**************************************/
  /*  best-fit rotation matric r[3][3]  */
  /**************************************/
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      rot[i+1][j+1] = (float)(  sign[1] * u[i*3 + 0]*v[j*3 + 0] 
				+ sign[2] * u[i*3 + 1]*v[j*3 + 1] 
				+ sign[3] * u[i*3 + 2]*v[j*3 + 2]);
    }
  }

  /****************************************/
  /*  best-fit translation vector vec[3]  */
  /****************************************/
  for(i=1; i<=3; i++){
    vec[i] = xc1[i];
    for(j=1; j<=3; j++){ vec[i] += -rot[i][j] * xc2[j]; }
  }
  
  /******************/
  /*      free      */
  /******************/
  for(iaa=0; iaa<=naa; iaa++){
    free(co1[iaa]);
    free(co2[iaa]);
  }
  free(co1);
  free(co2);
  free(xc1);
  free(xc2);
  free(sign);
  free(a);
  free(w);
  free(u);
  free(v);
  free(vt);

  return;
}/* end of function */




/**************************************************/
/*           calcurate  determinant               */
/**************************************************/
PRIVATE float determinant(double *a) {
  double dt;

  dt = a[0] * a[4] * a[8] + a[1] * a[5] * a[6] + a[2] * a[3] * a[7] -
       a[6] * a[4] * a[2] - a[7] * a[5] * a[0] - a[8] * a[3] * a[1];

  return ((float)dt);
} /* end of function */
