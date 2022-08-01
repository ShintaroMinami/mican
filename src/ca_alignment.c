#include "mican.h"
#include "main.h"

int ca_alignment(int naa_q, int naa_t, RESDAT *resdat_q, RESDAT *resdat_t,
                 ALIGN align, float ***mat, float **score, float **score_r) {
  int i;
  int icoord;
  int iaa, jaa;
  float dist_2, dmax_2;
  float sse_w;
  float angle1, angle2;
  float angle_w1, angle_w2;
  float mat_value;
  float *d0;
  float **co_q, **co_t;

  dmax_2 = input.dmax * input.dmax;
  
  /****  d0  ****/
  d0 = (float *)calloc((size_t)(Nite_loc + 1), sizeof(float));
  d0[0] = input.TM_d0 * input.TM_d0; // for TM-score d0 value
  d0[1] = input.greedy_cut1;
  d0[2] = input.greedy_cut2;
  d0[3] = input.greedy_cut3;

  /******************************/
  /**        initialize        **/
  /******************************/
  for (iaa = 1; iaa <= naa_q; iaa++) {
    align.align_q[iaa] = OFF;
    align.segmode[iaa] = OFF;
    align.termini[iaa] = OFF;
    align.local_score[iaa] = 0.0;
  }
  for (iaa = 1; iaa <= naa_t; iaa++) {
    align.align_t[iaa] = 0;
  }

  /******************************/
  /**     rotate structure     **/
  /******************************/
  co_q = (float **)calloc((size_t)(naa_q + 1), sizeof(float *));
  co_t = (float **)calloc((size_t)(naa_t + 1), sizeof(float *));

  for (iaa = 1; iaa <= naa_q; iaa++) {
    co_q[iaa] = (float *)calloc((size_t)(3 + 1), sizeof(float));
    for (icoord = 1; icoord <= Ncoord; icoord++) {
      co_q[iaa][icoord] = resdat_q[iaa].cacoord[icoord];
    }
  }
  for (iaa = 1; iaa <= naa_t; iaa++) {
    co_t[iaa] = (float *)calloc((size_t)(3 + 1), sizeof(float));
    for (i = 1; i <= Ncoord; i++) {
      co_t[iaa][i] = 0;
      for (icoord = 1; icoord <= Ncoord; icoord++) {
        co_t[iaa][i] += align.rot[i][icoord] * resdat_t[iaa].cacoord[icoord];
      }
      co_t[iaa][i] += align.vec[i];
    }
  }

  /********************************************/
  /**          initialize  matrix            **/
  /********************************************/
  for (i = 0; i <= Nite_loc; i++) {
    for (iaa = 0; iaa <= naa_q + 1; iaa++) {
      for (jaa = 0; jaa <= naa_t + 1; jaa++) {
	mat[i][iaa][jaa] = -100.0;
      }
    }
  }

  /*********************************************/
  /**           initialize score              **/
  /*********************************************/
  for (iaa = 1; iaa <= naa_q; iaa++) {
    for (jaa = 1; jaa <= naa_t; jaa++) {
      score[iaa][jaa] = 0;
      score_r[iaa][jaa] = 0;
    }
  }

  /*****************************/
  /**    calc input matrix    **/
  /*****************************/
  for (iaa = 1; iaa <= naa_q; iaa++) {
    for (jaa = 1; jaa <= naa_t; jaa++) {
      /** distance **/
      dist_2 = distance_2(co_q[iaa], co_t[jaa]);
      if(dist_2 > dmax_2){ continue; }
      /** angle weight **/
      if (resdat_q[iaa].chain_break == ON || resdat_t[jaa].chain_break == ON) {
        angle_w1 = DEF_angle_w_null;
	angle_w2 = DEF_angle_w_null;
      } else {
	angle1 = (float)( angle(co_q[iaa - 1], co_q[iaa], co_q[iaa + 1],
				co_t[jaa - 1], co_t[jaa], co_t[jaa + 1]) );
	angle2 = (float)( angle(co_q[iaa - 1], co_q[iaa + 1], co_q[iaa - 1],
				co_t[jaa - 1], co_t[jaa + 1], co_t[jaa - 1]) );
	if(input.mode == FWandRV){
	  if(angle2 > PI/2) angle2 = PI - angle2;
	}
	else if(input.mode == REVERSE){
	  angle2 = PI - angle2;
	}
	angle_w1 = exp( - angle1 * angle1 / (input.angle_sigma * input.angle_sigma) );
	angle_w2 = exp( - angle2 * angle2 / (input.angle_sigma * input.angle_sigma) );
      }
      /** sse weight **/
      if (resdat_q[iaa].ssetype == resdat_t[jaa].ssetype) {
        sse_w = 1.0;
      } else {
	sse_w = input.ssetype_w / (1 + input.ssetype_w);
      }
      /** matrix **/
      mat_value = angle_w1 * angle_w2 * sse_w / (1 + dist_2 / d0[0]);
      for (i = 1; i <= Nite_loc; i++) {
        if (mat_value >= d0[i]) {
          mat[i][iaa][jaa] = mat_value;
        }
      }
    }
  }

  /*********************************/
  /**   global alignment for CA   **/
  /*********************************/
  global_alignment(naa_q, naa_t, mat, align.align_q, align.segmode, score, score_r,
		   resdat_q, resdat_t);

  /**  extend alignment termini  **/
  extend_termini(naa_q, naa_t, align, co_q, co_t);

  /*********************************/
  /**   local score calcuration   **/
  /*********************************/
  align.naa_align = 0;
  for (iaa = 1; iaa <= naa_q; iaa++) {
    jaa = align.align_q[iaa];
    if (jaa > 0) {
      align.naa_align ++;
      /** distance **/
      dist_2 = distance_2(co_q[iaa], co_t[jaa]);
      /** angle weight **/
      if (resdat_q[iaa].chain_break == ON || resdat_t[jaa].chain_break == ON){
	angle_w1 = DEF_angle_w_null;
	angle_w2 = DEF_angle_w_null;
      } else {
	angle1 = angle(co_q[iaa - 1], co_q[iaa], co_q[iaa + 1],
				co_t[jaa - 1], co_t[jaa], co_t[jaa + 1]);
	angle2 = angle(co_q[iaa - 1], co_q[iaa + 1], co_q[iaa - 1],
				co_t[jaa - 1], co_t[jaa + 1], co_t[jaa - 1]);
	if(input.mode == FWandRV){
	  if(angle2 > PI/2) angle2 = PI - angle2;
	}
	else if(input.mode == REVERSE){
	  angle2 = PI - angle2;
	}
	angle_w1 = exp( - angle1 * angle1 / (input.angle_sigma * input.angle_sigma) );
	angle_w2 = exp( - angle2 * angle2 / (input.angle_sigma * input.angle_sigma) );
	if(align.termini[iaa] == ON){
	  if(angle_w1 < DEF_angle_w_null) angle_w1 = DEF_angle_w_null;
	  if(angle_w2 < DEF_angle_w_null) angle_w2 = DEF_angle_w_null;
	}
      }
      /** sse weight **/
      if (resdat_q[iaa].ssetype == resdat_t[jaa].ssetype) {
        sse_w = 1.0;
      } else {
	sse_w = input.ssetype_w / (1 + input.ssetype_w);
      }
      /** local score **/
      align.local_dist[iaa] = sqrt( dist_2 );
      align.local_TM[iaa]   = 1.0 / (1 + dist_2 / d0[0]);
      align.local_mTM[iaa]  = sse_w / (1 + dist_2 / d0[0]);
      align.local_aTM[iaa]  = sse_w * angle_w1 * angle_w2 / (1 + dist_2 / d0[0]);
      /** score used for refinement **/
      align.local_score[iaa] = align.local_aTM[iaa];
    }
  }

  /***************************/
  /**  align_q --> align_t  **/
  /***************************/
  for (iaa = 1; iaa <= naa_q; iaa++) {
    if (align.align_q[iaa] > 0) {
      align.align_t[align.align_q[iaa]] = iaa;
    }
  }

  /******************/
  /**     free      */
  /******************/
  for (iaa = 1; iaa <= naa_q; iaa++) {
    free(co_q[iaa]);
  }
  free(co_q);
  for (iaa = 1; iaa <= naa_t; iaa++) {
    free(co_t[iaa]);
  }
  free(co_t);
  free(d0);

  if (align.naa_align > 2) {
    return (TRUE);
  } else {
    return (FALSE);
  }
}


