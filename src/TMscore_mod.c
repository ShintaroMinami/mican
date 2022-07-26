#include "mican.h"
#include "main.h"

float TMscore_mod(int naa0, int naa1, int naa2, RESDAT *resdat1, RESDAT *resdat2, ALIGN align,
                  int angle_switch, int local_switch) {
  int iaa, jaa;
  int ico;
  int i;
  float **co1;
  float **co2;
  float TMscore;
  float d, d0;
  float weight;
  float angle1, angle2;
  float angle_w1, angle_w2;

  /*******************/
  /**     calloc    **/
  /*******************/
  co1 = (float **)calloc((size_t)(naa1 + 1), sizeof(float *));
  co2 = (float **)calloc((size_t)(naa2 + 1), sizeof(float *));
  for(iaa = 1; iaa <= naa1; iaa++) {
    co1[iaa] = (float *)calloc((size_t)(3 + 1), sizeof(float));    
  }
  for(iaa = 1; iaa <= naa2; iaa++) {
    co2[iaa] = (float *)calloc((size_t)(3 + 1), sizeof(float));    
  }

  /*************************/
  /**   copy coordinate   **/
  /*************************/
  /* for CA */
  for(iaa = 1; iaa <= naa1; iaa++) {
    for (ico = 1; ico <= Ncoord; ico++) {
      co1[iaa][ico] = resdat1[iaa].cacoord[ico];
    }
  }
  for(iaa = 1; iaa <= naa2; iaa++) {
    for (ico = 1; ico <= Ncoord; ico++) {
      co2[iaa][ico] = 0;
      for (i = 1; i <= Ncoord; i++) {
        co2[iaa][ico] += align.rot[ico][i] * resdat2[iaa].cacoord[i];
      }
      co2[iaa][ico] += align.vec[ico];
    }
  }

  /**************************/
  /**       TM score       **/
  /**************************/
  if (local_switch == ON) {
    if ( align.naa_align <= TM_len_MIN) {
      d0 = TM_d0_MIN;
    } else {
      d0 = (float)(1.24 * pow(align.naa_align - TM_len_MIN, 1.0 / 3.0) - 1.8);
    }
    if (d0 < TM_d0_MIN) {
      d0 = TM_d0_MIN;
    }
    if (d0 > TM_d0_MAX) {
      d0 = TM_d0_MAX;
    }
    d0 = d0 * d0;
  } else {
    d0 = input.TM_d0;
    d0 = d0 * d0;
  }
  /****  calc  ****/
  TMscore = 0;

  for (iaa = 1; iaa <= naa1; iaa++) {
    jaa = align.align_q[iaa];
    if(jaa < 1) continue;
    /** SSE weight **/
    if(resdat1[iaa].ssetype == resdat2[jaa].ssetype) {
      weight = 1.0;
    }
    else {
      weight = input.ssetype_w / (1 + input.ssetype_w);
    }
    /** local angle weight **/
    if (angle_switch == ON) { // if angle switch ON
      if (resdat1[iaa].chain_break == ON || resdat2[jaa].chain_break == ON){
        angle_w1 = DEF_angle_w_null;
        angle_w2 = DEF_angle_w_null;
      } else {
	angle1 = (float)( angle(co1[iaa - 1], co1[iaa], co1[iaa + 1],
				co2[jaa - 1], co2[jaa], co2[jaa + 1]) );
	angle2 = (float)( angle(co1[iaa - 1], co1[iaa + 1], co1[iaa - 1],
				co2[jaa - 1], co2[jaa + 1], co2[jaa - 1]) );
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
    } else { // angle switch OFF
      angle_w1 = 1.0;
      angle_w2 = 1.0;
    }
    /** distance **/
    d = 0;
    for (ico = 1; ico <= Ncoord; ico++) {
      d += (co1[iaa][ico] - co2[jaa][ico]) * (co1[iaa][ico] - co2[jaa][ico]);
    }
    /** TMscore_mod **/
    TMscore += weight * angle_w1 * angle_w2 / (1 + d / d0);
  }

  /** normarize **/
  TMscore /= (float)naa0;

  /**************/
  /**   free   **/
  /**************/
  for (iaa = 1; iaa <= naa1; iaa++) {
    free(co1[iaa]);
  }
  for(iaa = 1; iaa <= naa2; iaa++) {
    free(co2[iaa]);
  }
  free(co1);
  free(co2);

  return (TMscore);
}
