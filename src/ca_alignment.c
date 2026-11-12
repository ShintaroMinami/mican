#include "mican.h"
#include "main.h"

/* Grow an array to hold at least n elements (doubling); exits on failure. */
PRIVATE void *grow_buf(void *p, int *cap, int n, size_t sz) {
  if (n <= *cap) {
    return p;
  }
  *cap = (*cap > 0) ? *cap : 1024;
  while (*cap < n) {
    *cap *= 2;
  }
  p = realloc(p, (size_t)(*cap) * sz);
  if (p == NULL) {
    fprintf(stderr, "!! Error   : memory allocation failed\n");
    exit(1);
  }
  return p;
}

/* Per-residue part of angle(): out[0..2] = v1, out[3] = |v1|^2 for
   angle(p, c, n, ...) and out[4..6] = v2, out[7] = |v2|^2 for angle(p, n, p, ...),
   computed with exactly the same float operations as angle() in utils.c. */
PRIVATE void angle_vec(const float *p, const float *c, const float *n,
                       float *out) {
  int i;
  float s1 = 0.0, s2 = 0.0;
  for (i = 1; i <= 3; i++) {
    out[i - 1] = (c[i] - p[i]) + (c[i] - n[i]);
    out[i + 3] = (n[i] - p[i]) + (n[i] - p[i]);
    s1 += out[i - 1] * out[i - 1];
    s2 += out[i + 3] * out[i + 3];
  }
  out[3] = s1;
  out[7] = s2;
}

/* angle() between two precomputed vectors (see angle_vec). */
static inline float angle_pre(const float *u, const float *w) {
  int i;
  float ip = 0;
  for (i = 0; i < 3; i++) {
    ip += u[i] * w[i];
  }
  ip /= (float)sqrt(u[3] * w[3]);

  if(ip > 1.0) ip = 1.0;
  if(ip < -1.0) ip = -1.0;

  return ((float)acos(ip));
}

int ca_alignment(int naa_q, int naa_t, RESDAT *resdat_q, RESDAT *resdat_t,
                 ALIGN align) {
  int i;
  int icoord;
  int iaa, jaa;
  float dist_2, dmax_2;
  float sse_w;
  float angle1, angle2;
  float angle_w1, angle_w2;
  float mat_value;
  float d0[Nite_loc + 1];
  float **co_q, **co_t;
  int ncons;
  /* Buffers reused across calls (one call per GH candidate / refinement). */
  static float **co_q_buf, **co_t_buf, *co_q_mem, *co_t_mem;
  static int cap_q = 0, cap_t = 0;
  static float *t_xyz, *d2row; /* protein 2 coords as x[], y[], z[]; row dist^2 */
  static int cap_txyz = 0, cap_d2 = 0;
  float *tx, *ty, *tz, qx, qy, qz;
  static float *av_q, *av_t; /* angle vectors, 8 floats per residue */
  static int cap_avq = 0, cap_avt = 0;
  float *aq, *at;
  static int *cons_i, *cons_j;
  static float *cons_v;
  static int cap_ci = 0, cap_cj = 0, cap_cv = 0;

  dmax_2 = input.dmax * input.dmax;
  
  /****  d0  ****/
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
  /* rows 0 and naa+1 are zero padding (never read: termini are chain breaks) */
  if (naa_q > cap_q) {
    free(co_q_buf);
    free(co_q_mem);
    co_q_buf = (float **)malloc((size_t)(naa_q + 2) * sizeof(float *));
    co_q_mem = (float *)calloc((size_t)(naa_q + 2) * 4, sizeof(float));
    cap_q = naa_q;
  }
  if (naa_t > cap_t) {
    free(co_t_buf);
    free(co_t_mem);
    co_t_buf = (float **)malloc((size_t)(naa_t + 2) * sizeof(float *));
    co_t_mem = (float *)calloc((size_t)(naa_t + 2) * 4, sizeof(float));
    cap_t = naa_t;
  }
  co_q = co_q_buf;
  co_t = co_t_buf;
  for (iaa = 0; iaa <= naa_q + 1; iaa++) {
    co_q[iaa] = co_q_mem + 4 * iaa;
  }
  for (iaa = 0; iaa <= naa_t + 1; iaa++) {
    co_t[iaa] = co_t_mem + 4 * iaa;
  }

  for (iaa = 1; iaa <= naa_q; iaa++) {
    for (icoord = 1; icoord <= Ncoord; icoord++) {
      co_q[iaa][icoord] = resdat_q[iaa].cacoord[icoord];
    }
  }
  for (iaa = 1; iaa <= naa_t; iaa++) {
    for (i = 1; i <= Ncoord; i++) {
      co_t[iaa][i] = 0;
      for (icoord = 1; icoord <= Ncoord; icoord++) {
        co_t[iaa][i] += align.rot[i][icoord] * resdat_t[iaa].cacoord[icoord];
      }
      co_t[iaa][i] += align.vec[i];
    }
  }

  /**  per-residue angle vectors (rows 0 and naa+1 unused)  **/
  av_q = grow_buf(av_q, &cap_avq, 8 * (naa_q + 2), sizeof(float));
  av_t = grow_buf(av_t, &cap_avt, 8 * (naa_t + 2), sizeof(float));
  for (iaa = 1; iaa <= naa_q; iaa++) {
    angle_vec(co_q[iaa - 1], co_q[iaa], co_q[iaa + 1], av_q + 8 * iaa);
  }
  for (jaa = 1; jaa <= naa_t; jaa++) {
    angle_vec(co_t[jaa - 1], co_t[jaa], co_t[jaa + 1], av_t + 8 * jaa);
  }

  /*****************************/
  /**    calc input matrix    **/
  /*****************************/
  /* Only pairs that pass the loosest cut (mat[Nite_loc] > 0) are kept, as a
     list in (iaa, jaa) order; see global_alignment.c. */
  /* Distances of one row are computed first in a branch-free (vectorizable)
     loop over x[], y[], z[] arrays, with the same expression as distance_2. */
  t_xyz = grow_buf(t_xyz, &cap_txyz, 3 * (naa_t + 1), sizeof(float));
  d2row = grow_buf(d2row, &cap_d2, naa_t + 1, sizeof(float));
  tx = t_xyz;
  ty = t_xyz + (naa_t + 1);
  tz = t_xyz + 2 * (naa_t + 1);
  for (jaa = 1; jaa <= naa_t; jaa++) {
    tx[jaa] = co_t[jaa][1];
    ty[jaa] = co_t[jaa][2];
    tz[jaa] = co_t[jaa][3];
  }

  ncons = 0;
  for (iaa = 1; iaa <= naa_q; iaa++) {
    qx = co_q[iaa][1];
    qy = co_q[iaa][2];
    qz = co_q[iaa][3];
    for (jaa = 1; jaa <= naa_t; jaa++) {
      d2row[jaa] = (qx - tx[jaa]) * (qx - tx[jaa])
        + (qy - ty[jaa]) * (qy - ty[jaa])
        + (qz - tz[jaa]) * (qz - tz[jaa]);
    }
    for (jaa = 1; jaa <= naa_t; jaa++) {
      /** distance **/
      dist_2 = d2row[jaa];
      if(dist_2 > dmax_2){ continue; }
      /** angle weight **/
      if (resdat_q[iaa].chain_break == ON || resdat_t[jaa].chain_break == ON) {
        angle_w1 = DEF_angle_w_null;
	angle_w2 = DEF_angle_w_null;
      } else {
	aq = av_q + 8 * iaa;
	at = av_t + 8 * jaa;
	angle1 = angle_pre(aq, at);
	angle2 = angle_pre(aq + 4, at + 4);
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
      if (mat_value >= d0[Nite_loc] && mat_value > 0) {
        cons_i = grow_buf(cons_i, &cap_ci, ncons + 1, sizeof(int));
        cons_j = grow_buf(cons_j, &cap_cj, ncons + 1, sizeof(int));
        cons_v = grow_buf(cons_v, &cap_cv, ncons + 1, sizeof(float));
        cons_i[ncons] = iaa;
        cons_j[ncons] = jaa;
        cons_v[ncons] = mat_value;
        ncons++;
      }
    }
  }

  /*********************************/
  /**   global alignment for CA   **/
  /*********************************/
  global_alignment(naa_q, naa_t, ncons, cons_i, cons_j, cons_v, d0,
                   align.align_q, align.segmode, resdat_q, resdat_t);

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
	aq = av_q + 8 * iaa;
	at = av_t + 8 * jaa;
	angle1 = angle_pre(aq, at);
	angle2 = angle_pre(aq + 4, at + 4);
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
  if (align.naa_align > 2) {
    return (TRUE);
  } else {
    return (FALSE);
  }
}


