#include "mican.h"
#include "main.h"

PRIVATE void search(int *ialign, int nsse_q, int nsse_t, int nhash,
                    HASH_TABLE ***hash, SSEDAT *ssedat, HASH_RESULT *hres,
                    int mode);
PRIVATE void search_list(int id_x, int id_y, int id_z, HASH_TABLE ***hash,
                         SSEDAT *ssedat, int isse, int ksse, float *para,
                         float *perp, float *coor, float *count);

/*******************************************/
/**     FUNCTION  search_hash_table       **/
/*******************************************/
void search_hash(int nsse_q, int nsse_t, int nalign, int nhash,
                 HASH_TABLE ***hash, SSEDAT *ssedat, HASH_RESULT *hres) {

  int ialign;

  /*******************************/
  /**        initialize         **/
  /*******************************/
  for (ialign = 0; ialign <= nalign; ialign++) {
    hres[ialign].count = 0;
  }

  /********************************/
  /**      search hash table     **/
  /********************************/

  ialign = -1;
  /****  reverse only  ****/
  if (input.mode == REVERSE) {
    search(&ialign, nsse_q, nsse_t, nhash, hash, ssedat, hres, REVERSE);
  }
  /****  forward only  ****/
  else if (input.mode == FORWARD) {
    search(&ialign, nsse_q, nsse_t, nhash, hash, ssedat, hres, FORWARD);
  } else if (input.mode == SEQUENT) {
    search(&ialign, nsse_q, nsse_t, nhash, hash, ssedat, hres, FORWARD);
  }
  /****  forward/reverse mixed  ****/
  else {
    search(&ialign, nsse_q, nsse_t, nhash, hash, ssedat, hres, FORWARD);
    search(&ialign, nsse_q, nsse_t, nhash, hash, ssedat, hres, REVERSE);
  }

  return;
}

PRIVATE void search(int *ialign, int nsse_q, int nsse_t, int nhash,
                    HASH_TABLE ***hash, SSEDAT *ssedat, HASH_RESULT *hres,
                    int mode) {
  int i;
  int id_x, id_y, id_z;
  int isse, jsse;
  int isse_q;
  float *count;
  float c_a[3 + 1];
  float x[3 + 1];
  float y[3 + 1];
  float z[3 + 1];
  float para[3];
  float perp[3];
  float coor[3];

  count = (float *)calloc((size_t)(nsse_q + 1), sizeof(float));

  /********************************/
  /**      search hash table     **/
  /********************************/
  for (isse = 1; isse <= nsse_t; isse++) {

    /**  initialize  **/
    for (isse_q = 1; isse_q <= nsse_q; isse_q++) {
      count[isse_q] = 0.0;
    }

    /**  define xyz  **/
    if (mode == FORWARD) {
      coord_definition(x, y, z, ssedat, isse);
    } else {
      coord_definition_r(x, y, z, ssedat, isse);
    }

    /**  search start  **/
    for (jsse = 1; jsse <= nsse_t; jsse++) {
      if (isse == jsse) {
        continue;
      }
      for (i = 1; i <= 3; i++) {
        c_a[i] = ssedat[jsse].mid[i] - ssedat[isse].mid[i];
      }

      /** coordinate **/
      coor[0] = c_a[1] * x[1] + c_a[2] * x[2] + c_a[3] * x[3];
      coor[1] = c_a[1] * y[1] + c_a[2] * y[2] + c_a[3] * y[3];
      coor[2] = c_a[1] * z[1] + c_a[2] * z[2] + c_a[3] * z[3];
      id_x = (int)(coor[0] / input.d_hash + (float)nhash / 2);
      id_y = (int)(coor[1] / input.d_hash + (float)nhash / 2);
      id_z = (int)(coor[2] / input.d_hash + (float)nhash / 2);
      if (id_x > nhash) {
        continue;
      }
      if (id_y > nhash) {
        continue;
      }
      if (id_z > nhash) {
        continue;
      }
      if (id_x < 1) {
        continue;
      }
      if (id_y < 1) {
        continue;
      }
      if (id_z < 1) {
        continue;
      }
      if (hash[id_x][id_y][id_z].nlist < 1) {
        continue;
      }
      /**  sse parallel vector  **/
      para[0] = ssedat[jsse].para[1] * x[1] + ssedat[jsse].para[2] * x[2] +
                ssedat[jsse].para[3] * x[3];
      para[1] = ssedat[jsse].para[1] * y[1] + ssedat[jsse].para[2] * y[2] +
                ssedat[jsse].para[3] * y[3];
      para[2] = ssedat[jsse].para[1] * z[1] + ssedat[jsse].para[2] * z[2] +
                ssedat[jsse].para[3] * z[3];
      /**  sse perpendicular vector  **/
      perp[0] = ssedat[jsse].perp[1] * x[1] + ssedat[jsse].perp[2] * x[2] +
                ssedat[jsse].perp[3] * x[3];
      perp[1] = ssedat[jsse].perp[1] * y[1] + ssedat[jsse].perp[2] * y[2] +
                ssedat[jsse].perp[3] * y[3];
      perp[2] = ssedat[jsse].perp[1] * z[1] + ssedat[jsse].perp[2] * z[2] +
                ssedat[jsse].perp[3] * z[3];

      /**  search  **/
      search_list(id_x, id_y, id_z, hash, ssedat, isse, jsse, para, perp, coor,
                  count);
    }

    for (isse_q = 1; isse_q <= nsse_q; isse_q++) {
      if (count[isse_q] < input.min_GHcount) {
        continue;
      }
      *ialign += 1;
      hres[*ialign].count = count[isse_q];
      hres[*ialign].isse_q = isse_q;
      hres[*ialign].isse_t = isse;
      hres[*ialign].coord_mode = mode;
    }

    if (input.progress == ON) {
      print_progress("GE Matrix Search", isse, nsse_t);
    }
  }

  /***************/
  /**   free    **/
  /***************/
  free(count);

  return;
}

PRIVATE void search_list(int id_x, int id_y, int id_z, HASH_TABLE ***hash,
                         SSEDAT *ssedat, int isse, int jsse, float *para,
                         float *perp, float *coor, float *count) {
  int i;
  int ilist;
  int isse_q;
  float sum_para, sum_perp;

  for (ilist = 0; ilist <= hash[id_x][id_y][id_z].nlist - 1; ilist++) {

    /**  check ssetype  **/
    if (hash[id_x][id_y][id_z].ssetype[ilist] != ssedat[jsse].ssetype) {
      continue;
    }

    /**  check sse parallel vector  **/
    sum_para = 0.0;
    for (i = 0; i <= 2; i++) {
      sum_para += hash[id_x][id_y][id_z].para[ilist][i] * para[i];
    }
    if (input.mode == REVERSE) {
      sum_para *= -1;
    } else if (input.mode == FWandRV) {
      sum_para = (float)fabs((double)(sum_para));
    }
    if (sum_para < input.cutoff_para) {
      continue;
    }

    /**  check sse perpendicular vector  **/
    sum_perp = 0.0;
    for (i = 0; i <= 2; i++) {
      sum_perp += hash[id_x][id_y][id_z].perp[ilist][i] * perp[i];
    }
    if (sum_perp < input.cutoff_perp) {
      continue;
    }

    /**  count  **/
    isse_q = hash[id_x][id_y][id_z].isse[ilist];

    count[isse_q] += (sum_para - input.cutoff_para) / (1 - input.cutoff_para) +
                     (sum_perp - input.cutoff_perp) / (1 - input.cutoff_perp);
  }

  return;
}
