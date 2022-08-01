#include "mican.h"
#include "main.h"

/*************************************/
/**     FUNCTION  pre_make_hash     **/
/*************************************/
void pre_make_hash(int nsse, int nhash, HASH_TABLE ***hash, SSEDAT *ssedat) {
  int i, j, k;
  int id_x, id_y, id_z;
  int isse, ksse;
  int dx, dy, dz;
  float c_a[3 + 1];
  float x[3 + 1];
  float y[3 + 1];
  float z[3 + 1];

  /********************************/
  /**         initialize         **/
  /********************************/
  for (i = 0; i <= nhash; i++) {
    for (j = 0; j <= nhash; j++) {
      for (k = 0; k <= nhash; k++) {
        hash[i][j][k].nlist = 0;
      }
    }
  }

  /********************************/
  /**      check mem used        **/
  /********************************/
  for (isse = 1; isse <= nsse; isse++) {

    coord_definition(x, y, z, ssedat, isse);

    for (ksse = 1; ksse <= nsse; ksse++) {
      if (ksse == isse) {
        continue;
      }
      for (i = 1; i <= 3; i++) {
        c_a[i] = ssedat[ksse].mid[i] - ssedat[isse].mid[i];
      }
      id_x =
          (int)((c_a[1] * x[1] + c_a[2] * x[2] + c_a[3] * x[3]) / input.d_hash +
                (float)nhash / 2);
      id_y =
          (int)((c_a[1] * y[1] + c_a[2] * y[2] + c_a[3] * y[3]) / input.d_hash +
                (float)nhash / 2);
      id_z =
          (int)((c_a[1] * z[1] + c_a[2] * z[2] + c_a[3] * z[3]) / input.d_hash +
                (float)nhash / 2);

      for (dx = -input.max_nghbr; dx <= input.max_nghbr; dx++) {
        if (id_x + dx < 1 || nhash < id_x + dx) {
          continue;
        }
        for (dy = -input.max_nghbr; dy <= input.max_nghbr; dy++) {
          if (id_y + dy < 1 || nhash < id_y + dy) {
            continue;
          }
          for (dz = -input.max_nghbr; dz <= input.max_nghbr; dz++) {
            if (id_z + dz < 1 || nhash < id_z + dz) {
              continue;
            }

            hash[id_x + dx][id_y + dy][id_z + dz].nlist++;
          }
        }
      }
    }
  }

  return;
}

/*****************************************/
/**     FUNCTION  make_hash_table       **/
/*****************************************/
void make_hash_table(int nsse, int nhash, HASH_TABLE ***hash, SSEDAT *ssedat) {
  int i;
  int id_x, id_y, id_z;
  int isse, ksse;
  int ili;
  int dx, dy, dz;
  float c_a[3 + 1];
  float x[3 + 1];
  float y[3 + 1];
  float z[3 + 1];
  float para[3 + 1];
  float perp[3 + 1];
  float coor[3 + 1];

  /*******************************/
  /**        initialize         **/
  /*******************************/
  for (id_x = 0; id_x <= nhash; id_x++) {
    for (id_y = 0; id_y <= nhash; id_y++) {
      for (id_z = 0; id_z <= nhash; id_z++) {
        for (ili = 0; ili <= hash[id_x][id_y][id_z].nlist - 1; ili++) {
          hash[id_x][id_y][id_z].isse[ili] = 0;
          hash[id_x][id_y][id_z].ssetype[ili] = SSE_COIL;
          for (i = 0; i <= 2; i++) {
            hash[id_x][id_y][id_z].para[ili][i] = 0;
            hash[id_x][id_y][id_z].perp[ili][i] = 0;
          }
        }
      }
    }
  }
  for (id_x = 1; id_x <= nhash; id_x++) {
    for (id_y = 1; id_y <= nhash; id_y++) {
      for (id_z = 1; id_z <= nhash; id_z++) {
        hash[id_x][id_y][id_z].nlist = 0;
      }
    }
  }

  /********************************/
  /**       make hash table      **/
  /********************************/
  for (isse = 1; isse <= nsse; isse++) {

    coord_definition(x, y, z, ssedat, isse);

    for (ksse = 1; ksse <= nsse; ksse++) {
      if (ksse == isse) {
        continue;
      }
      for (i = 1; i <= 3; i++) {
        c_a[i] = ssedat[ksse].mid[i] - ssedat[isse].mid[i];
      }

      /**  coordinate  **/
      coor[1] = c_a[1] * x[1] + c_a[2] * x[2] + c_a[3] * x[3];
      coor[2] = c_a[1] * y[1] + c_a[2] * y[2] + c_a[3] * y[3];
      coor[3] = c_a[1] * z[1] + c_a[2] * z[2] + c_a[3] * z[3];
      id_x = (int)(coor[1] / input.d_hash + (float)nhash / 2);
      id_y = (int)(coor[2] / input.d_hash + (float)nhash / 2);
      id_z = (int)(coor[3] / input.d_hash + (float)nhash / 2);

      /**  sse vecter  **/
      para[1] = ssedat[ksse].para[1] * x[1] + ssedat[ksse].para[2] * x[2] +
                ssedat[ksse].para[3] * x[3];
      para[2] = ssedat[ksse].para[1] * y[1] + ssedat[ksse].para[2] * y[2] +
                ssedat[ksse].para[3] * y[3];
      para[3] = ssedat[ksse].para[1] * z[1] + ssedat[ksse].para[2] * z[2] +
                ssedat[ksse].para[3] * z[3];
      /**  sse perpendicular vecter  **/
      perp[1] = ssedat[ksse].perp[1] * x[1] + ssedat[ksse].perp[2] * x[2] +
                ssedat[ksse].perp[3] * x[3];
      perp[2] = ssedat[ksse].perp[1] * y[1] + ssedat[ksse].perp[2] * y[2] +
                ssedat[ksse].perp[3] * y[3];
      perp[3] = ssedat[ksse].perp[1] * z[1] + ssedat[ksse].perp[2] * z[2] +
                ssedat[ksse].perp[3] * z[3];

      for (dx = -input.max_nghbr; dx <= input.max_nghbr; dx++) {
        if (id_x + dx < 1 || nhash < id_x + dx) {
          continue;
        }
        for (dy = -input.max_nghbr; dy <= input.max_nghbr; dy++) {
          if (id_y + dy < 1 || nhash < id_y + dy) {
            continue;
          }
          for (dz = -input.max_nghbr; dz <= input.max_nghbr; dz++) {
            if (id_z + dz < 1 || nhash < id_z + dz) {
              continue;
            }
            hash[id_x + dx][id_y + dy][id_z + dz].nlist++;
            hash[id_x + dx][id_y + dy][id_z + dz]
                .isse[hash[id_x + dx][id_y + dy][id_z + dz].nlist - 1] = isse;
            hash[id_x + dx][id_y + dy][id_z + dz]
                .ssetype[hash[id_x + dx][id_y + dy][id_z + dz].nlist - 1] =
                ssedat[ksse].ssetype;
            for (i = 1; i <= 3; i++) {
              hash[id_x + dx][id_y + dy][id_z + dz]
                  .para[hash[id_x + dx][id_y + dy][id_z + dz].nlist -
                        1][i - 1] = para[i];
              hash[id_x + dx][id_y + dy][id_z + dz]
                  .perp[hash[id_x + dx][id_y + dy][id_z + dz].nlist -
                        1][i - 1] = perp[i];
            }
          }
        }
      }
    }
  }

  return;
}
