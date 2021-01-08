#include "mican.h"

void chain_break_check(int naa, RESDAT *resdat) {

  int iaa;

  for (iaa = 2; iaa <= naa - 1; iaa++) {
    /* chain break check */
    if (distance_2(resdat[iaa - 1].cacoord, resdat[iaa].cacoord) > 16.0) {
      resdat[iaa].chain_break = ON;
    } else if (distance_2(resdat[iaa + 1].cacoord, resdat[iaa].cacoord) >
               16.0) {
      resdat[iaa].chain_break = ON;
    } else {
      resdat[iaa].chain_break = OFF;
    }
  }
  resdat[1].chain_break = ON;
  resdat[naa].chain_break = ON;

  return;
}
