#include "mican.h"

void make_pdbdat(int naa, RESDAT *resdat, PDBDAT *pdbdat) {

  int iaa;
  int nchain;

  /* get chainID_org */
  /* get naa */
  nchain = 0;
  for (iaa = 1; iaa <= naa; iaa++) {
    if (nchain < resdat[iaa].chainID) {
      nchain = resdat[iaa].chainID;
      strcpy(pdbdat->chainID_org[nchain], resdat[iaa].chainID_org);
      pdbdat->chainID[nchain] = nchain;
      pdbdat->naa[nchain] = 0;
    }
    pdbdat->naa[nchain]++;
  }

  return;
}

int count_chain(int naa, RESDAT *resdat) {

  int iaa;
  int nchain;

  nchain = 0;
  for (iaa = 1; iaa <= naa; iaa++) {
    if (nchain < resdat[iaa].chainID) {
      nchain = resdat[iaa].chainID;
    }
  }
  return (nchain);
}
