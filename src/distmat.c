#include "mican.h"

void calc_distmat(int naa, RESDAT *resdat, float **distmat){
  int   iaa, jaa;
  float dist;

  for(iaa=1; iaa<=naa; iaa++){
    for(jaa=iaa+1; jaa<=naa; jaa++){
      dist =
	  (resdat[iaa].cacoord[1] - resdat[jaa].cacoord[1])
	* (resdat[iaa].cacoord[1] - resdat[jaa].cacoord[1])
	+ (resdat[iaa].cacoord[2] - resdat[jaa].cacoord[2])
	* (resdat[iaa].cacoord[2] - resdat[jaa].cacoord[2])
	+ (resdat[iaa].cacoord[3] - resdat[jaa].cacoord[3])
   	* (resdat[iaa].cacoord[3] - resdat[jaa].cacoord[3]);
      dist = sqrt(dist);
      distmat[iaa][jaa] = dist;
      distmat[jaa][iaa] = dist;
    }
  }

  return;
}
