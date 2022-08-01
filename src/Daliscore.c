#include "mican.h"
#include "main.h"

float Daliscore(int naa1, float **distmat1, float **distmat2, ALIGN *align) {

  int iaa, jaa;
  int iaa1, iaa2, jaa1, jaa2;
  int ialn, jaln, naln;
  int *ialn2iaa, *ialn2jaa;
  float Daliscore;
  float Dali_Zscore;
  float dist1, dist2, dist_ave;
  float L, mean, sigma;

  naln = align->naa_align;
  ialn2iaa = (int *)calloc((size_t)(naln+1), sizeof(int *));
  ialn2jaa = (int *)calloc((size_t)(naln+1), sizeof(int *));

  /**  aligned residues  **/
  ialn = 0;
  for(iaa=1; iaa<=naa1; iaa++){
    align->local_Dali[iaa] = 0.0;
    jaa = align->align_q[iaa];
    if(jaa < 1) continue;
    ialn ++;
    ialn2iaa[ialn] = iaa;
    ialn2jaa[ialn] = jaa;
  }

  /**  calc Daliscore  **/
  Daliscore = 0.0;
  for(ialn=1; ialn<=naln; ialn++){
    iaa1 = ialn2iaa[ialn];
    jaa1 = ialn2jaa[ialn];
    align->local_Dali[iaa1] = 0.0;
    for(jaln=1; jaln<=naln; jaln++){
      if(ialn == jaln){
	align->local_Dali[iaa1] += Dali_theteE;
	continue;
      }
      iaa2 = ialn2iaa[jaln];
      jaa2 = ialn2jaa[jaln];

      dist1 = distmat1[iaa1][iaa2];
      dist2 = distmat2[jaa1][jaa2];
      dist_ave = (dist1 + dist2)/2.0;

      align->local_Dali[iaa1] += 
	(Dali_theteE - fabs(dist1-dist2)/dist_ave)
	* exp(-dist_ave*dist_ave/Dali_alpha);
    }
    Daliscore += align->local_Dali[iaa1];
    if(align->local_Dali[iaa1] < 0) align->local_Dali[iaa1] = EPS;
  }

  /**  Dali score  **/
  align->Daliscore = Daliscore;

  /**  Dali Zscore  **/
  L = input.naa_eff;
  if(L > Dali_Lmax) L = Dali_Lmax;

  mean = Dali_c1 + Dali_c2*L + Dali_c3*L*L + Dali_c4*L*L*L;

  /** This is found at l519 in dp.f of DaliLite_3.3 software **/
  if(input.naa_eff > Dali_Lmax) mean += input.naa_eff - L;
  /************************************************************/

  sigma = 0.50 * mean;

  Dali_Zscore = (Daliscore - mean)/sigma;

  /** free **/
  free(ialn2iaa);
  free(ialn2jaa);

  return(Dali_Zscore);
}
