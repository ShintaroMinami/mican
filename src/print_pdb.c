#include "mican.h"
#define UNDEF -1000

void printrascript(FILE *fp, int naa1, int naa2, RESDAT *resdat1,
                   RESDAT *resdat2, int *align1, int *align2) {
  int iaa, iaa_ini, iaa_end, iaa_old;
  char cID_old;

  /******************************/
  /****     print FORMAT     ****/
  /******************************/
  /****  basic format  ****/
  fprintf(fp, "load inline\n");
  fprintf(fp, "background white\n");
  fprintf(fp, "set stereo off\n");
  fprintf(fp, "set ambient 60\n\n");

  fprintf(fp, "select all\n");
  fprintf(fp, "color bonds none\n");
  fprintf(fp, "color backbone none\n");
  fprintf(fp, "color hbonds none\n");
  fprintf(fp, "color ssbonds none\n");
  fprintf(fp, "color ribbons none\n");
  fprintf(fp, "spacefill off\n");
  fprintf(fp, "wireframe off\n");
  fprintf(fp, "ribbons off\n");
  fprintf(fp, "backbone 10\n\n");

  fprintf(fp, "select */1\n");
  fprintf(fp, "backbone 10\n");
  fprintf(fp, "color red\n");
  fprintf(fp, "select */2\n");
  fprintf(fp, "backbone 10\n");
  fprintf(fp, "color cyan\n\n");

  /****  aligned region  ****/
  /**  model 1  **/
  iaa_old = UNDEF;
  iaa_ini = UNDEF;
  iaa_end = UNDEF;
  cID_old = '\0';
  for (iaa = 1; iaa <= naa1; iaa++) {
    if (align1[iaa] != 0) {
      if (iaa_ini == UNDEF) {
        iaa_ini = resdat1[iaa].iaa_org;
        iaa_end = resdat1[iaa].iaa_org;
      } else if (resdat1[iaa].iaa_org != iaa_old + 1) {
        fprintf(fp, "select %d-%d:%c & */1\n", iaa_ini, iaa_end, cID_old);
        fprintf(fp, "backbone 80\n");
        iaa_ini = resdat1[iaa].iaa_org;
        iaa_end = resdat1[iaa].iaa_org;
      } else if (resdat1[iaa].chainID_org[0] != cID_old) {
        fprintf(fp, "select %d-%d:%c & */1\n", iaa_ini, iaa_end, cID_old);
        fprintf(fp, "backbone 80\n");
        iaa_ini = resdat1[iaa].iaa_org;
        iaa_end = resdat1[iaa].iaa_org;
      } else {
        iaa_end = resdat1[iaa].iaa_org;
      }
      iaa_old = resdat1[iaa].iaa_org;
      cID_old = resdat1[iaa].chainID_org[0];
    } else {
      if (iaa_ini != UNDEF) {
        fprintf(fp, "select %d-%d:%c & */1\n", iaa_ini, iaa_end, cID_old);
        fprintf(fp, "backbone 80\n");
        iaa_ini = UNDEF;
        iaa_end = UNDEF;
      }
    }
  }
  if (iaa_ini != UNDEF) {
    fprintf(fp, "select %d-%d:%c & */1\n", iaa_ini, iaa_end, cID_old);
    fprintf(fp, "backbone 80\n");
  }

  /**  model 2  **/
  iaa_old = UNDEF;
  iaa_ini = UNDEF;
  cID_old = '\0';
  for (iaa = 1; iaa <= naa2; iaa++) {
    if (align2[iaa] != 0) {
      if (iaa_ini == UNDEF) {
        iaa_ini = resdat2[iaa].iaa_org;
        iaa_end = resdat2[iaa].iaa_org;
      } else if (resdat2[iaa].iaa_org != iaa_old + 1) {
        fprintf(fp, "select %d-%d:%c & */2\n", iaa_ini, iaa_end, cID_old);
        fprintf(fp, "backbone 80\n");
        iaa_ini = resdat2[iaa].iaa_org;
        iaa_end = resdat2[iaa].iaa_org;
      } else if (resdat2[iaa].chainID_org[0] != cID_old) {
        fprintf(fp, "select %d-%d:%c & */2\n", iaa_ini, iaa_end, cID_old);
        fprintf(fp, "backbone 80\n");
        iaa_ini = resdat2[iaa].iaa_org;
        iaa_end = resdat2[iaa].iaa_org;
      } else {
        iaa_end = resdat2[iaa].iaa_org;
      }
      iaa_old = resdat2[iaa].iaa_org;
      cID_old = resdat2[iaa].chainID_org[0];
    } else {
      if (iaa_ini != UNDEF) {
        fprintf(fp, "select %d-%d:%c & */2\n", iaa_ini, iaa_end, cID_old);
        fprintf(fp, "backbone 80\n");
        iaa_ini = UNDEF;
        iaa_end = UNDEF;
      }
    }
  }
  if (iaa_ini != UNDEF) {
    fprintf(fp, "select %d-%d:%c & */2\n", iaa_ini, iaa_end, cID_old);
    fprintf(fp, "backbone 80\n");
  }

  fprintf(fp, "exit\n\n");

  return;
}

void printallatm(FILE *fp, int natom, ALLATM *allatm, int imodel) {
  int iatom;

  fprintf(fp, "MODEL        %1d\n", imodel);
  iatom = 0;
  for (iatom = 1; iatom <= natom; iatom++) {
    fprintf(
        fp, "%6s%5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f  %4.2f%6.2f\n", allatm[iatom].record,
        allatm[iatom].iatom_org, allatm[iatom].atmname, allatm[iatom].resname,
        allatm[iatom].chainID_org, allatm[iatom].iaa_org,
        allatm[iatom].coord[1], allatm[iatom].coord[2], allatm[iatom].coord[3],
        allatm[iatom].occupancy, allatm[iatom].tempfactor);
  }
  fprintf(fp, "ENDMDL\n");

  return;
}
