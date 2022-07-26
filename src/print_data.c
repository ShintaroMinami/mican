#include "mican.h"
#include "main.h"

/**  function print INPUT protein data  **/
void printinp(int naa_q, int naa_t, PDBDAT pdbdat_q, PDBDAT pdbdat_t) {
  char filename1[STRSIZE], filename2[STRSIZE];
  int naa1, naa2, nchain1, nchain2;
  char chainuse1[STRSIZE], chainuse2[STRSIZE];

  if (input.qtchange == 0) {
    strcpy(filename1, input.file_t);
    strcpy(filename2, input.file_q);
    naa1 = naa_t;
    naa2 = naa_q;
    nchain1 = pdbdat_t.nchain;
    nchain2 = pdbdat_q.nchain;
    strcpy(chainuse1, input.chain_t);
    strcpy(chainuse2, input.chain_q);
  } else {
    strcpy(filename1, input.file_q);
    strcpy(filename2, input.file_t);
    naa1 = naa_q;
    naa2 = naa_t;
    nchain1 = pdbdat_q.nchain;
    nchain2 = pdbdat_t.nchain;
    strcpy(chainuse1, input.chain_q);
    strcpy(chainuse2, input.chain_t);
  }

  /****  output input protein data  ****/
  if (strncmp(chainuse1, "0", 1) == 0) {
    printf(" Protein1 (%4d residues,%2d chain) = %s\n", naa1, nchain1,
           filename1);
  } else {
    printf(" Protein1 (%4d residues,%2d chain) = %s (chain %s)\n", naa1,
           nchain1, filename1, chainuse1);
  }
  if (strncmp(chainuse2, "0", 1) == 0) {
    printf(" Protein2 (%4d residues,%2d chain) = %s\n", naa2, nchain2,
           filename2);
  } else {
    printf(" Protein2 (%4d residues,%2d chain) = %s (chain %s)\n", naa2,
           nchain2, filename2, chainuse2);
  }
  printf("\n");

  return;
}

/**  function print superposition pdb  **/
void printsup(FILE *fp, int naa_q, int naa_t, int natom_q, int natom_t,
              RESDAT *resdat_q, RESDAT *resdat_t, ALLATM *allatm_q,
              ALLATM *allatm_t, ALIGN *align) {
  int i, iatom, isub_out;
  float co[4];

  isub_out = input.isub_out - 1;

  /**********************************/
  /**   output all atom  PDB file  **/
  /**********************************/
  /****  rotation  ****/
  if (input.qtchange == ON) {
    for (iatom = 1; iatom <= natom_q; iatom++) {
      for (i = 1; i <= 3; i++) {
        co[i] = align[isub_out].rot[i][1] * allatm_q[iatom].coord[1] +
                align[isub_out].rot[i][2] * allatm_q[iatom].coord[2] +
                align[isub_out].rot[i][3] * allatm_q[iatom].coord[3];
      }
      allatm_q[iatom].coord[1] = co[1] + align[isub_out].vec[1];
      allatm_q[iatom].coord[2] = co[2] + align[isub_out].vec[2];
      allatm_q[iatom].coord[3] = co[3] + align[isub_out].vec[3];
    }
  } else {
    for (iatom = 1; iatom <= natom_t; iatom++) {
      for (i = 1; i <= 3; i++) {
        co[i] = align[isub_out].rot[i][1] * allatm_t[iatom].coord[1] +
                align[isub_out].rot[i][2] * allatm_t[iatom].coord[2] +
                align[isub_out].rot[i][3] * allatm_t[iatom].coord[3];
      }
      allatm_t[iatom].coord[1] = co[1] + align[isub_out].vec[1];
      allatm_t[iatom].coord[2] = co[2] + align[isub_out].vec[2];
      allatm_t[iatom].coord[3] = co[3] + align[isub_out].vec[3];
    }
  }
  /****  output  ****/
  if (input.qtchange == OFF) {
    printrascript(fp, naa_t, naa_q, resdat_t, resdat_q, align[isub_out].align_t,
                  align[isub_out].align_q);
    printallatm(fp, natom_t, allatm_t, 1);
    printallatm(fp, natom_q, allatm_q, 2);
  } else {
    printrascript(fp, naa_q, naa_t, resdat_q, resdat_t, align[isub_out].align_q,
                  align[isub_out].align_t);
    printallatm(fp, natom_q, allatm_q, 1);
    printallatm(fp, natom_t, allatm_t, 2);
  }

  return;
}

/**  function print alignment  **/
void printali(FILE *fp, int naa_q, int naa_t, RESDAT *resdat_q,
                 RESDAT *resdat_t, ALIGN *align, int header) {
  int iaa, jaa;
  int isub_out;
  
  isub_out = input.isub_out - 1;
  
  if (input.qtchange == OFF) {
    fprintf(fp, "# PDB file1 = %s\n", input.file_t);
    fprintf(fp, "# PDB file2 = %s\n", input.file_q);
  } else {
    fprintf(fp, "# PDB file1 = %s\n", input.file_q);
    fprintf(fp, "# PDB file2 = %s\n", input.file_t);
  }

  fprintf(fp, "#\n");
  fprintf(fp, "# Alignment No. %02d\n", input.isub_out);
  fprintf(fp, "#\n");
  /** RMSD **/
  fprintf(fp, "# RMSD           %6.3f\n", align[isub_out].rmsd);
  /** Aligned length **/
  fprintf(fp, "# Length           %4d\n", align[isub_out].naa_align);
  /** mTM-score **/
  fprintf(fp, "# sTM-score      %6.4f\n", align[isub_out].TMscore_mod);
  /** C3-score **/
  fprintf(fp, "# C3-score       %6.4f\n", align[isub_out].C3score);
  /** SP-score **/
  fprintf(fp, "# SP-score       %6.4f\n", align[isub_out].SPscore);
  /** Dali Zscore **/
  fprintf(fp, "# Dali Zscore    %6.3f\n", align[isub_out].DaliZ);
  /** TM-score **/
  if (input.qtchange == OFF) {
    fprintf(fp, "# TM-score: P1   %6.4f\n", align[isub_out].TMscore_t);
    fprintf(fp, "# TM-score: P2   %6.4f\n", align[isub_out].TMscore_q);
  } else {
    fprintf(fp, "# TM-score: P1   %6.4f\n", align[isub_out].TMscore_q);
    fprintf(fp, "# TM-score: P2   %6.4f\n", align[isub_out].TMscore_t);
  }
  fprintf(fp, "#\n");
  fprintf(fp, "# Protein1  Protein2 Closeness_Score Distance\n#\n");
  /** Alignment **/
  if (input.qtchange == OFF) {
    for (iaa = 1; iaa <= naa_t; iaa++) {
      if(header==ON){ fprintf(fp, "ALIGN"); }
      fprintf(fp, "  %4d %1c %1s ", resdat_t[iaa].iaa_org,
	      resdat_t[iaa].reschar, resdat_t[iaa].chainID_org);
      jaa = align[isub_out].align_t[iaa];
      if (jaa == 0) {
	      fprintf(fp, "    . . .               .        .\n");
      } else {
	      fprintf(fp, " %4d %1c %1s           %5.2f    %5.2f\n", resdat_q[jaa].iaa_org,
		      resdat_q[jaa].reschar, resdat_q[jaa].chainID_org,
		      align[isub_out].local_mTM[jaa], align[isub_out].local_dist[jaa]);
      }
    }
  } else {
    for (iaa = 1; iaa <= naa_q; iaa++) {
      if(header==ON){ fprintf(fp, "ALIGN"); }
      fprintf(fp, "  %4d %1c %1s ", resdat_q[iaa].iaa_org,
	      resdat_q[iaa].reschar, resdat_q[iaa].chainID_org);
      jaa = align[isub_out].align_q[iaa];
      if (jaa == 0) {
	      fprintf(fp, "    . . .               .        .\n");
      } else {
	      fprintf(fp, " %4d %1c %1s           %5.2f    %5.2f\n", resdat_t[jaa].iaa_org,
		      resdat_t[jaa].reschar, resdat_t[jaa].chainID_org,
		      align[isub_out].local_mTM[iaa], align[isub_out].local_dist[iaa]);
      }
    }
  }

  return;
}

/**  function print translation matrix  **/
void printmat(FILE *fp, ALIGN *align, int header) {
  int i, j, isub_out;

  isub_out = input.isub_out - 1;

  /**  small value  **/
  for (i = 1; i <= 3; i++) {
    for (j = 1; j <= 3; j++) {
      if (fabs(align[isub_out].rot[i][j]) < EPS) {
	      align[isub_out].rot[i][j] = 0.0;
      }
    }
    if (fabs(align[isub_out].vec[i]) < EPS) {
      align[isub_out].vec[i] = 0.0;
    }
  }

  /**  output  **/
  if (header == ON){
    if (input.qtchange == ON) {
      fprintf(fp, " Alignment No. %02d  ( TM-score=%5.3f, SP-score=%5.3f )\n", input.isub_out,
	      align[isub_out].TMscore_q, align[isub_out].SPscore);
    } else {
      fprintf(fp, " Alignment No. %02d  ( TM-score=%5.3f, SP-score=%5.3f )\n", input.isub_out,
	      align[isub_out].TMscore_t, align[isub_out].SPscore);
    }
  }
  fprintf(fp, " ----- Rotation matrix to rotate Chain1 to Chain2 -----\n");
  fprintf(fp, " i       t(i)           u(i,1)      u(i,2)      u(i,3)\n");
  fprintf(fp, " 1   %11.6f      %9.6f   %9.6f   %9.6f\n", align[isub_out].vec[1],
	  align[isub_out].rot[1][1], align[isub_out].rot[1][2], align[isub_out].rot[1][3]);
  fprintf(fp, " 2   %11.6f      %9.6f   %9.6f   %9.6f\n", align[isub_out].vec[2],
	  align[isub_out].rot[2][1], align[isub_out].rot[2][2], align[isub_out].rot[2][3]);
  fprintf(fp, " 3   %11.6f      %9.6f   %9.6f   %9.6f\n", align[isub_out].vec[3],
	  align[isub_out].rot[3][1], align[isub_out].rot[3][2], align[isub_out].rot[3][3]);
  fprintf(fp, " ------------------------------------------------------\n");
  fprintf(fp, " t(i):Translation vector,  u(i,j):Rotation matrix\n");

  return;
}

void inverse_mat(float **rot, float *vec) {
  int i, j;
  float rot_tmp[4][4];
  float vec_tmp[4];

  for (i = 1; i <= 3; i++) {
    vec_tmp[i] = vec[i];
    for (j = 1; j <= 3; j++) {
      rot_tmp[i][j] = rot[j][i];
    }
  }

  for (i = 1; i <= 3; i++) {
    vec[i] = 0;
    for (j = 1; j <= 3; j++) {
      rot[i][j] = rot_tmp[i][j];
      vec[i] -= rot_tmp[i][j] * vec_tmp[j];
    }
  }

  return;
}
