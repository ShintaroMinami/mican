#include "mican.h"
#include "fsscanf.h"

PRIVATE void get_xyz(char *buffer, int naa, int *iaa, int *iatom,
		     RESDAT *resdat, ALLATM *allatm, int chainID);
PRIVATE int column_chainID = 21;

/****************************/
/*    FUNCTION  get_naa     */
/****************************/
void get_naa_natom(char *filename, char *chainID, int *naa, int *natom) {
  int iaa, i;
  int iatom;
  int sw, len;
  FILE *fp;
  char buffer[BUFSIZE + 1];
  char res[4 + 1];
  char ch;

  len = (int)strlen(chainID);

  /** initialize **/
  iaa = 0;
  iatom = 0;

  /** open file **/
  if ((fp = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "!! Error   : No such file (%s)\n", filename);
    exit(1);
  }

  /** count CA line **/
  while (fgets(buffer, BUFSIZE, fp) != NULL) {
    if (strncmp(buffer, "ENDMDL", 6) == 0) {
      break;
    }
    if (strncmp(buffer, "ATOM ", 5) == 0 || strncmp(buffer, "HETATM", 6) == 0) {
      fsscanf(buffer, "%12x%4s%5x%c", res, &ch);
      if (strcmp(chainID, "0") != 0) {
        sw = OFF;
        for (i = 0; i <= len; i++) {
          if (ch == chainID[i]) {
            sw = ON;
          }
        }
        if (sw == OFF) {
          continue;
        }
      }
      iatom++;
      if (strncmp(res, " CA ", 4) == 0) {
        iaa++;
      }
    }
  }

  /** Error message (No CA atoms) **/
  if (iaa == 0) {
    if (strcmp(chainID, "0") == 0) {
      fprintf(stderr, "!! Error   : No CA atoms in file (%s)\n", filename);
    } else {
      fprintf(stderr, "!! Error   : No CA atoms of chain %s in file (%s)\n",
              chainID, filename);
    }
    exit(1);
  }

  *naa = iaa;
  *natom = iatom;

  fclose(fp);

  return;
}

/*****************************/
/*                           */
/*    FUNCTION  check_CA     */
/*                           */
/*****************************/
void check_CA(char *filename, char *chainID, RESDAT *resdat) {
  int iaa, i;
  FILE *fp;
  char buffer[BUFSIZE + 1];
  char res[4 + 1], ch;
  int org;
  int sw, len;

  len = (int)strlen(chainID);

  /** initialize **/
  iaa = 0;

  /** open file **/
  if ((fp = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "!! Error   : No such file (%s)\n", filename);
    exit(1);
  }

  /** get iaa_org **/
  while (fgets(buffer, BUFSIZE, fp) != NULL) {
    if (strncmp(buffer, "ENDMDL", 6) == 0) {
      break;
    }
    if (strncmp(buffer, "ATOM ", 5) == 0 || strncmp(buffer, "HETATM", 6) == 0) {
      fsscanf(buffer, "%12x%4s%5x%c%4d", res, &ch, &org);
      if (strcmp(chainID, "0") != 0) {
        sw = OFF;
        for (i = 0; i <= len; i++) {
          if (ch == chainID[i]) {
            sw = ON;
          }
        }
        if (sw == OFF) {
          continue;
        }
      }
      if (strncmp(res, " CA ", 4) == 0) {
        iaa++;
        resdat[iaa].iaa_org = org;
      }
    }
  }

  fclose(fp);

  return;
}

/*****************************/
/*                           */
/*    FUNCTION  read_pdb     */
/*                           */
/*****************************/
void read_pdb(char *filename, char *chainID, int naa, int natom, RESDAT *resdat,
              ALLATM *allatm) {
  int iaa, iatom;
  FILE *fp;
  char buffer[BUFSIZE + 1];
  char list_chainID[MX_CHAIN + 1];
  char ch;
  int i, ichain, nchain;
  int sw, len;

  len = (int)strlen(chainID);

  /**********************/
  /**    initialize    **/
  /**********************/
  for (iaa = 1; iaa <= naa; iaa++) {
    for (iatom = 1; iatom <= MX_ATOM_TYPE; iatom++) {
      resdat[iaa].atom[iatom].exist = OFF;
    }
  }
  iaa = 0;
  iatom = 0;

  /***********************/
  /*      open file      */
  /***********************/
  if ((fp = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "!! Error   : No such file (%s)\n", filename);
    exit(1);
  }

  /***********************/
  /*   read ATOM record  */
  /***********************/
  nchain = 0;
  ichain = 0;
  while (fgets(buffer, BUFSIZE, fp) != NULL) {
    if (strncmp(buffer, "ENDMDL", 6) == 0) {
      break;
    }
    ch = '\0';
    if (strncmp(buffer, "ATOM ", 5) == 0 || strncmp(buffer, "HETATM", 6) == 0) {
      fsscanf(buffer, "%21x%c", &ch);
      if (strcmp(chainID, "0") != 0) {
        sw = OFF;
        for (i = 0; i <= len; i++) {
          if (ch == chainID[i]) {
            sw = ON;
          }
        }
        if (sw == OFF) {
          continue;
        }
      }
      if (nchain == 0) {
        nchain = 1;
        list_chainID[nchain] = buffer[column_chainID];
      } else {
        ichain = 0;
        for (i = 1; i <= nchain; i++) {
          if (buffer[column_chainID] == list_chainID[i]) {
            ichain = i;
            break;
          }
        }
        if (ichain == 0) {
          nchain++;
          list_chainID[nchain] = buffer[column_chainID];
          ichain = nchain;
        }
      }
      get_xyz(buffer, naa, &iaa, &iatom, resdat, allatm, ichain);
    }
  }
  fclose(fp);

  return;
}

/*****************************/
/*                           */
/*        SUBROUTINE         */
/*                           */
/*****************************/
PRIVATE void get_xyz(char *buffer, int naa, int *iaa, int *iatom,
		     RESDAT *resdat, ALLATM *allatm, int chainID) {
  /*   for atom lines  */
  char record[6 + 1];
  int serial;
  char name[4 + 1];
  char altLoc[1 + 1];
  char resName[3 + 1];
  char chainID_org[1 + 1];
  int resSeq;
  char iCode[1 + 1];
  float x;
  float y;
  float z;
  float occupancy;
  float tempFactor;
  char segID[4 + 1];
  char element[2 + 1];
  char charge[2 + 1];

  /*               1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7  8  9  0 */
  fsscanf(buffer,
          "%6s%5d%1x%4s%1s%3s%1x%1s%4d%1s%3x%8f%8f%8f%6f%6f%6x%4s%2s%2s",
          record /* %6s (1) */
          ,
          &serial /* %5d (2) */
                  /* %1x (3) */
          ,
          name /* %4s (4) */
          ,
          altLoc /* %1s (5) */
          ,
          resName /* %3s (6) */
                  /* %1x (7) */
          ,
          chainID_org /* %1s (8) */
          ,
          &resSeq /* %4d (9) */
          ,
          iCode /* %8f (0) */
                /* %1x (1) */
          ,
          &x /* %8f (2) */
          ,
          &y /* %8f (3) */
          ,
          &z /* %6f (4) */
          ,
          &occupancy /* %6f (5) */
          ,
          &tempFactor /* %6x (6) */
                      /* %6x (7) */
          ,
          segID /* %4s (8) */
          ,
          element /* %2s (9) */
          ,
          charge /* %2s (0) */
          );

  /**  check chainID  **/
  if (strncmp(chainID_org, " ", 1) == 0) {
    strcpy(chainID_org, "A");
  }
  /**************************************/
  /**         read  CA  atoms          **/
  /**************************************/
  if (strncmp(name, " CA ", 4) == 0) {
    *iaa = *iaa + 1;
    /** Residue Name **/
    strcpy(resdat[*iaa].resname, resName);
    /** Residue Character **/
    resdat[*iaa].reschar = resname2reschar(resName);
    /** Chain ID **/
    strcpy(resdat[*iaa].chainID_org, chainID_org);
    resdat[*iaa].chainID = chainID;
    /** Original Residue Number **/
    resdat[*iaa].iaa_org = resSeq;
    /** CA coordinate **/
    resdat[*iaa].cacoord[1] = x;
    resdat[*iaa].cacoord[2] = y;
    resdat[*iaa].cacoord[3] = z;
    /** For main chain output **/
    resdat[*iaa].atom[atom_CA].x = x;
    resdat[*iaa].atom[atom_CA].y = y;
    resdat[*iaa].atom[atom_CA].z = z;
    strcpy(resdat[*iaa].atom[atom_CA].name, name);
    resdat[*iaa].atom[atom_CA].exist = ON;
  }
  /***************************************/
  /**        For core output mode       **/
  /***************************************/
  else if (strncmp(name, " N  ", 4) == 0) {
    if(*iaa + 1 <= naa){
      if (resdat[*iaa + 1].iaa_org == resSeq) {
	resdat[*iaa + 1].atom[atom_N_].x = x;
	resdat[*iaa + 1].atom[atom_N_].y = y;
	resdat[*iaa + 1].atom[atom_N_].z = z;
	strcpy(resdat[*iaa + 1].atom[atom_N_].name, name);
	resdat[*iaa + 1].atom[atom_N_].exist = ON;
      }
    }
  } else if (strncmp(name, " C  ", 4) == 0) {
    if (resdat[*iaa].iaa_org == resSeq) {
      resdat[*iaa].atom[atom_C_].x = x;
      resdat[*iaa].atom[atom_C_].y = y;
      resdat[*iaa].atom[atom_C_].z = z;
      strcpy(resdat[*iaa].atom[atom_C_].name, name);
      resdat[*iaa].atom[atom_C_].exist = ON;
    }
  } else if (strncmp(name, " O  ", 4) == 0) {
    if (resdat[*iaa].iaa_org == resSeq) {
      resdat[*iaa].atom[atom_O_].x = x;
      resdat[*iaa].atom[atom_O_].y = y;
      resdat[*iaa].atom[atom_O_].z = z;
      strcpy(resdat[*iaa].atom[atom_O_].name, name);
      resdat[*iaa].atom[atom_O_].exist = ON;
    }
  } else if (strncmp(name, " CB ", 4) == 0) {
    if (resdat[*iaa].iaa_org == resSeq) {
      resdat[*iaa].atom[atom_CB].x = x;
      resdat[*iaa].atom[atom_CB].y = y;
      resdat[*iaa].atom[atom_CB].z = z;
      strcpy(resdat[*iaa].atom[atom_CB].name, name);
      resdat[*iaa].atom[atom_CB].exist = ON;
    }
  }

  /**************************************************/
  /**          For all atom output mode            **/
  /**************************************************/
  *iatom = *iatom + 1;
  /** Record **/
  strcpy(allatm[*iatom].record, record);
  /** iatom original **/
  allatm[*iatom].iatom_org = serial;
  /** Atom Name **/
  strcpy(allatm[*iatom].atmname, name);
  /** Residue Name **/
  strcpy(allatm[*iatom].resname, resName);
  /** Chain ID **/
  strcpy(allatm[*iatom].chainID_org, chainID_org);
  /** Original Residue Number **/
  allatm[*iatom].iaa_org = resSeq;
  /** Coordinate **/
  allatm[*iatom].coord[1] = x;
  allatm[*iatom].coord[2] = y;
  allatm[*iatom].coord[3] = z;
  /** Occupancy **/
  allatm[*iatom].occupancy = occupancy;
  /** Temp Factor **/
  allatm[*iatom].tempfactor = tempFactor;
  return;
}
