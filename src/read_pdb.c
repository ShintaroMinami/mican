#include "mican.h"
#include "fsscanf.h"

PRIVATE void get_xyz(char *buffer, int naa, int *iaa, int *iatom,
		     RESDAT *resdat, ALLATM *allatm, int chainID);
PRIVATE int column_chainID = 21;

/*
 * Fixed-column readers for PDB lines, equivalent to the fsscanf() formats
 * they replace: a line ends at its first '\n' or NUL, and a field is the part
 * of columns [start, start + width) before that end.
 *   field_str : like "%<width>s" (padded with spaces to width)
 *   field_chr : like "%c"        (' ' when the line is shorter)
 *   field_int / field_flt : like "%<width>d" / "%<width>f"
 *                (0 for a blank field, unchanged if not a number)
 * Example: line "ATOM      1  CA  ALA A   1 ..." -> field_str(.., 12, 4) = " CA ".
 */
PRIVATE int line_end(const char *s) { return (int)strcspn(s, "\n"); }

PRIVATE int get_field(const char *s, int end, int start, int width,
                      char *buf) {
  int n = 0;
  if (start < end) {
    n = (start + width <= end) ? width : end - start;
    memcpy(buf, s + start, (size_t)n);
  }
  buf[n] = '\0';
  return n;
}

PRIVATE void field_str(const char *s, int end, int start, int width,
                       char *out) {
  int n = get_field(s, end, start, width, out);
  for (; n < width; n++) {
    out[n] = ' ';
  }
  out[width] = '\0';
}

PRIVATE char field_chr(const char *s, int end, int start) {
  return (start < end && s[start] != '\0') ? s[start] : ' ';
}

PRIVATE void field_int(const char *s, int end, int start, int width,
                       int *out) {
  char buf[16];
  get_field(s, end, start, width, buf);
  if (scan_int(buf, out) == (-1)) {
    *out = 0;
  }
}

PRIVATE void field_flt(const char *s, int end, int start, int width,
                       float *out) {
  char buf[16];
  get_field(s, end, start, width, buf);
  if (scan_float(buf, out) == (-1)) {
    *out = (float)0.0;
  }
}

/****************************/
/*    FUNCTION  get_naa     */
/****************************/
void get_naa_natom(char *filename, char *chainID, int *naa, int *natom) {
  int iaa, i;
  int iatom;
  int sw, len, end;
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
      end = line_end(buffer);
      field_str(buffer, end, 12, 4, res);
      ch = field_chr(buffer, end, 21);
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
  int sw, len, end;

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
      end = line_end(buffer);
      field_str(buffer, end, 12, 4, res);
      ch = field_chr(buffer, end, 21);
      field_int(buffer, end, 22, 4, &org);
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
      ch = field_chr(buffer, line_end(buffer), 21);
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

  /* Same fields as the fsscanf() format
     "%6s%5d%1x%4s%1s%3s%1x%1s%4d%1s%3x%8f%8f%8f%6f%6f%6x%4s%2s%2s" */
  int end = line_end(buffer);
  field_str(buffer, end, 0, 6, record);
  field_int(buffer, end, 6, 5, &serial);
  field_str(buffer, end, 12, 4, name);
  field_str(buffer, end, 16, 1, altLoc);
  field_str(buffer, end, 17, 3, resName);
  field_str(buffer, end, 21, 1, chainID_org);
  field_int(buffer, end, 22, 4, &resSeq);
  field_str(buffer, end, 26, 1, iCode);
  field_flt(buffer, end, 30, 8, &x);
  field_flt(buffer, end, 38, 8, &y);
  field_flt(buffer, end, 46, 8, &z);
  field_flt(buffer, end, 54, 6, &occupancy);
  field_flt(buffer, end, 60, 6, &tempFactor);
  field_str(buffer, end, 72, 4, segID);
  field_str(buffer, end, 76, 2, element);
  field_str(buffer, end, 78, 2, charge);

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
