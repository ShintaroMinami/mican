#include "mican.h"
#include "fsscanf.h"

PRIVATE void get_xyz(const ATOM_REC *a, int naa, int *iaa, int *iatom,
                    RESDAT *resdat, ALLATM *allatm, int chainID);

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

/* Parses an ATOM/HETATM line; the same fields as the fsscanf() format
   "%6s%5d%1x%4s%1s%3s%1x%1s%4d%1s%3x%8f%8f%8f%6f%6f%6x%4s%2s%2s". */
PRIVATE void parse_atom_line(const char *buffer, ATOM_REC *a) {
  int end = line_end(buffer);
  memset(a, 0, sizeof(ATOM_REC));
  field_str(buffer, end, 0, 6, a->record);
  field_int(buffer, end, 6, 5, &a->serial);
  field_str(buffer, end, 12, 4, a->name);
  field_str(buffer, end, 17, 3, a->resName);
  a->chain[0] = field_chr(buffer, end, 21);
  field_int(buffer, end, 22, 4, &a->resSeq);
  field_flt(buffer, end, 30, 8, &a->x);
  field_flt(buffer, end, 38, 8, &a->y);
  field_flt(buffer, end, 46, 8, &a->z);
  field_flt(buffer, end, 54, 6, &a->occupancy);
  field_flt(buffer, end, 60, 6, &a->tempFactor);
}

/* Reads the ATOM/HETATM records up to the first ENDMDL into a new array
   *atoms (0-based) and returns their number. */
PRIVATE int read_pdb_atoms(const char *filename, ATOM_REC **atoms) {
  FILE *fp;
  char buffer[BUFSIZE + 1];
  int natom = 0, cap = 0;
  ATOM_REC *a = NULL;

  if ((fp = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "!! Error   : No such file (%s)\n", filename);
    exit(1);
  }
  while (fgets(buffer, BUFSIZE, fp) != NULL) {
    if (strncmp(buffer, "ENDMDL", 6) == 0) {
      break;
    }
    if (strncmp(buffer, "ATOM ", 5) == 0 || strncmp(buffer, "HETATM", 6) == 0) {
      if (natom >= cap) {
        cap = 2 * cap + 1024;
        a = (ATOM_REC *)realloc(a, (size_t)cap * sizeof(ATOM_REC));
        if (a == NULL) {
          fprintf(stderr, "!! Error   : memory allocation failed\n");
          exit(1);
        }
      }
      parse_atom_line(buffer, &a[natom]);
      natom++;
    }
  }
  fclose(fp);
  *atoms = a;
  return (natom);
}

/*
 * Atoms of a structure file (PDB or mmCIF), parsed once and kept for the
 * three passes below (get_naa_natom, check_CA, read_pdb) over both inputs.
 */
typedef struct {
  char *filename;
  int natom;
  ATOM_REC *atom;
} ATOM_FILE;

#define N_ATOM_CACHE 2
PRIVATE ATOM_FILE atom_cache[N_ATOM_CACHE];

PRIVATE const ATOM_FILE *load_atoms(const char *filename) {
  static int next = 0;
  ATOM_FILE *af;
  int i;

  for (i = 0; i < N_ATOM_CACHE; i++) {
    if (atom_cache[i].filename != NULL &&
        strcmp(atom_cache[i].filename, filename) == 0) {
      return (&atom_cache[i]);
    }
  }
  af = &atom_cache[next];
  next = (next + 1) % N_ATOM_CACHE;
  free(af->filename);
  free(af->atom);
  af->filename = (char *)malloc(strlen(filename) + 1);
  if (af->filename == NULL) {
    fprintf(stderr, "!! Error   : memory allocation failed\n");
    exit(1);
  }
  strcpy(af->filename, filename);
  if (is_mmcif_file(filename)) {
    af->natom = read_mmcif_atoms(filename, &af->atom);
  } else {
    af->natom = read_pdb_atoms(filename, &af->atom);
  }
  return (af);
}

/* Whether atom a belongs to the selected chains (chainID "0": all chains).
   Chain IDs longer than one character (mmCIF only) cannot be selected and
   stop the program when all chains are requested. */
PRIVATE int chain_selected(const ATOM_REC *a, const char *chainID,
                           const char *filename) {
  int i, len;

  if (a->chain[1] != '\0') {
    if (strcmp(chainID, "0") == 0) {
      fprintf(stderr,
              "!! Error   : Multi-character chain ID (%s) is not supported "
              "(%s); select chains with -c1/-c2\n",
              a->chain, filename);
      exit(1);
    }
    return (FALSE);
  }
  if (strcmp(chainID, "0") == 0) {
    return (TRUE);
  }
  len = (int)strlen(chainID);
  for (i = 0; i <= len; i++) {
    if (a->chain[0] == chainID[i]) {
      return (TRUE);
    }
  }
  return (FALSE);
}

/****************************/
/*    FUNCTION  get_naa     */
/****************************/
void get_naa_natom(char *filename, char *chainID, int *naa, int *natom) {
  const ATOM_FILE *af = load_atoms(filename);
  const ATOM_REC *a;
  int iaa = 0, iatom = 0, k;

  /** count CA atoms **/
  for (k = 0; k < af->natom; k++) {
    a = &af->atom[k];
    if (!chain_selected(a, chainID, filename)) {
      continue;
    }
    iatom++;
    if (strncmp(a->name, " CA ", 4) == 0) {
      iaa++;
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

  return;
}

/*****************************/
/*                           */
/*    FUNCTION  check_CA     */
/*                           */
/*****************************/
void check_CA(char *filename, char *chainID, RESDAT *resdat) {
  const ATOM_FILE *af = load_atoms(filename);
  const ATOM_REC *a;
  int iaa = 0, k;

  /** get iaa_org **/
  for (k = 0; k < af->natom; k++) {
    a = &af->atom[k];
    if (!chain_selected(a, chainID, filename)) {
      continue;
    }
    if (strncmp(a->name, " CA ", 4) == 0) {
      iaa++;
      resdat[iaa].iaa_org = a->resSeq;
    }
  }

  return;
}

/*****************************/
/*                           */
/*    FUNCTION  read_pdb     */
/*                           */
/*****************************/
void read_pdb(char *filename, char *chainID, int naa, int natom, RESDAT *resdat,
              ALLATM *allatm) {
  const ATOM_FILE *af = load_atoms(filename);
  const ATOM_REC *a;
  int iaa, iatom, k;
  char list_chainID[MX_CHAIN + 1];
  int i, ichain, nchain;

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
  /*   read ATOM record  */
  /***********************/
  nchain = 0;
  ichain = 0;
  for (k = 0; k < af->natom; k++) {
    a = &af->atom[k];
    if (!chain_selected(a, chainID, filename)) {
      continue;
    }
    if (nchain == 0) {
      nchain = 1;
      list_chainID[nchain] = a->chain[0];
    } else {
      ichain = 0;
      for (i = 1; i <= nchain; i++) {
        if (a->chain[0] == list_chainID[i]) {
          ichain = i;
          break;
        }
      }
      if (ichain == 0) {
        nchain++;
        list_chainID[nchain] = a->chain[0];
        ichain = nchain;
      }
    }
    get_xyz(a, naa, &iaa, &iatom, resdat, allatm, ichain);
  }

  return;
}

/*****************************/
/*                           */
/*        SUBROUTINE         */
/*                           */
/*****************************/
PRIVATE void get_xyz(const ATOM_REC *a, int naa, int *iaa, int *iatom,
                    RESDAT *resdat, ALLATM *allatm, int chainID) {
  const char *name = a->name;
  const char *resName = a->resName;
  int resSeq = a->resSeq;
  float x = a->x;
  float y = a->y;
  float z = a->z;
  char chainID_org[1 + 1];

  /**  check chainID  **/
  chainID_org[0] = (a->chain[0] == ' ') ? 'A' : a->chain[0];
  chainID_org[1] = '\0';
  /**************************************/
  /**         read  CA  atoms          **/
  /**************************************/
  if (strncmp(name, " CA ", 4) == 0) {
    *iaa = *iaa + 1;
    /** Residue Name **/
    strcpy(resdat[*iaa].resname, resName);
    /** Residue Character **/
    resdat[*iaa].reschar = resname2reschar((char *)resName);
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
  strcpy(allatm[*iatom].record, a->record);
  /** iatom original **/
  allatm[*iatom].iatom_org = a->serial;
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
  allatm[*iatom].occupancy = a->occupancy;
  /** Temp Factor **/
  allatm[*iatom].tempfactor = a->tempFactor;
  return;
}
