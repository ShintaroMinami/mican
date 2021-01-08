#include "mican.h"

/**  function resname2reschar  **/

char resname2reschar(char *resname) {
  if (strcmp(resname, "ALA") == 0) {
    return 'A';
  } else if (strcmp(resname, "CYS") == 0) {
    return 'C';
  } else if (strcmp(resname, "ASP") == 0) {
    return 'D';
  } else if (strcmp(resname, "GLU") == 0) {
    return 'E';
  } else if (strcmp(resname, "PHE") == 0) {
    return 'F';
  } else if (strcmp(resname, "GLY") == 0) {
    return 'G';
  } else if (strcmp(resname, "HIS") == 0) {
    return 'H';
  } else if (strcmp(resname, "ILE") == 0) {
    return 'I';
  } else if (strcmp(resname, "LYS") == 0) {
    return 'K';
  } else if (strcmp(resname, "LEU") == 0) {
    return 'L';
  } else if (strcmp(resname, "MET") == 0) {
    return 'M';
  } else if (strcmp(resname, "ASN") == 0) {
    return 'N';
  } else if (strcmp(resname, "PRO") == 0) {
    return 'P';
  } else if (strcmp(resname, "GLN") == 0) {
    return 'Q';
  } else if (strcmp(resname, "ARG") == 0) {
    return 'R';
  } else if (strcmp(resname, "SER") == 0) {
    return 'S';
  } else if (strcmp(resname, "THR") == 0) {
    return 'T';
  } else if (strcmp(resname, "VAL") == 0) {
    return 'V';
  } else if (strcmp(resname, "TRP") == 0) {
    return 'W';
  } else if (strcmp(resname, "TYR") == 0) {
    return 'Y';
  } else {
    return 'X';
  }
}
