#include "mican.h"
#include "main.h"

PRIVATE void read_opt(int argc, char *argv[], int iarg, int ich);
PRIVATE int check_argument(int argc, char *argv[], int iarg, int ich);

int getopt(int argc, char *argv[]) {
  int iarg;
  int ich;
  int iexist;

  /**  get option  **/
  for (iarg = 1; iarg <= argc - 1; iarg++) {
    if (argv[iarg][0] == '-') {
      ich = 1;
      while (argv[iarg][ich]) {
        read_opt(argc, argv, iarg, ich);
        ich++;
      }
      argv[iarg][0] = '\0'; // delete
    }
  }

  /**  read PDB filename  **/
  iexist = 0;
  for (iarg = 1; iarg <= argc - 1; iarg++) {
    if (argv[iarg][0] != '\0') {
      iexist++;
      if (iexist == 1) {
        strcpy(input.file_t, argv[iarg]);
      }
      if (iexist == 2) {
        strcpy(input.file_q, argv[iarg]);
      }
    }
  }

  if (iexist != 2) {
    return (FALSE);
  }
  return (TRUE);
} /** end of function **/

/** function read_opt **/
PRIVATE void read_opt(int argc, char *argv[], int iarg, int ich) {
  switch (argv[iarg][ich]) {
  case 'w': /**  rewiring mode  **/
    input.mode = FORWARD;
    break;
  case 'r': /**  rewiring & reverse mode  **/
    input.mode = FWandRV;
    break;
  case 's': /**  sequential mode  **/
    input.mode = SEQUENT;
    break;
  case 'R': /**  reverse only mode  **/
    input.mode = REVERSE;
    break;
  case 'f': /**  fast mode  **/
    input.GHcutnum = FAST_GHcutnum;
    break;
  case 'x': /**  silent mode  **/
    input.silent = ON;
    break;
  case 'p': /**  print progress  **/
    input.progress = ON;
    break;
  case 'c': /**  chain ID  **/
    if (check_argument(argc, argv, iarg, ich) == TRUE) {
      if (strncmp(argv[iarg], "-c1", 3)) {
        strcpy(input.chain_q, argv[iarg + 1]);
      } else if (strncmp(argv[iarg], "-c2", 3)) {
        strcpy(input.chain_t, argv[iarg + 1]);
      }
      argv[iarg][2] = '\0';
      argv[iarg + 1][0] = '\0';
    }
    break;
  case 'o': /**  output superposition filename  **/
    if (check_argument(argc, argv, iarg, ich) == TRUE) {
      strcpy(input.pdbout, argv[iarg + 1]);
      argv[iarg + 1][0] = '\0';
    }
    break;
  case 'a': /**  output alignment filename  **/
    if (check_argument(argc, argv, iarg, ich) == TRUE) {
      strcpy(input.aliout, argv[iarg + 1]);
      argv[iarg + 1][0] = '\0';
    }
    break;
  case 'm': /**  output translation matrix filename  **/
    if (check_argument(argc, argv, iarg, ich) == TRUE) {
      strcpy(input.matout, argv[iarg + 1]);
      argv[iarg + 1][0] = '\0';
    }
    break;
  case 'n': /**  nsub  **/
    if (check_argument(argc, argv, iarg, ich) == TRUE) {
      input.nsub = atoi(argv[iarg + 1]);
      argv[iarg + 1][0] = '\0';
    }
    break;
  case 'i': /**  output i-th sub-optimal alignment  **/
    if (check_argument(argc, argv, iarg, ich) == TRUE) {
      input.isub_out = atoi(argv[iarg + 1]);
      argv[iarg + 1][0] = '\0';
    }
    break;
  case 'g': /**  number of GH candidates  **/
    if (check_argument(argc, argv, iarg, ich) == TRUE) {
      input.GHcutnum = atoi(argv[iarg + 1]);
      argv[iarg + 1][0] = '\0';
    }
    break;
  case 'd': /**  TM_d0  **/
    if (check_argument(argc, argv, iarg, ich) == TRUE) {
      input.TM_d0 = (float)atof(argv[iarg + 1]);
      input.d0fix = ON;
      argv[iarg + 1][0] = '\0';
    }
    break;
  case 'h': /**  hash resolution  **/
    if (check_argument(argc, argv, iarg, ich) == TRUE) {
      input.d_hash = (float)atof(argv[iarg + 1]);
      argv[iarg + 1][0] = '\0';
    }
    break;
  case 'H': /**  max hash size  **/
    if (check_argument(argc, argv, iarg, ich) == TRUE) {
      input.max_hash = atoi(argv[iarg + 1]);
      argv[iarg + 1][0] = '\0';
    }
    break;
  case 'l': /**  minimum segment length  **/
    if (check_argument(argc, argv, iarg, ich) == TRUE) {
      input.lmin = atoi(argv[iarg + 1]);
      argv[iarg + 1][0] = '\0';
    }
    break;
  case 't': /**  selection score type  **/
    if (check_argument(argc, argv, iarg, ich) == TRUE) {
      input.type_select = atoi(argv[iarg + 1]);
      argv[iarg + 1][0] = '\0';
    }
    break;
  case 'u': /**  clustering cutoff for output alignments  **/
    if (check_argument(argc, argv, iarg, ich) == TRUE) {
      input.clst_cutoff = atof(argv[iarg + 1]);
      argv[iarg + 1][0] = '\0';
    }
    break;
  case 'q': /**  maximum distance of aligned Ca-Ca pairs  **/
    if (check_argument(argc, argv, iarg, ich) == TRUE) {
      input.dmax = atof(argv[iarg + 1]);
      argv[iarg + 1][0] = '\0';
    }
    break;
  /****  invisible options  ****/
  case 'z': /**  print alignment on STDOUT  **/
    input.printalign = ON;
    break;
  /****  for parameter tuning  ****/
  case '0': /**  Smin  **/
    if (check_argument(argc, argv, iarg, ich) == TRUE) {
      input.smin = atof(argv[iarg + 1]);
      argv[iarg + 1][0] = '\0';
    }
    break;
  case '1': /**  greedy cut1  **/
    if (check_argument(argc, argv, iarg, ich) == TRUE) {
      input.greedy_cut1 = atof(argv[iarg + 1]);
      argv[iarg + 1][0] = '\0';
    }
    break;
  case '2': /**  greedy cut2  **/
    if (check_argument(argc, argv, iarg, ich) == TRUE) {
      input.greedy_cut2 = atof(argv[iarg + 1]);
      argv[iarg + 1][0] = '\0';
    }
    break;
  case '3': /**  greedy cut3  **/
    if (check_argument(argc, argv, iarg, ich) == TRUE) {
      input.greedy_cut3 = atof(argv[iarg + 1]);
      argv[iarg + 1][0] = '\0';
    }
    break;
  default: /**  unknown option  **/
    fprintf(stderr, "!! Warning : Invalid option --'%c'\n", argv[iarg][ich]);
    break;
  }

  return;
} /** end of function **/

/** function check argument  **/
PRIVATE int check_argument(int argc, char *argv[], int iarg, int ich) {
  if (iarg == argc - 1) {
    fprintf(stderr, "!! Warning : The option requires an argument -- '%c'\n",
            argv[iarg][ich]);
    return (FALSE);
  } else if (argv[iarg + 1][0] == '-') {
    fprintf(stderr, "!! Warning : The option requires an argument -- '%c'\n",
            argv[iarg][ich]);
    return (FALSE);
  }
  return (TRUE);
} /** end of function **/
