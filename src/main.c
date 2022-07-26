#include "mican.h"
#include "main.h"

/** Define global variables**/

INPUT_DATA input;
clock_t time_init, time_end, time_tmp;

void usage(void);

int main(int argc, char *argv[]) {

  /**  time count  **/
  time_init = clock();

  /*******************************/
  /*     default parameters      */
  /*******************************/
  strcpy(input.pdbout, "OFF"); //  default value for superposition file
  strcpy(input.aliout, "OFF"); //  default value for alignment file
  strcpy(input.matout, "OFF"); //  default value for translation matrix file
  input.mode = DEF_mode;       //  default value for alignment mode
  input.bstft_weight = OFF;    //  RMSD weight for superposition
  strcpy(input.chain_t, "0");  //  considering all chain
  strcpy(input.chain_q, "0");  //  considering all chain
  /****  for GH search  ****/
  input.max_hash = DEF_max_hash;       //  hash matrix size
  input.d_hash = DEF_d_hash;           //  hash block resolution
  input.max_nghbr = DEF_max_nghbr;     //  neighbor hash search
  input.cutoff_para = DEF_cutoff_para; //  SSE para vector cutoff
  input.cutoff_perp = DEF_cutoff_perp; //  SSE perp vector cutoff
  input.min_GHcount = DEF_min_GHcount; //  minimum voting score
  input.GHcutnum = DEF_GHcutnum;       //  number of candidates
  /****  for Greedy refinement  ****/
  input.type_select = SELE_mTMscore;   //  controll selection score
  input.smin = DEF_smin;               //  alignment score cutoff
  input.lmin = DEF_lmin;               //  segment length cutoff
  input.dmax = DEF_dmax;               //  maximum distance
  input.greedy_cut1 = DEF_greedy_cut1; //  greedy cutoff 1
  input.greedy_cut2 = DEF_greedy_cut2; //  greedy cutoff 2
  input.greedy_cut3 = DEF_greedy_cut3; //  greedy cutoff 3
  input.ssetype_w = DEF_ssetype_w;     //  sse type weight
  input.angle_sigma = DEF_angle_sigma; //  angle gaussian sigma
  input.clst_red = DEF_clst_red;       //  clustering cutoff (1st)
  input.clst_cutoff = DEF_clst_cutoff; //  clustering cutoff (last)
  input.max_refine = DEF_max_refine;   //  max iteration of refinement
  /****  for sub-optimal solutions  ****/
  input.nsub = DEF_nsub; //  number of sub-optimal alignments considered
  input.isub_out = 1;    //  output superposition of suboptimal alignment
  /****  others  ****/
  input.TM_d0 = DEF_TM_d0;             //  d0 value for TM-score
  input.d0fix = OFF;                   //  controll fixed TM_d0 value mode
  input.qtchange = OFF;  //  exchange : query <-> template
  input.silent = OFF;    //  silent mode
  input.progress = OFF;  //  output progress
  /****  invisible options  ****/
  input.printalign = OFF;
  
  /****************************/
  /*  GET OPTION & FILENAMES  */
  /****************************/
  if (getopt(argc, argv) == FALSE) {
    usage();
  }

  /***************************/
  /*     CHECK ARGUMENTS     */
  /***************************/
  if (input.isub_out < 1) {
    fprintf(stderr, "!! Error   : Invarid argument for option -- '-i'\n");
    exit(0);
  }
  if(input.nsub < input.isub_out){ input.nsub = input.isub_out; }
  input.nsub_org = input.nsub;

  /***************************/
  /**        MESSAGE        **/
  /***************************/
  if (input.silent == OFF) {
    printf("\n");
    if (input.mode == SEQUENT) {
      printf(" Alignment mode : sequential (SQ)\n\n");
    }
    if (input.mode == FORWARD) {
      printf(" Alignment mode : rewiring (RW)\n\n");
    }
    if (input.mode == FWandRV) {
      printf(" Alignment mode : rewiring & reverse (RR)\n\n");
    }
    if (input.mode == REVERSE) {
      printf(" Alignment mode : reverse constrained\n\n");
    }
  }

  /***************************/
  /**  STRUCTURE ALIGNMENT  **/
  /***************************/
  main_align();

  /***************************/
  /**      TOTAL CLOCK      **/
  /***************************/
  time_end = clock(); // time count
  if (input.silent == OFF) {
    printf("\n Finished successfully ( %4.3f sec )\n\n",
           (double)(time_end - time_init) / CLOCKS_PER_SEC);
  }

  return (0);
}

/********************************/
/**      OUTPUT    USAGE       **/
/********************************/
void usage(void) {
  printf(
      "\n"	 
      " USAGE: %% mican protein1 protein2 [OPTION]\n\n"
      " Description:\n"
      "  -f             fast mode (same as \"-g %d\")\n"
      "  -s             sequential (SQ) alignment mode\n"
      "  -w             rewiring (RW) alignment mode\n"
      "  -r             rewiring & reverse (RR) alignment mode\n"
      "  -R             reverse constrained alignment mode\n"
      "  -x             silent mode (without any output on the console)\n"
      "  -p             print alignment progress\n"
      "  -c1 ChainIDs   chain ID specifier for protein1 (e.g. -c1 A, -c1 ABC)\n"
      "  -c2 ChainIDs   chain ID specifier for protein2\n"
      "  -o  Filename   superposition file (rasmol-script)\n"
      "  -a  Filename   alignment file\n"
      "  -m  Filename   translation matrix file\n"
      "  -n  Integer    number of solutions output (default=%d)\n"
      "  -i  Integer    output i-th solution on stdout & superposition file\n"
      "  -t  Integer    selection score ([0]:sTMscore, 1:TMscore, 2:Dali-Z)\n"
      "  -g  Integer    number of GH candidates used (default=%d)\n"
      "  -l  Integer    minimum segment length (default=%d)\n"
      "  -d  Real       fix TM-score scaling factor d0\n"
      "  -q  Real       maximum distance between Ca atoms to be aligned (default=%3.1f)\n\n"
      " Simple usage (SQ):\n"
      "   %% mican protein1 protein2\n"
      "   %% mican protein1 protein2 -a align.aln -o sup.pdb\n\n"
      " Rewiring mode alignment (RW):\n"
      "   %% mican protein1 protein2 -w\n\n"
      " Rewiring & reverse mode alignment (RR):\n"
      "   %% mican protein1 protein2 -r\n\n"
      " To visualize superposition:\n"
      "   %% mican protein1 protein2 -o sup.pdb\n"
      "   %% rasmol -script sup.pdb\n\n",
      FAST_GHcutnum, DEF_nsub, DEF_GHcutnum, DEF_lmin, DEF_dmax);
  exit(1);
}

