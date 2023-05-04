#ifndef INCLUDED_mican_h_
#define INCLUDED_mican_h_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <time.h>

#define PRIVATE static

/****  numerical constants  ****/
#define EPS 0.000001
#define PI 3.141592

/****  controll  ****/
#define TRUE 1
#define FALSE 0
#define ON 1
#define OFF 0
#define SEQUENT 2
#define FORWARD 1
#define FWandRV 0
#define REVERSE -1
#define ALICHAIN 0
#define ALLCHAIN 1

/****  string size  ****/
#define BUFSIZE 256
#define STRSIZE 256
#define MX_FILENAME 1024

/****  protein parameters  ****/
#define MX_AA_TYPE 20
#define MX_ATOM_TYPE 5
#define MX_CHAIN 500

/****  atom type  ****/
#define atom_N_ 1
#define atom_CA 2
#define atom_C_ 3
#define atom_O_ 4
#define atom_CB 5

/****  coordinate dimention  ****/
#define Ncoord 3

/****  SSE parameters  ****/
#define MX_SSE 10000
#define MX_seg 6
#define SSE_COIL 0
#define SSE_HELIX_ 1
#define SSE_STRAND 2
#define seg_STRAND 3
#define seg_HELIX_ 6

/****  max iteration number  ****/
#define Nite_loc 3

/****  parameters of TM-score  ****/
#define TM_len_MIN 15
#define TM_d0_MAX 10.0f
#define TM_d0_MIN 1.0f

/****  parameters of SP-score  ****/
#define SP_d0 4.0f
#define SP_dcore 8.0f
#define SP_dsurr 12.0f
#define SP_const 0.2f
#define SP_power 0.7f

/****  parameters of Daliscore  ****/
#define Dali_theteE   0.2f
#define Dali_alpha    400.0f
#define Dali_Lmax     400
/** Parameters can be found in DaliLite_3.3 software (at l518 in dp.f) **/
#define Dali_c1       7.9494f         // 
#define Dali_c2       0.70852f        //
#define Dali_c3       2.5895f/10000   // -2.59f/10000 in Holm & Sander (1998)
#define Dali_c4      -1.9156f/1000000 //

/****  selection score  ****/
#define SELE_mTMscore  0
#define SELE_TMscore   1
#define SELE_Daliscore 2
#define SELE_SPscore   3
#define SELE_aTMscore  4



/******************************/
/**    default parameters    **/
/******************************/
/****  alignment mode  ****/
#define DEF_mode 2 //  default value for alignment mode
/****  for GH search  ****/
#define DEF_max_hash 200     //  hash matrix size
#define DEF_d_hash 3.2f      //  hash block resolution
#define DEF_max_nghbr 1      //  neighbor hash search
#define DEF_cutoff_para 0.5f //  SSE para vector cutoff cos( PI/3 RAD )
#define DEF_cutoff_perp 0.5f //  SSE perp vector cutoff cos( PI/3 RAD )
#define DEF_min_GHcount 0.0f //  minimum GH boting score
#define DEF_GHcutnum 50      //  number of GH candidates
#define FAST_GHcutnum 15     //  number of GH candidates for fast mode
/****  for DP refinement  ****/
#define DEF_smin 0.88f        //  DP alignment score cutoff
#define DEF_lmin 3            //  segment length cutoff
#define DEF_dmax 10.0f        //  maximum distance
#define DEF_greedy_cut1 0.41f //  parameter: greedy cutoff 1
#define DEF_greedy_cut2 0.070f//  parameter: greedy cutoff 2
#define DEF_greedy_cut3 0.007f//  parameter: greedy cutoff 3
#define DEF_clst_red 0.95f    //  clustering cutoff : redundancy check
#define DEF_clst_cutoff 0.6f  //  clustering cutoff : result
#define DEF_max_refine 10     //  maximum iteration number of refinement
/****  for score  ****/
#define DEF_TM_d0 3.6f        //  TMscore d0 value
#define DEF_ssetype_w 1.0f    //  score weight for SSE type matching
#define DEF_angle_sigma 1.570796f //  PI/2
#define DEF_angle_w_null 0.641180f //  score with angle=PI/3
/****  for sub-optimal  ****/
#define DEF_nsub 5 //  number of sub-optimal alignment

typedef struct {
  /****  for file name  ****/
  char file_q[MX_FILENAME];
  char file_t[MX_FILENAME];
  char pdbout[MX_FILENAME];
  char aliout[MX_FILENAME];
  char matout[MX_FILENAME];
  /****  for chain identification  ****/
  char chain_t[STRSIZE];
  char chain_q[STRSIZE];
  /****  for controlling output alignment  ****/
  int qtchange;
  /****  for controling alignment mode  ****/
  int mode;
  int bstft_weight;
  /****  for GH hash search  ****/
  int max_hash;
  int max_nghbr;
  float d_hash;
  /****  for SSE comparison  ****/
  float cutoff_para;
  float cutoff_perp;
  float min_GHcount;
  int GHcutnum;
  /****  for alignment  ****/
  float smin;
  float lmin;
  float dmax;
  float greedy_cut1;
  float greedy_cut2;
  float greedy_cut3;
  float angle_sigma;
  float ssetype_w;
  float clst_red;
  float clst_cutoff;
  float max_refine;
  /****  for score  ****/
  int type_select;
  int d0fix;
  float TM_d0;
  float naa_eff;
  /****  for sub-optimal  ****/
  int nsub, nsub_org;
  int isub_out;
  /****  output option  ****/
  int silent;
  int progress;
  /****  invisible option  ****/
  int printalign;
} INPUT_DATA;


typedef struct {
  int exist;
  float x;
  float y;
  float z;
  char name[4 + 1];
} ATOM;

typedef struct {
  char resname[3 + 1];
  char chainID_org[1 + 1];
  int chainID;
  char reschar;
  char ssechar;
  int aa_type;
  int ssetype;
  int iaa_org;
  float cacoord[3 + 1];
  int chain_break;
  ATOM atom[MX_ATOM_TYPE + 1];
} RESDAT;

typedef struct {
  char record[6 + 1];
  char resname[3 + 1];
  char chainID_org[1 + 1];
  char atmname[4 + 1];
  int iatom_org;
  int iaa_org;
  float coord[3 + 1];
  float occupancy;
  float tempfactor;
} ALLATM;

typedef struct {
  int nchain;
  char **chainID_org;
  int *chainID;
  int *naa;
} PDBDAT;

typedef struct {
  char ssechar;
  int ssetype;
  float mid[3 + 1];
  float para[3 + 1];
  float perp[3 + 1];
  int naa;
  int iaa[MX_seg + 1];
} SSEDAT;

typedef struct {
  int nlist;
  int *isse;
  int *ssetype;
  float **para;
  float **perp;
} HASH_TABLE;

typedef struct {
  int isse_q;
  int isse_t;
  float count;
  int coord_mode;
} HASH_RESULT;

typedef struct {
  int isse_q;
  int isse_t;
  float count;
  int check;
  int coord_mode;

  float **rot;
  float *vec;

  int *align_q;
  int *align_t;
  int *segmode;
  int *termini;
  float *local_score;
  float *local_TM;
  float *local_mTM;
  float *local_aTM;
  float *local_Dali;
  float *local_dist;

  int naa_align;
  int naa_align_SSE;

  float cover_q;
  float cover_t;
  float cover_SSE_q;
  float cover_SSE_t;

  float rmsd;
  float TMscore_q;
  float TMscore_t;
  float TMscore_mod;
  float TMscore_ang;
  float Daliscore;
  float DaliZ;
  float SPscore;
  float C3score;
  float TMloc_mod;
  float score_cmp;

  float seqID;
} ALIGN;


/**********************/
/**     FUNCTION      */
/**********************/

/**  Main Controll Function  **/
void main_align(void);
int getopt(int argc, char *argv[]);

/**  Read PDB  **/
void get_naa_natom(char *filename, char *chainID, int *naa, int *natom);
void check_CA(char *filename, char *chainID, RESDAT *resdat);
void read_pdb(char *filename, char *chainID, int naa, int natom, RESDAT *resdat,
              ALLATM *allatm);
char resname2reschar(char *resname);
/**  Distance Matrix  **/
void calc_distmat(int naa, RESDAT *resdat, float **distmat);
/**  Check PDB  **/
void make_pdbdat(int naa, RESDAT *resdat, PDBDAT *pdbdat);
int count_chain(int naa, RESDAT *resdat);
void chain_break_check(int naa, RESDAT *resdat);

/**  SSE definition  **/
int assign_sse(int naa, SSEDAT *ssedat, RESDAT *resdat, int coil_as_strand);
void xyz2vec(int nsse, SSEDAT *ssedat, RESDAT *resdat);
void coord_definition(float *x, float *y, float *z, SSEDAT *ssedat, int isse);
void coord_definition_r(float *x, float *y, float *z, SSEDAT *ssedat, int isse);
float calc_maxdist(int nsse, SSEDAT *ssedat);

/**  Geometric Hashing  **/
void pre_make_hash(int nsse, int nhash, HASH_TABLE ***hash, SSEDAT *ssedat);
void make_hash_table(int nsse, int nhash, HASH_TABLE ***hash, SSEDAT *ssedat);
void search_hash(int nsse_q, int nsse_t, int nalign, int nhash,
                 HASH_TABLE ***hash, SSEDAT *ssedat, HASH_RESULT *hres);
void check_GHresult(int nsse, HASH_RESULT *hres, ALIGN *align, SSEDAT *ssedat_q,
                    SSEDAT *ssedat_t);
int hres2trans(SSEDAT *ssedat_q, SSEDAT *ssedat_t, ALIGN align);

/**  Alignment  **/
int ca_alignment(int naa_q, int naa_t, RESDAT *resdat_q, RESDAT *resdat_t,
                 ALIGN align, float ***mat, float **score, float **score_r);
void global_alignment(int naa1, int naa2, float ***mat, int *align, int *segmode,
                      float **score, float **score_r, RESDAT *resdat_q,
                      RESDAT *resdat_t);
void extend_termini(int naa1, int naa2, ALIGN align, float **co1, float **co2);
int superposition(int naa_q, RESDAT *resdat_q, RESDAT *resdat_t, ALIGN align, int sw);
void clustering_align(int naa_q, int naa_t, int nalign, ALIGN *align,
                      float cutoff);
int check_align(int naa_q, int naa_t, ALIGN *align, RESDAT *resdat_q, RESDAT *resdat_t);

/**  Score Function  **/
void score_calculation(int naa_q, int naa_t, int nsse_t, ALIGN *align,
                       RESDAT *resdat1, RESDAT *resdat2);
float TMscore(int naa0, int naa1, RESDAT *resdat1, RESDAT *resdat2, ALIGN align);
float TMscore_mod(int naa0, int naa1, int naa2, RESDAT *resdat1, RESDAT *resdat2,
		  ALIGN align, int angle_switch, int local_switch);
float SPscore(int naa1, int naa2, RESDAT *resdat1, RESDAT *resdat2,
              ALIGN align);
float Daliscore(int naa1, float **distmat1, float **distmat2, ALIGN *align);
float rmsd(int naa1, RESDAT *resdat1, RESDAT *resdat2, ALIGN align);
float C3score(int nch1, int nch2, int nalign, ALIGN **align_chain);
void calc_seq_identity(int naa_q, int naa_t, int nalign, ALIGN *align,
                       RESDAT *resdat_q, RESDAT *resdat_t);
void score_calc_chain(int naa_q, int naa_t, ALIGN align, ALIGN **align_chain,
                      RESDAT *resdat_q, RESDAT *resdat_t, PDBDAT pdbdat_q,
                      PDBDAT pdbdat_t);

/**  Output Function  **/
void print_progress(char *progress_name, int numerator, int denominator);
void output(int naa_q, int naa_t, int natom_q, int natom_t, RESDAT *resdat_q,
            RESDAT *resdat_t, ALLATM *allatm_q, ALLATM *allatm_t, ALIGN *align,
            ALIGN ***align_chain, PDBDAT pdbdat_q, PDBDAT pdbdat_t);
void printinp(int naa_q, int naa_t, PDBDAT pdbdat_q, PDBDAT pdbdat_t);
void printsup(FILE *fp, int naa_q, int naa_t, int natom_q, int natom_t,
              RESDAT *resdat_q, RESDAT *resdat_t, ALLATM *allatm_q,
              ALLATM *allatm_t, ALIGN *align);
void printrascript(FILE *fp, int naa1, int naa2, RESDAT *resdat1,
                   RESDAT *resdat2, int *align1, int *align2);
void printallatm(FILE *fp, int natom, ALLATM *allatm, int imodel);
void printali(FILE *fp, int naa_q, int naa_t, RESDAT *resdat_q,
                 RESDAT *resdat_t, ALIGN *align, int header);
void printmat(FILE *fp, ALIGN *align, int header);
void inverse_mat(float **rot, float *vec);

/**  Utils  **/
void bstft(int n, RESDAT *resdat1, RESDAT *resdat2, float *vec, float **rot);
void wbstft(int n, RESDAT *resdat1, RESDAT *resdat2, float *vec, float **rot, float *weight);
int svd(int m, int n, int withu, int withv, double eps, double tol, double *a,
        double *q, double *u, double *v, double *vt);
float angle(float *a1, float *a2, float *a3, float *b1, float *b2, float *b3);
float distance_2(float *a, float *b);

#endif
