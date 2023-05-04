#include "mican.h"
#include "main.h"

void copy_align(int naa_q, int naa_t, ALIGN *align1, ALIGN *align2);

void main_align(void) {
  int i, j, k;
  int iaa;
  int naa_q, naa_t, naa_tmp, naa_max;
  int natom_q, natom_t, natom_tmp;
  int ichain;
  int nsse_q, nsse_t;
  int ialign, nalign, maxalign;
  int ilist;
  int ntrue, itrue;
  int irefine;
  int nhash;
  float max_dist;
  char fname_tmp[STRSIZE];

  RESDAT *resdat_q;
  RESDAT *resdat_t;
  ALLATM *allatm_q;
  ALLATM *allatm_t;
  SSEDAT *ssedat_q;
  SSEDAT *ssedat_t;

  HASH_TABLE ***hash;
  HASH_RESULT *hres;
  ALIGN *align;
  ALIGN align_tmp;
  float ***mat;
  float **score, **score_r;

  PDBDAT pdbdat_q;
  PDBDAT pdbdat_t;
  ALIGN ***align_chain;

  float **distmat_q, **distmat_t;

  /********************************/
  /**      coordinate setup      **/
  /********************************/
  /**  count CA atom  **/
  get_naa_natom(input.file_q, input.chain_q, &naa_q, &natom_q);
  get_naa_natom(input.file_t, input.chain_t, &naa_t, &natom_t);

  /* define naa_max */
  if (naa_q > naa_t) {
    naa_max = naa_q;
  } else {
    naa_max = naa_t;
  }

  /**  exchange  query & template  **/
  if (naa_t < naa_q) {
    naa_tmp = naa_q;
    naa_q = naa_t;
    naa_t = naa_tmp;
    natom_tmp = natom_q;
    natom_q = natom_t;
    natom_t = natom_tmp;
    strcpy(fname_tmp, input.file_q);
    strcpy(input.file_q, input.file_t);
    strcpy(input.file_t, fname_tmp);
    strcpy(fname_tmp, input.chain_q);
    strcpy(input.chain_q, input.chain_t);
    strcpy(input.chain_t, fname_tmp);
    input.qtchange = ON;
  }

  /**  calculate TM_d0 value  **/
  if (input.d0fix == OFF) {
    if (naa_q <= TM_len_MIN) {
      input.TM_d0 = TM_d0_MIN;
    } else {
      input.TM_d0 = (float)(1.24 * pow(naa_q - TM_len_MIN, 1.0 / 3.0) - 1.8);
    }
    if (input.TM_d0 < TM_d0_MIN) {
      input.TM_d0 = TM_d0_MIN;
    }
    if (input.TM_d0 > TM_d0_MAX) {
      input.TM_d0 = TM_d0_MAX;
    }
  }

  /**  naa effective (for Dali Zscore)  **/
  input.naa_eff = sqrt((double)naa_q * (double)naa_t);

  /**  calloc RESDAT  **/
  resdat_q = (RESDAT *)calloc((size_t)(naa_q + 1), sizeof(RESDAT));
  resdat_t = (RESDAT *)calloc((size_t)(naa_t + 1), sizeof(RESDAT));

  /**  calloc all ATOM  **/
  allatm_q = (ALLATM *)calloc((size_t)(natom_q + 1), sizeof(ALLATM));
  allatm_t = (ALLATM *)calloc((size_t)(natom_t + 1), sizeof(ALLATM));

  /**  check CA iaa_org  **/
  check_CA(input.file_q, input.chain_q, resdat_q);
  check_CA(input.file_t, input.chain_t, resdat_t);

  /**  read pdb file  **/
  read_pdb(input.file_q, input.chain_q, naa_q, natom_q, resdat_q, allatm_q);
  read_pdb(input.file_t, input.chain_t, naa_t, natom_t, resdat_t, allatm_t);

  /**  read chain data  **/
  pdbdat_q.nchain = count_chain(naa_q, resdat_q);
  pdbdat_t.nchain = count_chain(naa_t, resdat_t);

  pdbdat_q.chainID = (int *)calloc((size_t)(pdbdat_q.nchain + 1), sizeof(int));
  pdbdat_t.chainID = (int *)calloc((size_t)(pdbdat_t.nchain + 1), sizeof(int));

  pdbdat_q.chainID_org =
      (char **)calloc((size_t)(pdbdat_q.nchain + 1), sizeof(char *));
  pdbdat_t.chainID_org =
      (char **)calloc((size_t)(pdbdat_t.nchain + 1), sizeof(char *));

  for (ichain = 0; ichain <= pdbdat_q.nchain; ichain++) {
    pdbdat_q.chainID_org[ichain] = (char *)calloc(2, sizeof(char));
  }
  for (ichain = 0; ichain <= pdbdat_t.nchain; ichain++) {
    pdbdat_t.chainID_org[ichain] = (char *)calloc(2, sizeof(char));
  }

  pdbdat_q.naa = (int *)calloc((size_t)(pdbdat_q.nchain + 1), sizeof(int));
  pdbdat_t.naa = (int *)calloc((size_t)(pdbdat_t.nchain + 1), sizeof(int));

  for (ichain = 0; ichain <= pdbdat_q.nchain; ichain++) {
    pdbdat_q.naa[ichain] = 0;
  }
  for (ichain = 0; ichain <= pdbdat_t.nchain; ichain++) {
    pdbdat_t.naa[ichain] = 0;
  }

  make_pdbdat(naa_q, resdat_q, &pdbdat_q);
  make_pdbdat(naa_t, resdat_t, &pdbdat_t);

  /*******************************/
  /**     print INPUT data      **/
  /*******************************/
  if (input.silent == OFF) {
    printinp(naa_q, naa_t, pdbdat_q, pdbdat_t);
  }

  /*******************************/
  /**      SSE assignment       **/
  /*******************************/
  ssedat_q = (SSEDAT *)calloc((size_t)(naa_max + 1), sizeof(SSEDAT));
  ssedat_t = (SSEDAT *)calloc((size_t)(naa_max + 1), sizeof(SSEDAT));

  /**  assign sse  **/
  nsse_q = assign_sse(naa_q, ssedat_q, resdat_q, FALSE);
  nsse_t = assign_sse(naa_t, ssedat_t, resdat_t, FALSE);

  if (nsse_q < 2 || nsse_t < 2) {
    if (nsse_q < 2) {
      fprintf(stderr, "! Worning   : nsse < 2 in %s\n", input.file_q);
    }
    if (nsse_t < 2) {
      fprintf(stderr, "! Worning   : nsse < 2 in %s\n", input.file_t);
    }
    fprintf(stderr, "! Worning   : Every COIL will be considered as STRAND\n");
    nsse_q = assign_sse(naa_q, ssedat_q, resdat_q, TRUE);
    nsse_t = assign_sse(naa_t, ssedat_t, resdat_t, TRUE);
  }

  /*******************************/
  /**  calc distance matrices   **/
  /*******************************/
  distmat_q = (float **)calloc((size_t)(naa_q+1), sizeof(float *));
  distmat_t = (float **)calloc((size_t)(naa_t+1), sizeof(float *));
  for(iaa=0; iaa<=naa_q; iaa++){
    distmat_q[iaa] = (float *)calloc((size_t)(naa_q+1), sizeof(float));
  }
  for(iaa=0; iaa<=naa_t; iaa++){
    distmat_t[iaa] = (float *)calloc((size_t)(naa_t+1), sizeof(float));
  }
  calc_distmat(naa_q, resdat_q, distmat_q);
  calc_distmat(naa_t, resdat_t, distmat_t);

  /*******************************/
  /**     xyz to sse vector     **/
  /*******************************/
  xyz2vec(nsse_q, ssedat_q, resdat_q);
  xyz2vec(nsse_t, ssedat_t, resdat_t);

  /*******************************/
  /**      make HASH_TABLE      **/
  /*******************************/
  /**  hash matrix size  **/
  max_dist = calc_maxdist(nsse_q, ssedat_q);
  nhash = 2 * (int)(max_dist / input.d_hash + 1);
  if (nhash > input.max_hash) {
    nhash = input.max_hash;
  }

  /**  calloc HASH_TABLE 1  **/
  hash = (HASH_TABLE ***)calloc((size_t)(nhash + 1), sizeof(HASH_TABLE **));
  for (i = 0; i <= nhash; i++) {
    hash[i] = (HASH_TABLE **)calloc((size_t)(nhash + 1), sizeof(HASH_TABLE *));
    for (j = 0; j <= nhash; j++) {
      hash[i][j] =
          (HASH_TABLE *)calloc((size_t)(nhash + 1), sizeof(HASH_TABLE));
    }
  }
  /**  pre make hash  **/
  pre_make_hash(nsse_q, nhash, hash, ssedat_q);

  /**  calloc HASH_TABLE 2  **/
  for (i = 0; i <= nhash; i++) {
    for (j = 0; j <= nhash; j++) {
      for (k = 0; k <= nhash; k++) {
        hash[i][j][k].isse =
            (int *)calloc((size_t)(hash[i][j][k].nlist + 1), sizeof(int));
        hash[i][j][k].ssetype =
            (int *)calloc((size_t)(hash[i][j][k].nlist + 1), sizeof(int));
        hash[i][j][k].para = (float **)calloc((size_t)(hash[i][j][k].nlist + 1),
                                              sizeof(float *));
        hash[i][j][k].perp = (float **)calloc((size_t)(hash[i][j][k].nlist + 1),
                                              sizeof(float *));
        for (ilist = 0; ilist <= hash[i][j][k].nlist; ilist++) {
          hash[i][j][k].para[ilist] = (float *)calloc(3, sizeof(float));
          hash[i][j][k].perp[ilist] = (float *)calloc(3, sizeof(float));
        }
      }
    }
  }
  /**  make hash table  **/
  make_hash_table(nsse_q, nhash, hash, ssedat_q);

  /**  calloc HASH_RESULT  **/
  maxalign = nsse_q * nsse_t * 2;
  hres = (HASH_RESULT *)calloc((size_t)(maxalign + 1), sizeof(HASH_RESULT));
  for (ialign = 0; ialign <= maxalign; ialign++) {
    hres[ialign].isse_q = 0;
    hres[ialign].isse_t = 0;
    hres[ialign].count = 0.0;
    hres[ialign].coord_mode = 0;
  }

  /*******************************/
  /**     search hash table     **/
  /*******************************/
  time_tmp = clock();
  search_hash(nsse_q, nsse_t, maxalign, nhash, hash, ssedat_t, hres);

  /***************************/
  /**    free hash table    **/
  /***************************/
  for (i = 0; i <= nhash; i++) {
    for (j = 0; j <= nhash; j++) {
      for (k = 0; k <= nhash; k++) {
        free(hash[i][j][k].isse);
        free(hash[i][j][k].ssetype);
        for (ilist = 0; ilist <= hash[i][j][k].nlist; ilist++) {
          free(hash[i][j][k].para[ilist]);
          free(hash[i][j][k].perp[ilist]);
        }
        free(hash[i][j][k].para);
        free(hash[i][j][k].perp);
      }
      free(hash[i][j]);
    }
    free(hash[i]);
  }
  free(hash);

  /*******************************/
  /**       calloc ALIGN        **/
  /*******************************/
  nalign = input.GHcutnum;

  align = (ALIGN *)calloc((size_t)(nalign + 1), sizeof(ALIGN));
  for (ialign = 0; ialign <= nalign; ialign++) {
    align[ialign].isse_q = 0;
    align[ialign].isse_t = 0;
    align[ialign].count = 0.0;
    align[ialign].coord_mode = 0;
    align[ialign].vec = (float *)calloc(3 + 1, sizeof(float));
    align[ialign].rot = (float **)calloc(3 + 1, sizeof(float *));
    for (i = 0; i <= 3; i++) {
      align[ialign].rot[i] = (float *)calloc(3 + 1, sizeof(float));
    }
    align[ialign].align_q = (int *)calloc((size_t)(naa_q + 1), sizeof(int));
    align[ialign].align_t = (int *)calloc((size_t)(naa_t + 1), sizeof(int));
    align[ialign].segmode = (int *)calloc((size_t)(naa_q + 1), sizeof(int));
    align[ialign].termini = (int *)calloc((size_t)(naa_q + 1), sizeof(int));
    align[ialign].local_score = (float *)calloc((size_t)(naa_q+1), sizeof(float));
    align[ialign].local_TM = (float *)calloc((size_t)(naa_q+1), sizeof(float));
    align[ialign].local_mTM = (float *)calloc((size_t)(naa_q+1), sizeof(float));
    align[ialign].local_aTM = (float *)calloc((size_t)(naa_q+1), sizeof(float));
    align[ialign].local_Dali = (float *)calloc((size_t)(naa_q+1), sizeof(float));
    align[ialign].local_dist = (float *)calloc((size_t)(naa_q+1), sizeof(float));
    align[ialign].naa_align = 0;
    align[ialign].rmsd = 0.0;
    align[ialign].TMscore_q = 0.0;
    align[ialign].TMscore_t = 0.0;
    align[ialign].TMscore_mod = 0.0;
    align[ialign].TMscore_ang = 0.0;
    align[ialign].C3score = 0.0;
    align[ialign].Daliscore = 0.0;
    align[ialign].DaliZ = 0.0;
    align[ialign].SPscore = 0.0;
    align[ialign].TMloc_mod = 0.0;
    align[ialign].score_cmp = 0.0;
    align[ialign].check = FALSE;
  }
  align_tmp.vec = (float *)calloc(3 + 1, sizeof(float));
  align_tmp.rot = (float **)calloc(3 + 1, sizeof(float *));
  for (i = 0; i <= 3; i++) {
    align_tmp.rot[i] = (float *)calloc(3 + 1, sizeof(float));
  }
  align_tmp.align_q = (int *)calloc((size_t)(naa_q + 1), sizeof(int));
  align_tmp.align_t = (int *)calloc((size_t)(naa_t + 1), sizeof(int));
  align_tmp.segmode = (int *)calloc((size_t)(naa_q + 1), sizeof(int));
  align_tmp.termini = (int *)calloc((size_t)(naa_q + 1), sizeof(int));
  align_tmp.local_score = (float *)calloc((size_t)(naa_q+1), sizeof(float));
  align_tmp.local_TM = (float *)calloc((size_t)(naa_q+1), sizeof(float));
  align_tmp.local_mTM = (float *)calloc((size_t)(naa_q+1), sizeof(float));
  align_tmp.local_aTM = (float *)calloc((size_t)(naa_q+1), sizeof(float));
  align_tmp.local_Dali = (float *)calloc((size_t)(naa_q+1), sizeof(float));
  align_tmp.local_dist = (float *)calloc((size_t)(naa_q+1), sizeof(float));

  /*********************************/
  /**       check GH result       **/
  /*********************************/
  check_GHresult(maxalign, hres, align, ssedat_q, ssedat_t);

  /**  computation time for GH search  **/
  time_end = clock();
  if (input.progress == ON) {
    fprintf(stderr, "\033[40C in %5.2f sec.\n",
            (double)(time_end - time_tmp) / CLOCKS_PER_SEC);
  }

  /*********************/
  /* chain break check */
  /*********************/
  chain_break_check(naa_q, resdat_q);
  chain_break_check(naa_t, resdat_t);

  /********************************/
  /**      calloc  matrix        **/
  /********************************/
  mat = (float ***)calloc((size_t)(Nite_loc + 1), sizeof(float **));
  for (i = 0; i <= Nite_loc; i++) {
    mat[i] = (float **)calloc((size_t)(naa_q + 2), sizeof(float *));
    for (iaa = 0; iaa <= naa_q + 1; iaa++) {
      mat[i][iaa] = (float *)calloc((size_t)(naa_t + 2), sizeof(float));
    }
  }
  score = (float **)calloc((size_t)(naa_q + 2), sizeof(float *));
  score_r = (float **)calloc((size_t)(naa_q + 2), sizeof(float *));
  for (iaa = 0; iaa <= naa_q + 1; iaa++) {
    score[iaa] = (float *)calloc((size_t)(naa_t + 2), sizeof(float));
    score_r[iaa] = (float *)calloc((size_t)(naa_t + 2), sizeof(float));
  }

  /**********************************/
  /**     first alignment step     **/
  /**********************************/
  time_tmp = clock();
  for (ialign = 0; ialign < nalign; ialign++) {

    if (input.progress == ON) {
      print_progress("Initial Alignment", ialign + 1, nalign);
    }
    if (align[ialign].check == FALSE) {
      continue;
    }
    /****   alignment   ****/
    align[ialign].check = ca_alignment(naa_q, naa_t, resdat_q, resdat_t,
                                       align[ialign], mat, score, score_r);
  }
  time_end = clock();
  if (input.progress == ON) {
    fprintf(stderr, "\033[40C in %5.2f sec.\n",
            (double)(time_end - time_tmp) / CLOCKS_PER_SEC);
  }

  /***************************/
  /**   score calculation   **/
  /***************************/
  score_calculation(naa_q, naa_t, nalign, align, resdat_q, resdat_t);

  /**************************/
  /**   check redundancy   **/
  /**************************/
  clustering_align(naa_q, naa_t, nalign, align, input.clst_red);

  /**  count non-redundant alignments  **/
  ntrue = 0;
  for (ialign = 0; ialign < nalign; ialign++) {
    if (align[ialign].check == TRUE) {
      ntrue++;
    }
  }

  /***********************************/
  /**     iterative refinement      **/
  /***********************************/
  time_tmp = clock();
  itrue = 0;
  for (ialign = 0; ialign < nalign; ialign++) {
    if (align[ialign].check == FALSE) {
      continue;
    }
    itrue++;
    if (input.progress == ON) {
      print_progress("Refine Alignment", itrue, ntrue);
    }
    irefine = 0;
    while (1) {
      irefine++;
      copy_align(naa_q, naa_t, &align_tmp, &(align[ialign]));
      /****   refine alignment   ****/
      superposition(naa_q, resdat_q, resdat_t, align[ialign], OFF);
      align[ialign].check = ca_alignment(naa_q, naa_t, resdat_q, resdat_t,
                                         align[ialign], mat, score, score_r);
      /****   check alignment   ****/
      if (align[ialign].check == FALSE) {
        copy_align(naa_q, naa_t, &(align[ialign]), &align_tmp);
        break;
      }
      if (check_align(naa_q, naa_t, &(align[ialign]), resdat_q, resdat_t) == FALSE) {
        copy_align(naa_q, naa_t, &(align[ialign]), &align_tmp);
        break;
      }
      if (irefine >= input.max_refine) {
        break;
      }
    }
  }
  time_end = clock();
  if (input.progress == ON) {
    fprintf(stderr, "\033[40C in %5.2f sec.\n\n",
            (double)(time_end - time_tmp) / CLOCKS_PER_SEC);
  }

  /*****************************/
  /**    score calculation    **/
  /*****************************/
  /**  Dali score & weight  **/
  for(ialign=0; ialign<nalign; ialign++){
    if(align[ialign].check == ON){
      align[ialign].DaliZ = Daliscore(naa_q, distmat_q, distmat_t, &(align[ialign]));
      for(iaa=1; iaa<=naa_q; iaa++){
	align[ialign].local_score[iaa] = align[ialign].local_Dali[iaa];
      }
      superposition(naa_q, resdat_q, resdat_t, align[ialign], ON);
    }
  }

  /**  other scores  **/
  score_calculation(naa_q, naa_t, nalign, align, resdat_q, resdat_t);
  calc_seq_identity(naa_q, naa_t, nalign, align, resdat_q, resdat_t);

  /*****************************/
  /**       clustering        **/
  /*****************************/
  clustering_align(naa_q, naa_t, nalign, align, input.clst_cutoff);

  /*****************************/
  /**     check alignment     **/
  /*****************************/
  for (ialign = 0; ialign <= input.nsub - 1; ialign++) {
    if (align[ialign].check == FALSE) {
      input.nsub = ialign;
      break;
    }
  }

  /**********************************************/
  /*  make align structures for oligomer output */
  /**********************************************/
  /* initialize align structure */
  align_chain = (ALIGN ***)calloc((size_t)(input.nsub + 1), sizeof(ALIGN **));
  for (i = 0; i <= input.nsub; i++) {
    align_chain[i] =
        (ALIGN **)calloc((size_t)(pdbdat_q.nchain + 1), sizeof(ALIGN *));
    for (ichain = 0; ichain <= pdbdat_q.nchain; ichain++) {
      align_chain[i][ichain] =
          (ALIGN *)calloc((size_t)(pdbdat_t.nchain + 1), sizeof(ALIGN));
    }
  }

  /* calc chain info */
  for (i = 0; i <= input.nsub - 1; i++) {
    score_calc_chain(naa_q, naa_t, align[i], align_chain[i], resdat_q, resdat_t,
		     pdbdat_q, pdbdat_t);
    align[i].C3score = C3score(pdbdat_q.nchain, pdbdat_t.nchain,
    			       align[i].naa_align, align_chain[i]);
  }

  /***********************************/
  /**      OUTPUT  structures       **/
  /***********************************/
  output(naa_q, naa_t, natom_q, natom_t, resdat_q, resdat_t, allatm_q,
	 allatm_t, align, align_chain, pdbdat_q, pdbdat_t);

  /********************************/
  /**           free             **/
  /********************************/
  /**  for resdat, ssedat  **/
  free(resdat_q);
  free(resdat_t);
  free(allatm_q);
  free(allatm_t);
  free(ssedat_q);
  free(ssedat_t);
  for(iaa=0; iaa<=naa_q; iaa++) free(distmat_q[iaa]);
  for(iaa=0; iaa<=naa_t; iaa++) free(distmat_t[iaa]);
  free(distmat_q);
  free(distmat_t);

  /**  for alignment  **/
  for (ialign = 0; ialign <= nalign; ialign++) {
    free(align[ialign].align_q);
    free(align[ialign].align_t);
    free(align[ialign].segmode);
    free(align[ialign].termini);
    free(align[ialign].local_score);
    free(align[ialign].local_TM);
    free(align[ialign].local_mTM);
    free(align[ialign].local_aTM);
    free(align[ialign].local_Dali);
    free(align[ialign].local_dist);
    free(align[ialign].vec);
    for (i = 0; i <= 3; i++) {
      free(align[ialign].rot[i]);
    }
    free(align[ialign].rot);
  }
  free(align);
  free(align_tmp.align_q);
  free(align_tmp.align_t);
  free(align_tmp.segmode);
  free(align_tmp.termini);
  free(align_tmp.local_score);
  free(align_tmp.local_TM);
  free(align_tmp.local_mTM);
  free(align_tmp.local_aTM);
  free(align_tmp.local_Dali);
  free(align_tmp.local_dist);
  free(align_tmp.vec);
  for (i = 0; i <= 3; i++) {
    free(align_tmp.rot[i]);
  }
  free(align_tmp.rot);
  /**  for hash result  **/
  free(hres);
  /**  for matrix  **/
  for (i = 0; i <= Nite_loc; i++) {
    for (iaa = 0; iaa <= naa_q + 1; iaa++)
      free(mat[i][iaa]);
    free(mat[i]);
  }
  free(mat);
  for (iaa = 0; iaa <= naa_q + 1; iaa++) {
    free(score[iaa]);
    free(score_r[iaa]);
  }
  free(score);
  free(score_r);

  /** for pdbdat **/
  for (ichain = 0; ichain <= pdbdat_q.nchain; ichain++) {
    free(pdbdat_q.chainID_org[ichain]);
  }
  free(pdbdat_q.chainID_org);
  free(pdbdat_q.chainID);
  for (ichain = 0; ichain <= pdbdat_t.nchain; ichain++) {
    free(pdbdat_t.chainID_org[ichain]);
  }
  free(pdbdat_t.chainID_org);
  free(pdbdat_t.chainID);

  free(pdbdat_q.naa);
  free(pdbdat_t.naa);

  /** for align_chain **/
  for (i = 0; i <= input.nsub; i++) {
    for (ichain = 0; ichain <= pdbdat_q.nchain; ichain++) {
      free(align_chain[i][ichain]);
    }
    free(align_chain[i]);
  }
  free(align_chain);

  return;
}



void copy_align(int naa_q, int naa_t, ALIGN *align1, ALIGN *align2) {

  int i, j;

  align1->isse_q = align2->isse_q;
  align1->isse_t = align2->isse_t;
  align1->count = align2->count;
  align1->coord_mode = align2->coord_mode;
  for (i = 0; i <= 3; i++) {
    align1->vec[i] = align2->vec[i];
  }
  for (i = 0; i <= 3; i++) {
    for (j = 0; j <= 3; j++) {
      align1->rot[i][j] = align2->rot[i][j];
    }
  }
  align1->naa_align = align2->naa_align;
  for (i = 0; i <= naa_q; i++) {
    align1->align_q[i] = align2->align_q[i];
    align1->segmode[i] = align2->segmode[i];
    align1->termini[i] = align2->termini[i];
    align1->local_score[i] = align2->local_score[i];
    align1->local_TM[i] = align2->local_TM[i];
    align1->local_mTM[i] = align2->local_mTM[i];
    align1->local_aTM[i] = align2->local_aTM[i];
    align1->local_Dali[i] = align2->local_Dali[i];
    align1->local_dist[i] = align2->local_dist[i];
  }
  for (i = 0; i <= naa_t; i++) {
    align1->align_t[i] = align2->align_t[i];
  }
  align1->rmsd = align2->rmsd;
  align1->TMscore_mod = align2->TMscore_mod;
  align1->TMscore_ang = align2->TMscore_ang;
  align1->TMscore_q = align2->TMscore_q;
  align1->TMscore_t = align2->TMscore_t;
  align1->Daliscore = align2->Daliscore;
  align1->DaliZ = align2->DaliZ;
  align1->SPscore = align2->SPscore;
  align1->TMloc_mod = align2->TMloc_mod;
  align1->score_cmp = align2->score_cmp;
  align1->check = align2->check;

  return;
}
