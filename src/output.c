#include "mican.h"
#include "main.h"

void output(int naa_q, int naa_t, int natom_q, int natom_t, RESDAT *resdat_q,
            RESDAT *resdat_t, ALLATM *allatm_q, ALLATM *allatm_t, ALIGN *align,
            ALIGN ***align_chain, PDBDAT pdbdat_q, PDBDAT pdbdat_t) {
  int i;
  FILE *fp;
  int ichain, jchain;

  char filename1[STRSIZE], filename2[STRSIZE];
  int naa1, nchain1, nchain2;
  char *chainID1, *chainID2;
  float TMscore1, TMscore2, TMloc_mod;
  float coverage1, coverage2;

  /*************************/
  /**    main output    ****/
  /*************************/
  /****  output protein data  ****/
  if (input.qtchange == OFF) {
    strcpy(filename1, input.file_t);
    strcpy(filename2, input.file_q);
    naa1 = naa_t;
    nchain1 = pdbdat_t.nchain;
    nchain2 = pdbdat_q.nchain;
  } else {
    strcpy(filename1, input.file_q);
    strcpy(filename2, input.file_t);
    naa1 = naa_q;
    nchain1 = pdbdat_q.nchain;
    nchain2 = pdbdat_t.nchain;
  }

  /*****************************/
  /**    brief description    **/
  /*****************************/
  if(input.silent == OFF){
    printf(" ===== Brief description of top %d alignments =====\n", input.nsub);
    printf(" Rank   sTMscore  TMscore  Dali.Z   SPscore  Length   RMSD   Seq.Id.\n");
    for (i = 0; i <= input.nsub - 1; i++) {
      if (align[i].check == OFF) {
	break;
      }
      if (input.qtchange == OFF) {
	TMscore1 = align[i].TMscore_t;
      } else {
	TMscore1 = align[i].TMscore_q;
      }
      printf(" %4d    %6.3f   %6.3f   %6.2f   %6.3f    %4d  %6.3f   %5.1f\n",
	     i + 1, align[i].TMscore_mod, TMscore1, align[i].DaliZ,
	     align[i].SPscore, align[i].naa_align,
	     align[i].rmsd, align[i].seqID);
    }
    printf(" (TMscore was normalized by size of Protein1. 'Dali.Z'=Dali"
	   " Zscore.\n  'Length'=number of aligned residues. 'Seq.Id.'=Sequence Identity.)\n");
    if (input.nsub_org != DEF_nsub && input.nsub_org > input.nsub) {
      fprintf(stderr,
	      "!! Warning: Mican detected only %d alignments, though\n",
	      input.nsub);
      fprintf(stderr,
	      "!! maximum number of solutions were set to %d.\n",
	      input.nsub_org);
      fprintf(stderr,
	      "!! For more alignments, please set a number of GH candidates\n");
      fprintf(stderr, "!! with -g option (E.g. '-g 100', '-g 500', or more.)\n");
    }
    printf("\n");

  }

  if (align[input.isub_out - 1].check != TRUE) {
    printf(" Warning: %d-th solution can not be detected.\n\n",
	   input.isub_out);
    exit(0);
  }

  /********************************/
  /**    detailed description    **/
  /********************************/
  if(input.silent == OFF){
    printf(" ===== Detailed description of alignment %d =====\n", input.isub_out);
    
    if (input.qtchange == OFF) {
      TMscore1 = align[input.isub_out - 1].TMscore_t;
      TMscore2 = align[input.isub_out - 1].TMscore_q;
      coverage1 = align[input.isub_out - 1].cover_t;
      coverage2 = align[input.isub_out - 1].cover_q;
    } else {
      TMscore1 = align[input.isub_out - 1].TMscore_q;
      TMscore2 = align[input.isub_out - 1].TMscore_t;
      coverage1 = align[input.isub_out - 1].cover_q;
      coverage2 = align[input.isub_out - 1].cover_t;
    }
    TMloc_mod = align[input.isub_out - 1].TMloc_mod;

    printf(" TM-score=%5.3f, Coverage=%5.1f%% (if normalized by size of Protein1)\n",
	    TMscore1, coverage1);
    printf(" TM-score=%5.3f, Coverage=%5.1f%% (if normalized by size of Protein2)\n",
	    TMscore2, coverage2);
    printf(" sTMscore(Nali)=%5.3f (sTMscore normalized by aligned length)\n",
	    TMloc_mod);
    printf(" C3score=%5.3f (Chain-to-Chain Correspondence score)\n",
      align[input.isub_out-1].C3score);
    printf(" Dali-score=%7.2f, Dali Zscore=%5.2f\n",
	    align[input.isub_out-1].Daliscore, align[input.isub_out-1].DaliZ);
    
    /****  chain results  ****/
    printf(" ----- Results for each chain -----\n");
    for (jchain = 1; jchain <= nchain1; jchain++) {
      if (input.qtchange == OFF) {
	      coverage1 = align_chain[input.isub_out - 1][0][jchain].cover_t;
      } else {
	      coverage1 = align_chain[input.isub_out - 1][jchain][0].cover_q;
      }
      if (coverage1 > 0.0) {
	      if (input.qtchange == OFF) {
	        chainID1 = pdbdat_t.chainID_org[jchain];
	        naa1 = pdbdat_t.naa[jchain];
	      } else {
	        chainID1 = pdbdat_q.chainID_org[jchain];
	        naa1 = pdbdat_q.naa[jchain];
	      }
	      printf(" [P1:%1s (size %4d)] Coverage=%5.1f%% (",
	        chainID1, naa1, coverage1);
	
	      for (ichain = 1; ichain <= nchain2; ichain++) {
	        if (input.qtchange == OFF) {
	          coverage2 = align_chain[input.isub_out - 1][ichain][jchain].cover_t;
	        } else {
	          coverage2 = align_chain[input.isub_out - 1][jchain][ichain].cover_q;
	        }
    	    if (coverage2 > 0.0) {
	          if (input.qtchange == OFF) {
	            chainID2 = pdbdat_q.chainID_org[ichain];
	          } else {
	            chainID2 = pdbdat_t.chainID_org[ichain];
	          }
	          printf("P2:%1s %5.1f%% ", chainID2, coverage2);
	        }
	      }
	      printf(")\n");
      }
    }
  }

  /***************************************/
  /**    output translation matrix      **/
  /***************************************/
  if(input.qtchange == ON){
    inverse_mat(align[input.isub_out-1].rot, align[input.isub_out-1].vec);
  }
  if(input.silent == OFF){
    printmat(stdout, align, OFF);
  }    

  /***************************************/
  /**      output result file name       **/
  /***************************************/
  if(input.silent == OFF){
    if (strcmp(input.aliout, "OFF") != 0) {
      printf(" alignment file          = %s\n", input.aliout);
    }
    if (strcmp(input.pdbout, "OFF") != 0) {
      printf(" superposition pdb file  = %s\n", input.pdbout);
    }
    if (strcmp(input.matout, "OFF") != 0) {
      printf(" translation matrix file = %s\n", input.matout);
    }
  }

  /***************************************/
  /**     print alignment on STDOUT     **/
  /***************************************/
  if(input.printalign == ON && input.silent == OFF){
    printali(stdout, naa_q, naa_t, resdat_q, resdat_t, align, ON);
  }

  /***************************************/
  /**       output alignment file       **/
  /***************************************/
  if (strncmp(input.aliout, "OFF", 3) != 0) {
    fp = fopen(input.aliout, "w");
    printali(fp, naa_q, naa_t, resdat_q, resdat_t, align, OFF);
    fclose(fp);
  }

  /****************************************/
  /**    output superposition pdbfile    **/
  /****************************************/
  if (strncmp(input.pdbout, "OFF", 3) != 0) {
    fp = fopen(input.pdbout, "w");
    printsup(fp, naa_q, naa_t, natom_q, natom_t, resdat_q, resdat_t, allatm_q,
             allatm_t, align);
    fclose(fp);
  }

  /****************************************/
  /**   output translation matrix file   **/
  /****************************************/
  if (strncmp(input.matout, "OFF", 3) != 0) {
    fp = fopen(input.matout, "w");
    printmat(fp, align, ON);
    fclose(fp);
  }

  return;
}
