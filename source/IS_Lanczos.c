/**********************************************************************
  IS_Lanczos.c:

     IS_Lanczos.c is a subroutine to calculate the inverse of 
     an overlap matrix using a recursion method in O(N) operation.
     (based on PRB 64, 195110 (2001) and a modified recurrence relation)

  Log of IS_Lanczos.c:

     22/Nov/2001  Released by T.Ozaki
     06/Jul/2005  Modified by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "openmx_common.h"
#include "lapack_prototypes.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

#include "tran_prototypes.h"



static void Construct_Solp(int Mc_AN, double ****OLP, double **Solp, int *MP);
static void Lanczos(int Mc_AN, double **Solp, int *MP, int rlmax_IS);
static void S_Re_Green(int Mc_AN, int T_switch, double ****IOLP, int rlmax_IS);
static void Inverse_Mat(int n, double **a);

static int *Msize,*MP;
static int *RCN_Rank;

static double **Solp;
static double ***U;
static double ***al;
static double ***be;
static double ***ibe;
static double **Cpre;
static double **Ccnt;
static double **Rcnt;
static double **SqBe,*ko,*iko;

static double *be2,*ter;
static double **inv;
static double **tp,**tp1;
static double **tp2,**tp3;
static double **G00,***G0;

static INTEGER *ipiv;
static double *A;
static double *work;

static int Msize_max;






void IS_Lanczos(double ****OLP, double ****IOLP, int rlmax_IS)
{
  static int firsttime=1;
  static int Mc_AN,Gc_AN,n2,Anum,Gi;
  static int wanA,NUM;
  static int size_U,size_al,size_Solp,size_Cpre;
  static int i,j,k,l,m,tno0,tno1;
  static int Hwan,Cwan,Gh_AN;

  static int numprocs,myid,ID,tag=999;

  MPI_Status stat;
  MPI_Request request; 

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  printf("A1\n");

  /****************************************************
     allocation of arrays:
  ****************************************************/

  /**********************
    find the matrix size
  ***********************/

  Msize = (int*)malloc(sizeof(int)*(Matomnum+1));
  MP = (int*)malloc(sizeof(int)*List_YOUSO[2]);

  Msize_max = 0;
  
  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];

    if (Mc_AN==0){
      Gc_AN = 0;
      FNAN[0] = 1;
      SNAN[0] = 0;
      n2 = 1;
      Msize[Mc_AN] = 1;
    }
    else{
      Anum = 1;
      for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
	Gi = natn[Gc_AN][i];
	wanA = WhatSpecies[Gi];
	Anum += Spe_Total_CNO[wanA];
      }
      NUM = Anum - 1;
      Msize[Mc_AN] = NUM;
    }
 
    if (Msize_max<Msize[Mc_AN]) Msize_max = Msize[Mc_AN] + 4;
  }

  /* RCN_Rank */

  RCN_Rank = (int*)malloc(sizeof(int)*(List_YOUSO[9]+3)); 

  /* Solp */

  size_Solp = Msize_max*Msize_max;

  Solp = (double**)malloc(sizeof(double*)*Msize_max);
  for (i=0; i<Msize_max; i++){
    Solp[i] = (double*)malloc(sizeof(double)*Msize_max);
  }

  /* U */
  
  size_U = (List_YOUSO[9]+1)*List_YOUSO[7]*Msize_max;

  U = (double***)malloc(sizeof(double**)*(List_YOUSO[9]+1));
  for (i=0; i<(List_YOUSO[9]+1); i++){
    U[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
    for (j=0; j<List_YOUSO[7]; j++){
      U[i][j] = (double*)malloc(sizeof(double)*Msize_max);
    }
  }

  /* al */

  size_al = List_YOUSO[9]*List_YOUSO[7]*List_YOUSO[7]; 

  al = (double***)malloc(sizeof(double**)*List_YOUSO[9]);
  for (i=0; i<List_YOUSO[9]; i++){
    al[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
    for (j=0; j<List_YOUSO[7]; j++){
      al[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    }
  }

  /* be */

  be = (double***)malloc(sizeof(double**)*List_YOUSO[9]);
  for (i=0; i<List_YOUSO[9]; i++){
    be[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
    for (j=0; j<List_YOUSO[7]; j++){
      be[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    }
  }

  /* ibe */

  ibe = (double***)malloc(sizeof(double**)*List_YOUSO[9]);
  for (i=0; i<List_YOUSO[9]; i++){
    ibe[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
    for (j=0; j<List_YOUSO[7]; j++){
      ibe[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    }
  }

  /* Cpre */

  size_Cpre = List_YOUSO[7]*Msize_max;

  Cpre = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    Cpre[i] = (double*)malloc(sizeof(double)*Msize_max);
  }

  /* Ccnt */

  Ccnt = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    Ccnt[i] = (double*)malloc(sizeof(double)*Msize_max);
  }

  /* Rcnt */

  Rcnt = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    Rcnt[i] = (double*)malloc(sizeof(double)*Msize_max);
  }

  /* SqBe */

  SqBe = (double**)malloc(sizeof(double*)*(List_YOUSO[7]+2));
  for (i=0; i<(List_YOUSO[7]+2); i++){
    SqBe[i] = (double*)malloc(sizeof(double)*(List_YOUSO[7]+2));
  }

  /* ko and iko */

  ko = (double*)malloc(sizeof(double)*(List_YOUSO[7]+2));
  iko = (double*)malloc(sizeof(double)*(List_YOUSO[7]+2));

  /* be2 */

  be2 = (double*)malloc(sizeof(double)*(List_YOUSO[7]+2));

  /* ter */

  ter = (double*)malloc(sizeof(double)*List_YOUSO[7]);

  /* inv */

  inv = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    inv[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
  } 

  /* tp */

  tp = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    tp[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
  } 

  /* tp1 */

  tp1 = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    tp1[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
  } 

  /* tp2 */

  tp2 = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    tp2[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
  } 

  /* tp3 */

  tp3 = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    tp3[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
  } 

  /* G00 */

  G00 = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    G00[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
  } 

  /* G0 */

  G0 = (double***)malloc(sizeof(double**)*List_YOUSO[9]);
  for (i=0; i<List_YOUSO[9]; i++){
    G0[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
    for (j=0; j<List_YOUSO[7]; j++){
      G0[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    }
  } 

  ipiv = (INTEGER*)malloc(sizeof(INTEGER)*List_YOUSO[7]);
  A = (double*)malloc(sizeof(double)*List_YOUSO[7]*List_YOUSO[7]);
  work = (double*)malloc(sizeof(double)*List_YOUSO[7]);

  /* PrintMemory */
  if (firsttime) {
    PrintMemory("IS_Lanczos: U",sizeof(double)*size_U,    NULL);
    PrintMemory("IS_Lanczos: al",sizeof(double)*size_al,  NULL);
    PrintMemory("IS_Lanczos: be",sizeof(double)*size_al,  NULL);
    PrintMemory("IS_Lanczos: ibe",sizeof(double)*size_al, NULL);
    PrintMemory("IS_Lanczos: Cpre",sizeof(double)*size_Cpre, NULL);
    PrintMemory("IS_Lanczos: Ccnt",sizeof(double)*size_Cpre, NULL);
    PrintMemory("IS_Lanczos: Rcnt",sizeof(double)*size_Cpre, NULL);
    firsttime=0;
  }

  /****************************************************
     start the main part....
  ****************************************************/

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];

    /***********************************************
     MP array

     Note:
         MP indicates the starting position of
              atom i in arraies H and S
    ***********************************************/
    
    Anum = 1;
    for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
      MP[i] = Anum;
      Gi = natn[Gc_AN][i];
      wanA = WhatSpecies[Gi];
      Anum += Spe_Total_CNO[wanA];
    }

    /****************************************************
                      construct Solp
    ****************************************************/

    Construct_Solp( Mc_AN, OLP, Solp, MP );

    /****************************************************
               block Lanczos transformation
    ****************************************************/

    Lanczos( Mc_AN, Solp, MP, rlmax_IS );

    /****************************************************
      evaluation of real resolvent (Green's functions)
    ****************************************************/

    S_Re_Green( Mc_AN, T_switch, IOLP, rlmax_IS );

  }

  /****************************************************
     freeing of arrays:
  ****************************************************/

  /* Msize */

  free(Msize);

  /* MP */

  free(MP);

  /* RCN_Rank */

  free(RCN_Rank);

  /* Solp */

  for (i=0; i<Msize_max; i++){
    free(Solp[i]);
  }
  free(Solp);

  /* U */
  
  for (i=0; i<(List_YOUSO[9]+1); i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(U[i][j]);
    }
    free(U[i]);
  }
  free(U);

  /* al */

  for (i=0; i<List_YOUSO[9]; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(al[i][j]);
    }
    free(al[i]);
  }
  free(al);

  /* be */

  for (i=0; i<List_YOUSO[9]; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(be[i][j]);
    }
    free(be[i]);
  }
  free(be);

  /* ibe */

  for (i=0; i<List_YOUSO[9]; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(ibe[i][j]);
    }
    free(ibe[i]);
  }
  free(ibe);

  /* Cpre */

  for (i=0; i<List_YOUSO[7]; i++){
    free(Cpre[i]);
  }
  free(Cpre);

  /* Ccnt */

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ccnt[i]);
  }
  free(Ccnt);

  /* Rcnt */

  for (i=0; i<List_YOUSO[7]; i++){
    free(Rcnt[i]);
  }
  free(Rcnt);

  /* SqBe */

  for (i=0; i<(List_YOUSO[7]+2); i++){
    free(SqBe[i]);
  }
  free(SqBe);

  /* ko and iko */

  free(ko);
  free(iko);

  /* be2 */

  free(be2);

  /* ter */

  free(ter);

  /* inv */

  for (i=0; i<List_YOUSO[7]; i++){
    free(inv[i]);
  } 
  free(inv);

  /* tp */

  for (i=0; i<List_YOUSO[7]; i++){
    free(tp[i]);
  } 
  free(tp);

  /* tp1 */

  for (i=0; i<List_YOUSO[7]; i++){
    free(tp1[i]);
  } 
  free(tp1);

  /* tp2 */

  for (i=0; i<List_YOUSO[7]; i++){
    free(tp2[i]);
  } 
  free(tp2);

  /* tp3 */

  for (i=0; i<List_YOUSO[7]; i++){
    free(tp3[i]);
  } 
  free(tp3);

  /* G00 */

  for (i=0; i<List_YOUSO[7]; i++){
    free(G00[i]);
  } 
  free(G00);

  /* G0 */

  for (i=0; i<List_YOUSO[9]; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(G0[i][j]);
    }
    free(G0[i]);
  } 
  free(G0);

  free(ipiv);
  free(A);
  free(work);

}








void Construct_Solp(int Mc_AN, double ****OLP, double **Solp, int *MP)
{
  static int Gc_AN,i,j,m,n,ig,jg,ian;
  static int Anum,Bnum,ih,kl,kg,jan;

  Gc_AN = M2G[Mc_AN];

  for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){

    ig = natn[Gc_AN][i];
    ian = Spe_Total_CNO[WhatSpecies[ig]];
    Anum = MP[i];
    ih = S_G2M[ig];

    for (j=0; j<=(FNAN[Gc_AN]+SNAN[Gc_AN]); j++){

      kl = RMI1[Mc_AN][i][j];
      jg = natn[Gc_AN][j];
      jan = Spe_Total_CNO[WhatSpecies[jg]];
      Bnum = MP[j];

      if (0<=kl){
	for (m=0; m<ian; m++){
	  for (n=0; n<jan; n++){
	    Solp[Anum+m][Bnum+n] = OLP[ih][kl][m][n];
	  }
	}
      }

      else{
	for (m=0; m<ian; m++){
	  for (n=0; n<jan; n++){
	    Solp[Anum+m][Bnum+n] = 0.0;
	  }
	}
      }
    }
  }
}








void Lanczos(int Mc_AN, double **Solp, int *MP, int rlmax_IS)
{
  static int firsttime=1;
  static int i,j,k,l,m,n,n2,ie,i1,j1,po,ig,jg,Rni,Rnj,Rn;
  static int kl,l1,l2,l3,ian,jan,count,Gc_AN,rl0;
  static int m1,m2,m3,ZeroNum,Rloop_END,Rloop_Term;
  static int fan,san,can,rl,wan,ct_on;
  static double sum,sum1,sum2,dum,dum1,sumx,sumy,sumz,xn,xa;

  Gc_AN = M2G[Mc_AN];

  fan = FNAN[Gc_AN];
  san = SNAN[Gc_AN];
  can = fan + san;
  wan = WhatSpecies[Gc_AN];
  ct_on = Spe_Total_CNO[wan];

  for (rl=0; rl<List_YOUSO[9]; rl++){
    for (i=0; i<ct_on; i++){
      for (j=0; j<ct_on; j++){
	al[rl][i][j] = 0.0;
	be[rl][i][j] = 0.0;
      }
    }
  }

  for (i=0; i<List_YOUSO[7]; i++){
    for (j=0; j<Msize_max; j++){
      Cpre[i][j] = 0.0;
      Ccnt[i][j] = 0.0;
      Rcnt[i][j] = 0.0;
    }
  }

  for (rl=0; rl<(List_YOUSO[9]+1); rl++){
    for (i=0; i<List_YOUSO[7]; i++){
      for (j=0; j<Msize_max; j++){
        U[rl][i][j] = 0.0;
      }
    }
  }

  for (i=0; i<ct_on; i++){
    Ccnt[i][i+1] = 1.0;
    U[0][i][i+1] = 1.0;
  }

  /****************************************************
       recursion in the block Lanczos algorithm
  ****************************************************/

  Rloop_END = 0;
  Rloop_Term = 0;
  RNUM[Gc_AN] = rlmax_IS;
  RCN_Rank[0] = ct_on;

  rl = 0;

  do {

    /****************************************************
                             H|Wn}
    ****************************************************/

    for (n=0; n<RCN_Rank[rl]; n++){
      for (i=1; i<=Msize[Mc_AN]; i++){
	sum = 0.0;
	for (j=1; j<=Msize[Mc_AN]; j++){
	  sum += Solp[i][j]*Ccnt[n][j];
	}
	Rcnt[n][i] = sum;
      } 
    }

    /****************************************************
                       Alpha_n = {Wn|H|Wn}
    ****************************************************/

    for (m=0; m<RCN_Rank[rl]; m++){
      for (n=0; n<RCN_Rank[rl]; n++){
        sum = 0.0; 
        for (i=1; i<=Msize[Mc_AN]; i++){
          sum += Ccnt[m][i]*Rcnt[n][i];
	}
        al[rl][m][n] = sum;
      }
    }


    if (rl!=rlmax_IS){

      /****************************************************
                              |Rcnt}
      ****************************************************/

      for (n=0; n<RCN_Rank[rl]; n++){
        for (i=1; i<=Msize[Mc_AN]; i++){

          sum = 0.0;  

	  if (1<=rl){ 
	   
	    /* note: transpotision */ 
	    for (k=0; k<RCN_Rank[rl-1]; k++){
	      sum += Cpre[k][i]*be[rl][n][k];
	    }
	  } 

	  for (k=0; k<RCN_Rank[rl]; k++){
	    sum += Ccnt[k][i]*al[rl][k][n];  
	  }

          Rcnt[n][i] -= sum; 
        }
      }         

      /****************************************************
           reorthogonalization by a modified block
                    Gram-Schmidt method
      ****************************************************/

      /* |Rcnt) := (I - \sum_{rl0} |U_rl0)(U_rl0|)|Rcnt) */

      for (rl0=0; rl0<=rl; rl0++){

        /* (U_rl0|Rcnt) */

	for (m=0; m<RCN_Rank[rl0]; m++){
	  for (n=0; n<RCN_Rank[rl]; n++){

	    sum = 0.0;
	    for (i=1; i<=Msize[Mc_AN]; i++){
	      sum += U[rl0][m][i]*Rcnt[n][i];
	    }

	    SqBe[m][n] = sum;
	  }
	}

        /* |Rcnt) := |Rcnt) - |U_rl0)(U_rl0|Rcnt) */

        for (n=0; n<RCN_Rank[rl]; n++){
          for (i=1; i<=Msize[Mc_AN]; i++){

            sum = 0.0;
	    for (k=0; k<RCN_Rank[rl0]; k++){
	      sum += U[rl0][k][i]*SqBe[k][n];  
	    }

            Rcnt[n][i] -= sum;
	  }
	}


      }



      /*
      for (n=0; n<ct_on; n++){
        for (i=1; i<=Msize[Mc_AN]; i++){
          printf("B rl=%2d n=%2d i=%2d Rcnt=%15.12f\n",rl,n,i,Rcnt[n][i]);
        }
      }
      */




      /****************************************************
                              Beta_(n+1)
      ****************************************************/

      for (m=0; m<RCN_Rank[rl]; m++){
	for (n=0; n<RCN_Rank[rl]; n++){
	  sum = 0.0; 
	  for (i=1; i<=Msize[Mc_AN]; i++){
	    sum += Rcnt[m][i]*Rcnt[n][i];
	  }
	  SqBe[m+1][n+1] = sum;
	}
      }

      /****************************************************
                     diagonalization of SqBe
      ****************************************************/

      printf("B7 Mc_AN=%2d rl=%2d\n",Mc_AN,rl);
      for (m=0; m<RCN_Rank[rl]; m++){
	for (n=0; n<RCN_Rank[rl]; n++){
	  printf("%18.15f ",SqBe[m+1][n+1]);
	}
        printf("\n");
      }


      if (RCN_Rank[rl]==1){
	ko[1] = SqBe[1][1];
	SqBe[1][1] = 1.0;
      }
      else{
	Eigen_lapack(SqBe,ko,RCN_Rank[rl],RCN_Rank[rl]);
      }

      /*
      for (m=0; m<RCN_Rank[rl]; m++){
        printf("A m=%2d ko=%18.15f\n",m,ko[m+1]);
      }
      */


      ZeroNum = ct_on - RCN_Rank[rl];

      for (i=1; i<=RCN_Rank[rl]; i++){

	/*
	if (1.0e-10<ko[i] && Rloop_Term==0){
	*/

	if (1.0e-10<ko[i]){
	  ko[i]  = sqrt(ko[i]);
	  iko[i] = 1.0/ko[i];
	}
	else{
	  ZeroNum++;
	  ko[i]  = 0.0;
	  iko[i] = 0.0;
	}
      }

      RCN_Rank[rl+1] = ct_on - ZeroNum;

      /*
      for (m=0; m<ct_on; m++){
        printf("B m=%2d ko=%18.15f\n",m,ko[m+1]);
      }
      */

  printf("B8 Mc_AN=%2d rl=%2d ZeroNum=%2d\n",Mc_AN,rl,ZeroNum);

      if (ZeroNum<ct_on){ 

        /* be = \lamda * U^t */

	if (ZeroNum!=0) Rloop_Term = 1;
	for (m=0; m<RCN_Rank[rl+1]; m++){
	  for (n=0; n<RCN_Rank[rl]; n++){
             /* note: transpotision */
 	     be[rl+1][m][n] =  ko[m+1+ZeroNum]*SqBe[n+1][m+1+ZeroNum];
	    ibe[rl+1][n][m] = iko[m+1+ZeroNum]*SqBe[n+1][m+1+ZeroNum];
	  }
	}

  printf("B9 Mc_AN=%2d rl=%2d ZeroNum=%2d\n",Mc_AN,rl,ZeroNum);

	/****************************************************
                              |Wn+1>
	****************************************************/


	for (n=0; n<RCN_Rank[rl]; n++){
	  for (i=1; i<=Msize[Mc_AN]; i++){
	    Cpre[n][i] = Ccnt[n][i];
	  }
	}

	for (n=0; n<RCN_Rank[rl+1]; n++){
	  for (i=1; i<=Msize[Mc_AN]; i++){

	    sum = 0.0;
	    for (k=0; k<RCN_Rank[rl]; k++){
	      sum += Rcnt[k][i]*ibe[rl+1][k][n];
	    }

	    Ccnt[n][i]    = sum;
	    U[rl+1][n][i] = sum;   
	  }
	}

  printf("B10 Mc_AN=%2d rl=%2d ZeroNum=%2d\n",Mc_AN,rl,ZeroNum);

      }
      else{
	RNUM[Gc_AN] = rl;
	Rloop_END = 1;
      }
    } /* end of if (rl!=rlmax_IS) */

    rl++;

  } while (rl<=rlmax_IS && Rloop_END==0);

}

















void S_Re_Green( int Mc_AN, int T_switch, double ****IOLP, int rlmax_IS )
{
  /****************************************************
    T_switch is a parameter which selects a terminator
    in the Green's function.

    T_switch = 1 is no terminator.
    T_switch = 2 is a square root terminator.
  ****************************************************/

  static int firsttime=1;
  static int i,j,k,n,Gc_AN,wan,tno,tno1,rl;
  static int itnum,dtnum,rl_loop,Anum,num0;
  static int rlm,rlm2,ig,ino1,p,q,m;
  static double d,sum,sum3,dum,dum1,dum2,Av_al;

  Gc_AN = M2G[Mc_AN];

  /****************************************************
                      start calc.
  ****************************************************/

  itnum = Av_num;
  dtnum = Av_num;

  wan = WhatSpecies[Gc_AN];
  tno = Spe_Total_CNO[wan];
  tno1 = tno - 1;




  for (rl=0; rl<=RNUM[Gc_AN]; rl++){   
    printf("al rl=%2d\n",rl);
    for (i=0; i<RCN_Rank[rl]; i++){
      for (j=0; j<RCN_Rank[rl]; j++){
        printf("%15.12f ",al[rl][i][j]);
      }
      printf("\n");
    }
  }

  for (rl=1; rl<=RNUM[Gc_AN]; rl++){   
    printf("be rl=%2d\n",rl);
    for (i=0; i<RCN_Rank[rl]; i++){
      for (j=0; j<RCN_Rank[rl-1]; j++){
        printf("%15.12f ",be[rl][i][j]);
      }
      printf("\n");
    }
  }

  /* initialize G0 */

  for (rl=0; rl<=RNUM[Gc_AN]; rl++){
    for (i=0; i<=tno1; i++){
      for (j=0; j<=tno1; j++){
	G0[rl][i][j] = 0.0;
      }
    }     
  }

  /****************************************************
                     without any terminator
  ****************************************************/

  if (T_switch==1){
    for (i=0; i<=tno1; i++){
      ter[i] = 0.0;
    }
  }

  /****************************************************
                   square root terminator
  ****************************************************/

  else{

    num0 = 0;
    for (i=0; i<=tno1; i++) be2[i] = 0.0;
    for (rl=RNUM[Gc_AN]; (RNUM[Gc_AN]-itnum+1)<=rl; rl--){   
      num0 += RCN_Rank[rl-1];
      for (i=0; i<RCN_Rank[rl-1]; i++){
	sum = 0.0;
	for (k=0; k<RCN_Rank[rl]; k++){
          /* note: transposition  */ 
	  sum += be[rl][k][i]*be[rl][k][i];
	}
	be2[i] += sum;
      }
    }

    /****************************************************
                           averaging
    ****************************************************/

    sum = 0.0;
    for (i=0; i<=tno1; i++){
      sum += be2[i];
    }
    sum = sum/(double)num0;
    for (i=0; i<=tno1; i++){
      be2[i] = sum;
    }

    sum = 0.0;
    for (rl=RNUM[Gc_AN]; (RNUM[Gc_AN]-itnum+1)<=rl; rl--){
      for (i=0; i<RCN_Rank[rl]; i++){
	sum += al[rl][i][i];
      }
    }

    Av_al = sum/(double)num0;

    for (i=0; i<=tno1; i++){
      d = Av_al*Av_al - 4.0*be2[i];
      if (0.0<=d){
	ter[i] = (Av_al-sqrt(d))/(2.0*be2[i]);
      }
      else{
	if (i==1){
	  printf("B Gc_AN Av_al be %i %10.9f %10.9f\n",
		  Gc_AN,Av_al,sqrt(be2[i]));
	}
	ter[i] = 0.0;
      }
    }
  }

  /****************************************************
            calculation of multiple inverse
  ****************************************************/

  for (i=0; i<=tno1; i++){
    for (j=0; j<=tno1; j++){
      if (i==j)
	inv[i][i] = ter[i];
      else{
	inv[i][j] = 0.0;
      }
    }
  }

  if (T_switch==1) rl_loop = RNUM[Gc_AN];
  else             rl_loop = RNUM[Gc_AN] - 1;

  for (rl=rl_loop; 0<=rl; rl--){

    if (rl==RNUM[Gc_AN]){

      /****************************************************
                     in case of no terminator
      ****************************************************/

      for (i=0; i<RCN_Rank[rl]; i++){
	for (j=0; j<RCN_Rank[rl]; j++){
	  inv[i][j] = al[rl][i][j];
	}
      }
    }
    else{

      /****************************************************
                           Beta^T * []^-1
      ****************************************************/

      for (i=0; i<RCN_Rank[rl]; i++){
	for (j=0; j<RCN_Rank[rl+1]; j++){

	  sum3 = 0.0;
	  for (k=0; k<RCN_Rank[rl+1]; k++){

            /* note: transpotision */

	    sum3 += be[rl+1][k][i]*inv[k][j]; 
	  }
	  tp1[i][j] = sum3;
	}
      }

      /****************************************************
                         Beta^T * []^-1 *Beta
      ****************************************************/

      for (i=0; i<RCN_Rank[rl]; i++){
	for (j=0; j<RCN_Rank[rl]; j++){
	  sum3 = 0.0;
	  for (k=0; k<RCN_Rank[rl+1]; k++){
	    dum = tp1[i][k]*be[rl+1][k][j];
	    sum3 = sum3 + dum;
	  }
	  tp2[i][j] = sum3;
	}
      }

      /****************************************************
                 Alpha - Beta * []^-1 * Beta
      ****************************************************/

      for (i=0; i<RCN_Rank[rl]; i++){
	for (j=0; j<RCN_Rank[rl]; j++){
	  inv[i][j] = al[rl][i][j] - tp2[i][j];
	}
      }

    } /* else */

    /****************************************************
            [ Alpha - Beta * []^-1 * Beta ]^-1
    ****************************************************/

  printf("C1 Mc_AN=%2d Gc_AN=%2d RNUM=%2d\n",Mc_AN,Gc_AN,RNUM[Gc_AN]);

    for (i=0; i<=tno1; i++){
      for (j=0; j<=tno1; j++){
        printf("%15.12f ",inv[i][j]);
      }
      printf("\n");
    }

    Inverse_Mat(RCN_Rank[rl]-1, inv);

  } /* rl */


  printf("G00\n");

  for (i=0; i<=tno1; i++){
    for (j=0; j<=tno1; j++){
      G00[i][j] = inv[i][j];

      printf("%18.15f ",G00[i][j]);

    }
    printf("\n");
  }

  /****************************************************
            Start of RecurG(Gc_AN,0.0,G00,G0)
  ****************************************************/

  for (i=0; i<=tno1; i++){
    for (j=0; j<=tno1; j++){
      G0[0][i][j] = G00[i][j];
    }
  }

  for (rl=1; rl<=RNUM[Gc_AN]; rl++){
    for (i=0; i<=tno1; i++){
      for (j=0; j<=tno1; j++){
        G0[rl][i][j] = 0.0;
      }
    }
  }

  /****************************************************
    calculate G01
  ****************************************************/

  /* I - G00*A0 */

  for (i=0; i<=tno1; i++){
    for (j=0; j<=tno1; j++){
      tp1[i][j] = al[0][i][j];
    }
  }

  for (i=0; i<=tno1; i++){
    for (j=0; j<=tno1; j++){
      sum = 0.0;
      for (k=0; k<=tno1; k++){
	sum += G0[0][i][k]*tp1[k][j];
      }
      tp2[i][j] = sum;
    }
  }

  for (i=0; i<=tno1; i++){
    for (j=0; j<=tno1; j++){
      if (i==j)	tp3[i][j] = -tp2[i][j] + 1.0;
      else	tp3[i][j] = -tp2[i][j];
    }
  }

  /* G01 = [ I - G00*A0 ] * B1^{-1} */

  printf("G01\n");

  for (i=0; i<=tno1; i++){
    for (j=0; j<RCN_Rank[1]; j++){
      sum = 0.0;
      for (k=0; k<=tno1; k++){
	sum += tp3[i][k]*ibe[1][k][j]; 
      }
      G0[1][i][j] = sum;

      printf("%18.15f ",G0[1][i][j]);

    }

    printf("\n");

  }

  /****************************************************
        calculate a series of ratio functions 
        by the following recurrence relation

    R_n = - Bn^{t} * [ An + R_{n+1}*B_{n+1} ]^{-1}
  ****************************************************/

  printf("C2 Mc_AN=%2d Gc_AN=%2d RNUM=%2d\n",Mc_AN,Gc_AN,RNUM[Gc_AN]);

  for (rl=RNUM[Gc_AN]; 2<=rl; rl--){

    /* set R_{RNUM[Gc_AN]+1} */
    /* This part can be modified by SRT. */

    if (rl==RNUM[Gc_AN]){
      for (i=0; i<=tno1; i++){
	for (j=0; j<=tno1; j++){
          G0[rl+1][i][j] = 0.0;
	}
      }     
    }

    /****************************************************
        R_n = - Bn^{t} * [ An + R_{n+1}*B_{n+1} ]^{-1}
    ****************************************************/

    /* R_{n+1}*B_{n+1} */    
    /* This part can be modified by SRT. */

    if (rl==RNUM[Gc_AN]){
      for (i=0; i<=tno1; i++){
	for (j=0; j<=tno1; j++){
	  tp[i][j] = 0.0;
	}
      }
    }

    else{
      for (i=0; i<RCN_Rank[rl]; i++){
	for (j=0; j<RCN_Rank[rl]; j++){
	  sum = 0.0;
	  for (k=0; k<RCN_Rank[rl+1]; k++){
	    sum += G0[rl+1][i][k]*be[rl+1][k][j];
	  }
	  tp[i][j] = sum;
	}
      }
    }

    /*  An + R_{n+1}*B_{n+1} */

    for (i=0; i<RCN_Rank[rl]; i++){
      for (j=0; j<RCN_Rank[rl]; j++){
        tp1[i][j] = al[rl][i][j] + tp[i][j];
      }
    }

    /* [ An + R_{n+1}*B_{n+1} ]^{-1} */    
    
    for (i=0; i<RCN_Rank[rl]; i++){
      for (j=0; j<RCN_Rank[rl]; j++){
        printf("%15.12f ",tp1[i][j]);
      }
      printf("\n");
    }

    Inverse_Mat(RCN_Rank[rl]-1,tp1);

    /* - Bn^{t} * [Z - An - R_{n+1}*B_{n+1}]^{-1} */    

    for (i=0; i<RCN_Rank[rl-1]; i++){
      for (j=0; j<RCN_Rank[rl]; j++){
	sum = 0.0;
	for (k=0; k<RCN_Rank[rl]; k++){
                     /* !! */
	  sum += be[rl][k][i]*tp1[k][j];
	}
	G0[rl][i][j] = -sum;
      }
    }

  } /* rl */

  /****************************************************
        calculate off-diagonal Green's functions
             by G_{0,n} = G_{0,n-1} * R_{n}
  ****************************************************/

  for (rl=2; rl<=RNUM[Gc_AN]; rl++){

    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[rl]; j++){
	sum = 0.0;
	for (k=0; k<RCN_Rank[rl-1]; k++){
	  sum += G0[rl-1][i][k]*G0[rl][k][j];
	}
	tp[i][j] = sum;
      }
    }

    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[rl]; j++){
	G0[rl][i][j] = tp[i][j];
      }
    }
  }

  /****************************************************
             inverse Lanczos transformation
  ****************************************************/

  for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){         /*   i     */

    ig = natn[Gc_AN][i];
    ino1 = Spe_Total_CNO[WhatSpecies[ig]] - 1;
    Anum = MP[i];

    for (p=0; p<=ino1; p++){                            /*  beta   */
      for (q=0; q<=tno1; q++){                          /*  alpha  */
	sum = 0.0;
	for (rl=0; rl<=RNUM[Gc_AN]; rl++){              /*   n     */
	  for (m=0; m<RCN_Rank[rl]; m++){               /*   mu    */
	    sum += G0[rl][q][m]*U[rl][m][Anum+p];
	  }
	}
	IOLP[Mc_AN][i][q][p] = sum;
      }
    }
  }
}















void Inverse_Mat(int n, double **A0)
{
  static int i,j,N;
  static INTEGER lda,info,lwork;

  N = n + 1;

  lda = N;
  lwork = N;

  /****************************************************
      A0 -> A
  ****************************************************/

  for (i=0;i<=n;i++) {
    for (j=0;j<=n;j++) {
       A[i*(n+1)+j]= A0[i][j];
    }
  }

  /****************************************************
                call zgetrf_() in clapack
  ****************************************************/

  F77_NAME(dgetrf,DGETRF)(&N, &N, A, &lda, ipiv, &info);

  if (info!=0){
    printf("error in dgetrf_() which is called from IS_Lanczos info=%2d\n",info);
  }

  /****************************************************
                Call dgetri_() in clapack
  ****************************************************/

  F77_NAME(dgetri,DGETRI)(&N, A, &lda, ipiv, work, &lwork, &info);
  if (info!=0){
    printf("error in dgetri_() which is called from IS_Lanczos info=%2d\n",info);
  }

  /****************************************************
               A -> A0
  ****************************************************/

  for (i=0; i<=n; i++) {
    for (j=0; j<=n; j++) {
      A0[i][j] = A[i*(n+1)+j];
    }
  }

}

