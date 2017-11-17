/**********************************************************************

  jx.c:

  This program code calculates spin-spin interaction 
  coupling constant J between the selected two atoms

  Log of jx.c:
     30/Aug/2003  Released by Myung Joon Han (supervised by Prof. J. Yu)
      7/Dec/2003  Modified by Taisuke Ozaki
     03/Mar/2011  Modified by Fumiyuki Ishii for MPI 
***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h> 
#include "read_scfout.h"
#include "lapack_prototypes.h"
#include "f77func.h"
#include "mpi.h"


#define Host_ID       0         /* ID of the host CPU in MPI */

#define printout  0    /* 0:off, 1:on */
#define PI   3.1415926535897932384626

static void Eigen_lapack(double **a, double *ko, int n);
static void EigenBand_lapack(dcomplex **A, double *W, int N);

static void jx_cluster1(int argc, char *argv[]);
static void calc_J_cluster1(int First_Atom, int Second_Atom, int *MPF,
                            double ***C, double **ko,
                            int F_TNumOrbs, int F_TNumOrbs3);

static void k_inversion(int i,  int j,  int k, 
                        int mi, int mj, int mk, 
                        int *ii, int *ij, int *ik );
static void Overlap_Band(double ****OLP,
                         dcomplex **S,int *MP,
                         double k1, double k2, double k3);
static void Hamiltonian_Band(double ****RH, dcomplex **H, int *MP,
                             double k1, double k2, double k3);
static void calc_J_band1(int First_Atom, int Second_Atom, int *MPF,
                         dcomplex ***C, double **ko,
                         int F_TNumOrbs, int F_TNumOrbs3, double Jmat[2]);
static void dtime(double *t);


MPI_Comm comm1;



int main(int argc, char *argv[]) 
{
  static int ct_AN,h_AN,Gh_AN,i,j,k,TNO1,TNO2;  
  static int spin,Rn;

  int spin1,spin2;   /* dummy variables for loop    */
  int ii,jj,ll;      /* dummy variables for loop    */
  int II,JJ ;        /* dummy variables for loop    */
  int Anum;          /* temporary variable          */
  int Number_Choo ;  /* the number of Choosen atom  */
  int TNumOrbs,TNumOrbs3;  
  int F_TNumOrbs,F_TNumOrbs3;
  int SPI,SPJ;
  int optEV ;        /* variable to be used for option of eigenvector-printing */
  int *Full_atom ;   /* array to store the glbal number of Full atoms */
  int *MPF ;         /* Full system version of MP */
  int First_Atom;    /* the first atom i in Jij  */
  int Second_Atom;   /* the second atom j in Jij */
  int end_switch;    /* switch in the evalusion loop of J */
  double ***FullHks; /* full Hamiltonian */
  double **FullOLP;  /* full overlap     */

  /* for MPI */
  int numprocs,myid,ID,ID1, *arpo;
  int T_knum,S_knum,E_knum,kloop0,kloop,num_kloop0;
  double TStime,TEtime;
  int Nk[4],atmij[3];
  int *T_k_op; 
  double *T_KGrids1, *T_KGrids2, *T_KGrids3, tmp4;
  int *Tkmesh2ID;
  double *JrTk;
  double *JiTk;
  /* for MPI end */

  static double snum_i,snum_j,snum_k;
  static double *KGrids1,*KGrids2,*KGrids3;
  static double Jmat[2],Jr,Ji;
  static int ij,ik,P_min;
  static int knum_i,knum_j,knum_k;
  static int knum_switch;
  static int ***k_op;

  static dcomplex **H,**S,**C,***Coes;
  static double **ko, *M1;
  static double k1,k2,k3; 
  static double OLP_eigen_cut = 1.0e-10;

  /* variable & arrays for PART-2; same with that of Cluster_DFT.c */
  static int l,n,n2,n1,i1,j1,l1;
  static double sum,sumi,sum1;

  /* MPI initialize */

  MPI_Status stat;
  MPI_Request request;
  MPI_Init(&argc, &argv);
  comm1 = MPI_COMM_WORLD;
  MPI_Comm_size(comm1,&numprocs);
  MPI_Comm_rank(comm1,&myid);

  /* MPI initialize end */

  char *s_vec[20];
  dtime(&TStime); 

  if (myid==Host_ID){
    printf("\n********************************************************************\n"); 
    printf("********************************************************************\n"); 
    printf(" jx: code for calculating an effective exchange coupling constant J \n"); 
    printf(" Copyright (C), 2003, Myung Joon Han, Jaejun Yu, and Taisuke Ozaki \n"); 
    printf(" This is free software, and you are welcome to         \n"); 
    printf(" redistribute it under the constitution of the GNU-GPL.\n");
    printf("********************************************************************\n"); 
    printf("********************************************************************\n"); 
  }
  read_scfout(argv);

  s_vec[0]="Recursion"; s_vec[1]="Cluster"; s_vec[2]="Band";
  s_vec[3]="NEGF";      s_vec[4]="DC";      s_vec[5]="GDC";
  s_vec[6]="Cluster2";

  if (myid==Host_ID){
    printf(" Previous eigenvalue solver = %s\n",s_vec[Solver-1]);
    printf(" atomnum                    = %i\n",atomnum);
    printf(" ChemP                      = %15.12f (Hartree)\n",ChemP);
    printf(" E_Temp                     = %15.12f (K)\n",E_Temp);
  }
  if      (Solver==2 || Solver==7) jx_cluster1(argc, argv);
  else if (Solver==3) { 

    /**** Calculation of J for bulk system  ****/

    /************************************************************************* 
     Calculation flow in the evaluation of J:

     PART-1 : set full Hamiltonian and overlap
     PART-2 : the generalized eigenvalue problem HC = ESC
     PART-3 : Calculation of J
    **************************************************************************/

    /************************************************************************* 
     PART-1 :  set full Hamiltonian and overlap
    **************************************************************************/

    if (myid==Host_ID){
      printf(" \nEvaluation of J based on band calculation\n");
    }
    Full_atom = (int*)malloc(sizeof(int)*atomnum); 
    MPF = (int*)malloc(sizeof(int)*(atomnum+1)); 

    /* read Full_atom */
    for (i=0; i<atomnum; i++){
      Full_atom[i] = i+1;
    }

    /*********************************************
   make an array MPF which specifies the starting
   position of atom II in the martix such as
   a full but small Hamiltonian
    *********************************************/

    Anum = 1;
    for (i=0; i<atomnum; i++){
      ct_AN = Full_atom[i];
      Anum = Anum + Total_NumOrbs[ct_AN];    
    }
    F_TNumOrbs = Anum - 1;
    F_TNumOrbs3 = F_TNumOrbs + 3;  

    /*********************************************
         get knum_i, knum_j, and knum_k
    *********************************************/

    do {
      if (myid==Host_ID){
	printf(" Specify the number of k-grids (e.g 4 4 4)  \n");
	scanf("%d %d %d",&Nk[1],&Nk[2],&Nk[3]);
	printf(" Nk1, Nk2, Nk3 = %d %d %d  \n", Nk[1],Nk[2], Nk[3]);
      }
      MPI_Bcast(Nk, 4, MPI_INT, 0, comm1); 
      knum_i=Nk[1];
      knum_j=Nk[2];
      knum_k=Nk[3];


      if ( knum_i<=0 || knum_j<=0 || knum_k<=0 ){
	printf("invalid number\n");
	knum_switch = 0;
      }
      else{
	knum_switch = 1;
      }
    } while(knum_switch==0);

    /****************************************************
   allocation of arrays:

   double    KGrids1[knum_i]
   double    KGrids2[knum_j]
   double    KGrids3[knum_k]
   double    k_op[knum_i][knum_j][knum_k]
    ****************************************************/

    ko = (double**)malloc(sizeof(double*)*(SpinP_switch+1));
    for (i=0; i<=SpinP_switch; i++){
      ko[i] = (double*)malloc(sizeof(double)*F_TNumOrbs3);
    }

    H = (dcomplex**)malloc(sizeof(dcomplex*)*F_TNumOrbs3);
    for (j=0; j<F_TNumOrbs3; j++){
      H[j] = (dcomplex*)malloc(sizeof(dcomplex)*F_TNumOrbs3);
    }

    S = (dcomplex**)malloc(sizeof(dcomplex*)*F_TNumOrbs3);
    for (i=0; i<F_TNumOrbs3; i++){
      S[i] = (dcomplex*)malloc(sizeof(dcomplex)*F_TNumOrbs3);
    }

    M1 = (double*)malloc(sizeof(double)*F_TNumOrbs3);

    C = (dcomplex**)malloc(sizeof(dcomplex*)*F_TNumOrbs3);
    for (j=0; j<F_TNumOrbs3; j++){
      C[j] = (dcomplex*)malloc(sizeof(dcomplex)*F_TNumOrbs3);
    }

    Coes = (dcomplex***)malloc(sizeof(dcomplex**)*(SpinP_switch+1));
    for (spin=0; spin<=SpinP_switch; spin++){
      Coes[spin] = (dcomplex**)malloc(sizeof(dcomplex*)*F_TNumOrbs3);
      for (j=0; j<F_TNumOrbs3; j++){
	Coes[spin][j] = (dcomplex*)malloc(sizeof(dcomplex)*F_TNumOrbs3);
      }
    }

    KGrids1 = (double*)malloc(sizeof(double)*knum_i);
    KGrids2 = (double*)malloc(sizeof(double)*knum_j);
    KGrids3 = (double*)malloc(sizeof(double)*knum_k);

    k_op=(int***)malloc(sizeof(int**)*(knum_i));
    for (i=0;i<knum_i; i++) {
      k_op[i]=(int**)malloc(sizeof(int*)*(knum_j));
      for (j=0;j<knum_j; j++) {
	k_op[i][j]=(int*)malloc(sizeof(int)*(knum_k));
      }
    }

    /*********************************************
   make an array MPF which specifies the starting
   position of atom II in the martix such as
   a full but small Hamiltonian
    *********************************************/

    snum_i = (double)knum_i;
    snum_j = (double)knum_j;
    snum_k = (double)knum_k;

    for (i=0; i<=(knum_i-1); i++){
      if (knum_i==1) k1 = 0.0;
      else           k1 = -0.5 + (2.0*(double)i+1.0)/(2.0*snum_i);
      KGrids1[i]=k1;
    }
    for (i=0; i<=(knum_j-1); i++){
      if (knum_j==1) k1 = 0.0;
      else           k1 = -0.5 + (2.0*(double)i+1.0)/(2.0*snum_j);
      KGrids2[i]=k1;
    }
    for (i=0; i<=(knum_k-1); i++){
      if (knum_k==1) k1 = 0.0;
      else           k1 = -0.5 + (2.0*(double)i+1.0)/(2.0*snum_k);
      KGrids3[i]=k1;
    }




    if (printout==1){
      printf("   KGrids1: ");
      for (i=0;i<=knum_i-1;i++) printf("%f ",KGrids1[i]);
      printf("\n");
      printf("   KGrids2: ");
      for (i=0;i<=knum_j-1;i++) printf("%f ",KGrids2[i]);
      printf("\n");
      printf("   KGrids3: ");
      for (i=0;i<=knum_k-1;i++) printf("%f ",KGrids3[i]);
      printf("\n");
    }

    /**************************************************************
   k_op[i][j][k]: weight of DOS 
                 =0   no calc.
                 =1   G-point
                 =2   which has k<->-k point
   Now , only the relation, E(k)=E(-k), is used. 
    *************************************************************/

    for (i=0;i<=knum_i-1;i++) {
      for (j=0;j<=knum_j-1;j++) {
	for (k=0;k<=knum_k-1;k++) {
	  k_op[i][j][k]=-999;
	}
      }
    }
    for (i=0;i<=knum_i-1;i++) {
      for (j=0;j<=knum_j-1;j++) {
	for (k=0;k<=knum_k-1;k++) {
	  if ( k_op[i][j][k]==-999 ) {
	    k_inversion(i,j,k,knum_i,knum_j,knum_k,&ii,&ij,&ik);
	    if ( i==ii && j==ij && k==ik ) {
	      /*  printf("* : %d %d %d (%f %f %f)\n",i,j,k,
		  KGrids1[i], KGrids2[j], KGrids3[k]);
	      */

	      k_op[i][j][k]=1;
	    }
	    else {
	      k_op[i][j][k]=2;
	      k_op[ii][ij][ik]=0;
	    }
	  }
	} /* k */
      } /* j */
    } /* i */


    /***********************************
       one-dimentionalize for MPI
    ************************************/


    T_knum = 0;
    for (i=0; i<=(knum_i-1); i++){
      for (j=0; j<=(knum_j-1); j++){
	for (k=0; k<=(knum_k-1); k++){
          T_knum++;
	}
      }
    }

    Tkmesh2ID= (int*)malloc(sizeof(int)*T_knum);
    JrTk = (double*)malloc(sizeof(double)*T_knum);
    JiTk = (double*)malloc(sizeof(double)*T_knum);
    T_KGrids1 = (double*)malloc(sizeof(double)*T_knum);
    T_KGrids2 = (double*)malloc(sizeof(double)*T_knum);
    T_KGrids3 = (double*)malloc(sizeof(double)*T_knum);
    arpo = (int*)malloc(sizeof(int)*numprocs);

    /* set T_KGrids1,2,3 and T_k_op */

    T_knum = 0;
    for (i=0; i<=(knum_i-1); i++){
      for (j=0; j<=(knum_j-1); j++){
	for (k=0; k<=(knum_k-1); k++){

	  T_KGrids1[T_knum] = KGrids1[i];
	  T_KGrids2[T_knum] = KGrids2[j];
	  T_KGrids3[T_knum] = KGrids3[k];
	  JrTk[T_knum]=0.0; 
	  JiTk[T_knum]=0.0; 

          T_knum++;
	}
      }
    }

    /****  allocate k-point  into proccessors *****/

    if (T_knum<=myid){
      S_knum = -10;
      E_knum = -100;
      num_kloop0 = 1;
    }
    else if (T_knum<numprocs) {
      S_knum = myid;
      E_knum = myid;
      num_kloop0 = 1;
    }
    else {
      tmp4 = (double)T_knum/(double)numprocs; 
      num_kloop0 = (int)tmp4 + 1;
      S_knum = (int)((double)myid*(tmp4+0.0001)); 
      E_knum = (int)((double)(myid+1)*(tmp4+0.0001)) - 1;
      if (myid==(numprocs-1)) E_knum = T_knum - 1;
      if (E_knum<0)           E_knum = 0;
    }

    /****************************************************
                         start calc J
    ****************************************************/

    /** printf("\n"); **/

    end_switch = 0;
    do {

      /* specify two atoms */ 
      if (myid==Host_ID){

	printf("\n Specify two atoms (e.g 1 2, quit: 0 0) \n");
	scanf("%d %d",&atmij[1],&atmij[2]);

	if (atmij[1] > 0){
	  printf(" Atom:i,j = %d, %d \n",atmij[1], atmij[2]);
	} 
	else{ 
          printf(" ******* exit \n");
        }
      }

      MPI_Bcast(atmij, 3, MPI_INT, 0, comm1); 
      MPI_Barrier(comm1); 
      First_Atom  = atmij[1];
      Second_Atom = atmij[2];

      if      (First_Atom==0 && Second_Atom==0) {
        end_switch = 1;       
        MPI_Finalize();
        exit(1);
      }
      else if (First_Atom==0 && Second_Atom!=0 || 
	       First_Atom!=0 && Second_Atom==0){
	printf(" Invalid atom\n");
      }
      else if (atomnum<First_Atom || atomnum<Second_Atom){
	printf(" Invalid atom\n");
      }
      else{

	/****************************************************
                          calculation of J
	****************************************************/

	/* loop for k-grids */

	Jr = 0.0;
	Ji = 0.0;

	for (kloop0=0; kloop0<num_kloop0; kloop0++){

	  kloop = S_knum + kloop0;

	  arpo[myid] = -1;
	  if (S_knum<=kloop && kloop<=E_knum) arpo[myid] = kloop;
	  for (ID=0; ID<numprocs; ID++){
	    MPI_Bcast(&arpo[ID], 1, MPI_INT, ID, comm1);
	  }
	  if (myid==Host_ID){
          }

	  kloop = arpo[myid];

	  for (ID=0; ID<numprocs; ID++){
	    Tkmesh2ID[kloop] = ID;
	    MPI_Bcast(Tkmesh2ID, T_knum, MPI_INT, ID, comm1);
	  }

	  if (0<=kloop){

	    k1 = T_KGrids1[kloop];
	    k2 = T_KGrids2[kloop];
	    k3 = T_KGrids3[kloop];
	    JrTk[kloop]=0.0;           
	    JiTk[kloop]=0.0;           

	    /****************************************************
		 MPF indicates the starting position of
			 atom i in arraies H and S
	    ****************************************************/

	    Overlap_Band(OLP,S,MPF,k1,k2,k3);

	    n = S[0][0].r;
	    EigenBand_lapack(S,ko[0],n);

	    P_min = 1;
	    for (l=1; l<=n; l++) if (ko[0][l]<OLP_eigen_cut) P_min = l + 1;
	    for (l=P_min; l<=n; l++) M1[l] = 1.0/sqrt(ko[0][l]);

	    for (spin=0; spin<=SpinP_switch; spin++){

	      Hamiltonian_Band(Hks[spin], H, MPF, k1, k2, k3);

	      /****************************************************
			   M1 * U^t * H * U * M1
	      ****************************************************/
 
	      for (i1=1; i1<=n; i1++){
		for (j1=P_min; j1<=n; j1++){
		  sum = 0.0; sumi=0.0;
		  for (l=1; l<=n; l++){
		    sum  += ( H[i1][l].r*S[l][j1].r
			    - H[i1][l].i*S[l][j1].i)*M1[j1];
		    sumi += ( H[i1][l].r*S[l][j1].i
			    + H[i1][l].i*S[l][j1].r)*M1[j1];
		  }
		  C[i1][j1].r = sum;
		  C[i1][j1].i = sumi;
		}
	      }     

	      for (i1=P_min; i1<=n; i1++){
		for (j1=1; j1<=n; j1++){
		  sum = 0.0; sumi=0.0;
		  for (l=1; l<=n; l++){
		    sum  +=  M1[i1]*( S[l][i1].r*C[l][j1].r +
				      S[l][i1].i*C[l][j1].i );
		    sumi +=  M1[i1]*( S[l][i1].r*C[l][j1].i -
				      S[l][i1].i*C[l][j1].r );
		  }
		  H[i1][j1].r = sum;
		  H[i1][j1].i = sumi;
		}
	      }     

	      for (i1=P_min; i1<=n; i1++){
		for (j1=P_min; j1<=n; j1++){
		  C[i1-(P_min-1)][j1-(P_min-1)] = H[i1][j1];
		}
	      }

	      n1 = n - (P_min - 1);
	      EigenBand_lapack(C,ko[spin],n1);

	      /****************************************************
		   Transformation to the original eigenvectors.
			   NOTE JRCAT-244p and JAIST-2122p 
	      ****************************************************/

	      for (i1=1; i1<=n; i1++){
		for (j1=1; j1<=n; j1++){
		  H[i1][j1].r = 0.0;
		  H[i1][j1].i = 0.0;
		}
	      }

	      for (i1=1; i1<=n; i1++){
		for (j1=1; j1<=n1; j1++){
		  sum = 0.0; sumi=0.0;
		  for (l=P_min; l<=n; l++){
		    sum  +=  S[i1][l].r*M1[l]*C[l-(P_min-1)][j1].r
		           - S[i1][l].i*M1[l]*C[l-(P_min-1)][j1].i;
		    sumi +=  S[i1][l].r*M1[l]*C[l-(P_min-1)][j1].i
		           + S[i1][l].i*M1[l]*C[l-(P_min-1)][j1].r;
		  }
		  Coes[spin][i1][j1].r = sum;
		  Coes[spin][i1][j1].i = sumi;
		}
	      }


	    } /* spin */

	    /************************************************************************* 
	       PART-3 : Calculation of J
  
		1) V_i,alpha,beta = < i,alpha | V_i | j,beta >
			          = 0.5 * ( H_i,alpha,i,beta(0) - H_i,alpha,i,beta(1) )

		2) V'_i = SUM_alpha,beta{ C_i,alpha * V_alpha,beta * C_i,beta }

		3) calculation of J from V'
	    **************************************************************************/

	    /* calculation of J */ 

	    calc_J_band1(First_Atom, Second_Atom, MPF, Coes, ko,
			 n1, F_TNumOrbs3, Jmat);

	    JrTk[kloop] += Jmat[0];
	    JiTk[kloop] += Jmat[1];

	    /*  printf("\n calculate kloop, Jr=%d, %f10.5 \n",kloop, JrTk[kloop]); */
	  } /* if(kloop =>0) */
        } /* kloop0  for MPI */


  
	MPI_Barrier(comm1); 
	kloop=0;
	for (i=0; i<=(knum_i-1); i++){
	  for (j=0; j<=(knum_j-1); j++){
	    for (k=0; k<=(knum_k-1); k++){
	      ID1 = Tkmesh2ID[kloop];
	      MPI_Bcast(&JrTk[kloop], 1, MPI_DOUBLE, ID1, comm1);
	      MPI_Bcast(&JiTk[kloop], 1, MPI_DOUBLE, ID1, comm1);
	      kloop++;
	    }
	  }
	}
     

	for (i=0; i<T_knum; i++){
          Jr += JrTk[i];
	}
	Jr = Jr/(double)(knum_i*knum_j*knum_k);
	if (myid==Host_ID){
	  printf("\n J_ij between %ith atom and %ith atom is %15.12f cm^{-1}\n",
		 First_Atom, Second_Atom, Jr);
	}
      } /* else */

    } while(end_switch==0);
    /*********************************************
    freeing of arrays:

    int *Full_atom;
    int *MPF;
    double ***FullHks ;  
    double **FullOLP ;  
    double **ko;
    double *M1;
    double **B;
    double ***C;
    double **D;
    *********************************************/

    free(Full_atom);
    free(MPF);

    for (spin=0; spin<=SpinP_switch; spin++){
      free(ko[spin]);
    }
    free(ko);

    for (j=0; j<F_TNumOrbs3; j++){
      free(H[j]);
    }
    free(H);

    for (i=0; i<F_TNumOrbs3; i++){
      free(S[i]);
    }
    free(S);

    free(M1);

    for (j=0; j<F_TNumOrbs3; j++){
      free(C[j]);
    }
    free(C);

    for (spin=0; spin<=SpinP_switch; spin++){
      for (j=0; j<F_TNumOrbs3; j++){
	free(Coes[spin][j]);
      }
      free(Coes[spin]);
    }
    free(Coes);

    free(KGrids1);
    free(KGrids2);
    free(KGrids3);

    for (i=0;i<knum_i; i++) {
      for (j=0;j<knum_j; j++) {
	free(k_op[i][j]);
      }
      free(k_op[i]);
    }
    free(k_op);
  }

  /* print message */

  MPI_Barrier(comm1); 

  dtime(&TEtime); 
  /* if (myid==Host_ID){ */
  printf(" \nElapsed time = %lf (s) for myid=%3d\n",TEtime-TStime,myid);fflush(stdout);
  /* printf(" \nElapsed time = %lf (s) with processors=%3d\n",TEtime-TStime,numprocs);fflush(stdout); */
  /*  } */
  MPI_Barrier(comm1); 
  /** printf(" \nThe calculation was finished normally in myid=%2d.\n",myid);fflush(stdout); **/

  /* MPI_Finalize */

  MPI_Finalize();

  exit(1);
}


void jx_cluster1(int argc, char *argv[]) 
{

  static int ct_AN,h_AN,Gh_AN,i,j,TNO1,TNO2;  
  static int spin,Rn;

  int spin1,spin2;   /* dummy variables for loop    */
  int ii,jj,ll;      /* dummy variables for loop    */
  int II,JJ ;        /* dummy variables for loop    */
  int Anum;          /* temporary variable          */
  int Number_Choo ;  /* the number of Choosen atom  */
  int TNumOrbs,TNumOrbs3;  
  int F_TNumOrbs,F_TNumOrbs3;
  int SPI,SPJ;
  int optEV ;        /* variable to be used for option of eigenvector-printing */
  int *Full_atom ;   /* array to store the glbal number of Full atoms */
  int *MPF ;         /* Full system version of MP */
  int First_Atom;    /* the first atom i in Jij  */
  int Second_Atom;   /* the second atom j in Jij */
  int end_switch;    /* switch in the evalusion loop of J */
  double ***FullHks; /* full Hamiltonian */
  double **FullOLP;  /* full overlap     */

  /* variable & arrays for PART-2; same with that of Cluster_DFT.c */
  static int l,n,n2,n1,i1,j1,k1,l1;
  static double **ko, *M1;
  static double **B, ***C, **D;
  static double sum,sum1;                                                      

  int myid,numprocs;

  MPI_Comm_size(comm1,&numprocs);
  MPI_Comm_rank(comm1,&myid);

  if (1<numprocs){
    if (myid==Host_ID){
      printf("\nMPI parallelization is not supported for the cluster mode.\n\n");
    }
    MPI_Finalize();
    exit(1);
  }

  /************************************************************************* 
     Calculation flow in the evaluation of J:

     PART-1 : set full Hamiltonian and overlap
     PART-2 : the generalized eigenvalue problem HC = ESC
     PART-3 : Calculation of J
  **************************************************************************/

  /************************************************************************* 
     PART-1 :  set full Hamiltonian and overlap
  **************************************************************************/

  printf(" \nEvaluation of J based on cluster calculation\n");

  Full_atom = (int*)malloc(sizeof(int)*atomnum); 
  MPF = (int*)malloc(sizeof(int)*(atomnum+1)); 

  /* read Full_atom */
  for (i=0; i<atomnum; i++){
    Full_atom[i] = i+1;
  }

  /*********************************************
   make an array MPF which specifies the starting
   position of atom II in the martix such as
   a full but small Hamiltonian
  *********************************************/

  Anum = 1;
  for (i=0; i<atomnum; i++){
    MPF[i+1] = Anum;
    ct_AN = Full_atom[i];
    Anum = Anum + Total_NumOrbs[ct_AN];    
  }
  F_TNumOrbs = Anum - 1;
  F_TNumOrbs3 = F_TNumOrbs + 3;  

  /*****************************************
    allocation of arrays:

    double ***FullHks ;  
    double **FullOLP ;  
  ******************************************/

  FullHks = (double***)malloc(sizeof(double**)*(SpinP_switch+1)); 
  for (spin=0; spin<=SpinP_switch; spin++){
    FullHks[spin] = (double**)malloc(sizeof(double*)*F_TNumOrbs3); 
    for (i=0; i<F_TNumOrbs3; i++){
      FullHks[spin][i] = (double*)malloc(sizeof(double)*F_TNumOrbs3); 
    }
  }

  FullOLP = (double**)malloc(sizeof(double*)*F_TNumOrbs3); 
  for (i=0; i<F_TNumOrbs3; i++){
    FullOLP[i] = (double*)malloc(sizeof(double)*F_TNumOrbs3); 
  }

  /*****************************************
   rearragement of data:
   Hks to FullHks
   OLP to FullOLP

   FullHks [1 to F_TNumOrbs] [1 to F_TNumOrbs]
   FullOLP [1 to F_TNumOrbs] [1 to F_TNumOrbs]
  ******************************************/

  for (spin=0; spin<=SpinP_switch; spin++){
    for (i=0; i<F_TNumOrbs3; i++){
      for (j=0; j<F_TNumOrbs3; j++){
        FullHks[spin][i][j] = 0.0;
      }
    }
  }

  for (i=0; i<F_TNumOrbs3; i++){
    for (j=0; j<F_TNumOrbs3; j++){
      FullOLP[i][j] = 0.0;
    }
  }

  for (spin=0; spin<=SpinP_switch; spin++){
    for (II=0; II<atomnum; II++){
      SPI = MPF[II+1];  
      for (JJ=0; JJ<atomnum; JJ++){
        SPJ = MPF[JJ+1];  
	ct_AN = Full_atom[II] ;
	TNO1 = Total_NumOrbs[ct_AN];
	for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
	  Gh_AN = natn[ct_AN][h_AN];
	  if (Gh_AN == Full_atom[JJ]){
	    Rn = ncn[ct_AN][h_AN];
	    TNO2 = Total_NumOrbs[Gh_AN];
	    for (i=0; i<TNO1; i++){
              for (j=0; j<TNO2; j++){
                FullHks[spin][i+SPI][j+SPJ] = Hks[spin][ct_AN][h_AN][i][j]; 
                FullOLP[i+SPI][j+SPJ] = OLP[ct_AN][h_AN][i][j]; 
	      }
            }
	  }
	}
      }
    }
  }
  
  /************************************************************************* 
     PART-2 : solve the generalized eigenvalue problem HC = ESC
  **************************************************************************/

  /*******************************************
   allocation of arrays:

   double ko[SpinP_switch+1][F_TNumOrbs3];
   double M1[F_TNumOrbs3];
   double B[F_TNumOrbs3][F_TNumOrbs3];
   double C[SpinP_switch+1][F_TNumOrbs3][F_TNumOrbs3];
   double D[F_TNumOrbs3][F_TNumOrbs3];
  ********************************************/

  ko = (double**)malloc(sizeof(double*)*(SpinP_switch+1));
  for (spin=0; spin<=SpinP_switch; spin++){
    ko[spin] = (double*)malloc(sizeof(double)*F_TNumOrbs3); 
  }
  M1 = (double*)malloc(sizeof(double)*F_TNumOrbs3); 

  B = (double**)malloc(sizeof(double*)*F_TNumOrbs3); 
  for (i=0; i<F_TNumOrbs3; i++){
    B[i] = (double*)malloc(sizeof(double)*F_TNumOrbs3); 
  }

  C = (double***)malloc(sizeof(double**)*(SpinP_switch+1)); 
  for (spin=0; spin<=SpinP_switch; spin++){
    C[spin] = (double**)malloc(sizeof(double*)*F_TNumOrbs3); 
    for (i=0; i<F_TNumOrbs3; i++){
      C[spin][i] = (double*)malloc(sizeof(double)*F_TNumOrbs3); 
    }
  }

  D = (double**)malloc(sizeof(double*)*F_TNumOrbs3); 
  for (i=0; i<F_TNumOrbs3; i++){
    D[i] = (double*)malloc(sizeof(double)*F_TNumOrbs3); 
  }

  /*******************************************
   diagonalize the overlap matrix

   first 
   FullOLP -> OLP matrix
   after call Eigen_lapack
   FullOLP -> eigenvectors of OLP matrix
  ********************************************/

  printf(" Diagonalize the overlap matrix\n");

  Eigen_lapack(FullOLP,ko[0],F_TNumOrbs);

  /*
    for (l=1; l<=F_TNumOrbs; l++){
    printf("l ko %2d %15.12f\n",l,ko[0][l]);
    }
  */

  /* check ill-conditioned eigenvalues */

  for (l=1; l<=F_TNumOrbs; l++){
    if (ko[0][l]<1.0e-14) {
      printf("Found ill-conditioned eigenvalues\n");
      printf("Stopped calculation\n");
      exit(1);
    }
  }

  for (l=1; l<=F_TNumOrbs; l++){
    M1[l] = 1.0/sqrt(ko[0][l]);
  }

  /****************************************************
   Calculations of eigenvalues for up and down spins
  ****************************************************/

  n = F_TNumOrbs;

  for (spin=0; spin<=SpinP_switch; spin++){

    printf(" Diagonalize the Hamiltonian for spin=%2d\n",spin);

    for (i1=1; i1<=n; i1++){
      for (j1=1; j1<=n; j1++){
        sum = 0.0;
        for (l=1; l<=n; l++){
	  sum = sum + FullHks[spin][i1][l]*FullOLP[l][j1]*M1[j1]; 
        }
        C[spin][i1][j1] = sum;
      }
    }

    for (i1=1; i1<=n; i1++){
      for (j1=1; j1<=n; j1++){
        sum = 0.0;
        for (l=1; l<=n; l++){
	  sum = sum + M1[i1]*FullOLP[l][i1]*C[spin][l][j1];
        }
        B[i1][j1] = sum;
      }
    }

    for (i1=1; i1<=n; i1++){
      for (j1=1; j1<=n; j1++){
        D[i1][j1] = B[i1][j1];    
      }
    }

    Eigen_lapack(D,ko[spin],n);

    /****************************************************
        Transformation to the original eigen vectors.
                          NOTE 244P
    ****************************************************/

    for (i1=1; i1<=n; i1++){
      for (j1=1; j1<=n; j1++){
        C[spin][i1][j1] = 0.0;
      }
    }

    for (i1=1; i1<=n; i1++){
      for (j1=1; j1<=n; j1++){
        sum = 0.0;
        for (l=1; l<=n; l++){
          sum = sum + FullOLP[i1][l]*M1[l]*D[l][j1];       
        }
        C[spin][i1][j1] = sum;
      }
    }
  } /* spin */

  /* printing out eigenvalues and the eigenvectors */

  if (printout==1){
    for (spin=0; spin<=SpinP_switch; spin++){
      printf("\nspin=%i \n",spin);                                       
      for (i=1; i<=F_TNumOrbs; i++){
	printf("%ith eigenvalue of HC=eSC: %15.12f\n",i,ko[spin][i]);
      }
    }
    printf("\nDo you want eigenvectors also? (yes:1 / no:0)");
    scanf("%i",&optEV);
    if (optEV == 1){
      for (spin=0; spin<=SpinP_switch; spin++){
	printf("\nspin=%i \n",spin);                                     
	for (i=1; i<=F_TNumOrbs; i++){
	  printf("%ith eigenvector: ",i);
	  printf("{");
	  for (j=1; j<=F_TNumOrbs; j++){
	    printf("%15.12f,",C[spin][i][j]);
	  }
	  printf("}\n");
	}
      }
    }
  }

  /************************************************************************* 

     PART-3 : Calculation of J
     
     1) V_i,alpha,beta = < i,alpha | V_i | j,beta >
                       = 0.5 * ( H_i,alpha,i,beta(0) - H_i,alpha,i,beta(1) )

     2) V'_i = SUM_alpha,beta{ C_i,alpha * V_alpha,beta * C_i,beta }

     3) calculation of J from V'
  **************************************************************************/

  printf("\n");

  end_switch = 0;
  do {

    /* specify two atoms */ 
    printf(" Specify two atoms (e.g 1 2, quit: 0 0)  ");
    scanf("%d %d",&First_Atom,&Second_Atom);

    if      (First_Atom==0 && Second_Atom==0) end_switch = 1;       
    else if (First_Atom==0 && Second_Atom!=0 || 
             First_Atom!=0 && Second_Atom==0){
      printf(" Invalid atom\n");
    }
    else if (atomnum<First_Atom || atomnum<Second_Atom){
      printf(" Invalid atom\n");
    }
    else{
      /* calculation of J */ 
      calc_J_cluster1(First_Atom, Second_Atom, MPF, C, ko, F_TNumOrbs, F_TNumOrbs3);
    } 

  } while(end_switch==0);

  /*********************************************
    freeing of arrays:

    int *Full_atom;
    int *MPF;
    double ***FullHks ;  
    double **FullOLP ;  
    double **ko;
    double *M1;
    double **B;
    double ***C;
    double **D;
  *********************************************/

  free(Full_atom);
  free(MPF);

  for (spin=0; spin<=SpinP_switch; spin++){
    for (i=0; i<F_TNumOrbs3; i++){
      free(FullHks[spin][i]);
    }
    free(FullHks[spin]);
  }
  free(FullHks);

  for (i=0; i<F_TNumOrbs3; i++){
    free(FullOLP[i]);
  }
  free(FullOLP);

  for (spin=0; spin<=SpinP_switch; spin++){
    free(ko[spin]);
  }
  free(ko);

  free(M1);

  for (i=0; i<F_TNumOrbs3; i++){
    free(B[i]);
  }
  free(B);

  for (spin=0; spin<=SpinP_switch; spin++){
    for (i=0; i<F_TNumOrbs3; i++){
      free(C[spin][i]);
    }
    free(C[spin]);
  }
  free(C);

  for (i=0; i<F_TNumOrbs3; i++){
    free(D[i]);
  }
  free(D);
}  





void calc_J_band1(int First_Atom, int Second_Atom, int *MPF, dcomplex ***C, double **ko,
                  int F_TNumOrbs, int F_TNumOrbs3, double Jmat[2])
{
  static int ct_AN,h_AN,Gh_AN,i,j,TNO1,TNO2;  
  static int spin,Rn;
  static int spin1,spin2;   /* dummy variables for loop                               */
  static int ii,jj,ll;      /* dummy variables for loop                               */
  static int II,JJ ;        /* dummy variables for loop                               */
  static int Number_Choo ;  /* the number of Choosen atom                             */
  static int *MP ;          /* array which specify a head position of a full matrix   */
  static int Anum;          /* temporary variable                                     */
  static int *Choo_atom ;   /* array to store the glbal number of choosen atoms       */
  static int TNumOrbs,TNumOrbs3;  
  static int SPI,SPJ;
  static int NOJ,NOI;

  static double sum_r,sum_i,sum1_r,sum1_i;
  static double ***SmallHks;
  static double **SmallOLP;
  static double **Vi,**Vj;         /* array for V_i and V_j */
  static double ****VVi_r,****VVi_i,****VVj_r,****VVj_i ; /* array for V' */
  static double *IVi_r,*IVi_i,*IVj_r,*IVj_i ; /* intermediate memory for calculating VVi, VVj */
  static double **Fftn;        /* Fermi function */
  static double J_ij_r,J_ij_i;         /* final exchange-interaction value  */
  static double kB=0.000003166813628;  /* Boltzman constant (Hatree/K) */
  static double dFftn;
  static double dko;

  /*********************************************
    allocation of arrays:

    int Choo_atom[Number_Choo];
    int MP[Number_Choo];
  *********************************************/

  Number_Choo = 2;

  Choo_atom = (int*)malloc(sizeof(int)*Number_Choo); 
  MP = (int*)malloc(sizeof(int)*Number_Choo); 

  Choo_atom[0] = First_Atom;
  Choo_atom[1] = Second_Atom;

  /*********************************************
   make an array MP which specify the starting
   position of atom II in the martix such as
   a full but small Hamiltonian
  *********************************************/

  Anum = 1;
  for (i=0; i<Number_Choo; i++){
    MP[i] = Anum;
    ct_AN = Choo_atom[i];
    Anum = Anum + Total_NumOrbs[ct_AN];    
  }
  TNumOrbs = Anum - 1;
  TNumOrbs3 = TNumOrbs + 3;  

  II = Choo_atom[0];
  NOI = Total_NumOrbs[II] ; /* this gives orbital number of i-atom (or 1st atom) */
  JJ = Choo_atom[1];
  NOJ = Total_NumOrbs[JJ] ; /* this gives orbital number of j-atom (or 2nd atom) */

  /*********************************************
    allocation of arrays:

    double ***SmallHks ;  
    double **SmallOLP ;  
  *********************************************/

  SmallHks = (double***)malloc(sizeof(double**)*(SpinP_switch+1)); 
  for (spin=0; spin<=SpinP_switch; spin++){
    SmallHks[spin] = (double**)malloc(sizeof(double*)*TNumOrbs3); 
    for (i=0; i<TNumOrbs3; i++){
      SmallHks[spin][i] = (double*)malloc(sizeof(double)*TNumOrbs3); 
    }
  }

  SmallOLP = (double**)malloc(sizeof(double*)*TNumOrbs3); 
  for (i=0; i<TNumOrbs3; i++){
    SmallOLP[i] = (double*)malloc(sizeof(double)*TNumOrbs3); 
  }

  /*********************************************
   Hamiltonian and overlap for selected two atom
  *********************************************/

  for (spin=0; spin<=SpinP_switch; spin++){
    for (II=0; II<Number_Choo; II++){
      SPI = MP[II];
      for (JJ=0; JJ<Number_Choo; JJ++){
        SPJ = MP[JJ];
	ct_AN = Choo_atom[II] ;
	TNO1 = Total_NumOrbs[ct_AN];
	for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
	  Gh_AN = natn[ct_AN][h_AN];
          Rn = ncn[ct_AN][h_AN];
	  if (Gh_AN == Choo_atom[JJ] && Rn==0){
	    TNO2 = Total_NumOrbs[Gh_AN];
	    for (i=0; i<TNO1; i++){
	      for (j=0; j<TNO2; j++){
		SmallHks[spin][i+SPI][j+SPJ] = Hks[spin][ct_AN][h_AN][i][j]; 
	      }
	    }
	  }
	}
      }
    }
  }

  if (printout==1){
    for (spin=0; spin<=SpinP_switch; spin++){
      printf("spin=%i Full Hamiltonian matrix of selected atoms\n",spin);
      for (i=1; i<=TNumOrbs; i++){
        for (j=1; j<=TNumOrbs; j++){
          printf("%7.3f ",SmallHks[spin][i][j] );           
        }
        printf("\n");
      }
    }

    for (spin=0; spin<=SpinP_switch; spin++){
      printf("spin=%i Hamiltonian matrix of 1st atom\n",spin);
      for (i=1; i<=NOI; i++){
	for (j=1; j<=NOI; j++){
	  printf("%7.3f ",SmallHks[spin][i][j]);
	}
	printf("\n");
      }
    }

    for (spin=0; spin<=SpinP_switch; spin++){
      printf("spin=%i Hamiltonian matrix of 2nd atom\n",spin);
      for (i=(NOI+1); i<=TNumOrbs; i++){
	for (j=(NOI+1); j<=TNumOrbs; j++){
	  printf("%7.3f ",SmallHks[spin][i][j]);              
	}
	printf("\n");
      }
    }
  }

  /*********************************************
        calculation of V from H(0)-H(1)
  *********************************************/

  /* memory allocation of Vi and Vj  */
  Vi = (double**)malloc(sizeof(double*)*(NOI+3)); 
  for (II=0; II<(NOI+3); II++){
    Vi[II] = (double*)malloc(sizeof(double)*(NOI+3)); 
  }

  Vj = (double**)malloc(sizeof(double*)*(NOJ+3)); 
  for (II=0; II<(NOJ+3); II++){
    Vj[II] = (double*)malloc(sizeof(double)*(NOJ+3)); 
  }

  /*  V[i] = 0.5*( H_i,up - H_i,down )  */
  for (i=1; i<=NOI; i++){
    for (j=1; j<=NOI; j++){
      Vi[i][j] = 0.5 * ( SmallHks[0][i][j] - SmallHks[1][i][j]);    
    }
  }

  /*  V[j] = 0.5*( H_j,up - H_j,down )  */
  for (i=1; i<=NOJ; i++){
    for (j=1; j<=NOJ; j++){
      Vj[i][j] = 0.5*(SmallHks[0][i+NOI][j+NOI] - SmallHks[1][i+NOI][j+NOI]);
    }
  }

  if (printout==1){
    printf("Vi for 1st atom\n");
    for (i=1; i<=NOI; i++){
      for (j=1; j<=NOI; j++){
        printf("%9.6f ",Vi[i][j]);
      }
      printf("\n");
    }

    printf("Vj  2nd atom\n");
    for (i=1; i<=NOJ; i++){
      for (j=1; j<=NOJ; j++){
        printf("%9.6f ",Vj[i][j]);
      }
      printf("\n");
    }
  }

  /*********************************************
           calculation of VVi and VVj
  *********************************************/

  /* memory allocation */
  VVi_r = (double****)malloc(sizeof(double***)*(SpinP_switch+1)); 
  for (spin1=0; spin1<=SpinP_switch; spin1++){
    VVi_r[spin1] = (double***)malloc(sizeof(double**)*(SpinP_switch+1)); 
    for (spin2=0; spin2<=1; spin2++){
      VVi_r[spin1][spin2] = (double**)malloc(sizeof(double*)*F_TNumOrbs3); 
      for (i=0; i<F_TNumOrbs3; i++){
	VVi_r[spin1][spin2][i] = (double*)malloc(sizeof(double)*F_TNumOrbs3); 
      }
    }
  }

  VVi_i = (double****)malloc(sizeof(double***)*(SpinP_switch+1)); 
  for (spin1=0; spin1<=SpinP_switch; spin1++){
    VVi_i[spin1] = (double***)malloc(sizeof(double**)*(SpinP_switch+1)); 
    for (spin2=0; spin2<=1; spin2++){
      VVi_i[spin1][spin2] = (double**)malloc(sizeof(double*)*F_TNumOrbs3); 
      for (i=0; i<F_TNumOrbs3; i++){
	VVi_i[spin1][spin2][i] = (double*)malloc(sizeof(double)*F_TNumOrbs3); 
      }
    }
  }

  VVj_r = (double****)malloc(sizeof(double***)*(SpinP_switch+1)); 
  for (spin1=0; spin1<=SpinP_switch; spin1++){
    VVj_r[spin1] = (double***)malloc(sizeof(double**)*(SpinP_switch+1)); 
    for (spin2=0; spin2<=SpinP_switch; spin2++){
      VVj_r[spin1][spin2] = (double**)malloc(sizeof(double*)*F_TNumOrbs3); 
      for (i=0; i<F_TNumOrbs3; i++){
	VVj_r[spin1][spin2][i] = (double*)malloc(sizeof(double)*F_TNumOrbs3); 
      }
    }
  }

  VVj_i = (double****)malloc(sizeof(double***)*(SpinP_switch+1)); 
  for (spin1=0; spin1<=SpinP_switch; spin1++){
    VVj_i[spin1] = (double***)malloc(sizeof(double**)*(SpinP_switch+1)); 
    for (spin2=0; spin2<=SpinP_switch; spin2++){
      VVj_i[spin1][spin2] = (double**)malloc(sizeof(double*)*F_TNumOrbs3); 
      for (i=0; i<F_TNumOrbs3; i++){
	VVj_i[spin1][spin2][i] = (double*)malloc(sizeof(double)*F_TNumOrbs3); 
      }
    }
  }

  IVi_r = (double*)malloc(sizeof(double)*(NOI+3)); 
  IVi_i = (double*)malloc(sizeof(double)*(NOI+3)); 
  IVj_r = (double*)malloc(sizeof(double)*(NOJ+3)); 
  IVj_i = (double*)malloc(sizeof(double)*(NOJ+3)); 

  /* calculation VVi */    
  for(spin1=0; spin1<=SpinP_switch; spin1++){
    for (spin2=0; spin2<=SpinP_switch; spin2++){
      for(i=1; i<=F_TNumOrbs; i++){
	for(j=1; j<=F_TNumOrbs; j++){
	  sum1_r = 0.0 ;
	  sum1_i = 0.0 ;

	  for(ll=1; ll<=NOI; ll++){	    
	    for(ii=1; ii<=NOI; ii++){	      
	      sum_r = 0.0;
	      sum_i = 0.0;

	      for(jj=1; jj<=NOI; jj++){
		sum_r +=  C[spin1][(jj-1+MPF[Choo_atom[0]])][i].r*Vi[jj][ll] ;
		sum_i += -C[spin1][(jj-1+MPF[Choo_atom[0]])][i].i*Vi[jj][ll] ;
	      }
	      IVi_r[ii] = sum_r ;
	      IVi_i[ii] = sum_i ;
	    }
	    sum1_r += IVi_r[ll]*C[spin2][(ll-1+MPF[Choo_atom[0]])][j].r
	      - IVi_i[ll]*C[spin2][(ll-1+MPF[Choo_atom[0]])][j].i ;

	    sum1_i += IVi_r[ll]*C[spin2][(ll-1+MPF[Choo_atom[0]])][j].i
	      + IVi_i[ll]*C[spin2][(ll-1+MPF[Choo_atom[0]])][j].r ;
	  }
	  VVi_r[spin1][spin2][i][j] = sum1_r ;
	  VVi_i[spin1][spin2][i][j] = sum1_i ;
	}
      }
    }
  }

  /* calculation VVj */    
  for(spin1=0; spin1<=SpinP_switch; spin1++){
    for (spin2=0; spin2<=SpinP_switch; spin2++){
      for(i=1; i<=F_TNumOrbs; i++){
	for(j=1; j<=F_TNumOrbs; j++){
	  sum1_r = 0.0 ;
	  sum1_i = 0.0 ;

	  for(ll=1; ll<=NOJ; ll++){
	    for(ii=1; ii<=NOJ; ii++){
	      sum_r = 0.0;
	      sum_i = 0.0;

	      for(jj=1; jj<=NOJ; jj++){
		sum_r +=  C[spin1][(jj-1+MPF[Choo_atom[1]])][i].r*Vj[jj][ll] ;
		sum_i += -C[spin1][(jj-1+MPF[Choo_atom[1]])][i].i*Vj[jj][ll] ;
	      }
	      IVj_r[ii] = sum_r ;
	      IVj_i[ii] = sum_i ;
	    }

	    sum1_r = sum1_r
	      + IVj_r[ll]*C[spin2][(ll-1+MPF[Choo_atom[1]])][j].r 
	      - IVj_i[ll]*C[spin2][(ll-1+MPF[Choo_atom[1]])][j].i ; 

	    sum1_i = sum1_i
	      + IVj_r[ll]*C[spin2][(ll-1+MPF[Choo_atom[1]])][j].i 
	      + IVj_i[ll]*C[spin2][(ll-1+MPF[Choo_atom[1]])][j].r ;
	  }
	  VVj_r[spin1][spin2][i][j] = sum1_r ;
	  VVj_i[spin1][spin2][i][j] = sum1_i ;

	}
      }
    }
  }

  /*********************************************
              calculation of J_ij
  *********************************************/
    
  /* allocation of memory */
  Fftn = (double**)malloc(sizeof(double*)*(SpinP_switch+1)); 
  for (spin=0; spin<=SpinP_switch; spin++){
    Fftn[spin] = (double*)malloc(sizeof(double)*F_TNumOrbs3); 
  }

  /* the Fftn-array for Fermi function */
  for (spin=0; spin<=SpinP_switch; spin++){
    for (i=1; i<=F_TNumOrbs; i++){
      Fftn[spin][i] = 1.0/( exp( (ko[spin][i]- ChemP)/(kB*E_Temp) ) + 1.0 )  ;
    }
  }

  /* J_ij calculation from VVi & VVj */

  J_ij_r = 0.0; 
  J_ij_i = 0.0; 
  for (i=1; i<=F_TNumOrbs; i++){    
    for (j=1; j<=F_TNumOrbs; j++){    
      dFftn = Fftn[0][i] - Fftn[1][j] ;
      dko = ko[1][j] - ko[0][i] ;
      J_ij_r = J_ij_r + 0.5*dFftn *
	(  VVj_r[0][1][i][j] * VVi_r[1][0][j][i]
	   - VVj_i[0][1][i][j] * VVi_i[1][0][j][i] ) / dko ;

      J_ij_i = J_ij_i + 0.5*dFftn *
	(  VVj_r[0][1][i][j] * VVi_i[1][0][j][i]
	   + VVj_i[0][1][i][j] * VVi_r[1][0][j][i] ) / dko ;
    }
  }

  /* unit conversion: Hartree to cm^{-1} */

  J_ij_r = 2.194746*100000.0*J_ij_r;
  J_ij_i = 2.194746*100000.0*J_ij_i;

  Jmat[0] = J_ij_r;
  Jmat[1] = J_ij_i;

  /*********************************************
    freeing of arrays:

    double ***SmallHks ;  
    double **SmallOLP ;  
    double **Vi;
    double **Vj;
    double ****VVi;
    double ****VVj;
    double *IVi;
    double *IVj;
    double **Fftn;
  *********************************************/

  free(Choo_atom);
  free(MP);

  for (spin=0; spin<=SpinP_switch; spin++){
    for (i=0; i<TNumOrbs3; i++){
      free(SmallHks[spin][i]);
    }
    free(SmallHks[spin]);
  }
  free(SmallHks);

  for (i=0; i<TNumOrbs3; i++){
    free(SmallOLP[i]);
  }
  free(SmallOLP);


  for (II=0; II<(NOI+3); II++){
    free(Vi[II]);
  }
  free(Vi);

  for (II=0; II<(NOJ+3); II++){
    free(Vj[II]);
  }
  free(Vj);

  for (spin1=0; spin1<=SpinP_switch; spin1++){
    for (spin2=0; spin2<=1; spin2++){
      for (i=0; i<F_TNumOrbs3; i++){
        free(VVi_r[spin1][spin2][i]);
      }
      free(VVi_r[spin1][spin2]);
    }
    free(VVi_r[spin1]);
  }
  free(VVi_r);

  for (spin1=0; spin1<=SpinP_switch; spin1++){
    for (spin2=0; spin2<=1; spin2++){
      for (i=0; i<F_TNumOrbs3; i++){
        free(VVi_i[spin1][spin2][i]);
      }
      free(VVi_i[spin1][spin2]);
    }
    free(VVi_i[spin1]);
  }
  free(VVi_i);

  for (spin1=0; spin1<=SpinP_switch; spin1++){
    for (spin2=0; spin2<=SpinP_switch; spin2++){
      for (i=0; i<(F_TNumOrbs3); i++){
	free(VVj_r[spin1][spin2][i]);
      }
      free(VVj_r[spin1][spin2]);
    }
    free(VVj_r[spin1]);
  }
  free(VVj_r);

  for (spin1=0; spin1<=SpinP_switch; spin1++){
    for (spin2=0; spin2<=SpinP_switch; spin2++){
      for (i=0; i<(F_TNumOrbs3); i++){
	free(VVj_i[spin1][spin2][i]);
      }
      free(VVj_i[spin1][spin2]);
    }
    free(VVj_i[spin1]);
  }
  free(VVj_i);

  free(IVi_r);
  free(IVi_i);
  free(IVj_r);
  free(IVj_i);


  for (spin=0; spin<=SpinP_switch; spin++){
    free(Fftn[spin]);
  }
  free(Fftn);
}    



void calc_J_cluster1(int First_Atom, int Second_Atom, int *MPF, double ***C, double **ko,
                     int F_TNumOrbs, int F_TNumOrbs3)
{
  static int ct_AN,h_AN,Gh_AN,i,j,TNO1,TNO2;  
  static int spin,Rn;
  static int spin1,spin2;   /* dummy variables for loop                               */
  static int ii,jj,ll;      /* dummy variables for loop                               */
  static int II,JJ ;        /* dummy variables for loop                               */
  static int Number_Choo ;  /* the number of Choosen atom                             */
  static int *MP ;          /* array which specify a head position of a full matrix   */
  static int Anum;          /* temporary variable                                     */
  static int *Choo_atom ;   /* array to store the glbal number of choosen atoms       */
  static int TNumOrbs,TNumOrbs3;  
  static int SPI,SPJ;
  static int NOJ,NOI;

  static double sum,sum1;                                                      
  static double ***SmallHks;
  static double **SmallOLP;
  static double **Vi,**Vj, ****VVi, ****VVj ;   /* array for V_i,V_j and V' */
  static double *IVi, *IVj ;   /* intermediate memory for calculating VVi, VVj */
  static double **Fftn;        /* Fermi function */
  static double J_ij;          /* final exchange-interaction value  */
  static double kB=0.000003166813628;  /* Boltzman constant (Hatree/K) */
  static double dFftn;
  static double dko;

  /*********************************************
    allocation of arrays:

    int Choo_atom[Number_Choo];
    int MP[Number_Choo];
  *********************************************/

  Number_Choo = 2;

  Choo_atom = (int*)malloc(sizeof(int)*Number_Choo); 
  MP = (int*)malloc(sizeof(int)*Number_Choo); 

  Choo_atom[0] = First_Atom;
  Choo_atom[1] = Second_Atom;

  /*********************************************
   make an array MP which specify the starting
   position of atom II in the martix such as
   a full but small Hamiltonian
  *********************************************/

  Anum = 1;
  for (i=0; i<Number_Choo; i++){
    MP[i] = Anum;
    ct_AN = Choo_atom[i];
    Anum = Anum + Total_NumOrbs[ct_AN];    
  }
  TNumOrbs = Anum - 1;
  TNumOrbs3 = TNumOrbs + 3;  

  II = Choo_atom[0];
  NOI = Total_NumOrbs[II] ; /* this gives orbital number of i-atom (or 1st atom) */
  JJ = Choo_atom[1];
  NOJ = Total_NumOrbs[JJ] ; /* this gives orbital number of j-atom (or 2nd atom) */

  /*********************************************
    allocation of arrays:

    double ***SmallHks ;  
    double **SmallOLP ;  
  *********************************************/

  SmallHks = (double***)malloc(sizeof(double**)*(SpinP_switch+1)); 
  for (spin=0; spin<=SpinP_switch; spin++){
    SmallHks[spin] = (double**)malloc(sizeof(double*)*TNumOrbs3); 
    for (i=0; i<TNumOrbs3; i++){
      SmallHks[spin][i] = (double*)malloc(sizeof(double)*TNumOrbs3); 
    }
  }

  SmallOLP = (double**)malloc(sizeof(double*)*TNumOrbs3); 
  for (i=0; i<TNumOrbs3; i++){
    SmallOLP[i] = (double*)malloc(sizeof(double)*TNumOrbs3); 
  }

  /*********************************************
   Hamiltonian and overlap for selected two atom
  *********************************************/

  for (spin=0; spin<=SpinP_switch; spin++){
    for (II=0; II<Number_Choo; II++){
      SPI = MP[II];  
      for (JJ=0; JJ<Number_Choo; JJ++){
        SPJ = MP[JJ];  
	ct_AN = Choo_atom[II] ;
	TNO1 = Total_NumOrbs[ct_AN];
	for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
	  Gh_AN = natn[ct_AN][h_AN];
	  if (Gh_AN == Choo_atom[JJ]){
	    Rn = ncn[ct_AN][h_AN];
	    TNO2 = Total_NumOrbs[Gh_AN];
	    for (i=0; i<TNO1; i++){
	      for (j=0; j<TNO2; j++){
		SmallHks[spin][i+SPI][j+SPJ] = Hks[spin][ct_AN][h_AN][i][j]; 
	      }
	    }
	  }
	}
      }
    }
  }

  if (printout==1){
    for (spin=0; spin<=SpinP_switch; spin++){
      printf("spin=%i Full Hamiltonian matrix of selected atoms\n",spin);
      for (i=1; i<=TNumOrbs; i++){
        for (j=1; j<=TNumOrbs; j++){
          printf("%7.3f ",SmallHks[spin][i][j] );           
        }
        printf("\n");
      }
    }

    for (spin=0; spin<=SpinP_switch; spin++){
      printf("spin=%i Hamiltonian matrix of 1st atom\n",spin);
      for (i=1; i<=NOI; i++){
	for (j=1; j<=NOI; j++){
	  printf("%7.3f ",SmallHks[spin][i][j]);
	}
	printf("\n");
      }
    }

    for (spin=0; spin<=SpinP_switch; spin++){
      printf("spin=%i Hamiltonian matrix of 2nd atom\n",spin);
      for (i=(NOI+1); i<=TNumOrbs; i++){
	for (j=(NOI+1); j<=TNumOrbs; j++){
	  printf("%7.3f ",SmallHks[spin][i][j]);              
	}
	printf("\n");
      }
    }
  }

  /*********************************************
        calculation of V from H(0)-H(1)
  *********************************************/

  /* memory allocation of Vi and Vj  */
  Vi = (double**)malloc(sizeof(double*)*(NOI+3)); 
  for (II=0; II<(NOI+3); II++){
    Vi[II] = (double*)malloc(sizeof(double)*(NOI+3)); 
  }

  Vj = (double**)malloc(sizeof(double*)*(NOJ+3)); 
  for (II=0; II<(NOJ+3); II++){
    Vj[II] = (double*)malloc(sizeof(double)*(NOJ+3)); 
  }

  /*  V[i] = 0.5*( H_i,up - H_i,down )  */
  for (i=1; i<=NOI; i++){
    for (j=1; j<=NOI; j++){
      Vi[i][j] = 0.5 * ( SmallHks[0][i][j] - SmallHks[1][i][j]);    
    }
  }

  /*  V[j] = 0.5*( H_j,up - H_j,down )  */
  for (i=1; i<=NOJ; i++){
    for (j=1; j<=NOJ; j++){
      Vj[i][j] = 0.5*(SmallHks[0][i+NOI][j+NOI] - SmallHks[1][i+NOI][j+NOI]);
    }
  }

  if (printout==1){
    printf("Vi for 1st atom\n");
    for (i=1; i<=NOI; i++){
      for (j=1; j<=NOI; j++){
        printf("%9.6f ",Vi[i][j]);
      }
      printf("\n");
    }

    printf("Vj  2nd atom\n");
    for (i=1; i<=NOJ; i++){
      for (j=1; j<=NOJ; j++){
        printf("%9.6f ",Vj[i][j]);
      }
      printf("\n");
    }
  }

  /*********************************************
           calculation of VVi and VVj
  *********************************************/

  /* memory allocation */
  VVi = (double****)malloc(sizeof(double***)*(SpinP_switch+1)); 
  for (spin1=0; spin1<=SpinP_switch; spin1++){
    VVi[spin1] = (double***)malloc(sizeof(double**)*(SpinP_switch+1)); 
    for (spin2=0; spin2<=1; spin2++){
      VVi[spin1][spin2] = (double**)malloc(sizeof(double*)*(F_TNumOrbs3)); 
      for (i=0; i<F_TNumOrbs3; i++){
	VVi[spin1][spin2][i] = (double*)malloc(sizeof(double)*(F_TNumOrbs3)); 
      }
    }
  }

  VVj = (double****)malloc(sizeof(double***)*(SpinP_switch+1)); 
  for (spin1=0; spin1<=SpinP_switch; spin1++){
    VVj[spin1] = (double***)malloc(sizeof(double**)*(SpinP_switch+1)); 
    for (spin2=0; spin2<=SpinP_switch; spin2++){
      VVj[spin1][spin2] = (double**)malloc(sizeof(double*)*(F_TNumOrbs3)); 
      for (i=0; i<(F_TNumOrbs3); i++){
	VVj[spin1][spin2][i] = (double*)malloc(sizeof(double)*(F_TNumOrbs3)); 
      }
    }
  }

  IVi = (double*)malloc(sizeof(double)*(NOI+3)); 
  IVj = (double*)malloc(sizeof(double)*(NOJ+3)); 

  /* calculation VVi */    
  for(spin1=0; spin1<=SpinP_switch; spin1++){
    for (spin2=0; spin2<=SpinP_switch; spin2++){
      for(i=1; i<=F_TNumOrbs; i++){
	for(j=1; j<=F_TNumOrbs; j++){
	  sum1 = 0.0 ;
	  for(ll=1; ll<=NOI; ll++){	    
	    for(ii=1; ii<=NOI; ii++){	      
	      sum = 0.0;
	      for(jj=1; jj<=NOI; jj++){
		sum = sum + C[spin1][(jj-1+MPF[Choo_atom[0]])][i]*Vi[jj][ll] ;
	      }
	      IVi[ii] = sum ;
	    }
	    sum1 = sum1 + IVi[ll]*C[spin2][(ll-1+MPF[Choo_atom[0]])][j] ;
	  }
	  VVi[spin1][spin2][i][j] = sum1 ;
	}
      }
    }
  }

  /* calculation VVj */    
  for(spin1=0; spin1<=SpinP_switch; spin1++){
    for (spin2=0; spin2<=SpinP_switch; spin2++){
      for(i=1; i<=F_TNumOrbs; i++){
	for(j=1; j<=F_TNumOrbs; j++){
	  sum1 = 0.0 ;
	  for(ll=1; ll<=NOJ; ll++){	    
	    for(ii=1; ii<=NOJ; ii++){	      
	      sum = 0.0;
	      for(jj=1; jj<=NOJ; jj++){
		sum = sum + C[spin1][(jj-1+MPF[Choo_atom[1]])][i]*Vj[jj][ll] ;
	      }
	      IVj[ii] = sum ;
	    }
	    sum1 = sum1 + IVj[ll]*C[spin2][(ll-1+MPF[Choo_atom[1]])][j] ;
	  }
	  VVj[spin1][spin2][i][j] = sum1 ;
	}
      }
    }
  }

  /*********************************************
              calculation of J_ij
  *********************************************/
    
  /* allocation of memory */
  Fftn = (double**)malloc(sizeof(double*)*(SpinP_switch+1)); 
  for (spin=0; spin<=SpinP_switch; spin++){
    Fftn[spin] = (double*)malloc(sizeof(double)*F_TNumOrbs3); 
  }

  J_ij = 0.0; 

  /* the Fftn-array for Fermi function */
  for (spin=0; spin<=SpinP_switch; spin++){
    for (i=1; i<=F_TNumOrbs; i++){
      Fftn[spin][i] = 1.0/( exp( (ko[spin][i]- ChemP)/(kB*E_Temp) ) + 1.0 )  ;
    }
  }

  /* J_ij calculation from VVi & VVj */
  for (i=1; i<=F_TNumOrbs; i++){    
    for (j=1; j<=F_TNumOrbs; j++){    
      dFftn = Fftn[0][i] - Fftn[1][j] ;
      dko = ko[1][j] - ko[0][i] ;
      J_ij = J_ij + 0.5*dFftn * VVj[0][1][i][j] * VVi[1][0][j][i] / dko ;
    }
  }

  /* unit conversion: Hartree to cm^{-1} */

  J_ij = 2.194746*100000.0*J_ij;

  /* printing the final result : J_ij*/

  printf(" J_ij between %ith atom and %ith atom is %15.12f cm^{-1}\n",
	 Choo_atom[0], Choo_atom[1], J_ij);

  /*********************************************
    freeing of arrays:

    double ***SmallHks ;  
    double **SmallOLP ;  
    double **Vi;
    double **Vj;
    double ****VVi;
    double ****VVj;
    double *IVi;
    double *IVj;
    double **Fftn;
  *********************************************/

  free(MP);
  free(Choo_atom);

  for (spin=0; spin<=SpinP_switch; spin++){
    for (i=0; i<TNumOrbs3; i++){
      free(SmallHks[spin][i]);
    }
    free(SmallHks[spin]);
  }
  free(SmallHks);

  for (i=0; i<TNumOrbs3; i++){
    free(SmallOLP[i]);
  }
  free(SmallOLP);


  for (II=0; II<(NOI+3); II++){
    free(Vi[II]);
  }
  free(Vi);

  for (II=0; II<(NOJ+3); II++){
    free(Vj[II]);
  }
  free(Vj);

  for (spin1=0; spin1<=SpinP_switch; spin1++){
    for (spin2=0; spin2<=1; spin2++){
      for (i=0; i<F_TNumOrbs3; i++){
        free(VVi[spin1][spin2][i]);
      }
      free(VVi[spin1][spin2]);
    }
    free(VVi[spin1]);
  }
  free(VVi);


  for (spin1=0; spin1<=SpinP_switch; spin1++){
    for (spin2=0; spin2<=SpinP_switch; spin2++){
      for (i=0; i<(F_TNumOrbs3); i++){
	free(VVj[spin1][spin2][i]);
      }
      free(VVj[spin1][spin2]);
    }
    free(VVj[spin1]);
  }
  free(VVj);


  free(IVi);
  free(IVj);


  for (spin=0; spin<=SpinP_switch; spin++){
    free(Fftn[spin]);
  }
  free(Fftn);
}    











































void k_inversion(int i,  int j,  int k, 
                 int mi, int mj, int mk, 
                 int *ii, int *ij, int *ik )
{
  *ii= mi-i-1;
  *ij= mj-j-1;
  *ik= mk-k-1;
}




void Overlap_Band(double ****OLP,
                  dcomplex **S,int *MP,
                  double k1, double k2, double k3)
{
  static int i,j,wanA,wanB,tnoA,tnoB,Anum,Bnum,NUM,GA_AN,LB_AN,GB_AN;
  static int l1,l2,l3,Rn,n2;
  static double **S1,**S2;
  static double kRn,si,co,s;

  Anum = 1;
  for (i=1; i<=atomnum; i++){
    MP[i] = Anum;
    Anum += Total_NumOrbs[i];
  }
  NUM = Anum - 1;

  /****************************************************
                       Allocation
  ****************************************************/

  n2 = NUM + 2;

  S1 = (double**)malloc(sizeof(double*)*n2);
  for (i=0; i<n2; i++){
    S1[i] = (double*)malloc(sizeof(double)*n2);
  }

  S2 = (double**)malloc(sizeof(double*)*n2);
  for (i=0; i<n2; i++){
    S2[i] = (double*)malloc(sizeof(double)*n2);
  }

  /****************************************************
                       set overlap
  ****************************************************/

  S[0][0].r = NUM;

  for (i=1; i<=NUM; i++){
    for (j=1; j<=NUM; j++){
      S1[i][j] = 0.0;
      S2[i][j] = 0.0;
    }
  }

  for (GA_AN=1; GA_AN<=atomnum; GA_AN++){
    tnoA = Total_NumOrbs[GA_AN];
    Anum = MP[GA_AN];

    for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
      GB_AN = natn[GA_AN][LB_AN];
      Rn = ncn[GA_AN][LB_AN];
      tnoB = Total_NumOrbs[GB_AN];

      l1 = atv_ijk[Rn][1];
      l2 = atv_ijk[Rn][2];
      l3 = atv_ijk[Rn][3];
      kRn = k1*(double)l1 + k2*(double)l2 + k3*(double)l3;

      si = sin(2.0*PI*kRn);
      co = cos(2.0*PI*kRn);
      Bnum = MP[GB_AN];
      for (i=0; i<tnoA; i++){
	for (j=0; j<tnoB; j++){
	  s = OLP[GA_AN][LB_AN][i][j];
	  S1[Anum+i][Bnum+j] = S1[Anum+i][Bnum+j] + s*co;
	  S2[Anum+i][Bnum+j] = S2[Anum+i][Bnum+j] + s*si;
	}
      }
    }
  }

  for (i=1; i<=NUM; i++){
    for (j=1; j<=NUM; j++){
      S[i][j].r =  S1[i][j];
      S[i][j].i =  S2[i][j];
    }
  }

  /****************************************************
                       free arrays
  ****************************************************/

  for (i=0; i<n2; i++){
    free(S1[i]);
    free(S2[i]);
  }
  free(S1);
  free(S2);

}


void Hamiltonian_Band(double ****RH, dcomplex **H, int *MP,
                      double k1, double k2, double k3)
{
  static int i,j,wanA,wanB,tnoA,tnoB,Anum,Bnum,NUM,GA_AN,LB_AN,GB_AN;
  static int l1,l2,l3,Rn,n2;
  static double **H1,**H2;
  static double kRn,si,co,h;

  Anum = 1;
  for (i=1; i<=atomnum; i++){
    MP[i] = Anum;
    Anum += Total_NumOrbs[i];
  }
  NUM = Anum - 1;

  /****************************************************
                       Allocation
  ****************************************************/

  n2 = NUM + 2;

  H1 = (double**)malloc(sizeof(double*)*n2);
  for (i=0; i<n2; i++){
    H1[i] = (double*)malloc(sizeof(double)*n2);
  }

  H2 = (double**)malloc(sizeof(double*)*n2);
  for (i=0; i<n2; i++){
    H2[i] = (double*)malloc(sizeof(double)*n2);
  }

  /****************************************************
                    set Hamiltonian
  ****************************************************/

  H[0][0].r = 2.0*NUM;
  for (i=1; i<=NUM; i++){
    for (j=1; j<=NUM; j++){
      H1[i][j] = 0.0;
      H2[i][j] = 0.0;
    }
  }

  for (GA_AN=1; GA_AN<=atomnum; GA_AN++){
    tnoA = Total_NumOrbs[GA_AN];
    Anum = MP[GA_AN];

    for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
      GB_AN = natn[GA_AN][LB_AN];
      Rn = ncn[GA_AN][LB_AN];
      tnoB = Total_NumOrbs[GB_AN];

      l1 = atv_ijk[Rn][1];
      l2 = atv_ijk[Rn][2];
      l3 = atv_ijk[Rn][3];
      kRn = k1*(double)l1 + k2*(double)l2 + k3*(double)l3;

      si = sin(2.0*PI*kRn);
      co = cos(2.0*PI*kRn);
      Bnum = MP[GB_AN];
      for (i=0; i<tnoA; i++){
	for (j=0; j<tnoB; j++){
	  h = RH[GA_AN][LB_AN][i][j];
	  H1[Anum+i][Bnum+j] = H1[Anum+i][Bnum+j] + h*co;
	  H2[Anum+i][Bnum+j] = H2[Anum+i][Bnum+j] + h*si;
	}
      }
    }
  }

  for (i=1; i<=NUM; i++){
    for (j=1; j<=NUM; j++){
      H[i][j].r = H1[i][j];
      H[i][j].i = H2[i][j];
    }
  }

  /****************************************************
                       free arrays
  ****************************************************/

  for (i=0; i<n2; i++){
    free(H1[i]);
    free(H2[i]);
  }
  free(H1);
  free(H2);
}



void Eigen_lapack(double **a, double *ko, int n)
{
  /* input:  n;
     input:  a[n][n];  matrix A

     output: a[n][n]; eigevectors
     output: ko[n];   eigenvalues  */
    
  static char *name="Eigen_lapack";

  char  *JOBZ="V";
  char  *RANGE="A";
  char  *UPLO="L";

  int LDA=n;
  double VL,VU; /* dummy */
  int IL,IU; /* dummy */
  double ABSTOL=1.0e-10;
  int M;

  double *A,*Z;
  int LDZ=n;
  int LWORK;
  double *WORK;
  int *IWORK;

  int *IFAIL, INFO;

  int i,j;

  A=(double*)malloc(sizeof(double)*n*n);
  Z=(double*)malloc(sizeof(double)*n*n);

  LWORK=n*8;
  WORK=(double*)malloc(sizeof(double)*LWORK);
  IWORK=(int*)malloc(sizeof(int)*n*5);
  IFAIL=(int*)malloc(sizeof(int)*n);

  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      A[i*n+j]= a[i+1][j+1];
    }
  }

#if 0
  printf("A=\n");
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      printf("%f ",A[i*n+j]);
    }
    printf("\n");
  }
  fflush(stdout);
#endif

  F77_NAME(dsyevx,DSYEVX)( JOBZ, RANGE, UPLO, &n, A, &LDA, &VL, &VU, &IL, &IU,
			   &ABSTOL, &M, ko, Z, &LDZ, WORK, &LWORK, IWORK,
			   IFAIL, &INFO ); 

  /* store eigenvectors */
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      /*  a[i+1][j+1]= Z[i*n+j]; */
      a[j+1][i+1]= Z[i*n+j];

    }
  }

  /* shift ko by 1 */
  for (i=n;i>=1;i--){
    ko[i]= ko[i-1];
  }

  if (INFO>0) {
    printf("\n%s: error in dsyevx_, info=%d\n\n",name,INFO);
  }
  if (INFO<0) {
    printf("%s: info=%d\n",name,INFO);
    exit(10);
  }
   
  free(IFAIL); free(IWORK); free(WORK); free(Z); free(A);

}





/*
  void EigenBand_lapack(dcomplex **A, double *W, int N)
  {
  static char *JOBZ="V";
  static char *RANGE="A";
  static char *UPLO="L";

  int LDA=N;
  double VL,VU; 
  int IL,IU; 
  double ABSTOL=1.0e-10;
  int M;

  int LWORK;
  dcomplex *Z;
  dcomplex *A0;
  dcomplex *WORK;
  double *RWORK;
  int *IFAIL, INFO;
  int LDZ=N;
  int i,j;
  int *IWORK;

  A0=(dcomplex*)malloc(sizeof(dcomplex)*N*N);

  for (i=1;i<=N;i++) {
  for (j=1;j<=N;j++) {
  A0[(j-1)*N+i-1] = A[i][j];
  }
  }

  LWORK=3*N; 
  WORK=(dcomplex*)malloc(sizeof(dcomplex)*LWORK);
  Z=(dcomplex*)malloc(sizeof(dcomplex)*N*N);
  RWORK=(double*)malloc(sizeof(double)*7*N);
  IWORK=(int*)malloc(sizeof(int)*5*N);
  IFAIL=(int*)malloc(sizeof(int)*N);


  F77_NAME(zheevx,ZHEEVX)( JOBZ, RANGE, UPLO, &N, A0, &LDA, &VL, &VU, &IL, &IU,
  &ABSTOL, &M, W, Z, &LDZ, WORK, &LWORK, RWORK, IWORK,
  IFAIL, &INFO );

  if (INFO!=0) {
  printf("************************************************************\n");
  printf("  EigenBand_lapack: zheevx_()=%d\n",INFO);
  printf("************************************************************\n");
  exit(10);
  }

  for (i=1;i<=N;i++) {
  for (j=1;j<=N;j++) {
  A[i][j].r = Z[(j-1)*N+i-1].r;
  A[i][j].i = Z[(j-1)*N+i-1].i;
  }
  }

  for (i=N;i>=1;i--) {
  W[i] =W[i-1];
  }

  for (i=1; i<=N; i++) {
  printf("i=%2d W=%15.12f\n",i,W[i]);
  }

  free(A0); 
  free(WORK);
  free(Z); 
  free(RWORK); 
  free(IWORK); 
  free(IFAIL); 
  }
*/




void EigenBand_lapack(dcomplex **A, double *W, int N)
{
  static char *JOBZ="V";
  static char *UPLO="L";
  int LWORK;
  dcomplex *A0;
  dcomplex *WORK;
  double *RWORK;
  int INFO;
  int i,j;

  A0=(dcomplex*)malloc(sizeof(dcomplex)*N*N);

  for (i=1;i<=N;i++) {
    for (j=1;j<=N;j++) {
      A0[(j-1)*N+i-1] = A[i][j];
    }
  }

  LWORK=3*N; 
  WORK=(dcomplex*)malloc(sizeof(dcomplex)*LWORK);
  RWORK=(double*)malloc(sizeof(double)*(3*N-2));
  F77_NAME(zheev,ZHEEV)(JOBZ,UPLO, &N, A0, &N, W, WORK, &LWORK, RWORK, &INFO  );

  if (INFO!=0) {
    printf("************************************************************\n");
    printf("  EigenBand_lapack: cheev_()=%d\n",INFO);
    printf("************************************************************\n");
    exit(10);
  }

  for (i=1;i<=N;i++) {
    for (j=1;j<=N;j++) {
      A[i][j].r = A0[(j-1)*N+i-1].r;
      A[i][j].i = A0[(j-1)*N+i-1].i;
    }
  }

  for (i=N;i>=1;i--) {
    W[i] =W[i-1];
  }

  free(A0); free(RWORK); free(WORK); 
}

void dtime(double *t)
{
  /* real time */
  struct timeval timev;
  gettimeofday(&timev, NULL);
  *t = timev.tv_sec + (double)timev.tv_usec*1e-6;

  /* user time + system time */
  /*
    float tarray[2];
    clock_t times(), wall;
    struct tms tbuf;
    wall = times(&tbuf);
    tarray[0] = (float) (tbuf.tms_utime / (float)CLOCKS_PER_SEC);
    tarray[1] = (float) (tbuf.tms_stime / (float)CLOCKS_PER_SEC);
    *t = (double) (tarray[0]+tarray[1]);
    printf("dtime: %lf\n",*t);
  */
}

