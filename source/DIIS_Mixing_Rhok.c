/**********************************************************************
  DIIS_Mixing_Rhok.c:

     DIIS_Mixing_Rhok.c is a subroutine to achieve self-consistent field
     using the direct inversion in the iterative subspace in k-space.

  Log of DIIS_Mixing_Rhok.c:

     3/Jan/2005  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "openmx_common.h"
#include "tran_prototypes.h"
#include "tran_variables.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

#ifdef noomp
#include "mimic_omp.h"
#else
#include <omp.h>
#endif

#define  measure_time   0
#define  maxima_step  1000000.0


static void Inverse(int n, double **a, double **ia);
static void Complex_Inverse(int n, double **a, double **ia, double **b, double **ib);

static void DIIS_Mixing_Rhok0(int SCF_iter,
			      double Mix_wgt,
			      double *****ReRhok,
			      double *****ImRhok,
			      double ****ReBestRhok,
			      double ****ImBestRhok,
			      double *****Residual_ReRhok,
			      double *****Residual_ImRhok,
			      double ***ReV1,
			      double ***ImV1,
			      double ***ReV2,
			      double ***ImV2,
			      double ***ReRhoAtomk,
			      double ***ImRhoAtomk);

static void DIIS_Mixing_Rhok_Normal(int SCF_iter,
			      double Mix_wgt,
			      double *****ReRhok,
			      double *****ImRhok,
			      double ****ReBestRhok,
			      double ****ImBestRhok,
			      double *****Residual_ReRhok,
			      double *****Residual_ImRhok,
			      double ***ReV1,
			      double ***ImV1,
			      double ***ReV2,
			      double ***ImV2,
			      double ***ReRhoAtomk,
			      double ***ImRhoAtomk);

static void DIIS_Mixing_Rhok_NEGF(int SCF_iter,
				  double Mix_wgt,
				  double *****ReRhok,
				  double *****ImRhok,
				  double ****ReBestRhok,
				  double ****ImBestRhok,
				  double *****Residual_ReRhok,
				  double *****Residual_ImRhok,
				  double ***ReV1,
				  double ***ImV1,
				  double ***ReV2,
				  double ***ImV2,
				  double ***ReRhoAtomk,
				  double ***ImRhoAtomk);
	


static void Broyden1(int SCF_iter,
		     double Mix_wgt,
		     double *****ReRhok,
		     double *****ImRhok,
		     double ****ReBestRhok,
		     double ****ImBestRhok,
		     double *****Residual_ReRhok,
		     double *****Residual_ImRhok,
		     double ***ReV1,
		     double ***ImV1,
		     double ***ReV2,
		     double ***ImV2,
		     double ***ReRhoAtomk,
		     double ***ImRhoAtomk);

static void Broyden2(int SCF_iter,
		     double Mix_wgt,
		     double *****ReRhok,
		     double *****ImRhok,
		     double ****ReBestRhok,
		     double ****ImBestRhok,
		     double *****Residual_ReRhok,
		     double *****Residual_ImRhok,
		     double ***ReV1,
		     double ***ImV1,
		     double ***ReV2,
		     double ***ImV2,
		     double ***ReRhoAtomk,
		     double ***ImRhoAtomk);



#ifdef MIN
#undef MIN
#endif
#define MIN(a,b) ((a)<(b))?  (a):(b)


/*
static double Re_DiaG[3][130][130][130];
static double Im_DiaG[3][130][130][130];
*/



void DIIS_Mixing_Rhok(int SCF_iter,
                      double Mix_wgt,
                      double *****ReRhok,
                      double *****ImRhok,
                      double ****ReBestRhok,
                      double ****ImBestRhok,
                      double *****Residual_ReRhok,
                      double *****Residual_ImRhok,
                      double ***ReV1,
                      double ***ImV1,
                      double ***ReV2,
                      double ***ImV2,
                      double ***ReRhoAtomk,
                      double ***ImRhoAtomk)
{

  if (Solver!=4 || TRAN_Poisson_flag==4){

    DIIS_Mixing_Rhok_Normal(
			    SCF_iter,
			    Mix_wgt,
			    ReRhok,
			    ImRhok,
			    ReBestRhok,
			    ImBestRhok,
			    Residual_ReRhok,
			    Residual_ImRhok,
			    ReV1,
			    ImV1,
			    ReV2,
			    ImV2,
			    ReRhoAtomk,
			    ImRhoAtomk);
  }

  else if (Solver==4){

    DIIS_Mixing_Rhok_NEGF(
			  SCF_iter,
			  Mix_wgt,
			  ReRhok,
			  ImRhok,
			  ReBestRhok,
			  ImBestRhok,
			  Residual_ReRhok,
			  Residual_ImRhok,
			  ReV1,
			  ImV1,
			  ReV2,
			  ImV2,
			  ReRhoAtomk,
			  ImRhoAtomk);
  }

}



void DIIS_Mixing_Rhok_Normal(int SCF_iter,
			     double Mix_wgt,
			     double *****ReRhok,
			     double *****ImRhok,
			     double ****ReBestRhok,
			     double ****ImBestRhok,
			     double *****Residual_ReRhok,
			     double *****Residual_ImRhok,
			     double ***ReV1,
			     double ***ImV1,
			     double ***ReV2,
			     double ***ImV2,
			     double ***ReRhoAtomk,
			     double ***ImRhoAtomk)
{
  static int firsttime=1;
  int spin,Mc_AN,Gc_AN,wan1,TNO1,h_AN,Gh_AN;
  int wan2,TNO2,i,j,k,NumMix,NumSlide;
  int SCFi,SCFj,tno1,tno0,Cwan,Hwan,k1,k2,k3,kk2; 
  int pSCF_iter,size_Re_OptRhok,size_Kerker_weight;
  int spinmax,MN,flag_nan,refiter,po3,po4,refiter2;
  int *One2SCFi,*One2SCFj;
  double time1,time2,time3,time4,time5; 
  double Stime1,Etime1,scale;
  double *alden,*alden0,*gradF;
  double **A,**My_A,**IA,*My_A_thrds;
  double tmp0,itmp0,tmp1,im1,im2,re1,re2;
  double OptRNorm,My_OptRNorm,sum0,sum1;
  double dum1,dum2,bunsi,bunbo,sum,coef_OptRho;
  double Gx,Gy,Gz,G2,G12,G22,G32,F,F0;
  double Min_Weight;
  double Max_Weight;
  double G0,G02,G02p,weight,wgt0,wgt1;
  double sk1,sk2,sk3;
  double ****Re_OptRhok;
  double ****Im_OptRhok;
  double ***Kerker_weight;
  int numprocs,myid;
  char nanchar[300];
  /* for OpenMP */
  int OMPID,Nthrds,Nprocs,Nloop,Nthrds0;
  int knum,knum_full;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /* find an optimum G0 */

  G12 = rtv[1][1]*rtv[1][1] + rtv[1][2]*rtv[1][2] + rtv[1][3]*rtv[1][3]; 
  G22 = rtv[2][1]*rtv[2][1] + rtv[2][2]*rtv[2][2] + rtv[2][3]*rtv[2][3]; 
  G32 = rtv[3][1]*rtv[3][1] + rtv[3][2]*rtv[3][2] + rtv[3][3]*rtv[3][3]; 

  if (G12<G22) G0 = G12;
  else         G0 = G22;
  if (G32<G0)  G0 = G32;

  G0 = sqrt(G0);
  G02 = Kerker_factor*Kerker_factor*G0*G0;
  G02p = (0.01*Kerker_factor*G0)*(0.01*Kerker_factor*G0);

  if      (SpinP_switch==0)  spinmax = 1;
  else if (SpinP_switch==1)  spinmax = 2;
  else if (SpinP_switch==3)  spinmax = 3;

  /****************************************************
    allocation of arrays:

    double alden[List_YOUSO[38]];
    double A[List_YOUSO[38]][List_YOUSO[38]];
    double IA[List_YOUSO[38]][List_YOUSO[38]];
    double My_A[List_YOUSO[38]][List_YOUSO[38]];
    double Re_OptRhok[spinmax]
                            [My_NGrid2_Poisson]
                            [Ngrid1]
                            [Ngrid3];
    double Im_OptRhok[spinmax]
                            [My_NGrid2_Poisson]
                            [Ngrid1]
                            [Ngrid3];
  ****************************************************/
  
  alden = (double*)malloc(sizeof(double)*List_YOUSO[38]);
  alden0 = (double*)malloc(sizeof(double)*List_YOUSO[38]);
  gradF = (double*)malloc(sizeof(double)*List_YOUSO[38]);

  A = (double**)malloc(sizeof(double*)*List_YOUSO[38]);
  for (i=0; i<List_YOUSO[38]; i++){
    A[i] = (double*)malloc(sizeof(double)*List_YOUSO[38]);
  }

  My_A = (double**)malloc(sizeof(double*)*List_YOUSO[38]);
  for (i=0; i<List_YOUSO[38]; i++){
    My_A[i] = (double*)malloc(sizeof(double)*List_YOUSO[38]);
  }

  IA = (double**)malloc(sizeof(double*)*List_YOUSO[38]);
  for (i=0; i<List_YOUSO[38]; i++){
    IA[i] = (double*)malloc(sizeof(double)*List_YOUSO[38]);
  }

  /* Re_OptRhok and Im_OptRhok */

  size_Re_OptRhok = spinmax*My_NGrid2_Poisson*Ngrid1*Ngrid3;

  Re_OptRhok = (double****)malloc(sizeof(double***)*spinmax); 
  for (spin=0; spin<spinmax; spin++){
    Re_OptRhok[spin] = (double***)malloc(sizeof(double**)*My_NGrid2_Poisson); 
    for (i=0; i<My_NGrid2_Poisson; i++){
      Re_OptRhok[spin][i] = (double**)malloc(sizeof(double*)*Ngrid1); 
      for (j=0; j<Ngrid1; j++){
	Re_OptRhok[spin][i][j] = (double*)malloc(sizeof(double)*Ngrid3); 
	for (k=0; k<Ngrid3; k++) Re_OptRhok[spin][i][j][k] = 0.0;
      }
    }
  }

  Im_OptRhok = (double****)malloc(sizeof(double***)*spinmax); 
  for (spin=0; spin<spinmax; spin++){
    Im_OptRhok[spin] = (double***)malloc(sizeof(double**)*My_NGrid2_Poisson); 
    for (i=0; i<My_NGrid2_Poisson; i++){
      Im_OptRhok[spin][i] = (double**)malloc(sizeof(double*)*Ngrid1); 
      for (j=0; j<Ngrid1; j++){
	Im_OptRhok[spin][i][j] = (double*)malloc(sizeof(double)*Ngrid3); 
	for (k=0; k<Ngrid3; k++) Im_OptRhok[spin][i][j][k] = 0.0;
      }
    }
  }

  /***********************************
            set Kerker_weight 
  ************************************/

  size_Kerker_weight = My_NGrid2_Poisson*Ngrid1*Ngrid3;
  knum = My_NGrid2_Poisson*Ngrid1*Ngrid3;
  knum_full = spinmax*My_NGrid2_Poisson*Ngrid1*Ngrid3;

  Kerker_weight = (double***)malloc(sizeof(double**)*My_NGrid2_Poisson); 
  for (i=0; i<My_NGrid2_Poisson; i++){
    Kerker_weight[i] = (double**)malloc(sizeof(double*)*Ngrid1); 
    for (j=0; j<Ngrid1; j++){
      Kerker_weight[i][j] = (double*)malloc(sizeof(double)*Ngrid3); 
    }
  }

  if (measure_time==1) dtime(&Stime1);

#pragma omp parallel shared(myid,Start_Grid2,knum,G02,G02p,Kerker_weight,rtv,Ngrid3,Ngrid2,Ngrid1) private(OMPID,Nthrds,Nprocs,k2,kk2,sk2,k,k1,k3,sk1,sk3,Gx,Gy,Gz,G2,weight)
  {

    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    /* one-dimensionalized loop */

    for (k=OMPID*knum/Nthrds; k<(OMPID+1)*knum/Nthrds; k++){

      /* get k2, k1, and k3 */
 
      k2 = k/(Ngrid1*Ngrid3);
      k1 = (k - k2*(Ngrid1*Ngrid3))/Ngrid3;
      k3 = k - k2*(Ngrid1*Ngrid3) - k1*Ngrid3;

      kk2 = k2 + Start_Grid2[myid];

      if (kk2<Ngrid2/2) sk2 = (double)kk2;
      else              sk2 = (double)(kk2 - Ngrid2);

      if (k1<Ngrid1/2) sk1 = (double)k1;
      else             sk1 = (double)(k1 - Ngrid1);

      if (k3<Ngrid3/2) sk3 = (double)k3;
      else             sk3 = (double)(k3 - Ngrid3);

      Gx = sk1*rtv[1][1] + sk2*rtv[2][1] + sk3*rtv[3][1];
      Gy = sk1*rtv[1][2] + sk2*rtv[2][2] + sk3*rtv[3][2]; 
      Gz = sk1*rtv[1][3] + sk2*rtv[2][3] + sk3*rtv[3][3];

      G2 = Gx*Gx + Gy*Gy + Gz*Gz;

      if (k1==0 && kk2==0 && k3==0)  weight = 1.0;
      else                           weight = (G2 + G02)/(G2 + G02p);

      Kerker_weight[k2][k1][k3] = weight;

    } /* k */
  } /* #pragma omp parallel */

  if (measure_time==1){ 
    dtime(&Etime1);
    time1 = Etime1 - Stime1;      
  }

  /****************************************************
     start calc.
  ****************************************************/
 
  if (SCF_iter==1){
    /* memory calc. done, only 1st iteration */
    if (firsttime) {
      PrintMemory("DIIS_Mixing_Rhok: Re_OptRhok",sizeof(double)*size_Re_OptRhok,NULL);
      PrintMemory("DIIS_Mixing_Rhok: Im_OptRhok",sizeof(double)*size_Re_OptRhok,NULL);
      PrintMemory("DIIS_Mixing_Rhok: Kerker_weight",sizeof(double)*size_Kerker_weight,NULL);
      firsttime=0;
    }

    Kerker_Mixing_Rhok(0,1.00,ReRhok,ImRhok,ReBestRhok,ImBestRhok,
                       Residual_ReRhok,Residual_ImRhok,
                       ReV1,ImV1,ReV2,ImV2,ReRhoAtomk,ImRhoAtomk);
  } 

  else if (SCF_iter==2){

    if (measure_time==1) dtime(&Stime1);

    /* construct residual rho */

#pragma omp parallel shared(Residual_ReRhok,Residual_ImRhok,ImRhok,ReRhok,knum_full,My_NGrid2_Poisson,Ngrid1,Ngrid3) private(OMPID,Nthrds,Nprocs,k,spin,k2,k1,k3,tmp0,tmp1)
    {

      /* get info. on OpenMP */ 

      OMPID = omp_get_thread_num();
      Nthrds = omp_get_num_threads();
      Nprocs = omp_get_num_procs();

      /* one-dimensionalized loop */

      for (k=OMPID*knum_full/Nthrds; k<(OMPID+1)*knum_full/Nthrds; k++){

	/* get spin, k2, k1, and k3 */

	spin = k/(My_NGrid2_Poisson*Ngrid1*Ngrid3);
	k2 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3))/(Ngrid1*Ngrid3);
	k1 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3))/Ngrid3;  
	k3 = k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3) - k1*Ngrid3;

	tmp0 = ReRhok[0][spin][k2][k1][k3] - ReRhok[1][spin][k2][k1][k3];
	tmp1 = ImRhok[0][spin][k2][k1][k3] - ImRhok[1][spin][k2][k1][k3];
	Residual_ReRhok[2][spin][k2][k1][k3] = tmp0;
	Residual_ImRhok[2][spin][k2][k1][k3] = tmp1;

      } /* k */
    } /* #pragma omp parallel */

    if (measure_time==1){ 
      dtime(&Etime1);
      time2 = Etime1 - Stime1;      
    }

    /* rho1 to rho2 */

    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
        for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){
	    ReRhok[2][spin][k2][k1][k3] = ReRhok[1][spin][k2][k1][k3];
	    ImRhok[2][spin][k2][k1][k3] = ImRhok[1][spin][k2][k1][k3];
	  }
	}
      }
    }

    Kerker_Mixing_Rhok(0,Mixing_weight,ReRhok,ImRhok,ReBestRhok,ImBestRhok,
                       Residual_ReRhok,Residual_ImRhok,
                       ReV1,ImV1,ReV2,ImV2,ReRhoAtomk,ImRhoAtomk);
  } 

  else {

    /****************************************************
                   construct residual rho1
    ****************************************************/

    if (measure_time==1) dtime(&Stime1);

#pragma omp parallel shared(Residual_ReRhok,Residual_ImRhok,ImRhok,ReRhok,knum_full,My_NGrid2_Poisson,Ngrid1,Ngrid3) private(OMPID,Nthrds,Nprocs,k,spin,k2,k1,k3,tmp0,tmp1)
    {
      /* get info. on OpenMP */ 

      OMPID = omp_get_thread_num();
      Nthrds = omp_get_num_threads();
      Nprocs = omp_get_num_procs();

      /* one-dimensionalized loop */

      for (k=OMPID*knum_full/Nthrds; k<(OMPID+1)*knum_full/Nthrds; k++){

	/* get spin, k2, k1, and k3 */

	spin = k/(My_NGrid2_Poisson*Ngrid1*Ngrid3);
	k2 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3))/(Ngrid1*Ngrid3);
	k1 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3))/Ngrid3;  
	k3 = k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3) - k1*Ngrid3;

	tmp0 = ReRhok[0][spin][k2][k1][k3] - ReRhok[1][spin][k2][k1][k3];
	tmp1 = ImRhok[0][spin][k2][k1][k3] - ImRhok[1][spin][k2][k1][k3];
	Residual_ReRhok[1][spin][k2][k1][k3] = tmp0;
	Residual_ImRhok[1][spin][k2][k1][k3] = tmp1;

      } /* k */
    } /* #pragma omp parallel */

    if (measure_time==1){ 
      dtime(&Etime1);
      time3 = Etime1 - Stime1;      
    }

    if ((SCF_iter-1)<Num_Mixing_pDM){
      NumMix   = SCF_iter - 1;
      NumSlide = NumMix + 1;
    }
    else{
      NumMix   = Num_Mixing_pDM;
      NumSlide = NumMix;
    }

    /****************************************************
                    alpha from residual rho
    ****************************************************/

    if (measure_time==1) dtime(&Stime1);

    /* allocation of arrrays */

    One2SCFi = (int*)malloc(sizeof(int)*(NumMix*(NumMix+1))/2); 
    One2SCFj = (int*)malloc(sizeof(int)*(NumMix*(NumMix+1))/2); 

    i = 0;
    for (SCFi=1; SCFi<=NumMix; SCFi++){
      for (SCFj=SCFi; SCFj<=NumMix; SCFj++){
        One2SCFi[i] = SCFi;
        One2SCFj[i] = SCFj;
        i++; 
      }
    }
    
    /* calculation of the norm matrix for the residual vectors */

#pragma omp parallel shared(NumMix,One2SCFi,One2SCFj,My_A,spinmax,My_NGrid2_Poisson,Ngrid1,Ngrid3,Kerker_weight,Residual_ReRhok,Residual_ImRhok) private(OMPID,Nthrds,Nprocs,i,SCFi,SCFj,spin,k2,k1,k3,weight,re1,im1,re2,im2,tmp0)
    {

      /* get info. on OpenMP */ 

      OMPID = omp_get_thread_num();
      Nthrds = omp_get_num_threads();
      Nprocs = omp_get_num_procs();

      for ( i=OMPID; i<(NumMix*(NumMix+1))/2; i+=Nthrds ){

	SCFi = One2SCFi[i]; 
	SCFj = One2SCFj[i]; 

        tmp0 = 0.0; 

	for (spin=0; spin<spinmax; spin++){
	  for (k2=0; k2<My_NGrid2_Poisson; k2++){
	    for (k1=0; k1<Ngrid1; k1++){
	      for (k3=0; k3<Ngrid3; k3++){

		weight = Kerker_weight[k2][k1][k3]; 

		re1 = Residual_ReRhok[SCFi][spin][k2][k1][k3];
		im1 = Residual_ImRhok[SCFi][spin][k2][k1][k3];

		re2 = Residual_ReRhok[SCFj][spin][k2][k1][k3];
		im2 = Residual_ImRhok[SCFj][spin][k2][k1][k3];

		tmp0 += (re1*re2 + im1*im2)*weight;
	      }
	    }
	  }
	}

	My_A[SCFi][SCFj] = tmp0;

      }
    } /* pragma omp parallel */

    /* freeing of arrrays */
    free(One2SCFj);
    free(One2SCFi);

    if (measure_time==1){ 
      dtime(&Etime1);
      time4 = Etime1 - Stime1;      
    }

    /* MPI My_A */
    for (SCFi=1; SCFi<=NumMix; SCFi++){
      for (SCFj=SCFi; SCFj<=NumMix; SCFj++){
        MPI_Allreduce(&My_A[SCFi][SCFj], &A[SCFi-1][SCFj-1],
                      1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
        A[SCFj-1][SCFi-1] = A[SCFi-1][SCFj-1];
      }
    }

    /* store NormRD */

    for (i=4; 1<=i; i--){
      NormRD[i] = NormRD[i-1];
      History_Uele[i] = History_Uele[i-1];
    }
    NormRD[0] = A[0][0]/(double)Ngrid1/(double)Ngrid2/(double)Ngrid3/(double)spinmax;
    History_Uele[0] = Uele;

    /* solve the linear equation */

    for (SCFi=1; SCFi<=NumMix; SCFi++){
       A[SCFi-1][NumMix] = -1.0;
       A[NumMix][SCFi-1] = -1.0;
    }
    A[NumMix][NumMix] = 0.0;

    Inverse(NumMix,A,IA);

    for (SCFi=1; SCFi<=NumMix; SCFi++){
       alden[SCFi] = -IA[SCFi-1][NumMix];
    }

    /*************************************
      check "nan", "NaN", "inf" or "Inf"
    *************************************/

    flag_nan = 0;
    for (SCFi=1; SCFi<=NumMix; SCFi++){

      sprintf(nanchar,"%8.4f",alden[SCFi]);
      if (strstr(nanchar,"nan")!=NULL || strstr(nanchar,"NaN")!=NULL 
       || strstr(nanchar,"inf")!=NULL || strstr(nanchar,"Inf")!=NULL){

        flag_nan = 1;
      }
    }

    if (flag_nan==1){
      for (SCFi=1; SCFi<=NumMix; SCFi++){
        alden[SCFi] = 0.0;
      }
      alden[1] = 0.1;
      alden[2] = 0.9;
    }

    /****************************************************
        refinement of alpha using an iterative method
    ****************************************************/

    if (0){

    po3 = 0;
    refiter = 0;

    do {
    
      F = 0.0;
      for (SCFi=0; SCFi<NumMix; SCFi++){
	for (SCFj=0; SCFj<NumMix; SCFj++){
	  F += alden[SCFi+1]*alden[SCFj+1]*A[SCFi][SCFj];
	}
      }

      /*
      printf("A refiter=%2d F=%15.12f\n",refiter,F);  
      */

      sum0 = 0.0; 
      for (SCFi=0; SCFi<NumMix; SCFi++){
	sum0 += 2.0*alden[SCFi+1]*A[SCFi][0];
      }

      for (SCFi=1; SCFi<NumMix; SCFi++){

        sum1 = 0.0;
	for (SCFj=0; SCFj<NumMix; SCFj++){
          sum1 += 2.0*alden[SCFj+1]*A[SCFj][SCFi];
	}

	gradF[SCFi+1] = -sum0 + sum1;
      }

      scale = 0.01;
      refiter2 = 1;
      po4 = 0;

      do {      

	for (SCFi=1; SCFi<NumMix; SCFi++){
	  alden0[SCFi+1] = alden[SCFi+1] - scale*gradF[SCFi+1];
	}

	sum = 0.0;
	for (SCFi=1; SCFi<NumMix; SCFi++){
	  sum += alden0[SCFi+1]; 
	}
	alden0[1] = 1.0 - sum;

	F0 = 0.0;
	for (SCFi=0; SCFi<NumMix; SCFi++){
	  for (SCFj=0; SCFj<NumMix; SCFj++){
	    F0 += alden0[SCFi+1]*alden0[SCFj+1]*A[SCFi][SCFj];
	  }
	}

	/*
        printf("B refiter=%2d F0=%15.12f scale=%15.12f\n",refiter,F,scale);  
	*/

	if (F0<=F) po4 = 1;
	else       scale /= 2.0; 

	refiter2++;

      } while (po4==0 && refiter2<20);

      for (SCFi=0; SCFi<NumMix; SCFi++){
        alden[SCFi+1] = alden0[SCFi+1];
      }

      refiter++;

    } while (refiter<10);

    }

    /****************************************************
              calculate optimized residual rho
    ****************************************************/

    if (measure_time==1) dtime(&Stime1);

    for (pSCF_iter=1; pSCF_iter<=NumMix; pSCF_iter++){

      tmp0 = alden[pSCF_iter]; 

#pragma omp parallel shared(pSCF_iter,tmp0,Im_OptRhok,Re_OptRhok,Residual_ReRhok,Residual_ImRhok,knum_full,My_NGrid2_Poisson,Ngrid1,Ngrid3) private(OMPID,Nthrds,Nprocs,k,spin,k2,k1,k3)
      {
	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	/* one-dimensionalized loop */

	for (k=OMPID*knum_full/Nthrds; k<(OMPID+1)*knum_full/Nthrds; k++){

	  /* get spin, k2, k1, and k3 */

	  spin = k/(My_NGrid2_Poisson*Ngrid1*Ngrid3);
	  k2 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3))/(Ngrid1*Ngrid3);
	  k1 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3))/Ngrid3;  
	  k3 = k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3) - k1*Ngrid3;

	  /* Pulay mixing */

	  Re_OptRhok[spin][k2][k1][k3] += tmp0*Residual_ReRhok[pSCF_iter][spin][k2][k1][k3];
	  Im_OptRhok[spin][k2][k1][k3] += tmp0*Residual_ImRhok[pSCF_iter][spin][k2][k1][k3];

	} /* k */
      } /* #pragma omp parallel */
    } /* pSCF_iter */

    /* the norm of the optimized residual rho */

    /*
    My_OptRNorm = 0.0;
    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
	for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){

	    weight = Kerker_weight[k2][k1][k3]; 

	    re1 = Re_OptRhok[spin][k2][k1][k3];
	    im1 = Im_OptRhok[spin][k2][k1][k3];

	    My_OptRNorm += (re1*re1 + im1*im1)*weight;
	  }
	}
      }
    }

    MPI_Allreduce(&My_OptRNorm, &OptRNorm, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
    OptRNorm = OptRNorm/(double)Ngrid1/(double)Ngrid2/(double)Ngrid3/(double)spinmax;
    printf("OptRNorm=%15.12f\n",sqrt(OptRNorm)); 
    */

    if      (1.0e-4<=NormRD[0])                        coef_OptRho = 0.5;
    else if (1.0e-6<=NormRD[0]  && NormRD[0]<1.0e-4)   coef_OptRho = 0.6;
    else if (1.0e-8<=NormRD[0]  && NormRD[0]<1.0e-6)   coef_OptRho = 0.7;
    else if (1.0e-10<=NormRD[0] && NormRD[0]<1.0e-8)   coef_OptRho = 0.8;
    else if (1.0e-12<=NormRD[0] && NormRD[0]<1.0e-10)  coef_OptRho = 0.9;
    else                                               coef_OptRho = 1.0;

    /****************************************************
            store the best input charge density 
    ****************************************************/

    if (NormRD[0]<BestNormRD ){

      BestNormRD = NormRD[0];

      for (spin=0; spin<spinmax; spin++){
        for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){
              ReBestRhok[spin][k2][k1][k3] = ReRhok[1][spin][k2][k1][k3];
              ImBestRhok[spin][k2][k1][k3] = ImRhok[1][spin][k2][k1][k3];
	    }
	  }
        }
      }
    }

    /****************************************************
                        mixing of rho 
    ****************************************************/

    /****************************************************
       Pulay mixing
    ****************************************************/

    if ( SCF_iter%EveryPulay_SCF==0 ) {

      /* reduce Mixing_weight so that charge sloshing can be avoided at the next SCF */

      Mixing_weight = 0.1*Max_Mixing_weight;
      if (Mixing_weight<Min_Mixing_weight) Mixing_weight = Min_Mixing_weight;

      /* initial ReRhok and ImRhok */

      for (spin=0; spin<spinmax; spin++){
	for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){
	      ReRhok[0][spin][k2][k1][k3] = 0.0;
	      ImRhok[0][spin][k2][k1][k3] = 0.0;
	    }
	  }
	}
      }

      /* Pulay mixing */

      for (pSCF_iter=1; pSCF_iter<=NumMix; pSCF_iter++){

	tmp0 =  alden[pSCF_iter];

#pragma omp parallel shared(pSCF_iter,tmp0,ReRhok,ImRhok,knum_full,My_NGrid2_Poisson,Ngrid1,Ngrid3) private(OMPID,Nthrds,Nprocs,k,spin,k2,k1,k3)
	{
	  /* get info. on OpenMP */ 

	  OMPID = omp_get_thread_num();
	  Nthrds = omp_get_num_threads();
	  Nprocs = omp_get_num_procs();

	  /* one-dimensionalized loop */

	  for (k=OMPID*knum_full/Nthrds; k<(OMPID+1)*knum_full/Nthrds; k++){

	    /* get spin, k2, k1, and k3 */

	    spin = k/(My_NGrid2_Poisson*Ngrid1*Ngrid3);
	    k2 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3))/(Ngrid1*Ngrid3);
	    k1 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3))/Ngrid3;  
	    k3 = k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3) - k1*Ngrid3;

	    ReRhok[0][spin][k2][k1][k3] += tmp0*ReRhok[pSCF_iter][spin][k2][k1][k3];
	    ImRhok[0][spin][k2][k1][k3] += tmp0*ImRhok[pSCF_iter][spin][k2][k1][k3];

	  } /* k */
	} /* #pragma omp parallel */
      } /* pSCF_iter */

      /* Correction by optimized residual rho */

#pragma omp parallel shared(coef_OptRho,Im_OptRhok,Re_OptRhok,ReRhok,ImRhok,knum_full,My_NGrid2_Poisson,Ngrid1,Ngrid3) private(OMPID,Nthrds,Nprocs,k,spin,k2,k1,k3)
      {
	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	/* one-dimensionalized loop */

	for (k=OMPID*knum_full/Nthrds; k<(OMPID+1)*knum_full/Nthrds; k++){

	  /* get spin, k2, k1, and k3 */

	  spin = k/(My_NGrid2_Poisson*Ngrid1*Ngrid3);
	  k2 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3))/(Ngrid1*Ngrid3);
	  k1 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3))/Ngrid3;  
	  k3 = k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3) - k1*Ngrid3;

	  ReRhok[0][spin][k2][k1][k3] += coef_OptRho*Re_OptRhok[spin][k2][k1][k3];
	  ImRhok[0][spin][k2][k1][k3] += coef_OptRho*Im_OptRhok[spin][k2][k1][k3];

	  /* correction to large changing components */

          tmp0 = ReRhok[0][spin][k2][k1][k3] - ReRhok[1][spin][k2][k1][k3];  
          tmp1 = ImRhok[0][spin][k2][k1][k3] - ImRhok[1][spin][k2][k1][k3];  

          if ( maxima_step<(fabs(tmp0)+fabs(tmp1)) ){
            ReRhok[0][spin][k2][k1][k3] = sgn(tmp0)*maxima_step + ReRhok[1][spin][k2][k1][k3]; 
            ImRhok[0][spin][k2][k1][k3] = sgn(tmp1)*maxima_step + ImRhok[1][spin][k2][k1][k3]; 
          }

	} /* k */
      } /* #pragma omp parallel */

      if (measure_time==1){ 
        dtime(&Etime1);
        time5 = Etime1 - Stime1;      
      }

    }

    /****************************************************
       Kerker mixing
    ****************************************************/

    else {

      /* find an optimum mixing weight */

      Min_Weight = Min_Mixing_weight;
      Max_Weight = Max_Mixing_weight;

      if ((int)sgn(History_Uele[0]-History_Uele[1])
	  ==(int)sgn(History_Uele[1]-History_Uele[2])
	  && NormRD[0]<NormRD[1]){

	/* tmp0 = 1.6*Mixing_weight; */

	tmp0 = NormRD[1]/(largest(NormRD[1]-NormRD[0],10e-10))*Mixing_weight;

	if (tmp0<Max_Weight){
	  if (Min_Weight<tmp0){
	    Mixing_weight = tmp0;
	  }
	  else{ 
	    Mixing_weight = Min_Weight;
	  }
	}
	else{ 
	  Mixing_weight = Max_Weight;
	  SCF_RENZOKU++;  
	}
      }
   
      else if ((int)sgn(History_Uele[0]-History_Uele[1])
	       ==(int)sgn(History_Uele[1]-History_Uele[2])
	       && NormRD[1]<NormRD[0]){

	tmp0 = NormRD[1]/(largest(NormRD[1]+NormRD[0],10e-10))*Mixing_weight;

	/* tmp0 = Mixing_weight/1.6; */

	if (tmp0<Max_Weight){
	  if (Min_Weight<tmp0)
	    Mixing_weight = tmp0;
	  else 
	    Mixing_weight = Min_Weight;
	}
	else 
	  Mixing_weight = Max_Weight;

	SCF_RENZOKU = -1;  
      }

      else if ((int)sgn(History_Uele[0]-History_Uele[1])
	       !=(int)sgn(History_Uele[1]-History_Uele[2])
	       && NormRD[0]<NormRD[1]){

	/* tmp0 = Mixing_weight*1.2; */

	tmp0 = NormRD[1]/(largest(NormRD[1]-NormRD[0],10e-10))*Mixing_weight;

	if (tmp0<Max_Weight){
	  if (Min_Weight<tmp0)
	    Mixing_weight = tmp0;
	  else 
	    Mixing_weight = Min_Weight;
	}
	else{ 
	  Mixing_weight = Max_Weight;
	  SCF_RENZOKU++;
	}
      }

      else if ((int)sgn(History_Uele[0]-History_Uele[1])
	       !=(int)sgn(History_Uele[1]-History_Uele[2])
	       && NormRD[1]<NormRD[0]){

	/* tmp0 = Mixing_weight/2.0; */

	tmp0 = NormRD[1]/(largest(NormRD[1]+NormRD[0],10e-10))*Mixing_weight;

	if (tmp0<Max_Weight){
	  if (Min_Weight<tmp0)
	    Mixing_weight = tmp0;
	  else 
	    Mixing_weight = Min_Weight;
	}
	else 
	  Mixing_weight = Max_Weight;

	SCF_RENZOKU = -1;
      }

      Mix_wgt = Mixing_weight;

      /* Kerer mixing */

#pragma omp parallel shared(Mix_wgt,Kerker_weight,ReRhok,ImRhok,knum_full,My_NGrid2_Poisson,Ngrid1,Ngrid3) private(OMPID,Nthrds,Nprocs,k,spin,k2,k1,k3,weight,wgt0,wgt1)
      {
	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	/* one-dimensionalized loop */

	for (k=OMPID*knum_full/Nthrds; k<(OMPID+1)*knum_full/Nthrds; k++){

	  /* get spin, k2, k1, and k3 */

	  spin = k/(My_NGrid2_Poisson*Ngrid1*Ngrid3);
	  k2 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3))/(Ngrid1*Ngrid3);
	  k1 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3))/Ngrid3;  
	  k3 = k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3) - k1*Ngrid3;

	  weight = 1.0/Kerker_weight[k2][k1][k3];
	  wgt0  = Mix_wgt*weight;
	  wgt1 =  1.0 - wgt0;

	  ReRhok[0][spin][k2][k1][k3] = wgt0*ReRhok[0][spin][k2][k1][k3]
	                              + wgt1*ReRhok[1][spin][k2][k1][k3];

	  ImRhok[0][spin][k2][k1][k3] = wgt0*ImRhok[0][spin][k2][k1][k3]
	                              + wgt1*ImRhok[1][spin][k2][k1][k3];

	} /* k */
      } /* #pragma omp parallel */      
    } /* else */

    /****************************************************
                         shift of rho
    ****************************************************/

    for (pSCF_iter=(List_YOUSO[38]-1); 0<pSCF_iter; pSCF_iter--){
      for (spin=0; spin<spinmax; spin++){
	for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){
	      ReRhok[pSCF_iter][spin][k2][k1][k3] = ReRhok[pSCF_iter-1][spin][k2][k1][k3]; 
	      ImRhok[pSCF_iter][spin][k2][k1][k3] = ImRhok[pSCF_iter-1][spin][k2][k1][k3]; 
	    }
	  }
	}
      }
    }

    /****************************************************
                    shift of residual rho
    ****************************************************/

    for (pSCF_iter=(List_YOUSO[38]-1); 0<pSCF_iter; pSCF_iter--){
      for (spin=0; spin<spinmax; spin++){
	for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){
	      Residual_ReRhok[pSCF_iter][spin][k2][k1][k3] = Residual_ReRhok[pSCF_iter-1][spin][k2][k1][k3];
	      Residual_ImRhok[pSCF_iter][spin][k2][k1][k3] = Residual_ImRhok[pSCF_iter-1][spin][k2][k1][k3]; 
	    }
	  }
	}
      }
    }

  }

  /****************************************************
        find the charge density in real space 
  ****************************************************/

  tmp0 = 1.0/(double)Ngrid1/(double)Ngrid2/(double)Ngrid3;

  for (spin=0; spin<spinmax; spin++){

    for (k2=0; k2<My_NGrid2_Poisson; k2++){
      for (k1=0; k1<Ngrid1; k1++){
        for (k3=0; k3<Ngrid3; k3++){
          ReV2[k2][k1][k3] = ReRhok[0][spin][k2][k1][k3];
          ImV2[k2][k1][k3] = ImRhok[0][spin][k2][k1][k3];
        }
      }
    }

    if (spin==0 || spin==1){
      Get_Value_inReal(0,ReV2,ImV2,ReV1,ImV1,Density_Grid[spin],Density_Grid[spin]);

      for (MN=0; MN<My_NumGrid1; MN++){
        Density_Grid[spin][MN] = Density_Grid[spin][MN]*tmp0;
      }
    }

    else if (spin==2){
      Get_Value_inReal(1,ReV2,ImV2,ReV1,ImV1,Density_Grid[2],Density_Grid[3]);
      for (MN=0; MN<My_NumGrid1; MN++){
        Density_Grid[2][MN] = Density_Grid[2][MN]*tmp0;
        Density_Grid[3][MN] = Density_Grid[3][MN]*tmp0;
      }
    }

  }

  if (SpinP_switch==0){
    for (MN=0; MN<My_NumGrid1; MN++){
      Density_Grid[1][MN] = Density_Grid[0][MN];
    }
  }


  /****************************************************
      set ReV2 and ImV2 which are used in Poisson.c
  ****************************************************/

  if (SpinP_switch==0){

#pragma omp parallel shared(ImRhoAtomk,ReRhoAtomk,ImRhok,ReRhok,ImV2,ReV2,Ngrid1,Ngrid3,knum) private(OMPID,Nthrds,Nprocs,k,k2,k1,k3)
    {

      /* get info. on OpenMP */ 

      OMPID = omp_get_thread_num();
      Nthrds = omp_get_num_threads();
      Nprocs = omp_get_num_procs();

      /* one-dimensionalized loop */

      for (k=OMPID*knum/Nthrds; k<(OMPID+1)*knum/Nthrds; k++){

	/* get k2, k1, and k3 */
 
	k2 = k/(Ngrid1*Ngrid3);
	k1 = (k - k2*(Ngrid1*Ngrid3))/Ngrid3;
	k3 = k - k2*(Ngrid1*Ngrid3) - k1*Ngrid3;
	ReV2[k2][k1][k3] = 2.0*ReRhok[0][0][k2][k1][k3] - ReRhoAtomk[k2][k1][k3];
	ImV2[k2][k1][k3] = 2.0*ImRhok[0][0][k2][k1][k3] - ImRhoAtomk[k2][k1][k3];

      } /* k */
    } /* #pragma omp parallel */
  }
  
  else {

#pragma omp parallel shared(ImRhoAtomk,ReRhoAtomk,ImRhok,ReRhok,ImV2,ReV2,Ngrid1,Ngrid3,knum) private(OMPID,Nthrds,Nprocs,k,k2,k1,k3)
    {

      /* get info. on OpenMP */ 

      OMPID = omp_get_thread_num();
      Nthrds = omp_get_num_threads();
      Nprocs = omp_get_num_procs();

      /* one-dimensionalized loop */

      for (k=OMPID*knum/Nthrds; k<(OMPID+1)*knum/Nthrds; k++){

	/* get k2, k1, and k3 */
 
	k2 = k/(Ngrid1*Ngrid3);
	k1 = (k - k2*(Ngrid1*Ngrid3))/Ngrid3;
	k3 = k - k2*(Ngrid1*Ngrid3) - k1*Ngrid3;

	ReV2[k2][k1][k3] = ReRhok[0][0][k2][k1][k3] + ReRhok[0][1][k2][k1][k3] - ReRhoAtomk[k2][k1][k3];
	ImV2[k2][k1][k3] = ImRhok[0][0][k2][k1][k3] + ImRhok[0][1][k2][k1][k3] - ImRhoAtomk[k2][k1][k3];

      } /* k */
    } /* #pragma omp parallel */
  }

  if (measure_time==1){ 
    printf("time1=%12.6f time2=%12.6f time3=%12.6f time4=%12.6f time5=%12.6f\n",time1,time2,time3,time4,time5);
  }

  /****************************************************
    freeing of arrays:

    double alden[List_YOUSO[38]];
    double A[List_YOUSO[38]][List_YOUSO[38]];
    double IA[List_YOUSO[38]][List_YOUSO[38]];
    double My_A[List_YOUSO[38]][List_YOUSO[38]];
    double Re_OptRhok[spinmax]
                            [My_NGrid2_Poisson]
                            [Ngrid1]
                            [Ngrid3];
    double Im_OptRhok[spinmax]
                            [My_NGrid2_Poisson]
                            [Ngrid1]
                            [Ngrid3];
  ****************************************************/
  
  free(alden);
  free(alden0);
  free(gradF);

  for (i=0; i<List_YOUSO[38]; i++){
    free(A[i]);
  }
  free(A);

  for (i=0; i<List_YOUSO[38]; i++){
    free(My_A[i]);
  }
  free(My_A);

  for (i=0; i<List_YOUSO[38]; i++){
    free(IA[i]);
  }
  free(IA);

  /* Re_OptRhok and Im_OptRhok */

  for (spin=0; spin<spinmax; spin++){
    for (i=0; i<My_NGrid2_Poisson; i++){
      for (j=0; j<Ngrid1; j++){
	free(Re_OptRhok[spin][i][j]);
      }
      free(Re_OptRhok[spin][i]);
    }
    free(Re_OptRhok[spin]);
  }
  free(Re_OptRhok);

  for (spin=0; spin<spinmax; spin++){
    for (i=0; i<My_NGrid2_Poisson; i++){
      for (j=0; j<Ngrid1; j++){
	free(Im_OptRhok[spin][i][j]);
      }
      free(Im_OptRhok[spin][i]);
    }
    free(Im_OptRhok[spin]);
  }
  free(Im_OptRhok);

  /* Kerker_weight */

  for (i=0; i<My_NGrid2_Poisson; i++){
    for (j=0; j<Ngrid1; j++){
      free(Kerker_weight[i][j]);
    }
    free(Kerker_weight[i]);
  }
  free(Kerker_weight);

}



void DIIS_Mixing_Rhok_NEGF(int SCF_iter,
			   double Mix_wgt,
			   double *****ReRhok,
			   double *****ImRhok,
			   double ****ReBestRhok,
			   double ****ImBestRhok,
			   double *****Residual_ReRhok,
			   double *****Residual_ImRhok,
			   double ***ReV1,
			   double ***ImV1,
			   double ***ReV2,
			   double ***ImV2,
			   double ***ReRhoAtomk,
			   double ***ImRhoAtomk)
{
  static int firsttime=1;
  int spin,Mc_AN,Gc_AN,wan1,TNO1,h_AN,Gh_AN;
  int wan2,TNO2,i,j,k,NumMix,NumSlide;
  int SCFi,SCFj,tno1,tno0,Cwan,Hwan,n1,k1,k2,k3,kk2; 
  int pSCF_iter,size_Re_OptRhok,size_Kerker_weight;
  int spinmax,MN,flag_nan;
  int *One2SCFi,*One2SCFj;
  double time1,time2,time3,time4,time5; 
  double Stime1,Etime1;
  double *alden;
  double **A,**My_A,**IA,*My_A_thrds;
  double tmp0,itmp0,tmp1,im1,im2,re1,re2;
  double OptRNorm,My_OptRNorm,x,wt;
  double dum1,dum2,bunsi,bunbo,sum,coef_OptRho;
  double Gx,Gy,Gz,G2,G12,G22,G32;
  double Min_Weight,Q,sumL,sumR;
  double Max_Weight;
  double G0,G02,G02p,weight,wgt0,wgt1;
  double sk1,sk2,sk3;
  double ****Re_OptRhok;
  double ****Im_OptRhok;
  double ***Kerker_weight;
  int numprocs,myid;
  char nanchar[300];
  /* for OpenMP */
  int OMPID,Nthrds,Nprocs,Nloop,Nthrds0;
  int knum,knum_full;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /* find an optimum G0 */

  G22 = rtv[2][1]*rtv[2][1] + rtv[2][2]*rtv[2][2] + rtv[2][3]*rtv[2][3]; 
  G32 = rtv[3][1]*rtv[3][1] + rtv[3][2]*rtv[3][2] + rtv[3][3]*rtv[3][3]; 

  if (G22<G32) G0 = G22;
  else         G0 = G32;

  G0 = sqrt(G0);
  G02 = Kerker_factor*Kerker_factor*G0*G0;
  G02p = (0.3*Kerker_factor*G0)*(0.3*Kerker_factor*G0);

  if      (SpinP_switch==0)  spinmax = 1;
  else if (SpinP_switch==1)  spinmax = 2;
  else if (SpinP_switch==3)  spinmax = 3;

  /****************************************************
    allocation of arrays:

    double alden[List_YOUSO[38]];
    double A[List_YOUSO[38]][List_YOUSO[38]];
    double IA[List_YOUSO[38]][List_YOUSO[38]];
    double My_A[List_YOUSO[38]][List_YOUSO[38]];
    double Re_OptRhok[spinmax]
                            [My_NGrid2_Poisson]
                            [Ngrid1]
                            [Ngrid3];
    double Im_OptRhok[spinmax]
                            [My_NGrid2_Poisson]
                            [Ngrid1]
                            [Ngrid3];
  ****************************************************/
  
  alden = (double*)malloc(sizeof(double)*List_YOUSO[38]);

  A = (double**)malloc(sizeof(double*)*List_YOUSO[38]);
  for (i=0; i<List_YOUSO[38]; i++){
    A[i] = (double*)malloc(sizeof(double)*List_YOUSO[38]);
  }

  My_A = (double**)malloc(sizeof(double*)*List_YOUSO[38]);
  for (i=0; i<List_YOUSO[38]; i++){
    My_A[i] = (double*)malloc(sizeof(double)*List_YOUSO[38]);
  }

  IA = (double**)malloc(sizeof(double*)*List_YOUSO[38]);
  for (i=0; i<List_YOUSO[38]; i++){
    IA[i] = (double*)malloc(sizeof(double)*List_YOUSO[38]);
  }

  /* Re_OptRhok and Im_OptRhok */

  size_Re_OptRhok = spinmax*My_NGrid2_Poisson*Ngrid1*Ngrid3;

  Re_OptRhok = (double****)malloc(sizeof(double***)*spinmax); 
  for (spin=0; spin<spinmax; spin++){
    Re_OptRhok[spin] = (double***)malloc(sizeof(double**)*My_NGrid2_Poisson); 
    for (i=0; i<My_NGrid2_Poisson; i++){
      Re_OptRhok[spin][i] = (double**)malloc(sizeof(double*)*Ngrid1); 
      for (j=0; j<Ngrid1; j++){
	Re_OptRhok[spin][i][j] = (double*)malloc(sizeof(double)*Ngrid3); 
	for (k=0; k<Ngrid3; k++) Re_OptRhok[spin][i][j][k] = 0.0;
      }
    }
  }

  Im_OptRhok = (double****)malloc(sizeof(double***)*spinmax); 
  for (spin=0; spin<spinmax; spin++){
    Im_OptRhok[spin] = (double***)malloc(sizeof(double**)*My_NGrid2_Poisson); 
    for (i=0; i<My_NGrid2_Poisson; i++){
      Im_OptRhok[spin][i] = (double**)malloc(sizeof(double*)*Ngrid1); 
      for (j=0; j<Ngrid1; j++){
	Im_OptRhok[spin][i][j] = (double*)malloc(sizeof(double)*Ngrid3); 
	for (k=0; k<Ngrid3; k++) Im_OptRhok[spin][i][j][k] = 0.0;
      }
    }
  }

  /***********************************
            set Kerker_weight 
  ************************************/

  size_Kerker_weight = My_NGrid2_Poisson*Ngrid1*Ngrid3;
  knum = My_NGrid2_Poisson*Ngrid1*Ngrid3;
  knum_full = spinmax*My_NGrid2_Poisson*Ngrid1*Ngrid3;

  Kerker_weight = (double***)malloc(sizeof(double**)*My_NGrid2_Poisson); 
  for (i=0; i<My_NGrid2_Poisson; i++){
    Kerker_weight[i] = (double**)malloc(sizeof(double*)*Ngrid1); 
    for (j=0; j<Ngrid1; j++){
      Kerker_weight[i][j] = (double*)malloc(sizeof(double)*Ngrid3); 
    }
  }

  if (measure_time==1) dtime(&Stime1);

#pragma omp parallel shared(myid,Start_Grid2,knum,G02,G02p,Kerker_weight,rtv,Ngrid3,Ngrid2,Ngrid1) private(OMPID,Nthrds,Nprocs,k2,kk2,sk2,k,k1,k3,sk1,sk3,Gx,Gy,Gz,G2,weight)
  {

    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    /* one-dimensionalized loop */

    for (k=OMPID*knum/Nthrds; k<(OMPID+1)*knum/Nthrds; k++){

      /* get k2, k1, and k3 */
 
      k2 = k/(Ngrid1*Ngrid3);
      k1 = (k - k2*(Ngrid1*Ngrid3))/Ngrid3;
      k3 = k - k2*(Ngrid1*Ngrid3) - k1*Ngrid3;

      kk2 = k2 + Start_Grid2[myid];

      if (kk2<Ngrid2/2) sk2 = (double)kk2;
      else              sk2 = (double)(kk2 - Ngrid2);

      if (k3<Ngrid3/2) sk3 = (double)k3;
      else             sk3 = (double)(k3 - Ngrid3);

      Gx = sk2*rtv[2][1] + sk3*rtv[3][1];
      Gy = sk2*rtv[2][2] + sk3*rtv[3][2]; 
      Gz = sk2*rtv[2][3] + sk3*rtv[3][3];

      G2 = Gx*Gx + Gy*Gy + Gz*Gz;

      if (kk2==0 && k3==0){

	sumL = 0.0;
	sumR = 0.0; 

	if (SpinP_switch==0){

	  for (n1=0; n1<k1; n1++){
	    sumL += 2.0*ReRhok[0][0][k2][n1][k3];
	  }

	  for (n1=k1+1; n1<Ngrid1; n1++){
	    sumR += 2.0*ReRhok[0][0][k2][n1][k3];
	  }
	}
        
	else if (SpinP_switch==1){

	  for (n1=0; n1<k1; n1++){
	    sumL += ReRhok[0][0][k2][n1][k3] + ReRhok[0][1][k2][n1][k3];
	  }

	  for (n1=k1+1; n1<Ngrid1; n1++){
	    sumR += ReRhok[0][0][k2][n1][k3] + ReRhok[0][1][k2][n1][k3];
	  }
	}

	Q = 4.0*fabs(sumL - sumR)*GridVol + 1.0;

      }    
      else {
	Q = 1.0;
      }

      weight = (G2 + G02)/(G2 + G02p)*Q;

      Kerker_weight[k2][k1][k3] = weight;

    } /* k */
  } /* #pragma omp parallel */

  if (measure_time==1){ 
    dtime(&Etime1);
    time1 = Etime1 - Stime1;      
  }

  /****************************************************
     start calc.
  ****************************************************/
 
  if (SCF_iter==1){
    /* memory calc. done, only 1st iteration */
    if (firsttime) {
      PrintMemory("DIIS_Mixing_Rhok: Re_OptRhok",sizeof(double)*size_Re_OptRhok,NULL);
      PrintMemory("DIIS_Mixing_Rhok: Im_OptRhok",sizeof(double)*size_Re_OptRhok,NULL);
      PrintMemory("DIIS_Mixing_Rhok: Kerker_weight",sizeof(double)*size_Kerker_weight,NULL);
      firsttime=0;
    }

    Kerker_Mixing_Rhok(0,1.00,ReRhok,ImRhok,ReBestRhok,ImBestRhok,
                       Residual_ReRhok,Residual_ImRhok,
                       ReV1,ImV1,ReV2,ImV2,ReRhoAtomk,ImRhoAtomk);
  } 

  else if (SCF_iter==2){

    if (measure_time==1) dtime(&Stime1);

    /* construct residual rho */

#pragma omp parallel shared(Residual_ReRhok,Residual_ImRhok,ImRhok,ReRhok,knum_full,My_NGrid2_Poisson,Ngrid1,Ngrid3) private(OMPID,Nthrds,Nprocs,k,spin,k2,k1,k3,tmp0,tmp1)
    {

      /* get info. on OpenMP */ 

      OMPID = omp_get_thread_num();
      Nthrds = omp_get_num_threads();
      Nprocs = omp_get_num_procs();

      /* one-dimensionalized loop */

      for (k=OMPID*knum_full/Nthrds; k<(OMPID+1)*knum_full/Nthrds; k++){

	/* get spin, k2, k1, and k3 */

	spin = k/(My_NGrid2_Poisson*Ngrid1*Ngrid3);
	k2 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3))/(Ngrid1*Ngrid3);
	k1 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3))/Ngrid3;  
	k3 = k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3) - k1*Ngrid3;

	tmp0 = ReRhok[0][spin][k2][k1][k3] - ReRhok[1][spin][k2][k1][k3];
	tmp1 = ImRhok[0][spin][k2][k1][k3] - ImRhok[1][spin][k2][k1][k3];
	Residual_ReRhok[2][spin][k2][k1][k3] = tmp0;
	Residual_ImRhok[2][spin][k2][k1][k3] = tmp1;

      } /* k */
    } /* #pragma omp parallel */

    if (measure_time==1){ 
      dtime(&Etime1);
      time2 = Etime1 - Stime1;      
    }

    /* rho1 to rho2 */

    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
        for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){
	    ReRhok[2][spin][k2][k1][k3] = ReRhok[1][spin][k2][k1][k3];
	    ImRhok[2][spin][k2][k1][k3] = ImRhok[1][spin][k2][k1][k3];
	  }
	}
      }
    }

    Kerker_Mixing_Rhok(0,Mixing_weight,ReRhok,ImRhok,ReBestRhok,ImBestRhok,
                       Residual_ReRhok,Residual_ImRhok,
                       ReV1,ImV1,ReV2,ImV2,ReRhoAtomk,ImRhoAtomk);
  } 

  else {

    /****************************************************
                   construct residual rho1
    ****************************************************/

    if (measure_time==1) dtime(&Stime1);

#pragma omp parallel shared(Residual_ReRhok,Residual_ImRhok,ImRhok,ReRhok,knum_full,My_NGrid2_Poisson,Ngrid1,Ngrid3) private(OMPID,Nthrds,Nprocs,k,spin,k2,k1,k3,tmp0,tmp1)
    {
      /* get info. on OpenMP */ 

      OMPID = omp_get_thread_num();
      Nthrds = omp_get_num_threads();
      Nprocs = omp_get_num_procs();

      /* one-dimensionalized loop */

      for (k=OMPID*knum_full/Nthrds; k<(OMPID+1)*knum_full/Nthrds; k++){

	/* get spin, k2, k1, and k3 */

	spin = k/(My_NGrid2_Poisson*Ngrid1*Ngrid3);
	k2 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3))/(Ngrid1*Ngrid3);
	k1 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3))/Ngrid3;  
	k3 = k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3) - k1*Ngrid3;

	tmp0 = ReRhok[0][spin][k2][k1][k3] - ReRhok[1][spin][k2][k1][k3];
	tmp1 = ImRhok[0][spin][k2][k1][k3] - ImRhok[1][spin][k2][k1][k3];
	Residual_ReRhok[1][spin][k2][k1][k3] = tmp0;
	Residual_ImRhok[1][spin][k2][k1][k3] = tmp1;

      } /* k */
    } /* #pragma omp parallel */

    if (measure_time==1){ 
      dtime(&Etime1);
      time3 = Etime1 - Stime1;      
    }

    if ((SCF_iter-1)<Num_Mixing_pDM){
      NumMix   = SCF_iter - 1;
      NumSlide = NumMix + 1;
    }
    else{
      NumMix   = Num_Mixing_pDM;
      NumSlide = NumMix;
    }

    /****************************************************
                    alpha from residual rho
    ****************************************************/

    if (measure_time==1) dtime(&Stime1);

    /* allocation of arrrays */

    One2SCFi = (int*)malloc(sizeof(int)*(NumMix*(NumMix+1))/2); 
    One2SCFj = (int*)malloc(sizeof(int)*(NumMix*(NumMix+1))/2); 

    i = 0;
    for (SCFi=1; SCFi<=NumMix; SCFi++){
      for (SCFj=SCFi; SCFj<=NumMix; SCFj++){
        One2SCFi[i] = SCFi;
        One2SCFj[i] = SCFj;
        i++; 
      }
    }
    
    /* calculation of the norm matrix for the residual vectors */

#pragma omp parallel shared(NumMix,One2SCFi,One2SCFj,My_A,spinmax,My_NGrid2_Poisson,Ngrid1,Ngrid3,Kerker_weight,Residual_ReRhok,Residual_ImRhok) private(OMPID,Nthrds,Nprocs,i,SCFi,SCFj,spin,k2,k1,k3,weight,re1,im1,re2,im2,tmp0)
    {

      /* get info. on OpenMP */ 

      OMPID = omp_get_thread_num();
      Nthrds = omp_get_num_threads();
      Nprocs = omp_get_num_procs();

      for ( i=OMPID; i<(NumMix*(NumMix+1))/2; i+=Nthrds ){

	SCFi = One2SCFi[i]; 
	SCFj = One2SCFj[i]; 

        tmp0 = 0.0; 

	for (spin=0; spin<spinmax; spin++){
	  for (k2=0; k2<My_NGrid2_Poisson; k2++){
	    for (k1=0; k1<Ngrid1; k1++){
	      for (k3=0; k3<Ngrid3; k3++){

		weight = Kerker_weight[k2][k1][k3]; 

		re1 = Residual_ReRhok[SCFi][spin][k2][k1][k3];
		im1 = Residual_ImRhok[SCFi][spin][k2][k1][k3];

		re2 = Residual_ReRhok[SCFj][spin][k2][k1][k3];
		im2 = Residual_ImRhok[SCFj][spin][k2][k1][k3];

		tmp0 += (re1*re2 + im1*im2)*weight;
	      }
	    }
	  }
	}

	My_A[SCFi][SCFj] = tmp0;

      }
    } /* pragma omp parallel */

    /* freeing of arrrays */
    free(One2SCFj);
    free(One2SCFi);

    if (measure_time==1){ 
      dtime(&Etime1);
      time4 = Etime1 - Stime1;      
    }

    /* MPI My_A */
    for (SCFi=1; SCFi<=NumMix; SCFi++){
      for (SCFj=SCFi; SCFj<=NumMix; SCFj++){
        MPI_Allreduce(&My_A[SCFi][SCFj], &A[SCFi-1][SCFj-1],
                      1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
        A[SCFj-1][SCFi-1] = A[SCFi-1][SCFj-1];
      }
    }

    /* store NormRD */

    for (i=4; 1<=i; i--){
      NormRD[i] = NormRD[i-1];
      History_Uele[i] = History_Uele[i-1];
    }
    NormRD[0] = A[0][0];
    History_Uele[0] = Uele;

    /* solve the linear equation */

    for (SCFi=1; SCFi<=NumMix; SCFi++){
       A[SCFi-1][NumMix] = -1.0;
       A[NumMix][SCFi-1] = -1.0;
    }
    A[NumMix][NumMix] = 0.0;

    Inverse(NumMix,A,IA);

    for (SCFi=1; SCFi<=NumMix; SCFi++){
       alden[SCFi] = -IA[SCFi-1][NumMix];
    }

    /*************************************
      check "nan", "NaN", "inf" or "Inf"
    *************************************/

    flag_nan = 0;
    for (SCFi=1; SCFi<=NumMix; SCFi++){

      sprintf(nanchar,"%8.4f",alden[SCFi]);
      if (strstr(nanchar,"nan")!=NULL || strstr(nanchar,"NaN")!=NULL 
       || strstr(nanchar,"inf")!=NULL || strstr(nanchar,"Inf")!=NULL){

        flag_nan = 1;
      }
    }

    if (flag_nan==1){
      for (SCFi=1; SCFi<=NumMix; SCFi++){
        alden[SCFi] = 0.0;
      }
      alden[1] = 0.1;
      alden[2] = 0.9;
    }

    /****************************************************
              calculate optimized residual rho
    ****************************************************/

    if (measure_time==1) dtime(&Stime1);

    for (pSCF_iter=1; pSCF_iter<=NumMix; pSCF_iter++){

      tmp0 =  alden[pSCF_iter]; 

#pragma omp parallel shared(pSCF_iter,tmp0,Im_OptRhok,Re_OptRhok,Residual_ReRhok,Residual_ImRhok,knum_full,My_NGrid2_Poisson,Ngrid1,Ngrid3) private(OMPID,Nthrds,Nprocs,k,spin,k2,k1,k3)
      {
	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	/* one-dimensionalized loop */

	for (k=OMPID*knum_full/Nthrds; k<(OMPID+1)*knum_full/Nthrds; k++){

	  /* get spin, k2, k1, and k3 */

	  spin = k/(My_NGrid2_Poisson*Ngrid1*Ngrid3);
	  k2 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3))/(Ngrid1*Ngrid3);
	  k1 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3))/Ngrid3;  
	  k3 = k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3) - k1*Ngrid3;

	  /* Pulay mixing */

	  Re_OptRhok[spin][k2][k1][k3] += tmp0*Residual_ReRhok[pSCF_iter][spin][k2][k1][k3];
	  Im_OptRhok[spin][k2][k1][k3] += tmp0*Residual_ImRhok[pSCF_iter][spin][k2][k1][k3];

	} /* k */
      } /* #pragma omp parallel */
    } /* pSCF_iter */

    /* the norm of the optimized residual rho */

    /*
    My_OptRNorm = 0.0;
    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
	for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){

	    weight = Kerker_weight[k2][k1][k3]; 

	    re1 = Re_OptRhok[spin][k2][k1][k3];
	    im1 = Im_OptRhok[spin][k2][k1][k3];

	    My_OptRNorm += (re1*re1 + im1*im1)*weight;
	  }
	}
      }
    }

    MPI_Allreduce(&My_OptRNorm, &OptRNorm, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
    OptRNorm = OptRNorm/(double)Ngrid1/(double)Ngrid2/(double)Ngrid3/(double)spinmax;
    printf("OptRNorm=%15.12f\n",sqrt(OptRNorm)); 
    */

    if      (1.0e-4<=NormRD[0])                        coef_OptRho = 0.5;
    else if (1.0e-6<=NormRD[0]  && NormRD[0]<1.0e-4)   coef_OptRho = 0.6;
    else if (1.0e-8<=NormRD[0]  && NormRD[0]<1.0e-6)   coef_OptRho = 0.7;
    else if (1.0e-10<=NormRD[0] && NormRD[0]<1.0e-8)   coef_OptRho = 0.8;
    else if (1.0e-12<=NormRD[0] && NormRD[0]<1.0e-10)  coef_OptRho = 0.9;
    else                                               coef_OptRho = 1.0;

    /****************************************************
            store the best input charge density 
    ****************************************************/

    if (NormRD[0]<BestNormRD ){

      BestNormRD = NormRD[0];

      for (spin=0; spin<spinmax; spin++){
        for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){
              ReBestRhok[spin][k2][k1][k3] = ReRhok[1][spin][k2][k1][k3];
              ImBestRhok[spin][k2][k1][k3] = ImRhok[1][spin][k2][k1][k3];
	    }
	  }
        }
      }
    }

    /****************************************************
                        mixing of rho 
    ****************************************************/

    /****************************************************
       Pulay mixing
    ****************************************************/

    if ( SCF_iter%EveryPulay_SCF==0 ) {

      /* reduce Mixing_weight so that charge sloshing can be avoided at the next SCF */

      Mixing_weight = 0.1*Max_Mixing_weight;
      if (Mixing_weight<Min_Mixing_weight) Mixing_weight = Min_Mixing_weight;

      /* initial ReRhok and ImRhok */

      for (spin=0; spin<spinmax; spin++){
	for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){
	      ReRhok[0][spin][k2][k1][k3] = 0.0;
	      ImRhok[0][spin][k2][k1][k3] = 0.0;
	    }
	  }
	}
      }

      /* Pulay mixing */

      for (pSCF_iter=1; pSCF_iter<=NumMix; pSCF_iter++){

	tmp0 =  alden[pSCF_iter];

#pragma omp parallel shared(pSCF_iter,tmp0,ReRhok,ImRhok,knum_full,My_NGrid2_Poisson,Ngrid1,Ngrid3) private(OMPID,Nthrds,Nprocs,k,spin,k2,k1,k3)
	{
	  /* get info. on OpenMP */ 

	  OMPID = omp_get_thread_num();
	  Nthrds = omp_get_num_threads();
	  Nprocs = omp_get_num_procs();

	  /* one-dimensionalized loop */

	  for (k=OMPID*knum_full/Nthrds; k<(OMPID+1)*knum_full/Nthrds; k++){

	    /* get spin, k2, k1, and k3 */

	    spin = k/(My_NGrid2_Poisson*Ngrid1*Ngrid3);
	    k2 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3))/(Ngrid1*Ngrid3);
	    k1 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3))/Ngrid3;  
	    k3 = k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3) - k1*Ngrid3;

	    ReRhok[0][spin][k2][k1][k3] += tmp0*ReRhok[pSCF_iter][spin][k2][k1][k3];
	    ImRhok[0][spin][k2][k1][k3] += tmp0*ImRhok[pSCF_iter][spin][k2][k1][k3];

	  } /* k */
	} /* #pragma omp parallel */
      } /* pSCF_iter */

      /* Correction by optimized residual rho */

#pragma omp parallel shared(coef_OptRho,Im_OptRhok,Re_OptRhok,ReRhok,ImRhok,knum_full,My_NGrid2_Poisson,Ngrid1,Ngrid3) private(OMPID,Nthrds,Nprocs,k,spin,k2,k1,k3)
      {
	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	/* one-dimensionalized loop */

	for (k=OMPID*knum_full/Nthrds; k<(OMPID+1)*knum_full/Nthrds; k++){

	  /* get spin, k2, k1, and k3 */

	  spin = k/(My_NGrid2_Poisson*Ngrid1*Ngrid3);
	  k2 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3))/(Ngrid1*Ngrid3);
	  k1 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3))/Ngrid3;  
	  k3 = k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3) - k1*Ngrid3;

	  ReRhok[0][spin][k2][k1][k3] += coef_OptRho*Re_OptRhok[spin][k2][k1][k3];
	  ImRhok[0][spin][k2][k1][k3] += coef_OptRho*Im_OptRhok[spin][k2][k1][k3];

	  /* correction to large changing components */

          tmp0 = ReRhok[0][spin][k2][k1][k3] - ReRhok[1][spin][k2][k1][k3];  
          tmp1 = ImRhok[0][spin][k2][k1][k3] - ImRhok[1][spin][k2][k1][k3];  

          if ( maxima_step<(fabs(tmp0)+fabs(tmp1)) ){
            ReRhok[0][spin][k2][k1][k3] = sgn(tmp0)*maxima_step + ReRhok[1][spin][k2][k1][k3]; 
            ImRhok[0][spin][k2][k1][k3] = sgn(tmp1)*maxima_step + ImRhok[1][spin][k2][k1][k3]; 
          }

	} /* k */
      } /* #pragma omp parallel */

      if (measure_time==1){ 
        dtime(&Etime1);
        time5 = Etime1 - Stime1;      
      }

    }

    /****************************************************
       Kerker mixing
    ****************************************************/

    else {

      /* find an optimum mixing weight */

      Min_Weight = Min_Mixing_weight;
      Max_Weight = Max_Mixing_weight;

      if ((int)sgn(History_Uele[0]-History_Uele[1])
	  ==(int)sgn(History_Uele[1]-History_Uele[2])
	  && NormRD[0]<NormRD[1]){

	/* tmp0 = 1.6*Mixing_weight; */

	tmp0 = NormRD[1]/(largest(NormRD[1]-NormRD[0],10e-10))*Mixing_weight;

	if (tmp0<Max_Weight){
	  if (Min_Weight<tmp0){
	    Mixing_weight = tmp0;
	  }
	  else{ 
	    Mixing_weight = Min_Weight;
	  }
	}
	else{ 
	  Mixing_weight = Max_Weight;
	  SCF_RENZOKU++;  
	}
      }
   
      else if ((int)sgn(History_Uele[0]-History_Uele[1])
	       ==(int)sgn(History_Uele[1]-History_Uele[2])
	       && NormRD[1]<NormRD[0]){

	tmp0 = NormRD[1]/(largest(NormRD[1]+NormRD[0],10e-10))*Mixing_weight;

	/* tmp0 = Mixing_weight/1.6; */

	if (tmp0<Max_Weight){
	  if (Min_Weight<tmp0)
	    Mixing_weight = tmp0;
	  else 
	    Mixing_weight = Min_Weight;
	}
	else 
	  Mixing_weight = Max_Weight;

	SCF_RENZOKU = -1;  
      }

      else if ((int)sgn(History_Uele[0]-History_Uele[1])
	       !=(int)sgn(History_Uele[1]-History_Uele[2])
	       && NormRD[0]<NormRD[1]){

	/* tmp0 = Mixing_weight*1.2; */

	tmp0 = NormRD[1]/(largest(NormRD[1]-NormRD[0],10e-10))*Mixing_weight;

	if (tmp0<Max_Weight){
	  if (Min_Weight<tmp0)
	    Mixing_weight = tmp0;
	  else 
	    Mixing_weight = Min_Weight;
	}
	else{ 
	  Mixing_weight = Max_Weight;
	  SCF_RENZOKU++;
	}
      }

      else if ((int)sgn(History_Uele[0]-History_Uele[1])
	       !=(int)sgn(History_Uele[1]-History_Uele[2])
	       && NormRD[1]<NormRD[0]){

	/* tmp0 = Mixing_weight/2.0; */

	tmp0 = NormRD[1]/(largest(NormRD[1]+NormRD[0],10e-10))*Mixing_weight;

	if (tmp0<Max_Weight){
	  if (Min_Weight<tmp0)
	    Mixing_weight = tmp0;
	  else 
	    Mixing_weight = Min_Weight;
	}
	else 
	  Mixing_weight = Max_Weight;

	SCF_RENZOKU = -1;
      }

      Mix_wgt = Mixing_weight;

      /* Kerer mixing */

#pragma omp parallel shared(Mix_wgt,Kerker_weight,ReRhok,ImRhok,knum_full,My_NGrid2_Poisson,Ngrid1,Ngrid3) private(OMPID,Nthrds,Nprocs,k,spin,k2,k1,k3,weight,wgt0,wgt1)
      {
	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	/* one-dimensionalized loop */

	for (k=OMPID*knum_full/Nthrds; k<(OMPID+1)*knum_full/Nthrds; k++){

	  /* get spin, k2, k1, and k3 */

	  spin = k/(My_NGrid2_Poisson*Ngrid1*Ngrid3);
	  k2 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3))/(Ngrid1*Ngrid3);
	  k1 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3))/Ngrid3;  
	  k3 = k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3) - k1*Ngrid3;

	  weight = 1.0/Kerker_weight[k2][k1][k3];
	  wgt0  = Mix_wgt*weight;
	  wgt1 =  1.0 - wgt0;

	  ReRhok[0][spin][k2][k1][k3] = wgt0*ReRhok[0][spin][k2][k1][k3]
	                              + wgt1*ReRhok[1][spin][k2][k1][k3];

	  ImRhok[0][spin][k2][k1][k3] = wgt0*ImRhok[0][spin][k2][k1][k3]
	                              + wgt1*ImRhok[1][spin][k2][k1][k3];

	} /* k */
      } /* #pragma omp parallel */      
    } /* else */

    /****************************************************
                        shift of rho
    ****************************************************/

    for (pSCF_iter=(List_YOUSO[38]-1); 0<pSCF_iter; pSCF_iter--){
      for (spin=0; spin<spinmax; spin++){
	for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){
	      ReRhok[pSCF_iter][spin][k2][k1][k3] = ReRhok[pSCF_iter-1][spin][k2][k1][k3]; 
	      ImRhok[pSCF_iter][spin][k2][k1][k3] = ImRhok[pSCF_iter-1][spin][k2][k1][k3]; 
	    }
	  }
	}
      }
    }

    /****************************************************
                    shift of residual rho
    ****************************************************/

    for (pSCF_iter=(List_YOUSO[38]-1); 0<pSCF_iter; pSCF_iter--){
      for (spin=0; spin<spinmax; spin++){
	for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){
	      Residual_ReRhok[pSCF_iter][spin][k2][k1][k3] = Residual_ReRhok[pSCF_iter-1][spin][k2][k1][k3];
	      Residual_ImRhok[pSCF_iter][spin][k2][k1][k3] = Residual_ImRhok[pSCF_iter-1][spin][k2][k1][k3]; 
	    }
	  }
	}
      }
    }

  }

  /****************************************************
        find the charge density in real space 
  ****************************************************/

  for (spin=0; spin<spinmax; spin++){

    for (k2=0; k2<My_NGrid2_Poisson; k2++){
      for (k1=0; k1<Ngrid1; k1++){
        for (k3=0; k3<Ngrid3; k3++){
          ReV2[k2][k1][k3] = ReRhok[0][spin][k2][k1][k3];
          ImV2[k2][k1][k3] = ImRhok[0][spin][k2][k1][k3];
        }
      }
    }

    Get_Value_inReal2D(ReV2,ImV2,ReV1,ImV1,Density_Grid[spin]);
  }

  if (SpinP_switch==0){
    for (MN=0; MN<My_NumGrid1; MN++){
      Density_Grid[1][MN] = Density_Grid[0][MN];
    }
  }

  if (measure_time==1){ 
    printf("time1=%12.6f time2=%12.6f time3=%12.6f time4=%12.6f time5=%12.6f\n",time1,time2,time3,time4,time5);
  }

  /****************************************************
    freeing of arrays:

    double alden[List_YOUSO[38]];
    double A[List_YOUSO[38]][List_YOUSO[38]];
    double IA[List_YOUSO[38]][List_YOUSO[38]];
    double My_A[List_YOUSO[38]][List_YOUSO[38]];
    double Re_OptRhok[spinmax]
                            [My_NGrid2_Poisson]
                            [Ngrid1]
                            [Ngrid3];
    double Im_OptRhok[spinmax]
                            [My_NGrid2_Poisson]
                            [Ngrid1]
                            [Ngrid3];
  ****************************************************/
  
  free(alden);

  for (i=0; i<List_YOUSO[38]; i++){
    free(A[i]);
  }
  free(A);

  for (i=0; i<List_YOUSO[38]; i++){
    free(My_A[i]);
  }
  free(My_A);

  for (i=0; i<List_YOUSO[38]; i++){
    free(IA[i]);
  }
  free(IA);

  /* Re_OptRhok and Im_OptRhok */

  for (spin=0; spin<spinmax; spin++){
    for (i=0; i<My_NGrid2_Poisson; i++){
      for (j=0; j<Ngrid1; j++){
	free(Re_OptRhok[spin][i][j]);
      }
      free(Re_OptRhok[spin][i]);
    }
    free(Re_OptRhok[spin]);
  }
  free(Re_OptRhok);

  for (spin=0; spin<spinmax; spin++){
    for (i=0; i<My_NGrid2_Poisson; i++){
      for (j=0; j<Ngrid1; j++){
	free(Im_OptRhok[spin][i][j]);
      }
      free(Im_OptRhok[spin][i]);
    }
    free(Im_OptRhok[spin]);
  }
  free(Im_OptRhok);

  /* Kerker_weight */

  for (i=0; i<My_NGrid2_Poisson; i++){
    for (j=0; j<Ngrid1; j++){
      free(Kerker_weight[i][j]);
    }
    free(Kerker_weight[i]);
  }
  free(Kerker_weight);

}




void Broyden2(int SCF_iter,
	      double Mix_wgt,
	      double *****ReRhok,
	      double *****ImRhok,
	      double ****ReBestRhok,
	      double ****ImBestRhok,
	      double *****Residual_ReRhok,
	      double *****Residual_ImRhok,
	      double ***ReV1,
	      double ***ImV1,
	      double ***ReV2,
	      double ***ImV2,
	      double ***ReRhoAtomk,
	      double ***ImRhoAtomk)
{
  static int firsttime=1;
  int spin,Mc_AN,Gc_AN,wan1,TNO1,h_AN,Gh_AN;
  int wan2,TNO2,i,j,k,NumMix;
  int SCFi,SCFj,tno1,tno0,Cwan,Hwan,k1,k2,k3,kk2; 
  int pSCF_iter,size_Re_OptRhok,size_Kerker_weight;
  int spinmax,MN,flag_nan;
  double *alden;
  double My_Norm,My_OptNorm,OptNorm;
  double **A,**My_A,**IA,*B,*My_B,*Gamma;
  double *B1,*My_B1,*Alpha;
  double tmp0,itmp0,tmp1,tmp2,im1,im2,im3,re1,re2,re3;
  double OptRNorm,My_OptRNorm;
  double dum1,dum2,bunsi,bunbo,sum,coef_OptRho;
  double Gx,Gy,Gz,G2,G12,G22,G32;
  double Min_Weight;
  double Max_Weight;
  double G0,G02,G02p,weight,wgt0,wgt1;
  double sk1,sk2,sk3;
  double ****Re_OptRhok;
  double ****Im_OptRhok;
  double ***Kerker_weight;
  int numprocs,myid;
  char nanchar[300];

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /* find an optimum G0 */

  G12 = rtv[1][1]*rtv[1][1] + rtv[1][2]*rtv[1][2] + rtv[1][3]*rtv[1][3]; 
  G22 = rtv[2][1]*rtv[2][1] + rtv[2][2]*rtv[2][2] + rtv[2][3]*rtv[2][3]; 
  G32 = rtv[3][1]*rtv[3][1] + rtv[3][2]*rtv[3][2] + rtv[3][3]*rtv[3][3]; 

  if (G12<G22) G0 = G12;
  else         G0 = G22;
  if (G32<G0)  G0 = G32;

  G0 = Kerker_factor*sqrt(G0);
  G02 = G0*G0;
  G02p = (0.1*G0)*(0.1*G0);

  if      (SpinP_switch==0)  spinmax = 1;
  else if (SpinP_switch==1)  spinmax = 2;
  else if (SpinP_switch==3)  spinmax = 3;

  /****************************************************
    allocation of arrays:

    double alden[List_YOUSO[38]];
    double A[List_YOUSO[38]][List_YOUSO[38]];
    double IA[List_YOUSO[38]][List_YOUSO[38]];
    double My_A[List_YOUSO[38]][List_YOUSO[38]];
    double Re_OptRhok[spinmax]
                            [My_NGrid2_Poisson]
                            [Ngrid1]
                            [Ngrid3];
    double Im_OptRhok[spinmax]
                            [My_NGrid2_Poisson]
                            [Ngrid1]
                            [Ngrid3];
  ****************************************************/
  
  alden = (double*)malloc(sizeof(double)*List_YOUSO[38]);

  A = (double**)malloc(sizeof(double*)*List_YOUSO[38]);
  for (i=0; i<List_YOUSO[38]; i++){
    A[i] = (double*)malloc(sizeof(double)*List_YOUSO[38]);
  }

  My_A = (double**)malloc(sizeof(double*)*List_YOUSO[38]);
  for (i=0; i<List_YOUSO[38]; i++){
    My_A[i] = (double*)malloc(sizeof(double)*List_YOUSO[38]);
  }

  IA = (double**)malloc(sizeof(double*)*List_YOUSO[38]);
  for (i=0; i<List_YOUSO[38]; i++){
    IA[i] = (double*)malloc(sizeof(double)*List_YOUSO[38]);
    for (j=0; j<List_YOUSO[38]; j++){
      IA[i][j] = 0.0;
    }
  }

  My_B = (double*)malloc(sizeof(double)*List_YOUSO[38]);
  B = (double*)malloc(sizeof(double)*List_YOUSO[38]);
  Gamma = (double*)malloc(sizeof(double)*List_YOUSO[38]);

  B1 = (double*)malloc(sizeof(double)*List_YOUSO[38]);
  My_B1 = (double*)malloc(sizeof(double)*List_YOUSO[38]);
  Alpha = (double*)malloc(sizeof(double)*List_YOUSO[38]);

  /* Re_OptRhok and Im_OptRhok */

  size_Re_OptRhok = spinmax*My_NGrid2_Poisson*Ngrid1*Ngrid3;

  Re_OptRhok = (double****)malloc(sizeof(double***)*spinmax); 
  for (spin=0; spin<spinmax; spin++){
    Re_OptRhok[spin] = (double***)malloc(sizeof(double**)*My_NGrid2_Poisson); 
    for (i=0; i<My_NGrid2_Poisson; i++){
      Re_OptRhok[spin][i] = (double**)malloc(sizeof(double*)*Ngrid1); 
      for (j=0; j<Ngrid1; j++){
	Re_OptRhok[spin][i][j] = (double*)malloc(sizeof(double)*Ngrid3); 
	for (k=0; k<Ngrid3; k++) Re_OptRhok[spin][i][j][k] = 0.0;
      }
    }
  }

  Im_OptRhok = (double****)malloc(sizeof(double***)*spinmax); 
  for (spin=0; spin<spinmax; spin++){
    Im_OptRhok[spin] = (double***)malloc(sizeof(double**)*My_NGrid2_Poisson); 
    for (i=0; i<My_NGrid2_Poisson; i++){
      Im_OptRhok[spin][i] = (double**)malloc(sizeof(double*)*Ngrid1); 
      for (j=0; j<Ngrid1; j++){
	Im_OptRhok[spin][i][j] = (double*)malloc(sizeof(double)*Ngrid3); 
	for (k=0; k<Ngrid3; k++) Im_OptRhok[spin][i][j][k] = 0.0;
      }
    }
  }

  /* Kerker_weight */

  size_Kerker_weight = My_NGrid2_Poisson*Ngrid1*Ngrid3;
  
  Kerker_weight = (double***)malloc(sizeof(double**)*My_NGrid2_Poisson); 
  for (i=0; i<My_NGrid2_Poisson; i++){
    Kerker_weight[i] = (double**)malloc(sizeof(double*)*Ngrid1); 
    for (j=0; j<Ngrid1; j++){
      Kerker_weight[i][j] = (double*)malloc(sizeof(double)*Ngrid3); 
    }
  }

  for (k2=0; k2<My_NGrid2_Poisson; k2++){

    kk2 = k2 + Start_Grid2[myid];

    if (kk2<Ngrid2/2) sk2 = (double)kk2;
    else              sk2 = (double)(kk2 - Ngrid2);

    for (k1=0; k1<Ngrid1; k1++){

      if (k1<Ngrid1/2) sk1 = (double)k1;
      else             sk1 = (double)(k1 - Ngrid1);

      for (k3=0; k3<Ngrid3; k3++){

	if (k3<Ngrid3/2) sk3 = (double)k3;
	else             sk3 = (double)(k3 - Ngrid3);

	Gx = sk1*rtv[1][1] + sk2*rtv[2][1] + sk3*rtv[3][1];
	Gy = sk1*rtv[1][2] + sk2*rtv[2][2] + sk3*rtv[3][2]; 
	Gz = sk1*rtv[1][3] + sk2*rtv[2][3] + sk3*rtv[3][3];
	G2 = Gx*Gx + Gy*Gy + Gz*Gz;

	if (k1==0 && kk2==0 && k3==0)  weight = 1.0;
	else                           weight = (G2 + G02)/(G2 + G02p);

        Kerker_weight[k2][k1][k3] = weight;
      }
    }
  }

  /* memory calc. done, only 1st iteration */
  if (firsttime) {
    PrintMemory("DIIS_Mixing_Rhok: Re_OptRhok",sizeof(double)*size_Re_OptRhok,NULL);
    PrintMemory("DIIS_Mixing_Rhok: Im_OptRhok",sizeof(double)*size_Re_OptRhok,NULL);
    PrintMemory("DIIS_Mixing_Rhok: Kerker_weight",sizeof(double)*size_Kerker_weight,NULL);
    firsttime=0;
  }

  /****************************************************
     start calc.
  ****************************************************/

  /*
  printf("SCF_iter=%2d\n",SCF_iter);fflush(stdout);
  */

  SCF_iter--;
  NumMix = SCF_iter%(Num_Mixing_pDM+1);
  if (NumMix<2) NumMix = 2;

  /*
  if (Num_Mixing_pDM<NumMix){
    NumMix = Num_Mixing_pDM%(Num_Mixing_pDM+1);
  }
  */

  if (SCF_iter<=1){

    /* construct residual rho */

    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
        for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){
	    tmp0 = ReRhok[0][spin][k2][k1][k3] - ReRhok[1][spin][k2][k1][k3];
	    tmp1 = ImRhok[0][spin][k2][k1][k3] - ImRhok[1][spin][k2][k1][k3];
	    Residual_ReRhok[2][spin][k2][k1][k3] = tmp0;
	    Residual_ImRhok[2][spin][k2][k1][k3] = tmp1;
	  }
	}
      }
    }

    /* rho1 to rho2 */

    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
        for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){
	    ReRhok[2][spin][k2][k1][k3] = ReRhok[1][spin][k2][k1][k3];
	    ImRhok[2][spin][k2][k1][k3] = ImRhok[1][spin][k2][k1][k3];
	  }
	}
      }
    }

    Kerker_Mixing_Rhok(1,Mixing_weight,ReRhok,ImRhok,ReBestRhok,ImBestRhok,
                       Residual_ReRhok,Residual_ImRhok,
                       ReV1,ImV1,ReV2,ImV2,ReRhoAtomk,ImRhoAtomk);
  } 

  else {

    /* initialize DiaG */

    /*
    if (SCF_iter==2){
      for (spin=0; spin<spinmax; spin++){
	for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){
	      Re_DiaG[spin][k2][k1][k3] = -1.0/Kerker_weight[k2][k1][k3];
	      Im_DiaG[spin][k2][k1][k3] = -1.0/Kerker_weight[k2][k1][k3];
	    }
	  }
	}
      }
    } 
    */
 
    /****************************************************
       construct residual rho1
       residual = output - input  
    ****************************************************/

    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
	for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){
	    tmp0 = ReRhok[0][spin][k2][k1][k3] - ReRhok[1][spin][k2][k1][k3];
	    tmp1 = ImRhok[0][spin][k2][k1][k3] - ImRhok[1][spin][k2][k1][k3];
	    Residual_ReRhok[1][spin][k2][k1][k3] = tmp0;
	    Residual_ImRhok[1][spin][k2][k1][k3] = tmp1;
	  }
	}
      }
    }

    /* calculate NormRD */
  
    My_Norm = 0.0;
    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
	for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){

	    weight = Kerker_weight[k2][k1][k3]; 

	    re1 = Residual_ReRhok[1][spin][k2][k1][k3];
	    im1 = Residual_ImRhok[1][spin][k2][k1][k3];

	    My_Norm += (re1*re1 + im1*im1)*weight;
	  }
	}
      }
    }

    MPI_Allreduce(&My_Norm, &NormRD[0], 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

    for (i=4; 1<=i; i--){
      NormRD[i] = NormRD[i-1];
      History_Uele[i] = History_Uele[i-1];
    }
    NormRD[0] = NormRD[0]/(double)Ngrid1/(double)Ngrid2/(double)Ngrid3/(double)spinmax;
    History_Uele[0] = Uele;

    /****************************************************
     The matrix A calculated by difference residual rho
    ****************************************************/
      
    for (SCFi=1; SCFi<NumMix; SCFi++){
      for (SCFj=SCFi; SCFj<NumMix; SCFj++){

	My_A[SCFi][SCFj] = 0.0;

	for (spin=0; spin<spinmax; spin++){
	  for (k2=0; k2<My_NGrid2_Poisson; k2++){
	    for (k1=0; k1<Ngrid1; k1++){
	      for (k3=0; k3<Ngrid3; k3++){

		weight = Kerker_weight[k2][k1][k3]; 

		re1 = Residual_ReRhok[SCFi][spin][k2][k1][k3] - Residual_ReRhok[SCFi+1][spin][k2][k1][k3];
		im1 = Residual_ImRhok[SCFi][spin][k2][k1][k3] - Residual_ImRhok[SCFi+1][spin][k2][k1][k3];

		re2 = Residual_ReRhok[SCFj][spin][k2][k1][k3] - Residual_ReRhok[SCFj+1][spin][k2][k1][k3];
		im2 = Residual_ImRhok[SCFj][spin][k2][k1][k3] - Residual_ImRhok[SCFj+1][spin][k2][k1][k3];

		My_A[SCFi][SCFj] += (re1*re2 + im1*im2)*weight;
	      }
	    }
	  }
	}
      }
    }

    /* MPI My_A */

    for (SCFi=1; SCFi<NumMix; SCFi++){
      for (SCFj=SCFi; SCFj<NumMix; SCFj++){
	MPI_Allreduce(&My_A[SCFi][SCFj], &A[SCFi-1][SCFj-1],
		      1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
	A[SCFj-1][SCFi-1] = A[SCFi-1][SCFj-1];
      }
    }

    /*
    printf("A myid=%2d SCF_iter=%2d\n",myid,SCF_iter);fflush(stdout);
    for (SCFi=1; SCFi<NumMix; SCFi++){
      for (SCFj=1; SCFj<NumMix; SCFj++){
	printf("%18.13f ",A[SCFi-1][SCFj-1]);fflush(stdout);  
      }
      printf("\n");fflush(stdout);
    }
    */

    /* calculate <delta F(k)|F(l)> */

    for (SCFi=1; SCFi<NumMix; SCFi++){

      My_B[SCFi-1] = 0.0;

      for (spin=0; spin<spinmax; spin++){
	for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){

	      weight = Kerker_weight[k2][k1][k3]; 

	      re1 = Residual_ReRhok[SCFi][spin][k2][k1][k3] - Residual_ReRhok[SCFi+1][spin][k2][k1][k3];
	      im1 = Residual_ImRhok[SCFi][spin][k2][k1][k3] - Residual_ImRhok[SCFi+1][spin][k2][k1][k3];

	      re2 = Residual_ReRhok[1][spin][k2][k1][k3];
	      im2 = Residual_ImRhok[1][spin][k2][k1][k3];

	      My_B[SCFi-1] += (re1*re2 + im1*im2)*weight; 
	    }
	  }
	}
      }      
    }

    /* MPI My_B */

    for (SCFi=1; SCFi<NumMix; SCFi++){
      MPI_Allreduce(&My_B[SCFi-1], &B[SCFi-1], 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
    }

    /*
    printf("B myid=%2d SCF_iter=%2d\n",myid,SCF_iter);fflush(stdout);
    for (SCFi=1; SCFi<NumMix; SCFi++){
      printf("%18.13f\n",B[SCFi-1]); fflush(stdout);   
    }
    */

    /* calculate the inverse of A */

    if (2<=NumMix) Inverse(NumMix-2,A,IA);

    /*
    printf("IA myid=%2d SCF_iter=%2d\n",myid,SCF_iter);fflush(stdout);
    for (SCFi=1; SCFi<NumMix; SCFi++){
      for (SCFj=1; SCFj<NumMix; SCFj++){
	printf("%18.13f ",IA[SCFi-1][SCFj-1]); fflush(stdout); 
      }
      printf("\n");fflush(stdout);
    }
    */

    /* calculate gamma */
    
    for (SCFi=1; SCFi<NumMix; SCFi++){
      Gamma[SCFi] = 0.0;
      for (SCFj=1; SCFj<NumMix; SCFj++){
	Gamma[SCFi] += IA[SCFi-1][SCFj-1]*B[SCFj-1];
      }
    }

    printf("Gamma myid=%2d SCF_iter=%2d\n",myid,SCF_iter);fflush(stdout);
    for (SCFi=1; SCFi<NumMix; SCFi++){
      printf("%18.13f\n",Gamma[SCFi]);  fflush(stdout);  
    }

    /* calculate the optimum residual vector */

    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
	for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){
	    ReBestRhok[spin][k2][k1][k3] = Residual_ReRhok[1][spin][k2][k1][k3];
	    ImBestRhok[spin][k2][k1][k3] = Residual_ImRhok[1][spin][k2][k1][k3];
	  }
	}
      }
    }

    for (SCFi=1; SCFi<NumMix; SCFi++){
      for (spin=0; spin<spinmax; spin++){
	for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){

	      re2 = Residual_ReRhok[SCFi][spin][k2][k1][k3] - Residual_ReRhok[SCFi+1][spin][k2][k1][k3];
	      im2 = Residual_ImRhok[SCFi][spin][k2][k1][k3] - Residual_ImRhok[SCFi+1][spin][k2][k1][k3];

	      ReBestRhok[spin][k2][k1][k3] -= Gamma[SCFi]*re2;
	      ImBestRhok[spin][k2][k1][k3] -= Gamma[SCFi]*im2;
	    }
	  }
	}
      }
    }

    /* calculate the optimum norm */

    My_OptNorm = 0.0;
    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
	for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){

	    weight = Kerker_weight[k2][k1][k3]; 

	    re1 = ReBestRhok[spin][k2][k1][k3];
	    im1 = ImBestRhok[spin][k2][k1][k3];

	    My_OptNorm += (re1*re1 + im1*im1)*weight;
	  }
	}
      }
    }

    MPI_Allreduce(&My_OptNorm, &OptNorm, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

    /*
    printf("myid=%2d SCF_iter=%2d OptNorm=%15.12f\n",myid,SCF_iter,OptNorm);fflush(stdout);
    */

    /* calculate <delta F(k)| \overbar{F}(l)> */

    for (SCFi=1; SCFi<NumMix; SCFi++){

      My_B1[SCFi-1] = 0.0;

      for (spin=0; spin<spinmax; spin++){
	for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){

	      weight = Kerker_weight[k2][k1][k3]; 

	      re1 = Residual_ReRhok[SCFi][spin][k2][k1][k3] - Residual_ReRhok[SCFi+1][spin][k2][k1][k3];
	      im1 = Residual_ImRhok[SCFi][spin][k2][k1][k3] - Residual_ImRhok[SCFi+1][spin][k2][k1][k3];

	      re2 = ReBestRhok[spin][k2][k1][k3];
	      im2 = ImBestRhok[spin][k2][k1][k3];

	      My_B1[SCFi-1] += (re1*re2 + im1*im2)*weight; 
	    }
	  }
	}
      }      
    }

    /* MPI My_B1 */

    for (SCFi=1; SCFi<NumMix; SCFi++){
      MPI_Allreduce(&My_B1[SCFi-1], &B1[SCFi-1], 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
    }

    /*
    if (myid==0){
    printf("B1 myid=%2d SCF_iter=%2d\n",myid,SCF_iter);fflush(stdout);
    for (SCFi=1; SCFi<NumMix; SCFi++){
      printf("%15.12f\n",B1[SCFi-1]); fflush(stdout);   
    }
    }
    */
    
    /* calculate alpha */
    
    for (SCFi=1; SCFi<NumMix; SCFi++){
      Alpha[SCFi] = 0.0;
      for (SCFj=1; SCFj<NumMix; SCFj++){
	Alpha[SCFi] += IA[SCFi-1][SCFj-1]*B1[SCFj-1];
      }
    }
    
    /*
    printf("Alpha myid=%2d SCF_iter=%2d\n",myid,SCF_iter);fflush(stdout);
    for (SCFi=1; SCFi<NumMix; SCFi++){
      printf("%15.12f\n",Alpha[SCFi]);  fflush(stdout);
    }
    */

    /****************************************************
                       update DiaG
    ****************************************************/

    /*
    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
	for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){

            tmp1 = Re_DiaG[spin][k2][k1][k3];
            tmp2 = Im_DiaG[spin][k2][k1][k3];

	    for (SCFi=1; SCFi<NumMix; SCFi++){
	      for (SCFj=1; SCFj<NumMix; SCFj++){

	        re1 = ReRhok[SCFi][spin][k2][k1][k3] - ReRhok[SCFi+1][spin][k2][k1][k3];
	        im1 = ImRhok[SCFi][spin][k2][k1][k3] - ImRhok[SCFi+1][spin][k2][k1][k3];

	        re2 = Residual_ReRhok[SCFj][spin][k2][k1][k3] - Residual_ReRhok[SCFj+1][spin][k2][k1][k3];
	        im2 = Residual_ImRhok[SCFj][spin][k2][k1][k3] - Residual_ImRhok[SCFj+1][spin][k2][k1][k3];

	        re3 = Residual_ReRhok[SCFi][spin][k2][k1][k3] - Residual_ReRhok[SCFi+1][spin][k2][k1][k3];
	        im3 = Residual_ImRhok[SCFi][spin][k2][k1][k3] - Residual_ImRhok[SCFi+1][spin][k2][k1][k3];

		tmp1 += IA[SCFi-1][SCFj-1]*(re1*re2 - Re_DiaG[spin][k2][k1][k3]*re3*re2);
		tmp2 += IA[SCFi-1][SCFj-1]*(im1*im2 - Im_DiaG[spin][k2][k1][k3]*im3*im2);
	      }
	    }

            Re_DiaG[spin][k2][k1][k3] = tmp1;
            Im_DiaG[spin][k2][k1][k3] = tmp2;

            if ( 1.1<Kerker_weight[k2][k1][k3] ){ 
              printf("spin=%2d k2=%2d k1=%2d k3=%2d %10.6f %10.6f %10.6f\n",
                      spin,k2,k1,k3,-1.0/Kerker_weight[k2][k1][k3],Re_DiaG[spin][k2][k1][k3],Im_DiaG[spin][k2][k1][k3]);
	    }

	  }
	}
      }
    }
    */

    /****************************************************
                 calculate the next input
    ****************************************************/

    /* |x(l)> -> |x(l+1)>  */

    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
	for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){
	    ReRhok[0][spin][k2][k1][k3] = ReRhok[1][spin][k2][k1][k3];
	    ImRhok[0][spin][k2][k1][k3] = ImRhok[1][spin][k2][k1][k3];
	  }
	}
      }
    }
    
    /* |x(l+1)> - G|F(l)>, it looks the Kerker mixing */

    if      (1.0e-4<=NormRD[0])                        Mix_wgt = 0.5;
    else if (1.0e-6<=NormRD[0]  && NormRD[0]<1.0e-4)   Mix_wgt = 0.6;
    else if (1.0e-8<=NormRD[0]  && NormRD[0]<1.0e-6)   Mix_wgt = 0.7;
    else if (1.0e-10<=NormRD[0] && NormRD[0]<1.0e-8)   Mix_wgt = 0.8;
    else if (1.0e-12<=NormRD[0] && NormRD[0]<1.0e-10)  Mix_wgt = 0.9;
    else                                               Mix_wgt = 1.0;
    
    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
	for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){

	    weight = 1.0/Kerker_weight[k2][k1][k3];
	    wgt0  = Mix_wgt*weight;

	    ReRhok[0][spin][k2][k1][k3] += wgt0*Residual_ReRhok[1][spin][k2][k1][k3];
	    ImRhok[0][spin][k2][k1][k3] += wgt0*Residual_ImRhok[1][spin][k2][k1][k3];


	    /*
            if (0){
	    weight = 1.0/Kerker_weight[k2][k1][k3];
	    wgt0  = Mix_wgt*weight;

	    ReRhok[0][spin][k2][k1][k3] += wgt0*Residual_ReRhok[1][spin][k2][k1][k3];
	    ImRhok[0][spin][k2][k1][k3] += wgt0*Residual_ImRhok[1][spin][k2][k1][k3];
	    }
            else{
	    ReRhok[0][spin][k2][k1][k3] += -Re_DiaG[spin][k2][k1][k3]*Residual_ReRhok[1][spin][k2][k1][k3];
	    ImRhok[0][spin][k2][k1][k3] += -Im_DiaG[spin][k2][k1][k3]*Residual_ImRhok[1][spin][k2][k1][k3];
	    }
	    */

	  }
	}
      }
    }
  
    /* sum of the remaining contributions */
  
    for (SCFi=1; SCFi<NumMix; SCFi++){
      for (spin=0; spin<spinmax; spin++){
	for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){

	      weight = 1.0/Kerker_weight[k2][k1][k3];
	      wgt0  = Mix_wgt*weight;
              
	      re1 = ReRhok[SCFi][spin][k2][k1][k3] - ReRhok[SCFi+1][spin][k2][k1][k3];
	      im1 = ImRhok[SCFi][spin][k2][k1][k3] - ImRhok[SCFi+1][spin][k2][k1][k3];

	      re2 = Residual_ReRhok[SCFi][spin][k2][k1][k3] - Residual_ReRhok[SCFi+1][spin][k2][k1][k3];
	      im2 = Residual_ImRhok[SCFi][spin][k2][k1][k3] - Residual_ImRhok[SCFi+1][spin][k2][k1][k3];

	      ReRhok[0][spin][k2][k1][k3] -= (Gamma[SCFi] + Alpha[SCFi])*(re1 + wgt0*re2);
	      ImRhok[0][spin][k2][k1][k3] -= (Gamma[SCFi] + Alpha[SCFi])*(im1 + wgt0*im2);

	      /*
	      if(0){
	      weight = 1.0/Kerker_weight[k2][k1][k3];
	      wgt0  = Mix_wgt*weight;
              
	      re1 = ReRhok[SCFi][spin][k2][k1][k3] - ReRhok[SCFi+1][spin][k2][k1][k3];
	      im1 = ImRhok[SCFi][spin][k2][k1][k3] - ImRhok[SCFi+1][spin][k2][k1][k3];

	      re2 = Residual_ReRhok[SCFi][spin][k2][k1][k3] - Residual_ReRhok[SCFi+1][spin][k2][k1][k3];
	      im2 = Residual_ImRhok[SCFi][spin][k2][k1][k3] - Residual_ImRhok[SCFi+1][spin][k2][k1][k3];

	      ReRhok[0][spin][k2][k1][k3] -= (Gamma[SCFi] + Alpha[SCFi])*(re1 + wgt0*re2);
	      ImRhok[0][spin][k2][k1][k3] -= (Gamma[SCFi] + Alpha[SCFi])*(im1 + wgt0*im2);
	      }
              else{
	      re1 = ReRhok[SCFi][spin][k2][k1][k3] - ReRhok[SCFi+1][spin][k2][k1][k3];
	      im1 = ImRhok[SCFi][spin][k2][k1][k3] - ImRhok[SCFi+1][spin][k2][k1][k3];

	      re2 = Residual_ReRhok[SCFi][spin][k2][k1][k3] - Residual_ReRhok[SCFi+1][spin][k2][k1][k3];
	      im2 = Residual_ImRhok[SCFi][spin][k2][k1][k3] - Residual_ImRhok[SCFi+1][spin][k2][k1][k3];

	      ReRhok[0][spin][k2][k1][k3] -= Gamma[SCFi]*(re1 - Re_DiaG[spin][k2][k1][k3]*re2);
	      ImRhok[0][spin][k2][k1][k3] -= Gamma[SCFi]*(im1 - Im_DiaG[spin][k2][k1][k3]*im2);
	      }
	      */

	    }
	  }
	}
      }
    }
  
    /****************************************************
        shift of rho and residual rho
    ****************************************************/
  
    for (pSCF_iter=(NumMix+1); 0<pSCF_iter; pSCF_iter--){
      for (spin=0; spin<spinmax; spin++){
	for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){
	      ReRhok[pSCF_iter][spin][k2][k1][k3] = ReRhok[pSCF_iter-1][spin][k2][k1][k3]; 
	      ImRhok[pSCF_iter][spin][k2][k1][k3] = ImRhok[pSCF_iter-1][spin][k2][k1][k3]; 
	      Residual_ReRhok[pSCF_iter][spin][k2][k1][k3] = Residual_ReRhok[pSCF_iter-1][spin][k2][k1][k3];
	      Residual_ImRhok[pSCF_iter][spin][k2][k1][k3] = Residual_ImRhok[pSCF_iter-1][spin][k2][k1][k3]; 
	    }
	  }
	}
      }
    }
  }

  /****************************************************
        find the charge density in real space 
  ****************************************************/

  tmp0 = 1.0/(double)Ngrid1/(double)Ngrid2/(double)Ngrid3;

  for (spin=0; spin<spinmax; spin++){

    for (k2=0; k2<My_NGrid2_Poisson; k2++){
      for (k1=0; k1<Ngrid1; k1++){
        for (k3=0; k3<Ngrid3; k3++){
          ReV2[k2][k1][k3] = ReRhok[0][spin][k2][k1][k3];
          ImV2[k2][k1][k3] = ImRhok[0][spin][k2][k1][k3];
        }
      }
    }

    if (spin==0 || spin==1){
      Get_Value_inReal(0,ReV2,ImV2,ReV1,ImV1,Density_Grid[spin],Density_Grid[spin]);

      for (MN=0; MN<My_NumGrid1; MN++){
        Density_Grid[spin][MN] = Density_Grid[spin][MN]*tmp0;
      }
    }

    else if (spin==2){
      Get_Value_inReal(1,ReV2,ImV2,ReV1,ImV1,Density_Grid[2],Density_Grid[3]);
      for (MN=0; MN<My_NumGrid1; MN++){
        Density_Grid[2][MN] = Density_Grid[2][MN]*tmp0;
        Density_Grid[3][MN] = Density_Grid[3][MN]*tmp0;
      }
    }

  }

  if (SpinP_switch==0){
    for (MN=0; MN<My_NumGrid1; MN++){
      Density_Grid[1][MN] = Density_Grid[0][MN];
    }
  }

  /****************************************************
     set ReV2 and ImV2 which are used in Poisson.c
  ****************************************************/

  if (SpinP_switch==0){
    for (k2=0; k2<My_NGrid2_Poisson; k2++){
      for (k1=0; k1<Ngrid1; k1++){
        for (k3=0; k3<Ngrid3; k3++){
          ReV2[k2][k1][k3] = 2.0*ReRhok[0][0][k2][k1][k3] - ReRhoAtomk[k2][k1][k3];
          ImV2[k2][k1][k3] = 2.0*ImRhok[0][0][k2][k1][k3] - ImRhoAtomk[k2][k1][k3];
        }
      }
    }
  }
  
  else {
    for (k2=0; k2<My_NGrid2_Poisson; k2++){
      for (k1=0; k1<Ngrid1; k1++){
        for (k3=0; k3<Ngrid3; k3++){
          ReV2[k2][k1][k3] = ReRhok[0][0][k2][k1][k3] + ReRhok[0][1][k2][k1][k3] - ReRhoAtomk[k2][k1][k3];
          ImV2[k2][k1][k3] = ImRhok[0][0][k2][k1][k3] + ImRhok[0][1][k2][k1][k3] - ImRhoAtomk[k2][k1][k3];
        }
      }
    }
  }

  /****************************************************
    freeing of arrays:

    double alden[List_YOUSO[38]];
    double A[List_YOUSO[38]][List_YOUSO[38]];
    double IA[List_YOUSO[38]][List_YOUSO[38]];
    double My_A[List_YOUSO[38]][List_YOUSO[38]];
    double Re_OptRhok[spinmax]
                            [My_NGrid2_Poisson]
                            [Ngrid1]
                            [Ngrid3];
    double Im_OptRhok[spinmax]
                            [My_NGrid2_Poisson]
                            [Ngrid1]
                            [Ngrid3];
  ****************************************************/
  
  free(alden);

  for (i=0; i<List_YOUSO[38]; i++){
    free(A[i]);
  }
  free(A);

  for (i=0; i<List_YOUSO[38]; i++){
    free(My_A[i]);
  }
  free(My_A);

  for (i=0; i<List_YOUSO[38]; i++){
    free(IA[i]);
  }
  free(IA);

  free(My_B);
  free(B);
  free(Gamma);

  free(B1);
  free(My_B1);
  free(Alpha);

  /* Re_OptRhok and Im_OptRhok */

  for (spin=0; spin<spinmax; spin++){
    for (i=0; i<My_NGrid2_Poisson; i++){
      for (j=0; j<Ngrid1; j++){
	free(Re_OptRhok[spin][i][j]);
      }
      free(Re_OptRhok[spin][i]);
    }
    free(Re_OptRhok[spin]);
  }
  free(Re_OptRhok);

  for (spin=0; spin<spinmax; spin++){
    for (i=0; i<My_NGrid2_Poisson; i++){
      for (j=0; j<Ngrid1; j++){
	free(Im_OptRhok[spin][i][j]);
      }
      free(Im_OptRhok[spin][i]);
    }
    free(Im_OptRhok[spin]);
  }
  free(Im_OptRhok);

  /* Kerker_weight */

  for (i=0; i<My_NGrid2_Poisson; i++){
    for (j=0; j<Ngrid1; j++){
      free(Kerker_weight[i][j]);
    }
    free(Kerker_weight[i]);
  }
  free(Kerker_weight);

}






void Broyden1(int SCF_iter,
	     double Mix_wgt,
	     double *****ReRhok,
	     double *****ImRhok,
	     double ****ReBestRhok,
	     double ****ImBestRhok,
	     double *****Residual_ReRhok,
	     double *****Residual_ImRhok,
	     double ***ReV1,
	     double ***ImV1,
	     double ***ReV2,
	     double ***ImV2,
	     double ***ReRhoAtomk,
	     double ***ImRhoAtomk)
{
  static int firsttime=1;
  int spin,Mc_AN,Gc_AN,wan1,TNO1,h_AN,Gh_AN;
  int wan2,TNO2,i,j,k,NumMix;
  int SCFi,SCFj,tno1,tno0,Cwan,Hwan,k1,k2,k3,kk2; 
  int pSCF_iter,size_Re_OptRhok,size_Kerker_weight;
  int spinmax,MN,flag_nan;
  double *alden;
  double My_Norm,My_OptNorm,OptNorm;
  double **A,**My_A,**IA,*B,*My_B,*Gamma;
  double *B1,*My_B1,*Alpha;
  double tmp0,itmp0,tmp1,im1,im2,re1,re2;
  double OptRNorm,My_OptRNorm;
  double dum1,dum2,bunsi,bunbo,sum,coef_OptRho;
  double Gx,Gy,Gz,G2,G12,G22,G32;
  double Min_Weight;
  double Max_Weight;
  double G0,G02,G02p,weight,wgt0,wgt1;
  double sk1,sk2,sk3;
  double ****Re_OptRhok;
  double ****Im_OptRhok;
  double ***Kerker_weight;
  int numprocs,myid;
  char nanchar[300];

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /* find an optimum G0 */

  G12 = rtv[1][1]*rtv[1][1] + rtv[1][2]*rtv[1][2] + rtv[1][3]*rtv[1][3]; 
  G22 = rtv[2][1]*rtv[2][1] + rtv[2][2]*rtv[2][2] + rtv[2][3]*rtv[2][3]; 
  G32 = rtv[3][1]*rtv[3][1] + rtv[3][2]*rtv[3][2] + rtv[3][3]*rtv[3][3]; 

  if (G12<G22) G0 = G12;
  else         G0 = G22;
  if (G32<G0)  G0 = G32;

  G0 = Kerker_factor*sqrt(G0);
  G02 = G0*G0;
  G02p = (0.1*G0)*(0.1*G0);

  if      (SpinP_switch==0)  spinmax = 1;
  else if (SpinP_switch==1)  spinmax = 2;
  else if (SpinP_switch==3)  spinmax = 3;

  /****************************************************
    allocation of arrays:

    double alden[List_YOUSO[38]];
    double A[List_YOUSO[38]][List_YOUSO[38]];
    double IA[List_YOUSO[38]][List_YOUSO[38]];
    double My_A[List_YOUSO[38]][List_YOUSO[38]];
    double Re_OptRhok[spinmax]
                            [My_NGrid2_Poisson]
                            [Ngrid1]
                            [Ngrid3];
    double Im_OptRhok[spinmax]
                            [My_NGrid2_Poisson]
                            [Ngrid1]
                            [Ngrid3];
  ****************************************************/
  
  alden = (double*)malloc(sizeof(double)*List_YOUSO[38]);

  A = (double**)malloc(sizeof(double*)*List_YOUSO[38]);
  for (i=0; i<List_YOUSO[38]; i++){
    A[i] = (double*)malloc(sizeof(double)*List_YOUSO[38]);
  }

  My_A = (double**)malloc(sizeof(double*)*List_YOUSO[38]);
  for (i=0; i<List_YOUSO[38]; i++){
    My_A[i] = (double*)malloc(sizeof(double)*List_YOUSO[38]);
  }

  IA = (double**)malloc(sizeof(double*)*List_YOUSO[38]);
  for (i=0; i<List_YOUSO[38]; i++){
    IA[i] = (double*)malloc(sizeof(double)*List_YOUSO[38]);
    for (j=0; j<List_YOUSO[38]; j++){
      IA[i][j] = 0.0;
    }
  }

  My_B = (double*)malloc(sizeof(double)*List_YOUSO[38]);
  B = (double*)malloc(sizeof(double)*List_YOUSO[38]);
  Gamma = (double*)malloc(sizeof(double)*List_YOUSO[38]);

  B1 = (double*)malloc(sizeof(double)*List_YOUSO[38]);
  My_B1 = (double*)malloc(sizeof(double)*List_YOUSO[38]);
  Alpha = (double*)malloc(sizeof(double)*List_YOUSO[38]);

  /* Re_OptRhok and Im_OptRhok */

  size_Re_OptRhok = spinmax*My_NGrid2_Poisson*Ngrid1*Ngrid3;

  Re_OptRhok = (double****)malloc(sizeof(double***)*spinmax); 
  for (spin=0; spin<spinmax; spin++){
    Re_OptRhok[spin] = (double***)malloc(sizeof(double**)*My_NGrid2_Poisson); 
    for (i=0; i<My_NGrid2_Poisson; i++){
      Re_OptRhok[spin][i] = (double**)malloc(sizeof(double*)*Ngrid1); 
      for (j=0; j<Ngrid1; j++){
	Re_OptRhok[spin][i][j] = (double*)malloc(sizeof(double)*Ngrid3); 
	for (k=0; k<Ngrid3; k++) Re_OptRhok[spin][i][j][k] = 0.0;
      }
    }
  }

  Im_OptRhok = (double****)malloc(sizeof(double***)*spinmax); 
  for (spin=0; spin<spinmax; spin++){
    Im_OptRhok[spin] = (double***)malloc(sizeof(double**)*My_NGrid2_Poisson); 
    for (i=0; i<My_NGrid2_Poisson; i++){
      Im_OptRhok[spin][i] = (double**)malloc(sizeof(double*)*Ngrid1); 
      for (j=0; j<Ngrid1; j++){
	Im_OptRhok[spin][i][j] = (double*)malloc(sizeof(double)*Ngrid3); 
	for (k=0; k<Ngrid3; k++) Im_OptRhok[spin][i][j][k] = 0.0;
      }
    }
  }

  /* Kerker_weight */

  size_Kerker_weight = My_NGrid2_Poisson*Ngrid1*Ngrid3;
  
  Kerker_weight = (double***)malloc(sizeof(double**)*My_NGrid2_Poisson); 
  for (i=0; i<My_NGrid2_Poisson; i++){
    Kerker_weight[i] = (double**)malloc(sizeof(double*)*Ngrid1); 
    for (j=0; j<Ngrid1; j++){
      Kerker_weight[i][j] = (double*)malloc(sizeof(double)*Ngrid3); 
    }
  }

  for (k2=0; k2<My_NGrid2_Poisson; k2++){

    kk2 = k2 + Start_Grid2[myid];

    if (kk2<Ngrid2/2) sk2 = (double)kk2;
    else              sk2 = (double)(kk2 - Ngrid2);

    for (k1=0; k1<Ngrid1; k1++){

      if (k1<Ngrid1/2) sk1 = (double)k1;
      else             sk1 = (double)(k1 - Ngrid1);

      for (k3=0; k3<Ngrid3; k3++){

	if (k3<Ngrid3/2) sk3 = (double)k3;
	else             sk3 = (double)(k3 - Ngrid3);

	Gx = sk1*rtv[1][1] + sk2*rtv[2][1] + sk3*rtv[3][1];
	Gy = sk1*rtv[1][2] + sk2*rtv[2][2] + sk3*rtv[3][2]; 
	Gz = sk1*rtv[1][3] + sk2*rtv[2][3] + sk3*rtv[3][3];
	G2 = Gx*Gx + Gy*Gy + Gz*Gz;

	if (k1==0 && kk2==0 && k3==0)  weight = 1.0;
	else                           weight = (G2 + G02)/(G2 + G02p);

        Kerker_weight[k2][k1][k3] = weight;
      }
    }
  }

  /* memory calc. done, only 1st iteration */
  if (firsttime) {
    PrintMemory("DIIS_Mixing_Rhok: Re_OptRhok",sizeof(double)*size_Re_OptRhok,NULL);
    PrintMemory("DIIS_Mixing_Rhok: Im_OptRhok",sizeof(double)*size_Re_OptRhok,NULL);
    PrintMemory("DIIS_Mixing_Rhok: Kerker_weight",sizeof(double)*size_Kerker_weight,NULL);
    firsttime=0;
  }

  /****************************************************
     start calc.
  ****************************************************/

  /*
  printf("SCF_iter=%2d\n",SCF_iter);fflush(stdout);
  */

  SCF_iter--;
  NumMix = SCF_iter%(Num_Mixing_pDM+1);
  if (NumMix<2) NumMix = 2;

  /*
  if (Num_Mixing_pDM<NumMix){
    NumMix = Num_Mixing_pDM%(Num_Mixing_pDM+1);
  }
  */

  if (SCF_iter<=1){

    /* construct residual rho */

    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
        for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){
	    tmp0 = ReRhok[0][spin][k2][k1][k3] - ReRhok[1][spin][k2][k1][k3];
	    tmp1 = ImRhok[0][spin][k2][k1][k3] - ImRhok[1][spin][k2][k1][k3];
	    Residual_ReRhok[2][spin][k2][k1][k3] = tmp0;
	    Residual_ImRhok[2][spin][k2][k1][k3] = tmp1;
	  }
	}
      }
    }

    /* rho1 to rho2 */

    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
        for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){
	    ReRhok[2][spin][k2][k1][k3] = ReRhok[1][spin][k2][k1][k3];
	    ImRhok[2][spin][k2][k1][k3] = ImRhok[1][spin][k2][k1][k3];
	  }
	}
      }
    }

    Kerker_Mixing_Rhok(1,Mixing_weight,ReRhok,ImRhok,ReBestRhok,ImBestRhok,
                       Residual_ReRhok,Residual_ImRhok,
                       ReV1,ImV1,ReV2,ImV2,ReRhoAtomk,ImRhoAtomk);
  } 

  else {
 
    /****************************************************
       construct residual rho1
       residual = output - input  
    ****************************************************/

    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
	for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){
	    tmp0 = ReRhok[0][spin][k2][k1][k3] - ReRhok[1][spin][k2][k1][k3];
	    tmp1 = ImRhok[0][spin][k2][k1][k3] - ImRhok[1][spin][k2][k1][k3];
	    Residual_ReRhok[1][spin][k2][k1][k3] = tmp0;
	    Residual_ImRhok[1][spin][k2][k1][k3] = tmp1;
	  }
	}
      }
    }

    /* calculate NormRD */
  
    My_Norm = 0.0;
    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
	for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){

	    weight = Kerker_weight[k2][k1][k3]; 

	    re1 = Residual_ReRhok[1][spin][k2][k1][k3];
	    im1 = Residual_ImRhok[1][spin][k2][k1][k3];

	    My_Norm += (re1*re1 + im1*im1)*weight;
	  }
	}
      }
    }

    MPI_Allreduce(&My_Norm, &NormRD[0], 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

    for (i=4; 1<=i; i--){
      NormRD[i] = NormRD[i-1];
      History_Uele[i] = History_Uele[i-1];
    }
    NormRD[0] = NormRD[0]/(double)Ngrid1/(double)Ngrid2/(double)Ngrid3/(double)spinmax;
    History_Uele[0] = Uele;

    /****************************************************
     The matrix A calculated by difference residual rho
    ****************************************************/
      
    for (SCFi=1; SCFi<NumMix; SCFi++){
      for (SCFj=SCFi; SCFj<NumMix; SCFj++){

	My_A[SCFi][SCFj] = 0.0;

	for (spin=0; spin<spinmax; spin++){
	  for (k2=0; k2<My_NGrid2_Poisson; k2++){
	    for (k1=0; k1<Ngrid1; k1++){
	      for (k3=0; k3<Ngrid3; k3++){

		weight = Kerker_weight[k2][k1][k3]; 

		re1 = Residual_ReRhok[SCFi][spin][k2][k1][k3] - Residual_ReRhok[SCFi+1][spin][k2][k1][k3];
		im1 = Residual_ImRhok[SCFi][spin][k2][k1][k3] - Residual_ImRhok[SCFi+1][spin][k2][k1][k3];

		re2 = Residual_ReRhok[SCFj][spin][k2][k1][k3] - Residual_ReRhok[SCFj+1][spin][k2][k1][k3];
		im2 = Residual_ImRhok[SCFj][spin][k2][k1][k3] - Residual_ImRhok[SCFj+1][spin][k2][k1][k3];

		My_A[SCFi][SCFj] += (re1*re2 + im1*im2)*weight;
	      }
	    }
	  }
	}
      }
    }

    /* MPI My_A */

    for (SCFi=1; SCFi<NumMix; SCFi++){
      for (SCFj=SCFi; SCFj<NumMix; SCFj++){
	MPI_Allreduce(&My_A[SCFi][SCFj], &A[SCFi-1][SCFj-1],
		      1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
	A[SCFj-1][SCFi-1] = A[SCFi-1][SCFj-1];
      }
    }

    /*
    printf("A myid=%2d SCF_iter=%2d\n",myid,SCF_iter);fflush(stdout);
    for (SCFi=1; SCFi<NumMix; SCFi++){
      for (SCFj=1; SCFj<NumMix; SCFj++){
	printf("%18.13f ",A[SCFi-1][SCFj-1]);fflush(stdout);  
      }
      printf("\n");fflush(stdout);
    }
    */

    /* calculate <delta F(k)|F(l)> */

    for (SCFi=1; SCFi<NumMix; SCFi++){

      My_B[SCFi-1] = 0.0;

      for (spin=0; spin<spinmax; spin++){
	for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){

	      weight = Kerker_weight[k2][k1][k3]; 

	      re1 = Residual_ReRhok[SCFi][spin][k2][k1][k3] - Residual_ReRhok[SCFi+1][spin][k2][k1][k3];
	      im1 = Residual_ImRhok[SCFi][spin][k2][k1][k3] - Residual_ImRhok[SCFi+1][spin][k2][k1][k3];

	      re2 = Residual_ReRhok[1][spin][k2][k1][k3];
	      im2 = Residual_ImRhok[1][spin][k2][k1][k3];

	      My_B[SCFi-1] += (re1*re2 + im1*im2)*weight; 
	    }
	  }
	}
      }      
    }

    /* MPI My_B */

    for (SCFi=1; SCFi<NumMix; SCFi++){
      MPI_Allreduce(&My_B[SCFi-1], &B[SCFi-1], 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
    }

    /*
    printf("B myid=%2d SCF_iter=%2d\n",myid,SCF_iter);fflush(stdout);
    for (SCFi=1; SCFi<NumMix; SCFi++){
      printf("%18.13f\n",B[SCFi-1]); fflush(stdout);   
    }
    */

    /* calculate the inverse of A */

    if (2<=NumMix) Inverse(NumMix-2,A,IA);

    /*
    printf("IA myid=%2d SCF_iter=%2d\n",myid,SCF_iter);fflush(stdout);
    for (SCFi=1; SCFi<NumMix; SCFi++){
      for (SCFj=1; SCFj<NumMix; SCFj++){
	printf("%18.13f ",IA[SCFi-1][SCFj-1]); fflush(stdout); 
      }
      printf("\n");fflush(stdout);
    }
    */

    /* calculate gamma */
    
    for (SCFi=1; SCFi<NumMix; SCFi++){
      Gamma[SCFi] = 0.0;
      for (SCFj=1; SCFj<NumMix; SCFj++){
	Gamma[SCFi] += IA[SCFi-1][SCFj-1]*B[SCFj-1];
      }
    }

    printf("Gamma myid=%2d SCF_iter=%2d\n",myid,SCF_iter);fflush(stdout);
    for (SCFi=1; SCFi<NumMix; SCFi++){
      printf("%18.13f\n",Gamma[SCFi]);  fflush(stdout);  
    }

    /* calculate the optimum residual vector */

    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
	for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){
	    ReBestRhok[spin][k2][k1][k3] = Residual_ReRhok[1][spin][k2][k1][k3];
	    ImBestRhok[spin][k2][k1][k3] = Residual_ImRhok[1][spin][k2][k1][k3];
	  }
	}
      }
    }

    for (SCFi=1; SCFi<NumMix; SCFi++){
      for (spin=0; spin<spinmax; spin++){
	for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){

	      re2 = Residual_ReRhok[SCFi][spin][k2][k1][k3] - Residual_ReRhok[SCFi+1][spin][k2][k1][k3];
	      im2 = Residual_ImRhok[SCFi][spin][k2][k1][k3] - Residual_ImRhok[SCFi+1][spin][k2][k1][k3];

	      ReBestRhok[spin][k2][k1][k3] -= Gamma[SCFi]*re2;
	      ImBestRhok[spin][k2][k1][k3] -= Gamma[SCFi]*im2;
	    }
	  }
	}
      }
    }

    /* calculate the optimum norm */

    My_OptNorm = 0.0;
    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
	for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){

	    weight = Kerker_weight[k2][k1][k3]; 

	    re1 = ReBestRhok[spin][k2][k1][k3];
	    im1 = ImBestRhok[spin][k2][k1][k3];

	    My_OptNorm += (re1*re1 + im1*im1)*weight;
	  }
	}
      }
    }

    MPI_Allreduce(&My_OptNorm, &OptNorm, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

    /*
    printf("myid=%2d SCF_iter=%2d OptNorm=%15.12f\n",myid,SCF_iter,OptNorm);fflush(stdout);
    */

    /* calculate <delta F(k)| \overbar{F}(l)> */

    for (SCFi=1; SCFi<NumMix; SCFi++){

      My_B1[SCFi-1] = 0.0;

      for (spin=0; spin<spinmax; spin++){
	for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){

	      weight = Kerker_weight[k2][k1][k3]; 

	      re1 = Residual_ReRhok[SCFi][spin][k2][k1][k3] - Residual_ReRhok[SCFi+1][spin][k2][k1][k3];
	      im1 = Residual_ImRhok[SCFi][spin][k2][k1][k3] - Residual_ImRhok[SCFi+1][spin][k2][k1][k3];

	      re2 = ReBestRhok[spin][k2][k1][k3];
	      im2 = ImBestRhok[spin][k2][k1][k3];

	      My_B1[SCFi-1] += (re1*re2 + im1*im2)*weight; 
	    }
	  }
	}
      }      
    }

    /* MPI My_B1 */

    for (SCFi=1; SCFi<NumMix; SCFi++){
      MPI_Allreduce(&My_B1[SCFi-1], &B1[SCFi-1], 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
    }

    /*
    if (myid==0){
    printf("B1 myid=%2d SCF_iter=%2d\n",myid,SCF_iter);fflush(stdout);
    for (SCFi=1; SCFi<NumMix; SCFi++){
      printf("%15.12f\n",B1[SCFi-1]); fflush(stdout);   
    }
    }
    */
    
    /* calculate alpha */
    
    for (SCFi=1; SCFi<NumMix; SCFi++){
      Alpha[SCFi] = 0.0;
      for (SCFj=1; SCFj<NumMix; SCFj++){
	Alpha[SCFi] += IA[SCFi-1][SCFj-1]*B1[SCFj-1];
      }
    }
    
    /*
    printf("Alpha myid=%2d SCF_iter=%2d\n",myid,SCF_iter);fflush(stdout);
    for (SCFi=1; SCFi<NumMix; SCFi++){
      printf("%15.12f\n",Alpha[SCFi]);  fflush(stdout);
    }
    */

    /****************************************************
                 calculate the next input
    ****************************************************/

    /* |x(l)> -> |x(l+1)>  */

    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
	for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){
	    ReRhok[0][spin][k2][k1][k3] = ReRhok[1][spin][k2][k1][k3];
	    ImRhok[0][spin][k2][k1][k3] = ImRhok[1][spin][k2][k1][k3];
	  }
	}
      }
    }
    
    /* |x(l+1)> - G|F(l)>, it looks the Kerker mixing */

    if      (1.0e-4<=NormRD[0])                        Mix_wgt = 0.5;
    else if (1.0e-6<=NormRD[0]  && NormRD[0]<1.0e-4)   Mix_wgt = 0.6;
    else if (1.0e-8<=NormRD[0]  && NormRD[0]<1.0e-6)   Mix_wgt = 0.7;
    else if (1.0e-10<=NormRD[0] && NormRD[0]<1.0e-8)   Mix_wgt = 0.8;
    else if (1.0e-12<=NormRD[0] && NormRD[0]<1.0e-10)  Mix_wgt = 0.9;
    else                                               Mix_wgt = 1.0;
    
    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
	for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){

	    weight = 1.0/Kerker_weight[k2][k1][k3];
	    wgt0  = Mix_wgt*weight;

	    ReRhok[0][spin][k2][k1][k3] += wgt0*Residual_ReRhok[1][spin][k2][k1][k3];
	    ImRhok[0][spin][k2][k1][k3] += wgt0*Residual_ImRhok[1][spin][k2][k1][k3];
	  }
	}
      }
    }
  
    /* sum of the remaining contributions */
  
    for (SCFi=1; SCFi<NumMix; SCFi++){
      for (spin=0; spin<spinmax; spin++){
	for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){

	      weight = 1.0/Kerker_weight[k2][k1][k3];
	      wgt0  = Mix_wgt*weight;
              
	      re1 = ReRhok[SCFi][spin][k2][k1][k3] - ReRhok[SCFi+1][spin][k2][k1][k3];
	      im1 = ImRhok[SCFi][spin][k2][k1][k3] - ImRhok[SCFi+1][spin][k2][k1][k3];

	      re2 = Residual_ReRhok[SCFi][spin][k2][k1][k3] - Residual_ReRhok[SCFi+1][spin][k2][k1][k3];
	      im2 = Residual_ImRhok[SCFi][spin][k2][k1][k3] - Residual_ImRhok[SCFi+1][spin][k2][k1][k3];

	      ReRhok[0][spin][k2][k1][k3] -= (Gamma[SCFi] + Alpha[SCFi])*(re1 + wgt0*re2);
	      ImRhok[0][spin][k2][k1][k3] -= (Gamma[SCFi] + Alpha[SCFi])*(im1 + wgt0*im2);
	    }
	  }
	}
      }
    }
  
    /****************************************************
        shift of rho and residual rho
    ****************************************************/
  
    for (pSCF_iter=(NumMix+1); 0<pSCF_iter; pSCF_iter--){
      for (spin=0; spin<spinmax; spin++){
	for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){
	      ReRhok[pSCF_iter][spin][k2][k1][k3] = ReRhok[pSCF_iter-1][spin][k2][k1][k3]; 
	      ImRhok[pSCF_iter][spin][k2][k1][k3] = ImRhok[pSCF_iter-1][spin][k2][k1][k3]; 
	      Residual_ReRhok[pSCF_iter][spin][k2][k1][k3] = Residual_ReRhok[pSCF_iter-1][spin][k2][k1][k3];
	      Residual_ImRhok[pSCF_iter][spin][k2][k1][k3] = Residual_ImRhok[pSCF_iter-1][spin][k2][k1][k3]; 
	    }
	  }
	}
      }
    }
  }

  /****************************************************
        find the charge density in real space 
  ****************************************************/

  tmp0 = 1.0/(double)Ngrid1/(double)Ngrid2/(double)Ngrid3;

  for (spin=0; spin<spinmax; spin++){

    for (k2=0; k2<My_NGrid2_Poisson; k2++){
      for (k1=0; k1<Ngrid1; k1++){
        for (k3=0; k3<Ngrid3; k3++){
          ReV2[k2][k1][k3] = ReRhok[0][spin][k2][k1][k3];
          ImV2[k2][k1][k3] = ImRhok[0][spin][k2][k1][k3];
        }
      }
    }

    if (spin==0 || spin==1){
      Get_Value_inReal(0,ReV2,ImV2,ReV1,ImV1,Density_Grid[spin],Density_Grid[spin]);

      for (MN=0; MN<My_NumGrid1; MN++){
        Density_Grid[spin][MN] = Density_Grid[spin][MN]*tmp0;
      }
    }

    else if (spin==2){
      Get_Value_inReal(1,ReV2,ImV2,ReV1,ImV1,Density_Grid[2],Density_Grid[3]);
      for (MN=0; MN<My_NumGrid1; MN++){
        Density_Grid[2][MN] = Density_Grid[2][MN]*tmp0;
        Density_Grid[3][MN] = Density_Grid[3][MN]*tmp0;
      }
    }

  }

  if (SpinP_switch==0){
    for (MN=0; MN<My_NumGrid1; MN++){
      Density_Grid[1][MN] = Density_Grid[0][MN];
    }
  }

  /****************************************************
     set ReV2 and ImV2 which are used in Poisson.c
  ****************************************************/

  if (SpinP_switch==0){
    for (k2=0; k2<My_NGrid2_Poisson; k2++){
      for (k1=0; k1<Ngrid1; k1++){
        for (k3=0; k3<Ngrid3; k3++){
          ReV2[k2][k1][k3] = 2.0*ReRhok[0][0][k2][k1][k3] - ReRhoAtomk[k2][k1][k3];
          ImV2[k2][k1][k3] = 2.0*ImRhok[0][0][k2][k1][k3] - ImRhoAtomk[k2][k1][k3];
        }
      }
    }
  }
  
  else {
    for (k2=0; k2<My_NGrid2_Poisson; k2++){
      for (k1=0; k1<Ngrid1; k1++){
        for (k3=0; k3<Ngrid3; k3++){
          ReV2[k2][k1][k3] = ReRhok[0][0][k2][k1][k3] + ReRhok[0][1][k2][k1][k3] - ReRhoAtomk[k2][k1][k3];
          ImV2[k2][k1][k3] = ImRhok[0][0][k2][k1][k3] + ImRhok[0][1][k2][k1][k3] - ImRhoAtomk[k2][k1][k3];
        }
      }
    }
  }

  /****************************************************
    freeing of arrays:

    double alden[List_YOUSO[38]];
    double A[List_YOUSO[38]][List_YOUSO[38]];
    double IA[List_YOUSO[38]][List_YOUSO[38]];
    double My_A[List_YOUSO[38]][List_YOUSO[38]];
    double Re_OptRhok[spinmax]
                            [My_NGrid2_Poisson]
                            [Ngrid1]
                            [Ngrid3];
    double Im_OptRhok[spinmax]
                            [My_NGrid2_Poisson]
                            [Ngrid1]
                            [Ngrid3];
  ****************************************************/
  
  free(alden);

  for (i=0; i<List_YOUSO[38]; i++){
    free(A[i]);
  }
  free(A);

  for (i=0; i<List_YOUSO[38]; i++){
    free(My_A[i]);
  }
  free(My_A);

  for (i=0; i<List_YOUSO[38]; i++){
    free(IA[i]);
  }
  free(IA);

  free(My_B);
  free(B);
  free(Gamma);

  free(B1);
  free(My_B1);
  free(Alpha);

  /* Re_OptRhok and Im_OptRhok */

  for (spin=0; spin<spinmax; spin++){
    for (i=0; i<My_NGrid2_Poisson; i++){
      for (j=0; j<Ngrid1; j++){
	free(Re_OptRhok[spin][i][j]);
      }
      free(Re_OptRhok[spin][i]);
    }
    free(Re_OptRhok[spin]);
  }
  free(Re_OptRhok);

  for (spin=0; spin<spinmax; spin++){
    for (i=0; i<My_NGrid2_Poisson; i++){
      for (j=0; j<Ngrid1; j++){
	free(Im_OptRhok[spin][i][j]);
      }
      free(Im_OptRhok[spin][i]);
    }
    free(Im_OptRhok[spin]);
  }
  free(Im_OptRhok);

  /* Kerker_weight */

  for (i=0; i<My_NGrid2_Poisson; i++){
    for (j=0; j<Ngrid1; j++){
      free(Kerker_weight[i][j]);
    }
    free(Kerker_weight[i]);
  }
  free(Kerker_weight);

}







void DIIS_Mixing_Rhok0(int SCF_iter,
                      double Mix_wgt,
                      double *****ReRhok,
                      double *****ImRhok,
                      double ****ReBestRhok,
                      double ****ImBestRhok,
                      double *****Residual_ReRhok,
                      double *****Residual_ImRhok,
                      double ***ReV1,
                      double ***ImV1,
                      double ***ReV2,
                      double ***ImV2,
                      double ***ReRhoAtomk,
                      double ***ImRhoAtomk)
{
  static int firsttime=1;
  int spin,Mc_AN,Gc_AN,wan1,TNO1,h_AN,Gh_AN;
  int wan2,TNO2,i,j,k,NumMix,NumSlide;
  int SCFi,SCFj,tno1,tno0,Cwan,Hwan,k1,k2,k3,kk2; 
  int pSCF_iter,size_Re_OptRhok,size_Kerker_weight;
  int spinmax,MN,flag_nan;
  double *alden;
  double **A,**My_A,**IA;
  double tmp0,itmp0,tmp1,im1,im2,re1,re2;
  double OptRNorm,My_OptRNorm;
  double dum1,dum2,bunsi,bunbo,sum,coef_OptRho;
  double Gx,Gy,Gz,G2,G12,G22,G32;
  double Min_Weight;
  double Max_Weight;
  double G0,G02,G02p,weight,wgt0,wgt1;
  double sk1,sk2,sk3;
  double ****Re_OptRhok;
  double ****Im_OptRhok;
  double ***Kerker_weight;
  int numprocs,myid;
  char nanchar[300];

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /* find an optimum G0 */

  G12 = rtv[1][1]*rtv[1][1] + rtv[1][2]*rtv[1][2] + rtv[1][3]*rtv[1][3]; 
  G22 = rtv[2][1]*rtv[2][1] + rtv[2][2]*rtv[2][2] + rtv[2][3]*rtv[2][3]; 
  G32 = rtv[3][1]*rtv[3][1] + rtv[3][2]*rtv[3][2] + rtv[3][3]*rtv[3][3]; 

  if (G12<G22) G0 = G12;
  else         G0 = G22;
  if (G32<G0)  G0 = G32;

  G0 = Kerker_factor*sqrt(G0);
  G02 = G0*G0;
  G02p = (0.01*G0)*(0.01*G0);

  if      (SpinP_switch==0)  spinmax = 1;
  else if (SpinP_switch==1)  spinmax = 2;
  else if (SpinP_switch==3)  spinmax = 3;

  /****************************************************
    allocation of arrays:

    double alden[List_YOUSO[38]];
    double A[List_YOUSO[38]][List_YOUSO[38]];
    double IA[List_YOUSO[38]][List_YOUSO[38]];
    double My_A[List_YOUSO[38]][List_YOUSO[38]];
    double Re_OptRhok[spinmax]
                            [My_NGrid2_Poisson]
                            [Ngrid1]
                            [Ngrid3];
    double Im_OptRhok[spinmax]
                            [My_NGrid2_Poisson]
                            [Ngrid1]
                            [Ngrid3];
  ****************************************************/
  
  alden = (double*)malloc(sizeof(double)*List_YOUSO[38]);

  A = (double**)malloc(sizeof(double*)*List_YOUSO[38]);
  for (i=0; i<List_YOUSO[38]; i++){
    A[i] = (double*)malloc(sizeof(double)*List_YOUSO[38]);
  }

  My_A = (double**)malloc(sizeof(double*)*List_YOUSO[38]);
  for (i=0; i<List_YOUSO[38]; i++){
    My_A[i] = (double*)malloc(sizeof(double)*List_YOUSO[38]);
  }

  IA = (double**)malloc(sizeof(double*)*List_YOUSO[38]);
  for (i=0; i<List_YOUSO[38]; i++){
    IA[i] = (double*)malloc(sizeof(double)*List_YOUSO[38]);
  }

  /* Re_OptRhok and Im_OptRhok */

  size_Re_OptRhok = spinmax*My_NGrid2_Poisson*Ngrid1*Ngrid3;

  Re_OptRhok = (double****)malloc(sizeof(double***)*spinmax); 
  for (spin=0; spin<spinmax; spin++){
    Re_OptRhok[spin] = (double***)malloc(sizeof(double**)*My_NGrid2_Poisson); 
    for (i=0; i<My_NGrid2_Poisson; i++){
      Re_OptRhok[spin][i] = (double**)malloc(sizeof(double*)*Ngrid1); 
      for (j=0; j<Ngrid1; j++){
	Re_OptRhok[spin][i][j] = (double*)malloc(sizeof(double)*Ngrid3); 
	for (k=0; k<Ngrid3; k++) Re_OptRhok[spin][i][j][k] = 0.0;
      }
    }
  }

  Im_OptRhok = (double****)malloc(sizeof(double***)*spinmax); 
  for (spin=0; spin<spinmax; spin++){
    Im_OptRhok[spin] = (double***)malloc(sizeof(double**)*My_NGrid2_Poisson); 
    for (i=0; i<My_NGrid2_Poisson; i++){
      Im_OptRhok[spin][i] = (double**)malloc(sizeof(double*)*Ngrid1); 
      for (j=0; j<Ngrid1; j++){
	Im_OptRhok[spin][i][j] = (double*)malloc(sizeof(double)*Ngrid3); 
	for (k=0; k<Ngrid3; k++) Im_OptRhok[spin][i][j][k] = 0.0;
      }
    }
  }

  /* Kerker_weight */

  size_Kerker_weight = My_NGrid2_Poisson*Ngrid1*Ngrid3;

  Kerker_weight = (double***)malloc(sizeof(double**)*My_NGrid2_Poisson); 
  for (i=0; i<My_NGrid2_Poisson; i++){
    Kerker_weight[i] = (double**)malloc(sizeof(double*)*Ngrid1); 
    for (j=0; j<Ngrid1; j++){
      Kerker_weight[i][j] = (double*)malloc(sizeof(double)*Ngrid3); 
    }
  }

  for (k2=0; k2<My_NGrid2_Poisson; k2++){

    kk2 = k2 + Start_Grid2[myid];

    if (kk2<Ngrid2/2) sk2 = (double)kk2;
    else              sk2 = (double)(kk2 - Ngrid2);

    for (k1=0; k1<Ngrid1; k1++){

      if (k1<Ngrid1/2) sk1 = (double)k1;
      else             sk1 = (double)(k1 - Ngrid1);

      for (k3=0; k3<Ngrid3; k3++){

	if (k3<Ngrid3/2) sk3 = (double)k3;
	else             sk3 = (double)(k3 - Ngrid3);

	Gx = sk1*rtv[1][1] + sk2*rtv[2][1] + sk3*rtv[3][1];
	Gy = sk1*rtv[1][2] + sk2*rtv[2][2] + sk3*rtv[3][2]; 
	Gz = sk1*rtv[1][3] + sk2*rtv[2][3] + sk3*rtv[3][3];
	G2 = Gx*Gx + Gy*Gy + Gz*Gz;

	if (k1==0 && kk2==0 && k3==0)  weight = 1.0;
	else                           weight = (G2 + G02)/(G2 + G02p);

        Kerker_weight[k2][k1][k3] = weight;
      }
    }
  }

  /****************************************************
     start calc.
  ****************************************************/
 
  if (SCF_iter==1){
    /* memory calc. done, only 1st iteration */
    if (firsttime) {
      PrintMemory("DIIS_Mixing_Rhok: Re_OptRhok",sizeof(double)*size_Re_OptRhok,NULL);
      PrintMemory("DIIS_Mixing_Rhok: Im_OptRhok",sizeof(double)*size_Re_OptRhok,NULL);
      PrintMemory("DIIS_Mixing_Rhok: Kerker_weight",sizeof(double)*size_Kerker_weight,NULL);
      firsttime=0;
    }

    Kerker_Mixing_Rhok(0,1.00,ReRhok,ImRhok,ReBestRhok,ImBestRhok,
                       Residual_ReRhok,Residual_ImRhok,
                       ReV1,ImV1,ReV2,ImV2,ReRhoAtomk,ImRhoAtomk);
  } 

  else if (SCF_iter==2){

    /* construct residual rho */

    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
        for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){
	    tmp0 = ReRhok[0][spin][k2][k1][k3] - ReRhok[1][spin][k2][k1][k3];
	    tmp1 = ImRhok[0][spin][k2][k1][k3] - ImRhok[1][spin][k2][k1][k3];
	    Residual_ReRhok[2][spin][k2][k1][k3] = tmp0;
	    Residual_ImRhok[2][spin][k2][k1][k3] = tmp1;
	  }
	}
      }
    }

    /* rho1 to rho2 */

    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
        for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){
	    ReRhok[2][spin][k2][k1][k3] = ReRhok[1][spin][k2][k1][k3];
	    ImRhok[2][spin][k2][k1][k3] = ImRhok[1][spin][k2][k1][k3];
	  }
	}
      }
    }

    Kerker_Mixing_Rhok(1,Mixing_weight,ReRhok,ImRhok,ReBestRhok,ImBestRhok,
                       Residual_ReRhok,Residual_ImRhok,
                       ReV1,ImV1,ReV2,ImV2,ReRhoAtomk,ImRhoAtomk);
  } 

  else {

    /****************************************************
                   construct residual rho1
    ****************************************************/

    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
        for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){
	    tmp0 = ReRhok[0][spin][k2][k1][k3] - ReRhok[1][spin][k2][k1][k3];
	    tmp1 = ImRhok[0][spin][k2][k1][k3] - ImRhok[1][spin][k2][k1][k3];
	    Residual_ReRhok[1][spin][k2][k1][k3] = tmp0;
	    Residual_ImRhok[1][spin][k2][k1][k3] = tmp1;
	  }
	}
      }
    }

    if ((SCF_iter-1)<Num_Mixing_pDM){
      NumMix   = SCF_iter - 1;
      NumSlide = NumMix + 1;
    }
    else{
      NumMix   = Num_Mixing_pDM;
      NumSlide = NumMix;
    }

    /****************************************************
                    alpha from residual rho
    ****************************************************/
      
    for (SCFi=1; SCFi<=NumMix; SCFi++){
      for (SCFj=SCFi; SCFj<=NumMix; SCFj++){

        My_A[SCFi][SCFj]  = 0.0;

        for (spin=0; spin<spinmax; spin++){
          for (k2=0; k2<My_NGrid2_Poisson; k2++){
	    for (k1=0; k1<Ngrid1; k1++){
	      for (k3=0; k3<Ngrid3; k3++){

		weight = Kerker_weight[k2][k1][k3]; 

		re1 = Residual_ReRhok[SCFi][spin][k2][k1][k3];
		im1 = Residual_ImRhok[SCFi][spin][k2][k1][k3];

		re2 = Residual_ReRhok[SCFj][spin][k2][k1][k3];
		im2 = Residual_ImRhok[SCFj][spin][k2][k1][k3];

		My_A[SCFi][SCFj] += (re1*re2 + im1*im2)*weight;
	      }
	    }
	  }
	}
      }
    }

    /* MPI My_A */
    for (SCFi=1; SCFi<=NumMix; SCFi++){
      for (SCFj=SCFi; SCFj<=NumMix; SCFj++){
        MPI_Allreduce(&My_A[SCFi][SCFj], &A[SCFi-1][SCFj-1],
                      1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
        A[SCFj-1][SCFi-1] = A[SCFi-1][SCFj-1];
      }
    }

    /* store NormRD */

    for (i=4; 1<=i; i--){
      NormRD[i] = NormRD[i-1];
      History_Uele[i] = History_Uele[i-1];
    }
    NormRD[0] = A[0][0]/(double)Ngrid1/(double)Ngrid2/(double)Ngrid3/(double)spinmax;
    History_Uele[0] = Uele;

    /* solve the linear equation */

    for (SCFi=1; SCFi<=NumMix; SCFi++){
       A[SCFi-1][NumMix] = -1.0;
       A[NumMix][SCFi-1] = -1.0;
    }
    A[NumMix][NumMix] = 0.0;

    for (SCFi=0; SCFi<=NumMix; SCFi++){
      for (SCFj=0; SCFj<=NumMix; SCFj++){
        printf("%12.8f ",A[SCFi][SCFj]);
      }
      printf("\n");
    }

    Inverse(NumMix,A,IA);

    for (SCFi=1; SCFi<=NumMix; SCFi++){
       alden[SCFi] = -IA[SCFi-1][NumMix];
    }

    for (SCFi=1; SCFi<=NumMix; SCFi++){
      printf("myid=%2d SCFi=%2d alden=%15.12f\n",myid,SCFi,alden[SCFi]);fflush(stdout);
    }

    /*************************************
      check "nan", "NaN", "inf" or "Inf"
    *************************************/

    flag_nan = 0;
    for (SCFi=1; SCFi<=NumMix; SCFi++){

      sprintf(nanchar,"%8.4f",alden[SCFi]);
      if (strstr(nanchar,"nan")!=NULL || strstr(nanchar,"NaN")!=NULL 
       || strstr(nanchar,"inf")!=NULL || strstr(nanchar,"Inf")!=NULL){

        flag_nan = 1;
      }
    }

    if (flag_nan==1){
      for (SCFi=1; SCFi<=NumMix; SCFi++){
        alden[SCFi] = 0.0;
      }
      alden[1] = 0.1;
      alden[2] = 0.9;
    }

    /****************************************************
              calculate optimized residual rho
    ****************************************************/

    for (pSCF_iter=1; pSCF_iter<=NumMix; pSCF_iter++){

      tmp0 =  alden[pSCF_iter]; 

      for (spin=0; spin<spinmax; spin++){
	for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){

	      /* Pulay mixing */

	      Re_OptRhok[spin][k2][k1][k3] += tmp0*Residual_ReRhok[pSCF_iter][spin][k2][k1][k3];
	      Im_OptRhok[spin][k2][k1][k3] += tmp0*Residual_ImRhok[pSCF_iter][spin][k2][k1][k3];
	    }
	  }
	}
      }
    }

    /* the norm of the optimized residual rho */

    /*
    My_OptRNorm = 0.0;
    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
	for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){

	    weight = Kerker_weight[k2][k1][k3]; 

	    re1 = Re_OptRhok[spin][k2][k1][k3];
	    im1 = Im_OptRhok[spin][k2][k1][k3];

	    My_OptRNorm += (re1*re1 + im1*im1)*weight;
	  }
	}
      }
    }

    MPI_Allreduce(&My_OptRNorm, &OptRNorm, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
    OptRNorm = OptRNorm/(double)Ngrid1/(double)Ngrid2/(double)Ngrid3/(double)spinmax;
    printf("OptRNorm=%15.12f\n",sqrt(OptRNorm)); 
    */

    if      (1.0e-4<=NormRD[0])                        coef_OptRho = 0.5;
    else if (1.0e-6<=NormRD[0]  && NormRD[0]<1.0e-4)   coef_OptRho = 0.6;
    else if (1.0e-8<=NormRD[0]  && NormRD[0]<1.0e-6)   coef_OptRho = 0.7;
    else if (1.0e-10<=NormRD[0] && NormRD[0]<1.0e-8)   coef_OptRho = 0.8;
    else if (1.0e-12<=NormRD[0] && NormRD[0]<1.0e-10)  coef_OptRho = 0.9;
    else                                               coef_OptRho = 1.0;

    /****************************************************
            store the best input charge density 
    ****************************************************/

    if (NormRD[0]<BestNormRD ){

      BestNormRD = NormRD[0];

      for (spin=0; spin<spinmax; spin++){
        for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){
              ReBestRhok[spin][k2][k1][k3] = ReRhok[1][spin][k2][k1][k3];
              ImBestRhok[spin][k2][k1][k3] = ImRhok[1][spin][k2][k1][k3];
	    }
	  }
        }
      }
    }

    /****************************************************
                        mixing of rho 
    ****************************************************/

    /****************************************************
       Pulay mixing
    ****************************************************/

    if ( SCF_iter%EveryPulay_SCF==0 ) {

      /* initial ReRhok and ImRhok */

      for (spin=0; spin<spinmax; spin++){
	for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){
	      ReRhok[0][spin][k2][k1][k3] = 0.0;
	      ImRhok[0][spin][k2][k1][k3] = 0.0;
	    }
	  }
	}
      }

      /* Pulay mixing */

      for (pSCF_iter=1; pSCF_iter<=NumMix; pSCF_iter++){

	tmp0 =  alden[pSCF_iter];

	for (spin=0; spin<spinmax; spin++){
	  for (k2=0; k2<My_NGrid2_Poisson; k2++){
	    for (k1=0; k1<Ngrid1; k1++){
	      for (k3=0; k3<Ngrid3; k3++){
		ReRhok[0][spin][k2][k1][k3] += tmp0*ReRhok[pSCF_iter][spin][k2][k1][k3];
		ImRhok[0][spin][k2][k1][k3] += tmp0*ImRhok[pSCF_iter][spin][k2][k1][k3];
	      }
	    }
	  }
	}
      }

      /* Correction by optimized residual rho */

      for (spin=0; spin<spinmax; spin++){
	for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){
	      ReRhok[0][spin][k2][k1][k3] += coef_OptRho*Re_OptRhok[spin][k2][k1][k3];
	      ImRhok[0][spin][k2][k1][k3] += coef_OptRho*Im_OptRhok[spin][k2][k1][k3];
	    }
	  }
	}
      }
    }

    /****************************************************
       Kerker mixing
    ****************************************************/

    else {

      /* find an optimum mixing weight */

      Min_Weight = Min_Mixing_weight;
      Max_Weight = Max_Mixing_weight;

      if (NormRD[0]<NormRD[1]){

	/* tmp0 = 1.8*Mixing_weight; */

	tmp0 = 2.0*Mixing_weight;

	if (tmp0<Max_Weight){
	  if (Min_Weight<tmp0){
	    Mixing_weight = tmp0;
	  }
	  else{ 
	    Mixing_weight = Min_Weight;
	  }
	}
	else{ 
	  Mixing_weight = Max_Weight;
	  SCF_RENZOKU++;  
	}
      }
   
      else if (NormRD[1]<NormRD[0]){

	/* tmp0 = Mixing_weight/8.0; */

        tmp0 = Mixing_weight/8.0;

	if (tmp0<Max_Weight){
	  if (Min_Weight<tmp0)
	    Mixing_weight = tmp0;
	  else 
	    Mixing_weight = Min_Weight;
	}
	else 
	  Mixing_weight = Max_Weight;

	SCF_RENZOKU = -1;  
      }

      /* reduce Mixing_weight so that charge sloshing can be avoided */
      if ( (SCF_iter-1)%EveryPulay_SCF==0 ){
        Mixing_weight = Min_Mixing_weight;
      }

      /* Kerer mixing */

      Mix_wgt = Mixing_weight;

      for (spin=0; spin<spinmax; spin++){
	for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){

	      weight = 1.0/Kerker_weight[k2][k1][k3];
	      wgt0  = Mix_wgt*weight;
	      wgt1 =  1.0 - wgt0;

	      ReRhok[0][spin][k2][k1][k3] = wgt0*ReRhok[0][spin][k2][k1][k3]
		                          + wgt1*ReRhok[1][spin][k2][k1][k3];

	      ImRhok[0][spin][k2][k1][k3] = wgt0*ImRhok[0][spin][k2][k1][k3]
		                          + wgt1*ImRhok[1][spin][k2][k1][k3];
	    }
	  }
	}
      }
      
    }

    /****************************************************
                         shift of rho
    ****************************************************/

    for (pSCF_iter=NumSlide; 0<pSCF_iter; pSCF_iter--){
      for (spin=0; spin<spinmax; spin++){
	for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){
	      ReRhok[pSCF_iter][spin][k2][k1][k3] = ReRhok[pSCF_iter-1][spin][k2][k1][k3]; 
	      ImRhok[pSCF_iter][spin][k2][k1][k3] = ImRhok[pSCF_iter-1][spin][k2][k1][k3]; 
	    }
	  }
	}
      }
    }

    /****************************************************
                    shift of residual rho
    ****************************************************/

    for (pSCF_iter=NumSlide; 0<pSCF_iter; pSCF_iter--){
      for (spin=0; spin<spinmax; spin++){
	for (k2=0; k2<My_NGrid2_Poisson; k2++){
	  for (k1=0; k1<Ngrid1; k1++){
	    for (k3=0; k3<Ngrid3; k3++){
	      Residual_ReRhok[pSCF_iter][spin][k2][k1][k3] = Residual_ReRhok[pSCF_iter-1][spin][k2][k1][k3];
	      Residual_ImRhok[pSCF_iter][spin][k2][k1][k3] = Residual_ImRhok[pSCF_iter-1][spin][k2][k1][k3]; 
	    }
	  }
	}
      }
    }

  }

  /****************************************************
        find the charge density in real space 
  ****************************************************/

  tmp0 = 1.0/(double)Ngrid1/(double)Ngrid2/(double)Ngrid3;

  for (spin=0; spin<spinmax; spin++){

    for (k2=0; k2<My_NGrid2_Poisson; k2++){
      for (k1=0; k1<Ngrid1; k1++){
        for (k3=0; k3<Ngrid3; k3++){
          ReV2[k2][k1][k3] = ReRhok[0][spin][k2][k1][k3];
          ImV2[k2][k1][k3] = ImRhok[0][spin][k2][k1][k3];
        }
      }
    }

    if (spin==0 || spin==1){
      Get_Value_inReal(0,ReV2,ImV2,ReV1,ImV1,Density_Grid[spin],Density_Grid[spin]);

      for (MN=0; MN<My_NumGrid1; MN++){
        Density_Grid[spin][MN] = Density_Grid[spin][MN]*tmp0;
      }
    }

    else if (spin==2){
      Get_Value_inReal(1,ReV2,ImV2,ReV1,ImV1,Density_Grid[2],Density_Grid[3]);
      for (MN=0; MN<My_NumGrid1; MN++){
        Density_Grid[2][MN] = Density_Grid[2][MN]*tmp0;
        Density_Grid[3][MN] = Density_Grid[3][MN]*tmp0;
      }
    }

  }

  if (SpinP_switch==0){
    for (MN=0; MN<My_NumGrid1; MN++){
      Density_Grid[1][MN] = Density_Grid[0][MN];
    }
  }

  /****************************************************
     set ReV2 and ImV2 which are used in Poisson.c
  ****************************************************/

  if (SpinP_switch==0){
    for (k2=0; k2<My_NGrid2_Poisson; k2++){
      for (k1=0; k1<Ngrid1; k1++){
        for (k3=0; k3<Ngrid3; k3++){
          ReV2[k2][k1][k3] = 2.0*ReRhok[0][0][k2][k1][k3] - ReRhoAtomk[k2][k1][k3];
          ImV2[k2][k1][k3] = 2.0*ImRhok[0][0][k2][k1][k3] - ImRhoAtomk[k2][k1][k3];
        }
      }
    }
  }
  
  else {
    for (k2=0; k2<My_NGrid2_Poisson; k2++){
      for (k1=0; k1<Ngrid1; k1++){
        for (k3=0; k3<Ngrid3; k3++){
          ReV2[k2][k1][k3] = ReRhok[0][0][k2][k1][k3] + ReRhok[0][1][k2][k1][k3] - ReRhoAtomk[k2][k1][k3];
          ImV2[k2][k1][k3] = ImRhok[0][0][k2][k1][k3] + ImRhok[0][1][k2][k1][k3] - ImRhoAtomk[k2][k1][k3];
        }
      }
    }
  }

  /****************************************************
    freeing of arrays:

    double alden[List_YOUSO[38]];
    double A[List_YOUSO[38]][List_YOUSO[38]];
    double IA[List_YOUSO[38]][List_YOUSO[38]];
    double My_A[List_YOUSO[38]][List_YOUSO[38]];
    double Re_OptRhok[spinmax]
                            [My_NGrid2_Poisson]
                            [Ngrid1]
                            [Ngrid3];
    double Im_OptRhok[spinmax]
                            [My_NGrid2_Poisson]
                            [Ngrid1]
                            [Ngrid3];
  ****************************************************/
  
  free(alden);

  for (i=0; i<List_YOUSO[38]; i++){
    free(A[i]);
  }
  free(A);

  for (i=0; i<List_YOUSO[38]; i++){
    free(My_A[i]);
  }
  free(My_A);

  for (i=0; i<List_YOUSO[38]; i++){
    free(IA[i]);
  }
  free(IA);

  /* Re_OptRhok and Im_OptRhok */

  for (spin=0; spin<spinmax; spin++){
    for (i=0; i<My_NGrid2_Poisson; i++){
      for (j=0; j<Ngrid1; j++){
	free(Re_OptRhok[spin][i][j]);
      }
      free(Re_OptRhok[spin][i]);
    }
    free(Re_OptRhok[spin]);
  }
  free(Re_OptRhok);

  for (spin=0; spin<spinmax; spin++){
    for (i=0; i<My_NGrid2_Poisson; i++){
      for (j=0; j<Ngrid1; j++){
	free(Im_OptRhok[spin][i][j]);
      }
      free(Im_OptRhok[spin][i]);
    }
    free(Im_OptRhok[spin]);
  }
  free(Im_OptRhok);

  /* Kerker_weight */

  for (i=0; i<My_NGrid2_Poisson; i++){
    for (j=0; j<Ngrid1; j++){
      free(Kerker_weight[i][j]);
    }
    free(Kerker_weight[i]);
  }
  free(Kerker_weight);

}



void Inverse(int n, double **a, double **ia)
{
  /****************************************************
                  LU decomposition
                      0 to n
   NOTE:
   This routine does not consider the reduction of rank
  ****************************************************/

  int i,j,k;
  double w;
  double *x,*y;
  double **da;

  /***************************************************
    allocation of arrays: 

     x[List_YOUSO[38]]
     y[List_YOUSO[38]]
     da[List_YOUSO[38]][List_YOUSO[38]]
  ***************************************************/

  x = (double*)malloc(sizeof(double)*List_YOUSO[38]);
  y = (double*)malloc(sizeof(double)*List_YOUSO[38]);
  for (i=0; i<List_YOUSO[38]; i++){
    x[i] = 0.0;
    y[i] = 0.0;
  }

  da = (double**)malloc(sizeof(double*)*List_YOUSO[38]);
  for (i=0; i<List_YOUSO[38]; i++){
    da[i] = (double*)malloc(sizeof(double)*List_YOUSO[38]);
    for (j=0; j<List_YOUSO[38]; j++){
      da[i][j] = 0.0;
    }
  }

  /* start calc. */

  if (n==-1){
    for (i=0; i<List_YOUSO[38]; i++){
      for (j=0; j<List_YOUSO[38]; j++){
	a[i][j] = 0.0;
      }
    }
  }
  else{
    for (i=0; i<=n; i++){
      for (j=0; j<=n; j++){
	da[i][j] = a[i][j];
      }
    }

    /****************************************************
                     LU factorization
    ****************************************************/

    for (k=0; k<=n-1; k++){
      w = 1.0/a[k][k];
      for (i=k+1; i<=n; i++){
	a[i][k] = w*a[i][k];
	for (j=k+1; j<=n; j++){
	  a[i][j] = a[i][j] - a[i][k]*a[k][j];
	}
      }
    }

    for (k=0; k<=n; k++){

      /****************************************************
                             Ly = b
      ****************************************************/

      for (i=0; i<=n; i++){
	if (i==k)
	  y[i] = 1.0;
	else
	  y[i] = 0.0;
	for (j=0; j<=i-1; j++){
	  y[i] = y[i] - a[i][j]*y[j];
	}
      }

      /****************************************************
                             Ux = y 
      ****************************************************/

      for (i=n; 0<=i; i--){
	x[i] = y[i];
	for (j=n; (i+1)<=j; j--){
	  x[i] = x[i] - a[i][j]*x[j];
	}
	x[i] = x[i]/a[i][i];
	ia[i][k] = x[i];
      }
    }

    for (i=0; i<=n; i++){
      for (j=0; j<=n; j++){
	a[i][j] = da[i][j];
      }
    }
  }

  /***************************************************
    freeing of arrays: 

     x[List_YOUSO[38]]
     y[List_YOUSO[38]]
     da[List_YOUSO[38]][List_YOUSO[38]]
  ***************************************************/

  free(x);
  free(y);

  for (i=0; i<List_YOUSO[38]; i++){
    free(da[i]);
  }
  free(da);
}






void Complex_Inverse(int n, double **a, double **ia, double **b, double **ib)
{
  /****************************************************
                  LU decomposition
                      0 to n
   NOTE:
   This routine does not consider the reduction of rank
  ****************************************************/

  int i,j,k;
  double w2,re,im;
  dcomplex w;
  dcomplex *x,*y;
  dcomplex **da;

  /***************************************************
    allocation of arrays: 

     x[List_YOUSO[38]]
     y[List_YOUSO[38]]
     da[List_YOUSO[38]][List_YOUSO[38]]
  ***************************************************/

  x = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[38]);
  y = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[38]);

  da = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[38]);
  for (i=0; i<List_YOUSO[38]; i++){
    da[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[38]);
  }

  /* start calc. */

  if (n==-1){
    for (i=0; i<List_YOUSO[38]; i++){
      for (j=0; j<List_YOUSO[38]; j++){
 	 a[i][j] = 0.0;
	ia[i][j] = 0.0;
      }
    }
  }

  else {

    for (i=0; i<=n; i++){
      for (j=0; j<=n; j++){
	da[i][j].r =  a[i][j];
	da[i][j].i = ia[i][j];
      }
    }

    /****************************************************
                     LU factorization
    ****************************************************/

    for (k=0; k<=n-1; k++){

      w.r =  a[k][k];
      w.i = ia[k][k];
      w2 = 1.0/(w.r*w.r + w.i*w.i); 

      for (i=k+1; i<=n; i++){

	re = w2*(w.r* a[i][k] + w.i*ia[i][k]);
        im = w2*(w.r*ia[i][k] - w.i* a[i][k]);
	a[i][k] = re;
       ia[i][k] = im;

	for (j=k+1; j<=n; j++){
 	   re =  a[i][j] - (a[i][k]* a[k][j] - ia[i][k]*ia[k][j]);
	   im = ia[i][j] - (a[i][k]*ia[k][j] + ia[i][k]* a[k][j]);
 	   a[i][j] = re;
	  ia[i][j] = im;
	}
      }
    }

    for (k=0; k<=n; k++){

      /****************************************************
                             Ly = b
      ****************************************************/

      for (i=0; i<=n; i++){
	if (i==k){
	  y[i].r = 1.0;
	  y[i].i = 0.0;
	}
	else{
	  y[i].r = 0.0;
	  y[i].i = 0.0;
	}

	for (j=0; j<=i-1; j++){
   	  y[i].r = y[i].r - (a[i][j]*y[j].r - ia[i][j]*y[j].i);
	  y[i].i = y[i].i - (a[i][j]*y[j].i + ia[i][j]*y[j].r);
	}
      }

      /****************************************************
                             Ux = y 
      ****************************************************/

      for (i=n; 0<=i; i--){
	x[i] = y[i];
	for (j=n; (i+1)<=j; j--){
	  x[i].r = x[i].r - (a[i][j]*x[j].r - ia[i][j]*x[j].i);
	  x[i].i = x[i].i - (a[i][j]*x[j].i + ia[i][j]*x[j].r);
	}

        w.r =  a[i][i];
        w.i = ia[i][i];
        w2 = 1.0/(w.r*w.r + w.i*w.i); 

	re = w2*(w.r*x[i].r + w.i*x[i].i);
        im = w2*(w.r*x[i].i - w.i*x[i].r);
	x[i].r = re;
	x[i].i = im;

         b[i][k] = x[i].r;
        ib[i][k] = x[i].i;
      }
    }

    for (i=0; i<=n; i++){
      for (j=0; j<=n; j++){
	 a[i][j] = da[i][j].r;
	ia[i][j] = da[i][j].i;
      }
    }
  }

  /***************************************************
    freeing of arrays: 

     x[List_YOUSO[38]]
     y[List_YOUSO[38]]
     da[List_YOUSO[38]][List_YOUSO[38]]
  ***************************************************/

  free(x);
  free(y);

  for (i=0; i<List_YOUSO[38]; i++){
    free(da[i]);
  }
  free(da);
}
