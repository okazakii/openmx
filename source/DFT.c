/**********************************************************************
  DFT.c:

     DFT.c is a subroutine to perform self-consistent calculations
     within LDA or GGA.

  Log of DFT.c:

     22/Nov/2001  Released by T.Ozaki

***********************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "openmx_common.h"
#include "Inputtools.h"
#include "mpi.h"


#include "tran_prototypes.h"
 
/*---------- added by TOYODA 15/FEB/2010 */
#include "exx_debug.h"
#include "exx_rhox.h"
/*---------- until here */

 
static void Output_Energies_Forces(FILE *fp);
static double dUele;
void Read_SCF_keywords();


double DFT(int MD_iter, int Cnt_Now)
{
  static int firsttime=1;
  double ECE[20],Eele0[2],Eele1[2];
  double pUele,ChemP_e0[2];
  double Norm1,Norm2,Norm3,Norm4,Norm5;
  double TotalE;
  double S_coordinate[3];
  double tmp,tmp0;
  int Cnt_kind,Calc_CntOrbital_ON,spin,spinmax,m;
  int SCF_iter,SCF_iter_shift,SCF_MAX;
  int i,j,k,fft_charge_flag,M1N;
  int orbitalOpt_iter,LSCF_iter,OrbOpt_end; 
  int n,n2,wanA,ETemp_controller;
  double time0,time1,time2,time3,time4;
  double time5,time6,time7,time8,time9;
  double time10,time11,time12,time13;
  double time14,time15,etime;
  double x,y,z;
  int po3,po,TRAN_Poisson_flag2;
  int  SucceedReadingDMfile,My_SucceedReadingDMfile;
  char file_DFTSCF[YOUSO10] = ".DFTSCF";
  char file_OrbOpt[YOUSO10] = ".OrbOpt";
  char operate[200];
  char *s_vec[20];
  double TStime,TEtime;
  FILE *fp_DFTSCF;
  FILE *fp_OrbOpt;

  double ***Cluster_ReCoes; 
  double **Cluster_ko;
  
  double *ReVk;
  double *ImVk;
  double *ReRhoAtomk;
  double *ImRhoAtomk;
  double ***ReRhok;
  double ***ImRhok;
  double **Residual_ReRhok;
  double **Residual_ImRhok;

  /* for Band_DFT_Col */

  int ii,ij,ik;
  int T_knum,size_H1;
  int *MP;
  int *order_GA;
  int *My_NZeros;
  int *SP_NZeros;
  int *SP_Atoms;
  double *ko;
  double *koS;
  dcomplex **H_Band_Col;
  dcomplex **S_Band;
  dcomplex **C_Band_Col;
  dcomplex *BLAS_S;
  double *H1_Band_Col;
  double *S1_Band_Col;
  double *CDM1_Band_Col;
  double *EDM1_Band_Col;

  int ***k_op;
  int *T_k_op;
  int **T_k_ID;
  double *T_KGrids1,*T_KGrids2,*T_KGrids3;
  double ***EIGEN_Band_Col;

  int numprocs0,myid0;
  int numprocs1,myid1;

  int Num_Comm_World1;
  int Num_Comm_World2;

  int myworld1;
  int myworld2;

  int *NPROCS_ID1;
  int *Comm_World1;
  int *NPROCS_WD1;
  int *Comm_World_StartID1;
  MPI_Comm *MPI_CommWD1;

  int *NPROCS_ID2;
  int *NPROCS_WD2;
  int *Comm_World2;
  int *Comm_World_StartID2;
  MPI_Comm *MPI_CommWD2;

  /* for orbital optimization */
  double ****His_CntCoes,****His_D_CntCoes;
  double ****His_CntCoes_Species,****His_D_CntCoes_Species;
  double **OrbOpt_Hessian,*His_OrbOpt_Etot;
  int dim_H,al,p,wan;

  dtime(&TStime);

  /* MPI */ 

  MPI_Comm_size(mpi_comm_level1,&numprocs0);
  MPI_Comm_rank(mpi_comm_level1,&myid0);

  /*******************************************************
   allocation of arrays for Cluster and Band methods
  *******************************************************/

  /* band and collinear calculation */

  if (Solver==3 && SpinP_switch<=1){

    n = 0;
    for (i=1; i<=atomnum; i++){
      wanA = WhatSpecies[i];
      n += Spe_Total_CNO[wanA];
    }
    n2 = n + 2;

    MP = (int*)malloc(sizeof(int)*List_YOUSO[1]);
    order_GA = (int*)malloc(sizeof(int)*(List_YOUSO[1]+1));

    My_NZeros = (int*)malloc(sizeof(int)*numprocs0);
    SP_NZeros = (int*)malloc(sizeof(int)*numprocs0);
    SP_Atoms = (int*)malloc(sizeof(int)*numprocs0);

    ko = (double*)malloc(sizeof(double)*(n+1));
    koS = (double*)malloc(sizeof(double)*(n+1));

    H_Band_Col = (dcomplex**)malloc(sizeof(dcomplex*)*(n+1));
    for (j=0; j<n+1; j++){
      H_Band_Col[j] = (dcomplex*)malloc(sizeof(dcomplex)*(n+1));
    }

    S_Band = (dcomplex**)malloc(sizeof(dcomplex*)*(n+1));
    for (i=0; i<n+1; i++){
      S_Band[i] = (dcomplex*)malloc(sizeof(dcomplex)*(n+1));
    }

    C_Band_Col = (dcomplex**)malloc(sizeof(dcomplex*)*(n+1));
    for (j=0; j<n+1; j++){
      C_Band_Col[j] = (dcomplex*)malloc(sizeof(dcomplex)*(n+1));
    }

    BLAS_S = (dcomplex*)malloc(sizeof(dcomplex)*n*n);

    /* find size_H1 */

    size_H1 = Get_OneD_HS_Col(0, H[0], &tmp, MP, order_GA, My_NZeros, SP_NZeros, SP_Atoms);

    H1_Band_Col   = (double*)malloc(sizeof(double)*size_H1);
    S1_Band_Col   = (double*)malloc(sizeof(double)*size_H1);
    CDM1_Band_Col = (double*)malloc(sizeof(double)*size_H1);
    EDM1_Band_Col = (double*)malloc(sizeof(double)*size_H1);

    k_op = (int***)malloc(sizeof(int**)*Kspace_grid1);
    for (i=0;i<Kspace_grid1; i++) {
      k_op[i] = (int**)malloc(sizeof(int*)*Kspace_grid2);
      for (j=0;j<Kspace_grid2; j++) {
	k_op[i][j] = (int*)malloc(sizeof(int)*Kspace_grid3);
      }
    }

    for (i=0; i<Kspace_grid1; i++) {
      for (j=0; j<Kspace_grid2; j++) {
	for (k=0; k<Kspace_grid3; k++) {
	  k_op[i][j][k] = -999;
	}
      }
    }

    for (i=0; i<Kspace_grid1; i++) {
      for (j=0; j<Kspace_grid2; j++) {
	for (k=0; k<Kspace_grid3; k++) {

	  if ( k_op[i][j][k]==-999 ) {
	    k_inversion(i,j,k,Kspace_grid1,Kspace_grid2,Kspace_grid3,&ii,&ij,&ik);
	    if ( i==ii && j==ij && k==ik ) {
	      k_op[i][j][k]    = 1;
	    }

	    else {
	      k_op[i][j][k]    = 2;
	      k_op[ii][ij][ik] = 0;
	    }
	  }
	} /* k */
      } /* j */
    } /* i */
  
    /* find T_knum */

    T_knum = 0;
    for (i=0; i<Kspace_grid1; i++) {
      for (j=0; j<Kspace_grid2; j++) {
	for (k=0; k<Kspace_grid3; k++) {
	  if (0<k_op[i][j][k]){
	    T_knum++;
	  }
	}
      }
    }

    /* Monkhorst-Pack k-points */ 
    if (way_of_kpoint==2) T_knum = num_non_eq_kpt; 

    T_KGrids1 = (double*)malloc(sizeof(double)*T_knum);
    T_KGrids2 = (double*)malloc(sizeof(double)*T_knum);
    T_KGrids3 = (double*)malloc(sizeof(double)*T_knum);
    T_k_op    = (int*)malloc(sizeof(int)*T_knum);

    T_k_ID    = (int**)malloc(sizeof(int*)*2);
    for (i=0; i<2; i++){
      T_k_ID[i] = (int*)malloc(sizeof(int)*T_knum);
    }

    EIGEN_Band_Col  = (double***)malloc(sizeof(double**)*List_YOUSO[23]);
    for (i=0; i<List_YOUSO[23]; i++){
      EIGEN_Band_Col[i] = (double**)malloc(sizeof(double*)*T_knum);
      for (j=0; j<T_knum; j++){
	EIGEN_Band_Col[i][j] = (double*)malloc(sizeof(double)*(n+1));
	for (k=0; k<(n+1); k++) EIGEN_Band_Col[i][j][k] = 1.0e+5;
      }
    }

    /***********************************************
      allocation of arrays for the first world 
      and 
      make the first level worlds 
    ***********************************************/

    Num_Comm_World1 = SpinP_switch + 1; 

    NPROCS_ID1 = (int*)malloc(sizeof(int)*numprocs0); 
    Comm_World1 = (int*)malloc(sizeof(int)*numprocs0); 
    NPROCS_WD1 = (int*)malloc(sizeof(int)*Num_Comm_World1); 
    Comm_World_StartID1 = (int*)malloc(sizeof(int)*Num_Comm_World1); 
    MPI_CommWD1 = (MPI_Comm*)malloc(sizeof(MPI_Comm)*Num_Comm_World1);

    Make_Comm_Worlds(mpi_comm_level1, myid0, numprocs0, Num_Comm_World1, &myworld1, MPI_CommWD1, 
		     NPROCS_ID1, Comm_World1, NPROCS_WD1, Comm_World_StartID1);

    MPI_Comm_size(MPI_CommWD1[myworld1],&numprocs1);
    MPI_Comm_rank(MPI_CommWD1[myworld1],&myid1);

    /***********************************************
        allocation of arrays for the second world 
        and 
        make the second level worlds 
    ***********************************************/

    if (T_knum<=numprocs1){

      Num_Comm_World2 = T_knum;

      NPROCS_ID2 = (int*)malloc(sizeof(int)*numprocs1);
      Comm_World2 = (int*)malloc(sizeof(int)*numprocs1);
      NPROCS_WD2 = (int*)malloc(sizeof(int)*Num_Comm_World2);
      Comm_World_StartID2 = (int*)malloc(sizeof(int)*Num_Comm_World2);
      MPI_CommWD2 = (MPI_Comm*)malloc(sizeof(MPI_Comm)*Num_Comm_World2);

      Make_Comm_Worlds(MPI_CommWD1[myworld1], myid1, numprocs1, Num_Comm_World2, &myworld2, MPI_CommWD2, 
		       NPROCS_ID2, Comm_World2, NPROCS_WD2, Comm_World_StartID2);
    }

  } /* if (Solver==3 && SpinP_switch<=1) */

  /* band and non-collinear calculation */

  else if (Solver==3 && SpinP_switch==3){

    n = 0;
    for (i=1; i<=atomnum; i++){
      wanA = WhatSpecies[i];
      n += Spe_Total_CNO[wanA];
    }
    n2 = n + 2;

    koS = (double*)malloc(sizeof(double)*(n+1));

    S_Band = (dcomplex**)malloc(sizeof(dcomplex*)*(n+1));
    for (i=0; i<n+1; i++){
      S_Band[i] = (dcomplex*)malloc(sizeof(dcomplex)*(n+1));
    }
  }

  /* revised by Y. Xiao for Noncollinear NEGF calculations*/
  if (Solver==4 && SpinP_switch==3 ){

    n = 0;
    for (i=1; i<=atomnum; i++){
      wanA = WhatSpecies[i];
      n += Spe_Total_CNO[wanA];
    }
    n2 = n + 2;

    koS = (double*)malloc(sizeof(double)*(n+1));

    S_Band = (dcomplex**)malloc(sizeof(dcomplex*)*(n+1));
    for (i=0; i<n+1; i++){
      S_Band[i] = (dcomplex*)malloc(sizeof(dcomplex)*(n+1));
    }
  }  /* until here by Y. Xiao for Noncollinear NEGF calculations*/

  /* cluster */

  if (Solver==2){
        
    n = 0;
    for (i=1; i<=atomnum; i++){
      wanA  = WhatSpecies[i];
      n += Spe_Total_CNO[wanA];
    }
    n2 = n + 2;

    if ( SpinP_switch==0 || SpinP_switch==1 ){

      Cluster_ReCoes = (double***)malloc(sizeof(double**)*List_YOUSO[23]);
      for (i=0; i<List_YOUSO[23]; i++){
	Cluster_ReCoes[i] = (double**)malloc(sizeof(double*)*n2);
	for (j=0; j<n2; j++){
	  Cluster_ReCoes[i][j] = (double*)malloc(sizeof(double)*n2);
	}
      }
    }

    Cluster_ko = (double**)malloc(sizeof(double*)*List_YOUSO[23]);
    for (i=0; i<List_YOUSO[23]; i++){
      Cluster_ko[i] = (double*)malloc(sizeof(double)*n2);
    }

  }

  /***************************************
   allocation of arrays for data on grid
  ***************************************/

  ReVk = (double*)malloc(sizeof(double)*My_Max_NumGridB); 
  for (i=0; i<My_Max_NumGridB; i++) ReVk[i] = 0.0;

  ImVk = (double*)malloc(sizeof(double)*My_Max_NumGridB); 
  for (i=0; i<My_Max_NumGridB; i++) ImVk[i] = 0.0;

  if ( Mixing_switch==3 || Mixing_switch==4 ){

    if      (SpinP_switch==0)  spinmax = 1;
    else if (SpinP_switch==1)  spinmax = 2;
    else if (SpinP_switch==3)  spinmax = 3;

    ReRhoAtomk = (double*)malloc(sizeof(double)*My_Max_NumGridB);
    for (i=0; i<My_Max_NumGridB; i++) ReRhoAtomk[i] = 0.0;

    ImRhoAtomk = (double*)malloc(sizeof(double)*My_Max_NumGridB);
    for (i=0; i<My_Max_NumGridB; i++) ImRhoAtomk[i] = 0.0;

    ReRhok = (double***)malloc(sizeof(double**)*List_YOUSO[38]); 
    for (m=0; m<List_YOUSO[38]; m++){
      ReRhok[m] = (double**)malloc(sizeof(double*)*spinmax); 
      for (spin=0; spin<spinmax; spin++){
	ReRhok[m][spin] = (double*)malloc(sizeof(double)*My_Max_NumGridB); 
	for (i=0; i<My_Max_NumGridB; i++) ReRhok[m][spin][i] = 0.0;
      }
    }

    ImRhok = (double***)malloc(sizeof(double**)*List_YOUSO[38]); 
    for (m=0; m<List_YOUSO[38]; m++){
      ImRhok[m] = (double**)malloc(sizeof(double*)*spinmax); 
      for (spin=0; spin<spinmax; spin++){
	ImRhok[m][spin] = (double*)malloc(sizeof(double)*My_Max_NumGridB); 
	for (i=0; i<My_Max_NumGridB; i++) ImRhok[m][spin][i] = 0.0;
      }
    }

    Residual_ReRhok = (double**)malloc(sizeof(double*)*spinmax); 
    for (spin=0; spin<spinmax; spin++){
      Residual_ReRhok[spin] = (double*)malloc(sizeof(double)*My_NumGridB_CB*List_YOUSO[38]); 
      for (i=0; i<My_NumGridB_CB*List_YOUSO[38]; i++) Residual_ReRhok[spin][i] = 0.0;
    }

    Residual_ImRhok = (double**)malloc(sizeof(double*)*spinmax); 
    for (spin=0; spin<spinmax; spin++){
      Residual_ImRhok[spin] = (double*)malloc(sizeof(double)*My_NumGridB_CB*List_YOUSO[38]); 
      for (i=0; i<My_NumGridB_CB*List_YOUSO[38]; i++) Residual_ImRhok[spin][i] = 0.0;
    }
  }

  /*************************************************
      allocation of arrays for orbital optimization
  *************************************************/

  if (Cnt_switch==1){

    His_CntCoes = (double****)malloc(sizeof(double***)*(orbitalOpt_History+1));
    for (k=0; k<(orbitalOpt_History+1); k++){
      His_CntCoes[k] = (double***)malloc(sizeof(double**)*(Matomnum+1));
      for (i=0; i<=Matomnum; i++){
	His_CntCoes[k][i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
	for (j=0; j<List_YOUSO[7]; j++){
	  His_CntCoes[k][i][j] = (double*)malloc(sizeof(double)*List_YOUSO[24]);
	}
      }
    }

    His_D_CntCoes = (double****)malloc(sizeof(double***)*(orbitalOpt_History+1));
    for (k=0; k<(orbitalOpt_History+1); k++){
      His_D_CntCoes[k] = (double***)malloc(sizeof(double**)*(Matomnum+1));
      for (i=0; i<=Matomnum; i++){
	His_D_CntCoes[k][i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
	for (j=0; j<List_YOUSO[7]; j++){
	  His_D_CntCoes[k][i][j] = (double*)malloc(sizeof(double)*List_YOUSO[24]);
	}
      }
    }

    His_CntCoes_Species = (double****)malloc(sizeof(double***)*(orbitalOpt_History+1));
    for (k=0; k<(orbitalOpt_History+1); k++){
      His_CntCoes_Species[k] = (double***)malloc(sizeof(double**)*(SpeciesNum+1));
      for (i=0; i<=SpeciesNum; i++){
	His_CntCoes_Species[k][i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
	for (j=0; j<List_YOUSO[7]; j++){
	  His_CntCoes_Species[k][i][j] = (double*)malloc(sizeof(double)*List_YOUSO[24]);
	}
      }
    }

    His_D_CntCoes_Species = (double****)malloc(sizeof(double***)*(orbitalOpt_History+1));
    for (k=0; k<(orbitalOpt_History+1); k++){
      His_D_CntCoes_Species[k] = (double***)malloc(sizeof(double**)*(SpeciesNum+1));
      for (i=0; i<=SpeciesNum; i++){
	His_D_CntCoes_Species[k][i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
	for (j=0; j<List_YOUSO[7]; j++){
	  His_D_CntCoes_Species[k][i][j] = (double*)malloc(sizeof(double)*List_YOUSO[24]);
	}
      }
    }

    dim_H = 0;
    for (wan=0; wan<SpeciesNum; wan++){
      for (al=0; al<Spe_Total_CNO[wan]; al++){
	for (p=0; p<Spe_Specified_Num[wan][al]; p++){
	  dim_H++;
	}
      }
    }

    OrbOpt_Hessian = (double**)malloc(sizeof(double*)*(dim_H+2));
    for (i=0; i<(dim_H+2); i++){
      OrbOpt_Hessian[i] = (double*)malloc(sizeof(double)*(dim_H+2));
    }      

    His_OrbOpt_Etot = (double*)malloc(sizeof(double)*(orbitalOpt_History+2));
  }

  /* PrintMemory */
  if (firsttime) {

    PrintMemory("DFT: ReVk",sizeof(double)*My_Max_NumGridB,NULL);
    PrintMemory("DFT: ImVk",sizeof(double)*My_Max_NumGridB,NULL);

    if (Mixing_switch==3 || Mixing_switch==4){
      PrintMemory("DFT: ReRhoAtomk",sizeof(double)*My_Max_NumGridB,NULL);
      PrintMemory("DFT: ImRhoAtomk",sizeof(double)*My_Max_NumGridB,NULL);
      PrintMemory("DFT: ReRhok",sizeof(double)*My_Max_NumGridB*List_YOUSO[38]*spinmax,NULL);
      PrintMemory("DFT: ImRhok",sizeof(double)*My_Max_NumGridB*List_YOUSO[38]*spinmax,NULL);
      PrintMemory("DFT: Residual_ReRhok", sizeof(double)*My_NumGridB_CB*List_YOUSO[38]*spinmax,NULL);
      PrintMemory("DFT: Residual_ImRhok", sizeof(double)*My_NumGridB_CB*List_YOUSO[38]*spinmax,NULL);
    }

    if (Cnt_switch==1){

      PrintMemory("DFT: His_CntCoes",
		  sizeof(double)*(orbitalOpt_History+1)*(Matomnum+1)*List_YOUSO[7]*List_YOUSO[24],NULL);
      PrintMemory("DFT: His_D_CntCoes",
		  sizeof(double)*(orbitalOpt_History+1)*(Matomnum+1)*List_YOUSO[7]*List_YOUSO[24],NULL);
      PrintMemory("DFT: His_CntCoes_Species",
		  sizeof(double)*(orbitalOpt_History+1)*(SpeciesNum+1)*List_YOUSO[7]*List_YOUSO[24],NULL);
      PrintMemory("DFT: His_D_CntCoes_Species",
		  sizeof(double)*(orbitalOpt_History+1)*(SpeciesNum+1)*List_YOUSO[7]*List_YOUSO[24],NULL);
    }

    /* turn off firsttime flag */
    firsttime=0;
  }

  /****************************************************
    print some informations to the standard output
    and initialize times
  ****************************************************/

  if (Cnt_switch==1 && Cnt_Now==1 && MYID_MPI_COMM_WORLD==Host_ID && 0<level_stdout){
    printf("\n*******************************************************\n");      fflush(stdout);
    printf("             Orbital optimization                      \n");        fflush(stdout);
    printf("             SCF calculation at MD =%2d                \n",MD_iter);fflush(stdout);
    printf("*******************************************************\n\n");      fflush(stdout);
  } 
  else if (MYID_MPI_COMM_WORLD==Host_ID  && 0<level_stdout){
    printf("\n*******************************************************\n");      fflush(stdout);
    printf("             SCF calculation at MD =%2d                \n",MD_iter);fflush(stdout);
    printf("*******************************************************\n\n");      fflush(stdout);
  }

  fnjoint(filepath,filename,file_DFTSCF);
  fnjoint(filepath,filename,file_OrbOpt);

  /* initialize */ 

  time3  = 0.0;
  time4  = 0.0; 
  time5  = 0.0;
  time6  = 0.0;
  time7  = 0.0;
  time8  = 0.0;
  time10 = 0.0; 
  time11 = 0.0;
  time12 = 0.0;
  time13 = 0.0;
  time14 = 0.0;
  time15 = 0.0;

  for (i=0; i<20; i++) ECE[i] = 0.0;

  fft_charge_flag = 0;
  TRAN_Poisson_flag2 = 0;

  /****************************************************
         Calculations of overlap and Hamiltonian
         matrices for DFTSCF
  ****************************************************/

  if (MYID_MPI_COMM_WORLD==Host_ID  && 0<level_stdout){
    printf("<MD=%2d>  Calculation of the overlap matrix\n",MD_iter);fflush(stdout);
  }

  time1 = Set_OLP_Kin(OLP,H0);

  if (MYID_MPI_COMM_WORLD==Host_ID && 0<level_stdout){
    printf("<MD=%2d>  Calculation of the nonlocal matrix\n",MD_iter);fflush(stdout);
  }

  time2 = Set_Nonlocal(HNL,DS_NL);

  if (ProExpn_VNA==1){
    if (MYID_MPI_COMM_WORLD==Host_ID && 0<level_stdout){
      printf("<MD=%2d>  Calculation of the VNA projector matrix\n",MD_iter);fflush(stdout);
    }
    time12 = Set_ProExpn_VNA(HVNA, HVNA2, DS_VNA);  
  }

  /*---------- added by TOYODA 06/Jan/2010 */
  g_exx_switch = 0;
  if (5==XC_switch) { 
    g_exx_switch = 1;

    EXX_on_OpenMX_Init(mpi_comm_level1); 
    EXX_Initial_DM(g_exx, g_exx_DM[0]);
  }
  /*---------- until here */

  /* SCF loop */
  
  if (Cnt_switch==1 && Cnt_Now==1){
    SCF_MAX = orbitalOpt_SCF*(orbitalOpt_MD+1)-1;
    Oopt_NormD[1] = 1.0e+5;
  }
  else{
    SCF_MAX = DFTSCF_loop;
  }
  
  orbitalOpt_iter = 1;
  OrbOpt_end = 0;
  SCF_iter = 0;
  LSCF_iter = 0;
  ETemp_controller = 0;
  SCF_iter_shift = 0;
  NormRD[0] = 100.0;

  SCF_RENZOKU = -1;
  po = 0;
  pUele  = 100.0;
  Norm1 = 100.0;
  Norm2 = 100.0;
  Norm3 = 100.0;
  Norm4 = 100.0;

  /*****************************************************
              read contraction coefficients
  *****************************************************/

  if ( Cnt_switch==1 && MD_iter!=1 ) File_CntCoes("read");

  /****************************************************
            Start self consistent calculation
  ****************************************************/

  do {
    
    SCF_iter++;
    LSCF_iter++;
    
    /*****************************************************
                         print stdout
    *****************************************************/
    
    if (MYID_MPI_COMM_WORLD==Host_ID && 0<level_stdout){     
      if (Cnt_switch==1 && Cnt_Now==1){
        printf("\n***************** Orbital optimization **************\n");fflush(stdout);
        printf("    MD=%2d  orbitalOpt_iter=%2d  G-SCF=%2d  L-SCF=%2d  \n",
               MD_iter,orbitalOpt_iter,SCF_iter,LSCF_iter);fflush(stdout);
      }
      else{ 
        printf("\n******************* MD=%2d  SCF=%2d *******************\n",
               MD_iter,SCF_iter);fflush(stdout);
      }
      fflush(stdout);
    }

    /*****************************************************
     setting of densities, potentials and KS Hamiltonian
    *****************************************************/
    /****************************************************
     if (SCF==1)
    ****************************************************/

    if (SCF_iter==1){

      if (MD_iter!=1) Mixing_weight = Min_Mixing_weight;

      if (Cnt_switch==0 || (Cnt_switch==1 && Cnt_Now==1)) Cnt_kind = 0;
      else if (Cnt_switch==1)                             Cnt_kind = 1;
      else                                                Cnt_kind = 0;

      time9 = Set_Aden_Grid();

      if (Solver==4) {

        /*****************************************************
                          set grid for NEGF
	*****************************************************/

	TRAN_Set_Electrode_Grid(mpi_comm_level1, 
                                &TRAN_Poisson_flag2,                                
				Grid_Origin, tv, Left_tv, Right_tv, gtv,
				Ngrid1, Ngrid2, Ngrid3);

	/* revised by Y. Xiao for Noncollinear NEGF calculations */
	if (SpinP_switch<2) {
	  TRAN_Allocate_Lead_Region(mpi_comm_level1);
	  TRAN_Allocate_Cregion(mpi_comm_level1, SpinP_switch,atomnum,WhatSpecies,Spe_Total_CNO);
	}
	else {
	  TRAN_Allocate_Lead_Region_NC(mpi_comm_level1);
	  TRAN_Allocate_Cregion_NC(mpi_comm_level1, SpinP_switch,atomnum,WhatSpecies,Spe_Total_CNO);
	} /* until here by Y. Xiao for Noncollinear NEGF calculations*/

        /*****************************************************
                         add density from Leads
 	*****************************************************/

	TRAN_Add_Density_Lead(mpi_comm_level1,
			      SpinP_switch, Ngrid1, Ngrid2, Ngrid3,
			      My_NumGridB_AB, Density_Grid_B); 

	TRAN_Add_ADensity_Lead(mpi_comm_level1,
 	 		       SpinP_switch, Ngrid1, Ngrid2, Ngrid3,
			       My_NumGridB_AB, ADensity_Grid_B); 
      }

      time10 += Set_Orbitals_Grid(Cnt_kind);

      /******************************************************
                        read restart files
      ******************************************************/

      SucceedReadingDMfile = 0;
      if (Scf_RestartFromFile && ((Cnt_switch==1 && Cnt_Now!=1) || Cnt_switch==0 ) ) {

	My_SucceedReadingDMfile = RestartFileDFT("read",MD_iter,&Uele,H,CntH,&etime);
        time13 += etime;
        MPI_Barrier(mpi_comm_level1);
	MPI_Allreduce(&My_SucceedReadingDMfile, &SucceedReadingDMfile,
		      1, MPI_INT, MPI_PROD, mpi_comm_level1);

        if (MYID_MPI_COMM_WORLD==Host_ID && 0<level_stdout){
          if (SucceedReadingDMfile==1){
            printf("<Restart>  Found restart files\n");          fflush(stdout);
	  }
          else{
            printf("<Restart>  Could not find restart files\n"); fflush(stdout);
	  }
	}

        /***********************************************************************
         If reading the restart files is terminated, the densities on grid are 
         partially overwritten. So, Set_Aden_Grid() is called once more.  
        ***********************************************************************/

        if (SucceedReadingDMfile==0){

          time9 += Set_Aden_Grid();

          if (Solver==4) {
	    TRAN_Add_ADensity_Lead(mpi_comm_level1,
				   SpinP_switch, Ngrid1, Ngrid2, Ngrid3,
 		  	           My_NumGridB_AB, ADensity_Grid_B); 
	  }
        } 
      }

      /*****************************************************
       FFT of the initial density for k-space charge mixing 
      *****************************************************/

      if ( (Mixing_switch==3 || Mixing_switch==4) && SCF_iter==1) {

        /* non-spin polarization */
	if (SpinP_switch==0){
          if (Solver!=4 || TRAN_Poisson_flag2==2){
            time15 += FFT_Density(1,ReRhok[1][0],ImRhok[1][0]);
	  }
          else{
            time15 += FFT2D_Density(1,ReRhok[1][0],ImRhok[1][0]);
	  }
	}

	/* collinear spin polarization */
	else if (SpinP_switch==1) {
	  if (Solver!=4 || TRAN_Poisson_flag2==2){
            time15 += FFT_Density(1,ReRhok[1][0],ImRhok[1][0]);
            time15 += FFT_Density(2,ReRhok[1][1],ImRhok[1][1]);
	  }
          else{
            time15 += FFT2D_Density(1,ReRhok[1][0],ImRhok[1][0]);
            time15 += FFT2D_Density(2,ReRhok[1][1],ImRhok[1][1]);
          }
	}

	/* non-collinear spin polarization */
	else if (SpinP_switch==3) {
          if (Solver!=4 || TRAN_Poisson_flag2==2){
  	    time15 += FFT_Density(1,ReRhok[1][0],ImRhok[1][0]);
	    time15 += FFT_Density(2,ReRhok[1][1],ImRhok[1][1]);
	    time15 += FFT_Density(4,ReRhok[1][2],ImRhok[1][2]);
	  }
          else{
  	    time15 += FFT2D_Density(1,ReRhok[1][0],ImRhok[1][0]);
	    time15 += FFT2D_Density(2,ReRhok[1][1],ImRhok[1][1]);
	    time15 += FFT2D_Density(4,ReRhok[1][2],ImRhok[1][2]);
          }
	}
      }

      if (SpinP_switch==3) diagonalize_nc_density();

      /* In case the restart file is found */

      if (SucceedReadingDMfile && Cnt_switch==0) {

	if (Solver!=4 && ESM_switch==0)      time4 += Poisson(1,ReVk,ImVk);
        else if (Solver!=4 && ESM_switch!=0) time4 += Poisson_ESM(1,ReVk,ImVk); /* added by Ohwaki */
        else                                 time4 += TRAN_Poisson(ReVk,ImVk); 

	/*
	  Although the matrix of Hamiltonian is read from files, it would be better to 
	  reconstruct the matrix using charge density read from files in order to avoid 
	  abrupt change of NormRD at the second SCF step.
	  Also, if (Correct_Position_flag), then the reconstruction has to be done 
	  regardless of the above reason.
	*/

	time3 += Set_Hamiltonian("nostdout",SCF_iter,SucceedReadingDMfile,Cnt_kind,H0,HNL,DM[0],H);
      }

      /* failure of reading restart files */
      else{

	if (Solver==4) time4 += TRAN_Poisson(ReVk,ImVk); 
	time3 += Set_Hamiltonian("nostdout",SCF_iter,SucceedReadingDMfile,Cnt_kind,H0,HNL,DM[0],H);
      }  

      if (Cnt_switch==1 && Cnt_Now==1){
        if (MD_iter==1) Initial_CntCoes2(H,OLP);

	Contract_Hamiltonian(H,CntH,OLP,CntOLP);
        if (SO_switch==1) Contract_iHNL(iHNL,iCntHNL);
      } 

      /* switch a restart flag on for proceeding MD steps */
      if (MD_switch!=0) Scf_RestartFromFile = 1;

    } /* end of if (SCF_iter==1) */

    /****************************************************
     if (SCF!=1)
    ****************************************************/

    else {

      if (Cnt_switch==0 || (Cnt_switch==1 && Cnt_Now==1)) Cnt_kind = 0;
      else if (Cnt_switch==1)                             Cnt_kind = 1;
      else                                                Cnt_kind = 0;

      if (Solver!=4 && ESM_switch==0)       time4 += Poisson(fft_charge_flag,ReVk,ImVk);
      else if (Solver!=4 && ESM_switch!=0)  time4 += Poisson_ESM(fft_charge_flag,ReVk,ImVk); /* added by Ohwaki */
      else                                  time4 += TRAN_Poisson(ReVk,ImVk); 

      /* construct matrix elements for LDA+U or Zeeman term */

      if (   Hub_U_switch==1 
	     || Constraint_NCS_switch==1
	     || Zeeman_NCS_switch==1 
	     || Zeeman_NCO_switch==1) { 

	Eff_Hub_Pot(SCF_iter, OLP[0]) ;  /* added by MJ */
      }

      time3  += Set_Hamiltonian("stdout",SCF_iter,SucceedReadingDMfile,Cnt_kind,H0,HNL,DM[0],H);

      if (Cnt_switch==1 && Cnt_Now==1){
	Contract_Hamiltonian(H,CntH,OLP,CntOLP);
	if (SO_switch==1) Contract_iHNL(iHNL,iCntHNL);
      }
    }

    /****************************************************
                Solve the eigenvalue problem
    ****************************************************/

    s_vec[0]="Recursion";     s_vec[1]="Cluster"; s_vec[2]="Band";
    s_vec[3]="NEGF";          s_vec[4]="DC";      s_vec[5]="GDC";
    s_vec[6]="Cluster-DIIS";  s_vec[7]="Krylov";  s_vec[8]="Cluster2";

    if (MYID_MPI_COMM_WORLD==Host_ID && 0<level_stdout){
      printf("<%s>  Solving the eigenvalue problem...\n",s_vec[Solver-1]);fflush(stdout);
    }

    if (Cnt_switch==0){

      switch (Solver) {

      case 1:
	/* not supported */ 
	break;

      case 2:

        /*---------- modified by TOYODA 08/JAN/2010 */
        if (g_exx_switch) {
          time5 += Cluster_DFT("scf",LSCF_iter,SpinP_switch,
			       Cluster_ReCoes,Cluster_ko, H,iHNL,OLP[0],DM[0],EDM,
			       g_exx,g_exx_DM[0],g_exx_U,Eele0,Eele1);
        } else {
 
	  time5 += Cluster_DFT("scf",LSCF_iter,SpinP_switch,
			       Cluster_ReCoes,Cluster_ko, H,iHNL,OLP[0],DM[0],EDM,
			       NULL,NULL,NULL,Eele0,Eele1);
        }
        /*---------- until here */

	break;

      case 3:

	if (SpinP_switch<=1)
          /*---------- modified by TOYODA 15/FEB/2010 */
          if (5==XC_switch) {
            time5 += Band_DFT_Col(LSCF_iter,
                                  Kspace_grid1,Kspace_grid2,Kspace_grid3,
	    			  SpinP_switch,H,iHNL,OLP[0],DM[0],EDM,Eele0,Eele1, 
                                  MP,order_GA,ko,koS,EIGEN_Band_Col,
                                  H1_Band_Col,S1_Band_Col,
                                  CDM1_Band_Col,EDM1_Band_Col,
                                  H_Band_Col,S_Band,C_Band_Col, 
                                  BLAS_S,
                                  k_op,T_k_op,T_k_ID,
                                  T_KGrids1,T_KGrids2,T_KGrids3,
				  myworld1,
				  NPROCS_ID1,
				  Comm_World1,
				  NPROCS_WD1,
				  Comm_World_StartID1,
				  MPI_CommWD1,
				  myworld2,
				  NPROCS_ID2,
 			          NPROCS_WD2,
 			          Comm_World2,
				  Comm_World_StartID2,
				  MPI_CommWD2,
                                  g_exx,
                                  g_exx_DM[0],
                                  g_exx_U
                                  );
          } else {

            time5 += Band_DFT_Col(LSCF_iter,
                                  Kspace_grid1,Kspace_grid2,Kspace_grid3,
	    			  SpinP_switch,H,iHNL,OLP[0],DM[0],EDM,Eele0,Eele1, 
                                  MP,order_GA,ko,koS,EIGEN_Band_Col,
                                  H1_Band_Col,S1_Band_Col,
                                  CDM1_Band_Col,EDM1_Band_Col,
                                  H_Band_Col,S_Band,C_Band_Col, 
                                  BLAS_S,
                                  k_op,T_k_op,T_k_ID,
                                  T_KGrids1,T_KGrids2,T_KGrids3,
				  myworld1,
				  NPROCS_ID1,
				  Comm_World1,
				  NPROCS_WD1,
				  Comm_World_StartID1,
				  MPI_CommWD1,
				  myworld2,
				  NPROCS_ID2,
				  NPROCS_WD2,
				  Comm_World2,
				  Comm_World_StartID2,
				  MPI_CommWD2,
                                  NULL, NULL, NULL);
          }
	/*---------- until here */

	else 
	  time5 += Band_DFT_NonCol(LSCF_iter,
				   koS,S_Band, 
				   Kspace_grid1,Kspace_grid2,Kspace_grid3,
				   SpinP_switch,H,iHNL,OLP[0],DM[0],EDM,Eele0,Eele1);
	break;

      case 4:
	/* revised by Y. Xiao for Noncollinear NEGF calculations*/
	if (SpinP_switch<2) {
	  time5 += TRAN_DFT(mpi_comm_level1, SucceedReadingDMfile,
			    level_stdout, LSCF_iter,SpinP_switch,H,iHNL,OLP[0],
			    atomnum,Matomnum,WhatSpecies, Spe_Total_CNO, FNAN, natn,ncn,
			    M2G, G2ID, F_G2M, atv_ijk, List_YOUSO,
			    DM[0],EDM,TRAN_DecMulP,Eele0,Eele1,ChemP_e0);
	}
	else {
	  time5 += TRAN_DFT_NC(mpi_comm_level1, SucceedReadingDMfile,
			       level_stdout, LSCF_iter,SpinP_switch,H,iHNL,OLP[0],
			       atomnum,Matomnum,WhatSpecies, Spe_Total_CNO, FNAN, natn,ncn,
			       M2G, G2ID, F_G2M, atv_ijk, List_YOUSO,koS,S_Band,
			       DM[0],iDM[0],EDM,TRAN_DecMulP,Eele0,Eele1,ChemP_e0);
	} /* until here by Y. Xiao for Noncollinear NEGF calculations*/

	break;

      case 5:    
	time5 += Divide_Conquer("scf",LSCF_iter,H,iHNL,OLP[0],DM[0],EDM,Eele0,Eele1);
	break;

      case 6:    
	/* not supported */
	break;       

      case 7:    
	break;       

      case 8:
	time5 += Krylov("scf",LSCF_iter,H,iHNL,OLP[0],DM[0],EDM,Eele0,Eele1);
	break;       

      case 9:
	time5 += Cluster_DFT_ON2("scf",LSCF_iter,SpinP_switch,Cluster_ReCoes,Cluster_ko,
				 H,iHNL,OLP[0],DM[0],EDM,Eele0,Eele1);
	break;       

      }

    }
    else{

      switch (Solver) {

      case 1:
	/* not supported */ 
	break;

      case 2:

        /*---------- modified by TOYODA 08/JAN/2010 */
        if (g_exx_switch) {
          time5 += Cluster_DFT("scf",LSCF_iter,SpinP_switch,
                               Cluster_ReCoes,Cluster_ko, CntH,iCntHNL,CntOLP[0],DM[0],EDM,
                               g_exx,g_exx_DM[0],g_exx_U,Eele0,Eele1);
        } 
        else {
          time5 += Cluster_DFT("scf",LSCF_iter,SpinP_switch,
                               Cluster_ReCoes,Cluster_ko, CntH,iCntHNL,CntOLP[0],DM[0],EDM,
                               NULL,NULL,NULL,Eele0,Eele1);
        }
        /*---------- until here */

	break;

      case 3:

	if (SpinP_switch<=1)
          /*---------- modified by TOYODA 15/FEB/2010 */
          if (5==XC_switch) {
            time5 += Band_DFT_Col(LSCF_iter,
                                  Kspace_grid1,Kspace_grid2,Kspace_grid3,
	                          SpinP_switch,CntH,iCntHNL,CntOLP[0],DM[0],EDM,Eele0,Eele1, 
                                  MP,order_GA,ko,koS,EIGEN_Band_Col,
                                  H1_Band_Col,S1_Band_Col,
                                  CDM1_Band_Col,EDM1_Band_Col,
                                  H_Band_Col,S_Band,C_Band_Col, 
                                  BLAS_S,
                                  k_op,T_k_op,T_k_ID,
                                  T_KGrids1,T_KGrids2,T_KGrids3,
		  	          myworld1,
		  	          NPROCS_ID1,
			          Comm_World1,
			          NPROCS_WD1,
  			          Comm_World_StartID1,
				  MPI_CommWD1,
				  myworld2,
				  NPROCS_ID2,
				  NPROCS_WD2,
				  Comm_World2,
				  Comm_World_StartID2,
				  MPI_CommWD2,
                                  g_exx,
                                  g_exx_DM[0],
                                  g_exx_U
                                  );
          } 
          else {

            time5 += Band_DFT_Col(LSCF_iter,
                                  Kspace_grid1,Kspace_grid2,Kspace_grid3,
	                          SpinP_switch,CntH,iCntHNL,CntOLP[0],DM[0],EDM,Eele0,Eele1, 
                                  MP,order_GA,ko,koS,EIGEN_Band_Col,
                                  H1_Band_Col,S1_Band_Col,
                                  CDM1_Band_Col,EDM1_Band_Col,
                                  H_Band_Col,S_Band,C_Band_Col, 
                                  BLAS_S,
                                  k_op,T_k_op,T_k_ID,
                                  T_KGrids1,T_KGrids2,T_KGrids3,
		  	          myworld1,
		  	          NPROCS_ID1,
			          Comm_World1,
			          NPROCS_WD1,
  			          Comm_World_StartID1,
				  MPI_CommWD1,
				  myworld2,
				  NPROCS_ID2,
				  NPROCS_WD2,
				  Comm_World2,
				  Comm_World_StartID2,
				  MPI_CommWD2,
                                  NULL,
                                  NULL,
                                  NULL
                                  );
 
          }
        /*---------- until here */

	else 
	  time5 += Band_DFT_NonCol(LSCF_iter,
				   koS,S_Band, 
				   Kspace_grid1,Kspace_grid2,Kspace_grid3,
				   SpinP_switch,CntH,iCntHNL,CntOLP[0],DM[0],EDM,Eele0,Eele1);
	break;

      case 4:
	/* revised by Y. Xiao for Noncollinear NEGF calculations*/
	if (SpinP_switch<2) {
	  time5 += TRAN_DFT(mpi_comm_level1, SucceedReadingDMfile,
			    level_stdout, LSCF_iter,SpinP_switch,H,iHNL,OLP[0],
			    atomnum,Matomnum,WhatSpecies, Spe_Total_CNO, FNAN, natn,ncn,
			    M2G, G2ID, F_G2M, atv_ijk, List_YOUSO,
			    DM[0],EDM,TRAN_DecMulP,Eele0,Eele1,ChemP_e0);
	} 
	else {
	  time5 += TRAN_DFT_NC(mpi_comm_level1, SucceedReadingDMfile,
			       level_stdout, LSCF_iter,SpinP_switch,H,iHNL,OLP[0],
			       atomnum,Matomnum,WhatSpecies, Spe_Total_CNO, FNAN, natn,ncn,
			       M2G, G2ID, F_G2M, atv_ijk, List_YOUSO,koS,S_Band,
			       DM[0],iDM[0],EDM,TRAN_DecMulP,Eele0,Eele1,ChemP_e0);
	}  /* until here by Y. Xiao for Noncollinear NEGF calculations*/

	break;

      case 5:    
	time5 += Divide_Conquer("scf",LSCF_iter,CntH,iCntHNL,CntOLP[0],DM[0],EDM,Eele0,Eele1);
	break;

      case 6: 
	/* not supported */
	break;       

      case 7:
	break;       

      case 8:
	time5 += Krylov("scf",LSCF_iter,CntH,iCntHNL,CntOLP[0],DM[0],EDM,Eele0,Eele1);
	break;       

      case 9:
	time5 += Cluster_DFT_ON2("scf",LSCF_iter,SpinP_switch,Cluster_ReCoes,Cluster_ko,
				 CntH,iCntHNL,CntOLP[0],DM[0],EDM,Eele0,Eele1);
	break;       

      }
    }

    Uele_OS0 = Eele0[0];
    Uele_OS1 = Eele0[1];
    Uele_IS0 = Eele1[0];
    Uele_IS1 = Eele1[1];
    Uele = Uele_IS0 + Uele_IS1;

    /*****************************************************
                   Orbital magnetic moment 
    *****************************************************/

    if (SpinP_switch==3) Orbital_Moment("non");

    /*****************************************************
                      Mulliken charge
    *****************************************************/

    time14 += Mulliken_Charge("stdout");

    /*****************************************************
                    check SCF-convergence 
    *****************************************************/

    if (Solver==4) {
      dUele = 1.0; /* do not calculate Uele */
    }
    else {

      if (SCF_iter==1)
        dUele = 1.0;
      else  
        dUele = fabs(Uele - pUele);
    }

    if (
	2<=SCF_iter &&        
	((dUele<SCF_Criterion && NormRD[0]<(10*SCF_Criterion) && Cnt_switch==0)   ||
	 (dUele<SCF_Criterion && Cnt_switch==1 && Cnt_Now!=1)                     ||
	 (dUele<SCF_Criterion && Cnt_switch==1 && Cnt_Now==1 && OrbOpt_end==1))
	) po = 1;

    /*****************************************************
                     orbital optimization
    *****************************************************/

    if ( Cnt_switch==1 && Cnt_Now==1 && OrbOpt_end!=1 && 
	 ( LSCF_iter%orbitalOpt_SCF==0 ||
	   ( dUele<SCF_Criterion && NormRD[0]<1.0e-7 )
	   )
	 ){

      /* calculate the total energy */

      /*
	time10 += Set_Orbitals_Grid(1);
	time11 += Set_Density_Grid(1, 0, DM[0]);
	time7 += Force(H0,DS_NL,OLP,DM[0],EDM);
	time8 += Total_Energy(MD_iter,DM[0],ECE);
	time10 += Set_Orbitals_Grid(0);

	TotalE = 0.0;
	for (i=0; i<=12; i++){
	TotalE += ECE[i];
	}
	printf("orbitalOpt_iter=%3d Etot=%15.12f\n",orbitalOpt_iter,TotalE);
      */

      /* optimization of contraction coefficients */

      if (Opt_Contraction(orbitalOpt_iter,
                          TotalE,
                          H,OLP,DM[0],EDM,
                          His_CntCoes,His_D_CntCoes,
                          His_CntCoes_Species,His_D_CntCoes_Species,
                          His_OrbOpt_Etot,OrbOpt_Hessian)<orbitalOpt_criterion){

        SCF_MAX = SCF_iter + orbitalOpt_SCF - 1;
        OrbOpt_end = 1;
      }

      /* save the information at iteration */
      outputfile1(3,MD_iter,orbitalOpt_iter,0,SCF_iter,file_OrbOpt,ChemP_e0); 

      orbitalOpt_iter++;
      LSCF_iter = 0;
      SCF_iter_shift = 0;
      ETemp_controller = 0;

      if ( (orbitalOpt_SCF*(orbitalOpt_MD+1)-1) < SCF_iter) po = 1;
      if (orbitalOpt_MD < orbitalOpt_iter){
        SCF_MAX = SCF_iter + orbitalOpt_SCF - 1;
        OrbOpt_end = 1;
      }
    }

    /*****************************************************
                   mixing of density matrices
    *****************************************************/

    if (SCF_iter==1 || LSCF_iter==0)          Calc_CntOrbital_ON = 1;
    else                                      Calc_CntOrbital_ON = 0;

    /********************************************************
      control the electric temperature 
      for accelerating SCF convergence
    ********************************************************/

    if (Solver!=4 && SCF_Control_Temp==1){

      Norm5 = Norm4;
      Norm4 = Norm3;
      Norm3 = Norm2;
      Norm2 = Norm1;
      Norm1 = NormRD[0];
      tmp = (Norm1+Norm2+Norm3+Norm4+Norm5)/5.0;
      tmp0 = sqrt(fabs(NormRD[0]));

      if      (0.01<tmp0 && tmp<Norm1 && 6<LSCF_iter && ETemp_controller==0){
        E_Temp = 10.0*Original_E_Temp;
        n = LSCF_iter - Pulay_SCF + 2; 
        if (0<=n && n<LSCF_iter) SCF_iter_shift = n;
        ETemp_controller = 1;
      }
      else if (tmp0<=0.01 && ETemp_controller<=1){
        E_Temp = 1.0*Original_E_Temp;
        n = LSCF_iter - Pulay_SCF + 2; 
        if (0<=n && n<LSCF_iter) SCF_iter_shift = n;
        ETemp_controller = 2;
      }

      /* update Beta */ 
      Beta = 1.0/kB/E_Temp;
    }

    /********************************************************
     simple, RMM-DIIS, or GR-Pulay mixing for density matrix
    ********************************************************/

    if ( Mixing_switch==0
	 || Mixing_switch==1
	 || Mixing_switch==2 ){

      time6  += Mixing_DM(1,
			  LSCF_iter-SCF_iter_shift,
			  SCF_iter-SCF_iter_shift,
			  SucceedReadingDMfile,
			  ReRhok,ImRhok,
			  Residual_ReRhok,Residual_ImRhok,
			  ReVk,ImVk,ReRhoAtomk,ImRhoAtomk);

      time11 += Set_Density_Grid(Cnt_kind,Calc_CntOrbital_ON,DM[0]);

      if (Solver==4) {
	TRAN_Add_Density_Lead(mpi_comm_level1,
			      SpinP_switch, Ngrid1, Ngrid2, Ngrid3,
			      My_NumGridB_AB, Density_Grid_B); 
      }

      if (SpinP_switch==3) diagonalize_nc_density();

      fft_charge_flag = 1;
    }

    /********************************************************
         Kerker or RMM-DIISK mixing for density in k-space
    ********************************************************/
    
    else if (Mixing_switch==3 || Mixing_switch==4) {

      time11 += Set_Density_Grid(Cnt_kind,Calc_CntOrbital_ON,DM[0]);

      if (Solver==4) {
	TRAN_Add_Density_Lead(mpi_comm_level1,
			      SpinP_switch, Ngrid1, Ngrid2, Ngrid3,
			      My_NumGridB_AB, Density_Grid_B); 
      }

      /* non-spin polarization */
      if (SpinP_switch==0){
	if (Solver!=4 || TRAN_Poisson_flag2==2) time15 += FFT_Density(1,ReRhok[0][0],ImRhok[0][0]);
	else                                    time15 += FFT2D_Density(1,ReRhok[0][0],ImRhok[0][0]);
      }

      /* collinear spin polarization */
      else if (SpinP_switch==1) {
	if (Solver!=4 || TRAN_Poisson_flag2==2){
	  time15 += FFT_Density(1,ReRhok[0][0],ImRhok[0][0]);
	  time15 += FFT_Density(2,ReRhok[0][1],ImRhok[0][1]);
	}
	else{
	  time15 += FFT2D_Density(1,ReRhok[0][0],ImRhok[0][0]);
	  time15 += FFT2D_Density(2,ReRhok[0][1],ImRhok[0][1]);
	}
      }

      /* non-collinear spin polarization */
      else if (SpinP_switch==3) {
	if (Solver!=4 || TRAN_Poisson_flag2==2){
	  time15 += FFT_Density(1,ReRhok[0][0],ImRhok[0][0]);
	  time15 += FFT_Density(2,ReRhok[0][1],ImRhok[0][1]);
	  time15 += FFT_Density(4,ReRhok[0][2],ImRhok[0][2]);
	}
	else{
	  time15 += FFT2D_Density(1,ReRhok[0][0],ImRhok[0][0]);
	  time15 += FFT2D_Density(2,ReRhok[0][1],ImRhok[0][1]);
	  time15 += FFT2D_Density(4,ReRhok[0][2],ImRhok[0][2]);
	}
      }

      if (SCF_iter==1){
	if (Solver!=4 || TRAN_Poisson_flag2==2)  time15 += FFT_Density(3,ReRhoAtomk,ImRhoAtomk);
	else                                     time15 += FFT2D_Density(3,ReRhoAtomk,ImRhoAtomk);  
      }

      time6 += Mixing_DM(1,
			 LSCF_iter-SCF_iter_shift,
			 SCF_iter-SCF_iter_shift,
			 SucceedReadingDMfile,
			 ReRhok,ImRhok,
			 Residual_ReRhok,Residual_ImRhok,
			 ReVk,ImVk,ReRhoAtomk,ImRhoAtomk);

      if (SpinP_switch==3) diagonalize_nc_density();

      fft_charge_flag = 0;

    }
     
    else{
      printf("unknown scf.Mixing.Type\n");fflush(stdout);
      MPI_Finalize();
      exit(0);
    }

    /*****************************************************
        calculate occupation number for LDA+U
        and calculate the effective potential 
        for LDA+U, constraint for spin, Zeeman term for 
        spin, Zeeman term for orbital 
    *****************************************************/

    if (  Hub_U_switch==1
	  || Constraint_NCS_switch==1 
	  || Zeeman_NCS_switch==1 
	  || Zeeman_NCO_switch==1){ 

      Occupation_Number_LDA_U(SCF_iter, SucceedReadingDMfile, dUele, ECE, "stdout"); 
    }

    /*****************************************************
                       print informations
    *****************************************************/

    if (MYID_MPI_COMM_WORLD==Host_ID){

      /* spin non-collinear */
      if (SpinP_switch==3){

        if (0<level_stdout){
          printf("<DFT>  Total Spin Moment    (muB) %12.9f   Angles %14.9f  %14.9f\n",
		 2.0*Total_SpinS,Total_SpinAngle0/PI*180.0,
		 Total_SpinAngle1/PI*180.0); fflush(stdout);
          printf("<DFT>  Total Orbital Moment (muB) %12.9f   Angles %14.9f  %14.9f\n",
		 Total_OrbitalMoment,
		 Total_OrbitalMomentAngle0/PI*180.0,
		 Total_OrbitalMomentAngle1/PI*180.0); fflush(stdout);
	}

        x = 2.0*Total_SpinS*sin(Total_SpinAngle0)*cos(Total_SpinAngle1);
        y = 2.0*Total_SpinS*sin(Total_SpinAngle0)*sin(Total_SpinAngle1);
        z = 2.0*Total_SpinS*cos(Total_SpinAngle0);

        x += Total_OrbitalMoment*sin(Total_OrbitalMomentAngle0)*cos(Total_OrbitalMomentAngle1);
        y += Total_OrbitalMoment*sin(Total_OrbitalMomentAngle0)*sin(Total_OrbitalMomentAngle1);
        z += Total_OrbitalMoment*cos(Total_OrbitalMomentAngle0);

        xyz2spherical( x,y,z, 0.0,0.0,0.0, S_coordinate ); 

        if (0<level_stdout){
          printf("<DFT>  Total Moment         (muB) %12.9f   Angles %14.9f  %14.9f\n",
                 S_coordinate[0],S_coordinate[1]/PI*180.0,S_coordinate[2]/PI*180.0);
	  fflush(stdout);
	}

      }
      /* spin collinear */
      else{
        if (0<level_stdout){
          printf("<DFT>  Total Spin Moment (muB) = %15.12f\n",2.0*Total_SpinS);fflush(stdout);
	}
      }

      if (Mixing_switch==0 && 0<level_stdout){
        printf("<DFT>  Mixing_weight=%15.12f SCF_RENZOKU=%2d\n",
	       Mixing_weight,SCF_RENZOKU);fflush(stdout);
      }
      else{
        if (0<level_stdout){
          printf("<DFT>  Mixing_weight=%15.12f\n",Mixing_weight);fflush(stdout);
	}
      }

      if (0<level_stdout){
        printf("<DFT>  Uele   =%18.12f  dUele     =%17.12f\n",Uele,dUele);fflush(stdout);
        printf("<DFT>  NormRD =%18.12f  Criterion =%17.12f\n",
               sqrt(fabs(NormRD[0])),SCF_Criterion);fflush(stdout);
      }
    }

    if (Solver==4) {

      /* check whether the SCF converges or not */

      if ( sqrt(fabs(NormRD[0]))< SCF_Criterion ) {

	po = 1;

	/* printf("TRAN  NormRD < SCF_Criterion, SCF is satisfied \n");	*/
      } 
    }

    outputfile1(1,MD_iter,orbitalOpt_iter,Cnt_Now,SCF_iter,file_DFTSCF,ChemP_e0); 

    /*****************************************************
                         Uele -> pUele
    *****************************************************/

    pUele = Uele;

    /*************************************************************
     If there is a proper file, change parameters controlling the 
     SCF iteration, based on keywords described in the file.
    **************************************************************/

    Read_SCF_keywords();
    if (Cnt_switch==0) SCF_MAX = DFTSCF_loop;
  
    /************************************************************************
                              end of SCF calculation
    ************************************************************************/

  } while (po==0 && SCF_iter<SCF_MAX);

  /*****************************************************
          making of the input data for TranMain
  *****************************************************/

  /* revised by Y. Xiao for Noncollinear NEGF calculations */ /* iHNL is outputed */
  TRAN_Output_Trans_HS( mpi_comm_level1, Solver, SpinP_switch, ChemP, H, iHNL, OLP,
                        atomnum, SpeciesNum, WhatSpecies, 
                        Spe_Total_CNO, FNAN, natn, ncn, G2ID, atv_ijk,
                        Max_FSNAN, ScaleSize, F_G2M, TCpyCell, List_YOUSO, 
                        filepath,filename,"tranb");

  /*****************************************************
              save contraction coefficients
  *****************************************************/

  if (Cnt_switch==1) File_CntCoes("write");

  /*****************************************************
              save orbital magnetic moment 
  *****************************************************/

  if (SpinP_switch==3) Orbital_Moment("write");

  /*****************************************************
              save Mulliken charge
  *****************************************************/

  time14 += Mulliken_Charge("write");

  /*****************************************************
            save occupation number in LDA+U
  *****************************************************/

  /* ---- added by MJ */
  if (  Hub_U_switch==1
	|| Constraint_NCS_switch==1 
	|| Zeeman_NCS_switch==1 
	|| Zeeman_NCO_switch==1){ 

    Occupation_Number_LDA_U(SCF_iter, SucceedReadingDMfile, dUele, ECE, "write"); 

  }

  /*****************************************************
                      Band dispersion 
  *****************************************************/

  if ( Band_disp_switch==1 && Band_Nkpath>0 && (Solver==3 || PeriodicGamma_flag==1) ){
    if (Cnt_switch==0){
      Band_DFT_kpath(Band_Nkpath, Band_N_perpath,
                     Band_kpath, Band_kname, 
                     SpinP_switch,H,iHNL,OLP[0]);
    }
    else {
      Band_DFT_kpath(Band_Nkpath, Band_N_perpath,
                     Band_kpath, Band_kname,
                     SpinP_switch,CntH,iCntHNL,CntOLP[0]);
    }
  }            

  /******************************************************
          calculation of density of states (DOS)
  ******************************************************/

  if ( Dos_fileout || DosGauss_fileout) {

    if ( Solver==2 && PeriodicGamma_flag==0 ){  /* cluster */

      if ( Cnt_switch==0 ){

        if (Opticalconductivity_fileout==1){
  	  Cluster_DFT_Dosout( SpinP_switch, H, iHNL, OLP[0] );
	}
        else {      

          /*---------- modified by TOYODA 08/JAN/2010 */
          if (g_exx_switch) {
            time5 += Cluster_DFT("dos",LSCF_iter,SpinP_switch,
				 Cluster_ReCoes,Cluster_ko, H,iHNL,OLP[0],DM[0],EDM,
				 g_exx,g_exx_DM[0],g_exx_U,Eele0,Eele1);
          } else {
            time5 += Cluster_DFT("dos",LSCF_iter,SpinP_switch,
				 Cluster_ReCoes,Cluster_ko, H,iHNL,OLP[0],DM[0],EDM,
				 NULL,NULL,NULL,Eele0,Eele1);
          }
          /*---------- until here */
	}
      }
      else {

        if (Opticalconductivity_fileout==1){
 	  Cluster_DFT_Dosout( SpinP_switch, CntH, iCntHNL, CntOLP[0] );
	}
        else {

          /*---------- modifined by TOYODA 08/JAN/2010 */
          if (g_exx_switch) {
            
            time5 += Cluster_DFT("dos",LSCF_iter,SpinP_switch,
				 Cluster_ReCoes,Cluster_ko,CntH,iCntHNL,CntOLP[0],DM[0],EDM,
				 g_exx,g_exx_DM[0],g_exx_U,Eele0,Eele1);
          } else {
            time5 += Cluster_DFT("dos",LSCF_iter,SpinP_switch,
				 Cluster_ReCoes,Cluster_ko,CntH,iCntHNL,CntOLP[0],DM[0],EDM,
				 NULL,NULL,NULL,Eele0,Eele1);
          }
          /*---------- until here */
        }
      }
    }

    if ( Solver==3 || PeriodicGamma_flag==1 ){  /* band */
      if (Cnt_switch==0){

	Band_DFT_Dosout( Dos_Kgrid[0], Dos_Kgrid[1], Dos_Kgrid[2],
			 SpinP_switch, H, iHNL, OLP[0] );
      }
      else {
	Band_DFT_Dosout( Dos_Kgrid[0], Dos_Kgrid[1], Dos_Kgrid[2],
			 SpinP_switch, CntH, iCntHNL, CntOLP[0] );
      }
    }

    if (Solver==4){ /* NEGF */

      TRAN_DFT_Dosout( mpi_comm_level1, level_stdout, LSCF_iter,SpinP_switch,H,iHNL,OLP[0],
		       atomnum,Matomnum,WhatSpecies, Spe_Total_CNO, FNAN, natn, ncn, 
		       M2G, G2ID, atv_ijk, List_YOUSO, Spe_Num_CBasis, SpeciesNum, filename, filepath,
		       DM[0],EDM,Eele0,Eele1 );

    }

    if ( Solver==5 ){  /* divide-conquer */
      if (Cnt_switch==0){
        time5 += Divide_Conquer("dos",LSCF_iter,H,iHNL,OLP[0],DM[0],EDM,Eele0,Eele1);
      }
      else {
        time5 += Divide_Conquer("dos",LSCF_iter,CntH,iCntHNL,CntOLP[0],DM[0],EDM,Eele0,Eele1);
      }
    }

    if ( Solver==8 ){  /* Krylov subspace method */
      if (Cnt_switch==0){
        time5 += Krylov("dos",LSCF_iter,H,iHNL,OLP[0],DM[0],EDM,Eele0,Eele1);
      }
      else {
        time5 += Krylov("dos",LSCF_iter,CntH,iCntHNL,CntOLP[0],DM[0],EDM,Eele0,Eele1);
      }
    }
  }

  /*****************************************************
         Calc. of Bloch waves at given k-points
  *****************************************************/

  if (MO_fileout==1 && MO_Nkpoint>0 && (Solver==3 || PeriodicGamma_flag==1) ) {
    if (Cnt_switch==0){
      Band_DFT_MO(MO_Nkpoint, MO_kpoint, SpinP_switch, H, iHNL, OLP[0]);
    }
    else {
      Band_DFT_MO(MO_Nkpoint, MO_kpoint, SpinP_switch, CntH, iCntHNL, CntOLP[0]);
    }
  }

  /*****************************************************
       write a restart file if SpinP_switch!=3
  *****************************************************/

  if (SpinP_switch!=3){
    RestartFileDFT("write",MD_iter,&Uele,H,CntH,&etime);
    time13 += etime;
    MPI_Barrier(mpi_comm_level1);
  }

  /****************************************************
                         force
      If you first call Total_Energy.c before
      Force.c, then the force calculation will fail. 
  ****************************************************/

  if (MYID_MPI_COMM_WORLD==Host_ID && 0<level_stdout){
    printf("<MD=%2d>  Force calculation\n",MD_iter);fflush(stdout);
  }

  if (Cnt_switch==1 && Cnt_Now==1) time10 += Set_Orbitals_Grid(1);
  time11 += Set_Density_Grid(1, 0, DM[0]);

  if (Solver==4) {
    TRAN_Add_Density_Lead(mpi_comm_level1,
			  SpinP_switch, Ngrid1, Ngrid2, Ngrid3,
 	                  My_NumGridB_AB, Density_Grid_B); 
  }

  /* to save mixing densities between the previous and current ones
     do charge mixing */

  if (Mixing_switch==3 || Mixing_switch==4) {

    /* non-spin polarization */
    if (SpinP_switch==0){
      if (Solver!=4 || TRAN_Poisson_flag2==2)  time15 += FFT_Density(1,ReRhok[0][0],ImRhok[0][0]);
      else                                     time15 += FFT2D_Density(1,ReRhok[0][0],ImRhok[0][0]);
    }

    /* collinear spin polarization */
    else if (SpinP_switch==1) {
      if (Solver!=4 || TRAN_Poisson_flag2==2){
        time15 += FFT_Density(1,ReRhok[0][0],ImRhok[0][0]);
        time15 += FFT_Density(2,ReRhok[0][1],ImRhok[0][1]);
      }
      else{
        time15 += FFT2D_Density(1,ReRhok[0][0],ImRhok[0][0]);
        time15 += FFT2D_Density(2,ReRhok[0][1],ImRhok[0][1]);
      }
    }
    /* non-collinear spin polarization */
    else if (SpinP_switch==3) {
      if (Solver!=4 || TRAN_Poisson_flag2==2){
        time15 += FFT_Density(1,ReRhok[0][0],ImRhok[0][0]);
        time15 += FFT_Density(2,ReRhok[0][1],ImRhok[0][1]);
        time15 += FFT_Density(4,ReRhok[0][2],ImRhok[0][2]);
      }
      else{
        time15 += FFT2D_Density(1,ReRhok[0][0],ImRhok[0][0]);
        time15 += FFT2D_Density(2,ReRhok[0][1],ImRhok[0][1]);
        time15 += FFT2D_Density(4,ReRhok[0][2],ImRhok[0][2]);
      }
    }

    time6 += Mixing_DM(1,
                       LSCF_iter-SCF_iter_shift,
                       SCF_iter-SCF_iter_shift,
                       SucceedReadingDMfile,
		       ReRhok,ImRhok,
		       Residual_ReRhok,Residual_ImRhok,
		       ReVk,ImVk,ReRhoAtomk,ImRhoAtomk);
  }

  /* write a restart file, it is important to save
     charge density before calling diagonalize_nc_density */

  if (SpinP_switch==3){
    RestartFileDFT("write",MD_iter,&Uele,H,CntH,&etime);
    time13 += etime;
    MPI_Barrier(mpi_comm_level1);
  }

  if (SpinP_switch==3) diagonalize_nc_density();

  /****************************************************
   if (Solver==9), calculate the energy density matrix 
  ****************************************************/

  if (Solver==9){

    if (Cnt_switch==0){
      time5 += Cluster_DFT_ON2("force",LSCF_iter,SpinP_switch,Cluster_ReCoes,Cluster_ko,
			       H,iHNL,OLP[0],DM[0],EDM,Eele0,Eele1);
    }
    else{

      time5 += Cluster_DFT_ON2("force",LSCF_iter,SpinP_switch,Cluster_ReCoes,Cluster_ko,
			       CntH,iCntHNL,CntOLP[0],DM[0],EDM,Eele0,Eele1);
    }
  }

  /****************************************************
               calculation of forces
  ****************************************************/

  if (!orbitalOpt_Force_Skip) time7 += Force(H0,DS_NL,OLP,DM[0],EDM);

  /****************************************************
     if the SCF iteration converged, change Pulay_SCF
  ****************************************************/

  if (po==1){
    Pulay_SCF = 3;
  }

  /****************************************************
               calculate the total energy 
  ****************************************************/

  if (MYID_MPI_COMM_WORLD==Host_ID && 0<level_stdout){
    printf("<MD=%2d>  Total Energy\n",MD_iter);fflush(stdout);
  }

  if (!orbitalOpt_Force_Skip) time8 += Total_Energy(MD_iter,DM[0],ECE);

  Uatom  = 0.0;
  Ucore  = ECE[0];
  UH0    = ECE[1];
  Ukin   = ECE[2];
  Una    = ECE[3];
  Unl    = ECE[4];
  UH1    = ECE[5];
  Uxc0   = ECE[6];
  Uxc1   = ECE[7];
  Uhub   = ECE[8];    
  Ucs    = ECE[9]; 
  Uzs    = ECE[10];
  Uzo    = ECE[11];
  Uef    = ECE[12];
  UvdW   = ECE[13];

  Utot  = Ucore + UH0 + Ukin + Una + Unl + UH1
    + Uxc0 + Uxc1 + Uhub + Ucs + Uzs + Uzo + Uef + UvdW; 

  /*---------- added by TOYODA 20/JAN/2010 */
  if (g_exx_switch) { Utot += g_exx_U[0] + g_exx_U[1]; }
  /*---------- until here */

  /* elapsed time */
  dtime(&TEtime);
  time0 = TEtime - TStime;

  if (MYID_MPI_COMM_WORLD==Host_ID && 0<level_stdout){

    printf("\n*******************************************************\n"); 
    printf("                Total Energy (Hartree) at MD =%2d        \n",MD_iter);
    printf("*******************************************************\n\n"); 

    printf("  Uele  = %20.12f\n\n",Uele);
    printf("  Ukin  = %20.12f\n",Ukin);
    printf("  UH0   = %20.12f\n",UH0);
    printf("  UH1   = %20.12f\n",UH1);
    printf("  Una   = %20.12f\n",Una);
    printf("  Unl   = %20.12f\n",Unl);
    printf("  Uxc0  = %20.12f\n",Uxc0);
    printf("  Uxc1  = %20.12f\n",Uxc1);
    /*---------- added by TOYODA 20/JAN/2010 */
    if (g_exx_switch) {
      printf("  Uexx0 = %20.12f\n",g_exx_U[0]);
      printf("  Uexx1 = %20.12f\n",g_exx_U[1]);
    }
    printf("  Ucore = %20.12f\n",Ucore);
    printf("  Uhub  = %20.12f\n",Uhub);
    printf("  Ucs   = %20.12f\n",Ucs);
    printf("  Uzs   = %20.12f\n",Uzs);
    printf("  Uzo   = %20.12f\n",Uzo);
    printf("  Uef   = %20.12f\n",Uef);
    printf("  UvdW  = %20.12f\n",UvdW);
    printf("  Utot  = %20.12f\n\n",Utot);

    printf("  Note:\n\n");
    /*---------- modified by TOYODA 20/JAN/2010 */
    if (g_exx_switch) {
      printf("  Utot = Ukin+UH0+UH1+Una+Unl+Uxc0+Uxc1+Uexx0+Uexx1");
      printf("+Ucore+Uhub+Ucs+Uzs+Uzo+Uef+UvdW\n\n");
    } else {
      printf("  Utot = Ukin+UH0+UH1+Una+Unl+Uxc0+Uxc1");
      printf("+Ucore+Uhub+Ucs+Uzs+Uzo+Uef+UvdW\n\n");
    }
    /*---------- until here */
    printf("  Uele:   band energy\n");
    printf("  Ukin:   kinetic energy\n");
    printf("  UH0:    electric part of screened Coulomb energy\n");
    printf("  UH1:    difference electron-electron Coulomb energy\n");
    printf("  Una:    neutral atom potential energy\n");
    printf("  Unl:    non-local potential energy\n");
    printf("  Uxc0:   exchange-correlation energy for alpha spin\n");
    printf("  Uxc1:   exchange-correlation energy for beta spin\n");
    printf("  Ucore:  core-core Coulomb energy\n");
    printf("  Uhub:   LDA+U energy\n");
    printf("  Ucs:    constraint energy for spin orientation\n");
    printf("  Uzs:    Zeeman term for spin magnetic moment\n");
    printf("  Uzo:    Zeeman term for orbital magnetic moment\n");
    printf("  Uef:    electric energy by electric field\n");
    printf("  UvdW:   semi-empirical vdW energy \n\n"); /* okuno */
    printf("  (see also PRB 72, 045121(2005) for the energy contributions)\n\n");

    printf("\n*******************************************************\n"); 
    printf("           Computational times (s) at MD =%2d            \n",MD_iter);
    printf("*******************************************************\n\n"); 

    printf("  DFT in total      = %10.5f\n\n",time0);
    printf("  Set_OLP_Kin       = %10.5f\n",time1);
    printf("  Set_Nonlocal      = %10.5f\n",time2);
    printf("  Set_ProExpn_VNA   = %10.5f\n",time12);
    printf("  Set_Hamiltonian   = %10.5f\n",time3);
    printf("  Poisson           = %10.5f\n",time4);
    printf("  diagonalization   = %10.5f\n",time5);
    printf("  Mixing_DM         = %10.5f\n",time6);
    printf("  Force             = %10.5f\n",time7);
    printf("  Total_Energy      = %10.5f\n",time8);
    printf("  Set_Aden_Grid     = %10.5f\n",time9);
    printf("  Set_Orbitals_Grid = %10.5f\n",time10);
    printf("  Set_Density_Grid  = %10.5f\n",time11);
    printf("  RestartFileDFT    = %10.5f\n",time13);
    printf("  Mulliken_Charge   = %10.5f\n",time14);
    printf("  FFT(2D)_Density   = %10.5f\n",time15);

    /*---------- added by TOYODA 18/FEB/2010 */
    if (g_exx_switch) { EXX_Time_Report(); }
    /*---------- until here */
  }

  outputfile1(2,MD_iter,0,0,SCF_iter,file_DFTSCF,ChemP_e0);

  /****************************************************
               Output energies and Forces 
  ****************************************************/

  /* The Sum of Atomic Energy */
  /* if (iter==1) Atomic_Energy();  */
  /* Output_Energies_Forces(fp_DFTSCF); */

  CompTime[myid0][5]  += time1;  /* Set_OLP_Kin       */
  CompTime[myid0][6]  += time2;  /* Set_Nonlocal      */
  CompTime[myid0][7]  += time3;  /* Set_Hamiltonian   */
  CompTime[myid0][8]  += time4;  /* Poisson           */
  CompTime[myid0][9]  += time5;  /* diagonalization   */
  CompTime[myid0][10] += time6;  /* Mixing_DM         */
  CompTime[myid0][11] += time7;  /* Force             */
  CompTime[myid0][12] += time8;  /* Total_Energy      */
  CompTime[myid0][13] += time9;  /* Set_Aden_Grid     */
  CompTime[myid0][14] += time10; /* Set_Orbitals_Grid */
  CompTime[myid0][15] += time11; /* Set_Density_Grid  */
  CompTime[myid0][16] += time12; /* Set_ProExpn_VNA   */
  CompTime[myid0][17] += time13; /* RestartFileDFT    */
  CompTime[myid0][18] += time14; /* Mulliken_Charge   */
  CompTime[myid0][19] += time15; /* FFT(2D)_Density   */

  /*********************************************************
   freeing of arrays 
  *********************************************************/

  /*---------- added by TOYODA 06/Jan/2010 */
  if (g_exx_switch) { 
    if (g_exx_switch_output_DM && Host_ID==myid1) { 
      EXX_Output_DM(g_exx, g_exx_DM[0]); 
    }
    EXX_on_OpenMX_Free(); 
  }
  /*---------- until here */

  free(ReVk);
  free(ImVk);

  if ( Mixing_switch==3 || Mixing_switch==4 ){

    if      (SpinP_switch==0)  spinmax = 1;
    else if (SpinP_switch==1)  spinmax = 2;
    else if (SpinP_switch==3)  spinmax = 3;

    free(ReRhoAtomk);
    free(ImRhoAtomk);

    for (m=0; m<List_YOUSO[38]; m++){
      for (spin=0; spin<spinmax; spin++){
	free(ReRhok[m][spin]);
      }
      free(ReRhok[m]);
    }
    free(ReRhok);

    for (m=0; m<List_YOUSO[38]; m++){
      for (spin=0; spin<spinmax; spin++){
	free(ImRhok[m][spin]);
      }
      free(ImRhok[m]);
    }
    free(ImRhok);

    for (spin=0; spin<spinmax; spin++){
      free(Residual_ReRhok[spin]);
    }
    free(Residual_ReRhok);

    for (spin=0; spin<spinmax; spin++){
      free(Residual_ImRhok[spin]);
    }
    free(Residual_ImRhok);
  }

  /*************************************************
      allocation of arrays for orbital optimization
  *************************************************/

  if (Cnt_switch==1){

    for (k=0; k<(orbitalOpt_History+1); k++){
      for (i=0; i<=Matomnum; i++){
	for (j=0; j<List_YOUSO[7]; j++){
	  free(His_CntCoes[k][i][j]);
	}
	free(His_CntCoes[k][i]);
      }
      free(His_CntCoes[k]);
    }
    free(His_CntCoes);

    for (k=0; k<(orbitalOpt_History+1); k++){
      for (i=0; i<=Matomnum; i++){
	for (j=0; j<List_YOUSO[7]; j++){
	  free(His_D_CntCoes[k][i][j]);
	}
	free(His_D_CntCoes[k][i]);
      }
      free(His_D_CntCoes[k]);
    }
    free(His_D_CntCoes);

    for (k=0; k<(orbitalOpt_History+1); k++){
      for (i=0; i<=SpeciesNum; i++){
	for (j=0; j<List_YOUSO[7]; j++){
	  free(His_CntCoes_Species[k][i][j]);
	}
	free(His_CntCoes_Species[k][i]);
      }
      free(His_CntCoes_Species[k]);
    }
    free(His_CntCoes_Species);

    for (k=0; k<(orbitalOpt_History+1); k++){
      for (i=0; i<=SpeciesNum; i++){
	for (j=0; j<List_YOUSO[7]; j++){
	  free(His_D_CntCoes_Species[k][i][j]);
	}
	free(His_D_CntCoes_Species[k][i]);
      }
      free(His_D_CntCoes_Species[k]);
    }
    free(His_D_CntCoes_Species);

    dim_H = 0;
    for (wan=0; wan<SpeciesNum; wan++){
      for (al=0; al<Spe_Total_CNO[wan]; al++){
	for (p=0; p<Spe_Specified_Num[wan][al]; p++){
	  dim_H++;
	}
      }
    }

    for (i=0; i<(dim_H+2); i++){
      free(OrbOpt_Hessian[i]);
    }      
    free(OrbOpt_Hessian);

    free(His_OrbOpt_Etot);
  }

  /* band and collinear calculation */

  if (Solver==3 && SpinP_switch<=1){

    n = 0;
    for (i=1; i<=atomnum; i++){
      wanA  = WhatSpecies[i];
      n += Spe_Total_CNO[wanA];
    }
    n2 = n + 2;

    free(MP);
    free(order_GA);

    free(My_NZeros);
    free(SP_NZeros);
    free(SP_Atoms);

    free(ko);
    free(koS);

    for (j=0; j<n+1; j++){
      free(H_Band_Col[j]);
    }
    free(H_Band_Col);

    for (i=0; i<n+1; i++){
      free(S_Band[i]);
    }
    free(S_Band);

    for (j=0; j<n+1; j++){
      free(C_Band_Col[j]);
    }
    free(C_Band_Col);

    free(BLAS_S);

    free(H1_Band_Col);
    free(S1_Band_Col);
    free(CDM1_Band_Col);
    free(EDM1_Band_Col);

    for (i=0; i<Kspace_grid1; i++) {
      for (j=0; j<Kspace_grid2; j++) {
	free(k_op[i][j]);
      }
      free(k_op[i]);
    }
    free(k_op);

    free(T_KGrids1);
    free(T_KGrids2);
    free(T_KGrids3);
    free(T_k_op);

    for (i=0; i<2; i++){
      free(T_k_ID[i]);
    }
    free(T_k_ID);

    for (i=0; i<List_YOUSO[23]; i++){
      for (j=0; j<T_knum; j++){
	free(EIGEN_Band_Col[i][j]);
      }
      free(EIGEN_Band_Col[i]);
    }
    free(EIGEN_Band_Col);

    /* freeing of arrays for the second world */

    if (T_knum<=numprocs1){

      if (Num_Comm_World2<=numprocs1){
        MPI_Comm_free(&MPI_CommWD2[myworld2]);
      }

      free(MPI_CommWD2);
      free(Comm_World_StartID2);
      free(NPROCS_WD2);
      free(Comm_World2);
      free(NPROCS_ID2);
    }

    /* freeing of arrays for the first world */

    if (Num_Comm_World1<=numprocs0){
      MPI_Comm_free(&MPI_CommWD1[myworld1]);
    }

    free(MPI_CommWD1);
    free(Comm_World_StartID1);
    free(NPROCS_WD1);
    free(Comm_World1);
    free(NPROCS_ID1);
  }

  /* band and non-collinear calculation */

  else if (Solver==3 && SpinP_switch==3){

    n = 0;
    for (i=1; i<=atomnum; i++){
      wanA = WhatSpecies[i];
      n += Spe_Total_CNO[wanA];
    }
    n2 = n + 2;

    free(koS);

    for (i=0; i<n+1; i++){
      free(S_Band[i]);
    }
    free(S_Band);
  }

  /* cluster */
 
  if (Solver==2){
        
    n = 0;
    for (i=1; i<=atomnum; i++){
      wanA  = WhatSpecies[i];
      n += Spe_Total_CNO[wanA];
    }
    n2 = n + 2;

    if ( SpinP_switch==0 || SpinP_switch==1 ){
      for (i=0; i<List_YOUSO[23]; i++){
	for (j=0; j<n2; j++){
	  free(Cluster_ReCoes[i][j]);
	}
        free(Cluster_ReCoes[i]);
      }
      free(Cluster_ReCoes);
    }

    for (i=0; i<List_YOUSO[23]; i++){
      free(Cluster_ko[i]);
    }
    free(Cluster_ko);
  }

  if (Solver==4)  {
    /* revised by Y. Xiao for Noncollinear NEGF calculations */
    if (SpinP_switch<2) {
      TRAN_Deallocate_Lead_Region();
      TRAN_Deallocate_Cregion( SpinP_switch );
    } else {
      TRAN_Deallocate_Lead_Region_NC();
      TRAN_Deallocate_Cregion_NC( SpinP_switch );
    }
    TRAN_Deallocate_Electrode_Grid(Ngrid2);
  } /* if (Solver==4) */

  /* freeing koS and S_Band */
  if (Solver==4 && SpinP_switch==3){
    n = 0;
    for (i=1; i<=atomnum; i++){
      wanA = WhatSpecies[i];
      n += Spe_Total_CNO[wanA];
    }
    n2 = n + 2;

    free(koS);

    for (i=0; i<n+1; i++){
      free(S_Band[i]);
    }
    free(S_Band);
  }
  /* until here by Y. Xiao for Noncollinear NEGF calculations */

  MPI_Barrier(MPI_COMM_WORLD1);

  return time0;
}




void Read_SCF_keywords()
{
  int po,myid;
  char fname[YOUSO10];
  FILE *fp;

  sprintf(fname,"%s%s_SCF_keywords",filepath,filename);

  if ((fp = fopen(fname,"r")) != NULL){
    po = 1;
  }
  else {
    po = 0;
  }
  
  if (po==1){

    MPI_Comm_rank(MPI_COMM_WORLD1,&myid);

    if (myid==Host_ID){
      printf("\n The keywords for SCF iteration are renewed by %s.\n",fname);fflush(stdout);
    }
    
    /* open the file */

    input_open(fname);

    /* scf.maxIter */

    input_int("scf.maxIter",&DFTSCF_loop,40);
    if (DFTSCF_loop<0) {
      MPI_Finalize();
      exit(0);
    }

    /* scf.Min.Mixing.Weight */

    input_double("scf.Min.Mixing.Weight",&Min_Mixing_weight,(double)0.001);

    /* scf.Max.Mixing.Weight */

    input_double("scf.Max.Mixing.Weight",&Max_Mixing_weight,(double)0.4);

    /* scf.Kerker.factor */

    input_double("scf.Kerker.factor",&Kerker_factor,(double)4.0);

    /* scf.Mixing.StartPulay */

    input_int("scf.Mixing.StartPulay",&Pulay_SCF,6);

    /* close the file */
      
    input_close();
  }

}








void Output_Energies_Forces(FILE *fp)
{ 
  int ct_AN;
  double sumx,sumy,sumz;
  double AEx,AEy,AEz;
  
  fprintf(fp,"\n");
  fprintf(fp,"***********************************************************\n");
  fprintf(fp,"***********************************************************\n");
  fprintf(fp,"                     Energies (Hartree)                    \n");
  fprintf(fp,"***********************************************************\n");
  fprintf(fp,"***********************************************************\n");
  fprintf(fp,"    Uatom                         = %20.14f\n",Uatom);
  fprintf(fp,"    Uele for   up-spin (OS)       = %20.14f\n",Uele_OS0);
  fprintf(fp,"    Uele for   up-spin (IS)       = %20.14f\n",Uele_IS0);
  fprintf(fp,"    Uele for down-spin (OS)       = %20.14f\n",Uele_OS1);
  fprintf(fp,"    Uele for down-spin (IS)       = %20.14f\n",Uele_IS1);
  fprintf(fp,"    Uxc for up-spin               = %20.14f\n",Uxc0);
  fprintf(fp,"    Uxc for down-spin             = %20.14f\n",Uxc1);
  fprintf(fp,"    UH0                           = %20.14f\n",UH0);
  fprintf(fp,"    UH1                           = %20.14f\n",UH1);
  fprintf(fp,"    UH2                           = %20.14f\n",UH2);
  fprintf(fp,"    Ucore                         = %20.14f\n",Ucore);
  fprintf(fp,"    Udc  (Uxc+UH0+UH1+UH2+Ucore)  = %20.14f\n",Udc);
  fprintf(fp,"    Utot (Uele+Udc)               = %20.14f\n",Utot);
  fprintf(fp,"    Ucoh (Utot-Uatom)             = %20.14f\n",Ucoh);

  fprintf(fp,"\n");
  fprintf(fp,"***********************************************************\n");
  fprintf(fp,"***********************************************************\n");
  fprintf(fp,"                  Energies/atom (Hartree)                  \n");
  fprintf(fp,"***********************************************************\n");
  fprintf(fp,"***********************************************************\n");
  fprintf(fp,"    Uatom                         = %20.14f\n",Uatom/atomnum);
  fprintf(fp,"    Uele for   up-spin (OS)       = %20.14f\n",Uele_OS0/atomnum);
  fprintf(fp,"    Uele for   up-spin (IS)       = %20.14f\n",Uele_IS0/atomnum);
  fprintf(fp,"    Uele for down-spin (OS)       = %20.14f\n",Uele_OS1/atomnum);
  fprintf(fp,"    Uele for down-spin (IS)       = %20.14f\n",Uele_IS1/atomnum);
  fprintf(fp,"    Uxc for up-spin               = %20.14f\n",Uxc0/atomnum);
  fprintf(fp,"    Uxc for down-spin             = %20.14f\n",Uxc1/atomnum);
  fprintf(fp,"    UH0                           = %20.14f\n",UH0/atomnum);
  fprintf(fp,"    UH1                           = %20.14f\n",UH1/atomnum);
  fprintf(fp,"    UH2                           = %20.14f\n",UH2/atomnum);
  fprintf(fp,"    Ucore                         = %20.14f\n",Ucore/atomnum);
  fprintf(fp,"    Udc  (Uxc+UH0+UH1+UH2+Ucore)  = %20.14f\n",Udc/atomnum);
  fprintf(fp,"    Utot (Uele+Udc)               = %20.14f\n",Utot/atomnum);
  fprintf(fp,"    Ucoh (Utot-Uatom)             = %20.14f\n",Ucoh/atomnum);

  fprintf(fp,"\n");
  fprintf(fp,"***********************************************************\n");
  fprintf(fp,"***********************************************************\n");
  fprintf(fp,"               Force on Atom (Hartree/bohr)                \n");
  fprintf(fp,"***********************************************************\n");
  fprintf(fp,"***********************************************************\n");
  fprintf(fp,"\n");

  fprintf(fp,"                   Fx            Fy            Fz\n");

  sumx = 0.0;
  sumy = 0.0;
  sumz = 0.0;

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    fprintf(fp,"   %i %i        %12.9f  %12.9f  %12.9f\n",
            ct_AN,WhatSpecies[ct_AN],
            -Gxyz[ct_AN][17],-Gxyz[ct_AN][18],-Gxyz[ct_AN][19]);

    sumx = sumx - Gxyz[ct_AN][17];
    sumy = sumy - Gxyz[ct_AN][18];
    sumz = sumz - Gxyz[ct_AN][19];
  }
  fprintf(fp,"\n");
  fprintf(fp,"   Sum of F   %12.9f  %12.9f  %12.9f\n",sumx,sumy,sumz);

  /****************************************************
                   Correction of Force
  ****************************************************/

  AEx = sumx/(double)atomnum;
  AEy = sumy/(double)atomnum;
  AEz = sumz/(double)atomnum;
 
  fprintf(fp,"\n");
  fprintf(fp,"***********************************************************\n");
  fprintf(fp,"***********************************************************\n");
  fprintf(fp,"            Corrected Force on Atom (Hartree/bohr)         \n");
  fprintf(fp,"***********************************************************\n");
  fprintf(fp,"***********************************************************\n");
  fprintf(fp,"\n");

  fprintf(fp,"                   Fx            Fy            Fz\n");

  sumx = 0.0;
  sumy = 0.0;
  sumz = 0.0;

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    Gxyz[ct_AN][17] = Gxyz[ct_AN][17] + AEx;
    Gxyz[ct_AN][18] = Gxyz[ct_AN][18] + AEy;
    Gxyz[ct_AN][19] = Gxyz[ct_AN][19] + AEz;

    fprintf(fp,"   %i %i        %12.9f  %12.9f  %12.9f\n",
            ct_AN,WhatSpecies[ct_AN],
            -Gxyz[ct_AN][17],-Gxyz[ct_AN][18],-Gxyz[ct_AN][19]);

    sumx = sumx - Gxyz[ct_AN][17];
    sumy = sumy - Gxyz[ct_AN][18];
    sumz = sumz - Gxyz[ct_AN][19];
  }
  fprintf(fp,"\n");
  fprintf(fp,"   Sum of F   %12.9f  %12.9f  %12.9f\n",sumx,sumy,sumz);

}
