/**********************************************************************
  truncation.c:

     truncation.c is a subrutine to divide a large system
     to small systems and set grid data.

  Log of truncation.c:

     22/Nov/2001  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "openmx_common.h" 

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

#ifdef noomp
#include "mimic_omp.h"
#else
#include "omp.h"
#endif

#include "tran_prototypes.h"

#define Scale_BCR_DPCC    2.1    /* a cutoff factor for BCR_DPCC */
#define rcut_buffer       0.0
#define Min_Bond_Length   1.0e-12

static void Fixed_FNAN_SNAN();
static void LT(int CorL);
static void Reflexive(int Iatom, int Catom, int a, int b, int c, int Chops,int Fhops);
static int Get_R_index(int R_switch,int i,int j,int Rn);
static void Set_R_index(int R_switch,int i,int j,int Rn);
static int Set_Periodic(int CpyN, int Allocate_switch);
static void Free_truncation(int CpyN, int TN, int Free_switch);
static void free_arrays_truncation0();
static void Trn_System(int MD_iter, int CpyCell, int TCpyCell, int flag_GDC, int flag_GDC2);
static void Estimate_Trn_System(int CpyCell, int TCpyCell, int flag_GDC);
static void Check_System();
static void Set_RMI();
static void Output_Connectivity(FILE *fp);
static void UCell_Box(int MD_iter, int estimate_switch, int CpyCell);
static void Set_Inf_SndRcv();
static void allocate_grids2atoms(int MD_iter);

int TFNAN,TFNAN2,TSNAN,TSNAN2;
double **CellDis;
unsigned long asize1;
long *R_index1,*R_index2;



double truncation(int MD_iter, int flag_GDC, int UCell_flag)
{
  static int firsttime=1;
  int i,j,m,l,ct_AN,h_AN,Gh_AN,Mc_AN,Gc_AN,s1,s2;
  int k,tno,tno0,tno1,Cwan,Hwan,N,Nmax,so,nc,ns,spin;
  int num,wan,n2,wanA,Gi,Mc_AN_GDC,Max_Num_Cells0;
  int vsize,Anum,Bnum,NUM,p,MAnum,fan;
  int size_Orbs_Grid,size_COrbs_Grid;
  int size_H0,size_CntH0,size_H,size_CntH;
  int size_HNL,size_Left_U0;
  int size_iHNL,size_iHNL0,size_iCntHNL;
  int size_DS_NL,size_CntDS_NL,size_IOLP;
  int size_NumOLG,size_OLP,size_CntOLP;
  int size_OLP_L;
  int size_HVNA,size_DS_VNA,size_CntDS_VNA;
  int size_HVNA2,size_CntHVNA2;
  int size_H_Hub,size_DM_onsite;
  int size_v_eff,size_NC_OcpN,size_NC_v_eff;
  int size_S12,size_iDM;
  int size_ResidualDM,size_iResidualDM;
  int size_Krylov_U,size_EC_matrix;
  int size_H_Zeeman_NCO;
  int size_TRAN_DecMulP;
  int My_Max_GridN_Atom,My_Max_OneD_Grids,My_Max_NumOLG;
  int rlmax2,wanB;
  int numprocs,myid,tag=999,ID;
  double rcutA,r0,scale_rc0;
  int po,po0;
  double time0;

  double My_KryDH,My_KryDS;
  double KryDH,KryDS;

  char file_TRN[YOUSO10] = ".TRN";
  double TStime,TEtime;
  FILE *fp_TRN;
  char buf[fp_bsize];          /* setvbuf */

  /* MPI */
  if (atomnum<=MYID_MPI_COMM_WORLD) return 0.0;
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (myid==Host_ID && 0<level_stdout){
    printf("\n*******************************************************\n"); 
    printf("            Truncation and setting of grids            \n");
    printf("*******************************************************\n\n"); 
    printf("<truncation>  Logically truncation of the whole system\n");
  }

  dtime(&TStime); 

  /****************************************************
     freeing of dynamically allocated arrays using
     FNAN and SNAN at the previous MD step
  ****************************************************/

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    ID = G2ID[Gc_AN];
    MPI_Bcast(&time_per_atom[Gc_AN], 1, MPI_DOUBLE, ID, mpi_comm_level1);
  }

   free_arrays_truncation0();
  
  /****************************************************
            allocation of atoms to processors
  ****************************************************/

  if (2<=MD_iter && flag_GDC==0){

    /*****************************
      final input
      0: atomnum, 
      1: the neighbor 
      2: elapsed time 
    *****************************/    

    /* DC or Krylov */
    if (Solver==5 || Solver==8)
      Set_Allocate_Atom2CPU(MD_iter,0,2);  /* 2: elapsed time */
    else 
      Set_Allocate_Atom2CPU(MD_iter,0,0);  /* 0: atomnum */
  }

  else if (2<=MD_iter && flag_GDC==1){ /* for GDC */

    /************************************************
    the way of allocation of atoms to cpus in the GDC

    (1) find FNAN, SNAN, natn, ncn, Dis
    (2) Set_Allocate_Atom2CPU in truncation.c
    (3) free arrays allocated by Set_Periodic 
    ************************************************/

    TCpyCell = Set_Periodic(CpyCell,0);
    Estimate_Trn_System(CpyCell,TCpyCell,0);
    Allocate_Arrays(3);
    Trn_System(MD_iter,CpyCell,TCpyCell,0,1);
    Set_Allocate_Atom2CPU(MD_iter,0,2);
    Free_truncation(CpyCell,TCpyCell,0);
  }

  /****************************************************
                       Truncation
  ****************************************************/

  if (MD_iter==1){

    /****************************************************
     find:

       List_YOUSO[4] 

     Note: 

       In Set_Periodic(CpyCell,0)
         allocation of arrays:
            ratv           (global)
            atv            (global
            atv_ijk        (global)
            R_index1       (local)
            R_index2       (local)
            CellDis        (local)
            and
        Generation_ATV(CpyN);

       In Set_Periodic(CpyCell,1)
         allocation of arrays:
            R_index1       (local)
            R_index2       (local) 
            CellDis        (local) 

       In Allocate_Arrays(3) 
         freeing and allocation of arrays:
            natn           (global)
            ncn            (global)
            Dis            (global)

       In Free_truncation(CpyCell,TCpyCell,0)
         freeing of arrays:
            ratv           (global)
            atv            (global)
            atv_ijk        (global)
            R_index1       (local)  
            R_index2       (local) 
            CellDis        (local)

       In Free_truncation(CpyCell,TCpyCell,1)
         freeing of arrays:
            R_index1       (local) 
            R_index2       (local) 
            CellDis        (local)
    ****************************************************/

    CpyCell = 0;
    po = 0;
    TFNAN = 0;
    TSNAN = 0;

    do{

      CpyCell++;

      /**********************************************
           allocation of arrays listed above and 
                 Generation_ATV(CpyN);
      **********************************************/

      TCpyCell = Set_Periodic(CpyCell,0);

      /**********************************************
        find Max_FSNAN by the physical truncation
          for allocation of natn, ncn, and Dis 
      **********************************************/

      Estimate_Trn_System(CpyCell,TCpyCell,0);

      /**********************************************
              allocation of natn, ncn, and Dis 
      **********************************************/

      Allocate_Arrays(3);

      /**********************
       find TFNAN and TSNAN
      **********************/

      TFNAN2 = TFNAN;
      TSNAN2 = TSNAN;

      Trn_System(MD_iter,CpyCell,TCpyCell,0,1);

      if ( TFNAN==TFNAN2 && TSNAN==TSNAN2 && (Solver==1 || Solver==5 || Solver==6 || Solver==8) ) po++;
      else if (TFNAN==TFNAN2 && (Solver==2 || Solver==3 || Solver==4 || Solver==7 || Solver==9) ) po++;
      else if (CellNN_flag==1)                                                                    po++;

      /**********************************************
         if (flag_GDC==1)
      **********************************************/

      if (flag_GDC==1)  Set_Allocate_Atom2CPU(MD_iter,0,2);

      /**********************************************
         freeing of arrays which are allocated in
                 Set_Periodic(CpyCell,0)
      **********************************************/

      Free_truncation(CpyCell,TCpyCell,0);

    } while (po==0);

    List_YOUSO[4] = CpyCell;
    TCpyCell = Set_Periodic(CpyCell,0);
    Estimate_Trn_System(CpyCell,TCpyCell,flag_GDC);
    Allocate_Arrays(3);
    Trn_System(MD_iter,CpyCell,TCpyCell,flag_GDC,1);

    /* for GDC */
    if (flag_GDC==1)  Set_Allocate_Atom2CPU(MD_iter,0,2);

    Set_Inf_SndRcv();
    Set_RMI();

  } /* if (MD_iter==1) */ 
    
  else{

    TCpyCell = Set_Periodic(CpyCell,0);
    Estimate_Trn_System(CpyCell,TCpyCell,flag_GDC);
    Allocate_Arrays(3);
    Trn_System(MD_iter,CpyCell,TCpyCell,flag_GDC,1);

    /* for GDC */
    if (flag_GDC==1)  Set_Allocate_Atom2CPU(MD_iter,0,2);

    Set_Inf_SndRcv();
    Set_RMI();
  }

  if (2<=level_stdout){
    printf("List_YOUSO[4]=%2d\n",List_YOUSO[4]);
  } 

  /****************************************************
    allocation of arrays:

    NumOLG
  ****************************************************/

  FNAN[0] = 0;
  NumOLG = (int**)malloc(sizeof(int*)*(Matomnum+1));
  size_NumOLG = 0;
  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
    if (Mc_AN==0) Gc_AN = 0;
    else          Gc_AN = M2G[Mc_AN];
    NumOLG[Mc_AN] = (int*)malloc(sizeof(int)*(FNAN[Gc_AN]+1));
    size_NumOLG += FNAN[Gc_AN] + 1;
  }
  alloc_first[5] = 0;
  
  /* PrintMemory */
  if (Solver==2){
    if (firsttime)
    PrintMemory("SetPara_DFT: S",
                sizeof(double)*(Size_Total_Matrix+2)*(Size_Total_Matrix+2),NULL);
  }

  if (firsttime)
  PrintMemory("truncation: NumOLG", sizeof(int)*size_NumOLG, NULL);

  /****************************************************
                  check the system type
  ****************************************************/

  if (MD_iter==1) Check_System();

  /****************************************************
                       UCell_Box
  ****************************************************/

  if (MD_iter==1 && UCell_flag==1){

    /*************************************
      find 
            Max_GridN_Atom
            Max_OneD_Grids
    *************************************/

    Max_GridN_Atom = 0;
    Max_OneD_Grids = 0;
    if (2<=level_stdout) printf("\n***** UCell_Box(MD_iter,1,CpyCell) *****\n");

    UCell_Box(MD_iter,1,CpyCell);

    My_Max_GridN_Atom = Max_GridN_Atom;

    MPI_Reduce(&My_Max_GridN_Atom, &Max_GridN_Atom, 1,
                MPI_INT, MPI_MAX, Host_ID, mpi_comm_level1);
    MPI_Bcast(&Max_GridN_Atom, 1, MPI_INT, Host_ID, mpi_comm_level1);

    My_Max_OneD_Grids = Max_OneD_Grids;
    MPI_Reduce(&My_Max_OneD_Grids, &Max_OneD_Grids, 1,
                MPI_INT, MPI_MAX, Host_ID, mpi_comm_level1);
    MPI_Bcast(&Max_OneD_Grids, 1, MPI_INT, Host_ID, mpi_comm_level1);
    List_YOUSO[11] = (int)(Max_GridN_Atom*ScaleSize) + 1;
    List_YOUSO[17] = (int)(Max_OneD_Grids*ScaleSize);

    if (2<=level_stdout){
      printf("Max_OneD_Grids=%2d\n",Max_OneD_Grids);
    }

    /*************************************
      find 
            Max_NumOLG
    *************************************/

    Max_NumOLG = 0;
    if (2<=level_stdout) printf("\n***** UCell_Box(MD_iter,2,CpyCell) *****\n");

    UCell_Box(MD_iter,2,CpyCell);

    My_Max_NumOLG = Max_NumOLG;
    MPI_Reduce(&My_Max_NumOLG, &Max_NumOLG, 1,
                MPI_INT, MPI_MAX, Host_ID, mpi_comm_level1);
    MPI_Bcast(&Max_NumOLG, 1, MPI_INT, Host_ID, mpi_comm_level1);
    List_YOUSO[12] = (int)(Max_NumOLG*ScaleSize) + 1;

    if (2<=level_stdout){
      printf("YOUSO11=%i YOUSO12=%i YOUSO17=%i\n",
             List_YOUSO[11],List_YOUSO[12],List_YOUSO[17]);
    }

    /*************************************
      setting of
                 GListTAtoms0
                 GListTAtoms1
                 GListTAtoms2
                 GListTAtoms3 
                 GListTCells0
    *************************************/

    if (myid==Host_ID && 0<level_stdout)  printf("<UCell_Box> Info. of cutoff energy and num. of grids\n");

    UCell_Box(MD_iter,0,CpyCell);

    My_Max_GridN_Atom = Max_GridN_Atom;
    MPI_Reduce(&My_Max_GridN_Atom, &Max_GridN_Atom, 1,
                MPI_INT, MPI_MAX, Host_ID, mpi_comm_level1);
    MPI_Bcast(&Max_GridN_Atom, 1, MPI_INT, Host_ID, mpi_comm_level1);

    My_Max_OneD_Grids = Max_OneD_Grids;
    MPI_Reduce(&My_Max_OneD_Grids, &Max_OneD_Grids, 1,
                MPI_INT, MPI_MAX, Host_ID, mpi_comm_level1);
    MPI_Bcast(&Max_OneD_Grids, 1, MPI_INT, Host_ID, mpi_comm_level1);
    List_YOUSO[11] = (int)(Max_GridN_Atom*ScaleSize) + 1;
    List_YOUSO[17] = (int)(Max_OneD_Grids*ScaleSize);
  }

  else if (UCell_flag==1) {

    if (myid==Host_ID && 0<level_stdout) printf("<UCell_Box> Info. of cutoff energy and num. of grids\n");

    UCell_Box(MD_iter,0,CpyCell);

    List_YOUSO[11] = (int)(Max_GridN_Atom*ScaleSize) + 1;
    List_YOUSO[17] = (int)(Max_OneD_Grids*ScaleSize);
  }

  /****************************************************
    allocation of arrays:

     H0
     CntH0
     OLP
     CntOLP
     H
     CntH
     iCntHNL
     DS_NL
     CntDS_NL
     DM
     ResidualDM
     EDM
     PDM
     IOLP 
     double CntCoes[Matomnum+MatomnumF+1][YOUSO7][YOUSO24];
     HVNA
     DS_VNA
     CntDS_VNA
     HVNA2
     CntHVNA2
     H_Hub
     DM_onsite
     v_eff
     NC_OcpN
     NC_v_eff
  ****************************************************/

  /* H0 */  

  size_H0 = 0;
  H0 = (double*****)malloc(sizeof(double****)*4);
  for (k=0; k<4; k++){
    H0[k] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
    FNAN[0] = 0;
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

      if (Mc_AN==0){
        Gc_AN = 0;
        tno0 = 1;
      }
      else{
        Gc_AN = M2G[Mc_AN];
        Cwan = WhatSpecies[Gc_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    

      H0[k][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

        if (Mc_AN==0){
          tno1 = 1;  
        }
        else{
          Gh_AN = natn[Gc_AN][h_AN];
          Hwan = WhatSpecies[Gh_AN];
          tno1 = Spe_Total_NO[Hwan];
        } 

        H0[k][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
        for (i=0; i<tno0; i++){
          H0[k][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
          size_H0 += tno1;
        }
      }
    }
  }

  /* CntH0 */  

  size_CntH0 = 0;

  if (Cnt_switch==1){

    CntH0 = (double*****)malloc(sizeof(double****)*4);
    for (k=0; k<4; k++){
      CntH0[k] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
      FNAN[0] = 0;
      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	if (Mc_AN==0){
	  Gc_AN = 0;
	  tno0 = 1;
	}
	else{
	  Gc_AN = M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_CNO[Cwan];  
	}    

	CntH0[k][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_CNO[Hwan];
	  } 

	  CntH0[k][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
	  for (i=0; i<tno0; i++){
	    CntH0[k][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
	    size_CntH0 += tno1;
	  }
	}
      }
    }
  }

  /* HNL */  

  size_HNL = 0;
  HNL = (double*****)malloc(sizeof(double****)*List_YOUSO[5]);
  for (k=0; k<List_YOUSO[5]; k++){
    HNL[k] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
    FNAN[0] = 0;
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

      if (Mc_AN==0){
        Gc_AN = 0;
        tno0 = 1;
      }
      else{
        Gc_AN = M2G[Mc_AN];
        Cwan = WhatSpecies[Gc_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    

      HNL[k][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

        if (Mc_AN==0){
          tno1 = 1;  
        }
        else{
          Gh_AN = natn[Gc_AN][h_AN];
          Hwan = WhatSpecies[Gh_AN];
          tno1 = Spe_Total_NO[Hwan];
        } 

        HNL[k][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
        for (i=0; i<tno0; i++){
          HNL[k][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
          size_HNL += tno1;
	  for (j=0; j<tno1; j++) HNL[k][Mc_AN][h_AN][i][j] = 0.0;
        }
      }
    }
  }

  /* iHNL */  

  if ( SpinP_switch==3 ){

    size_iHNL = 0;
    iHNL = (double*****)malloc(sizeof(double****)*List_YOUSO[5]);
    for (k=0; k<List_YOUSO[5]; k++){
      iHNL[k] = (double****)malloc(sizeof(double***)*(Matomnum+MatomnumF+MatomnumS+1)); 
      FNAN[0] = 0;
      for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF+MatomnumS); Mc_AN++){

	if (Mc_AN==0){
	  Gc_AN = 0;
	  tno0 = 1;
	}
	else{
	  Gc_AN = S_M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	iHNL[k][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_NO[Hwan];
	  } 

	  iHNL[k][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
	  for (i=0; i<tno0; i++){
	    iHNL[k][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
  	    for (j=0; j<tno1; j++) iHNL[k][Mc_AN][h_AN][i][j] = 0.0;
	  }
          size_iHNL += tno0*tno1;  
	}
      }
    }
  }

  /* iCntHNL */  

  if (SO_switch==1 && Cnt_switch==1){

    size_iCntHNL = 0;
    iCntHNL = (double*****)malloc(sizeof(double****)*List_YOUSO[5]);
    for (k=0; k<List_YOUSO[5]; k++){
      iCntHNL[k] = (double****)malloc(sizeof(double***)*(Matomnum+MatomnumF+MatomnumS+1)); 
      FNAN[0] = 0;
      for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF+MatomnumS); Mc_AN++){

	if (Mc_AN==0){
	  Gc_AN = 0;
	  tno0 = 1;
	}
	else{
	  Gc_AN = S_M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_CNO[Cwan];  
	}    

	iCntHNL[k][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_CNO[Hwan];
	  } 

	  iCntHNL[k][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
	  for (i=0; i<tno0; i++){
	    iCntHNL[k][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
            size_iCntHNL += tno1;  
	  }
	}
      }
    }
  }

  /* H_Hub  --- added by MJ */  

  if (Hub_U_switch==1 || Constraint_NCS_switch==1 || Zeeman_NCS_switch==1 || Zeeman_NCO_switch==1){

    size_H_Hub = 0;

    H_Hub = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
    for (k=0; k<=SpinP_switch; k++){
      H_Hub[k] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
      FNAN[0] = 0;
      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	if (Mc_AN==0){
	  Gc_AN = 0;
	  tno0 = 1;
	}
	else{
	  Gc_AN = M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	H_Hub[k][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_NO[Hwan];
	  } 

	  H_Hub[k][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
	  for (i=0; i<tno0; i++){
	    H_Hub[k][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
  	    for (j=0; j<tno1; j++) H_Hub[k][Mc_AN][h_AN][i][j] = 0.0;
	  }

          size_H_Hub += tno0*tno1;

	}
      }
    }
  }

  /* H_Zeeman_NCO */  

  if (Zeeman_NCO_switch==1){

    size_H_Zeeman_NCO = 0;

    H_Zeeman_NCO = (double***)malloc(sizeof(double**)*(Matomnum+1)); 
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

      if (Mc_AN==0){
	Gc_AN = 0;
	tno0 = 1;
      }
      else{
	Gc_AN = M2G[Mc_AN];
	Cwan = WhatSpecies[Gc_AN];
	tno0 = Spe_Total_NO[Cwan];  
      }    

      H_Zeeman_NCO[Mc_AN] = (double**)malloc(sizeof(double*)*tno0); 
      for (i=0; i<tno0; i++){
        H_Zeeman_NCO[Mc_AN][i] = (double*)malloc(sizeof(double)*tno0); 
	for (j=0; j<tno0; j++) H_Zeeman_NCO[Mc_AN][i][j] = 0.0;
      }

      size_H_Zeeman_NCO += tno0*tno0;
    }
  }

  /* iHNL0 */

  size_iHNL0 = 0;

  if (SpinP_switch==3){

    iHNL0 = (double*****)malloc(sizeof(double****)*List_YOUSO[5]);
    for (k=0; k<List_YOUSO[5]; k++){
      iHNL0[k] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
      FNAN[0] = 0;
      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	if (Mc_AN==0){
	  Gc_AN = 0;
	  tno0 = 1;
	}
	else{
	  Gc_AN = M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	iHNL0[k][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_NO[Hwan];
	  } 

	  iHNL0[k][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
	  for (i=0; i<tno0; i++){
	    iHNL0[k][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
  	    for (j=0; j<tno1; j++) iHNL0[k][Mc_AN][h_AN][i][j] = 0.0;
	  }
          size_iHNL0 += tno0*tno1;  
	}
      }
    }
  }

  /* OLP_L */  

  size_OLP_L = 0;

  OLP_L = (double*****)malloc(sizeof(double****)*3); 
  for (k=0; k<3; k++){

    OLP_L[k] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
    FNAN[0] = 0;
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

      if (Mc_AN==0){
	Gc_AN = 0;
	tno0 = 1;
      }
      else{
	Gc_AN = F_M2G[Mc_AN];
	Cwan = WhatSpecies[Gc_AN];
	tno0 = Spe_Total_NO[Cwan];  
      }    

      OLP_L[k][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	if (Mc_AN==0){
	  tno1 = 1;  
	}
	else{
	  Gh_AN = natn[Gc_AN][h_AN];
	  Hwan = WhatSpecies[Gh_AN];
	  tno1 = Spe_Total_NO[Hwan];
	} 

	OLP_L[k][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
	for (i=0; i<tno0; i++){
	  OLP_L[k][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
	  size_OLP_L += tno1;
	}
      }
    }
  }

  /* OLP */  

  OLP = (double*****)malloc(sizeof(double****)*4);
  size_OLP = 0;
  for (k=0; k<4; k++){
    OLP[k] = (double****)malloc(sizeof(double***)*(Matomnum+MatomnumF+MatomnumS+1)); 
    FNAN[0] = 0;
    for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF+MatomnumS); Mc_AN++){

      if (Mc_AN==0){
        Gc_AN = 0;
        tno0 = 1;
      }
      else{
        Gc_AN = S_M2G[Mc_AN];
        Cwan = WhatSpecies[Gc_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    

      OLP[k][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

        if (Mc_AN==0){
          tno1 = 1;  
        }
        else{
          Gh_AN = natn[Gc_AN][h_AN];
          Hwan = WhatSpecies[Gh_AN];
          tno1 = Spe_Total_NO[Hwan];
        } 

        OLP[k][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
        for (i=0; i<tno0; i++){
          OLP[k][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
          size_OLP += tno1;
        }
      }
    }
  }

  /* CntOLP */  

  size_CntOLP = 0;

  if (Cnt_switch==1){
 
    CntOLP = (double*****)malloc(sizeof(double****)*4);
    for (k=0; k<4; k++){
      CntOLP[k] = (double****)malloc(sizeof(double***)*(Matomnum+MatomnumF+MatomnumS+1)); 
      FNAN[0] = 0;
      for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF+MatomnumS); Mc_AN++){

	if (Mc_AN==0){
	  Gc_AN = 0;
	  tno0 = 1;
	}
	else{
	  Gc_AN = S_M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_CNO[Cwan];  
	}    

	CntOLP[k][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_CNO[Hwan];
	  } 

	  CntOLP[k][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
	  for (i=0; i<tno0; i++){
	    CntOLP[k][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
	    size_CntOLP += tno1;
	  }
	}
      }
    }
  }

  /* H */  

  size_H = 0;
  H = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
  for (k=0; k<=SpinP_switch; k++){
    H[k] = (double****)malloc(sizeof(double***)*(Matomnum+MatomnumF+MatomnumS+1)); 
    FNAN[0] = 0;
    for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF+MatomnumS); Mc_AN++){

      if (Mc_AN==0){
        Gc_AN = 0;
        tno0 = 1;
      }
      else{
        Gc_AN = S_M2G[Mc_AN];
        Cwan = WhatSpecies[Gc_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    

      H[k][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

        if (Mc_AN==0){
          tno1 = 1;  
        }
        else{
          Gh_AN = natn[Gc_AN][h_AN];
          Hwan = WhatSpecies[Gh_AN];
          tno1 = Spe_Total_NO[Hwan];
        } 

        H[k][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
        for (i=0; i<tno0; i++){
          H[k][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
          for (j=0; j<tno1; j++) H[k][Mc_AN][h_AN][i][j] = 0.0;
          size_H += tno1;  
        }
      }
    }
  }

  /* CntH */  

  size_CntH = 0;

  if (Cnt_switch==1){

    CntH = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
    for (k=0; k<=SpinP_switch; k++){
      CntH[k] = (double****)malloc(sizeof(double***)*(Matomnum+MatomnumF+MatomnumS+1)); 
      FNAN[0] = 0;
      for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF+MatomnumS); Mc_AN++){

	if (Mc_AN==0){
	  Gc_AN = 0;
	  tno0 = 1;
	}
	else{
	  Gc_AN = S_M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_CNO[Cwan];  
	}    

	CntH[k][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_CNO[Hwan];
	  } 

	  CntH[k][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
	  for (i=0; i<tno0; i++){
	    CntH[k][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
	    size_CntH += tno1;  
	  }
	}
      }
    }
  }

  /* DS_NL */  

  size_DS_NL = 0;
  DS_NL = (double******)malloc(sizeof(double*****)*(SO_switch+1));
  for (so=0; so<(SO_switch+1); so++){
    DS_NL[so] = (double*****)malloc(sizeof(double****)*4);
    for (k=0; k<4; k++){
      DS_NL[so][k] = (double****)malloc(sizeof(double***)*(Matomnum+MatomnumF+1)); 
      FNAN[0] = 0;
      for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){

	if (Mc_AN==0){
	  Gc_AN = 0;
	  tno0 = 1;
	}
	else{
	  Gc_AN = F_M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	DS_NL[so][k][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_VPS_Pro[Hwan] + 2;
	  } 

	  DS_NL[so][k][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
	  for (i=0; i<tno0; i++){
	    DS_NL[so][k][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
	    size_DS_NL += tno1;
	  }
	}
      }
    }
  }
 
  /* CntDS_NL */  

  size_CntDS_NL = 0;

  if (Cnt_switch==1){

    CntDS_NL = (double******)malloc(sizeof(double*****)*(SO_switch+1));
    for (so=0; so<(SO_switch+1); so++){
      CntDS_NL[so] = (double*****)malloc(sizeof(double****)*4); 
      for (k=0; k<4; k++){
	CntDS_NL[so][k] = (double****)malloc(sizeof(double***)*(Matomnum+MatomnumF+1)); 
	FNAN[0] = 0;
	for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	  }
	  else{
	    Gc_AN = F_M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_CNO[Cwan];  
	  }    

	  CntDS_NL[so][k][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	    if (Mc_AN==0){
	      tno1 = 1;  
	    }
	    else{
	      Gh_AN = natn[Gc_AN][h_AN];
	      Hwan = WhatSpecies[Gh_AN];
	      tno1 = Spe_Total_VPS_Pro[Hwan] + 2;
	    } 

	    CntDS_NL[so][k][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
	    for (i=0; i<tno0; i++){
	      CntDS_NL[so][k][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
	      size_CntDS_NL += tno1;
	    }
	  }
	}
      }
    }
  }

  /* DM */

  DM = (double******)malloc(sizeof(double*****)*List_YOUSO[16]); 
  for (m=0; m<List_YOUSO[16]; m++){
    DM[m] = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
    for (k=0; k<=SpinP_switch; k++){
      DM[m][k] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
      FNAN[0] = 0;
      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	if (Mc_AN==0){
          Gc_AN = 0;
	  tno0 = 1;
	}
	else{
          Gc_AN = M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	DM[m][k][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1));
	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_NO[Hwan];
	  } 

	  DM[m][k][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
	  for (i=0; i<tno0; i++){
	    DM[m][k][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
  	    for (j=0; j<tno1; j++) DM[m][k][Mc_AN][h_AN][i][j] = 0.0; 
	  }
	}
      }
    }
  }

  /* Partial_DM */

  if (cal_partial_charge==1){

    Partial_DM = (double*****)malloc(sizeof(double****)*2); 
    for (k=0; k<=1; k++){
      Partial_DM[k] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
      FNAN[0] = 0;
      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	if (Mc_AN==0){
          Gc_AN = 0;
	  tno0 = 1;
	}
	else{
          Gc_AN = M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	Partial_DM[k][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1));
	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_NO[Hwan];
	  } 

	  Partial_DM[k][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
	  for (i=0; i<tno0; i++){
	    Partial_DM[k][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
  	    for (j=0; j<tno1; j++) Partial_DM[k][Mc_AN][h_AN][i][j] = 0.0; 
	  }
	}
      }
    }
  }

  /* DM_onsite   --- MJ */  

  if (Hub_U_switch==1 || Constraint_NCS_switch==1 || Zeeman_NCS_switch==1 || Zeeman_NCO_switch==1){
  
    size_DM_onsite = 0;

    DM_onsite = (double*****)malloc(sizeof(double****)*2); 

    for (m=0; m<2; m++){

      DM_onsite[m] = (double****)malloc(sizeof(double***)*(SpinP_switch+1)); 
      for (k=0; k<=SpinP_switch; k++){

	DM_onsite[m][k] = (double***)malloc(sizeof(double**)*(Matomnum+1)); 
	FNAN[0] = 0;

	for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
	  
	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	  }
	  else{
	    Gc_AN = M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	  }  

	  h_AN = 0;
	      
	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_NO[Hwan];
	  } 

	  DM_onsite[m][k][Mc_AN] = (double**)malloc(sizeof(double*)*tno0); 

	  for (i=0; i<tno0; i++){
	    DM_onsite[m][k][Mc_AN][i] = (double*)malloc(sizeof(double)*tno1); 
	    for (j=0; j<tno1; j++) DM_onsite[m][k][Mc_AN][i][j] = 0.0;
	  }

	  size_DM_onsite += tno0*tno1;

	}
      }
    }
  }

  /*  v_eff  --- MJ */

  if (Hub_U_switch==1 || Constraint_NCS_switch==1 || Zeeman_NCS_switch==1 || Zeeman_NCO_switch==1){

    size_v_eff = 0;

    v_eff = (double****)malloc(sizeof(double***)*(SpinP_switch+1)); 
    for (k=0; k<=SpinP_switch; k++){

      v_eff[k] = (double***)malloc(sizeof(double**)*(Matomnum+MatomnumF+1)); 
      FNAN[0] = 0;

      for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
	  
	if (Mc_AN==0){
	  Gc_AN = 0;
	  tno0 = 1;
	}
	else{
	  Gc_AN = F_M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}  

	v_eff[k][Mc_AN] = (double**)malloc(sizeof(double*)*tno0); 

	for (i=0; i<tno0; i++){
	  v_eff[k][Mc_AN][i] = (double*)malloc(sizeof(double)*tno0); 
	  for (j=0; j<tno0; j++) v_eff[k][Mc_AN][i][j] = 0.0; 
	}

        size_v_eff += tno0*tno0;
      }
    }
  }

  /*  NC_OcpN */

  if ( (Hub_U_switch==1 || Constraint_NCS_switch==1 || Zeeman_NCS_switch==1 || Zeeman_NCO_switch==1) 
       && SpinP_switch==3 ){

    size_NC_OcpN = 0;

    NC_OcpN = (dcomplex******)malloc(sizeof(dcomplex*****)*2); 

    for (m=0; m<2; m++){

      NC_OcpN[m] = (dcomplex*****)malloc(sizeof(dcomplex****)*2); 
      for (s1=0; s1<2; s1++){
	NC_OcpN[m][s1] = (dcomplex****)malloc(sizeof(dcomplex***)*2);
	for (s2=0; s2<2; s2++){
	  NC_OcpN[m][s1][s2] = (dcomplex***)malloc(sizeof(dcomplex**)*(Matomnum+1));
	  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	    if (Mc_AN==0){
	      Gc_AN = 0;
	      tno0 = 1;
	    }
	    else{
	      Gc_AN = M2G[Mc_AN];
	      Cwan = WhatSpecies[Gc_AN];
	      tno0 = Spe_Total_NO[Cwan];  
	    }  

	    NC_OcpN[m][s1][s2][Mc_AN] = (dcomplex**)malloc(sizeof(dcomplex*)*tno0);

	    for (i=0; i<tno0; i++){
	      NC_OcpN[m][s1][s2][Mc_AN][i] = (dcomplex*)malloc(sizeof(dcomplex)*tno0);
	      for (j=0; j<tno0; j++)  NC_OcpN[m][s1][s2][Mc_AN][i][j] = Complex(0.0,0.0);
	    }

	    size_NC_OcpN += tno0*tno0;
	  }
	}
      }
    }           
  }

  /*  NC_v_eff */

  if ( (Hub_U_switch==1 
     || Constraint_NCS_switch==1 
     || Zeeman_NCS_switch==1 
     || Zeeman_NCO_switch==1) && SpinP_switch==3 ){

    size_NC_v_eff = 0;

    NC_v_eff = (dcomplex*****)malloc(sizeof(dcomplex****)*2); 
    for (s1=0; s1<2; s1++){

      NC_v_eff[s1] = (dcomplex****)malloc(sizeof(dcomplex***)*2);

      for (s2=0; s2<2; s2++){

	NC_v_eff[s1][s2] = (dcomplex***)malloc(sizeof(dcomplex**)*(Matomnum+MatomnumF+1));

	for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	  }
	  else{
	    Gc_AN = F_M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	  }  

	  NC_v_eff[s1][s2][Mc_AN] = (dcomplex**)malloc(sizeof(dcomplex*)*tno0);

	  for (i=0; i<tno0; i++){

	    NC_v_eff[s1][s2][Mc_AN][i] = (dcomplex*)malloc(sizeof(dcomplex)*tno0);
	    for (j=0; j<tno0; j++)  NC_v_eff[s1][s2][Mc_AN][i][j] = Complex(0.0,0.0);
	  }

	  size_NC_v_eff += tno0*tno0;
	}
      }
    }
  }

  /* ResidualDM */  

  size_ResidualDM = 0;
  if ( Mixing_switch==0 || Mixing_switch==1 || Mixing_switch==2 ){

    ResidualDM = (double******)malloc(sizeof(double*****)*List_YOUSO[16]); 
    for (m=0; m<List_YOUSO[16]; m++){
      ResidualDM[m] = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
      for (k=0; k<=SpinP_switch; k++){
	ResidualDM[m][k] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
	FNAN[0] = 0;
	for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	  }
	  else{
	    Gc_AN = M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	  }    

	  ResidualDM[m][k][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	    if (Mc_AN==0){
	      tno1 = 1;  
	    }
	    else{
	      Gh_AN = natn[Gc_AN][h_AN];
	      Hwan = WhatSpecies[Gh_AN];
	      tno1 = Spe_Total_NO[Hwan];
	    } 

	    ResidualDM[m][k][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
	    for (i=0; i<tno0; i++){
	      ResidualDM[m][k][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
  	      for (j=0; j<tno1; j++) ResidualDM[m][k][Mc_AN][h_AN][i][j] = 0.0;
	    }

            size_ResidualDM += tno0*tno1; 
	  }
	}
      }
    }
  }

  /* iResidualDM */

  size_iResidualDM = 0;

  if ( (Mixing_switch==0 || Mixing_switch==1 || Mixing_switch==2)
       && SpinP_switch==3 && ( SO_switch==1 || Hub_U_switch==1 || Constraint_NCS_switch==1 
        || Zeeman_NCS_switch==1 || Zeeman_NCO_switch==1) ){

    iResidualDM = (double******)malloc(sizeof(double*****)*List_YOUSO[16]); 
    for (m=0; m<List_YOUSO[16]; m++){
      iResidualDM[m] = (double*****)malloc(sizeof(double****)*2); 
      for (k=0; k<2; k++){
	iResidualDM[m][k] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
	FNAN[0] = 0;
	for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	  }
	  else{
	    Gc_AN = M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	  }    

	  iResidualDM[m][k][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	    if (Mc_AN==0){
	      tno1 = 1;  
	    }
	    else{
	      Gh_AN = natn[Gc_AN][h_AN];
	      Hwan = WhatSpecies[Gh_AN];
	      tno1 = Spe_Total_NO[Hwan];
	    } 

	    iResidualDM[m][k][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
	    for (i=0; i<tno0; i++){
	      iResidualDM[m][k][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
  	      for (j=0; j<tno1; j++) iResidualDM[m][k][Mc_AN][h_AN][i][j] = 0.0;
	    }

            size_iResidualDM += tno0*tno1; 
	  }
	}
      }
    }
  }

  else {
    iResidualDM = (double******)malloc(sizeof(double*****)*List_YOUSO[16]); 
    for (m=0; m<List_YOUSO[16]; m++){
      iResidualDM[m] = (double*****)malloc(sizeof(double****)*1); 
      iResidualDM[m][0] = (double****)malloc(sizeof(double***)*1); 
      iResidualDM[m][0][0] = (double***)malloc(sizeof(double**)*1); 
      iResidualDM[m][0][0][0] = (double**)malloc(sizeof(double*)*1); 
      iResidualDM[m][0][0][0][0] = (double*)malloc(sizeof(double)*1); 
    }   
  }

  /* EDM */  

  EDM = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
  for (k=0; k<=SpinP_switch; k++){
    EDM[k] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
    FNAN[0] = 0;
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

      if (Mc_AN==0){
        Gc_AN = 0;
        tno0 = 1;
      }
      else{
        Gc_AN = M2G[Mc_AN];
        Cwan = WhatSpecies[Gc_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    

      EDM[k][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

        if (Mc_AN==0){
          tno1 = 1;  
        }
        else{
          Gh_AN = natn[Gc_AN][h_AN];
          Hwan = WhatSpecies[Gh_AN];
          tno1 = Spe_Total_NO[Hwan];
        } 

        EDM[k][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
        for (i=0; i<tno0; i++){
          EDM[k][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
        }
      }
    }
  }

  /* PDM */  

  PDM = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
  for (k=0; k<=SpinP_switch; k++){
    PDM[k] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
    FNAN[0] = 0;
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

      if (Mc_AN==0){
        Gc_AN = 0;
        tno0 = 1;
      }
      else{
        Gc_AN = M2G[Mc_AN];
        Cwan = WhatSpecies[Gc_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    

      PDM[k][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

        if (Mc_AN==0){
          tno1 = 1;  
        }
        else{
          Gh_AN = natn[Gc_AN][h_AN];
          Hwan = WhatSpecies[Gh_AN];
          tno1 = Spe_Total_NO[Hwan];
        } 

        PDM[k][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
        for (i=0; i<tno0; i++){
          PDM[k][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
        }
      }
    }
  }

  /* iDM */  

  size_iDM = 0;
  iDM = (double******)malloc(sizeof(double*****)*List_YOUSO[16]);
  for (m=0; m<List_YOUSO[16]; m++){
    iDM[m] = (double*****)malloc(sizeof(double****)*2); 
    for (k=0; k<2; k++){
      iDM[m][k] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
      FNAN[0] = 0;
      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	if (Mc_AN==0){
	  Gc_AN = 0;
	  tno0 = 1;
	}
	else{
	  Gc_AN = M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	iDM[m][k][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_NO[Hwan];
	  } 

	  iDM[m][k][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
	  for (i=0; i<tno0; i++){
	    iDM[m][k][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
	    for (j=0; j<tno1; j++)  iDM[m][k][Mc_AN][h_AN][i][j] = 0.0; 
	  }
	  size_iDM += tno0*tno1;
	}
      }
    }
  }

  /* IOLP */  

  size_IOLP = 0;
  if (Solver==1 || Solver==7){
    IOLP = (double****)malloc(sizeof(double***)*(Matomnum+MatomnumF+MatomnumS+1)); 
    FNAN[0] = 0;
    SNAN[0] = 0;
    for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF+MatomnumS); Mc_AN++){

      if (Mc_AN==0){
	Gc_AN = 0;
	tno0 = 1;
      }
      else{
	Gc_AN = S_M2G[Mc_AN];
	Cwan = WhatSpecies[Gc_AN];
	tno0 = Spe_Total_NO[Cwan];  
      }

      IOLP[Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+SNAN[Gc_AN]+1)); 
      for (h_AN=0; h_AN<=(FNAN[Gc_AN]+SNAN[Gc_AN]); h_AN++){

	if (Mc_AN==0){
	  tno1 = 1;  
	}
	else{
	  Gh_AN = natn[Gc_AN][h_AN];
	  Hwan = WhatSpecies[Gh_AN];
	  tno1 = Spe_Total_NO[Hwan];
	} 

	IOLP[Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
	for (i=0; i<tno0; i++){
	  IOLP[Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
	  size_IOLP += tno1;
	}
      }
    }
  }

  /* S12 for recursion or DC */  

  if (Solver==1 || Solver==5){

    size_S12 = 0;

    S12 = (double***)malloc(sizeof(double**)*(Matomnum+1));

    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

      if (Mc_AN==0) n2 = 1;
      else{
        Gc_AN = M2G[Mc_AN];
        wan = WhatSpecies[Gc_AN];

	num = 1;
	for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
	  Gi = natn[Gc_AN][i];
	  wanA = WhatSpecies[Gi];
	  num += Spe_Total_CNO[wanA];
	}
	n2 = num + 2;
      }

      S12[Mc_AN] = (double**)malloc(sizeof(double*)*n2);
      for (i=0; i<n2; i++){
        S12[Mc_AN][i] = (double*)malloc(sizeof(double)*n2);
      }
      size_S12 += n2*n2;
    }
  }

  else if (Solver==6){ /* for GDC */

    size_S12 = 0;

    S12 = (double***)malloc(sizeof(double**)*(Matomnum_GDC+1));

    for (Mc_AN_GDC=0; Mc_AN_GDC<=Matomnum_GDC; Mc_AN_GDC++){

      if (Mc_AN_GDC==0) n2 = 1;
      else{

        Mc_AN = Mnatn_GDC[Mc_AN_GDC][0];
        Gc_AN = M2G[Mc_AN];
        wan = WhatSpecies[Gc_AN];

	num = 1;
	for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
	  Gi = natn[Gc_AN][i];
	  wanA = WhatSpecies[Gi];
	  num += Spe_Total_CNO[wanA];
	}
	n2 = num + 2;
      }

      S12[Mc_AN_GDC] = (double**)malloc(sizeof(double*)*n2);
      for (i=0; i<n2; i++){
        S12[Mc_AN_GDC][i] = (double*)malloc(sizeof(double)*n2);
      }
      size_S12 += n2*n2;
    }
  }

  if (Solver==1) { /* for recursion */

    /* Left_U0 */

    size_Left_U0 = 0; 

    Left_U0 = (double*****)malloc(sizeof(double****)*(SpinP_switch+1));
    for (i=0; i<=SpinP_switch; i++){
      Left_U0[i] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
      for (j=0; j<=Matomnum; j++){

	if (j==0){
          tno1 = 1;
          vsize = 1;
	}
        else{
          Gc_AN = M2G[j];
          wanA = WhatSpecies[Gc_AN];
          tno1 = Spe_Total_CNO[wanA];

          Anum = 1;
          for (p=0; p<=(FNAN[Gc_AN]+SNAN[Gc_AN]); p++){
  	    Gi = natn[Gc_AN][p];
	    wanA = WhatSpecies[Gi];
	    Anum += Spe_Total_CNO[wanA];
          }

          NUM = Anum - 1;
          vsize = NUM + 2;
	}

	Left_U0[i][j] = (double***)malloc(sizeof(double**)*List_YOUSO[3]); 
	for (k=0; k<List_YOUSO[3]; k++){
	  Left_U0[i][j][k] = (double**)malloc(sizeof(double*)*tno1); 
	  for (l=0; l<tno1; l++){
	    Left_U0[i][j][k][l] = (double*)malloc(sizeof(double)*vsize); 
	    for (m=0; m<vsize; m++) Left_U0[i][j][k][l][m] = 0.0; 
	  }
          size_Left_U0 += tno1*vsize; 
	}
      }
    }

    /* Right_U0 */

    Right_U0 = (double*****)malloc(sizeof(double****)*(SpinP_switch+1));
    for (i=0; i<=SpinP_switch; i++){
      Right_U0[i] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
      for (j=0; j<=Matomnum; j++){

	if (j==0){
          tno1 = 1;
          vsize = 1;
	}
        else{
          Gc_AN = M2G[j];
          wanA = WhatSpecies[Gc_AN];
          tno1 = Spe_Total_CNO[wanA];

          Anum = 1;
          for (p=0; p<=(FNAN[Gc_AN]+SNAN[Gc_AN]); p++){
  	    Gi = natn[Gc_AN][p];
	    wanA = WhatSpecies[Gi];
	    Anum += Spe_Total_CNO[wanA];
          }

          NUM = Anum - 1;
          vsize = NUM + 2;
	}

	Right_U0[i][j] = (double***)malloc(sizeof(double**)*List_YOUSO[3]); 
	for (k=0; k<List_YOUSO[3]; k++){
	  Right_U0[i][j][k] = (double**)malloc(sizeof(double*)*tno1); 
	  for (l=0; l<tno1; l++){
	    Right_U0[i][j][k][l] = (double*)malloc(sizeof(double)*vsize); 
	    for (m=0; m<vsize; m++) Right_U0[i][j][k][l][m] = 0.0;
	  }
	}
      }
    }
  }

  /* CntCoes */

  if (Cnt_switch==1){
    CntCoes = (double***)malloc(sizeof(double**)*(Matomnum+MatomnumF+1));
    for (i=0; i<=(Matomnum+MatomnumF); i++){
      CntCoes[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
      for (j=0; j<List_YOUSO[7]; j++){
	CntCoes[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[24]);
      }
    }

    CntCoes_Species = (double***)malloc(sizeof(double**)*(SpeciesNum+1));
    for (i=0; i<(SpeciesNum+1); i++){
      CntCoes_Species[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
      for (j=0; j<List_YOUSO[7]; j++){
	CntCoes_Species[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[24]);
      }
    }
  }

  if (ProExpn_VNA==1){

    /* HVNA */  

    size_HVNA = 0;
    HVNA = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
    FNAN[0] = 0;
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

      if (Mc_AN==0){
        Gc_AN = 0;
        tno0 = 1;
      }
      else{
        Gc_AN = M2G[Mc_AN];
        Cwan = WhatSpecies[Gc_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    

      HVNA[Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

        if (Mc_AN==0){
          tno1 = 1;  
        }
        else{
          Gh_AN = natn[Gc_AN][h_AN];
          Hwan = WhatSpecies[Gh_AN];
          tno1 = Spe_Total_NO[Hwan];
        } 

        HVNA[Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
        for (i=0; i<tno0; i++){
          HVNA[Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
          size_HVNA += tno1;
        }
      }
    }

    /* DS_VNA */

    size_DS_VNA = 0;
    DS_VNA = (Type_DS_VNA*****)malloc(sizeof(Type_DS_VNA****)*4); 
    for (k=0; k<4; k++){

      DS_VNA[k] = (Type_DS_VNA****)malloc(sizeof(Type_DS_VNA***)*(Matomnum+2)); 

      FNAN[0] = 0;
      for (Mc_AN=0; Mc_AN<(Matomnum+2); Mc_AN++){
          
	if (Mc_AN==0){
	  Gc_AN = 0;
	  tno0 = 1;
          fan = FNAN[Gc_AN];
	}
	else if ( (Matomnum+1)<=Mc_AN ){
          fan = List_YOUSO[8];
          tno0 = List_YOUSO[7];
        }
	else{
	  Gc_AN = F_M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
          fan = FNAN[Gc_AN];
	}
          
	DS_VNA[k][Mc_AN] = (Type_DS_VNA***)malloc(sizeof(Type_DS_VNA**)*(fan+1)); 

	for (h_AN=0; h_AN<(fan+1); h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    tno1 = (List_YOUSO[35]+1)*(List_YOUSO[35]+1)*List_YOUSO[34];
	  } 

	  DS_VNA[k][Mc_AN][h_AN] = (Type_DS_VNA**)malloc(sizeof(Type_DS_VNA*)*tno0); 
	  for (i=0; i<tno0; i++){
	    DS_VNA[k][Mc_AN][h_AN][i] = (Type_DS_VNA*)malloc(sizeof(Type_DS_VNA)*tno1);
	    size_DS_VNA += tno1;
	  }
	}
      }
    }

    /* CntDS_VNA */  

    if (Cnt_switch==1){

      size_CntDS_VNA = 0;
      CntDS_VNA = (Type_DS_VNA*****)malloc(sizeof(Type_DS_VNA****)*4); 
      for (k=0; k<4; k++){

	CntDS_VNA[k] = (Type_DS_VNA****)malloc(sizeof(Type_DS_VNA***)*(Matomnum+2));

	FNAN[0] = 0;
	for (Mc_AN=0; Mc_AN<(Matomnum+2); Mc_AN++){
          
	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
            fan = FNAN[Gc_AN];
	  }
	  else if ( Mc_AN==(Matomnum+1) ){
            fan = List_YOUSO[8];
            tno0 = List_YOUSO[7];
          }
	  else{
	    Gc_AN = F_M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_CNO[Cwan];  
            fan = FNAN[Gc_AN];
	  }
          
	  CntDS_VNA[k][Mc_AN] = (Type_DS_VNA***)malloc(sizeof(Type_DS_VNA**)*(fan+1)); 
	  for (h_AN=0; h_AN<(fan+1); h_AN++){

	    if (Mc_AN==0){
	      tno1 = 1;  
	    }
	    else{
	      tno1 = (List_YOUSO[35]+1)*(List_YOUSO[35]+1)*List_YOUSO[34];
	    } 

	    CntDS_VNA[k][Mc_AN][h_AN] = (Type_DS_VNA**)malloc(sizeof(Type_DS_VNA*)*tno0); 
	    for (i=0; i<tno0; i++){
	      CntDS_VNA[k][Mc_AN][h_AN][i] = (Type_DS_VNA*)malloc(sizeof(Type_DS_VNA)*tno1); 
	      size_CntDS_VNA += tno1;
	    }
	  }
	}
      }
    }

    /* HVNA2 */  

    size_HVNA2 = 0;
    HVNA2 = (double*****)malloc(sizeof(double****)*4);
    for (k=0; k<4; k++){
      HVNA2[k] = (double****)malloc(sizeof(double***)*(Matomnum+MatomnumF+1)); 
      FNAN[0] = 0;
      for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){

	if (Mc_AN==0){
	  Gc_AN = 0;
	  tno0 = 1;
	}
	else{
	  Gc_AN = F_M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	HVNA2[k][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  HVNA2[k][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
	  for (i=0; i<tno0; i++){
	    HVNA2[k][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno0); 
	    size_HVNA2 += tno0;
	  }
	}
      }
    }

    /* CntHVNA2 */  

    if (Cnt_switch==1){

      size_CntHVNA2 = 0;
      CntHVNA2 = (double*****)malloc(sizeof(double****)*4);
      for (k=0; k<4; k++){
	CntHVNA2[k] = (double****)malloc(sizeof(double***)*(Matomnum+MatomnumF+1)); 
	FNAN[0] = 0;
	for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	  }
	  else{
	    Gc_AN = F_M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_CNO[Cwan];  
	  }    

	  CntHVNA2[k][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	    CntHVNA2[k][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
	    for (i=0; i<tno0; i++){
	      CntHVNA2[k][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno0); 
	      size_CntHVNA2 += tno0;
	    }
	  }
	}
      }
    }

  }

  if (Solver==8) { /* Krylov subspace method */

    if (EKC_expand_core_flag==1)
      scale_rc0 = 1.10;
    else 
      scale_rc0 = 0.80;

    EKC_core_size = (int*)malloc(sizeof(int)*(Matomnum+1));
    scale_rc_EKC  = (double*)malloc(sizeof(double)*(Matomnum+1));

    for (i=1; i<=Matomnum; i++){

      scale_rc_EKC[i] = scale_rc0;

      ct_AN = M2G[i];
      wanA = WhatSpecies[ct_AN];
      rcutA = Spe_Atom_Cut1[wanA];

      /* find the nearest atom with distance of r0 */   

      r0 = 1.0e+10;
      for (h_AN=1; h_AN<=FNAN[ct_AN]; h_AN++){
        Gh_AN = natn[ct_AN][h_AN];
        wanB = WhatSpecies[Gh_AN];
        if (Dis[ct_AN][h_AN]<r0)  r0 = Dis[ct_AN][h_AN]; 
      }

      /* find atoms within scale_rc_EKC times r0 */

      po0 = 0;        

      do {
      
        Anum = 0;
        for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
          Gh_AN = natn[ct_AN][h_AN];
          wanB = WhatSpecies[Gh_AN];
          if ( Dis[ct_AN][h_AN]<(scale_rc_EKC[i]*r0) ){
            Anum += Spe_Total_CNO[wanB];
	  }
        }

        if (EKC_expand_core_flag==1){
          if (30<Anum || 1<=po0) po0 = 3;
          else                   scale_rc_EKC[i] *= 1.2;  
	}
        else{
     	  po0 = 3;
	}

	po0++;

      } while(po0<2);

      EKC_core_size[i] = Anum;  
    }

    rlmax_EC  = (int*)malloc(sizeof(int)*(Matomnum+1));
    rlmax_EC2 = (int*)malloc(sizeof(int)*(Matomnum+1));

    size_Krylov_U = 0;

    Krylov_U = (double*****)malloc(sizeof(double****)*(SpinP_switch+1));
    for (i=0; i<=SpinP_switch; i++){

      Krylov_U[i] = (double****)malloc(sizeof(double***)*(Matomnum+1));

      My_KryDH = 0.0;
      My_KryDS = 0.0;

      for (j=0; j<=Matomnum; j++){

	if (j==0){
	  Gc_AN = 0;
	  tno0 = 1;
          FNAN[0] = 1;
          rlmax_EC[0] = 1;
          Anum = 1; 
          Bnum = 0;
          EKC_core_size[0] = 1;
	}
	else{
	  Gc_AN = M2G[j];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  

          Anum = 0;  
          for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
            Gh_AN = natn[Gc_AN][h_AN];
            wanA = WhatSpecies[Gh_AN];
            Anum += Spe_Total_CNO[wanA];
	  }

          Bnum = 0;  
          for (h_AN=0; h_AN<=(FNAN[Gc_AN]+SNAN[Gc_AN]); h_AN++){
            Gh_AN = natn[Gc_AN][h_AN];
            wanA = WhatSpecies[Gh_AN];
            Bnum += Spe_Total_CNO[wanA];
	  }

          rlmax2 = (int)ceil((double)Bnum/(double)EKC_core_size[j]);

          if (KrylovH_order<Bnum){
            rlmax_EC[j] = (int)ceil((double)KrylovH_order/(double)EKC_core_size[j]);
	  }	
          else{
            rlmax_EC[j] = (int)ceil((double)Bnum/(double)EKC_core_size[j]);
          } 

          if (KrylovS_order<Bnum){
            rlmax_EC2[j] = (int)ceil((double)KrylovS_order/(double)EKC_core_size[j]);
	  }
          else{
            rlmax_EC2[j] = (int)ceil((double)Bnum/(double)EKC_core_size[j]);
          } 

          if (2<=level_stdout){
            printf("<Krylov parameters>  Gc_AN=%4d rlmax_EC=%3d rlmax_EC2=%3d EKC_core_size=%3d\n",
                    Gc_AN,rlmax_EC[j],rlmax_EC2[j],EKC_core_size[j]);
	  }

          My_KryDH += rlmax_EC[j]*EKC_core_size[j];
          My_KryDS += rlmax_EC2[j]*EKC_core_size[j];

	}    

        Krylov_U[i][j] = (double***)malloc(sizeof(double**)*rlmax_EC[j]);
        for (k=0; k<rlmax_EC[j]; k++){
          Krylov_U[i][j][k] = (double**)malloc(sizeof(double*)*EKC_core_size[j]);
          for (l=0; l<EKC_core_size[j]; l++){
            Krylov_U[i][j][k][l] = (double*)malloc(sizeof(double)*(Bnum+2));
	  }
        }

        size_Krylov_U += rlmax_EC[j]*EKC_core_size[j]*(Bnum+2);
      }

      if (i==0){
        MPI_Reduce(&My_KryDH, &KryDH, 1, MPI_DOUBLE, MPI_SUM, Host_ID, mpi_comm_level1);
        MPI_Reduce(&My_KryDS, &KryDS, 1, MPI_DOUBLE, MPI_SUM, Host_ID, mpi_comm_level1);

        if (myid==Host_ID && 2<=level_stdout){
	  printf("<Krylov parameters>  Av Krlov dimension H=%10.5f S=%10.5f\n",
                   KryDH/atomnum,KryDS/atomnum);
        }       
      }

    }

    size_EC_matrix = 0;

    EC_matrix = (double****)malloc(sizeof(double***)*(SpinP_switch+1));
    for (i=0; i<=SpinP_switch; i++){
      EC_matrix[i] = (double***)malloc(sizeof(double**)*(Matomnum+1));
      for (j=0; j<=Matomnum; j++){

	if (j==0){
	  tno0 = 1;
          EKC_core_size[0] = 1;
          rlmax_EC[0] = 1;
	}
	else{
	  Gc_AN = M2G[j];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

        EC_matrix[i][j] = (double**)malloc(sizeof(double*)*(rlmax_EC[j]*EKC_core_size[j]+1));

        for (k=0; k<(rlmax_EC[j]*EKC_core_size[j]+1); k++){
          EC_matrix[i][j][k] = (double*)malloc(sizeof(double)*(rlmax_EC[j]*EKC_core_size[j]+1));
	}

        size_EC_matrix += (rlmax_EC[j]*EKC_core_size[j]+1)*(rlmax_EC[j]*EKC_core_size[j]+1);
      }
    }

    /* find EKC_core_size_max */
        
    EKC_core_size_max = 0;
    for (j=1; j<=Matomnum; j++){
      if (EKC_core_size_max<EKC_core_size[j]) EKC_core_size_max = EKC_core_size[j];
    }
  }

  /* NEGF */

  if (Solver==4){
 
    size_TRAN_DecMulP = (SpinP_switch+1)*(Matomnum+1)*List_YOUSO[7];

    TRAN_DecMulP = (double***)malloc(sizeof(double**)*(SpinP_switch+1));
    for (spin=0; spin<(SpinP_switch+1); spin++){
      TRAN_DecMulP[spin] = (double**)malloc(sizeof(double*)*(Matomnum+1));
      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
	TRAN_DecMulP[spin][Mc_AN] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
	for (i=0; i<List_YOUSO[7]; i++) TRAN_DecMulP[spin][Mc_AN][i] = 0.0;
      }

    }
  }

  /* set zero */

  alloc_first[4] = 0;

  /* PrintMemory */

  if (firsttime){
  PrintMemory("truncation: H0",      sizeof(double)*size_H0,      NULL);
  PrintMemory("truncation: CntH0",   sizeof(double)*size_CntH0,   NULL);
  PrintMemory("truncation: HNL",     sizeof(double)*size_HNL,     NULL);
  PrintMemory("truncation: OLP",     sizeof(double)*size_OLP,     NULL);
  PrintMemory("truncation: CntOLP",  sizeof(double)*size_CntOLP,  NULL);
  PrintMemory("truncation: OLP_L",   sizeof(double)*size_OLP_L,   NULL);
  PrintMemory("truncation: H",       sizeof(double)*size_H,       NULL);
  PrintMemory("truncation: CntH",    sizeof(double)*size_CntH,    NULL);
  PrintMemory("truncation: DS_NL",   sizeof(double)*size_DS_NL,   NULL);
  PrintMemory("truncation: CntDS_NL",sizeof(double)*size_CntDS_NL,NULL);
  PrintMemory("truncation: DM",      sizeof(double)*List_YOUSO[16]*
                                                        size_H0,  NULL);
  PrintMemory("truncation: iDM",     sizeof(double)*size_iDM,     NULL);
  PrintMemory("truncation: ResidualDM",sizeof(double)*size_ResidualDM,  NULL);
  PrintMemory("truncation: EDM",     sizeof(double)*size_H0,      NULL);
  PrintMemory("truncation: PDM",     sizeof(double)*size_H0,      NULL);
  PrintMemory("truncation: IOLP",    sizeof(double)*size_IOLP,    NULL);
  if (SO_switch==1 || SpinP_switch==3){
  PrintMemory("truncation: iResidualDM",sizeof(double)*size_iResidualDM,  NULL);
  PrintMemory("truncation: iHNL",    sizeof(double)*size_iHNL,    NULL);
  PrintMemory("truncation: iHNL0",   sizeof(double)*size_iHNL0,   NULL);
  }
  if (SO_switch==1 && Cnt_switch==1){
  PrintMemory("truncation: iCntHNL", sizeof(double)*size_iCntHNL,     NULL);
  }
  if (ProExpn_VNA==1){
  PrintMemory("truncation: DS_VNA",  sizeof(Type_DS_VNA)*size_DS_VNA,      NULL);
  PrintMemory("truncation: HVNA",    sizeof(double)*size_HVNA,        NULL);
  PrintMemory("truncation: HVNA2",   sizeof(double)*size_HVNA2,       NULL);
  }
  if (Cnt_switch==1){
  PrintMemory("truncation: CntDS_VNA",sizeof(Type_DS_VNA)*size_CntDS_VNA,  NULL);
  PrintMemory("truncation: CntHVNA2", sizeof(double)*size_CntHVNA2,   NULL);
  }
  if (Hub_U_switch==1 || Constraint_NCS_switch==1 || Zeeman_NCS_switch==1 || Zeeman_NCO_switch==1){
  PrintMemory("truncation: H_Hub",     sizeof(double)*size_H_Hub,     NULL);
  PrintMemory("truncation: DM_onsite", sizeof(double)*size_DM_onsite, NULL);
  PrintMemory("truncation: v_eff",     sizeof(double)*size_v_eff,     NULL);
  }
  if (Zeeman_NCO_switch==1){
  PrintMemory("truncation: size_H_Zeeman_NCO",   sizeof(double)*size_H_Zeeman_NCO,     NULL);
  }
  if ( (Hub_U_switch==1 || Constraint_NCS_switch==1 || Zeeman_NCS_switch==1 || Zeeman_NCO_switch==1) 
        && SpinP_switch==3 ){
  PrintMemory("truncation: size_NC_OcpN",   sizeof(dcomplex)*size_NC_OcpN,     NULL);
  PrintMemory("truncation: size_NC_v_eff",  sizeof(dcomplex)*size_NC_v_eff,    NULL);
  }
  if (Solver==1 || Solver==5 || Solver==6){
  PrintMemory("truncation: S12",       sizeof(double)*size_S12,       NULL);
  }
  if (Solver==1){
  PrintMemory("truncation: Left_U0",       sizeof(double)*size_Left_U0,      NULL);
  PrintMemory("truncation: Right_U0",      sizeof(double)*size_Left_U0,      NULL);
  }
  if (Solver==8){
  PrintMemory("truncation: Krylov_U",      sizeof(double)*size_Krylov_U,      NULL);
  PrintMemory("truncation: EC_matrix",     sizeof(double)*size_EC_matrix,     NULL);
  }
  if (Solver==4){
  PrintMemory("truncation: TRAN_DecMulP",  sizeof(double)*(SpinP_switch+1)*(Matomnum+1)*List_YOUSO[7], NULL);
  }
  }

  /****************************************************
    allocation of arrays:

      Density_Grid
      ADensity_Grid
      PCCDensity_Grid
      Vxc_Grid
      RefVxc_Grid
      VNA_Grid
      dVHart_Grid
      Vpot_Grid
      Orbs_Grid
      COrbs_Grid
  ****************************************************/

  if (UCell_flag==1){

    MPI_Reduce(&Num_Cells0, &Max_Num_Cells0, 1, MPI_INT, MPI_MAX, Host_ID, mpi_comm_level1);
    MPI_Bcast(&Max_Num_Cells0, 1, MPI_INT, Host_ID, mpi_comm_level1);

    N = Num_Cells0*Ngrid2*Ngrid3;
    Nmax = Max_Num_Cells0*Ngrid2*Ngrid3;

    /******************************************************************************
     Since Density_Grid is used in TRAN_FFT_Dinterpolation3D.c as a work array,
     the size of the second index must be Nmax in the NEGF calculation (Solver==4). 
    ******************************************************************************/

    /*
    printf("myid=%2d Num_Cells0*Ngrid2*Ngrid3=%2d\n",myid,Num_Cells0*Ngrid2*Ngrid3);
    */

    if (SpinP_switch==3){ /* spin non-collinear */
      Density_Grid = (double**)malloc(sizeof(double*)*4); 
      for (k=0; k<=3; k++){
        if (Solver==4){
          Density_Grid[k] = (double*)malloc(sizeof(double)*Nmax); 
          for (i=0; i<Nmax; i++) Density_Grid[k][i] = 0.0;
	}
        else{
          Density_Grid[k] = (double*)malloc(sizeof(double)*N); 
          for (i=0; i<N; i++) Density_Grid[k][i] = 0.0;
        }
      }
    }
    else{
      Density_Grid = (double**)malloc(sizeof(double*)*2); 
      for (k=0; k<=1; k++){
        if (Solver==4){
          Density_Grid[k] = (double*)malloc(sizeof(double)*Nmax); 
          for (i=0; i<Nmax; i++) Density_Grid[k][i] = 0.0;
	}
        else{
          Density_Grid[k] = (double*)malloc(sizeof(double)*N); 
          for (i=0; i<N; i++) Density_Grid[k][i] = 0.0;
        }
      }
    }

    ADensity_Grid   = (double*)malloc(sizeof(double)*N); 
    PCCDensity_Grid = (double*)malloc(sizeof(double)*N); 

    if (SpinP_switch==3){ /* spin non-collinear */
      Vxc_Grid = (double**)malloc(sizeof(double*)*4); 
      for (k=0; k<=3; k++){
        Vxc_Grid[k] = (double*)malloc(sizeof(double)*N); 
      }
    }
    else{
      Vxc_Grid = (double**)malloc(sizeof(double*)*2); 
      for (k=0; k<=1; k++){
        Vxc_Grid[k] = (double*)malloc(sizeof(double)*N); 
      }
    }

    RefVxc_Grid = (double*)malloc(sizeof(double)*N); 

    VNA_Grid = (double*)malloc(sizeof(double)*N); 

    dVHart_Grid = (double*)malloc(sizeof(double)*N); 
    for (i=0; i<N; i++) dVHart_Grid[i] = 0.0;

    if (SpinP_switch==3){ /* spin non-collinear */
      Vpot_Grid = (double**)malloc(sizeof(double*)*4); 
      for (k=0; k<=3; k++){
        Vpot_Grid[k] = (double*)malloc(sizeof(double)*N); 
      }
    }
    else{
      Vpot_Grid = (double**)malloc(sizeof(double*)*2); 
      for (k=0; k<=1; k++){
        Vpot_Grid[k] = (double*)malloc(sizeof(double)*N); 
      }
    }

    /* external electric field */
    VEF_Grid = (double*)malloc(sizeof(double)*N); 

    /* Orbs_Grid */
    size_Orbs_Grid = 0;
    Orbs_Grid = (Type_Orbs_Grid***)malloc(sizeof(Type_Orbs_Grid**)*(Matomnum+MatomnumF+1)); 
    Orbs_Grid[0] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*1); 
    Orbs_Grid[0][0] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*1); 
    for (Mc_AN=1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
      Gc_AN = F_M2G[Mc_AN];
      Cwan = WhatSpecies[Gc_AN];
      Orbs_Grid[Mc_AN] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*Spe_Total_NO[Cwan]); 
      for (i=0; i<Spe_Total_NO[Cwan]; i++){
        Orbs_Grid[Mc_AN][i] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*GridN_Atom[Gc_AN]); 
        size_Orbs_Grid += GridN_Atom[Gc_AN];
      }     
    }

    /* COrbs_Grid */
    size_COrbs_Grid = 0;
    if (Cnt_switch!=0){
      COrbs_Grid = (Type_Orbs_Grid***)malloc(sizeof(Type_Orbs_Grid**)*(Matomnum+MatomnumF+1)); 
      COrbs_Grid[0] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*1); 
      COrbs_Grid[0][0] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*1); 
      for (Mc_AN=1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
        Gc_AN = F_M2G[Mc_AN];
        Cwan = WhatSpecies[Gc_AN];
        COrbs_Grid[Mc_AN] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*Spe_Total_CNO[Cwan]); 
        for (i=0; i<Spe_Total_CNO[Cwan]; i++){
          COrbs_Grid[Mc_AN][i] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*GridN_Atom[Gc_AN]); 
          size_COrbs_Grid += GridN_Atom[Gc_AN];
        }
      }
    }

    alloc_first[3] = 0;

    /* PrintMemory */

    if (firsttime){
    if (SpinP_switch==3){
    PrintMemory("truncation: Density_Grid",    sizeof(double)*N*4,              NULL);
    PrintMemory("truncation: ADensity_Grid",   sizeof(double)*N,                NULL);
    PrintMemory("truncation: PCCDensity_Grid", sizeof(double)*N,                NULL);
    PrintMemory("truncation: Vxc_Grid",        sizeof(double)*N*4,              NULL);
    PrintMemory("truncation: RefVxc_Grid",     sizeof(double)*N,                NULL);
    PrintMemory("truncation: VNA_Grid",        sizeof(double)*N,                NULL);
    PrintMemory("truncation: dVHart_Grid",     sizeof(double)*N,                NULL);
    PrintMemory("truncation: Vpot_Grid",       sizeof(double)*N*4,              NULL);
    PrintMemory("truncation: Orbs_Grid",       sizeof(Type_Orbs_Grid)*size_Orbs_Grid,   NULL);
    PrintMemory("truncation: COrbs_Grid",      sizeof(Type_Orbs_Grid)*size_COrbs_Grid,  NULL);
    }  
    else{
    PrintMemory("truncation: Density_Grid",    sizeof(double)*N*2,              NULL);
    PrintMemory("truncation: ADensity_Grid",   sizeof(double)*N,                NULL);
    PrintMemory("truncation: PCCDensity_Grid", sizeof(double)*N,                NULL);
    PrintMemory("truncation: Vxc_Grid",        sizeof(double)*N*2,              NULL);
    PrintMemory("truncation: RefVxc_Grid",     sizeof(double)*N,                NULL);
    PrintMemory("truncation: VNA_Grid",        sizeof(double)*N,                NULL);
    PrintMemory("truncation: dVHart_Grid",     sizeof(double)*N,                NULL);
    PrintMemory("truncation: Vpot_Grid",       sizeof(double)*N*2,              NULL);
    PrintMemory("truncation: Orbs_Grid",       sizeof(Type_Orbs_Grid)*size_Orbs_Grid,   NULL);
    PrintMemory("truncation: COrbs_Grid",      sizeof(Type_Orbs_Grid)*size_COrbs_Grid,  NULL);
    }

    /* external electric field */
    PrintMemory("truncation: VEF_Grid",        sizeof(double)*N,                NULL);
    }
  
  } /* if (UCell_flag==1) */

  /****************************************************
      Output the truncation data to filename.TRN
  ****************************************************/

  if (2<=level_fileout && myid==Host_ID){
    fnjoint(filepath,filename,file_TRN);
    if ((fp_TRN = fopen(file_TRN,"w")) != NULL){

#ifdef xt3
      setvbuf(fp_TRN,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      Output_Connectivity(fp_TRN);
      fclose(fp_TRN);
    }
    else
      printf("Failure of saving the TRN file.\n");

  }

  /****************************************************
    freeing of arrays:

            R_index1       (local) 
            R_index2       (local) 
            CellDis        (local)
  ****************************************************/

  Free_truncation(CpyCell,TCpyCell,1);

  /* for PrintMemory */
  firsttime = 0;
 
  /* for time */
  dtime(&TEtime);
  time0 = TEtime - TStime;
  return time0;
}







void Trn_System(int MD_iter, int CpyCell, int TCpyCell, int flag_GDC, int flag_GDC2)
{
  int i,j,k,l,m,Rn,fan,san,tan,wanA,wanB,po0;
  int ct_AN,h_AN,m2,m3,i0,size_RMI1,size_array;
  int My_TFNAN,My_TSNAN,Gh_AN,LT_switch,Nloop;
  double r,rcutA,rcutB,rcut,dx,dy,dz;
  double *fDis,*sDis;
  int *fnan2,*fncn;
  int *snan2,*sncn;
  int numprocs,myid,tag=999,ID;
  double Stime_atom, Etime_atom;
  /* for OpenMP */
  int OMPID,Nthrds,Nprocs;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /*******************************************************
                      start calc.
  *******************************************************/

#pragma omp parallel shared(ScaleSize,Max_FSNAN,level_stdout,BCR,CellDis,atv,Gxyz,Dis,ncn,natn,TCpyCell,CpyCell,atomnum,FNAN,SNAN,Spe_Atom_Cut1,WhatSpecies,M2G,Matomnum) private(OMPID,Nthrds,Nprocs,i,ct_AN,wanA,rcutA,j,wanB,rcutB,rcut,Rn,dx,dy,dz,r,l,k,fDis,fncn,fnan2,sDis,sncn,snan2,size_array)
  {

    /* allocation of arrays */
 
    size_array = (int)((Max_FSNAN+2)*ScaleSize);

    fDis = (double*)malloc(sizeof(double)*size_array); 
    sDis = (double*)malloc(sizeof(double)*size_array);

    fnan2 = (int*)malloc(sizeof(int)*size_array);
    fncn  = (int*)malloc(sizeof(int)*size_array); 
    snan2 = (int*)malloc(sizeof(int)*size_array); 
    sncn  = (int*)malloc(sizeof(int)*size_array); 

    /* get info. on OpenMP */ 
  
    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    for (i=(OMPID+1); i<=Matomnum; i+=Nthrds){

      ct_AN = M2G[i];
      wanA = WhatSpecies[ct_AN];
      rcutA = Spe_Atom_Cut1[wanA];

      FNAN[ct_AN] = 0;
      SNAN[ct_AN] = 0;

      for (j=1; j<=atomnum; j++){

	wanB = WhatSpecies[j];
	rcutB = Spe_Atom_Cut1[wanB];
	rcut = rcutA + rcutB + rcut_buffer;

	for (Rn=0; Rn<=TCpyCell; Rn++){

	  if ((ct_AN==j) && Rn==0){
	    natn[ct_AN][0] = ct_AN;
	    ncn[ct_AN][0]  = 0;
	    Dis[ct_AN][0]  = 0.0;
	  }
            
	  else{
            
	    if (((Rn==0)&&(ct_AN!=j)) || Rn!=0){
	      dx = fabs(Gxyz[ct_AN][1] - Gxyz[j][1] - atv[Rn][1]);
	      dy = fabs(Gxyz[ct_AN][2] - Gxyz[j][2] - atv[Rn][2]);
	      dz = fabs(Gxyz[ct_AN][3] - Gxyz[j][3] - atv[Rn][3]);
	      r = sqrt(dx*dx + dy*dy + dz*dz);
	    }

	    if (ct_AN==j && ct_AN==1){
	      CellDis[Rn][0] = r;
	      CellDis[Rn][1] = dx;
	      CellDis[Rn][2] = dy;
	      CellDis[Rn][3] = dz;
	    }

	    if (r<=rcut){
	      FNAN[ct_AN]++;
	      fnan2[FNAN[ct_AN]] = j;
	      fncn[FNAN[ct_AN]] = Rn;
	      fDis[FNAN[ct_AN]] = r;
	    }

	    else if (r<=BCR){
	      SNAN[ct_AN]++;
	      snan2[SNAN[ct_AN]] = j;
	      sncn[SNAN[ct_AN]] = Rn;
	      sDis[SNAN[ct_AN]] = r;
	    }

	  } /* else */

	} /* Rn */
      } /* j */

      for (k=1; k<=FNAN[ct_AN]; k++){
	natn[ct_AN][k] = fnan2[k];
	ncn[ct_AN][k] = fncn[k];
	Dis[ct_AN][k] = fDis[k];
      }
      for (k=1; k<=SNAN[ct_AN]; k++){
	l = FNAN[ct_AN] + k;
	natn[ct_AN][l] = snan2[k];
	ncn[ct_AN][l] = sncn[k];
	Dis[ct_AN][l] = sDis[k];
      }
      if (2<=level_stdout){
	printf("<physical truncation> myid=%4d CpyCell=%2d ct_AN=%4d FNAN SNAN %3d %3d\n",
	       myid,CpyCell,ct_AN,FNAN[ct_AN],SNAN[ct_AN]);fflush(stdout);
      }

    } /* i */  

    /* freeing of arrays */

    free(sncn);
    free(snan2);
    free(fncn);
    free(fnan2);
    free(sDis);
    free(fDis);

  } /* #pragma omp parallel */

  /****************************************************
   MPI: 

         FNAN
         SNAN 
         natn
         ncn
         Dis
  ****************************************************/

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){

    ID = G2ID[ct_AN];
    MPI_Bcast(&FNAN[ct_AN], 1, MPI_INT, ID, mpi_comm_level1);
    MPI_Bcast(&SNAN[ct_AN], 1, MPI_INT, ID, mpi_comm_level1);
    MPI_Bcast(&natn[ct_AN][0], FNAN[ct_AN]+SNAN[ct_AN]+1,
              MPI_INT, ID, mpi_comm_level1);
    MPI_Bcast(&ncn[ct_AN][0], FNAN[ct_AN]+SNAN[ct_AN]+1,
              MPI_INT, ID, mpi_comm_level1);
    MPI_Bcast(&Dis[ct_AN][0], FNAN[ct_AN]+SNAN[ct_AN]+1,
              MPI_DOUBLE, ID, mpi_comm_level1);
  }  

  if (myid==Host_ID && 2<=level_stdout){
    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      printf("   myid=%i ct_AN=%i  FNAN=%2d  SNAN=%2d\n",
              myid,ct_AN,FNAN[ct_AN],SNAN[ct_AN]);fflush(stdout);
    }
  }

  /****************************************************
                   logical truncation
  ****************************************************/

  if      (orderN_FNAN_SNAN_flag==0) LT_switch = 1;
  else if (orderN_FNAN_SNAN_flag==1) LT_switch = 0;

  if (LT_switch==1) LT(1);
  if (orderN_FNAN_SNAN_flag==1)  Fixed_FNAN_SNAN();

#pragma omp parallel shared(atv,Gxyz,ncn_GDC,natn_GDC,SNAN_GDC,flag_GDC2,flag_GDC,True_SNAN,Dis,ScaleSize,Max_FSNAN,LT_switch,ncn,natn,FNAN,SNAN,M2G,Matomnum) private(i,OMPID,Nthrds,Nprocs,ct_AN,m2,m3,j,l,k,Rn,size_array,fDis,sDis,fnan2,fncn,snan2,sncn,po0,Gh_AN,dx,dy,dz)
  {

    /* allocation of arrays */
 
    size_array = (int)((Max_FSNAN+2)*ScaleSize);

    fDis = (double*)malloc(sizeof(double)*size_array); 
    sDis = (double*)malloc(sizeof(double)*size_array);

    fnan2 = (int*)malloc(sizeof(int)*size_array);
    fncn  = (int*)malloc(sizeof(int)*size_array); 
    snan2 = (int*)malloc(sizeof(int)*size_array); 
    sncn  = (int*)malloc(sizeof(int)*size_array); 

    /* get info. on OpenMP */ 
  
    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    for (i=(OMPID+1); i<=Matomnum; i+=Nthrds){

      ct_AN = M2G[i];

      m2 = 0;
      m3 = 0;

      for (j=1; j<=SNAN[ct_AN]; j++){

	l = FNAN[ct_AN] + j;
	k = natn[ct_AN][l];
	Rn = ncn[ct_AN][l];

	if (LT_switch==1){  
	  if (Get_R_index(2,ct_AN,k,Rn)==1){
	    m2++;
	    snan2[m2] = k;
	    sncn[m2] = Rn;
	    sDis[m2] = Dis[ct_AN][l];
	  }
	}

	else{

	  m2++;
	  snan2[m2] = k;
	  sncn[m2]  = Rn;
	  sDis[m2]  = Dis[ct_AN][l];
	}  
      }

      SNAN[ct_AN]      = m2;
      True_SNAN[ct_AN] = m2;

      for (j=1; j<=SNAN[ct_AN]; j++){
	l = FNAN[ct_AN] + j;
	natn[ct_AN][l] = snan2[j];
	ncn[ct_AN][l] = sncn[j];
	Dis[ct_AN][l] = sDis[j];
      }

      /***********************************************************
       taking account of SNAN_GDC[ct_AN] when GDC method is used. 
      ***********************************************************/

      if (flag_GDC==1 && flag_GDC2==1){

	m2 = 0;

	for (k=0; k<SNAN_GDC[ct_AN]; k++){

	  po0 = 0;
	  for (j=0; j<=(FNAN[ct_AN]+SNAN[ct_AN]); j++){
	    if (natn[ct_AN][j]==natn_GDC[ct_AN][k] && 
		ncn[ct_AN][j]==ncn_GDC[ct_AN][k]) po0 = 1;                         
	  }

	  if (po0==0){

	    m2++;
	    l = FNAN[ct_AN] + SNAN[ct_AN] + m2;

	    Gh_AN = natn_GDC[ct_AN][k];
	    Rn = ncn_GDC[ct_AN][k];

	    natn[ct_AN][l] = Gh_AN;
	    ncn[ct_AN][l]  = Rn;

	    dx = fabs(Gxyz[ct_AN][1] - Gxyz[Gh_AN][1] - atv[Rn][1]);
	    dy = fabs(Gxyz[ct_AN][2] - Gxyz[Gh_AN][2] - atv[Rn][2]);
	    dz = fabs(Gxyz[ct_AN][3] - Gxyz[Gh_AN][3] - atv[Rn][3]);

	    Dis[ct_AN][l]  = sqrt(dx*dx + dy*dy + dz*dz);
	  } 

	}

	SNAN[ct_AN] += m2;
      }
    
    } /*  for (i=1; i<=Matomnum; i++){ */

    /* freeing of arrays */

    free(sncn);
    free(snan2);
    free(fncn);
    free(fnan2);
    free(sDis);
    free(fDis);

  } /* #pragma omp parallel */

  My_TFNAN = 0;
  My_TSNAN = 0;
  for (i=1; i<=Matomnum; i++){
    ct_AN = M2G[i];
    My_TFNAN = My_TFNAN + FNAN[ct_AN];
    My_TSNAN = My_TSNAN + SNAN[ct_AN];
  }

  /****************************************************
   MPI: 

         My_TFNAN -> TFNAN
         My_TSNAN -> TSNAN
         FNAN
         SNAN 
         natn
         ncn
         Dis
  ****************************************************/

  MPI_Reduce(&My_TFNAN, &TFNAN, 1, MPI_INT, MPI_SUM, Host_ID, mpi_comm_level1);
  MPI_Bcast(&TFNAN, 1, MPI_INT, Host_ID, mpi_comm_level1);
  MPI_Reduce(&My_TSNAN, &TSNAN, 1, MPI_INT, MPI_SUM, Host_ID, mpi_comm_level1);
  MPI_Bcast(&TSNAN, 1, MPI_INT, Host_ID, mpi_comm_level1);

  if (myid==Host_ID && 0<level_stdout){
    printf("TFNAN=%8d   Average FNAN=%10.5f\n",TFNAN,(double)TFNAN/(double)atomnum);fflush(stdout);
    printf("TSNAN=%8d   Average SNAN=%10.5f\n",TSNAN,(double)TSNAN/(double)atomnum);fflush(stdout);
  }

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    ID = G2ID[ct_AN];
    MPI_Bcast(&FNAN[ct_AN], 1, MPI_INT, ID, mpi_comm_level1);
    MPI_Bcast(&SNAN[ct_AN], 1, MPI_INT, ID, mpi_comm_level1);
    MPI_Bcast(&True_SNAN[ct_AN], 1, MPI_INT, ID, mpi_comm_level1);

    MPI_Bcast(&natn[ct_AN][0], FNAN[ct_AN]+SNAN[ct_AN]+1,
              MPI_INT, ID, mpi_comm_level1);
    MPI_Bcast(&ncn[ct_AN][0], FNAN[ct_AN]+SNAN[ct_AN]+1,
              MPI_INT, ID, mpi_comm_level1);
    MPI_Bcast(&Dis[ct_AN][0], FNAN[ct_AN]+SNAN[ct_AN]+1,
              MPI_DOUBLE, ID, mpi_comm_level1);
  }  

  if (myid==Host_ID && 0<level_stdout){
    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      printf("<logical truncation> myid=%4d CpyCell=%2d ct_AN=%4d FNAN SNAN %3d %3d\n",
              myid,CpyCell,ct_AN,FNAN[ct_AN],SNAN[ct_AN]);fflush(stdout);
    }
  }

  /****************************************************
                     Set_R_index
  ****************************************************/

  for (i=0; i<=asize1; i++){
    R_index1[i] = 0;
    R_index2[i] = 0;
  }

#pragma omp parallel shared(atomnum,FNAN,SNAN,natn,ncn) private(i,k,j,Rn,OMPID,Nthrds,Nprocs)
  {

    /* get info. on OpenMP */ 
  
    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    for (i=(OMPID+1); i<=atomnum; i+=Nthrds){
      for (k=0; k<=(FNAN[i]+SNAN[i]); k++){
	j = natn[i][k];
	Rn = ncn[i][k];
	if (k<=FNAN[i]) Set_R_index(1,i,j,Rn);
	Set_R_index(2,i,j,Rn);
      }
    }

  } /* #pragma omp parallel */

}



void Estimate_Trn_System(int CpyCell, int TCpyCell, int flag_GDC)
{
  /****************************************************
     FNAN, SNAN, Max_FNAN, Max_FSNAN are determined
               by the physical truncation
  ****************************************************/

  int i,j,ct_AN,Rn,wanA,wanB;
  double r,rcutA,rcutB,rcut;
  double dx,dy,dz;
  int numprocs,myid,tag=999,ID;
  int MFNAN,MFSNAN;
  int abnormal_bond,my_abnormal_bond;
  int spe1,spe2;

  /* for OpenMP */
  int OMPID,Nthrds,Nprocs;

  MPI_Status stat;
  MPI_Request request;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  abnormal_bond = 0;
  my_abnormal_bond = 0;

#pragma omp parallel shared(SNAN_GDC,flag_GDC,CpyCell,level_stdout,BCR,CellDis,Spe_WhatAtom,atv,Gxyz,TCpyCell,SNAN,FNAN,Spe_Atom_Cut1,WhatSpecies,M2G,Matomnum,atomnum) private(i,ct_AN,wanA,rcutA,j,wanB,rcutB,rcut,Rn,dx,dy,dz,r,spe1,spe2,OMPID,Nthrds,Nprocs)

  {

    /* get info. on OpenMP */ 
  
    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    for (i=(OMPID+1); i<=Matomnum; i+=Nthrds){
        
      ct_AN = M2G[i];
      wanA = WhatSpecies[ct_AN];
      rcutA = Spe_Atom_Cut1[wanA];

      FNAN[ct_AN] = 0;
      SNAN[ct_AN] = 0;

      for (j=1; j<=atomnum; j++){

	wanB = WhatSpecies[j];
	rcutB = Spe_Atom_Cut1[wanB];
	rcut = rcutA + rcutB + rcut_buffer;

	for (Rn=0; Rn<=TCpyCell; Rn++){

	  if ((ct_AN==j) && Rn==0){
	    /* Nothing to be done */
	  }
 
	  else{

	    if (((Rn==0)&&(ct_AN!=j)) || Rn!=0){
	      dx = fabs(Gxyz[ct_AN][1] - Gxyz[j][1] - atv[Rn][1]);
	      dy = fabs(Gxyz[ct_AN][2] - Gxyz[j][2] - atv[Rn][2]);
	      dz = fabs(Gxyz[ct_AN][3] - Gxyz[j][3] - atv[Rn][3]);
	      r = sqrt(dx*dx + dy*dy + dz*dz);

	      /* check unphysical structure */

	      spe1 = Spe_WhatAtom[wanA];
	      spe2 = Spe_WhatAtom[wanB];

	      /*
		if (r<0.5 && spe1!=0 && spe2!=0){
		my_abnormal_bond = 1;
		}
	      */

	    }

	    if (ct_AN==j && ct_AN==1){
	      CellDis[Rn][0] = r;
	      CellDis[Rn][1] = dx;
	      CellDis[Rn][2] = dy;
	      CellDis[Rn][3] = dz;
	    }

	    /*
	      if (Min_Bond_Length<r && r<=rcut)     FNAN[ct_AN] += 1;
	      else if (Min_Bond_Length<r && r<=BCR) SNAN[ct_AN] += 1;
	    */

	    if (r<=rcut)     FNAN[ct_AN]++;
	    else if (r<=BCR) SNAN[ct_AN]++;

	  }
	}
      }

      if (2<=level_stdout){
	printf("<physical truncation> CpyCell=%2d ct_AN=%2d FNAN SNAN %2d %2d\n",
	       CpyCell,i,FNAN[ct_AN],SNAN[ct_AN]);
      }

      if (flag_GDC==1) SNAN[ct_AN] += SNAN_GDC[ct_AN];
  
    } /* i */

  } /* #pragma omp parallel */

  Max_FNAN  = 0;
  Max_FSNAN = 0;

  for (i=1; i<=Matomnum; i++){
    ct_AN = M2G[i];

    if (Max_FNAN<FNAN[ct_AN])                Max_FNAN  = FNAN[ct_AN];
    if (Max_FSNAN<(FNAN[ct_AN]+SNAN[ct_AN])) Max_FSNAN = FNAN[ct_AN] + SNAN[ct_AN];
  }

  MFNAN  = Max_FNAN;
  MFSNAN = Max_FSNAN;

  /***************************************************
   MPI:
         MFNAN
         MFSNAN
  ***************************************************/ 

  /* printf("A myid=%i Max_FNAN=%i Max_FSNAN=%i\n",myid,Max_FNAN,Max_FSNAN); */

  MPI_Allreduce(&my_abnormal_bond, &abnormal_bond, 1, MPI_INT, MPI_MAX, mpi_comm_level1);

  /* check unphysical structure */

  if (abnormal_bond==1){
    if (myid==Host_ID && 0<level_stdout){
      printf("\nFound unphysical bond length: check your structure!\n\n");
    }
    MPI_Finalize();
    exit(0);
  }

  /*  MFNAN */

  MPI_Reduce(&MFNAN, &Max_FNAN, 1, MPI_INT, MPI_MAX, Host_ID, mpi_comm_level1);
  MPI_Bcast(&Max_FNAN, 1, MPI_INT, Host_ID, mpi_comm_level1);

  /*  MFSNAN */
  MPI_Reduce(&MFSNAN, &Max_FSNAN, 1, MPI_INT, MPI_MAX, Host_ID, mpi_comm_level1);
  MPI_Bcast(&Max_FSNAN, 1, MPI_INT, Host_ID, mpi_comm_level1);

  List_YOUSO[8] = Max_FNAN  + 7;
  if (Solver==1 || Solver==5 || Solver==6 || Solver==7 || Solver==8)
    List_YOUSO[2] = Max_FSNAN + 7;
  else
    List_YOUSO[2] = Max_FNAN  + 7;
}




void Check_System()
{

  char *s_vec[20];
  int i,ct_AN,h_AN,Rn,po[4],num;
  int myid;

  MPI_Comm_rank(mpi_comm_level1,&myid);

  po[1] = 0;
  po[2] = 0;
  po[3] = 0;

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    for (h_AN=1; h_AN<=FNAN[ct_AN]; h_AN++){
      Rn = ncn[ct_AN][h_AN];
      if (Rn!=0){
        if (atv_ijk[Rn][1]!=0) po[1] = 1;
        if (atv_ijk[Rn][2]!=0) po[2] = 1;
        if (atv_ijk[Rn][3]!=0) po[3] = 1;
      }
    }
  }

  num = 0;  
  for (i=1; i<=3; i++){
    if (po[i]==1) num++;
  }

  s_vec[0] = "molecule.";
  s_vec[1] = "chain.";
  s_vec[2] = "slab.";
  s_vec[3] = "bulk.";

  specified_system = num;
  
  if (myid==Host_ID && 0<level_stdout) printf("<Check_System> The system is %s\n",s_vec[num]);

  /* in case of solver==cluster */
  if (specified_system!=0 && Solver==2) PeriodicGamma_flag = 1;
}




void Set_RMI()
{
  /****************************************************

    What is RMI?

    RMI[Mc_AN][i][j] is a array which specifies
    the position of arraies storing hopping and
    overlap integrals between atoms i and j.

    If ig = natn[Gc_AN][i]
       mi = F_G2M[ig] or S_G2M[ig]
       k = RMI[Mc_AN][i][j],

    then, we can find the hopping and overlap integrals
    between atoms i and j from h[mi][k] and OLP[mi][k],
    respectively.
  ****************************************************/

  static int firsttime=1;
  int i,j,k,Mc_AN,Gc_AN,size_RMI1;
  int fan,san,can,wan,ig,Rni,jg,Rnj;
  int l1,l2,l3,m1,m2,m3,Rn;
  int i_rlt,po;
  int numprocs,myid,tag=999,ID;
  double Stime_atom, Etime_atom;
  /* for OpenMP */
  int OMPID,Nthrds,Nprocs;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
   allocation of arrays:

      RMI1[Matomnum+1]
          [FNAN[Gc_AN]+SNAN[Gc_AN]+1]
          [FNAN[Gc_AN]+SNAN[Gc_AN]+1]

      RMI2[Matomnum+1]
          [FNAN[Gc_AN]+SNAN[Gc_AN]+1]
          [FNAN[Gc_AN]+SNAN[Gc_AN]+1]
  ****************************************************/
  
  FNAN[0] = 0;
  SNAN[0] = 0;

  size_RMI1 = 0;
  RMI1 = (int***)malloc(sizeof(int**)*(Matomnum+1)); 
  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
    if (Mc_AN==0) Gc_AN = 0;
    else          Gc_AN = M2G[Mc_AN];
    RMI1[Mc_AN] = (int**)malloc(sizeof(int*)*(FNAN[Gc_AN]+SNAN[Gc_AN]+1)); 
    for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
      RMI1[Mc_AN][i] = (int*)malloc(sizeof(int)*(FNAN[Gc_AN]+SNAN[Gc_AN]+1)); 
      size_RMI1 += FNAN[Gc_AN]+SNAN[Gc_AN]+1;
    }      
  }
  
  RMI2 = (int***)malloc(sizeof(int**)*(Matomnum+1)); 
  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
    if (Mc_AN==0) Gc_AN = 0;
    else          Gc_AN = M2G[Mc_AN];
    RMI2[Mc_AN] = (int**)malloc(sizeof(int*)*(FNAN[Gc_AN]+SNAN[Gc_AN]+1)); 
    for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
      RMI2[Mc_AN][i] = (int*)malloc(sizeof(int)*(FNAN[Gc_AN]+SNAN[Gc_AN]+1)); 
    }      
  }

  alloc_first[6] = 0;
  
  if (firsttime){
  PrintMemory("truncation: RMI1", sizeof(int)*size_RMI1, NULL);
  PrintMemory("truncation: RMI2", sizeof(int)*size_RMI1, NULL);
  firsttime = 0; 
  }

  /****************************************************
                  setting of RMI matrix
  ****************************************************/

#pragma omp parallel shared(time_per_atom,RMI2,RMI1,ratv,CpyCell,atv_ijk,ncn,natn,WhatSpecies,SNAN,FNAN,M2G) private(OMPID,Nthrds,Nprocs,Mc_AN,Stime_atom,Gc_AN,fan,san,can,wan,i,ig,Rni,j,jg,Rnj,l1,l2,l3,m1,m2,m3,Rn,i_rlt,k,po,Etime_atom)
  {

    /* get info. on OpenMP */ 
  
    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();
  
    for (Mc_AN=(OMPID+1); Mc_AN<=Matomnum; Mc_AN+=Nthrds){

      dtime(&Stime_atom);

      Gc_AN = M2G[Mc_AN];
      fan = FNAN[Gc_AN];
      san = SNAN[Gc_AN];
      can = fan + san;
      wan = WhatSpecies[Gc_AN];

      for (i=0; i<=can; i++){

	ig = natn[Gc_AN][i];
	Rni = ncn[Gc_AN][i];

	for (j=0; j<=can; j++){

	  jg = natn[Gc_AN][j];
	  Rnj = ncn[Gc_AN][j];
	  l1 = atv_ijk[Rnj][1] - atv_ijk[Rni][1];
	  l2 = atv_ijk[Rnj][2] - atv_ijk[Rni][2];
	  l3 = atv_ijk[Rnj][3] - atv_ijk[Rni][3];
	  if (l1<0) m1=-l1; else m1 = l1;
	  if (l2<0) m2=-l2; else m2 = l2;
	  if (l3<0) m3=-l3; else m3 = l3;

	  if (m1<=CpyCell && m2<=CpyCell && m3<=CpyCell){

	    Rn = ratv[l1+CpyCell][l2+CpyCell][l3+CpyCell];

	    /* FNAN */

	    i_rlt = Get_R_index(1,ig,jg,Rn);

	    if (i_rlt==1){
	      k = -1;
	      po = 0;
	      while(po==0){
		k++;
		if (natn[ig][k]==jg && ncn[ig][k]==Rn){
		  RMI1[Mc_AN][i][j] = k;
		  po = 1;
		}
	      }
	    }
	    else{
	      RMI1[Mc_AN][i][j] = -1;
	    }

	    /* FNAN + SNAN */

	    i_rlt = Get_R_index(2,ig,jg,Rn);

	    if (i_rlt==1){
	      k = -1;
	      po = 0;
	      while(po==0){
		k++;
		if (natn[ig][k]==jg && ncn[ig][k]==Rn){
		  RMI2[Mc_AN][i][j] = k;
		  po = 1;
		}
	      }
	    }
	    else{
	      RMI2[Mc_AN][i][j] = -1;
	    }

	  }
	  else{
	    RMI1[Mc_AN][i][j] = -1;
	    RMI2[Mc_AN][i][j] = -1;
	  }
	}
      }

      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;

    } /* Mc_AN */
  } /* #pragma omp parallel */


}




void Fixed_FNAN_SNAN()
{
  int i,j,Gc_AN;
  int hops;
  int numprocs,myid;
  double Stime_atom, Etime_atom;
  /* for OpenMP */
  int OMPID,Nthrds,Nprocs;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  for (i=1; i<=Matomnum; i++){

    dtime(&Stime_atom);

    Gc_AN = M2G[i];

    if ( orderN_FNAN_SNAN[Gc_AN]<=(FNAN[Gc_AN]+SNAN[Gc_AN]) ){

      qsort_double3B(SNAN[Gc_AN], &Dis[Gc_AN][FNAN[Gc_AN]+1], 
                                  &natn[Gc_AN][FNAN[Gc_AN]+1], 
                                  &ncn[Gc_AN][FNAN[Gc_AN]+1]);

      SNAN[Gc_AN] = orderN_FNAN_SNAN[Gc_AN] - FNAN[Gc_AN];
      if (SNAN[Gc_AN]<0) SNAN[Gc_AN] = 0;
    }

    dtime(&Etime_atom);
    time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
  }

}



void LT(int CorL)
{
  int i,Gc_AN;
  int hops;
  int numprocs,myid;
  double Stime_atom, Etime_atom;
  /* for OpenMP */
  int OMPID,Nthrds,Nprocs;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  for (i=0; i<=asize1; i++){
    R_index2[i] = 0;
  }

#pragma omp parallel shared(atomnum) private(i,OMPID,Nthrds,Nprocs)
 {
  /* get info. on OpenMP */ 
  
  OMPID = omp_get_thread_num();
  Nthrds = omp_get_num_threads();
  Nprocs = omp_get_num_procs();

  for (i=(OMPID+1); i<=atomnum; i+=Nthrds) Set_R_index(2,i,i,0);

 } /* #pragma omp parallel */

#pragma omp parallel shared(time_per_atom,NOHS_L,CorL,M2G,Matomnum) private(i,OMPID,Nthrds,Nprocs,Stime_atom,Etime_atom,Gc_AN,hops)
  {
    /* get info. on OpenMP */ 
  
    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    for (i=(OMPID+1); i<=Matomnum; i+=Nthrds){

      dtime(&Stime_atom);

      Gc_AN = M2G[i];
      if (CorL==1) hops = NOHS_L;
      if (CorL==2) hops = NOHS_C;

      Reflexive(Gc_AN,Gc_AN,0,0,0,1,hops);

      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
    }

  } /* #pragma omp parallel */

}





void Reflexive(int Iatom,int Catom,int a,int b,int c,int Chops,int Fhops)
{
  int i,gi,l1,l2,l3;
  int m1,m2,m3,wanA,wanB;
  int Rn,IRn,GH_AN;
  double r,rcut,rcutA,rcutB;

  wanA = WhatSpecies[Catom];
  rcutA = Spe_Atom_Cut1[wanA];
  for (i=1; i<=FNAN[Catom]; i++){
    r = Dis[Catom][i];
    GH_AN = natn[Catom][i];
    wanB = WhatSpecies[GH_AN];
    rcutB = Spe_Atom_Cut1[wanB];
    rcut = rcutA + rcutB;

    if (0.1<r && r<rcut){
      gi = natn[Catom][i];
      Rn = ncn[Catom][i];
      l1 = atv_ijk[Rn][1] + a;
      l2 = atv_ijk[Rn][2] + b;
      l3 = atv_ijk[Rn][3] + c;
      if (l1<0) m1=-l1; else m1 = l1;
      if (l2<0) m2=-l2; else m2 = l2;
      if (l3<0) m3=-l3; else m3 = l3;
      if (m1<=CpyCell && m2<=CpyCell && m3<=CpyCell){
	IRn = ratv[l1+CpyCell][l2+CpyCell][l3+CpyCell];
	if (gi!=Iatom || IRn!=0){
	  if (Get_R_index(2,Iatom,gi,IRn)==0) Set_R_index(2,Iatom,gi,IRn);
	  if ((Chops+1)<=Fhops) Reflexive(Iatom,gi,l1,l2,l3,Chops+1,Fhops);
	}
      }
    }
  }
}





void Set_R_index(int R_switch,int i,int j,int Rn)
{
  int odn,odni,odnj;

  /* Warning for mixing int and long types */
  odn = (i-1)*atomnum*(TCpyCell+1) + (j-1)*(TCpyCell+1) + Rn;
  odnj = odn % 32;
  odni = (odn-odnj)/32;
  if (R_switch==1) R_index1[odni] = R_index1[odni]+(1<<odnj);
  if (R_switch==2) R_index2[odni] = R_index2[odni]+(1<<odnj);
}




int Get_R_index(int R_switch,int i,int j,int Rn)
{
  long ltmp;
  int odn,odni,odnj,itmp;
  int result,b;

  /* Warning for mixing int and long types */
  odn = (i-1)*atomnum*(TCpyCell+1) + (j-1)*(TCpyCell+1) + Rn;
  odnj = odn % 32;
  odni = (odn-odnj)/32;

  if (R_switch==1){

    if (sizeof(long)==8){
      itmp = R_index1[odni];
      b = itmp<<(31-odnj)&1<<BYTESIZE*sizeof(int)-1;
    }  
    else {
      b = R_index1[odni]<<(31-odnj)&1<<BYTESIZE*sizeof(int)-1;
    }

    if (b==0)
      result = 0;
    else 
      result = 1;
  }
  else if (R_switch==2){

    if (sizeof(long)==8){
      itmp = R_index2[odni];
      b = itmp<<(31-odnj)&1<<BYTESIZE*sizeof(int)-1;
    }
    else{
      b = R_index2[odni]<<(31-odnj)&1<<BYTESIZE*sizeof(int)-1;
    }

    if (b==0)
      result = 0;
    else 
      result = 1;
  }
  return result;
}



void Output_Connectivity(FILE *fp)
{
  int ct_AN,h_AN,i,j,can;

  fprintf(fp,"\n");
  fprintf(fp,"***********************************************************\n");
  fprintf(fp,"***********************************************************\n");
  fprintf(fp,"                       Connectivity                        \n");
  fprintf(fp,"***********************************************************\n");
  fprintf(fp,"***********************************************************\n");

  fprintf(fp,"   FNAN SNAN\n");
  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    fprintf(fp,"      %i   ",ct_AN);
    fprintf(fp,"%i %i\n",FNAN[ct_AN],SNAN[ct_AN]);
  }

  fprintf(fp,"   natn\n");
  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    fprintf(fp,"      %i   ",ct_AN);
    for (h_AN=0; h_AN<=(FNAN[ct_AN]+SNAN[ct_AN]); h_AN++){
      fprintf(fp,"%i ",natn[ct_AN][h_AN]);
    }
    fprintf(fp,"\n");
  }

  fprintf(fp,"   ncn\n");
  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    fprintf(fp,"      %i   ",ct_AN);
    for (h_AN=0; h_AN<=(FNAN[ct_AN]+SNAN[ct_AN]); h_AN++){
      fprintf(fp,"%i ",ncn[ct_AN][h_AN]);
    }
    fprintf(fp,"\n");
  }

  fprintf(fp,"   Dis\n");
  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    fprintf(fp,"      %i   ",ct_AN);
    for (h_AN=0; h_AN<=(FNAN[ct_AN]+SNAN[ct_AN]); h_AN++){
      fprintf(fp,"%7.4f ",0.529177*Dis[ct_AN][h_AN]);
    }
    fprintf(fp,"\n");
  }

  /*
  fprintf(fp,"   RMI1\n");
  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    fprintf(fp,"      %i\n",ct_AN);
    can = FNAN[ct_AN] + SNAN[ct_AN];
    for (i=0; i<=can; i++){
      for (j=0; j<=can; j++){
        if (j==0)
          fprintf(fp,"      %i ",RMI1[ct_AN][i][j]);
        else 
          fprintf(fp,"%i ",RMI1[ct_AN][i][j]);
      }
      fprintf(fp,"\n");
    }
  }

  fprintf(fp,"   RMI2\n");
  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    fprintf(fp,"      %i\n",ct_AN);
    can = FNAN[ct_AN] + SNAN[ct_AN];
    for (i=0; i<=can; i++){
      for (j=0; j<=can; j++){
        if (j==0)
          fprintf(fp,"      %i ",RMI2[ct_AN][i][j]);
        else 
          fprintf(fp,"%i ",RMI2[ct_AN][i][j]);
      }
      fprintf(fp,"\n");
    }
  }
  */

}



  

void UCell_Box(int MD_iter, int estimate_switch, int CpyCell)
{
  static int firsttime=1;
  int size_GListTAtoms0;
  int po,N3[4],i,j;
  int size_GridListAtom,size_MGridListAtom;
  int NOC[4],nn1,nn2,nn3,l1,l2,l3,lmax;
  int ct_AN,N,p,n1,n2,n3,Cwan,Rn,Rn1;
  int Nct,MinN,Scale,ScaleA,ScaleB,ScaleC,Nm[4];
  int Mc_AN,Mc_AN0,Gc_AN,h_AN,h_AN0,Mh_AN,Gh_AN,Nh,Rh,Nog,GNh,GRh;
  int ll1,ll2,ll3,Nnb,My_Max_NumOLG;
  int lll1,lll2,lll3,GRh1,Nc,GNc,GRc,size_array;
  int *TAtoms0,*TCells0,*TAtoms1,*TAtoms2;
  int **Tmp_GridListAtom,**Tmp_CellListAtom;
  double MinR,CutR2,r2,rws,tmp0;
  double sa,sa_cri,tmp[4],Cxyz[4];
  double xc,yc,zc,xm,ym,zm;
  double dx,dy,dz,r,sn1,sn2,sn3;
  double A2,B2,C2,CellV;
  double S_Lng,L_Lng,LngA,LngB,LngC,x,y,z;
  double GVolume,buffer_scale,GridV;
  int *TempGrid,*TempCell;
  int numprocs,myid,ID,IDS,IDR,tag=999;
  double Stime_atom, Etime_atom;
  char file_UCell[YOUSO10];
  FILE *fp;
  char buf[fp_bsize];          /* setvbuf */
  /* for OpenMP */
  int OMPID,Nthrds,Nprocs;

  MPI_Status stat;
  MPI_Request request;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
                Reciprocal lattice vectors
  ****************************************************/

  if (estimate_switch<=1){
  
    Cross_Product(tv[2],tv[3],tmp);
    CellV = Dot_Product(tv[1],tmp); 
    Cell_Volume = fabs(CellV);
  
    Cross_Product(tv[2],tv[3],tmp);
    rtv[1][1] = 2.0*PI*tmp[1]/CellV;
    rtv[1][2] = 2.0*PI*tmp[2]/CellV;
    rtv[1][3] = 2.0*PI*tmp[3]/CellV;
  
    Cross_Product(tv[3],tv[1],tmp);
    rtv[2][1] = 2.0*PI*tmp[1]/CellV;
    rtv[2][2] = 2.0*PI*tmp[2]/CellV;
    rtv[2][3] = 2.0*PI*tmp[3]/CellV;
  
    Cross_Product(tv[1],tv[2],tmp);
    rtv[3][1] = 2.0*PI*tmp[1]/CellV;
    rtv[3][2] = 2.0*PI*tmp[2]/CellV;
    rtv[3][3] = 2.0*PI*tmp[3]/CellV;

    if (myid==Host_ID && 0<level_stdout){    

      printf("lattice vectors (bohr)\n");
      printf("A  = %15.12f, %15.12f, %15.12f\n",tv[1][1],tv[1][2],tv[1][3]);
      printf("B  = %15.12f, %15.12f, %15.12f\n",tv[2][1],tv[2][2],tv[2][3]);
      printf("C  = %15.12f, %15.12f, %15.12f\n",tv[3][1],tv[3][2],tv[3][3]);

      printf("reciprocal lattice vectors (bohr^-1)\n");
      printf("RA = %15.12f, %15.12f, %15.12f\n",rtv[1][1],rtv[1][2],rtv[1][3]);
      printf("RB = %15.12f, %15.12f, %15.12f\n",rtv[2][1],rtv[2][2],rtv[2][3]);
      printf("RB = %15.12f, %15.12f, %15.12f\n",rtv[3][1],rtv[3][2],rtv[3][3]);

    }  

    if (Ngrid_fixed_flag==0){

      /****************************************************
                        set real space grids
      ****************************************************/
  
      LngA = sqrt( Dot_Product(tv[1],tv[1]) );
      LngB = sqrt( Dot_Product(tv[2],tv[2]) );
      LngC = sqrt( Dot_Product(tv[3],tv[3]) );
  
      S_Lng = smallest(LngA, LngB);
      S_Lng = smallest(S_Lng, LngC);
      L_Lng = largest(LngA, LngB);
      L_Lng = largest(L_Lng, LngC);

      if      ((L_Lng/S_Lng)<1.1) Scale =16;
      else if ((L_Lng/S_Lng)<1.3) Scale = 8;
      else if ((L_Lng/S_Lng)<2.0) Scale = 4;
      else if ((L_Lng/S_Lng)<3.0) Scale = 2;
      else                        Scale = 1;

      ScaleA = (int)(Scale*LngA/S_Lng + 0.4999);
      ScaleB = (int)(Scale*LngB/S_Lng + 0.4999);
      ScaleC = (int)(Scale*LngC/S_Lng + 0.4999);

      l1 = 0;
      po = 0;
  
      do {
	l1++;
      
	gtv[1][1] = tv[1][1]/((double)(ScaleA*l1));
	gtv[1][2] = tv[1][2]/((double)(ScaleA*l1));
	gtv[1][3] = tv[1][3]/((double)(ScaleA*l1));
    
	gtv[2][1] = tv[2][1]/((double)(ScaleB*l1));
	gtv[2][2] = tv[2][2]/((double)(ScaleB*l1));
	gtv[2][3] = tv[2][3]/((double)(ScaleB*l1));
    
	gtv[3][1] = tv[3][1]/((double)(ScaleC*l1));
	gtv[3][2] = tv[3][2]/((double)(ScaleC*l1));
	gtv[3][3] = tv[3][3]/((double)(ScaleC*l1));
    
	Cross_Product(gtv[2],gtv[3],tmp);

	GridV = Dot_Product(gtv[1],tmp); 
	GVolume = fabs( GridV );

	Cross_Product(gtv[2],gtv[3],tmp);
	rgtv[1][1] = 2.0*PI*tmp[1]/GridV;
	rgtv[1][2] = 2.0*PI*tmp[2]/GridV;
	rgtv[1][3] = 2.0*PI*tmp[3]/GridV;

	Cross_Product(gtv[3],gtv[1],tmp);
	rgtv[2][1] = 2.0*PI*tmp[1]/GridV;
	rgtv[2][2] = 2.0*PI*tmp[2]/GridV;
	rgtv[2][3] = 2.0*PI*tmp[3]/GridV;
    
	Cross_Product(gtv[1],gtv[2],tmp);
	rgtv[3][1] = 2.0*PI*tmp[1]/GridV;
	rgtv[3][2] = 2.0*PI*tmp[2]/GridV;
	rgtv[3][3] = 2.0*PI*tmp[3]/GridV;

	A2 = rgtv[1][1]*rgtv[1][1] + rgtv[1][2]*rgtv[1][2] + rgtv[1][3]*rgtv[1][3];
	B2 = rgtv[2][1]*rgtv[2][1] + rgtv[2][2]*rgtv[2][2] + rgtv[2][3]*rgtv[2][3];
	C2 = rgtv[3][1]*rgtv[3][1] + rgtv[3][2]*rgtv[3][2] + rgtv[3][3]*rgtv[3][3];

	A2 = A2/4.0;
	B2 = B2/4.0;
	C2 = C2/4.0;

	if (estimate_switch==0 || 2<=level_stdout){ 
	  if (myid==Host_ID && 0<level_stdout){ 
	    printf("Trial cutoff energies (a,b,c) = %7.3f (%i), %7.3f (%i), %7.3f (%i)\n",
		   A2,l1*ScaleA,B2,l1*ScaleB,C2,l1*ScaleC);
	  }
	}

	if (Grid_Ecut<=A2 && Grid_Ecut<=B2 && Grid_Ecut<=C2){
	  po = 1;
	  lmax = l1;
	}

      }while (po==0);
   
      Ngrid1 = lmax*ScaleA;
      Ngrid2 = lmax*ScaleB;
      Ngrid3 = lmax*ScaleC;

      if (myid==Host_ID && 0<level_stdout){    
	printf("UCell_Box: Cutoff=%lf(%i) %lf(%i) %lf(%i)\n",
	       A2,Ngrid1,B2,Ngrid2,C2,Ngrid3);
      }


      /* min. and max. is found.  max=Ngrid?, min=Ngrid/2 */

      Find_ApproxFactN(tv,&Grid_Ecut,&Ngrid1,&Ngrid2,&Ngrid3,&A2,&B2,&C2);

      if (Solver==4) {
        TRAN_adjust_Ngrid(mpi_comm_level1, &Ngrid1, &Ngrid2, &Ngrid3);
      }

      /* recalculate gtv, rgtv, A2, B2, C2 */

      gtv[1][1] = tv[1][1]/(double)Ngrid1;
      gtv[1][2] = tv[1][2]/(double)Ngrid1;
      gtv[1][3] = tv[1][3]/(double)Ngrid1;
    
      gtv[2][1] = tv[2][1]/(double)Ngrid2;
      gtv[2][2] = tv[2][2]/(double)Ngrid2;
      gtv[2][3] = tv[2][3]/(double)Ngrid2;
    
      gtv[3][1] = tv[3][1]/(double)Ngrid3;
      gtv[3][2] = tv[3][2]/(double)Ngrid3;
      gtv[3][3] = tv[3][3]/(double)Ngrid3;

      Cross_Product(gtv[2],gtv[3],tmp);

      GridV = Dot_Product(gtv[1],tmp); 
      GVolume = fabs( GridV );

      Cross_Product(gtv[2],gtv[3],tmp);
      rgtv[1][1] = 2.0*PI*tmp[1]/GridV;
      rgtv[1][2] = 2.0*PI*tmp[2]/GridV;
      rgtv[1][3] = 2.0*PI*tmp[3]/GridV;

      Cross_Product(gtv[3],gtv[1],tmp);
      rgtv[2][1] = 2.0*PI*tmp[1]/GridV;
      rgtv[2][2] = 2.0*PI*tmp[2]/GridV;
      rgtv[2][3] = 2.0*PI*tmp[3]/GridV;
    
      Cross_Product(gtv[1],gtv[2],tmp);
      rgtv[3][1] = 2.0*PI*tmp[1]/GridV;
      rgtv[3][2] = 2.0*PI*tmp[2]/GridV;
      rgtv[3][3] = 2.0*PI*tmp[3]/GridV;

      A2 = rgtv[1][1]*rgtv[1][1] + rgtv[1][2]*rgtv[1][2] + rgtv[1][3]*rgtv[1][3];
      B2 = rgtv[2][1]*rgtv[2][1] + rgtv[2][2]*rgtv[2][2] + rgtv[2][3]*rgtv[2][3];
      C2 = rgtv[3][1]*rgtv[3][1] + rgtv[3][2]*rgtv[3][2] + rgtv[3][3]*rgtv[3][3];

      A2 = A2/4.0;
      B2 = B2/4.0;
      C2 = C2/4.0;

      if (myid==Host_ID && 0<level_stdout){    
        printf("UCell_Box: (tuned) Cutoff=%lf(%i) %lf(%i) %lf(%i)\n",
                  A2,Ngrid1,B2,Ngrid2,C2,Ngrid3);
      }   

      /* for mixed basis set */

      if (Mixed_Basis_flag==1){
        
        if ( ((Ngrid1/Ngrid1_FE)!=NG_Mixed_Basis)
          || ((Ngrid2/Ngrid2_FE)!=NG_Mixed_Basis)
	  || ((Ngrid3/Ngrid3_FE)!=NG_Mixed_Basis)
          || ((Ngrid1%Ngrid1_FE)!=0)    
          || ((Ngrid1%Ngrid2_FE)!=0)    
          || ((Ngrid1%Ngrid3_FE)!=0) ){

          if (myid==Host_ID && 0<level_stdout){
            printf("Failed setting up grids for the mixed basis scheme.\n");
          }
          MPI_Finalize();
          exit(0);
        }    
      } 

    } /* if (Ngrid_fixed_flag==0) */

    else{

      gtv[1][1] = tv[1][1]/(double)Ngrid1;
      gtv[1][2] = tv[1][2]/(double)Ngrid1;
      gtv[1][3] = tv[1][3]/(double)Ngrid1;
    
      gtv[2][1] = tv[2][1]/(double)Ngrid2;
      gtv[2][2] = tv[2][2]/(double)Ngrid2;
      gtv[2][3] = tv[2][3]/(double)Ngrid2;
    
      gtv[3][1] = tv[3][1]/(double)Ngrid3;
      gtv[3][2] = tv[3][2]/(double)Ngrid3;
      gtv[3][3] = tv[3][3]/(double)Ngrid3;
       
      Cross_Product(gtv[2],gtv[3],tmp);

      GridV = Dot_Product(gtv[1],tmp); 
      GVolume = fabs( GridV );

      Cross_Product(gtv[2],gtv[3],tmp);
      rgtv[1][1] = 2.0*PI*tmp[1]/GridV;
      rgtv[1][2] = 2.0*PI*tmp[2]/GridV;
      rgtv[1][3] = 2.0*PI*tmp[3]/GridV;

      Cross_Product(gtv[3],gtv[1],tmp);
      rgtv[2][1] = 2.0*PI*tmp[1]/GridV;
      rgtv[2][2] = 2.0*PI*tmp[2]/GridV;
      rgtv[2][3] = 2.0*PI*tmp[3]/GridV;
    
      Cross_Product(gtv[1],gtv[2],tmp);
      rgtv[3][1] = 2.0*PI*tmp[1]/GridV;
      rgtv[3][2] = 2.0*PI*tmp[2]/GridV;
      rgtv[3][3] = 2.0*PI*tmp[3]/GridV;

      A2 = rgtv[1][1]*rgtv[1][1] + rgtv[1][2]*rgtv[1][2] + rgtv[1][3]*rgtv[1][3];
      B2 = rgtv[2][1]*rgtv[2][1] + rgtv[2][2]*rgtv[2][2] + rgtv[2][3]*rgtv[2][3];
      C2 = rgtv[3][1]*rgtv[3][1] + rgtv[3][2]*rgtv[3][2] + rgtv[3][3]*rgtv[3][3];

      A2 = A2/4.0;
      B2 = B2/4.0;
      C2 = C2/4.0;

      Grid_Ecut = (A2 + B2 + C2)/3.0;
    }

    if (estimate_switch==0 || 2<=level_stdout){
      if (myid==Host_ID && 0<level_stdout) {
        printf("Required cutoff energy (Ryd) for 3D-grids = %7.4f\n",Grid_Ecut);
        printf("    Used cutoff energy (Ryd) for 3D-grids = %7.4f, %7.4f, %7.4f\n",
                A2,B2,C2);
        printf("Num. of grids of a-, b-, and c-axes = %2d, %2d, %2d\n",
                Ngrid1,Ngrid2,Ngrid3);

        /***********************************************
            output informations on grids to a file
        ***********************************************/

        sprintf(file_UCell,"%s%s.UCell",filepath,filename);
        fp = fopen(file_UCell, "w");

#ifdef xt3
        setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

        if (fp!=NULL && estimate_switch==0){
          fprintf(fp,"\n\n***********************************************************\n");
          fprintf(fp,"***********************************************************\n\n");
          fprintf(fp,"  Required cutoff energy (Ryd) for 3D-grids = %7.4f\n",Grid_Ecut);
          fprintf(fp,"      Used cutoff energy (Ryd) for 3D-grids = %7.4f, %7.4f, %7.4f\n",
                  A2,B2,C2);
          fprintf(fp,"  Num. of grids of a-, b-, and c-axes = %2d, %2d, %2d\n\n",
                  Ngrid1,Ngrid2,Ngrid3);

          fprintf(fp,"  Num.Grid1. %5d\n",Ngrid1);
          fprintf(fp,"  Num.Grid2. %5d\n",Ngrid2);
          fprintf(fp,"  Num.Grid3. %5d\n",Ngrid3);
          fprintf(fp,"\n\n");

          fclose(fp); 
        }
      }
    }

    Max_OneD_Grids = Ngrid1;
    if (Max_OneD_Grids<Ngrid2) Max_OneD_Grids = Ngrid2;
    if (Max_OneD_Grids<Ngrid3) Max_OneD_Grids = Ngrid3;

  } /* if (estimate_switch<=1) */

  /****************************************************
      Setting of the center of unit cell and grids
  ****************************************************/

  /* the center of the system */

  xc = 0.0;
  yc = 0.0;
  zc = 0.0;

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    xc += Gxyz[ct_AN][1];
    yc += Gxyz[ct_AN][2];
    zc += Gxyz[ct_AN][3];
  }

  xc = xc/(double)atomnum;
  yc = yc/(double)atomnum;
  zc = zc/(double)atomnum;

  /* added by T.Ohwaki */
  if(length_gtv[1] != 0.0 && ESM_switch!=0){

    double xc_tmp; /* added by T.Ohwaki */

    xc_tmp = (int)(xc / length_gtv[1]);
    xc = ((double)xc_tmp * length_gtv[1]);
    if(myid==Host_ID && 0<level_stdout){
      printf("<ESM> length_gtv[1],xc/length_gtv[1] = %12.9f,%12.9f \n",length_gtv[1],xc/length_gtv[1]);
    }
  } /* added by T.Ohwaki */

  dx = xc - LastBoxCenterX;
  dy = yc - LastBoxCenterY;
  dz = zc - LastBoxCenterZ;
  r = sqrt(dx*dx + dy*dy + dz*dz);
  tmp0 = 3.0*Cell_Volume/(4.0*PI);
  rws = pow(tmp0,0.3333333333333333333);

  if (MD_iter==1 || Last_TNumGrid<Ngrid1*Ngrid2*Ngrid3){

    /************
     start calc.
    ************/

    /* gtv */

    gtv[1][1] = tv[1][1]/(double)Ngrid1;
    gtv[1][2] = tv[1][2]/(double)Ngrid1;
    gtv[1][3] = tv[1][3]/(double)Ngrid1;

    gtv[2][1] = tv[2][1]/(double)Ngrid2;
    gtv[2][2] = tv[2][2]/(double)Ngrid2;
    gtv[2][3] = tv[2][3]/(double)Ngrid2;

    gtv[3][1] = tv[3][1]/(double)Ngrid3;
    gtv[3][2] = tv[3][2]/(double)Ngrid3;
    gtv[3][3] = tv[3][3]/(double)Ngrid3;

    sn1 = 0.5*( (Ngrid1+1) % 2 ); 
    sn2 = 0.5*( (Ngrid2+1) % 2 ); 
    sn3 = 0.5*( (Ngrid3+1) % 2 ); 

    xm = ( (double)(Ngrid1/2) - sn1 )*gtv[1][1]
       + ( (double)(Ngrid2/2) - sn2 )*gtv[2][1]
       + ( (double)(Ngrid3/2) - sn3 )*gtv[3][1];

    ym = ( (double)(Ngrid1/2) - sn1 )*gtv[1][2]
       + ( (double)(Ngrid2/2) - sn2 )*gtv[2][2]
       + ( (double)(Ngrid3/2) - sn3 )*gtv[3][2];

    zm = ( (double)(Ngrid1/2) - sn1 )*gtv[1][3]
       + ( (double)(Ngrid2/2) - sn2 )*gtv[2][3]
       + ( (double)(Ngrid3/2) - sn3 )*gtv[3][3];

    if ( 1.0e+8<scf_fixed_origin[0] &&
         1.0e+8<scf_fixed_origin[1] &&
	 1.0e+8<scf_fixed_origin[2] ){

      Grid_Origin[1] = xc - xm;
      Grid_Origin[2] = yc - ym;
      Grid_Origin[3] = zc - zm;

      /* added by T.Ohwaki */
      if (myid==Host_ID && ESM_switch!=0 && 0<level_stdout){
        printf("xc=%15.12f yc=%15.12f zc=%15.12f\n",xc,yc,zc);
        printf("xm=%15.12f ym=%15.12f zm=%15.12f\n",xm,ym,zm);
      }

    }
    else{
      Grid_Origin[1] = scf_fixed_origin[0];
      Grid_Origin[2] = scf_fixed_origin[1];
      Grid_Origin[3] = scf_fixed_origin[2];
    }
    
    if (myid==Host_ID && 0<level_stdout){
      printf("Grid_Origin %15.12f %15.12f %15.12f\n",
              Grid_Origin[1],Grid_Origin[2],Grid_Origin[3]);
    }    

    TNumGrid = Ngrid1*Ngrid2*Ngrid3;
    Last_TNumGrid = (int)(ScaleSize*Ngrid1*Ngrid2*Ngrid3);

    LastBoxCenterX = xc;
    LastBoxCenterY = yc;
    LastBoxCenterZ = zc;

  }

  /****************************************************
            xyz-coordinate to cell-coordinate
  ****************************************************/

  if (estimate_switch<=1){

#pragma omp parallel shared(level_stdout,rtv,Cell_Gxyz,Grid_Origin,Gxyz,M2G,Matomnum) private(Mc_AN,Gc_AN,Cxyz,OMPID,Nthrds,Nprocs)
    {

      /* get info. on OpenMP */ 
  
      OMPID = omp_get_thread_num();
      Nthrds = omp_get_num_threads();
      Nprocs = omp_get_num_procs();

      for (Mc_AN=(OMPID+1); Mc_AN<=Matomnum; Mc_AN+=Nthrds){

	Gc_AN = M2G[Mc_AN];
	Cxyz[1] = Gxyz[Gc_AN][1] - Grid_Origin[1];
	Cxyz[2] = Gxyz[Gc_AN][2] - Grid_Origin[2];
	Cxyz[3] = Gxyz[Gc_AN][3] - Grid_Origin[3];
	Cell_Gxyz[Gc_AN][1] = Dot_Product(Cxyz,rtv[1])*0.5/PI;
	Cell_Gxyz[Gc_AN][2] = Dot_Product(Cxyz,rtv[2])*0.5/PI;
	Cell_Gxyz[Gc_AN][3] = Dot_Product(Cxyz,rtv[3])*0.5/PI;

	if (2<=level_stdout){
	  printf("Cell_Gxyz %3d  %15.12f %15.12f %15.12f\n",
		 Gc_AN,
		 Cell_Gxyz[Gc_AN][1],
		 Cell_Gxyz[Gc_AN][2],
		 Cell_Gxyz[Gc_AN][3]); 
	}
      }

    } /* #pragma omp parallel */

    /****************
     MPI: 
         Cell_Gxyz
    *****************/

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      ID = G2ID[ct_AN];
      MPI_Bcast(&Cell_Gxyz[ct_AN][0], 4, MPI_DOUBLE, ID, mpi_comm_level1);
    }
  }

  /****************************************************
            Find grids overlaping to each atom
  ****************************************************/

  /* gtv */

  gtv[1][1] = tv[1][1]/(double)Ngrid1;
  gtv[1][2] = tv[1][2]/(double)Ngrid1;
  gtv[1][3] = tv[1][3]/(double)Ngrid1;

  gtv[2][1] = tv[2][1]/(double)Ngrid2;
  gtv[2][2] = tv[2][2]/(double)Ngrid2;
  gtv[2][3] = tv[2][3]/(double)Ngrid2;

  gtv[3][1] = tv[3][1]/(double)Ngrid3;
  gtv[3][2] = tv[3][2]/(double)Ngrid3;
  gtv[3][3] = tv[3][3]/(double)Ngrid3;

  Cross_Product(gtv[2],gtv[3],tmp);
  GridV = Dot_Product(gtv[1],tmp);
  GridVol = fabs( GridV );

  length_gtv[1] = sqrt( Dot_Product(gtv[1], gtv[1]) );
  length_gtv[2] = sqrt( Dot_Product(gtv[2], gtv[2]) );
  length_gtv[3] = sqrt( Dot_Product(gtv[3], gtv[3]) );

  /* rgtv */

  Cross_Product(gtv[2],gtv[3],tmp);
  rgtv[1][1] = 2.0*PI*tmp[1]/GridV;
  rgtv[1][2] = 2.0*PI*tmp[2]/GridV;
  rgtv[1][3] = 2.0*PI*tmp[3]/GridV;

  Cross_Product(gtv[3],gtv[1],tmp);
  rgtv[2][1] = 2.0*PI*tmp[1]/GridV;
  rgtv[2][2] = 2.0*PI*tmp[2]/GridV;
  rgtv[2][3] = 2.0*PI*tmp[3]/GridV;

  Cross_Product(gtv[1],gtv[2],tmp);
  rgtv[3][1] = 2.0*PI*tmp[1]/GridV;
  rgtv[3][2] = 2.0*PI*tmp[2]/GridV;
  rgtv[3][3] = 2.0*PI*tmp[3]/GridV;

  if (myid==Host_ID && 0<level_stdout){    
    printf("Cell_Volume = %19.12f (Bohr^3)\n",Cell_Volume); 
    printf("GridVol     = %19.12f (Bohr^3)\n",GridVol); 
  }

  if ( (estimate_switch==0 || 2<=level_stdout) && myid==Host_ID ){

    if (0<level_stdout){ 

      printf("Cell vectors (bohr) of the grid cell (gtv)\n");
      printf("  gtv_a = %15.12f, %15.12f, %15.12f\n",gtv[1][1],gtv[1][2],gtv[1][3]);
      printf("  gtv_b = %15.12f, %15.12f, %15.12f\n",gtv[2][1],gtv[2][2],gtv[2][3]);
      printf("  gtv_c = %15.12f, %15.12f, %15.12f\n",gtv[3][1],gtv[3][2],gtv[3][3]);
      printf("  |gtv_a| = %15.12f\n",length_gtv[1]);
      printf("  |gtv_b| = %15.12f\n",length_gtv[2]);
      printf("  |gtv_c| = %15.12f\n",length_gtv[3]);
    }

    /***********************************************
        output informations on grids to a file
    ***********************************************/

    sprintf(file_UCell,"%s%s.UCell",filepath,filename);
    fp = fopen(file_UCell, "a");

#ifdef xt3
    setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

    if (fp!=NULL && estimate_switch==0){

      fprintf(fp,"  Cell_Volume = %19.12f (Bohr^3)\n",Cell_Volume); 
      fprintf(fp,"  GridVol     = %19.12f (Bohr^3)\n",GridVol); 

      fprintf(fp,"  Cell vectors (bohr) of the grid cell (gtv)\n");
      fprintf(fp,"    gtv_a = %15.12f, %15.12f, %15.12f\n",gtv[1][1],gtv[1][2],gtv[1][3]);
      fprintf(fp,"    gtv_b = %15.12f, %15.12f, %15.12f\n",gtv[2][1],gtv[2][2],gtv[2][3]);
      fprintf(fp,"    gtv_c = %15.12f, %15.12f, %15.12f\n",gtv[3][1],gtv[3][2],gtv[3][3]);
      fprintf(fp,"    |gtv_a| = %15.12f\n",length_gtv[1]);
      fprintf(fp,"    |gtv_b| = %15.12f\n",length_gtv[2]);
      fprintf(fp,"    |gtv_c| = %15.12f\n",length_gtv[3]);

      fprintf(fp,"\n***********************************************************\n");
      fprintf(fp,"***********************************************************\n\n");

      fclose(fp); 
    }
  }

  /**********************************
    allocation of arrays: 

    Tmp_GridListAtom
    Tmp_CellListAtom
    MGridListAtom
  **********************************/
 
  Tmp_GridListAtom = (int**)malloc(sizeof(int*)*(Matomnum+MatomnumF+1));
  Tmp_CellListAtom = (int**)malloc(sizeof(int*)*(Matomnum+MatomnumF+1));
  MGridListAtom = (int**)malloc(sizeof(int*)*(Matomnum+MatomnumF+1));
  Tmp_GridListAtom[0] = (int*)malloc(sizeof(int)*1);
  Tmp_CellListAtom[0] = (int*)malloc(sizeof(int)*1);
  MGridListAtom[0] = (int*)malloc(sizeof(int)*1);
  alloc_first[2] = 0;

  /****************************************************
   1) find a neighbouring point of the atom Mc_AN
   2) the ranges which deterinie a box on the atom Mc_AN 
   3) determine whether overlap exists or not
  ****************************************************/

  Mc_AN0 = 1;

#pragma omp parallel shared(myid,ScaleSize,Max_GridN_Atom,time_per_atom,MGridListAtom,Tmp_CellListAtom,Tmp_GridListAtom,GridN_Atom,estimate_switch,Ngrid2,Ngrid3,gtv,Spe_Atom_Cut1,S_Lng,TNumGrid,Gxyz,WhatSpecies,M2G,Mc_AN0,Matomnum) private(OMPID,Nthrds,Nprocs,size_array,TempGrid,TempCell,Mc_AN,Stime_atom,Gc_AN,Cwan,sa_cri,po,MinR,MinN,N,Cxyz,dx,dy,dz,sa,buffer_scale,CutR2,p,r2,Nm,N3,nn1,nn2,nn3,Nc,Nct,n1,n2,n3,NOC,Rn,l1,l2,l3,i,j,Etime_atom)
  {

    /* get info. on OpenMP */ 
  
    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    /* allocation of arrays */

    size_array = (int)(Max_GridN_Atom*ScaleSize);
    TempGrid = (int*)malloc(sizeof(int)*size_array); 
    TempCell = (int*)malloc(sizeof(int)*size_array); 

    do {

#pragma omp barrier
      Mc_AN = Mc_AN0 + OMPID;

      if (Mc_AN<=Matomnum){ 

	dtime(&Stime_atom);
	Gc_AN = M2G[Mc_AN];
	Cwan = WhatSpecies[Gc_AN];

	/* find a neighbouring point of the atom Gc_AN */

	sa_cri = 1.0;
	po = 0;

	do {  

	  MinR = 10000000.0;
	  N = -1;

	  do {

	    N++;

	    Get_Grid_XYZ(N,Cxyz);
	    dx = fabs(Gxyz[Gc_AN][1] - Cxyz[1]);
	    dy = fabs(Gxyz[Gc_AN][2] - Cxyz[2]);
	    dz = fabs(Gxyz[Gc_AN][3] - Cxyz[3]);
	    sa = dx + dy + dz;

	    if (sa<MinR){
	      MinR = sa;
	      MinN = N;
	      if (sa<sa_cri) po = 1;
	    }
	  } while(po==0 && N<TNumGrid);
 
	  if (po==0) sa_cri = 2.0*sa_cri; 

	} while(po==0);

	/* the ranges which determine a box on the atom ct_AN  */

	Get_Grid_XYZ(MinN,Cxyz);
	Cxyz[1] = Cxyz[1] - Gxyz[Gc_AN][1];
	Cxyz[2] = Cxyz[2] - Gxyz[Gc_AN][2];
	Cxyz[3] = Cxyz[3] - Gxyz[Gc_AN][3];

	if      (S_Lng<2.0*Spe_Atom_Cut1[Cwan]) buffer_scale = 16.0;
	else if (S_Lng<4.0*Spe_Atom_Cut1[Cwan]) buffer_scale = 12.0;
	else                                    buffer_scale =  4.0;

	CutR2 = Spe_Atom_Cut1[Cwan]*Spe_Atom_Cut1[Cwan];

	for (p=1; p<=3; p++){

	  po = 0;
	  N = -1;

	  do {
	    N++;
	    dx = (double)N*gtv[p][1] + Cxyz[1];
	    dy = (double)N*gtv[p][2] + Cxyz[2];
	    dz = (double)N*gtv[p][3] + Cxyz[3];
	    r2 = dx*dx + dy*dy + dz*dz; 
	    if ( (buffer_scale*CutR2)<=r2){
	      Nm[p] = N;
	      po = 1;
	    }
	  } while(po==0);
	}

	/* determine whether overlap exists or not  */

	GN2N(MinN,N3);
	nn1 = N3[1];
	nn2 = N3[2];
	nn3 = N3[3];

	Nct = -1;
	for (n1=-Nm[1]+nn1; n1<=(Nm[1]+nn1); n1++){
	  for (n2=-Nm[2]+nn2; n2<=(Nm[2]+nn2); n2++){
	    for (n3=-Nm[3]+nn3; n3<=(Nm[3]+nn3); n3++){
    
	      Find_CGrids(1,n1,n2,n3,Cxyz,NOC);
	      Rn = NOC[0];
	      l1 = NOC[1];
	      l2 = NOC[2];
	      l3 = NOC[3];
	      N = l1*Ngrid2*Ngrid3 + l2*Ngrid3 + l3;
                          
	      dx = Cxyz[1] - Gxyz[Gc_AN][1];
	      dy = Cxyz[2] - Gxyz[Gc_AN][2];
	      dz = Cxyz[3] - Gxyz[Gc_AN][3];

	      r2 = dx*dx + dy*dy + dz*dz;
	      if (r2<=CutR2){
		Nct++;
		if (estimate_switch!=1){
		  TempGrid[Nct+1] = N;
		  TempCell[Nct+1] = Rn;
		}             
	      }

	    }    
	  }    
	}    

        GridN_Atom[Gc_AN] = Nct + 1;

      } /* if (Mc_AN<=Matomnum) */

#pragma omp barrier
#pragma omp flush(GridN_Atom)

      /* allocation of arrays */

      if (OMPID==0){
        for (i=Mc_AN0; i<(Mc_AN0+Nthrds); i++){

          if (i<=Matomnum){
     	    j = M2G[i];
            Tmp_GridListAtom[i] = (int*)malloc(sizeof(int)*GridN_Atom[j]);
            Tmp_CellListAtom[i] = (int*)malloc(sizeof(int)*GridN_Atom[j]);
            MGridListAtom[i] = (int*)malloc(sizeof(int)*GridN_Atom[j]);
	  }

	}
      }

#pragma omp barrier
#pragma omp flush(Tmp_GridListAtom,Tmp_CellListAtom,MGridListAtom)

      if (estimate_switch!=1 && Mc_AN<=Matomnum){

	/* sorting */

	qsort_int((long)(Nct+1),TempGrid,TempCell);

	for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
	  Tmp_GridListAtom[Mc_AN][Nc] = TempGrid[Nc+1];
	  Tmp_CellListAtom[Mc_AN][Nc] = TempCell[Nc+1];
	}
      }

      if (Mc_AN<=Matomnum){
        dtime(&Etime_atom);
        time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
      }

      /* increament of Mc_AN0 */

      if (OMPID==0) Mc_AN0 += Nthrds;
#pragma omp barrier
#pragma omp flush(Mc_AN0)

    } while (Mc_AN0<=Matomnum);

    /* freeing of arrays */

    free(TempCell);
    free(TempGrid);

  } /* #pragma omp parallel */

  /* calculate size_GridListAtom */

  size_GridListAtom = 0;
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    size_GridListAtom += GridN_Atom[Gc_AN];
    if (Max_GridN_Atom<GridN_Atom[Gc_AN]) Max_GridN_Atom = GridN_Atom[Gc_AN];
  }

  /****************************************************
   MPI: 

       GridN_Atom
  ****************************************************/

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    ID = G2ID[ct_AN];
    MPI_Bcast(&GridN_Atom[ct_AN], 1, MPI_INT, ID, mpi_comm_level1);

    if (myid==Host_ID && estimate_switch==0 && 0<level_stdout){
      printf("Num. of grids overlapping with atom %4d = %4d\n",
             ct_AN, GridN_Atom[ct_AN]);
    }
  }

  /****************************************************
    allocation of arrays:

       Tmp_GridListAtom
       Tmp_CellListAtom
  ****************************************************/
  
  size_MGridListAtom = size_GridListAtom;

  for (Mc_AN=Matomnum+1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
    Gc_AN = F_M2G[Mc_AN];
    Tmp_GridListAtom[Mc_AN] = (int*)malloc(sizeof(int)*GridN_Atom[Gc_AN]);
    Tmp_CellListAtom[Mc_AN] = (int*)malloc(sizeof(int)*GridN_Atom[Gc_AN]);
    MGridListAtom[Mc_AN] = (int*)malloc(sizeof(int)*GridN_Atom[Gc_AN]);
    size_MGridListAtom += GridN_Atom[Gc_AN];
  }
  /* PrintMemory */
  if (firsttime){
  PrintMemory("truncation: GridListAtom", sizeof(int)*size_GridListAtom, NULL);
  PrintMemory("truncation: CellListAtom", sizeof(int)*size_MGridListAtom, NULL);
  PrintMemory("truncation: MGridListAtom", sizeof(int)*size_MGridListAtom, NULL);
  }

  /****************************************************
   MPI: 

       Tmp_GridListAtom
       Tmp_CellListAtom
  ****************************************************/

  /* MPI_Barrier */
  MPI_Barrier(mpi_comm_level1);

  /* Tmp_GridListAtom */

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){
      tag = 999;
      /* Sending of data to IDS */
      if (F_Snd_Num[IDS]!=0){
        for (i=0; i<F_Snd_Num[IDS]; i++){
          Mc_AN = Snd_MAN[IDS][i];
          Gc_AN = Snd_GAN[IDS][i];
          MPI_Isend(&Tmp_GridListAtom[Mc_AN][0], GridN_Atom[Gc_AN], MPI_INT,
                    IDS, tag, mpi_comm_level1, &request);
	}
      }

      /* Receiving of data from IDR */
      if (F_Rcv_Num[IDR]!=0){
        tag = 999;
        for (i=0; i<F_Rcv_Num[IDR]; i++){
          Mc_AN = F_TopMAN[IDR] + i;
          Gc_AN = F_M2G[Mc_AN];
          MPI_Recv(&Tmp_GridListAtom[Mc_AN][0], GridN_Atom[Gc_AN], MPI_INT,
                   IDR, tag, mpi_comm_level1, &stat);
	}          
      }

      if (F_Snd_Num[IDS]!=0) MPI_Wait(&request,&stat);
    }     
  }

  /* CellListAtom */

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){
      tag = 999;
      /* Sending of data to IDS */
      if (F_Snd_Num[IDS]!=0){
        for (i=0; i<F_Snd_Num[IDS]; i++){
          Mc_AN = Snd_MAN[IDS][i];
          Gc_AN = Snd_GAN[IDS][i];
          MPI_Isend(&Tmp_CellListAtom[Mc_AN][0], GridN_Atom[Gc_AN], MPI_INT,
                    IDS, tag, mpi_comm_level1, &request);
	}
      }

      /* Receiving of data from IDR */
      if (F_Rcv_Num[IDR]!=0){
        tag = 999;
        for (i=0; i<F_Rcv_Num[IDR]; i++){
          Mc_AN = F_TopMAN[IDR] + i;
          Gc_AN = F_M2G[Mc_AN];
          MPI_Recv(&Tmp_CellListAtom[Mc_AN][0], GridN_Atom[Gc_AN], MPI_INT,
                   IDR, tag, mpi_comm_level1, &stat);
	}          
      }

      if (F_Snd_Num[IDS]!=0) MPI_Wait(&request,&stat);
    }     
  }

  /* MPI_Barrier */
  MPI_Barrier(mpi_comm_level1);

  if (estimate_switch!=1){

    /****************************************************
            Find overlap grids between two orbitals
    ****************************************************/
    
    if (estimate_switch==0){
      size_GListTAtoms0 = 0;

      if (Allocate_TAtoms0==1){
        GListTAtoms0 = (int***)malloc(sizeof(int**)*(Matomnum+1));
        GListTAtoms3 = (int***)malloc(sizeof(int**)*(Matomnum+1));
        GListTCells0 = (int***)malloc(sizeof(int**)*(Matomnum+1));
      }

      GListTAtoms1 = (int***)malloc(sizeof(int**)*(Matomnum+1));
      GListTAtoms2 = (int***)malloc(sizeof(int**)*(Matomnum+1));
      alloc_first[0] = 0;
    }

    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

      dtime(&Stime_atom);

      if (Mc_AN==0) Gc_AN = 0;
      else          Gc_AN = M2G[Mc_AN];

      if (Mc_AN==0){

	FNAN[0] = 0; 

        if (estimate_switch==0){

          if (Allocate_TAtoms0==1){
            GListTAtoms0[0] = (int**)malloc(sizeof(int*)*1);
            GListTAtoms3[0] = (int**)malloc(sizeof(int*)*1);
            GListTCells0[0] = (int**)malloc(sizeof(int*)*1);
	  }

          GListTAtoms1[0] = (int**)malloc(sizeof(int*)*1);
          GListTAtoms2[0] = (int**)malloc(sizeof(int*)*1);

          if (Allocate_TAtoms0==1){
            GListTAtoms0[0][0] = (int*)malloc(sizeof(int)*1);
            GListTAtoms3[0][0] = (int*)malloc(sizeof(int)*1);
            GListTCells0[0][0] = (int*)malloc(sizeof(int)*1);
	  }

          GListTAtoms1[0][0] = (int*)malloc(sizeof(int)*1);
          GListTAtoms2[0][0] = (int*)malloc(sizeof(int)*1);
        }
      }

      else{

        if (estimate_switch==0){

          if (Allocate_TAtoms0==1){
            GListTAtoms0[Mc_AN] = (int**)malloc(sizeof(int*)*(FNAN[Gc_AN]+1));
            GListTAtoms3[Mc_AN] = (int**)malloc(sizeof(int*)*(FNAN[Gc_AN]+1));
            GListTCells0[Mc_AN] = (int**)malloc(sizeof(int*)*(FNAN[Gc_AN]+1));
	  }

          GListTAtoms1[Mc_AN] = (int**)malloc(sizeof(int*)*(FNAN[Gc_AN]+1));
          GListTAtoms2[Mc_AN] = (int**)malloc(sizeof(int*)*(FNAN[Gc_AN]+1));
        }

        h_AN0 = 0;

#pragma omp parallel shared(List_YOUSO,GListTAtoms0,GListTAtoms3,GListTCells0,GListTAtoms1,GListTAtoms2,ScaleSize,Max_NumOLG,size_GListTAtoms0,level_stdout,NumOLG,estimate_switch,CpyCell,Mc_AN,Tmp_CellListAtom,Tmp_GridListAtom,GridN_Atom,atv_ijk,ncn,F_G2M,natn,h_AN0,FNAN,Gc_AN) private(OMPID,Nthrds,Nprocs,h_AN,Gh_AN,Mh_AN,Rh,l1,l2,l3,Nog,Nc,Nh,GNh,GRh,ll1,ll2,ll3,lll1,lll2,lll3,GRh1,po,GNc,GRc,TAtoms0,TCells0,TAtoms1,TAtoms2,size_array,i)
	{        

	  /*******************************************************
            allocation of temporal arrays
	  *******************************************************/

	  size_array = (int)((Max_NumOLG+2)*ScaleSize);
	  TAtoms0 = (int*)malloc(sizeof(int)*size_array);
	  TCells0 = (int*)malloc(sizeof(int)*size_array);
	  TAtoms1 = (int*)malloc(sizeof(int)*size_array);
	  TAtoms2 = (int*)malloc(sizeof(int)*size_array);

	  /* get info. on OpenMP */ 
  
	  OMPID = omp_get_thread_num();
	  Nthrds = omp_get_num_threads();
	  Nprocs = omp_get_num_procs();

          do {  

#pragma omp barrier
            h_AN = h_AN0 + OMPID;

            if (h_AN<=FNAN[Gc_AN]){

	      Gh_AN = natn[Gc_AN][h_AN];
	      Mh_AN = F_G2M[Gh_AN];
	      Rh = ncn[Gc_AN][h_AN];

	      l1 = atv_ijk[Rh][1];
	      l2 = atv_ijk[Rh][2];
	      l3 = atv_ijk[Rh][3];
  
	      Nog = -1;
	      Nc = 0;

	      for (Nh=0; Nh<GridN_Atom[Gh_AN]; Nh++){

		GNh = Tmp_GridListAtom[Mh_AN][Nh];
		GRh = Tmp_CellListAtom[Mh_AN][Nh];

		ll1 = atv_ijk[GRh][1];
		ll2 = atv_ijk[GRh][2];
		ll3 = atv_ijk[GRh][3];

		lll1 = l1 + ll1;
		lll2 = l2 + ll2;
		lll3 = l3 + ll3;

		if (Tmp_GridListAtom[Mc_AN][0]<=GNh) {

		  /* find the initial Nc */

		  if (GNh==0) {
		    Nc = 0;
		  }
		  else {
		    while ( GNh<=Tmp_GridListAtom[Mc_AN][Nc] && Nc!=0 ){
		      Nc = Nc - 10;
		      if (Nc<0) Nc = 0;
		    }
		  }

		  /*  find whether there is the overlapping or not. */

		  if (abs(lll1)<=CpyCell && abs(lll2)<=CpyCell && abs(lll3)<=CpyCell){

		    GRh1 = R_atv(CpyCell,lll1,lll2,lll3);

		    po = 0;

		    while (po==0 && Nc<GridN_Atom[Gc_AN]) {

		      GNc = Tmp_GridListAtom[Mc_AN][Nc];
		      GRc = Tmp_CellListAtom[Mc_AN][Nc];

		      if (GNc==GNh && GRc==GRh1){

			Nog++;

			if (estimate_switch==0){
			  TAtoms0[Nog] = GNc;
			  TCells0[Nog] = GRc;
			  TAtoms1[Nog] = Nc;
			  TAtoms2[Nog] = Nh;
			}

			po = 1;
		      }

		      else if (GNh<GNc){
			po = 1; 
		      } 

		      Nc++;

		    } /* while (...)                      */

		    /* for Nc==GridN_Atom[Gc_AN] */

		    Nc--;
		    if (Nc<0) Nc = 0;

		  } /* if (abs.... )                    */
		} /* if (Tmp_GridListAtom[Mc_AN][0]<=GNh) */
	      } /* Nh */

	      NumOLG[Mc_AN][h_AN] = Nog + 1;

	      if (List_YOUSO[12]<(Nog+1) && estimate_switch==0){
		printf("YOUSO12<(Nog+1)\n");
                MPI_Finalize();
		exit(1);
	      }

	      if (2<=level_stdout){
		printf("Num. of grids overlapping between atoms %2d (G) and %2d (L) = %2d\n",
		       Gc_AN,h_AN,Nog+1);
	      }
	    
	    } /* if (h_AN<=FNAN[Gc_AN]) */

#pragma omp barrier
#pragma omp flush(NumOLG)


	    if (estimate_switch==0){

              /* allocation of arrays */

              if (OMPID==0){

		for (i=h_AN0; i<(h_AN0+Nthrds); i++){

                  if (i<=FNAN[Gc_AN]){

		    size_GListTAtoms0 += NumOLG[Mc_AN][i];

		    if (Allocate_TAtoms0==1){
		      GListTAtoms0[Mc_AN][i] = (int*)malloc(sizeof(int)*NumOLG[Mc_AN][i]);
		      GListTAtoms3[Mc_AN][i] = (int*)malloc(sizeof(int)*NumOLG[Mc_AN][i]);
		      GListTCells0[Mc_AN][i] = (int*)malloc(sizeof(int)*NumOLG[Mc_AN][i]);
		    }

		    GListTAtoms1[Mc_AN][i] = (int*)malloc(sizeof(int)*NumOLG[Mc_AN][i]);
		    GListTAtoms2[Mc_AN][i] = (int*)malloc(sizeof(int)*NumOLG[Mc_AN][i]);

		  }
		}
	      }

#pragma omp barrier
#pragma omp flush(GListTAtoms1,GListTAtoms2)
 	      if (Allocate_TAtoms0==1){
#pragma omp flush(GListTAtoms0,GListTAtoms3,GListTCells0)
	      }

              if (h_AN<=FNAN[Gc_AN]){

		for (Nog=0; Nog<NumOLG[Mc_AN][h_AN]; Nog++){

		  if (Allocate_TAtoms0==1){
		    GListTAtoms0[Mc_AN][h_AN][Nog] = TAtoms0[Nog];
		    GListTCells0[Mc_AN][h_AN][Nog] = TCells0[Nog];
		  }

		  GListTAtoms1[Mc_AN][h_AN][Nog] = TAtoms1[Nog];
		  GListTAtoms2[Mc_AN][h_AN][Nog] = TAtoms2[Nog];

		}
	      }
	    }

            /* increament of h_AN0 */

            if (OMPID==0) h_AN0 += Nthrds;
#pragma omp barrier
#pragma omp flush(h_AN0)

	  } while (h_AN0<=FNAN[Gc_AN]);

          /* freeing of arrays */

	  free(TAtoms2);
	  free(TAtoms1);
	  free(TCells0);
	  free(TAtoms0);

	} /* #pragma omp parallel */

	dtime(&Etime_atom);
	time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
      }

    } /* Mc_AN */ 

    if (estimate_switch==0){

      if (Allocate_TAtoms0==1){

        if (firsttime){
  	PrintMemory("truncation: GListTAtoms0", sizeof(int)*size_GListTAtoms0, NULL);
	PrintMemory("truncation: GListTAtoms3", sizeof(int)*size_GListTAtoms0, NULL);
	PrintMemory("truncation: GListTCells0", sizeof(int)*size_GListTAtoms0, NULL);
	}
      }

      if (firsttime){
      PrintMemory("truncation: GListTAtoms1", sizeof(int)*size_GListTAtoms0, NULL);
      PrintMemory("truncation: GListTAtoms2", sizeof(int)*size_GListTAtoms0, NULL);
      }
    }

    /* find Max_NumOLG */

    My_Max_NumOLG = 0;
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];
      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
        if (My_Max_NumOLG<NumOLG[Mc_AN][h_AN]) My_Max_NumOLG = NumOLG[Mc_AN][h_AN];
      }
    }

    MPI_Allreduce(&My_Max_NumOLG, &Max_NumOLG, 1, MPI_INT, MPI_MAX, mpi_comm_level1);

  } /* if (estimate_switch!=1) */

  /* MPI_Barrier */
  MPI_Barrier(mpi_comm_level1);

  /****************************************************
       Tmp_GridListAtom -> GridListAtom
       Tmp_CellListAtom -> CellListAtom                     
  ****************************************************/

  GridListAtom = (int**)malloc(sizeof(int*)*(Matomnum+1));
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    GridListAtom[Mc_AN] = (int*)malloc(sizeof(int)*GridN_Atom[Gc_AN]);
    for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
      GridListAtom[Mc_AN][Nc] = Tmp_GridListAtom[Mc_AN][Nc];
    }
  }

  CellListAtom = (int**)malloc(sizeof(int*)*(Matomnum+MatomnumF+1));
  for (Mc_AN=1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
    Gc_AN = F_M2G[Mc_AN];
    CellListAtom[Mc_AN] = (int*)malloc(sizeof(int)*GridN_Atom[Gc_AN]);
    for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
      CellListAtom[Mc_AN][Nc] = Tmp_CellListAtom[Mc_AN][Nc];
    }
  }

  /****************************************************
      find grids that each processor has to know
       and 
      setting of grids (intermediate) 
  ****************************************************/

  if (estimate_switch==0){

    allocate_grids2atoms(MD_iter);

    if (Allocate_TAtoms0==1){

      for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

	dtime(&Stime_atom);

	Gc_AN = M2G[Mc_AN];

	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	  for (Nog=0; Nog<NumOLG[Mc_AN][h_AN]; Nog++){
	    GNc = GListTAtoms0[Mc_AN][h_AN][Nog];

	    GN2N(GNc,N3);
	    n1 = N3[1];  
	    n2 = N3[2];  
	    n3 = N3[3];  
	    nn1 = My_Cell0[n1];

	    if (nn1==-1){
	      printf("myid=%i error in truncation.c  n1=%i\n",myid,n1);fflush(stdout);
	      MPI_Finalize();
	      exit(1);
	    }

	    N = nn1*Ngrid2*Ngrid3 + n2*Ngrid3 + n3;
	    GListTAtoms3[Mc_AN][h_AN][Nog] = N;  /* medium grid number in myid */
	  }
	}

	dtime(&Etime_atom);
	time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
      }
    }

    for (Mc_AN=1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){

      dtime(&Stime_atom);

      Gc_AN = F_M2G[Mc_AN];

#pragma omp parallel shared(Ngrid2,Ngrid3,MGridListAtom,My_Cell0,Mc_AN,Tmp_GridListAtom,GridN_Atom,Gc_AN) private(Nc,GNc,N3,n1,n2,n3,nn1,N,OMPID,Nthrds,Nprocs)
      {

	/* get info. on OpenMP */ 
  
	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	for (Nc=OMPID; Nc<GridN_Atom[Gc_AN]; Nc+=Nthrds){

	  GNc = Tmp_GridListAtom[Mc_AN][Nc];
	  GN2N(GNc,N3);
	  n1 = N3[1];  
	  n2 = N3[2];  
	  n3 = N3[3];  
	  nn1 = My_Cell0[n1];
	  if (nn1==-1){
	    MGridListAtom[Mc_AN][Nc] = -1;
	  }
	  else{
	    N = nn1*Ngrid2*Ngrid3 + n2*Ngrid3 + n3;
	    MGridListAtom[Mc_AN][Nc] = N;
	  }
	}

      } /* #pragma omp parallel */

      dtime(&Etime_atom); 
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
    }
  }

  /* for PrintMemory */
  firsttime = 0;

  /****************************************************
                          Free
  ****************************************************/

  /* MPI_Barrier */
  MPI_Barrier(mpi_comm_level1);

  for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
    free(Tmp_CellListAtom[Mc_AN]);
    free(Tmp_GridListAtom[Mc_AN]);
  }
  free(Tmp_CellListAtom);
  free(Tmp_GridListAtom);

  if (alloc_first[2]==0 && estimate_switch!=0){

    for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
      free(MGridListAtom[Mc_AN]);
    }
    free(MGridListAtom);

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      free(GridListAtom[Mc_AN]);
    }
    free(GridListAtom);

    for (Mc_AN=1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
      free(CellListAtom[Mc_AN]);
    }
    free(CellListAtom);
  }

}








int Set_Periodic(int CpyN, int Allocate_switch)
{
  static int firsttime=1;
  int i,j,n,n2,TN;
  long n3;

  TN = (2*CpyN+1)*(2*CpyN+1)*(2*CpyN+1) - 1;

  if (Allocate_switch==0){

    /****************************************************
                         Allocation
    ****************************************************/

    atv = (double**)malloc(sizeof(double*)*(TN+1));
    for (i=0; i<(TN+1); i++){
      atv[i] = (double*)malloc(sizeof(double)*4);
    }

    n = 2*CpyN + 4;
    ratv = (int***)malloc(sizeof(int**)*n);
    for (i=0; i<n; i++){
      ratv[i] = (int**)malloc(sizeof(int*)*n);
      for (j=0; j<n; j++){
	ratv[i][j] = (int*)malloc(sizeof(int)*n);
      }
    }
    if (firsttime) 
    PrintMemory("Set_Periodic: ratv",sizeof(int)*n*n*n,NULL);

    atv_ijk = (int**)malloc(sizeof(int*)*(TN+1));
    for (i=0; i<(TN+1); i++){
      atv_ijk[i] = (int*)malloc(sizeof(int)*4);
    }

    alloc_first[7] = 0;

    /* local */
    n2 = atomnum + 2;
    n3 = ((long)n2*((long)TN+1))/32 + 1;
    asize1 = n3*(long)n2 + 10;
    R_index1 = (long*)malloc(sizeof(long)*(asize1+2));
    for (i=0; i<(asize1+2); i++)  R_index1[i] = 0;
    R_index2 = (long*)malloc(sizeof(long)*(asize1+2));
    for (i=0; i<(asize1+2); i++)  R_index2[i] = 0;

    CellDis = (double**)malloc(sizeof(double*)*(TN+3));
    for (i=0; i<(TN+3); i++){
      CellDis[i] = (double*)malloc(sizeof(double)*4);
      for (j=0; j<4; j++)  CellDis[i][j] = 0.0;
    }

    /* for PrintMemory */
    firsttime=0;

    /****************************************************
           setting of parameters of periodic cells
    ****************************************************/

    Generation_ATV(CpyN);
  }

  else if (Allocate_switch==1){
    /* local */
    n2 = atomnum + 2;
    n3 = ((long)n2*((long)TN+1))/32 + 1;
    asize1 = n3*(long)n2 + 10;
    R_index1 = (long*)malloc(sizeof(long)*(asize1+2));
    R_index2 = (long*)malloc(sizeof(long)*(asize1+2));

    CellDis = (double**)malloc(sizeof(double*)*(TN+3));
    for (i=0; i<(TN+3); i++){
      CellDis[i] = (double*)malloc(sizeof(double)*4);
    }
  }

  /* return */

  return TN;
}











void Free_truncation(int CpyN, int TN, int Free_switch)
{
  int i,j,n;

  if (Free_switch==0){

    if (alloc_first[7]==0){

      for (i=0; i<(TN+1); i++){
        free(atv[i]);
      }
      free(atv);

      n = 2*CpyN + 4;

      for (i=0; i<n; i++){
        for (j=0; j<n; j++){
	  free(ratv[i][j]);
        }
        free(ratv[i]);
      }
      free(ratv);

      for (i=0; i<(TN+1); i++){
        free(atv_ijk[i]);
      }
      free(atv_ijk);
    }

    /* local arrays */

    free(R_index1);
    free(R_index2);

    for (i=0; i<(TN+3); i++){
      free(CellDis[i]);
    }
    free(CellDis);

  }

  else if (Free_switch==1){

    free(R_index1);
    free(R_index2);

    for (i=0; i<(TN+3); i++){
      free(CellDis[i]);
    }
    free(CellDis);

  }
} 








void free_arrays_truncation0()
{
  int i,j,k,m,ct_AN,h_AN,Gh_AN,Hwan;
  int tno0,tno1,tno,Cwan,so,s1,s2;
  int num,wan,n2,wanA,Gi,Mc_AN_GDC;
  int Gc_AN,Mc_AN,nc,ns,spin,fan;
  int Anum,p,vsize,l,NUM,MAnum;
  int numprocs,myid;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
    freeing of arrays:

      H0
      CntH0
      HNL
      iCntHNL
      OLP
      CntOLP
      H
      CntH
      DS_NL
      CntDS_NL
      DM
      ResidualDM
      EDM
      PDM
      IOLP
      CntCoes
      HVNA
      DS_VNA
      HVNA2
      CntHVNA2
      DM_onsite
      v_eff
      NC_OcpN
      NC_v_eff
  ****************************************************/

  if (alloc_first[4]==0){

    /* H0 */

    for (k=0; k<4; k++){
      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	if (Mc_AN==0){
          Gc_AN = 0;
	  tno0 = 1;
          FNAN[0] = 0;
	}
	else{
          Gc_AN = M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	  for (i=0; i<tno0; i++){
	    free(H0[k][Mc_AN][h_AN][i]);
	  }
	  free(H0[k][Mc_AN][h_AN]);
	}
	free(H0[k][Mc_AN]);
      }
      free(H0[k]);
    }
    free(H0);

    /* CntH0 */

    if (Cnt_switch==1){
      for (k=0; k<4; k++){
	for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	    FNAN[0] = 0;
	  }
	  else{
	    Gc_AN = M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_CNO[Cwan];  
	  }    

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	    for (i=0; i<tno0; i++){
	      free(CntH0[k][Mc_AN][h_AN][i]);
	    }
	    free(CntH0[k][Mc_AN][h_AN]);
	  }
	  free(CntH0[k][Mc_AN]);
	}
	free(CntH0[k]);
      }
      free(CntH0);
    }

    /* HNL */

    for (k=0; k<List_YOUSO[5]; k++){
      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	if (Mc_AN==0){
          Gc_AN = 0;
	  tno0 = 1;
          FNAN[0] = 0;
	}
	else{
          Gc_AN = M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	  for (i=0; i<tno0; i++){
	    free(HNL[k][Mc_AN][h_AN][i]);
	  }
	  free(HNL[k][Mc_AN][h_AN]);
	}
	free(HNL[k][Mc_AN]);
      }
      free(HNL[k]);
    }
    free(HNL);

    /* iHNL */

    if ( SpinP_switch==3 ){

      for (k=0; k<List_YOUSO[5]; k++){
	for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF+MatomnumS); Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	    FNAN[0] = 0;
	  }
	  else{
	    Gc_AN = S_M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	  }    

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	    for (i=0; i<tno0; i++){
	      free(iHNL[k][Mc_AN][h_AN][i]);
	    }
	    free(iHNL[k][Mc_AN][h_AN]);
	  }
	  free(iHNL[k][Mc_AN]);
	}
	free(iHNL[k]);
      }
      free(iHNL);
    }

    /* iCntHNL */

    if (SO_switch==1 && Cnt_switch==1){

      for (k=0; k<List_YOUSO[5]; k++){
	for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF+MatomnumS); Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	    FNAN[0] = 0;
	  }
	  else{
	    Gc_AN = S_M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_CNO[Cwan];  
	  }    

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	    for (i=0; i<tno0; i++){
	      free(iCntHNL[k][Mc_AN][h_AN][i]);
	    }
	    free(iCntHNL[k][Mc_AN][h_AN]);
	  }
	  free(iCntHNL[k][Mc_AN]);
	}
	free(iCntHNL[k]);
      }
      free(iCntHNL);
    }

    /* H_Hub  --- added by MJ */  

    if (Hub_U_switch==1 || Constraint_NCS_switch==1 || Zeeman_NCS_switch==1 || Zeeman_NCO_switch==1){

      for (k=0; k<=SpinP_switch; k++){

	FNAN[0] = 0;
	for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	  }
	  else{
	    Gc_AN = M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	  }    

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	    if (Mc_AN==0){
	      tno1 = 1;  
	    }
	    else{
	      Gh_AN = natn[Gc_AN][h_AN];
	      Hwan = WhatSpecies[Gh_AN];
	      tno1 = Spe_Total_NO[Hwan];
	    } 

	    for (i=0; i<tno0; i++){
	      free(H_Hub[k][Mc_AN][h_AN][i]);
	    }
            free(H_Hub[k][Mc_AN][h_AN]);
	  }
          free(H_Hub[k][Mc_AN]);
	}
        free(H_Hub[k]);
      }
      free(H_Hub);
    }

    /* H_Zeeman_NCO */  

    if (Zeeman_NCO_switch==1){

      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	if (Mc_AN==0){
	  Gc_AN = 0;
	  tno0 = 1;
	}
	else{
	  Gc_AN = M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	for (i=0; i<tno0; i++){
	  free(H_Zeeman_NCO[Mc_AN][i]);
	}
        free(H_Zeeman_NCO[Mc_AN]);
      }
      free(H_Zeeman_NCO);
    }

    /* iHNL0 */

    if (SpinP_switch==3){

      for (k=0; k<List_YOUSO[5]; k++){

	FNAN[0] = 0;
	for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	  }
	  else{
	    Gc_AN = M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	  }    

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	    if (Mc_AN==0){
	      tno1 = 1;  
	    }
	    else{
	      Gh_AN = natn[Gc_AN][h_AN];
	      Hwan = WhatSpecies[Gh_AN];
	      tno1 = Spe_Total_NO[Hwan];
	    } 

	    for (i=0; i<tno0; i++){
	      free(iHNL0[k][Mc_AN][h_AN][i]);
	    }
	    free(iHNL0[k][Mc_AN][h_AN]);
	  }
	  free(iHNL0[k][Mc_AN]);
	}
        free(iHNL0[k]);
      }
      free(iHNL0);
    }

    /* OLP_L */  

    for (k=0; k<3; k++){

      FNAN[0] = 0;
      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	if (Mc_AN==0){
	  Gc_AN = 0;
	  tno0 = 1;
	}
	else{
	  Gc_AN = F_M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_NO[Hwan];
	  } 

	  for (i=0; i<tno0; i++){
	    free(OLP_L[k][Mc_AN][h_AN][i]);
	  }
          free(OLP_L[k][Mc_AN][h_AN]);
	}
        free(OLP_L[k][Mc_AN]);
      }
      free(OLP_L[k]);
    }
    free(OLP_L);

    /* OLP */

    for (k=0; k<4; k++){
      for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF+MatomnumS); Mc_AN++){

	if (Mc_AN==0){
          Gc_AN = 0;
	  tno0 = 1;
          FNAN[0] = 0;
	}
	else{
          Gc_AN = S_M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	  for (i=0; i<tno0; i++){
	    free(OLP[k][Mc_AN][h_AN][i]);
	  }
	  free(OLP[k][Mc_AN][h_AN]);
	}
	free(OLP[k][Mc_AN]);
      }
      free(OLP[k]);
    }
    free(OLP);

    /* CntOLP */

    if (Cnt_switch==1){

      for (k=0; k<4; k++){
	for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF+MatomnumS); Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	    FNAN[0] = 0;
	  }
	  else{
	    Gc_AN = S_M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_CNO[Cwan];  
	  }    

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	    for (i=0; i<tno0; i++){
	      free(CntOLP[k][Mc_AN][h_AN][i]);
	    }
	    free(CntOLP[k][Mc_AN][h_AN]);
	  }
	  free(CntOLP[k][Mc_AN]);
	}
	free(CntOLP[k]);
      }
      free(CntOLP);
    }

    /* H */

    for (k=0; k<=SpinP_switch; k++){
      for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF+MatomnumS); Mc_AN++){

	if (Mc_AN==0){
          Gc_AN = 0;
	  tno0 = 1;
	}
	else{
          Gc_AN = S_M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_NO[Hwan];
	  } 

	  for (i=0; i<tno0; i++){
	    free(H[k][Mc_AN][h_AN][i]);
	  }

	  free(H[k][Mc_AN][h_AN]);
	}
	free(H[k][Mc_AN]);
      }
      free(H[k]);
    }
    free(H);

    /* CntH */

    if (Cnt_switch==1){

      for (k=0; k<=SpinP_switch; k++){
	for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF+MatomnumS); Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	  }
	  else{
	    Gc_AN = S_M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_CNO[Cwan];  
	  }    

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	    if (Mc_AN==0){
	      tno1 = 1;  
	    }
	    else{
	      Gh_AN = natn[Gc_AN][h_AN];
	      Hwan = WhatSpecies[Gh_AN];
	      tno1 = Spe_Total_CNO[Hwan];
	    } 

	    for (i=0; i<tno0; i++){
	      free(CntH[k][Mc_AN][h_AN][i]);
	    }

	    free(CntH[k][Mc_AN][h_AN]);
	  }
	  free(CntH[k][Mc_AN]);
	}
	free(CntH[k]);
      }
      free(CntH);
    }

    /* DS_NL */  

    for (so=0; so<(SO_switch+1); so++){
      for (k=0; k<4; k++){
	FNAN[0] = 0;
	for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	  }
	  else{
	    Gc_AN = F_M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	  }    

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	    if (Mc_AN==0){
	      tno1 = 1;  
	    }
	    else{
	      Gh_AN = natn[Gc_AN][h_AN];
	      Hwan = WhatSpecies[Gh_AN];
	      tno1 = Spe_Total_VPS_Pro[Hwan] + 2;
	    } 

	    for (i=0; i<tno0; i++){
	      free(DS_NL[so][k][Mc_AN][h_AN][i]);
	    }
	    free(DS_NL[so][k][Mc_AN][h_AN]);
	  }
	  free(DS_NL[so][k][Mc_AN]);
	}
	free(DS_NL[so][k]);
      }
      free(DS_NL[so]);
    }
    free(DS_NL);

    /* CntDS_NL */  

    if (Cnt_switch==1){

      for (so=0; so<(SO_switch+1); so++){
	for (k=0; k<4; k++){
	  FNAN[0] = 0;
	  for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){

	    if (Mc_AN==0){
	      Gc_AN = 0;
	      tno0 = 1;
	      FNAN[0] = 0;
	    }
	    else{
	      Gc_AN = F_M2G[Mc_AN];
	      Cwan = WhatSpecies[Gc_AN];
	      tno0 = Spe_Total_CNO[Cwan];  
	    }    

	    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	      if (Mc_AN==0){
		tno1 = 1;  
	      }
	      else{
		Gh_AN = natn[Gc_AN][h_AN];
		Hwan = WhatSpecies[Gh_AN];
		tno1 = Spe_Total_VPS_Pro[Hwan] + 2;
	      } 

	      for (i=0; i<tno0; i++){
		free(CntDS_NL[so][k][Mc_AN][h_AN][i]);
	      }
	      free(CntDS_NL[so][k][Mc_AN][h_AN]);
	    }
	    free(CntDS_NL[so][k][Mc_AN]);
	  }
	  free(CntDS_NL[so][k]);
	}
	free(CntDS_NL[so]);
      }
      free(CntDS_NL);
    }

    /* DM */  

    for (m=0; m<List_YOUSO[16]; m++){
      for (k=0; k<=SpinP_switch; k++){
	FNAN[0] = 0;
	for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	  if (Mc_AN==0){
            Gc_AN = 0;
	    tno0 = 1;
            FNAN[0] = 0;
	  }
	  else{
            Gc_AN = M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	  }    

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	    if (Mc_AN==0){
	      tno1 = 1;  
	    }
	    else{
	      Gh_AN = natn[Gc_AN][h_AN];
	      Hwan = WhatSpecies[Gh_AN];
	      tno1 = Spe_Total_NO[Hwan];
	    } 

	    for (i=0; i<tno0; i++){
	      free(DM[m][k][Mc_AN][h_AN][i]);
	    }
            free(DM[m][k][Mc_AN][h_AN]);
	  }
          free(DM[m][k][Mc_AN]);
	}
        free(DM[m][k]);
      }
      free(DM[m]);
    }
    free(DM);

    /* Partial_DM */

    if (cal_partial_charge==1){

      for (k=0; k<=1; k++){

	FNAN[0] = 0;
	for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	  }
	  else{
	    Gc_AN = M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	  }    

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	    if (Mc_AN==0){
	      tno1 = 1;  
	    }
	    else{
	      Gh_AN = natn[Gc_AN][h_AN];
	      Hwan = WhatSpecies[Gh_AN];
	      tno1 = Spe_Total_NO[Hwan];
	    } 

	    for (i=0; i<tno0; i++){
	      free(Partial_DM[k][Mc_AN][h_AN][i]);
	    }
            free(Partial_DM[k][Mc_AN][h_AN]);
	  }
          free(Partial_DM[k][Mc_AN]);
	}
        free(Partial_DM[k]);
      }
      free(Partial_DM);
    }

    /* DM_onsite   --- MJ */  

    if (Hub_U_switch==1 || Constraint_NCS_switch==1 || Zeeman_NCS_switch==1 || Zeeman_NCO_switch==1){

      for (m=0; m<2; m++){
	for (k=0; k<=SpinP_switch; k++){

	  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
	  
	    if (Mc_AN==0){
	      Gc_AN = 0;
	      tno0 = 1;
	    }
	    else{
	      Gc_AN = M2G[Mc_AN];
	      Cwan = WhatSpecies[Gc_AN];
	      tno0 = Spe_Total_NO[Cwan];  
	    }  

	    h_AN = 0;
	      
	    if (Mc_AN==0){
	      tno1 = 1;  
	    }
	    else{
	      Gh_AN = natn[Gc_AN][h_AN];
	      Hwan = WhatSpecies[Gh_AN];
	      tno1 = Spe_Total_NO[Hwan];
	    } 

	    for (i=0; i<tno0; i++){
	      free(DM_onsite[m][k][Mc_AN][i]);
	    }

            free(DM_onsite[m][k][Mc_AN]);
	  }
          free(DM_onsite[m][k]);
	}
        free(DM_onsite[m]);
      }
      free(DM_onsite);
    }

    /* v_eff   --- MJ */  

    if (Hub_U_switch==1 || Constraint_NCS_switch==1 || Zeeman_NCS_switch==1 || Zeeman_NCO_switch==1){
   
      for (k=0; k<=SpinP_switch; k++){

	FNAN[0] = 0;

        for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	  }
	  else{
	    Gc_AN = F_M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	  }  

	  for (i=0; i<tno0; i++){
	    free(v_eff[k][Mc_AN][i]);
	  }
	  free(v_eff[k][Mc_AN]);

	}
	free(v_eff[k]);
      }
      free(v_eff);
    }

    /*  NC_OcpN */

    if ( (Hub_U_switch==1 || Constraint_NCS_switch==1 || Zeeman_NCS_switch==1 || Zeeman_NCO_switch==1) 
         && SpinP_switch==3 ){

      for (m=0; m<2; m++){
	for (s1=0; s1<2; s1++){
	  for (s2=0; s2<2; s2++){
	    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	      if (Mc_AN==0){
		Gc_AN = 0;
		tno0 = 1;
	      }
	      else{
		Gc_AN = M2G[Mc_AN];
		Cwan = WhatSpecies[Gc_AN];
		tno0 = Spe_Total_NO[Cwan];  
	      }  

	      for (i=0; i<tno0; i++){
                free(NC_OcpN[m][s1][s2][Mc_AN][i]);
	      }
              free(NC_OcpN[m][s1][s2][Mc_AN]);
	    }
            free(NC_OcpN[m][s1][s2]);
	  }
          free(NC_OcpN[m][s1]);
	}
        free(NC_OcpN[m]);
      }
      free(NC_OcpN);
    }

    /*  NC_v_eff */

    if ( (Hub_U_switch==1 || Constraint_NCS_switch==1 || Zeeman_NCS_switch==1 || Zeeman_NCO_switch==1)
         && SpinP_switch==3 ){

      for (s1=0; s1<2; s1++){
	for (s2=0; s2<2; s2++){
	  for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){

	    if (Mc_AN==0){
	      Gc_AN = 0;
	      tno0 = 1;
	    }
	    else{
	      Gc_AN = F_M2G[Mc_AN];
	      Cwan = WhatSpecies[Gc_AN];
	      tno0 = Spe_Total_NO[Cwan];  
	    }  

	    for (i=0; i<tno0; i++){
	      free(NC_v_eff[s1][s2][Mc_AN][i]);
	    }
            free(NC_v_eff[s1][s2][Mc_AN]);
	  }
          free(NC_v_eff[s1][s2]);
	}
        free(NC_v_eff[s1]);
      }           
      free(NC_v_eff);
    }

    /* ResidualDM */  

    if ( Mixing_switch==0 || Mixing_switch==1 || Mixing_switch==2 ){

      for (m=0; m<List_YOUSO[16]; m++){
	for (k=0; k<=SpinP_switch; k++){
	  FNAN[0] = 0;
	  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	    if (Mc_AN==0){
	      Gc_AN = 0;
	      tno0 = 1;
	    }
	    else{
	      Gc_AN = M2G[Mc_AN];
	      Cwan = WhatSpecies[Gc_AN];
	      tno0 = Spe_Total_NO[Cwan];  
	    }    

	    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	      if (Mc_AN==0){
		tno1 = 1;  
	      }
	      else{
		Gh_AN = natn[Gc_AN][h_AN];
		Hwan = WhatSpecies[Gh_AN];
		tno1 = Spe_Total_NO[Hwan];
	      } 

	      for (i=0; i<tno0; i++){
		free(ResidualDM[m][k][Mc_AN][h_AN][i]);
	      }
	      free(ResidualDM[m][k][Mc_AN][h_AN]);
	    }
	    free(ResidualDM[m][k][Mc_AN]);
	  }
	  free(ResidualDM[m][k]);
	}
	free(ResidualDM[m]);
      }
      free(ResidualDM);
    }

    /* iResidualDM */  

    if ( (Mixing_switch==0 || Mixing_switch==1 || Mixing_switch==2)
	 && SpinP_switch==3 && ( SO_switch==1 || Hub_U_switch==1 || Constraint_NCS_switch==1 
         || Zeeman_NCS_switch==1 || Zeeman_NCO_switch==1) ){

      for (m=0; m<List_YOUSO[16]; m++){
        for (k=0; k<2; k++){
	  FNAN[0] = 0;
	  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	    if (Mc_AN==0){
	      Gc_AN = 0;
	      tno0 = 1;
	    }
	    else{
	      Gc_AN = M2G[Mc_AN];
	      Cwan = WhatSpecies[Gc_AN];
	      tno0 = Spe_Total_NO[Cwan];  
	    }    

	    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	      if (Mc_AN==0){
		tno1 = 1;  
	      }
	      else{
		Gh_AN = natn[Gc_AN][h_AN];
		Hwan = WhatSpecies[Gh_AN];
		tno1 = Spe_Total_NO[Hwan];
	      } 

	      for (i=0; i<tno0; i++){
		free(iResidualDM[m][k][Mc_AN][h_AN][i]);
	      }
 	      free(iResidualDM[m][k][Mc_AN][h_AN]);
	    }
            free(iResidualDM[m][k][Mc_AN]);
	  }
          free(iResidualDM[m][k]);
	}
        free(iResidualDM[m]);
      }
      free(iResidualDM);
    }
    else{
      for (m=0; m<List_YOUSO[16]; m++){
        free(iResidualDM[m][0][0][0][0]);
        free(iResidualDM[m][0][0][0]);
        free(iResidualDM[m][0][0]);
        free(iResidualDM[m][0]);
        free(iResidualDM[m]);
      }   
      free(iResidualDM);
    }

    /* EDM */

    for (k=0; k<=SpinP_switch; k++){
      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	if (Mc_AN==0){
          Gc_AN = 0;
	  tno0 = 1;
	}
	else{
          Gc_AN = M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_NO[Hwan];
	  } 

	  for (i=0; i<tno0; i++){
	    free(EDM[k][Mc_AN][h_AN][i]);
	  }

	  free(EDM[k][Mc_AN][h_AN]);
	}
	free(EDM[k][Mc_AN]);
      }
      free(EDM[k]);
    }
    free(EDM);

    /* PDM */

    for (k=0; k<=SpinP_switch; k++){
      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	if (Mc_AN==0){
          Gc_AN = 0;
	  tno0 = 1;
	}
	else{
          Gc_AN = M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_NO[Hwan];
	  } 

	  for (i=0; i<tno0; i++){
	    free(PDM[k][Mc_AN][h_AN][i]);
	  }

	  free(PDM[k][Mc_AN][h_AN]);
	}
	free(PDM[k][Mc_AN]);
      }
      free(PDM[k]);
    }
    free(PDM);

    /* iDM */

    for (m=0; m<List_YOUSO[16]; m++){
      for (k=0; k<2; k++){

	FNAN[0] = 0;
	for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	  }
	  else{
	    Gc_AN = M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	  }    

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	    if (Mc_AN==0){
	      tno1 = 1;  
	    }
	    else{
	      Gh_AN = natn[Gc_AN][h_AN];
	      Hwan = WhatSpecies[Gh_AN];
	      tno1 = Spe_Total_NO[Hwan];
	    } 

	    for (i=0; i<tno0; i++){
	      free(iDM[m][k][Mc_AN][h_AN][i]);
	    }
	    free(iDM[m][k][Mc_AN][h_AN]);
	  }
	  free(iDM[m][k][Mc_AN]);
	}
	free(iDM[m][k]);
      }
      free(iDM[m]);
    }
    free(iDM);

    /* IOLP */  

    if (Solver==1 || Solver==7){
      FNAN[0] = 0;
      SNAN[0] = 0;
      for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF+MatomnumS); Mc_AN++){

	if (Mc_AN==0){
	  Gc_AN = 0;
	  tno0 = 1;
	}
	else{
	  Gc_AN = S_M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	for (h_AN=0; h_AN<=(FNAN[Gc_AN]+SNAN[Gc_AN]); h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_NO[Hwan];
	  } 

	  for (i=0; i<tno0; i++){
	    free(IOLP[Mc_AN][h_AN][i]);
	  }
	  free(IOLP[Mc_AN][h_AN]);
	}
	free(IOLP[Mc_AN]);
      }
      free(IOLP);
    }

    /* S12 for DC or recursion */  

    if (Solver==1 || Solver==5){

      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	if (Mc_AN==0) n2 = 1;
	else{
	  Gc_AN = M2G[Mc_AN];
	  wan = WhatSpecies[Gc_AN];

	  num = 1;
	  for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
	    Gi = natn[Gc_AN][i];
	    wanA = WhatSpecies[Gi];
	    num += Spe_Total_CNO[wanA];
	  }
	  n2 = num + 2;
	}

	for (i=0; i<n2; i++){
	  free(S12[Mc_AN][i]);
	}
        free(S12[Mc_AN]);
      }
      free(S12);
    }

    else if (Solver==6){ /* for GDC */

      for (Mc_AN_GDC=0; Mc_AN_GDC<=Matomnum_GDC; Mc_AN_GDC++){

	if (Mc_AN_GDC==0) n2 = 1;
	else{

	  Mc_AN = Mnatn_GDC[Mc_AN_GDC][0];
	  Gc_AN = M2G[Mc_AN];
	  wan = WhatSpecies[Gc_AN];

	  num = 1;
	  for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
	    Gi = natn[Gc_AN][i];
	    wanA = WhatSpecies[Gi];
	    num += Spe_Total_CNO[wanA];
	  }
	  n2 = num + 2;
	}

	for (i=0; i<n2; i++){
	  free(S12[Mc_AN_GDC][i]);
	}
        free(S12[Mc_AN_GDC]);
      }
      free(S12);
    }

    if (Solver==1) { /* for recursion */

      /* Left_U0 */

      for (i=0; i<=SpinP_switch; i++){
	for (j=0; j<=Matomnum; j++){

	  if (j==0){
	    tno1 = 1;
	    vsize = 1;
	  }
	  else{
	    Gc_AN = M2G[j];
	    wanA = WhatSpecies[Gc_AN];
	    tno1 = Spe_Total_CNO[wanA];

	    Anum = 1;
	    for (p=0; p<=(FNAN[Gc_AN]+SNAN[Gc_AN]); p++){
	      Gi = natn[Gc_AN][p];
	      wanA = WhatSpecies[Gi];
	      Anum += Spe_Total_CNO[wanA];
	    }

	    NUM = Anum - 1;
	    vsize = NUM + 2;
	  }

	  for (k=0; k<List_YOUSO[3]; k++){
	    for (l=0; l<tno1; l++){
	      free(Left_U0[i][j][k][l]);
	    }
            free(Left_U0[i][j][k]);
	  }
          free(Left_U0[i][j]);
	}
        free(Left_U0[i]);
      }
      free(Left_U0);

      /* Right_U0 */

      for (i=0; i<=SpinP_switch; i++){
	for (j=0; j<=Matomnum; j++){

	  if (j==0){
	    tno1 = 1;
	    vsize = 1;
	  }
	  else{
	    Gc_AN = M2G[j];
	    wanA = WhatSpecies[Gc_AN];
	    tno1 = Spe_Total_CNO[wanA];

	    Anum = 1;
	    for (p=0; p<=(FNAN[Gc_AN]+SNAN[Gc_AN]); p++){
	      Gi = natn[Gc_AN][p];
	      wanA = WhatSpecies[Gi];
	      Anum += Spe_Total_CNO[wanA];
	    }

	    NUM = Anum - 1;
	    vsize = NUM + 2;
	  }

	  for (k=0; k<List_YOUSO[3]; k++){
	    for (l=0; l<tno1; l++){
	      free(Right_U0[i][j][k][l]);
	    }
            free(Right_U0[i][j][k]);
	  }
          free(Right_U0[i][j]);
	}
        free(Right_U0[i]);
      }
      free(Right_U0);
    }

    /* CntCoes */

    if (Cnt_switch==1){
      for (i=0; i<=(Matomnum+MatomnumF); i++){
	for (j=0; j<List_YOUSO[7]; j++){
	  free(CntCoes[i][j]);
	}
	free(CntCoes[i]);
      }
      free(CntCoes);

      for (i=0; i<(SpeciesNum+1); i++){
	for (j=0; j<List_YOUSO[7]; j++){
	  free(CntCoes_Species[i][j]);
	}
	free(CntCoes_Species[i]);
      }
      free(CntCoes_Species);
    }

    if (ProExpn_VNA==1){

      /* HVNA */  

      FNAN[0] = 0;
      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	if (Mc_AN==0){
	  Gc_AN = 0;
	  tno0 = 1;
	}
	else{
	  Gc_AN = M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_NO[Cwan];  
	}    

	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_NO[Hwan];
	  } 

	  for (i=0; i<tno0; i++){
	    free(HVNA[Mc_AN][h_AN][i]);
	  }
          free(HVNA[Mc_AN][h_AN]);
	}
        free(HVNA[Mc_AN]);
      }
      free(HVNA);

      /* DS_VNA */  

      for (k=0; k<4; k++){

	FNAN[0] = 0;
	for (Mc_AN=0; Mc_AN<(Matomnum+2); Mc_AN++){
          
	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	    fan = FNAN[Gc_AN];
	  }
  	  else if ( (Matomnum+1)<=Mc_AN ){
	    fan = List_YOUSO[8];
	    tno0 = List_YOUSO[7];
	  }
	  else{
	    Gc_AN = F_M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	    fan = FNAN[Gc_AN];
	  }
          
	  for (h_AN=0; h_AN<(fan+1); h_AN++){

	    if (Mc_AN==0){
	      tno1 = 1;  
	    }
	    else{
	      tno1 = (List_YOUSO[35]+1)*(List_YOUSO[35]+1)*List_YOUSO[34];
	    } 

	    for (i=0; i<tno0; i++){
	      free(DS_VNA[k][Mc_AN][h_AN][i]);
	    }
            free(DS_VNA[k][Mc_AN][h_AN]);
	  }
          free(DS_VNA[k][Mc_AN]);
	}
        free(DS_VNA[k]);
      }
      free(DS_VNA);

      /* CntDS_VNA */  

      if (Cnt_switch==1){

	for (k=0; k<4; k++){

	  FNAN[0] = 0;
	  for (Mc_AN=0; Mc_AN<(Matomnum+2); Mc_AN++){
          
	    if (Mc_AN==0){
	      Gc_AN = 0;
	      tno0 = 1;
	      fan = FNAN[Gc_AN];
	    }
	    else if ( Mc_AN==(Matomnum+1) ){
	      fan = List_YOUSO[8];
	      tno0 = List_YOUSO[7];
	    }
	    else{
	      Gc_AN = F_M2G[Mc_AN];
	      Cwan = WhatSpecies[Gc_AN];
	      tno0 = Spe_Total_CNO[Cwan];  
              fan = FNAN[Gc_AN];
	    }

	    for (h_AN=0; h_AN<(fan+1); h_AN++){

	      if (Mc_AN==0){
		tno1 = 1;  
	      }
	      else{
		tno1 = (List_YOUSO[35]+1)*(List_YOUSO[35]+1)*List_YOUSO[34];
	      } 

	      for (i=0; i<tno0; i++){
		free(CntDS_VNA[k][Mc_AN][h_AN][i]);
	      }
 	      free(CntDS_VNA[k][Mc_AN][h_AN]);
	    }
            free(CntDS_VNA[k][Mc_AN]);
	  }
          free(CntDS_VNA[k]);
	}
        free(CntDS_VNA);
      }

      /* HVNA2 */

      for (k=0; k<4; k++){
	FNAN[0] = 0;
	for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){

	  if (Mc_AN==0){
	    Gc_AN = 0;
	    tno0 = 1;
	  }
	  else{
	    Gc_AN = F_M2G[Mc_AN];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	  }    

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	    for (i=0; i<tno0; i++){
	      free(HVNA2[k][Mc_AN][h_AN][i]);
	    }
	    free(HVNA2[k][Mc_AN][h_AN]);
	  }
	  free(HVNA2[k][Mc_AN]);
	}
	free(HVNA2[k]);
      }
      free(HVNA2);

      /* CntHVNA2 */

      if (Cnt_switch==1){

	for (k=0; k<4; k++){
	  FNAN[0] = 0;
	  for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){

	    if (Mc_AN==0){
	      Gc_AN = 0;
	      tno0 = 1;
	    }
	    else{
	      Gc_AN = F_M2G[Mc_AN];
	      Cwan = WhatSpecies[Gc_AN];
	      tno0 = Spe_Total_CNO[Cwan];  
	    }    

	    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	      for (i=0; i<tno0; i++){
		free(CntHVNA2[k][Mc_AN][h_AN][i]);
	      }
	      free(CntHVNA2[k][Mc_AN][h_AN]);
	    }
	    free(CntHVNA2[k][Mc_AN]);
	  }
	  free(CntHVNA2[k]);
	}
	free(CntHVNA2);
      }

    }

    if (Solver==8) { /* Krylov subspace method */

      for (i=0; i<=SpinP_switch; i++){
	for (j=0; j<=Matomnum; j++){

	  if (j==0){
	    tno0 = 1;
	  }
	  else{
	    Gc_AN = M2G[j];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	  }    

          for (k=0; k<rlmax_EC[j]; k++){
	    for (l=0; l<EKC_core_size[j]; l++){
	      free(Krylov_U[i][j][k][l]);
	    }
            free(Krylov_U[i][j][k]);
	  }
          free(Krylov_U[i][j]);
	}
        free(Krylov_U[i]);
      }
      free(Krylov_U);

      for (i=0; i<=SpinP_switch; i++){
	for (j=0; j<=Matomnum; j++){

	  if (j==0){
	    tno0 = 1;
	  }
	  else{
	    Gc_AN = M2G[j];
	    Cwan = WhatSpecies[Gc_AN];
	    tno0 = Spe_Total_NO[Cwan];  
	  }    

	  for (k=0; k<(rlmax_EC[j]*EKC_core_size[j]+1); k++){
	    free(EC_matrix[i][j][k]);
	  }
          free(EC_matrix[i][j]);
	}
        free(EC_matrix[i]);
      }
      free(EC_matrix);

      free(rlmax_EC);
      free(rlmax_EC2);
      free(EKC_core_size);
      free(scale_rc_EKC);

    } /* if (Solver==8) */

    /* NEGF */

    if (Solver==4){
      for (spin=0; spin<(SpinP_switch+1); spin++){
	for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
	  free(TRAN_DecMulP[spin][Mc_AN]);
	}
	free(TRAN_DecMulP[spin]);
      }
      free(TRAN_DecMulP);
    }

  } /*  if (alloc_first[4]==0){ */

  /****************************************************
    freeing of arrays:

      GListTAtoms0;
      GListTCells0;
      GListTAtoms1;
      GListTAtoms2;
      GListTAtoms3;

      Density_Grid
      ADensity_Grid
      PCCDensity_Grid
      Vxc_Grid
      RefVxc_Grid
      VNA_Grid
      dVHart_Grid
      Vpot_Grid
      Orbs_Grid
      COrbs_Grid
      VEF_Grid
  ****************************************************/

  if (alloc_first[0]==0){

    FNAN[0] = 0; 
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

      if (Mc_AN==0) Gc_AN = 0;
      else          Gc_AN = M2G[Mc_AN];

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
        free(GListTAtoms2[Mc_AN][h_AN]);
        free(GListTAtoms1[Mc_AN][h_AN]);

        if (Allocate_TAtoms0==1){
          free(GListTCells0[Mc_AN][h_AN]);
          free(GListTAtoms3[Mc_AN][h_AN]);
          free(GListTAtoms0[Mc_AN][h_AN]);
	}
      }
      free(GListTAtoms2[Mc_AN]);
      free(GListTAtoms1[Mc_AN]);

      if (Allocate_TAtoms0==1){
        free(GListTCells0[Mc_AN]);
        free(GListTAtoms3[Mc_AN]);
        free(GListTAtoms0[Mc_AN]);
      }
    }
    free(GListTAtoms2);
    free(GListTAtoms1);

    if (Allocate_TAtoms0==1){
      free(GListTCells0);
      free(GListTAtoms3);
      free(GListTAtoms0);
    }
  }

  if (alloc_first[3]==0){

    if (SpinP_switch==3){ /* spin non-collinear */
      for (k=0; k<=3; k++){
        free(Density_Grid[k]);
      }
      free(Density_Grid);
    }
    else{ 
      for (k=0; k<=1; k++){
        free(Density_Grid[k]);
      }
      free(Density_Grid);
    }

    free(ADensity_Grid);
    free(PCCDensity_Grid);

    if (SpinP_switch==3){ /* spin non-collinear */
      for (k=0; k<=3; k++){
        free(Vxc_Grid[k]);
      }
      free(Vxc_Grid);
    }
    else{
      for (k=0; k<=1; k++){
        free(Vxc_Grid[k]);
      }
      free(Vxc_Grid);
    }

    free(RefVxc_Grid);

    free(VNA_Grid);

    free(dVHart_Grid);

    if (SpinP_switch==3){ /* spin non-collinear */
      for (k=0; k<=3; k++){
        free(Vpot_Grid[k]);
      }
      free(Vpot_Grid);
    }
    else{
      for (k=0; k<=1; k++){
        free(Vpot_Grid[k]);
      }
      free(Vpot_Grid);
    }

    /* external electric field */
    free(VEF_Grid);

    /* Orbs_Grid */

    for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
      if (Mc_AN==0){
        tno = 1;
        Gc_AN = 0;
      }
      else{
        Gc_AN = F_M2G[Mc_AN];
        Cwan = WhatSpecies[Gc_AN];
        tno = Spe_Total_NO[Cwan];
      }

      for (i=0; i<tno; i++){
        free(Orbs_Grid[Mc_AN][i]); 
      }
      free(Orbs_Grid[Mc_AN]); 
    }
    free(Orbs_Grid); 

    /* COrbs_Grid */

    if (Cnt_switch!=0){
      for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
        if (Mc_AN==0){
          tno = 1;
          Gc_AN = 0;
        }
        else{
          Gc_AN = F_M2G[Mc_AN];
          Cwan = WhatSpecies[Gc_AN];
          tno = Spe_Total_CNO[Cwan];
        }

        for (i=0; i<tno; i++){
          free(COrbs_Grid[Mc_AN][i]); 
        }
        free(COrbs_Grid[Mc_AN]); 
      }
      free(COrbs_Grid);
    }

  }

  /****************************************************
    freeing of arrays:

     NumOLG
  ****************************************************/

  if (alloc_first[5]==0){

    /* NumOLG */
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      free(NumOLG[Mc_AN]);
    }
    free(NumOLG);

  } /* if (alloc_first[5]==0){ */

  /****************************************************
    freeing of arrays:

      RMI1
      RMI2
  ****************************************************/

  if (alloc_first[6]==0){

    FNAN[0] = 0;
    SNAN[0] = 0;

    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      if (Mc_AN==0) Gc_AN = 0;
      else          Gc_AN = M2G[Mc_AN];
      for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
        free(RMI1[Mc_AN][i]);
      }      
      free(RMI1[Mc_AN]);
    }
    free(RMI1);

    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      if (Mc_AN==0) Gc_AN = 0;
      else          Gc_AN = M2G[Mc_AN];
      for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
        free(RMI2[Mc_AN][i]);
      }      
      free(RMI2[Mc_AN]);
    }
    free(RMI2);

  } /* if (alloc_first[6]==0){ */

  /****************************************************
    freeing of arrays:

      ratv
      atv
      atv_ijk
  ****************************************************/

  if (alloc_first[7]==0){

    for (i=0; i<(TCpyCell+1); i++){
      free(atv[i]);
    }
    free(atv);

    for (i=0; i<(2*CpyCell+4); i++){
      for (j=0; j<(2*CpyCell+4); j++){
	free(ratv[i][j]);
      }
      free(ratv[i]);
    }
    free(ratv);

    for (i=0; i<(TCpyCell+1); i++){
      free(atv_ijk[i]);
    }
    free(atv_ijk);

  } /* if (alloc_first[7]==0){ */

  /**********************************
    freeing of arrays: 
      GridListAtom
      CellListAtom
      MGridListAtom
  **********************************/

  if (alloc_first[2]==0){

    for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
      free(MGridListAtom[Mc_AN]);
    }
    free(MGridListAtom);

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      free(GridListAtom[Mc_AN]);
    }
    free(GridListAtom);

    for (Mc_AN=1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
      free(CellListAtom[Mc_AN]);
    }
    free(CellListAtom);
  }

  /**********************************
    freeing of arrays: 
      F_M2G
      S_M2G
  **********************************/

  if (alloc_first[13]==0){
    free(F_M2G);
    free(S_M2G);
  }

}



void Set_Inf_SndRcv()
{ 
  int i,ID,IDS,IDR,Mc_AN,Gc_AN,Num,ID1,Lh_AN,Gh_AN;
  int myid,numprocs,tag=999;
  int **po_ID;
  int **Rcv_FGAN,**Rcv_SGAN;
  int **Snd_FGAN,**Snd_SGAN;
  int *Num_Pro_Snd;

  int Rn,Rn2,m1,m2,m3,n1,n2,n3,j,po,Gj_AN;
  MPI_Status stat;
  MPI_Request request;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /*********************************
   allocation of arrays:

   int po_ID[atomnum+1][numprocs]
  *********************************/

  po_ID = (int**)malloc(sizeof(int*)*(atomnum+1));
  for (Gc_AN=0; Gc_AN<(atomnum+1); Gc_AN++){
    po_ID[Gc_AN] = (int*)malloc(sizeof(int)*numprocs);
  }

  /*********************************
   initialize

      F_Snd_Num
      F_Rcv_Num
      S_Snd_Num
      S_Rcv_Num
  *********************************/

  for (ID=0; ID<numprocs; ID++){
    F_Snd_Num[ID] = 0;
    F_Rcv_Num[ID] = 0;
    S_Snd_Num[ID] = 0;
    S_Rcv_Num[ID] = 0;
  }
    
  /************************************************
      find F_Rcv_Num and S_Rcv_Num
  *************************************************/

  /* initialize po_ID */

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    for (ID=0; ID<numprocs; ID++){
      po_ID[Gc_AN][ID] = 0;
    }
  }

  /* find F_Rcv_Num */

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    for (Lh_AN=0; Lh_AN<=FNAN[Gc_AN]; Lh_AN++){
      Gh_AN = natn[Gc_AN][Lh_AN];
      ID1 = G2ID[Gh_AN];
      if (po_ID[Gh_AN][ID1]==0 && ID1!=myid){
        F_Rcv_Num[ID1]++;
        po_ID[Gh_AN][ID1] = 1;
      }
    }
  }

  /* find S_Rcv_Num */

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    for (Lh_AN=(FNAN[Gc_AN]+1); Lh_AN<=(FNAN[Gc_AN]+SNAN[Gc_AN]); Lh_AN++){
      Gh_AN = natn[Gc_AN][Lh_AN];
      ID1 = G2ID[Gh_AN];
      if (po_ID[Gh_AN][ID1]==0 && ID1!=myid){
        S_Rcv_Num[ID1]++;
        po_ID[Gh_AN][ID1] = 1;
      }
    }
  }

  /************************************************
   allocation of array:

   int Rcv_GAN[numprocs]
              [F_Rcv_Num[ID]+S_Rcv_Num[ID]]
  *************************************************/

  if (alloc_first[12]==0){
    for (ID=0; ID<numprocs; ID++){
      free(Rcv_GAN[ID]); 
    }
    free(Rcv_GAN); 
  }

  Rcv_GAN = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Rcv_GAN[ID] = (int*)malloc(sizeof(int)*(F_Rcv_Num[ID]+S_Rcv_Num[ID]));
  }

  alloc_first[12] = 0;

  Rcv_FGAN = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Rcv_FGAN[ID] = (int*)malloc(sizeof(int)*F_Rcv_Num[ID]);
  }

  Rcv_SGAN = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Rcv_SGAN[ID] = (int*)malloc(sizeof(int)*S_Rcv_Num[ID]);
  }

  /************************************************
             set Rcv_FGAN and Rcv_SGAN
  *************************************************/
  
  /* initialize F_Rcv_Num and S_Rcv_Num */

  for (ID=0; ID<numprocs; ID++){
    F_Rcv_Num[ID] = 0;
    S_Rcv_Num[ID] = 0;
  }

  /* initialized po_ID */

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    for (ID=0; ID<numprocs; ID++){
      po_ID[Gc_AN][ID] = 0;
    }
  }

  /* set Rcv_FGAN */

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    for (Lh_AN=0; Lh_AN<=FNAN[Gc_AN]; Lh_AN++){
      Gh_AN = natn[Gc_AN][Lh_AN];
      ID1 = G2ID[Gh_AN];
      if (po_ID[Gh_AN][ID1]==0 && ID1!=myid){

        Rcv_FGAN[ID1][F_Rcv_Num[ID1]] = Gh_AN;
        F_Rcv_Num[ID1]++;
        po_ID[Gh_AN][ID1] = 1;
      }
    }
  }

  /* set Rcv_SGAN */

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    for (Lh_AN=(FNAN[Gc_AN]+1); Lh_AN<=(FNAN[Gc_AN]+SNAN[Gc_AN]); Lh_AN++){
      Gh_AN = natn[Gc_AN][Lh_AN];
      ID1 = G2ID[Gh_AN];
      if (po_ID[Gh_AN][ID1]==0 && ID1!=myid){

        Rcv_SGAN[ID1][S_Rcv_Num[ID1]] = Gh_AN;
        S_Rcv_Num[ID1]++;
        po_ID[Gh_AN][ID1] = 1;
      }
    }
  }

  /*****************************************
       MPI:  F_Rcv_Num
  *****************************************/

  /* MPI_Barrier */
  MPI_Barrier(mpi_comm_level1);
  
  tag = 999;
  for (ID=0; ID<numprocs; ID++){
    if (ID!=0){
      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;
      MPI_Isend(&F_Rcv_Num[IDS],1,MPI_INT,IDS,tag,mpi_comm_level1,&request);
      MPI_Recv( &F_Snd_Num[IDR],1,MPI_INT,IDR,tag,mpi_comm_level1,&stat); 
      MPI_Wait(&request,&stat);
    }
  }

  Snd_FGAN = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Snd_FGAN[ID] = (int*)malloc(sizeof(int)*F_Snd_Num[ID]);
  }

  for (ID=0; ID<numprocs; ID++){
    if (ID!=0){
      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;
      MPI_Isend(&Rcv_FGAN[IDS][0],F_Rcv_Num[IDS],MPI_INT,IDS,tag,mpi_comm_level1,&request);
      MPI_Recv( &Snd_FGAN[IDR][0],F_Snd_Num[IDR],MPI_INT,IDR,tag,mpi_comm_level1,&stat); 
      MPI_Wait(&request,&stat);
    }
  }

  /*****************************************
       MPI:  S_Rcv_Num
  *****************************************/

  /* MPI_Barrier */
  MPI_Barrier(mpi_comm_level1);
  
  tag = 999;
  for (ID=0; ID<numprocs; ID++){
    if (ID!=0){
      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;
      MPI_Isend(&S_Rcv_Num[IDS],1,MPI_INT,IDS,tag,mpi_comm_level1,&request);
      MPI_Recv( &S_Snd_Num[IDR],1,MPI_INT,IDR,tag,mpi_comm_level1,&stat); 
      MPI_Wait(&request,&stat);
    }
  }

  Snd_SGAN = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Snd_SGAN[ID] = (int*)malloc(sizeof(int)*S_Snd_Num[ID]);
  }

  for (ID=0; ID<numprocs; ID++){
    if (ID!=0){
      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;
      MPI_Isend(&Rcv_SGAN[IDS][0],S_Rcv_Num[IDS],MPI_INT,IDS,tag,mpi_comm_level1,&request);
      MPI_Recv( &Snd_SGAN[IDR][0],S_Snd_Num[IDR],MPI_INT,IDR,tag,mpi_comm_level1,&stat); 
      MPI_Wait(&request,&stat);
    }
  }

  /********************************************************
    allocation of arrays:
  
    int Snd_MAN[numprocs][F_Snd_Num[ID]+S_Snd_Num[ID]+1]
    int Snd_GAN[numprocs][F_Snd_Num[ID]+S_Snd_Num[ID]+1]
  *********************************************************/

  if (alloc_first[11]==0){
    for (ID=0; ID<numprocs; ID++){
      free(Snd_MAN[ID]); 
    }
    free(Snd_MAN);

    for (ID=0; ID<numprocs; ID++){
      free(Snd_GAN[ID]); 
    }
    free(Snd_GAN); 
  }  

  Snd_MAN = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Snd_MAN[ID] = (int*)malloc(sizeof(int)*(F_Snd_Num[ID]+S_Snd_Num[ID]+1));
  }

  Snd_GAN = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Snd_GAN[ID] = (int*)malloc(sizeof(int)*(F_Snd_Num[ID]+S_Snd_Num[ID]+1));
  }

  alloc_first[11] = 0;

  /************************************************
      find data structures to send informations
      related to FNAN from myid to the other IDs 
  *************************************************/

  for (ID=0; ID<numprocs; ID++){

    Num = 0;

    for (i=0; i<F_Snd_Num[ID]; i++){ 

      for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

	Gc_AN = M2G[Mc_AN];

	if (Gc_AN==Snd_FGAN[ID][i]){

	  Snd_MAN[ID][Num] = Mc_AN;
	  Snd_GAN[ID][Num] = Gc_AN;
	  Num++;
	}      
      }
    }
  }

  /************************************************
      find data structures to send informations
      related to SNAN from myid to the other IDs 
  *************************************************/

  for (ID=0; ID<numprocs; ID++){

    Num = 0;

    for (i=0; i<S_Snd_Num[ID]; i++){ 

      for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

	Gc_AN = M2G[Mc_AN];

	if (Gc_AN==Snd_SGAN[ID][i]){

	  Snd_MAN[ID][F_Snd_Num[ID]+Num] = Mc_AN;
	  Snd_GAN[ID][F_Snd_Num[ID]+Num] = Gc_AN;
	  Num++;
	}      
      }
    }
  }

  /************************************************
   MPI:

     Snd_GAN
     Rcv_GAN
  *************************************************/

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){
      tag = 999;
      MPI_Isend(&Snd_GAN[IDS][0], F_Snd_Num[IDS]+S_Snd_Num[IDS],
                MPI_INT,IDS, tag, mpi_comm_level1, &request);
      MPI_Recv(&Rcv_GAN[IDR][0],  F_Rcv_Num[IDR]+S_Rcv_Num[IDR],
                MPI_INT, IDR, tag, mpi_comm_level1, &stat);
      MPI_Wait(&request,&stat);
    }
  }

  /************************************************
          setting of F_TopMAN and S_TopMAN

    F_TopMAN and S_TopMAN give the first intermediate
    atom number in atoms sent from ID in the size of
    F_Rcv_Num[ID] and F_Rcv_Num[ID] + S_Rcv_Num[ID],
    respectively.
  *************************************************/

  Num = Matomnum + 1;
  for (ID=0; ID<numprocs; ID++){
    if (F_Rcv_Num[ID]!=0 && ID!=myid){
      F_TopMAN[ID] = Num;
      Num = Num + F_Rcv_Num[ID];
    }
  }
  
  Num = Matomnum + 1;
  for (ID=0; ID<numprocs; ID++){
    if ((F_Rcv_Num[ID]!=0 || S_Rcv_Num[ID]!=0) && ID!=myid){
      S_TopMAN[ID] = Num;
      Num = Num + F_Rcv_Num[ID] + S_Rcv_Num[ID];
    }
  }

  /************************************************
       MatomnumF = the sum of F_Rcv_Num[ID]
       MatomnumS = the sum of S_Rcv_Num[ID]
  *************************************************/

  MatomnumF = 0;
  for (ID=0; ID<numprocs; ID++){
    if (ID!=myid) MatomnumF += F_Rcv_Num[ID];
  }

  MatomnumS = 0;
  for (ID=0; ID<numprocs; ID++){
    if (ID!=myid) MatomnumS += S_Rcv_Num[ID];
  }

  /************************************************
       allocation of arrays:
         
          F_M2G
          S_M2G
  *************************************************/

  F_M2G = (int*)malloc(sizeof(int)*(Matomnum+MatomnumF+1));
  S_M2G = (int*)malloc(sizeof(int)*(Matomnum+MatomnumF+MatomnumS+1));
  alloc_first[13] = 0;

  /************************************************
           setting of F_G2M and S_G2M

    F_G2M and S_G2M give a conversion from the
    global atom number to the medium atom number
    for atoms sent from ID in the size of
    F_Rcv_Num[ID] and F_Rcv_Num[ID] + S_Rcv_Num[ID],
    respectively. 
  *************************************************/
  
  /* initialization of F_G2M and S_G2M */
  
  for (Gc_AN=0; Gc_AN<=atomnum; Gc_AN++){
    F_G2M[Gc_AN] = -1;
    S_G2M[Gc_AN] = -1;
  } 
  
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    F_G2M[Gc_AN] = Mc_AN;
    S_G2M[Gc_AN] = Mc_AN;
    F_M2G[Mc_AN] = Gc_AN;
    S_M2G[Mc_AN] = Gc_AN;
  }

  for (ID=0; ID<numprocs; ID++){
    if (ID!=myid && F_Rcv_Num[ID]!=0){
      for (Num=0; Num<F_Rcv_Num[ID]; Num++){
        Gc_AN = Rcv_GAN[ID][Num];
        F_G2M[Gc_AN] = F_TopMAN[ID] + Num;
        F_M2G[F_TopMAN[ID] + Num] = Gc_AN;  
      }
    }
  }

  for (ID=0; ID<numprocs; ID++){
    if (ID!=myid && (F_Rcv_Num[ID]!=0 || S_Rcv_Num[ID]!=0)){
      for (Num=0; Num<(F_Rcv_Num[ID]+S_Rcv_Num[ID]); Num++){

        Gc_AN = Rcv_GAN[ID][Num];
        S_G2M[Gc_AN] = S_TopMAN[ID] + Num;
        S_M2G[S_TopMAN[ID] + Num] = Gc_AN;
      }
    }
  }

  /************************************************
      analysis of data structures  
      for MPI communication of DS_VNA
  *************************************************/

  Num_Pro_Snd = (int*)malloc(sizeof(int)*numprocs);

  for (ID=0; ID<numprocs; ID++)  Num_Pro_Snd[ID] = 0;

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    for (Lh_AN=0; Lh_AN<=FNAN[Gc_AN]; Lh_AN++){
      Gh_AN = natn[Gc_AN][Lh_AN];
      ID1 = G2ID[Gh_AN];
      Num_Pro_Snd[ID1]++;
    }
  }

  if (alloc_first[22]==0){

    for (ID=0; ID<numprocs; ID++){
      free(Pro_Snd_GAtom[ID]);
    }
    free(Pro_Snd_GAtom);

    for (ID=0; ID<numprocs; ID++){
      free(Pro_Snd_MAtom[ID]);
    }
    free(Pro_Snd_MAtom);

    for (ID=0; ID<numprocs; ID++){
      free(Pro_Snd_LAtom[ID]);
    }
    free(Pro_Snd_LAtom);

    for (ID=0; ID<numprocs; ID++){
      free(Pro_Snd_LAtom2[ID]);
    }
    free(Pro_Snd_LAtom2);
  }

  Pro_Snd_GAtom = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Pro_Snd_GAtom[ID] = (int*)malloc(sizeof(int)*(Num_Pro_Snd[ID]+1));
    for (i=0; i<(Num_Pro_Snd[ID]+1); i++){
      Pro_Snd_GAtom[ID][i] = 0;
    }   
  }

  Pro_Snd_MAtom = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Pro_Snd_MAtom[ID] = (int*)malloc(sizeof(int)*(Num_Pro_Snd[ID]+1));
    for (i=0; i<(Num_Pro_Snd[ID]+1); i++){
      Pro_Snd_MAtom[ID][i] = 0;
    }   
  }

  Pro_Snd_LAtom = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Pro_Snd_LAtom[ID] = (int*)malloc(sizeof(int)*(Num_Pro_Snd[ID]+1));
    for (i=0; i<(Num_Pro_Snd[ID]+1); i++){
      Pro_Snd_LAtom[ID][i] = 0;
    }   
  }

  Pro_Snd_LAtom2 = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Pro_Snd_LAtom2[ID] = (int*)malloc(sizeof(int)*(Num_Pro_Snd[ID]+1));
    for (i=0; i<(Num_Pro_Snd[ID]+1); i++){
      Pro_Snd_LAtom2[ID][i] = 0;
    }   
  }

  alloc_first[22] = 0;

  for (ID=0; ID<numprocs; ID++)  Num_Pro_Snd[ID] = 0;
  
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    for (Lh_AN=0; Lh_AN<=FNAN[Gc_AN]; Lh_AN++){
    
      Gh_AN = natn[Gc_AN][Lh_AN];
      ID1 = G2ID[Gh_AN];
    
      Pro_Snd_GAtom[ID1][Num_Pro_Snd[ID1]] = Gh_AN;
      Pro_Snd_MAtom[ID1][Num_Pro_Snd[ID1]] = Mc_AN;
      Pro_Snd_LAtom[ID1][Num_Pro_Snd[ID1]] = Lh_AN;
      
      Num_Pro_Snd[ID1]++;
    }
  }

  /* quick sorting of Pro_Snd_GAtom */ 

  for (ID=0; ID<numprocs; ID++){ 
    qsort_int3(Num_Pro_Snd[ID], Pro_Snd_GAtom[ID], Pro_Snd_MAtom[ID], Pro_Snd_LAtom[ID]);
  }

  /* setting of Pro_Snd_LAtom2 */

  for (ID=0; ID<numprocs; ID++){ 

    for (i=0; i<Num_Pro_Snd[ID]; i++){

       Gh_AN = Pro_Snd_GAtom[ID][i];
       Mc_AN = Pro_Snd_MAtom[ID][i];
       Lh_AN = Pro_Snd_LAtom[ID][i];

       Gc_AN = M2G[Mc_AN];
       Rn = ncn[Gc_AN][Lh_AN];

       m1 =-atv_ijk[Rn][1];
       m2 =-atv_ijk[Rn][2];
       m3 =-atv_ijk[Rn][3];

       j = 0;
       po = 0;

       do{ 

         Gj_AN = natn[Gh_AN][j];
         Rn2 = ncn[Gh_AN][j];

         n1 = atv_ijk[Rn2][1];
         n2 = atv_ijk[Rn2][2];
         n3 = atv_ijk[Rn2][3];
              
         if (m1==n1 && m2==n2 && m3==n3 &&  Gj_AN==Gc_AN){
           Pro_Snd_LAtom2[ID][i] = j;
	   po = 1;
         }

         j++;

       } while (po==0);

    }
  }

  /*********************************
   freeing of arrays:
  *********************************/

  free(Num_Pro_Snd);

  for (ID=0; ID<numprocs; ID++){
    free(Snd_SGAN[ID]);
  }
  free(Snd_SGAN);

  for (ID=0; ID<numprocs; ID++){
    free(Snd_FGAN[ID]);
  }
  free(Snd_FGAN);

  for (ID=0; ID<numprocs; ID++){
    free(Rcv_SGAN[ID]);
  }
  free(Rcv_SGAN);

  for (ID=0; ID<numprocs; ID++){
    free(Rcv_FGAN[ID]);
  }
  free(Rcv_FGAN);

  for (Gc_AN=0; Gc_AN<(atomnum+1); Gc_AN++){
    free(po_ID[Gc_AN]);
  }
  free(po_ID);

} 








void allocate_grids2atoms(int MD_iter)
{
  int i,j,po,Mc_AN,Gc_AN,wan;
  int MinC,MaxC,num,n1,gn[10];
  int *Cell_IDtmp;
  double Cxyz[4],b[4],c[4],v[4];
  double coef,rcut,Vec0,Vec1;
  double MinV,MaxV;
  double avg_Ngrid1,avg_Ngrid2;
  double atom_pos,a_unit;
  int cell,FNAN2;
  int *exist_switch;
  int *natn2;
  int *Num_Rcv_FNAN2, *Num_Snd_FNAN2;
  int **Rcv_FNAN2, **Snd_FNAN2;
  int *TopMAN2;
  int **GLAtom;
  int **GLCell;
  int **Rcv_FNAN2_At;
  int **Rcv_FNAN2_Nc;
  int n2,n3,nn1,N,N3[4],GNc,Nc;
  int myid,numprocs,ID,IDS,IDR,Geta,tag=999;

  MPI_Status stat;
  MPI_Request request;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /******************************************************
    find the unit vector perpendicular to the bc-plane
  ******************************************************/

  b[1] = tv[2][1];
  b[2] = tv[2][2];
  b[3] = tv[2][3];

  c[1] = tv[3][1];
  c[2] = tv[3][2];
  c[3] = tv[3][3];

  Cross_Product(b,c,v);
  coef = 1.0/sqrt(fabs( Dot_Product(v,v) ));

  v[1] = coef*v[1];
  v[2] = coef*v[2];
  v[3] = coef*v[3];

  /******************************************************
    find the minimun and maximum grid numbers of a-axis
  ******************************************************/

  MinV =  1.0e+10;
  MaxV = -1.0e+10;

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    wan  = WhatSpecies[Gc_AN];
    rcut = Spe_Atom_Cut1[wan] + 0.5;

    Cxyz[1] = Gxyz[Gc_AN][1] + rcut*v[1] - Grid_Origin[1];
    Cxyz[2] = Gxyz[Gc_AN][2] + rcut*v[2] - Grid_Origin[2];
    Cxyz[3] = Gxyz[Gc_AN][3] + rcut*v[3] - Grid_Origin[3];

    Vec0 = Dot_Product(Cxyz,rgtv[1])*0.5/PI;

    Cxyz[1] = Gxyz[Gc_AN][1] - rcut*v[1] - Grid_Origin[1];
    Cxyz[2] = Gxyz[Gc_AN][2] - rcut*v[2] - Grid_Origin[2];
    Cxyz[3] = Gxyz[Gc_AN][3] - rcut*v[3] - Grid_Origin[3];

    Vec1 = Dot_Product(Cxyz,rgtv[1])*0.5/PI;

    if (Vec0<MinV) MinV = Vec0;
    if (Vec1<MinV) MinV = Vec1;
    if (MaxV<Vec0) MaxV = Vec0;   
    if (MaxV<Vec1) MaxV = Vec1;
  }

  MinC = (int)MinV - 3;  /* buffer for GGA */
  MaxC = (int)MaxV + 3;  /* buffer for GGA */


  /*
  printf("MinC=%2d MaxC=%2d\n",MinC,MaxC);
  */

  /******************************************************
    set cells of a-axis overlapping to atoms in ID
  ******************************************************/

  if (alloc_first[15]==0){
    free(My_Cell0);
    free(Cell_ID0); 
    free(edge_block);
  }
  My_Cell0 = (int*)malloc(sizeof(int)*Ngrid1);
  Cell_ID0 = (int*)malloc(sizeof(int)*Ngrid1);
  edge_block = (int*)malloc(sizeof(int)*Ngrid1);
  alloc_first[15] = 0;

  for (i=0; i<Ngrid1; i++){
    edge_block[i] =  0;
    My_Cell0[i]   = -1;
  }

  num = -1;
  for (i=MinC; i<=MaxC; i++){
    if (0<=i) j = i%Ngrid1;
    else      j = (Ngrid1 + i%Ngrid1)%Ngrid1;

    if (My_Cell0[j]==-1){
      num++;
      My_Cell0[j] = num;
    }
  }

  Num_Cells0 = num + 1;

  /*
  for (i=0; i<Ngrid1; i++){
    printf("i=%3d  My_Cell0=%3d\n",i,My_Cell0[i]);
  }
  exit(0);
  */

  /******************************************************
                      set edge_block
  ******************************************************/

  for (i=0; i<Ngrid1; i++){

    gn[0] = i - 2;
    gn[1] = i - 1;
    gn[2] = i + 1;
    gn[3] = i + 2;
 
    po = 0;
    for (j=0; j<=3; j++){
      if (0<=gn[j]) gn[j] = gn[j]%Ngrid1;
      else          gn[j] = (Ngrid1 + gn[j]%Ngrid1)%Ngrid1;
      if (My_Cell0[ gn[j] ]==-1) po = 1; 
    }

    if ( po==1 ) edge_block[i] = 1;
  }

  /******************************************************
       set cells on a-axis overlapping to atoms in ID
  ******************************************************/

  if (alloc_first[14]==0){
    free(My_Cell1);
  }
  My_Cell1 = (int*)malloc(sizeof(int)*Num_Cells0);
  alloc_first[14] = 0;

  for (i=0; i<Ngrid1; i++){
    if (My_Cell0[i]!=-1) My_Cell1[My_Cell0[i]] = i;
  }

  My_NumGrid1 = Num_Cells0*Ngrid2*Ngrid3;

  /******************************************************
   ******************************************************
              set informations of grids for
               solving Poisson's equation 
   ******************************************************
  ******************************************************/

  /* Start_Grid1 and End_Grid1 
     Start_Grid2 and End_Grid2 */

  avg_Ngrid1 = (double)Ngrid1/(double)numprocs;  
  avg_Ngrid2 = (double)Ngrid2/(double)numprocs;  

  for (ID=0; ID<numprocs; ID++){

    if (ID==0){
       Start_Grid1[ID] = 0;
       Start_Grid2[ID] = 0;
    }
    else{ 
       Start_Grid1[ID] = End_Grid1[ID-1] + 1;
       Start_Grid2[ID] = End_Grid2[ID-1] + 1;
    }

    End_Grid1[ID] = (int)floor(avg_Ngrid1*(double)(ID+1)-0.01);
    End_Grid2[ID] = (int)floor(avg_Ngrid2*(double)(ID+1)-0.01);
  
    if (ID==(numprocs-1)){
      End_Grid1[ID] = Ngrid1 - 1;
      End_Grid2[ID] = Ngrid2 - 1;
    }
  }

  My_NGrid1_Poisson = End_Grid1[myid] - Start_Grid1[myid] + 1;
  My_NGrid2_Poisson = End_Grid2[myid] - Start_Grid2[myid] + 1;

  /* Cell_ID0 */
 
  for (i=0; i<Ngrid1; i++){
    Cell_ID0[i] = -1;
  }

  Geta = 1000000;
  for (i=0; i<Num_Cells0; i++){
    j = My_Cell1[i];
    if (edge_block[j]==0){
      if (Start_Grid1[myid]<=j && j<=End_Grid1[myid])
        Cell_ID0[j] = myid + Geta;   
      else         
        Cell_ID0[j] = myid;   
    }
  }

  Cell_IDtmp = (int*)malloc(sizeof(int)*Ngrid1);
  for (i=0; i<Ngrid1; i++){
    MPI_Reduce(&Cell_ID0[i], &Cell_IDtmp[i], 1, MPI_INT, MPI_MAX, Host_ID, mpi_comm_level1);

    if (myid==Host_ID){
      if (Geta<=Cell_IDtmp[i])
        Cell_ID0[i] = Cell_IDtmp[i] - Geta;  
      else
        Cell_ID0[i] = Cell_IDtmp[i];
    }
  }
  free(Cell_IDtmp);

  MPI_Bcast(&Cell_ID0[0], Ngrid1, MPI_INT, Host_ID, mpi_comm_level1);

  /* informations for sending and recieving */

  /************************************
  allocation of arrays

  Num_Rcv_Grid1[numprocs]
  Num_Snd_Grid1[numprocs]
  Rcv_Grid1[numprocs][Num_Rcv_Grid1[ID]]
  Snd_Grid1[numprocs][Num_Rcv_Grid1[ID]]
  ************************************/

  if (alloc_first[16]==0){

    free(Num_Rcv_Grid1);
    free(Num_Snd_Grid1);

    for (ID=0; ID<numprocs; ID++){
      free(Rcv_Grid1[ID]);
    }
    free(Rcv_Grid1);

    for (ID=0; ID<numprocs; ID++){
      free(Snd_Grid1[ID]);
    }
    free(Snd_Grid1);
  }

  Num_Rcv_Grid1 = (int*)malloc(sizeof(int)*numprocs);
  Num_Snd_Grid1 = (int*)malloc(sizeof(int)*numprocs);

  for (ID=0; ID<numprocs; ID++) Num_Rcv_Grid1[ID] = 0;
  for (i=Start_Grid1[myid]; i<=End_Grid1[myid]; i++){
    if (Cell_ID0[i]!=myid && Cell_ID0[i]!=-1){
      ID = Cell_ID0[i];
      Num_Rcv_Grid1[ID]++;
    }
  }
  
  /* allocation of Rcv_Grid1 */
  Rcv_Grid1 = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Rcv_Grid1[ID] = (int*)malloc(sizeof(int)*Num_Rcv_Grid1[ID]);
  }

  /* Rcv_Grid1 */   
  for (ID=0; ID<numprocs; ID++) Num_Rcv_Grid1[ID] = 0;
  for (i=Start_Grid1[myid]; i<=End_Grid1[myid]; i++){
    if (Cell_ID0[i]!=myid && Cell_ID0[i]!=-1){
      ID = Cell_ID0[i];
      Rcv_Grid1[ID][Num_Rcv_Grid1[ID]] = i;
      Num_Rcv_Grid1[ID]++;
    }    
  }

  /* MPI, Num_Rcv_Grid1 */
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){
      tag = 999;
      MPI_Isend(&Num_Rcv_Grid1[IDS], 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
      MPI_Recv( &Num_Snd_Grid1[IDR], 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
      MPI_Wait(&request,&stat);
    }
    else{
      Num_Snd_Grid1[IDR] = 0;
    }     
  }

  /* allocation of Snd_Grid1 */
  Snd_Grid1 = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Snd_Grid1[ID] = (int*)malloc(sizeof(int)*Num_Snd_Grid1[ID]);
  }
  alloc_first[16] = 0;

  /* MPI, Rcv_Grid1 */
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){
      tag = 999;
      MPI_Isend(&Rcv_Grid1[IDS][0], Num_Rcv_Grid1[IDS], MPI_INT, IDS, tag,
                mpi_comm_level1, &request);
      MPI_Recv(&Snd_Grid1[IDR][0], Num_Snd_Grid1[IDR], MPI_INT, IDR, tag,
               mpi_comm_level1, &stat);
      MPI_Wait(&request,&stat);
    }
  }

  /******************************************************
   ******************************************************
             set informations for converting
           from Poisson's grid to atom's grid
   ******************************************************
  ******************************************************/

  /*******************************************
   allocation of arrays

   Num_IRcv_Grid1[numprocs]
   Num_ISnd_Grid1[numprocs]
   IRcv_Grid1[numprocs][Num_IRcv_Grid1[ID]]
   ISnd_Grid1[numprocs][Num_IRcv_Grid1[ID]]
  ********************************************/

  if (alloc_first[17]==0){
    free(Num_IRcv_Grid1);
    free(Num_ISnd_Grid1);
    for (ID=0; ID<numprocs; ID++){
      free(IRcv_Grid1[ID]);
    }
    free(IRcv_Grid1);

    for (ID=0; ID<numprocs; ID++){
      free(ISnd_Grid1[ID]);
    }
    free(ISnd_Grid1);
  }

  /* allocation of Num_IRcv_Grid1 and Num_ISnd_Grid1 */ 
  Num_IRcv_Grid1 = (int*)malloc(sizeof(int)*numprocs);
  Num_ISnd_Grid1 = (int*)malloc(sizeof(int)*numprocs);

  /* Cell_IDtmp */
  Cell_IDtmp = (int*)malloc(sizeof(int)*Ngrid1);
  for (ID=0; ID<numprocs; ID++){
    for (n1=Start_Grid1[ID]; n1<=End_Grid1[ID]; n1++){
      Cell_IDtmp[n1] = ID;
    }  
  }

  /* Num_IRcv_Grid1 */
  for (ID=0; ID<numprocs; ID++) Num_IRcv_Grid1[ID] = 0;
  for (i=0; i<Num_Cells0; i++){
    n1 = My_Cell1[i];
    ID = Cell_IDtmp[n1];
    if (ID!=myid) Num_IRcv_Grid1[ID]++;
  }

  /* allocation of IRcv_Grid1 */
  IRcv_Grid1 = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    IRcv_Grid1[ID] = (int*)malloc(sizeof(int)*Num_IRcv_Grid1[ID]);
  }

  /* IRcv_Grid1 */   
  for (ID=0; ID<numprocs; ID++) Num_IRcv_Grid1[ID] = 0;
  for (i=0; i<Num_Cells0; i++){
    n1 = My_Cell1[i];
    ID = Cell_IDtmp[n1];
    if (ID!=myid){
      IRcv_Grid1[ID][Num_IRcv_Grid1[ID]] = n1;
      Num_IRcv_Grid1[ID]++;
    }
  }

  /* MPI, Num_IRcv_Grid1 */
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){
      tag = 999;
      MPI_Isend(&Num_IRcv_Grid1[IDS], 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
      MPI_Recv( &Num_ISnd_Grid1[IDR], 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
      MPI_Wait(&request,&stat);
    }
    else{
      Num_ISnd_Grid1[IDR] = 0;
    }
  }

  /* allocation of ISnd_Grid1 */
  ISnd_Grid1 = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    ISnd_Grid1[ID] = (int*)malloc(sizeof(int)*Num_ISnd_Grid1[ID]);
  }
  alloc_first[17] = 0;

  /* MPI, IRcv_Grid1 */
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){
      tag = 999;
      MPI_Isend(&IRcv_Grid1[IDS][0], Num_IRcv_Grid1[IDS], MPI_INT, IDS, tag,
                mpi_comm_level1, &request);
      MPI_Recv(&ISnd_Grid1[IDR][0], Num_ISnd_Grid1[IDR], MPI_INT, IDR, tag,
               mpi_comm_level1, &stat);
      MPI_Wait(&request,&stat);
    }
  }

  /* free */
  free(Cell_IDtmp);

  /******************************************************
              search FNAN2 atoms (not FNAN)
      which are used to calculate electron densities 
  ******************************************************/

  exist_switch = (int*)malloc(sizeof(int)*(atomnum+1));

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    exist_switch[Gc_AN] = 0;
  }

  for (Mc_AN=1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
    Gc_AN = F_M2G[Mc_AN];
    exist_switch[Gc_AN] = 1;
  }

  a_unit = sqrt(Dot_Product( gtv[1], gtv[1] ));

  num = 0;
  for (cell=-CpyCell; cell<=CpyCell; cell++){
    for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){

      if (exist_switch[Gc_AN]==0){
        wan = WhatSpecies[Gc_AN];
        rcut = Spe_Atom_Cut1[wan]/a_unit;
        atom_pos = (double)Ngrid1*(Cell_Gxyz[Gc_AN][1] + (double)cell);

        if (   (MinC<=atom_pos && atom_pos<=MaxC)
            ||
               (atom_pos<MinC && MinC<=(atom_pos+rcut))
            ||
               ((atom_pos-rcut)<=MaxC && MaxC<atom_pos)
           ){
            num++; 
            exist_switch[Gc_AN] = -1;
        }

      }
    }
  } 

  FNAN2 = num;

  natn2 = (int*)malloc(sizeof(int)*FNAN2);
  num = -1;  
  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    if (exist_switch[Gc_AN]<0){
      num++;
      natn2[num] = Gc_AN;
    } 
  }  

  free(exist_switch);

  /****************************************************
             data for sending and receiving
  ****************************************************/

  /* allocation of Num_Rcv_FNAN2 and  Num_Snd_FNAN2 */
  Num_Rcv_FNAN2 = (int*)malloc(sizeof(int)*numprocs);
  Num_Snd_FNAN2 = (int*)malloc(sizeof(int)*numprocs);

  /* count Num_Rcv_FNAN2 */
  for (ID=0; ID<numprocs; ID++) Num_Rcv_FNAN2[ID] = 0;
  for (i=0; i<FNAN2; i++){
    Gc_AN = natn2[i];
    ID = G2ID[Gc_AN];
    Num_Rcv_FNAN2[ID]++; 
  }

  /* allocation of Rcv_FNAN2 */
  Rcv_FNAN2 = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Rcv_FNAN2[ID] = (int*)malloc(sizeof(int)*Num_Rcv_FNAN2[ID]);
  }  

  /* set Rcv_FNAN2 */
  for (ID=0; ID<numprocs; ID++) Num_Rcv_FNAN2[ID] = 0;
  for (i=0; i<FNAN2; i++){
    Gc_AN = natn2[i];
    ID = G2ID[Gc_AN];
    Rcv_FNAN2[ID][Num_Rcv_FNAN2[ID]] = Gc_AN;
    Num_Rcv_FNAN2[ID]++; 
  }

  /* reset natn2 */
  num = 0;
  for (ID=0; ID<numprocs; ID++){
    if (Num_Rcv_FNAN2[ID]!=0){
      for (i=0; i<Num_Rcv_FNAN2[ID]; i++){
        natn2[num] = Rcv_FNAN2[ID][i];
        num++;
      }
    }
  }

  /* set TopMAN2 */
  TopMAN2 = (int*)malloc(sizeof(int)*numprocs);
  num = 0;
  for (ID=0; ID<numprocs; ID++){
    if (Num_Rcv_FNAN2[ID]!=0){
      TopMAN2[ID] = num;
      num += Num_Rcv_FNAN2[ID];
    }
  }

  /* MPI, Num_Rcv_FNAN2 */
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){
      tag = 999;
      MPI_Isend(&Num_Rcv_FNAN2[IDS], 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
      MPI_Recv( &Num_Snd_FNAN2[IDR], 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
      MPI_Wait(&request,&stat);
    }
    else{
      Num_Snd_FNAN2[IDR] = 0;
    }     
  }

  /* allocation of Snd_FNAN2 */
  Snd_FNAN2 = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Snd_FNAN2[ID] = (int*)malloc(sizeof(int)*Num_Snd_FNAN2[ID]);
  }  

  /* MPI, Rcv_FNAN2 to Snd_FNAN2 */
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){
      tag = 999;
      MPI_Isend(&Rcv_FNAN2[IDS][0], Num_Rcv_FNAN2[IDS], MPI_INT,
                IDS, tag, mpi_comm_level1, &request);
      MPI_Recv( &Snd_FNAN2[IDR][0], Num_Snd_FNAN2[IDR], MPI_INT,
                IDR, tag, mpi_comm_level1, &stat);
      MPI_Wait(&request,&stat);
    }
  }

  /**********************************
    allocation of arrays: 

    GLAtom
    GLCell
  **********************************/

  GLAtom = (int**)malloc(sizeof(int*)*FNAN2);
  for (i=0; i<FNAN2; i++){
    Gc_AN = natn2[i];
    GLAtom[i] = (int*)malloc(sizeof(int)*GridN_Atom[Gc_AN]);
  }

  GLCell = (int**)malloc(sizeof(int*)*FNAN2);
  for (i=0; i<FNAN2; i++){
    Gc_AN = natn2[i];
    GLCell[i] = (int*)malloc(sizeof(int)*GridN_Atom[Gc_AN]);
  }

  /* MPI, GLAtom and GLCell */

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){
      tag = 999;

      /****************************
                 GLAtom
      ****************************/

      /* Sending of data to IDS */
      if (Num_Snd_FNAN2[IDS]!=0){

        for (i=0; i<Num_Snd_FNAN2[IDS]; i++){
          Gc_AN = Snd_FNAN2[IDS][i];
          Mc_AN = F_G2M[Gc_AN];
          MPI_Isend(&GridListAtom[Mc_AN][0], GridN_Atom[Gc_AN], MPI_INT,
                    IDS, tag, mpi_comm_level1, &request);
	}
      }

      /* Receiving of data from IDR */
      if (Num_Rcv_FNAN2[IDR]!=0){

        for (i=0; i<Num_Rcv_FNAN2[IDR]; i++){
          Gc_AN = Rcv_FNAN2[IDR][i];
          Mc_AN = TopMAN2[IDR] + i;
          MPI_Recv(&GLAtom[Mc_AN][0], GridN_Atom[Gc_AN], MPI_INT,
                   IDR, tag, mpi_comm_level1, &stat);
	} 
      }

      if (Num_Snd_FNAN2[IDS]!=0) MPI_Wait(&request,&stat);

      /****************************
                 GLCell
      ****************************/

      /* Sending of data to IDS */
      if (Num_Snd_FNAN2[IDS]!=0){

        for (i=0; i<Num_Snd_FNAN2[IDS]; i++){
          Gc_AN = Snd_FNAN2[IDS][i];
          Mc_AN = F_G2M[Gc_AN];
          MPI_Isend(&CellListAtom[Mc_AN][0], GridN_Atom[Gc_AN], MPI_INT,
                    IDS, tag, mpi_comm_level1, &request);
	}
      }

      /* Receiving of data from IDR */
      if (Num_Rcv_FNAN2[IDR]!=0){

        for (i=0; i<Num_Rcv_FNAN2[IDR]; i++){
          Gc_AN = Rcv_FNAN2[IDR][i];
          Mc_AN = TopMAN2[IDR] + i;
          MPI_Recv(&GLCell[Mc_AN][0], GridN_Atom[Gc_AN], MPI_INT,
                   IDR, tag, mpi_comm_level1, &stat);
	} 
      }

      if (Num_Snd_FNAN2[IDS]!=0) MPI_Wait(&request,&stat);

    } 
  }

  /*****************************************************************
   use the following variables when data comunicatios
   in terms of FNAN2

   Num_Rcv_FNAN2_Grid
   Num_Snd_FNAN2_Grid
   Snd_FNAN2_At
   Snd_FNAN2_Nc
   Rcv_FNAN2_GA
   Rcv_FNAN2_MN
   Rcv_FNAN2_GRc
   FNAN2_Grid
   TopMAN2_Grid
  *****************************************************************/

  if (alloc_first[18]==0){
    free(Num_Rcv_FNAN2_Grid);
    free(Num_Snd_FNAN2_Grid);
  }

  Num_Rcv_FNAN2_Grid = (int*)malloc(sizeof(int)*numprocs);
  Num_Snd_FNAN2_Grid = (int*)malloc(sizeof(int)*numprocs);

  for (ID=0; ID<numprocs; ID++) Num_Rcv_FNAN2_Grid[ID] = 0;

  FNAN2_Grid = 0;
  for (ID=0; ID<numprocs; ID++){
    if (ID!=myid){
      if (Num_Rcv_FNAN2[ID]!=0){
        for (i=0; i<Num_Rcv_FNAN2[ID]; i++){
          Gc_AN = Rcv_FNAN2[ID][i];
          Mc_AN = TopMAN2[ID] + i;
          for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
            GNc = GLAtom[Mc_AN][Nc];
            GN2N(GNc,N3);
            n1 = N3[1];  
            n2 = N3[2];  
            n3 = N3[3];  
            nn1 = My_Cell0[n1];
            if (nn1!=-1){
              Num_Rcv_FNAN2_Grid[ID]++; 
              FNAN2_Grid++;
	    }
	  }
	}
      }
    }
  }

  /**********************************************
   allocation of local arrays

   Rcv_FNAN2_At, Rcv_FNAN2_Nc, and Rcv_FNAN2_MN
  **********************************************/

  if (alloc_first[18]==0){
    free(TopMAN2_Grid);
    free(Rcv_FNAN2_MN);
    free(Rcv_FNAN2_GA);
    free(Rcv_FNAN2_GRc);
  }

  Rcv_FNAN2_MN = (int*)malloc(sizeof(int)*FNAN2_Grid);
  Rcv_FNAN2_GA = (int*)malloc(sizeof(int)*FNAN2_Grid);
  Rcv_FNAN2_GRc = (int*)malloc(sizeof(int)*FNAN2_Grid);

  Rcv_FNAN2_At = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Rcv_FNAN2_At[ID] = (int*)malloc(sizeof(int)*Num_Rcv_FNAN2_Grid[ID]);
  }
  Rcv_FNAN2_Nc = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Rcv_FNAN2_Nc[ID] = (int*)malloc(sizeof(int)*Num_Rcv_FNAN2_Grid[ID]);
  }

  /* set TopMAN2_Grid */
  TopMAN2_Grid = (int*)malloc(sizeof(int)*numprocs);
  num = 0;
  for (ID=0; ID<numprocs; ID++){
    if (Num_Rcv_FNAN2_Grid[ID]!=0){
      TopMAN2_Grid[ID] = num;
      num += Num_Rcv_FNAN2_Grid[ID];
    }
  }

  for (ID=0; ID<numprocs; ID++) Num_Rcv_FNAN2_Grid[ID] = 0;
  FNAN2_Grid = 0;

  for (ID=0; ID<numprocs; ID++){

    if (ID!=myid){
      if (Num_Rcv_FNAN2[ID]!=0){

        for (i=0; i<Num_Rcv_FNAN2[ID]; i++){

          Gc_AN = Rcv_FNAN2[ID][i];
          Mc_AN = TopMAN2[ID] + i;

          for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){

            GNc = GLAtom[Mc_AN][Nc];
            GN2N(GNc,N3);
            n1 = N3[1];  
            n2 = N3[2];  
            n3 = N3[3];  
            nn1 = My_Cell0[n1];

            if (nn1!=-1){

              N = nn1*Ngrid2*Ngrid3 + n2*Ngrid3 + n3;

              Rcv_FNAN2_At[ID][Num_Rcv_FNAN2_Grid[ID]] = Gc_AN;
              Rcv_FNAN2_Nc[ID][Num_Rcv_FNAN2_Grid[ID]] = Nc;
              Rcv_FNAN2_MN[FNAN2_Grid] = N;
              Rcv_FNAN2_GA[FNAN2_Grid] = Gc_AN;

              /* for NEGF */
              Rcv_FNAN2_GRc[FNAN2_Grid] = GLCell[Mc_AN][Nc];

              Num_Rcv_FNAN2_Grid[ID]++; 
              FNAN2_Grid++;
	    }
	  }
	}
      }
    }
  }

  /* MPI, Num_Rcv_FNAN2_Grid */
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    tag = 999;
    if (ID!=0){
      MPI_Isend(&Num_Rcv_FNAN2_Grid[IDS], 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
      MPI_Recv( &Num_Snd_FNAN2_Grid[IDR], 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
      MPI_Wait(&request,&stat);
    }
    else{
      Num_Snd_FNAN2_Grid[IDR] = 0;
    }     
  }

  /* allocation of Snd_FNAN2_At */

  if (alloc_first[18]==0){
    for (ID=0; ID<numprocs; ID++){
      free(Snd_FNAN2_At[ID]);
    }  
    free(Snd_FNAN2_At);
  }

  Snd_FNAN2_At = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Snd_FNAN2_At[ID] = (int*)malloc(sizeof(int)*Num_Snd_FNAN2_Grid[ID]);
  }  

  /* allocation of Snd_FNAN2_Nc */

  if (alloc_first[18]==0){
    for (ID=0; ID<numprocs; ID++){
      free(Snd_FNAN2_Nc[ID]);
    }  
    free(Snd_FNAN2_Nc);
  }

  Snd_FNAN2_Nc = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Snd_FNAN2_Nc[ID] = (int*)malloc(sizeof(int)*Num_Snd_FNAN2_Grid[ID]);
  }

  /* MPI, Rcv_FNAN2_At to Snd_FNAN2_At */
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){
      tag = 999;
      MPI_Isend(&Rcv_FNAN2_At[IDS][0], Num_Rcv_FNAN2_Grid[IDS], MPI_INT,
                IDS, tag, mpi_comm_level1, &request);
      MPI_Recv( &Snd_FNAN2_At[IDR][0], Num_Snd_FNAN2_Grid[IDR], MPI_INT,
                IDR, tag, mpi_comm_level1, &stat);
      MPI_Wait(&request,&stat);
    }
  }

  /* MPI, Rcv_FNAN2_Nc to Snd_FNAN2_Nc */
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){
      tag = 999;
      MPI_Isend(&Rcv_FNAN2_Nc[IDS][0], Num_Rcv_FNAN2_Grid[IDS], MPI_INT,
                IDS, tag, mpi_comm_level1, &request);
      MPI_Recv( &Snd_FNAN2_Nc[IDR][0], Num_Snd_FNAN2_Grid[IDR], MPI_INT,
                IDR, tag, mpi_comm_level1, &stat);
      MPI_Wait(&request,&stat);
    }
  }

  /* alloc_first[18] = 0; */
  alloc_first[18] = 0;

  /* free */

  free(natn2);
 
  for (i=0; i<FNAN2; i++){
    free(GLAtom[i]);
  }
  free(GLAtom);

  for (i=0; i<FNAN2; i++){
    free(GLCell[i]);
  }
  free(GLCell);

  free(Num_Rcv_FNAN2);
  free(Num_Snd_FNAN2);

  for (ID=0; ID<numprocs; ID++){
    free(Rcv_FNAN2[ID]);
  }  
  free(Rcv_FNAN2);

  for (ID=0; ID<numprocs; ID++){
    free(Snd_FNAN2[ID]);
  }  
  free(Snd_FNAN2);

  free(TopMAN2);

  for (ID=0; ID<numprocs; ID++){
    free(Rcv_FNAN2_Nc[ID]);
  }
  free(Rcv_FNAN2_Nc);

  for (ID=0; ID<numprocs; ID++){
    free(Rcv_FNAN2_At[ID]);
  }  
  free(Rcv_FNAN2_At);
}

