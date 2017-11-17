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
#include "mpi.h"
#include "omp.h"
#include "tran_prototypes.h"

#define  measure_time   0


static void Fixed_FNAN_SNAN();
static int Set_Periodic(int CpyN, int Allocate_switch);
static void Free_truncation(int CpyN, int TN, int Free_switch);
static void free_arrays_truncation0();
static void Trn_System(int MD_iter, int CpyCell, int TCpyCell);
static void Estimate_Trn_System(int CpyCell, int TCpyCell);
static void Check_System();
static void Set_RMI();
static void Output_Connectivity(FILE *fp);
static void UCell_Box(int MD_iter, int estimate_switch, int CpyCell);
static void Set_Inf_SndRcv();
static void Construct_MPI_Data_Structure_Grid();

int TFNAN,TFNAN2,TSNAN,TSNAN2;




double truncation(int MD_iter,int UCell_flag)
{
  static int firsttime=1;
  int i,j,k,m,l,ct_AN,h_AN,Gh_AN,Mc_AN,Gc_AN,s1,s2;
  int tno,tno0,tno1,Cwan,Hwan,N,so,nc,ns,spin;
  int num,wan,n2,wanA,Gi,Max_Num_Cells0;
  int NO1,Mh_AN,Rnh;
  int vsize,Anum,Bnum,NUM,p,MAnum,fan,csize;
  int size_Orbs_Grid,size_COrbs_Grid;
  int size_Orbs_Grid_FNAN;
  int size_H0,size_CntH0,size_H,size_CntH;
  int size_HNL;
  int size_iHNL,size_iHNL0,size_iCntHNL;
  int size_DS_NL,size_CntDS_NL;
  int size_NumOLG,size_OLP,size_CntOLP;
  int size_OLP_L;
  int size_HVNA,size_DS_VNA,size_CntDS_VNA;
  int size_HVNA2,size_CntHVNA2;
  int size_HVNA3,size_CntHVNA3;
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
  double time0,time1,time2,time3,time4,time5,time6,time7;
  double time8,time9,time10,time11,time12,time13,time14;
  double time15,time16,time17,time18,time19,time20,time21;
  double stime,etime;

  double My_KryDH,My_KryDS;
  double KryDH,KryDS;

  char file_TRN[YOUSO10] = ".TRN";
  double TStime,TEtime;
  FILE *fp_TRN;
  char buf[fp_bsize];          /* setvbuf */

  /* MPI */

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  time0 = 0.0; time1 = 0.0; time2 = 0.0; time3 = 0.0; time4 = 0.0;
  time5 = 0.0; time6 = 0.0; time7 = 0.0; time8 = 0.0; time9 = 0.0;
  time10= 0.0; time11= 0.0; time12= 0.0; time13= 0.0; time14= 0.0;
  time15= 0.0; time16= 0.0; time17= 0.0; time18= 0.0; time19= 0.0;

  if (myid==Host_ID && 0<level_stdout){
    printf("\n*******************************************************\n"); 
    printf("        Analysis of neigbors and setting of grids        \n");
    printf("*******************************************************\n\n"); 
  }

  dtime(&TStime); 

  /****************************************************
     freeing of dynamically allocated arrays using
     FNAN and SNAN at the previous MD step
  ****************************************************/

  if (measure_time) dtime(&stime); 

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    ID = G2ID[Gc_AN];
    MPI_Bcast(&time_per_atom[Gc_AN], 1, MPI_DOUBLE, ID, mpi_comm_level1);
  }

  free_arrays_truncation0();

  if (measure_time){
    dtime(&etime); 
    time0 = etime - stime;
  }

  /****************************************************
            allocation of atoms to processors
  ****************************************************/

  if (measure_time) dtime(&stime); 

  if (2<=MD_iter){

    /*****************************
      last input
      0: atomnum, 
      1: elapsed time 
    *****************************/    

    if (Solver==5 || Solver==8)  /* DC or Krylov */
      Set_Allocate_Atom2CPU(MD_iter,0,1);  /* 1: elapsed time */
    else 
      Set_Allocate_Atom2CPU(MD_iter,0,0);  /* 0: atomnum */
  }

  if (measure_time){
    dtime(&etime); 
    time1 = etime - stime;
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
            and
        Generation_ATV(CpyN);

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

      if (measure_time) dtime(&stime); 

      TCpyCell = Set_Periodic(CpyCell,0);

      if (measure_time){
	dtime(&etime); 
	time2 += etime - stime;
      }

      /**********************************************
        find Max_FSNAN by the physical truncation
          for allocation of natn, ncn, and Dis 
      **********************************************/

      if (measure_time) dtime(&stime); 

      Estimate_Trn_System(CpyCell,TCpyCell);

      if (measure_time){
	dtime(&etime); 
	time3 += etime - stime;
      }

      /**********************************************
              allocation of natn, ncn, and Dis 
      **********************************************/

      if (measure_time) dtime(&stime); 

      Allocate_Arrays(3);

      if (measure_time){
	dtime(&etime); 
	time4 += etime - stime;
      }

      /**********************
       find TFNAN and TSNAN
      **********************/

      TFNAN2 = TFNAN;
      TSNAN2 = TSNAN;

      if (measure_time) dtime(&stime); 

      Trn_System(MD_iter,CpyCell,TCpyCell);

      if (measure_time){
	dtime(&etime); 
	time5 += etime - stime;
      }

      if ( TFNAN==TFNAN2 && TSNAN==TSNAN2 && (Solver==1 || Solver==5 || Solver==6 || Solver==8) ) po++;
      else if (TFNAN==TFNAN2 && (Solver==2 || Solver==3 || Solver==4 || Solver==7 || Solver==9) ) po++;
      else if (CellNN_flag==1)                                                                    po++;

      /**********************************************
         freeing of arrays which are allocated in
                 Set_Periodic(CpyCell,0)
      **********************************************/

      if (measure_time) dtime(&stime); 

      Free_truncation(CpyCell,TCpyCell,0);

      if (measure_time){
	dtime(&etime); 
	time6 += etime - stime;
      }

    } while (po==0);

    List_YOUSO[4] = CpyCell;

    if (measure_time) dtime(&stime); 

    TCpyCell = Set_Periodic(CpyCell,0);

    if (measure_time){
      dtime(&etime); 
      time2 += etime - stime;
    }

    if (measure_time) dtime(&stime); 

    Estimate_Trn_System(CpyCell,TCpyCell);

    if (measure_time){
      dtime(&etime); 
      time3 += etime - stime;
    }

    if (measure_time) dtime(&stime); 

    Allocate_Arrays(3);

    if (measure_time){
      dtime(&etime); 
      time7 += etime - stime;
    }

    if (measure_time) dtime(&stime); 

    Trn_System(MD_iter,CpyCell,TCpyCell);

    if (measure_time){
      dtime(&etime); 
      time5 += etime - stime;
    }

    if (measure_time) dtime(&stime); 

    Set_Inf_SndRcv();

    if (measure_time){
      dtime(&etime); 
      time8 += etime - stime;
    }

    if (measure_time) dtime(&stime); 

    Set_RMI();

    if (measure_time){
      dtime(&etime); 
      time9 += etime - stime;
    }

  } /* if (MD_iter==1) */ 
    
  else{

    TCpyCell = Set_Periodic(CpyCell,0);
    Estimate_Trn_System(CpyCell,TCpyCell);
    Allocate_Arrays(3);
    Trn_System(MD_iter,CpyCell,TCpyCell);

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

  if (measure_time) dtime(&stime); 

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

  if (measure_time){
    dtime(&etime); 
    time10 += etime - stime;
  }

  /****************************************************
                  check the system type
  ****************************************************/

  if (measure_time) dtime(&stime); 

  if (MD_iter==1) Check_System();

  if (measure_time){
    dtime(&etime); 
    time11 += etime - stime;
  }

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

    if (measure_time) dtime(&stime); 

    UCell_Box(MD_iter,1,CpyCell);

    if (measure_time){
      dtime(&etime); 
      time12 += etime - stime;
    }

    if (measure_time) dtime(&stime); 

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

    if (measure_time){
      dtime(&etime); 
      time13 += etime - stime;
    }

    /*************************************
      find 
            Max_NumOLG
    *************************************/

    Max_NumOLG = 0;
    if (2<=level_stdout) printf("\n***** UCell_Box(MD_iter,2,CpyCell) *****\n");

    if (measure_time) dtime(&stime); 

    UCell_Box(MD_iter,2,CpyCell);

    if (measure_time){
      dtime(&etime); 
      time14 += etime - stime;
    }

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
                 GListTAtoms1
                 GListTAtoms2
    *************************************/

    if (myid==Host_ID && 0<level_stdout)  printf("<UCell_Box> Info. of cutoff energy and num. of grids\n");

    if (measure_time) dtime(&stime); 

    UCell_Box(MD_iter,0,CpyCell);

    if (measure_time){
      dtime(&etime); 
      time15 += etime - stime;
    }

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

  if (measure_time) dtime(&stime); 

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
      else if ( (Hub_U_switch==0 || Hub_U_occupation!=1) && 0<k && Matomnum<Mc_AN){
        Gc_AN = S_M2G[Mc_AN];
        Cwan = WhatSpecies[Gc_AN];
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
        else if ( (Hub_U_switch==0 || Hub_U_occupation!=1) && 0<k && Matomnum<Mc_AN){
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
      DS_NL[so][k] = (double****)malloc(sizeof(double***)*(Matomnum+2)); 
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

	DS_NL[so][k][Mc_AN] = (double***)malloc(sizeof(double**)*(fan+1)); 
	for (h_AN=0; h_AN<(fan+1); h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
  	  else if ( (Matomnum+1)<=Mc_AN ){
	    tno1 = List_YOUSO[20];  
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
	CntDS_NL[so][k] = (double****)malloc(sizeof(double***)*(Matomnum+2)); 
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
	    tno0 = Spe_Total_CNO[Cwan];  
            fan = FNAN[Gc_AN];
	  }    

	  CntDS_NL[so][k][Mc_AN] = (double***)malloc(sizeof(double**)*(fan+1)); 
	  for (h_AN=0; h_AN<(fan+1); h_AN++){

	    if (Mc_AN==0){
	      tno1 = 1;  
	    }
  	    else if ( (Matomnum+1)<=Mc_AN ){
	      tno1 = List_YOUSO[20];  
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
      HVNA2[k] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
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

    /* HVNA3 */  

    size_HVNA3 = 0;
    HVNA3 = (double*****)malloc(sizeof(double****)*4);
    for (k=0; k<4; k++){
      HVNA3[k] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
      FNAN[0] = 0;
      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	if (Mc_AN==0)  Gc_AN = 0;
	else           Gc_AN = F_M2G[Mc_AN];

	HVNA3[k][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno0 = 1;
	  }
	  else{
            Gh_AN = natn[Gc_AN][h_AN];        
            Hwan = WhatSpecies[Gh_AN];
	    tno0 = Spe_Total_NO[Hwan];  
	  }    

	  HVNA3[k][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
	  for (i=0; i<tno0; i++){
	    HVNA3[k][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno0); 
	    size_HVNA3 += tno0;
	  }
	}
      }
    }

    /* CntHVNA2 */  

    if (Cnt_switch==1){

      size_CntHVNA2 = 0;
      CntHVNA2 = (double*****)malloc(sizeof(double****)*4);
      for (k=0; k<4; k++){
	CntHVNA2[k] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
	FNAN[0] = 0;
	for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

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

    /* CntHVNA3 */  

    if (Cnt_switch==1){

      size_CntHVNA3 = 0;
      CntHVNA3 = (double*****)malloc(sizeof(double****)*4);
      for (k=0; k<4; k++){
	CntHVNA3[k] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
	FNAN[0] = 0;
	for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	  if (Mc_AN==0) Gc_AN = 0;
	  else          Gc_AN = F_M2G[Mc_AN];

	  CntHVNA3[k][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	    if (Mc_AN==0){
	      tno0 = 1;
	    }
	    else{
	      Gh_AN = natn[Gc_AN][h_AN];        
	      Hwan = WhatSpecies[Gh_AN];
	      tno0 = Spe_Total_CNO[Hwan];  
	    }    

	    CntHVNA3[k][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
	    for (i=0; i<tno0; i++){
	      CntHVNA3[k][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno0); 
	      size_CntHVNA3 += tno0;
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

    Krylov_U = (double***)malloc(sizeof(double**)*(SpinP_switch+1));
    for (i=0; i<=SpinP_switch; i++){

      Krylov_U[i] = (double**)malloc(sizeof(double*)*(Matomnum+1));

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

	csize = Bnum + 2;
	
	Krylov_U[i][j] = (double*)malloc(sizeof(double)*rlmax_EC[j]*EKC_core_size[j]*csize);

        size_Krylov_U += rlmax_EC[j]*EKC_core_size[j]*csize;

      } /* j */

      if (i==0){
        MPI_Reduce(&My_KryDH, &KryDH, 1, MPI_DOUBLE, MPI_SUM, Host_ID, mpi_comm_level1);
        MPI_Reduce(&My_KryDS, &KryDS, 1, MPI_DOUBLE, MPI_SUM, Host_ID, mpi_comm_level1);

        if (myid==Host_ID && 2<=level_stdout){
	  printf("<Krylov parameters>  Av Krlov dimension H=%10.5f S=%10.5f\n",
                   KryDH/atomnum,KryDS/atomnum);
        }       
      }

    } /* i */

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
  PrintMemory("truncation: HVNA3",   sizeof(double)*size_HVNA3,       NULL);
  }
  if (Cnt_switch==1){
  PrintMemory("truncation: CntDS_VNA",sizeof(Type_DS_VNA)*size_CntDS_VNA,  NULL);
  PrintMemory("truncation: CntHVNA2", sizeof(double)*size_CntHVNA2,   NULL);
  PrintMemory("truncation: CntHVNA3", sizeof(double)*size_CntHVNA3,   NULL);
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
  ****************************************************/

  if (UCell_flag==1){

    N = My_NumGridC;

    if (SpinP_switch==3){ /* spin non-collinear */
      Density_Grid = (double**)malloc(sizeof(double*)*4); 
      for (k=0; k<=3; k++){
        Density_Grid[k] = (double*)malloc(sizeof(double)*N); 
        for (i=0; i<N; i++) Density_Grid[k][i] = 0.0;
      }
    }
    else{
      Density_Grid = (double**)malloc(sizeof(double*)*2); 
      for (k=0; k<=1; k++){
        Density_Grid[k] = (double*)malloc(sizeof(double)*N); 
        for (i=0; i<N; i++) Density_Grid[k][i] = 0.0;
      }
    }

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

    /* arrays for the partitions B and C */

    if (SpinP_switch==3){ /* spin non-collinear */
      Density_Grid_B = (double**)malloc(sizeof(double*)*4); 
      for (k=0; k<=3; k++){
        Density_Grid_B[k] = (double*)malloc(sizeof(double)*My_NumGridB_AB); 
        for (i=0; i<My_NumGridB_AB; i++) Density_Grid_B[k][i] = 0.0;
      }
    }
    else{
      Density_Grid_B = (double**)malloc(sizeof(double*)*2); 
      for (k=0; k<=1; k++){
        Density_Grid_B[k] = (double*)malloc(sizeof(double)*My_NumGridB_AB); 
        for (i=0; i<My_NumGridB_AB; i++) Density_Grid_B[k][i] = 0.0;
      }
    }

    ADensity_Grid_B = (double*)malloc(sizeof(double)*My_NumGridB_AB); 
    PCCDensity_Grid_B = (double*)malloc(sizeof(double)*My_NumGridB_AB); 
    for (i=0; i<My_NumGridB_AB; i++) PCCDensity_Grid_B[i] = 0.0;

    dVHart_Grid_B = (double*)malloc(sizeof(double)*My_Max_NumGridB); 
    for (i=0; i<My_Max_NumGridB; i++) dVHart_Grid_B[i] = 0.0;

    RefVxc_Grid_B = (double*)malloc(sizeof(double)*My_NumGridB_AB); 

    if (SpinP_switch==3){ /* spin non-collinear */
      Vxc_Grid_B = (double**)malloc(sizeof(double*)*4); 
      for (k=0; k<=3; k++){
        Vxc_Grid_B[k] = (double*)malloc(sizeof(double)*My_NumGridB_AB); 
      }
    }
    else{
      Vxc_Grid_B = (double**)malloc(sizeof(double*)*2); 
      for (k=0; k<=1; k++){
        Vxc_Grid_B[k] = (double*)malloc(sizeof(double)*My_NumGridB_AB);
      }
    }

    if (SpinP_switch==3){ /* spin non-collinear */
      Vpot_Grid_B = (double**)malloc(sizeof(double*)*4); 
      for (k=0; k<=3; k++){
        Vpot_Grid_B[k] = (double*)malloc(sizeof(double)*My_NumGridB_AB); 
      }
    }
    else{
      Vpot_Grid_B = (double**)malloc(sizeof(double*)*2); 
      for (k=0; k<=1; k++){
        Vpot_Grid_B[k] = (double*)malloc(sizeof(double)*My_NumGridB_AB);
      }
    }

    /* if (ProExpn_VNA==off) */
    if (ProExpn_VNA==0){
      VNA_Grid = (double*)malloc(sizeof(double)*N); 
      VNA_Grid_B = (double*)malloc(sizeof(double)*My_NumGridB_AB); 
    }

    /* electric energy by electric field */
    if (E_Field_switch==1){
      VEF_Grid = (double*)malloc(sizeof(double)*N); 
      VEF_Grid_B = (double*)malloc(sizeof(double)*My_NumGridB_AB); 
    }

    /* arrays for the partitions D */

    PCCDensity_Grid_D = (double*)malloc(sizeof(double)*My_NumGridD); 

    if (SpinP_switch==3){ /* spin non-collinear */
      Density_Grid_D = (double**)malloc(sizeof(double*)*4); 
      for (k=0; k<=3; k++){
        Density_Grid_D[k] = (double*)malloc(sizeof(double)*My_NumGridD); 
        for (i=0; i<My_NumGridD; i++) Density_Grid_D[k][i] = 0.0;
      }
    }
    else{
      Density_Grid_D = (double**)malloc(sizeof(double*)*2); 
      for (k=0; k<=1; k++){
        Density_Grid_D[k] = (double*)malloc(sizeof(double)*My_NumGridD); 
        for (i=0; i<My_NumGridD; i++) Density_Grid_D[k][i] = 0.0;
      }
    }

    if (SpinP_switch==3){ /* spin non-collinear */
      Vxc_Grid_D = (double**)malloc(sizeof(double*)*4); 
      for (k=0; k<=3; k++){
        Vxc_Grid_D[k] = (double*)malloc(sizeof(double)*My_NumGridD); 
      }
    }
    else{
      Vxc_Grid_D = (double**)malloc(sizeof(double*)*2); 
      for (k=0; k<=1; k++){
        Vxc_Grid_D[k] = (double*)malloc(sizeof(double)*My_NumGridD);
      }
    }

    /* Orbs_Grid */
    size_Orbs_Grid = 0;
    Orbs_Grid = (Type_Orbs_Grid***)malloc(sizeof(Type_Orbs_Grid**)*(Matomnum+1)); 
    Orbs_Grid[0] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*1); 
    Orbs_Grid[0][0] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*1); 
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = F_M2G[Mc_AN];
      Cwan = WhatSpecies[Gc_AN];
      /* AITUNE */
      Orbs_Grid[Mc_AN] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*GridN_Atom[Gc_AN]); 
      int Nc;
      for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
        Orbs_Grid[Mc_AN][Nc] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*Spe_Total_NO[Cwan]); 
	size_Orbs_Grid += Spe_Total_NO[Cwan];
      }
      /* AITUNE */
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

    /* Orbs_Grid_FNAN */
    size_Orbs_Grid_FNAN = 0;
    Orbs_Grid_FNAN = (Type_Orbs_Grid****)malloc(sizeof(Type_Orbs_Grid***)*(Matomnum+1)); 
    Orbs_Grid_FNAN[0] = (Type_Orbs_Grid***)malloc(sizeof(Type_Orbs_Grid**)*1); 
    Orbs_Grid_FNAN[0][0] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*1); 
    Orbs_Grid_FNAN[0][0][0] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*1); 

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

      Gc_AN = M2G[Mc_AN];    
      Orbs_Grid_FNAN[Mc_AN] = (Type_Orbs_Grid***)malloc(sizeof(Type_Orbs_Grid**)*(FNAN[Gc_AN]+1)); 

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

        Gh_AN = natn[Gc_AN][h_AN];

        if (G2ID[Gh_AN]!=myid){

	  Mh_AN = F_G2M[Gh_AN];
	  Rnh = ncn[Gc_AN][h_AN];
	  Hwan = WhatSpecies[Gh_AN];
          NO1 = Spe_Total_NO[Hwan];
          /* AITUNE */
          if (0<NumOLG[Mc_AN][h_AN]){
  	    Orbs_Grid_FNAN[Mc_AN][h_AN] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*NumOLG[Mc_AN][h_AN]); 
            int Nc;
	    for (Nc=0; Nc<NumOLG[Mc_AN][h_AN]; Nc++){
              Orbs_Grid_FNAN[Mc_AN][h_AN][Nc] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*NO1); 
              size_Orbs_Grid_FNAN += NO1;
	    }
	  }
          /* AITUNE */
	}
        else {
          Orbs_Grid_FNAN[Mc_AN][h_AN] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*1); 
          Orbs_Grid_FNAN[Mc_AN][h_AN][0] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*1); 
          size_Orbs_Grid_FNAN += 1;
        }
      }
    }

    alloc_first[3] = 0;

    /* PrintMemory */

    if (firsttime){

      int mul;
      if (SpinP_switch==3) mul = 4;
      else                 mul = 2;

      PrintMemory("truncation: Density_Grid",     sizeof(double)*N*mul,              NULL);
      PrintMemory("truncation: Vxc_Grid",         sizeof(double)*N*mul,              NULL);
      PrintMemory("truncation: RefVxc_Grid",      sizeof(double)*N,                  NULL);
      PrintMemory("truncation: Vpot_Grid",        sizeof(double)*N*mul,              NULL);
      PrintMemory("truncation: dVHart_Grid",      sizeof(double)*N,                  NULL);
      PrintMemory("truncation: Vpot_Grid_B",      sizeof(double)*My_NumGridB_AB*mul, NULL);
      PrintMemory("truncation: Density_Grid_B",   sizeof(double)*My_NumGridB_AB*mul, NULL);
      PrintMemory("truncation: ADensity_Grid_B",  sizeof(double)*My_NumGridB_AB,     NULL);
      PrintMemory("truncation: PCCDensity_Grid_B",sizeof(double)*My_NumGridB_AB,     NULL);
      PrintMemory("truncation: dVHart_Grid_B",    sizeof(double)*My_Max_NumGridB,    NULL);
      PrintMemory("truncation: RefVxc_Grid_B",    sizeof(double)*My_NumGridB_AB,     NULL);
      PrintMemory("truncation: Vxc_Grid_B",       sizeof(double)*My_NumGridB_AB*mul, NULL);
      PrintMemory("truncation: Density_Grid_D",   sizeof(double)*My_NumGridD*mul,    NULL);
      PrintMemory("truncation: Vxc_Grid_D",       sizeof(double)*My_NumGridD*mul,    NULL);
      PrintMemory("truncation: PCCDensity_Grid_D",sizeof(double)*My_NumGridD,        NULL);
      PrintMemory("truncation: Orbs_Grid",        sizeof(Type_Orbs_Grid)*size_Orbs_Grid,   NULL);
      PrintMemory("truncation: COrbs_Grid",       sizeof(Type_Orbs_Grid)*size_COrbs_Grid,  NULL);
      PrintMemory("truncation: Orbs_Grid_FNAN",   sizeof(Type_Orbs_Grid)*size_Orbs_Grid_FNAN,  NULL);

      if (ProExpn_VNA==0){
        PrintMemory("truncation: VNA_Grid",         sizeof(double)*N,                  NULL);
	PrintMemory("truncation: VNA_Grid_B",       sizeof(double)*My_NumGridB_AB,   NULL);
      }
      if (E_Field_switch==1){
	PrintMemory("truncation: VEF_Grid",         sizeof(double)*N,                NULL);
	PrintMemory("truncation: VEF_Grid_B",       sizeof(double)*My_NumGridB_AB,   NULL);
      }
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

  if (measure_time){
    dtime(&etime); 
    time16 += etime - stime;
  }

  if (measure_time){
    printf("myid=%5d time0 =%6.3f time1 =%6.3f time2 =%6.3f time3 =%6.3f time4 =%6.3f time5 =%6.3f\n",
            myid,time0,time1,time2,time3,time4,time5);
    printf("myid=%5d time6 =%6.3f time7 =%6.3f time8 =%6.3f time9 =%6.3f time10=%6.3f time11=%6.3f\n",
            myid,time6,time7,time8,time9,time10,time11);
    printf("myid=%5d time12=%6.3f time13=%6.3f time14=%6.3f time15=%6.3f time16=%6.3f\n",
            myid,time12,time13,time14,time15,time16);
  }

  /* for PrintMemory */
  firsttime = 0;
 
  /* for time */
  dtime(&TEtime);
  return (TEtime-TStime);
}







void Trn_System(int MD_iter, int CpyCell, int TCpyCell)
{
  int i,j,k,l,m,Rn,fan,san,tan,wanA,wanB,po0;
  int ct_AN,h_AN,m2,m3,i0,size_RMI1,size_array;
  int My_TFNAN,My_TSNAN,Gh_AN,LT_switch,Nloop;
  double r,rcutA,rcutB,rcut,dx,dy,dz,rcut_max;
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

#pragma omp parallel shared(myid,ScaleSize,Max_FSNAN,level_stdout,BCR,atv,Gxyz,Dis,ncn,natn,TCpyCell,CpyCell,atomnum,FNAN,SNAN,Spe_Atom_Cut1,WhatSpecies,M2G,Matomnum) private(OMPID,Nthrds,Nprocs,i,ct_AN,wanA,rcutA,j,wanB,rcutB,rcut,Rn,dx,dy,dz,r,l,k,fDis,fncn,fnan2,sDis,sncn,snan2,size_array,rcut_max)
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
	rcut = rcutA + rcutB;

        if (rcut<BCR) rcut_max = BCR;
        else          rcut_max = rcut; 

	for (Rn=0; Rn<=TCpyCell; Rn++){

	  if ((ct_AN==j) && Rn==0){
	    natn[ct_AN][0] = ct_AN;
	    ncn[ct_AN][0]  = 0;
	    Dis[ct_AN][0]  = 0.0;
	  }
            
	  else{

	    dx = fabs(Gxyz[ct_AN][1] - Gxyz[j][1] - atv[Rn][1]);
	    dy = fabs(Gxyz[ct_AN][2] - Gxyz[j][2] - atv[Rn][2]);
	    dz = fabs(Gxyz[ct_AN][3] - Gxyz[j][3] - atv[Rn][3]);

            if (dx<=rcut_max && dy<=rcut_max && dz<=rcut_max){

	      r = sqrt(dx*dx + dy*dy + dz*dz);

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
	printf("<truncation> myid=%4d CpyCell=%2d ct_AN=%4d FNAN SNAN %3d %3d\n",
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
                     Fixed_FNAN_SNAN
  ****************************************************/

  if (orderN_FNAN_SNAN_flag==1)  Fixed_FNAN_SNAN();

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

  My_TFNAN = 0;
  My_TSNAN = 0;
  for (i=1; i<=Matomnum; i++){
    ct_AN = M2G[i];
    My_TFNAN = My_TFNAN + FNAN[ct_AN];
    My_TSNAN = My_TSNAN + SNAN[ct_AN];
  }

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

    MPI_Bcast(&natn[ct_AN][0], FNAN[ct_AN]+SNAN[ct_AN]+1,
              MPI_INT, ID, mpi_comm_level1);
    MPI_Bcast(&ncn[ct_AN][0], FNAN[ct_AN]+SNAN[ct_AN]+1,
              MPI_INT, ID, mpi_comm_level1);
    MPI_Bcast(&Dis[ct_AN][0], FNAN[ct_AN]+SNAN[ct_AN]+1,
              MPI_DOUBLE, ID, mpi_comm_level1);
  }  

  if (myid==Host_ID && 0<level_stdout){
    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      if (ct_AN<=20 && level_stdout<=1){
        printf("<truncation> CpyCell=%2d ct_AN=%4d FNAN SNAN %3d %3d\n",
                 CpyCell,ct_AN,FNAN[ct_AN],SNAN[ct_AN]);fflush(stdout);
      }
    }

    if (20<atomnum && level_stdout<=1){
      printf("     ..........\n");
      printf("     ......\n\n");
    }
  }
}



void Estimate_Trn_System(int CpyCell, int TCpyCell)
{
  /****************************************************
     FNAN, SNAN, Max_FNAN, Max_FSNAN are determined
               by the physical truncation
  ****************************************************/

  int i,j,ct_AN,Rn,wanA,wanB;
  double r,rcutA,rcutB,rcut;
  double dx,dy,dz,rcut_max;
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

#pragma omp parallel shared(CpyCell,level_stdout,BCR,Spe_WhatAtom,atv,Gxyz,TCpyCell,SNAN,FNAN,Spe_Atom_Cut1,WhatSpecies,M2G,Matomnum,atomnum) private(i,ct_AN,wanA,rcutA,j,wanB,rcutB,rcut,Rn,dx,dy,dz,r,spe1,spe2,OMPID,Nthrds,Nprocs,rcut_max)

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
	rcut = rcutA + rcutB;

        if (rcut<BCR) rcut_max = BCR;
        else          rcut_max = rcut; 

	for (Rn=0; Rn<=TCpyCell; Rn++){

	  if ((ct_AN==j) && Rn==0){
	    /* Nothing to be done */
	  }
 
	  else{

	    dx = fabs(Gxyz[ct_AN][1] - Gxyz[j][1] - atv[Rn][1]);
	    dy = fabs(Gxyz[ct_AN][2] - Gxyz[j][2] - atv[Rn][2]);
	    dz = fabs(Gxyz[ct_AN][3] - Gxyz[j][3] - atv[Rn][3]);

            if (dx<=rcut_max && dy<=rcut_max && dz<=rcut_max){

	      r = sqrt(dx*dx + dy*dy + dz*dz);

	      if (r<=rcut)     FNAN[ct_AN]++;
	      else if (r<=BCR) SNAN[ct_AN]++;
	    }

	  }
	}
      }

      if (2<=level_stdout){
	printf("<truncation> CpyCell=%2d ct_AN=%2d FNAN SNAN %2d %2d\n",
	       CpyCell,i,FNAN[ct_AN],SNAN[ct_AN]);
      }
  
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

	    k = 0; po = 0;
            RMI1[Mc_AN][i][j] = -1;
	    do {
	      if (natn[ig][k]==jg && ncn[ig][k]==Rn){
		RMI1[Mc_AN][i][j] = k;
		po = 1;
	      }
	      k++;
	    } while (po==0 && k<=FNAN[ig]);

	    /* FNAN + SNAN */

	    k = 0; po = 0;
            RMI2[Mc_AN][i][j] = -1;
	    do {
	      if (natn[ig][k]==jg && ncn[ig][k]==Rn){
		RMI2[Mc_AN][i][j] = k;
		po = 1;
	      }
	      k++;
	    } while(po==0 && k<=(FNAN[ig]+SNAN[ig]));

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
  int size_GListTAtoms1;
  int po,N3[4],i,j,k;
  int size_GridListAtom,size_MGridListAtom;
  int NOC[4],nn1,nn2,nn3,l1,l2,l3,lmax;
  int ct_AN,N,n1,n2,n3,Cwan,Rn,Rn1,Ng1,Ng2,Ng3;
  int p,q,r,s,pmax,qmax,rmax,smax,popt,qopt,ropt,sopt;
  int Nct,MinN,Scale,ScaleA,ScaleB,ScaleC,Nm[4];
  int Mc_AN,Mc_AN0,Gc_AN,h_AN,h_AN0,Mh_AN,Gh_AN,Nh,Rh,Nog,GNh,GRh;
  int ll1,ll2,ll3,Nnb,My_Max_NumOLG;
  int lll1,lll2,lll3,GRh1,Nc,GNc,GRc,size_array;
  int *TAtoms0,*TCells0,*TAtoms1,*TAtoms2;
  int **Tmp_GridListAtom,**Tmp_CellListAtom;
  int nmin[4],nmax[4],Np;
  double LgTN,LgN,Lg2,Lg3,Lg5,Lg7,DouN[4],A1,A2,A3;
  double MinD,MinR,CutR2,r2,coef;
  double sa,sa_cri,tmp[4],Cxyz[4];
  double b[4],c[4],v[4],rcut;
  double xc,yc,zc,xm,ym,zm;
  double dx,dy,dz,sn1,sn2,sn3;
  double B2,C2,CellV;
  double S_Lng,L_Lng,LngA,LngB,LngC,x,y,z;
  double GVolume,buffer_scale,GridV;
  double stime,etime;
  double time0,time1,time2,time3,time4,time5,time6,time7;
  double time8,time9,time10,time11,time12,time13,time14;

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

  time0 = 0.0; time1 = 0.0; time2 = 0.0; time3 = 0.0; time4 = 0.0;
  time5 = 0.0; time6 = 0.0; time7 = 0.0; time8 = 0.0; time9 = 0.0;

  /****************************************************
                Reciprocal lattice vectors
  ****************************************************/

  if (measure_time) dtime(&stime); 

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

      /* find proper N1, N2, and N3 */ 

      Cross_Product(tv[2],tv[3],tmp);
      A1 = PI*PI*Dot_Product(tmp,tmp)/Cell_Volume/Cell_Volume; 
      
      Cross_Product(tv[3],tv[1],tmp);
      A2 = PI*PI*Dot_Product(tmp,tmp)/Cell_Volume/Cell_Volume; 

      Cross_Product(tv[1],tv[2],tmp);
      A3 = PI*PI*Dot_Product(tmp,tmp)/Cell_Volume/Cell_Volume; 

      DouN[1] = sqrt(Grid_Ecut/A1);
      DouN[2] = sqrt(Grid_Ecut/A2);
      DouN[3] = sqrt(Grid_Ecut/A3);

      Lg2 = log(2);
      Lg3 = log(3);
      Lg5 = log(5);
      Lg7 = log(7);

      for (i=1; i<=3; i++){

        LgN = log(DouN[i]);
        MinD = 10e+10;
        
        pmax = ceil(LgN/Lg2); 
        for (p=0; p<=pmax; p++){
          qmax = ceil((LgN-p*Lg2)/Lg3); 
          for (q=0; q<=qmax; q++){
            rmax = ceil((LgN-p*Lg2-q*Lg3)/Lg5); 
            for (r=0; r<=rmax; r++){
              smax = ceil((LgN-p*Lg2-q*Lg3-r*Lg5)/Lg7); 
              for (s=0; s<=smax; s++){
                
                LgTN = p*Lg2+q*Lg3+r*Lg5+s*Lg7;
                
                if (fabs(LgTN-LgN)<MinD){
                  MinD = fabs(LgTN-LgN);
                  popt = p;
                  qopt = q;
                  ropt = r;
                  sopt = s;
                }
	      }
	    }
	  }
        }
        
        k = 1;
        for (p=0; p<popt; p++) k *= 2;
        for (q=0; q<qopt; q++) k *= 3;
        for (r=0; r<ropt; r++) k *= 5;
        for (s=0; s<sopt; s++) k *= 7;

        if      (i==1) Ngrid1 = k;
        else if (i==2) Ngrid2 = k;
        else if (i==3) Ngrid3 = k;
      }

      /* adjust Ngrid for NEGF  */

      if (Solver==4) {
        TRAN_adjust_Ngrid(mpi_comm_level1, &Ngrid1, &Ngrid2, &Ngrid3);
      }

    } /* if (Ngrid_fixed_flag==0) */

    /* calculate gtv, rgtv, A2, B2, C2, and used cutoff energies */

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

    A2 = A2/4.0;  /* note: change the unit from Hatree to Rydberg by multiplying 1/2 */
    B2 = B2/4.0;
    C2 = C2/4.0;

    if (Ngrid_fixed_flag==1)  Grid_Ecut = (A2 + B2 + C2)/3.0;

    /* for calculation of the delta-factor */

    if (MD_switch==16) Ngrid_fixed_flag = 1;

    /* print information to std output */

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

  if (measure_time){
    dtime(&etime); 
    time0 += etime - stime;
  }

  /****************************************************
       Setting the center of unit cell and grids
  ****************************************************/

  if (measure_time) dtime(&stime); 

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

  if (measure_time){
    dtime(&etime); 
    time1 += etime - stime;
  }

  /**********************************
    allocation of arrays: 

    Tmp_GridListAtom
    Tmp_CellListAtom
    MGridListAtom
  **********************************/
 
  Tmp_GridListAtom = (int**)malloc(sizeof(int*)*(Matomnum+MatomnumF+1));
  Tmp_CellListAtom = (int**)malloc(sizeof(int*)*(Matomnum+MatomnumF+1));
  MGridListAtom = (int**)malloc(sizeof(int*)*(Matomnum+1));
  Tmp_GridListAtom[0] = (int*)malloc(sizeof(int)*1);
  Tmp_CellListAtom[0] = (int*)malloc(sizeof(int)*1);
  MGridListAtom[0] = (int*)malloc(sizeof(int)*1);
  alloc_first[2] = 0;

  /****************************************************
   1) find a neighbouring point of the atom Mc_AN
   2) the ranges which deterinie a box on the atom Mc_AN 
   3) determine whether overlap exists or not
  ****************************************************/

  if (measure_time) dtime(&stime); 

  /* for allocation of arrays */

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];
    Cwan = WhatSpecies[Gc_AN];
    rcut = Spe_Atom_Cut1[Cwan] + 0.5;

    for (k=1; k<=3; k++){

      if      (k==1){ i = 2; j = 3; }
      else if (k==2){ i = 3; j = 1; }
      else if (k==3){ i = 1; j = 2; }

      b[1] = tv[i][1];
      b[2] = tv[i][2];
      b[3] = tv[i][3];

      c[1] = tv[j][1];
      c[2] = tv[j][2];
      c[3] = tv[j][3];

      Cross_Product(b,c,v);
      coef = 1.0/sqrt(fabs( Dot_Product(v,v) ));

      v[1] = coef*v[1];
      v[2] = coef*v[2];
      v[3] = coef*v[3];

      Cxyz[1] = Gxyz[Gc_AN][1] + rcut*v[1] - Grid_Origin[1];
      Cxyz[2] = Gxyz[Gc_AN][2] + rcut*v[2] - Grid_Origin[2];
      Cxyz[3] = Gxyz[Gc_AN][3] + rcut*v[3] - Grid_Origin[3];

      /* find the maximum range of grids */
      nmax[k] = Dot_Product(Cxyz,rgtv[k])*0.5/PI;

      Cxyz[1] = Gxyz[Gc_AN][1] - rcut*v[1] - Grid_Origin[1];
      Cxyz[2] = Gxyz[Gc_AN][2] - rcut*v[2] - Grid_Origin[2];
      Cxyz[3] = Gxyz[Gc_AN][3] - rcut*v[3] - Grid_Origin[3];

      /* find the mimum range of grids */
      nmin[k] = Dot_Product(Cxyz,rgtv[k])*0.5/PI;

      if (nmax[k]<nmin[k]){
        i = nmin[k];
        j = nmax[k];
        nmin[k] = j;
        nmax[k] = i;
      } 
  
    } /* k */  

    /* allocation of arrays */ 

    Np = (nmax[1]-nmin[1]+1)*(nmax[2]-nmin[2]+1)*(nmax[3]-nmin[3]+1)*3/2;
    
    Tmp_GridListAtom[Mc_AN] = (int*)malloc(sizeof(int)*Np);
    Tmp_CellListAtom[Mc_AN] = (int*)malloc(sizeof(int)*Np);
    MGridListAtom[Mc_AN] = (int*)malloc(sizeof(int)*Np);

  } /* Mc_AN */

  /* store Tmp_GridListAtom and Tmp_CellListAtom */

  size_array = (int)(Max_GridN_Atom*ScaleSize);
  TempGrid = (int*)malloc(sizeof(int)*size_array); 
  TempCell = (int*)malloc(sizeof(int)*size_array); 

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];
    Cwan = WhatSpecies[Gc_AN];
    rcut = Spe_Atom_Cut1[Cwan] + 0.5;

    for (k=1; k<=3; k++){

      if      (k==1){ i = 2; j = 3; }
      else if (k==2){ i = 3; j = 1; }
      else if (k==3){ i = 1; j = 2; }

      b[1] = tv[i][1];
      b[2] = tv[i][2];
      b[3] = tv[i][3];

      c[1] = tv[j][1];
      c[2] = tv[j][2];
      c[3] = tv[j][3];

      Cross_Product(b,c,v);
      coef = 1.0/sqrt(fabs( Dot_Product(v,v) ));

      v[1] = coef*v[1];
      v[2] = coef*v[2];
      v[3] = coef*v[3];

      Cxyz[1] = Gxyz[Gc_AN][1] + rcut*v[1] - Grid_Origin[1];
      Cxyz[2] = Gxyz[Gc_AN][2] + rcut*v[2] - Grid_Origin[2];
      Cxyz[3] = Gxyz[Gc_AN][3] + rcut*v[3] - Grid_Origin[3];

      /* find the maximum range of grids */
      nmax[k] = Dot_Product(Cxyz,rgtv[k])*0.5/PI;

      Cxyz[1] = Gxyz[Gc_AN][1] - rcut*v[1] - Grid_Origin[1];
      Cxyz[2] = Gxyz[Gc_AN][2] - rcut*v[2] - Grid_Origin[2];
      Cxyz[3] = Gxyz[Gc_AN][3] - rcut*v[3] - Grid_Origin[3];

      /* find the mimum range of grids */
      nmin[k] = Dot_Product(Cxyz,rgtv[k])*0.5/PI;
  
      if (nmax[k]<nmin[k]){
        i = nmin[k];
        j = nmax[k];
        nmin[k] = j;
        nmax[k] = i;
      } 

    } /* k */  

    CutR2 = Spe_Atom_Cut1[Cwan]*Spe_Atom_Cut1[Cwan];

    Nct = 0;
    for (n1=nmin[1]; n1<=nmax[1]; n1++){
      for (n2=nmin[2]; n2<=nmax[2]; n2++){
	for (n3=nmin[3]; n3<=nmax[3]; n3++){

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
	    if (estimate_switch!=1){
	      TempGrid[Nct+1] = N;
	      TempCell[Nct+1] = Rn;
	    }             
	    Nct++;
	  }
	}
      }
    }

    Np = (nmax[1]-nmin[1]+1)*(nmax[2]-nmin[2]+1)*(nmax[3]-nmin[3]+1)*3/2;
    if (Np<Nct){
      printf("Invalid access in truncation.c\n"); 
      MPI_Finalize();
      exit(0); 
    }

    GridN_Atom[Gc_AN] = Nct;

    if (estimate_switch!=1){

      /* sorting */
      qsort_int((long)Nct,TempGrid,TempCell);

      for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
	Tmp_GridListAtom[Mc_AN][Nc] = TempGrid[Nc+1];
	Tmp_CellListAtom[Mc_AN][Nc] = TempCell[Nc+1];
      }
    }
  }

  free(TempCell);
  free(TempGrid);

  /* calculate size_GridListAtom */

  size_GridListAtom = 0;
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    size_GridListAtom += GridN_Atom[Gc_AN];
    if (Max_GridN_Atom<GridN_Atom[Gc_AN]) Max_GridN_Atom = GridN_Atom[Gc_AN];
  }

  if (measure_time){
    dtime(&etime); 
    time2 += etime - stime;
  }

  /****************************************************
   MPI: 

       GridN_Atom
  ****************************************************/

  if (measure_time) dtime(&stime); 

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    ID = G2ID[ct_AN];
    MPI_Bcast(&GridN_Atom[ct_AN], 1, MPI_INT, ID, mpi_comm_level1);
  }

  if (myid==Host_ID && estimate_switch==0 && 0<level_stdout){
    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      if (ct_AN<=20 && level_stdout<=1){
         printf("Num. of grids overlapping with atom %4d = %4d\n",
                 ct_AN, GridN_Atom[ct_AN]);
      }
    }

    if (20<atomnum && level_stdout<=1){
      printf("     ..........\n");
      printf("     ......\n\n");
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
    size_MGridListAtom += GridN_Atom[Gc_AN];
  }

  /* PrintMemory */
  if (firsttime){
  PrintMemory("truncation: Tmp_GridListAtom", sizeof(int)*size_MGridListAtom, NULL);
  PrintMemory("truncation: Tmp_CellListAtom", sizeof(int)*size_MGridListAtom, NULL);
  }

  if (measure_time){
    dtime(&etime); 
    time3 += etime - stime;
  }

  /****************************************************
   MPI: 

       Tmp_GridListAtom
       Tmp_CellListAtom
  ****************************************************/

  if (measure_time) dtime(&stime); 

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

  /* Tmp_CellListAtom */

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

  if (measure_time){
    dtime(&etime); 
    time4 += etime - stime;
  }

  if (measure_time) dtime(&stime); 

  if (estimate_switch!=1){

    /****************************************************
            Find overlap grids between two orbitals
    ****************************************************/
    
    if (estimate_switch==0){
      size_GListTAtoms1 = 0;

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

          GListTAtoms1[0] = (int**)malloc(sizeof(int*)*1);
          GListTAtoms2[0] = (int**)malloc(sizeof(int*)*1);

          GListTAtoms1[0][0] = (int*)malloc(sizeof(int)*1);
          GListTAtoms2[0][0] = (int*)malloc(sizeof(int)*1);
        }
      }

      else{

        if (estimate_switch==0){

          GListTAtoms1[Mc_AN] = (int**)malloc(sizeof(int*)*(FNAN[Gc_AN]+1));
          GListTAtoms2[Mc_AN] = (int**)malloc(sizeof(int*)*(FNAN[Gc_AN]+1));
        }

        h_AN0 = 0;

#pragma omp parallel shared(List_YOUSO,GListTAtoms1,GListTAtoms2,ScaleSize,Max_NumOLG,size_GListTAtoms1,level_stdout,NumOLG,estimate_switch,CpyCell,Mc_AN,Tmp_CellListAtom,Tmp_GridListAtom,GridN_Atom,atv_ijk,ncn,F_G2M,natn,h_AN0,FNAN,Gc_AN) private(OMPID,Nthrds,Nprocs,h_AN,Gh_AN,Mh_AN,Rh,l1,l2,l3,Nog,Nc,Nh,GNh,GRh,ll1,ll2,ll3,lll1,lll2,lll3,GRh1,po,GNc,GRc,TAtoms0,TCells0,TAtoms1,TAtoms2,size_array,i)
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

		    } /* while (...) */

		    /* for Nc==GridN_Atom[Gc_AN] */

		    Nc--;
		    if (Nc<0) Nc = 0;

		  } /* if (abs.... ) */
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

		    size_GListTAtoms1 += NumOLG[Mc_AN][i];

		    GListTAtoms1[Mc_AN][i] = (int*)malloc(sizeof(int)*NumOLG[Mc_AN][i]);
		    GListTAtoms2[Mc_AN][i] = (int*)malloc(sizeof(int)*NumOLG[Mc_AN][i]);

		  }
		}
	      }

#pragma omp barrier
#pragma omp flush(GListTAtoms1,GListTAtoms2)

              if (h_AN<=FNAN[Gc_AN]){

		for (Nog=0; Nog<NumOLG[Mc_AN][h_AN]; Nog++){

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

      if (firsttime){
      PrintMemory("truncation: GListTAtoms1", sizeof(int)*size_GListTAtoms1, NULL);
      PrintMemory("truncation: GListTAtoms2", sizeof(int)*size_GListTAtoms1, NULL);
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

  if (measure_time){
    dtime(&etime); 
    time5 += etime - stime;
  }

  /****************************************************
       Tmp_GridListAtom -> GridListAtom
       Tmp_CellListAtom -> CellListAtom                     
  ****************************************************/

  if (measure_time) dtime(&stime); 

  size_GridListAtom = 0;

  GridListAtom = (int**)malloc(sizeof(int*)*(Matomnum+1));
  GridListAtom[0] = (int*)malloc(sizeof(int)*1);
  size_GridListAtom++;
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    GridListAtom[Mc_AN] = (int*)malloc(sizeof(int)*GridN_Atom[Gc_AN]);
    size_GridListAtom += GridN_Atom[Gc_AN];
    for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
      GridListAtom[Mc_AN][Nc] = Tmp_GridListAtom[Mc_AN][Nc];
    }
  }

  CellListAtom = (int**)malloc(sizeof(int*)*(Matomnum+1));
  CellListAtom[0] = (int*)malloc(sizeof(int)*1);
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = F_M2G[Mc_AN];
    CellListAtom[Mc_AN] = (int*)malloc(sizeof(int)*GridN_Atom[Gc_AN]);
    for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
      CellListAtom[Mc_AN][Nc] = Tmp_CellListAtom[Mc_AN][Nc];
    }
  }

  /* PrintMemory */
  if (firsttime){
  PrintMemory("truncation: GridListAtom", sizeof(int)*size_GridListAtom, NULL);
  PrintMemory("truncation: CellListAtom", sizeof(int)*size_GridListAtom, NULL);
  PrintMemory("truncation: MGridListAtom", sizeof(int)*size_GridListAtom, NULL);
  }

  /****************************************************
   construct the data structure for MPI communications
   for grid data
  ****************************************************/

  if (estimate_switch==0){

    Construct_MPI_Data_Structure_Grid();

    Ng1 = Max_Grid_Index[1] - Min_Grid_Index[1] + 1;
    Ng2 = Max_Grid_Index[2] - Min_Grid_Index[2] + 1;
    Ng3 = Max_Grid_Index[3] - Min_Grid_Index[3] + 1;

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

      dtime(&Stime_atom);

      Gc_AN = F_M2G[Mc_AN];

#pragma omp parallel shared(Min_Grid_Index,atv_ijk,Ng1,Ng2,Ng3,Ngrid1,Ngrid2,Ngrid3,MGridListAtom,Mc_AN,Tmp_GridListAtom,Tmp_CellListAtom,GridN_Atom,Gc_AN) private(Nc,GNc,GRc,N3,n1,n2,n3,N,OMPID,Nthrds,Nprocs)
      {

	/* get info. on OpenMP */ 
  
	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	for (Nc=OMPID; Nc<GridN_Atom[Gc_AN]; Nc+=Nthrds){

	  GNc = Tmp_GridListAtom[Mc_AN][Nc];
          GRc = Tmp_CellListAtom[Mc_AN][Nc];

	  GN2N(GNc,N3);

	  n1 = N3[1] + Ngrid1*atv_ijk[GRc][1] - Min_Grid_Index[1];  
	  n2 = N3[2] + Ngrid2*atv_ijk[GRc][2] - Min_Grid_Index[2];  
	  n3 = N3[3] + Ngrid3*atv_ijk[GRc][3] - Min_Grid_Index[3];  

	  MGridListAtom[Mc_AN][Nc] = n1*Ng2*Ng3 + n2*Ng3 + n3;
	}

      } /* #pragma omp parallel */

      dtime(&Etime_atom); 
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
    }

  }

  if (measure_time){
    dtime(&etime); 
    time6 += etime - stime;
  }

  /* for PrintMemory */
  firsttime = 0;

  /****************************************************
                          Free
  ****************************************************/

  if (measure_time) dtime(&stime); 

  /* MPI_Barrier */
  MPI_Barrier(mpi_comm_level1);

  for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
    free(Tmp_CellListAtom[Mc_AN]);
    free(Tmp_GridListAtom[Mc_AN]);
  }
  free(Tmp_CellListAtom);
  free(Tmp_GridListAtom);

  if (alloc_first[2]==0 && estimate_switch!=0){

    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      free(MGridListAtom[Mc_AN]);
    }
    free(MGridListAtom);

    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      free(GridListAtom[Mc_AN]);
    }
    free(GridListAtom);

    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      free(CellListAtom[Mc_AN]);
    }
    free(CellListAtom);
  }

  if (measure_time){
    dtime(&etime); 
    time7 += etime - stime;
  }

  if (measure_time){
    printf("UCell_Box myid=%5d time0=%6.3f time1=%6.3f time2=%6.3f time3=%6.3f time4=%6.3f time5=%6.3f time6=%6.3f time7=%6.3f\n",
            myid,time0,time1,time2,time3,time4,time5,time6,time7);
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

    /* for PrintMemory */
    firsttime=0;

    /****************************************************
           setting of parameters of periodic cells
    ****************************************************/

    Generation_ATV(CpyN);
  }

  /* return */

  return TN;
}











void Free_truncation(int CpyN, int TN, int Free_switch)
{
  int i,j,n,Nc;

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
  }
} 








void free_arrays_truncation0()
{
  int i,j,k,m,ct_AN,h_AN,Gh_AN,Hwan;
  int tno0,tno1,tno,Cwan,so,s1,s2;
  int num,wan,n2,wanA,Gi;
  int NO1,Rnh,Mh_AN,Nc;
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
        else if ( (Hub_U_switch==0 || Hub_U_occupation!=1) && 0<k && Matomnum<Mc_AN){
          Gc_AN = S_M2G[Mc_AN];
          Cwan = WhatSpecies[Gc_AN];
          tno0 = 1;
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
	    else if ( (Matomnum+1)<=Mc_AN ){
	      tno1 = List_YOUSO[20];  
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
	      tno0 = Spe_Total_CNO[Cwan];  
	      fan = FNAN[Gc_AN];
	    }    

	    for (h_AN=0; h_AN<(fan+1); h_AN++){

	      if (Mc_AN==0){
		tno1 = 1;  
	      }
	      else if ( (Matomnum+1)<=Mc_AN ){
		tno1 = List_YOUSO[20];  
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

      /* HVNA3 */

      for (k=0; k<4; k++){
	FNAN[0] = 0;
	for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	  if (Mc_AN==0)  Gc_AN = 0;
	  else           Gc_AN = F_M2G[Mc_AN];

	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	    if (Mc_AN==0){
	      tno0 = 1;
	    }
	    else{
	      Gh_AN = natn[Gc_AN][h_AN];        
	      Hwan = WhatSpecies[Gh_AN];
	      tno0 = Spe_Total_NO[Hwan];  
	    }    

	    for (i=0; i<tno0; i++){
	      free(HVNA3[k][Mc_AN][h_AN][i]);
	    }
	    free(HVNA3[k][Mc_AN][h_AN]);
	  }
	  free(HVNA3[k][Mc_AN]);
	}
	free(HVNA3[k]);
      }
      free(HVNA3);

      /* CntHVNA2 */

      if (Cnt_switch==1){

	for (k=0; k<4; k++){
	  FNAN[0] = 0;
	  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

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

      /* CntHVNA3 */

      if (Cnt_switch==1){

	for (k=0; k<4; k++){
	  FNAN[0] = 0;
	  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	    if (Mc_AN==0) Gc_AN = 0;
	    else          Gc_AN = F_M2G[Mc_AN];

	    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	      if (Mc_AN==0){
		tno0 = 1;
	      }
	      else{
		Gh_AN = natn[Gc_AN][h_AN];        
		Hwan = WhatSpecies[Gh_AN];
		tno0 = Spe_Total_CNO[Hwan];  
	      }    

	      for (i=0; i<tno0; i++){
		free(CntHVNA3[k][Mc_AN][h_AN][i]);
	      }
	      free(CntHVNA3[k][Mc_AN][h_AN]);
	    }
	    free(CntHVNA3[k][Mc_AN]);
	  }
	  free(CntHVNA3[k]);
	}
	free(CntHVNA3);
      }
    }  

    if (Solver==8) { /* Krylov subspace method */

      for (i=0; i<=SpinP_switch; i++){
	for (j=0; j<=Matomnum; j++){
          free(Krylov_U[i][j]);
	}
        free(Krylov_U[i]);
      }
      free(Krylov_U);

      for (i=0; i<=SpinP_switch; i++){
	for (j=0; j<=Matomnum; j++){
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
  ****************************************************/

  if (alloc_first[0]==0){

    FNAN[0] = 0; 
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

      if (Mc_AN==0) Gc_AN = 0;
      else          Gc_AN = M2G[Mc_AN];

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
        free(GListTAtoms2[Mc_AN][h_AN]);
        free(GListTAtoms1[Mc_AN][h_AN]);
      }
      free(GListTAtoms2[Mc_AN]);
      free(GListTAtoms1[Mc_AN]);
    }
    free(GListTAtoms2);
    free(GListTAtoms1);
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
    free(dVHart_Grid);
    free(RefVxc_Grid_B);

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

    /* arrays for the partitions B and C */

    if (SpinP_switch==3){ /* spin non-collinear */
      for (k=0; k<=3; k++){
        free(Density_Grid_B[k]);
      }
      free(Density_Grid_B);
    }
    else{
      for (k=0; k<=1; k++){
        free(Density_Grid_B[k]);
      }
      free(Density_Grid_B);
    }

    free(ADensity_Grid_B);
    free(PCCDensity_Grid_B);
    free(dVHart_Grid_B);

    if (SpinP_switch==3){ /* spin non-collinear */
      for (k=0; k<=3; k++){
        free(Vxc_Grid_B[k]);
      }
      free(Vxc_Grid_B);
    }
    else{
      for (k=0; k<=1; k++){
        free(Vxc_Grid_B[k]);
      }
      free(Vxc_Grid_B);
    }

    if (SpinP_switch==3){ /* spin non-collinear */
      for (k=0; k<=3; k++){
        free(Vpot_Grid_B[k]);
      }
      free(Vpot_Grid_B);
    }
    else{
      for (k=0; k<=1; k++){
        free(Vpot_Grid_B[k]);
      }
      free(Vpot_Grid_B);
    }

    /* if (ProExpn_VNA==off) */
    if (ProExpn_VNA==0){
      free(VNA_Grid);
      free(VNA_Grid_B);
    }

    /* electric energy by electric field */
    if (E_Field_switch==1){
      free(VEF_Grid);
      free(VEF_Grid_B);
    }

    /* arrays for the partition D */

    free(PCCDensity_Grid_D);

    if (SpinP_switch==3){ /* spin non-collinear */
      for (k=0; k<=3; k++){
        free(Density_Grid_D[k]);
      }
      free(Density_Grid_D);
    }
    else{
      for (k=0; k<=1; k++){
        free(Density_Grid_D[k]);
      }
      free(Density_Grid_D);
    }

    if (SpinP_switch==3){ /* spin non-collinear */
      for (k=0; k<=3; k++){
        free(Vxc_Grid_D[k]);
      }
      free(Vxc_Grid_D);
    }
    else{
      for (k=0; k<=1; k++){
        free(Vxc_Grid_D[k]);
      }
      free(Vxc_Grid_D);
    }

    /* Orbs_Grid */

    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      if (Mc_AN==0){
        Gc_AN = 0;
        num = 1;
      }
      else{
        Gc_AN = F_M2G[Mc_AN];
        num = GridN_Atom[Gc_AN];
      }

      for (Nc=0; Nc<num; Nc++){
        free(Orbs_Grid[Mc_AN][Nc]);
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

    /* Orbs_Grid_FNAN */

    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

      if (Mc_AN==0){
        free(Orbs_Grid_FNAN[0][0][0]);
        free(Orbs_Grid_FNAN[0][0]);
      }
      else{

	Gc_AN = M2G[Mc_AN];    

	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  Gh_AN = natn[Gc_AN][h_AN];

	  if (G2ID[Gh_AN]!=myid){

	    Mh_AN = F_G2M[Gh_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    NO1 = Spe_Total_NO[Hwan];
            num = NumOLG[Mc_AN][h_AN];  
	  }
	  else {
            num = 1;
	  }

          if (0<NumOLG[Mc_AN][h_AN]){
            for (Nc=0; Nc<num; Nc++){
              free(Orbs_Grid_FNAN[Mc_AN][h_AN][Nc]);
            }
            free(Orbs_Grid_FNAN[Mc_AN][h_AN]);
	  }

	} /* h_AN */
      } /* else */

      free(Orbs_Grid_FNAN[Mc_AN]);
    }
    free(Orbs_Grid_FNAN);

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

    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      free(MGridListAtom[Mc_AN]);
    }
    free(MGridListAtom);

    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      free(GridListAtom[Mc_AN]);
    }
    free(GridListAtom);

    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
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
  int *flag_DoubleCounting;
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

   int flag_DoubleCounting[atomnum+1]
  *********************************/

  flag_DoubleCounting = (int*)malloc(sizeof(int)*(atomnum+1));

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

  /* initialize flag_DoubleCounting */

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    flag_DoubleCounting[Gc_AN] = 0;
  }

  /* find F_Rcv_Num */

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    for (Lh_AN=0; Lh_AN<=FNAN[Gc_AN]; Lh_AN++){
      Gh_AN = natn[Gc_AN][Lh_AN];
      ID1 = G2ID[Gh_AN];
      if (flag_DoubleCounting[Gh_AN]==0 && ID1!=myid){
        F_Rcv_Num[ID1]++;
        flag_DoubleCounting[Gh_AN] = 1;
      }
    }
  }

  /* find S_Rcv_Num */

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    for (Lh_AN=(FNAN[Gc_AN]+1); Lh_AN<=(FNAN[Gc_AN]+SNAN[Gc_AN]); Lh_AN++){
      Gh_AN = natn[Gc_AN][Lh_AN];
      ID1 = G2ID[Gh_AN];
      if (flag_DoubleCounting[Gh_AN]==0 && ID1!=myid){
        S_Rcv_Num[ID1]++;
        flag_DoubleCounting[Gh_AN] = 1;
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

  /* initialized flag_DoubleCounting */

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    flag_DoubleCounting[Gc_AN] = 0;
  }

  /* set Rcv_FGAN */

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    for (Lh_AN=0; Lh_AN<=FNAN[Gc_AN]; Lh_AN++){
      Gh_AN = natn[Gc_AN][Lh_AN];
      ID1 = G2ID[Gh_AN];
      if (flag_DoubleCounting[Gh_AN]==0 && ID1!=myid){

        Rcv_FGAN[ID1][F_Rcv_Num[ID1]] = Gh_AN;
        F_Rcv_Num[ID1]++;
        flag_DoubleCounting[Gh_AN] = 1;
      }
    }
  }

  /* set Rcv_SGAN */

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    for (Lh_AN=(FNAN[Gc_AN]+1); Lh_AN<=(FNAN[Gc_AN]+SNAN[Gc_AN]); Lh_AN++){
      Gh_AN = natn[Gc_AN][Lh_AN];
      ID1 = G2ID[Gh_AN];
      if (flag_DoubleCounting[Gh_AN]==0 && ID1!=myid){

        Rcv_SGAN[ID1][S_Rcv_Num[ID1]] = Gh_AN;
        S_Rcv_Num[ID1]++;
        flag_DoubleCounting[Gh_AN] = 1;
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
    global atom number to the intermediate atom number
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

  free(flag_DoubleCounting);
} 








void Construct_MPI_Data_Structure_Grid()
{
  static int firsttime=1;
  int i,j,k,Mc_AN,Gc_AN,wan,n1,n2,n3;
  int min_n1,max_n1,min_n2,max_n2,N3[4];
  unsigned long long int AN,BN,CN,DN;
  unsigned long long int B_AB2,Bs,Be;
  unsigned long long int BN_AB,BN_CB,BN_CA,GN_B_AB;
  unsigned long long int GN,GNs,GR,n2D,N2D;
  unsigned long long int GN_AB,GN_CB,GN_CA;
  int size_Index_Snd_Grid_A2B;
  int size_Index_Rcv_Grid_A2B;
  int size_Index_Snd_Grid_B2C;
  int size_Index_Rcv_Grid_B2C;
  int size_Index_Snd_Grid_B2D;
  int size_Index_Rcv_Grid_B2D;
  int size_Index_Snd_Grid_B_AB2CA;
  int size_Index_Rcv_Grid_B_AB2CA;
  int size_Index_Snd_Grid_B_CA2CB;
  int size_Index_Rcv_Grid_B_CA2CB;
  int myid,numprocs,ID,IDS,IDR,tag=999;
  double Vec0,Vec1,coef,MinV,MaxV,rcut;
  double Cxyz[4],b[4],c[4],v[4];

  MPI_Status stat;
  MPI_Request request;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /******************************************************
    find the smallest parallelepipedon which contains 
    atoms allocated to my ID under consideration of 
    cutoff radii of basis functions. 
  ******************************************************/
 
  for (k=1; k<=3; k++){

    if      (k==1){ i = 2; j = 3; }
    else if (k==2){ i = 3; j = 1; }
    else if (k==3){ i = 1; j = 2; }

    b[1] = tv[i][1];
    b[2] = tv[i][2];
    b[3] = tv[i][3];

    c[1] = tv[j][1];
    c[2] = tv[j][2];
    c[3] = tv[j][3];

    Cross_Product(b,c,v);
    coef = 1.0/sqrt(fabs( Dot_Product(v,v) ));
    
    v[1] = coef*v[1];
    v[2] = coef*v[2];
    v[3] = coef*v[3];

    MinV =  1.0e+10;
    MaxV = -1.0e+10;

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];
      wan  = WhatSpecies[Gc_AN];
      rcut = Spe_Atom_Cut1[wan];

      Cxyz[1] = Gxyz[Gc_AN][1] + rcut*v[1] - Grid_Origin[1];
      Cxyz[2] = Gxyz[Gc_AN][2] + rcut*v[2] - Grid_Origin[2];
      Cxyz[3] = Gxyz[Gc_AN][3] + rcut*v[3] - Grid_Origin[3];

      Vec0 = Dot_Product(Cxyz,rgtv[k])*0.5/PI;

      Cxyz[1] = Gxyz[Gc_AN][1] - rcut*v[1] - Grid_Origin[1];
      Cxyz[2] = Gxyz[Gc_AN][2] - rcut*v[2] - Grid_Origin[2];
      Cxyz[3] = Gxyz[Gc_AN][3] - rcut*v[3] - Grid_Origin[3];

      Vec1 = Dot_Product(Cxyz,rgtv[k])*0.5/PI;

      if (Vec0<MinV) MinV = Vec0;
      if (Vec1<MinV) MinV = Vec1;
      if (MaxV<Vec0) MaxV = Vec0;   
      if (MaxV<Vec1) MaxV = Vec1;
    }

    Min_Grid_Index[k] = (int)MinV;  /* buffer for GGA */
    Max_Grid_Index[k] = (int)MaxV;  /* buffer for GGA */

  } /* k */

  /******************************************************
    find the smallest parallelepipedon which contains 
    grids in the partition B_AB
    
    The parallelepipedon defines the partition D.
  ******************************************************/

  N2D = Ngrid1*Ngrid2;
  Bs = (myid*N2D+numprocs-1)/numprocs;
  Be = ((myid+1)*N2D+numprocs-1)/numprocs;
  
  min_n1 = 1000000;
  max_n1 =-1000000;
  min_n2 = 1000000;
  max_n2 =-1000000;

  for (B_AB2=Bs; B_AB2<Be; B_AB2++){

    n1 = B_AB2/Ngrid2;
    n2 = B_AB2 - n1*Ngrid2;

    if (n1<min_n1) min_n1 = n1;
    if (max_n1<n1) max_n1 = n1;
    if (n2<min_n2) min_n2 = n2;
    if (max_n2<n2) max_n2 = n2;
  }

  Min_Grid_Index_D[1] = min_n1 - 2;
  Max_Grid_Index_D[1] = max_n1 + 2;
  Min_Grid_Index_D[2] = min_n2 - 2;
  Max_Grid_Index_D[2] = max_n2 + 2;
  Min_Grid_Index_D[3] = -2;
  Max_Grid_Index_D[3] = (Ngrid3-1) + 2;

  /****************************************************************
      The partitions A to B

      construct the data structure for transfering rho_i from 
      the partitions A to B when rho is calculated 
      in the partition B using rho_i 
  ****************************************************************/

  /* find Num_Snd_Grid_A2B[ID] */

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_A2B[ID] = 0;

  N2D = Ngrid1*Ngrid2;

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];
    wan  = WhatSpecies[Gc_AN];

    for (AN=0; AN<GridN_Atom[Gc_AN]; AN++){

      GN = GridListAtom[Mc_AN][AN];

      /* get process ID and increment Num_Snd_Grid_A2B */

      GN2N(GN,N3);
      n2D = N3[1]*Ngrid2 + N3[2];
      ID = n2D*numprocs/N2D;
      Num_Snd_Grid_A2B[ID]++;
    }
  }    

  /* MPI: Num_Snd_Grid_A2B */  

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    tag = 999;
    MPI_Isend(&Num_Snd_Grid_A2B[IDS], 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
    MPI_Recv(&Num_Rcv_Grid_A2B[IDR], 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
    MPI_Wait(&request,&stat);
  }

  /* allocation of arrays */  

  if (alloc_first[26]==0){

    for (ID=0; ID<numprocs; ID++){
      free(Index_Snd_Grid_A2B[ID]);
    }  
    free(Index_Snd_Grid_A2B);
  
    for (ID=0; ID<numprocs; ID++){
      free(Index_Rcv_Grid_A2B[ID]);
    }  
    free(Index_Rcv_Grid_A2B);
  }

  size_Index_Snd_Grid_A2B = 0;
  Index_Snd_Grid_A2B = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Index_Snd_Grid_A2B[ID] = (int*)malloc(sizeof(int)*3*Num_Snd_Grid_A2B[ID]);
    size_Index_Snd_Grid_A2B += 3*Num_Snd_Grid_A2B[ID];
  }  
  
  size_Index_Rcv_Grid_A2B = 0; 
  Index_Rcv_Grid_A2B = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Index_Rcv_Grid_A2B[ID] = (int*)malloc(sizeof(int)*3*Num_Rcv_Grid_A2B[ID]);
    size_Index_Rcv_Grid_A2B += 3*Num_Rcv_Grid_A2B[ID]; 
  }  

  alloc_first[26] = 0;

  /* construct Index_Snd_Grid_A2B */  

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_A2B[ID] = 0;

  N2D = Ngrid1*Ngrid2;

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];
    wan  = WhatSpecies[Gc_AN];

    for (AN=0; AN<GridN_Atom[Gc_AN]; AN++){

      GN = GridListAtom[Mc_AN][AN];
      GR = CellListAtom[Mc_AN][AN];

      /* get process ID and grid index (BN) for the partition B */

      GN2N(GN,N3);
      n2D = N3[1]*Ngrid2 + N3[2];
      ID = n2D*numprocs/N2D;
      GN_B_AB = N3[1]*Ngrid2*Ngrid3 + N3[2]*Ngrid3 + N3[3];
      BN = GN_B_AB - ((ID*N2D+numprocs-1)/numprocs)*Ngrid3;

      /*
      printf("ABC1 myid=%2d ID=%2d\n",myid,ID);
      */

      Index_Snd_Grid_A2B[ID][3*Num_Snd_Grid_A2B[ID]+0] = BN; 
      Index_Snd_Grid_A2B[ID][3*Num_Snd_Grid_A2B[ID]+1] = Gc_AN; 
      Index_Snd_Grid_A2B[ID][3*Num_Snd_Grid_A2B[ID]+2] = GR; 

      Num_Snd_Grid_A2B[ID]++;
    }
  }    

  /* MPI: Index_Snd_Grid_A2B */
      
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    tag = 999;
   
    if (Num_Snd_Grid_A2B[IDS]!=0){
      MPI_Isend( &Index_Snd_Grid_A2B[IDS][0], 3*Num_Snd_Grid_A2B[IDS], 
                 MPI_INT, IDS, tag, mpi_comm_level1, &request);
    }

    if (Num_Rcv_Grid_A2B[IDR]!=0){
      MPI_Recv( &Index_Rcv_Grid_A2B[IDR][0], 3*Num_Rcv_Grid_A2B[IDR], 
                MPI_INT, IDR, tag, mpi_comm_level1, &stat);
    }

    if (Num_Snd_Grid_A2B[IDS]!=0)  MPI_Wait(&request,&stat);
  }

  /* count NN_A2B_S and NN_A2B_R */
  
  NN_A2B_S = 0;
  NN_A2B_R = 0;

  for (ID=0; ID<numprocs; ID++){
    if (Num_Snd_Grid_A2B[ID]!=0) NN_A2B_S++;
    if (Num_Rcv_Grid_A2B[ID]!=0) NN_A2B_R++;
  }

  /****************************************************************
      The partitions B to C

      construct the data structure for transfering rho from 
      the partitions B to C when rho is constructed in the 
      partition C using rho stored in the partition B.
  ****************************************************************/

  /* find Num_Rcv_Grid_B2C[ID] */

  for (ID=0; ID<numprocs; ID++) Num_Rcv_Grid_B2C[ID] = 0;

  N2D = Ngrid1*Ngrid2;

  for (n1=Min_Grid_Index[1]; n1<=Max_Grid_Index[1]; n1++){
    for (n2=Min_Grid_Index[2]; n2<=Max_Grid_Index[2]; n2++){
      for (n3=Min_Grid_Index[3]; n3<=Max_Grid_Index[3]; n3++){

         Find_CGrids(0,n1,n2,n3,Cxyz,N3);
         n2D = N3[1]*Ngrid2 + N3[2];
         ID = n2D*numprocs/N2D;
         Num_Rcv_Grid_B2C[ID]++;
      }
    }
  }

  /* MPI: Num_Rcv_Grid_B2C */  

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    tag = 999;
    MPI_Isend(&Num_Rcv_Grid_B2C[IDS], 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
    MPI_Recv(&Num_Snd_Grid_B2C[IDR], 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
    MPI_Wait(&request,&stat);
  }

  /* find Max_Num_Snd_Grid_B2C and Max_Num_Rcv_Grid_B2C */

  Max_Num_Snd_Grid_B2C = 0;
  Max_Num_Rcv_Grid_B2C = 0;
  for (ID=0; ID<numprocs; ID++){
    if ( Max_Num_Snd_Grid_B2C<Num_Snd_Grid_B2C[ID] ) Max_Num_Snd_Grid_B2C = Num_Snd_Grid_B2C[ID];
    if ( Max_Num_Rcv_Grid_B2C<Num_Rcv_Grid_B2C[ID] ) Max_Num_Rcv_Grid_B2C = Num_Rcv_Grid_B2C[ID];
  }

  /* allocation of arrays */  

  if (alloc_first[27]==0){
    for (ID=0; ID<numprocs; ID++){
      free(Index_Snd_Grid_B2C[ID]);
    }  
    free(Index_Snd_Grid_B2C);

    for (ID=0; ID<numprocs; ID++){
      free(Index_Rcv_Grid_B2C[ID]);
    }  
    free(Index_Rcv_Grid_B2C);

    free(ID_NN_B2C_S);
    free(ID_NN_B2C_R);
    free(GP_B2C_S);
    free(GP_B2C_R);
  }

  size_Index_Snd_Grid_B2C = 0;
  Index_Snd_Grid_B2C = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Index_Snd_Grid_B2C[ID] = (int*)malloc(sizeof(int)*Num_Snd_Grid_B2C[ID]);
    size_Index_Snd_Grid_B2C += Num_Snd_Grid_B2C[ID];
  }  

  size_Index_Rcv_Grid_B2C = 0;  
  Index_Rcv_Grid_B2C = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Index_Rcv_Grid_B2C[ID] = (int*)malloc(sizeof(int)*Num_Rcv_Grid_B2C[ID]);
    size_Index_Rcv_Grid_B2C += Num_Rcv_Grid_B2C[ID];
  }  
    
  alloc_first[27] = 0;

  /* construct Index_Rcv_Grid_B2C
     first BN is stored to Index_Rcv_Grid_B2C.  
     and after MPI communication, CN is stored 
     to Index_Rcv_Grid_B2C.  
     Also note that Index_Snd_Grid_B2C stores BN
  */

  for (ID=0; ID<numprocs; ID++) Num_Rcv_Grid_B2C[ID] = 0;

  N2D = Ngrid1*Ngrid2;

  for (n1=Min_Grid_Index[1]; n1<=Max_Grid_Index[1]; n1++){
    for (n2=Min_Grid_Index[2]; n2<=Max_Grid_Index[2]; n2++){
      for (n3=Min_Grid_Index[3]; n3<=Max_Grid_Index[3]; n3++){

         Find_CGrids(0,n1,n2,n3,Cxyz,N3);
         n2D = N3[1]*Ngrid2 + N3[2];
         ID = n2D*numprocs/N2D;
         GN = N3[1]*Ngrid2*Ngrid3 + N3[2]*Ngrid3 + N3[3];
         BN = GN - ((ID*N2D+numprocs-1)/numprocs)*Ngrid3;
         Index_Rcv_Grid_B2C[ID][Num_Rcv_Grid_B2C[ID]] = BN; 
         Num_Rcv_Grid_B2C[ID]++;
      }
    }
  }

  /* MPI: Index_Rcv_Grid_B2C */

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    tag = 999;

    if (Num_Rcv_Grid_B2C[IDS]!=0){
      MPI_Isend( &Index_Rcv_Grid_B2C[IDS][0], Num_Rcv_Grid_B2C[IDS],
                 MPI_INT, IDS, tag, mpi_comm_level1, &request);
    }

    if (Num_Snd_Grid_B2C[IDR]!=0){ 
      MPI_Recv( &Index_Snd_Grid_B2C[IDR][0], Num_Snd_Grid_B2C[IDR], 
                MPI_INT, IDR, tag, mpi_comm_level1, &stat);
    }

    if (Num_Rcv_Grid_B2C[IDS]!=0)  MPI_Wait(&request,&stat);
  }

  /* reconstruct Index_Rcv_Grid_B2C
     Index_Rcv_Grid_B2C stores CN.
  */

  for (ID=0; ID<numprocs; ID++) Num_Rcv_Grid_B2C[ID] = 0;

  N2D = Ngrid1*Ngrid2;
  CN = 0;
  for (n1=Min_Grid_Index[1]; n1<=Max_Grid_Index[1]; n1++){
    for (n2=Min_Grid_Index[2]; n2<=Max_Grid_Index[2]; n2++){
      for (n3=Min_Grid_Index[3]; n3<=Max_Grid_Index[3]; n3++){

         Find_CGrids(0,n1,n2,n3,Cxyz,N3);
         n2D = N3[1]*Ngrid2 + N3[2];
         ID = n2D*numprocs/N2D;
         Index_Rcv_Grid_B2C[ID][Num_Rcv_Grid_B2C[ID]] = CN; 
         Num_Rcv_Grid_B2C[ID]++;

         CN++;
      }
    }
  }

  /* find NN_B2C_S and NN_B2C_R
     and set ID_NN_B2C_S
             ID_NN_B2C_R
             GP_B2C_S
             GP_B2C_R 
  */ 

  NN_B2C_S = 0;
  NN_B2C_R = 0;

  for (ID=0; ID<numprocs; ID++){
    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;
    if (Num_Snd_Grid_B2C[IDS]!=0) NN_B2C_S++;
    if (Num_Rcv_Grid_B2C[IDR]!=0) NN_B2C_R++;
  }

  ID_NN_B2C_S = (int*)malloc(sizeof(int)*NN_B2C_S);
  ID_NN_B2C_R = (int*)malloc(sizeof(int)*NN_B2C_R);
  GP_B2C_S = (int*)malloc(sizeof(int)*(NN_B2C_S+1));
  GP_B2C_R = (int*)malloc(sizeof(int)*(NN_B2C_R+1));

  NN_B2C_S = 0;
  NN_B2C_R = 0;
  GP_B2C_S[0] = 0;
  GP_B2C_R[0] = 0;

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (Num_Snd_Grid_B2C[IDS]!=0){
      ID_NN_B2C_S[NN_B2C_S] = IDS;
      NN_B2C_S++;
      GP_B2C_S[NN_B2C_S] = GP_B2C_S[NN_B2C_S-1] + Num_Snd_Grid_B2C[IDS];
    }

    if (Num_Rcv_Grid_B2C[IDR]!=0){
      ID_NN_B2C_R[NN_B2C_R] = IDR;
      NN_B2C_R++;
      GP_B2C_R[NN_B2C_R] = GP_B2C_R[NN_B2C_R-1] + Num_Rcv_Grid_B2C[IDR];
    }
  }

  /* set the number of grids allocated to the partitions B and C for myid */

  N2D = Ngrid1*Ngrid2;
  My_NumGridB_AB =  (((myid+1)*N2D+numprocs-1)/numprocs)*Ngrid3
                  - ((myid*N2D+numprocs-1)/numprocs)*Ngrid3;
  My_NumGridC = CN;

  /****************************************************************
      The partitions B to D

      construct the data structure for transfering rho from 
      the partitions B to D for the case that rho is constructed
      in the partition D using rho stored in the partition B.
  ****************************************************************/

  /* find Num_Rcv_Grid_B2D[ID] */

  for (ID=0; ID<numprocs; ID++) Num_Rcv_Grid_B2D[ID] = 0;

  N2D = Ngrid1*Ngrid2;

  for (n1=Min_Grid_Index_D[1]; n1<=Max_Grid_Index_D[1]; n1++){
    for (n2=Min_Grid_Index_D[2]; n2<=Max_Grid_Index_D[2]; n2++){
      for (n3=Min_Grid_Index_D[3]; n3<=Max_Grid_Index_D[3]; n3++){

         Find_CGrids(0,n1,n2,n3,Cxyz,N3);
         n2D = N3[1]*Ngrid2 + N3[2];
         ID = n2D*numprocs/N2D;
         Num_Rcv_Grid_B2D[ID]++;
      }
    }
  }

  /* MPI: Num_Rcv_Grid_B2D */  

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    tag = 999;
    MPI_Isend(&Num_Rcv_Grid_B2D[IDS], 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
    MPI_Recv(&Num_Snd_Grid_B2D[IDR], 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
    MPI_Wait(&request,&stat);
  }

  /* find Max_Num_Snd_Grid_B2D and Max_Num_Rcv_Grid_B2D */

  Max_Num_Snd_Grid_B2D = 0;
  Max_Num_Rcv_Grid_B2D = 0;
  for (ID=0; ID<numprocs; ID++){
    if ( Max_Num_Snd_Grid_B2D<Num_Snd_Grid_B2D[ID] ) Max_Num_Snd_Grid_B2D = Num_Snd_Grid_B2D[ID];
    if ( Max_Num_Rcv_Grid_B2D<Num_Rcv_Grid_B2D[ID] ) Max_Num_Rcv_Grid_B2D = Num_Rcv_Grid_B2D[ID];
  }

  /* allocation of arrays */  

  if (alloc_first[30]==0){
    for (ID=0; ID<numprocs; ID++){
      free(Index_Snd_Grid_B2D[ID]);
    }  
    free(Index_Snd_Grid_B2D);

    for (ID=0; ID<numprocs; ID++){
      free(Index_Rcv_Grid_B2D[ID]);
    }  
    free(Index_Rcv_Grid_B2D);

    free(ID_NN_B2D_S);
    free(ID_NN_B2D_R);
    free(GP_B2D_S);
    free(GP_B2D_R);
  }

  size_Index_Snd_Grid_B2D = 0;
  Index_Snd_Grid_B2D = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Index_Snd_Grid_B2D[ID] = (int*)malloc(sizeof(int)*Num_Snd_Grid_B2D[ID]);
    size_Index_Snd_Grid_B2D += Num_Snd_Grid_B2D[ID];
  }  

  size_Index_Rcv_Grid_B2D = 0;  
  Index_Rcv_Grid_B2D = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Index_Rcv_Grid_B2D[ID] = (int*)malloc(sizeof(int)*Num_Rcv_Grid_B2D[ID]);
    size_Index_Rcv_Grid_B2D += Num_Rcv_Grid_B2D[ID];
  }  
    
  alloc_first[30] = 0;

  /* construct Index_Rcv_Grid_B2D
     first BN is stored to Index_Rcv_Grid_B2D.  
     and after MPI communication, DN is stored 
     to Index_Rcv_Grid_B2D.  
     Also note that Index_Snd_Grid_B2D stores BN
  */

  for (ID=0; ID<numprocs; ID++) Num_Rcv_Grid_B2D[ID] = 0;
  
  N2D = Ngrid1*Ngrid2;
  
  for (n1=Min_Grid_Index_D[1]; n1<=Max_Grid_Index_D[1]; n1++){
    for (n2=Min_Grid_Index_D[2]; n2<=Max_Grid_Index_D[2]; n2++){
      for (n3=Min_Grid_Index_D[3]; n3<=Max_Grid_Index_D[3]; n3++){
  
         Find_CGrids(0,n1,n2,n3,Cxyz,N3);
         n2D = N3[1]*Ngrid2 + N3[2];
         ID = n2D*numprocs/N2D;
         GN = N3[1]*Ngrid2*Ngrid3 + N3[2]*Ngrid3 + N3[3];
         BN = GN - ((ID*N2D+numprocs-1)/numprocs)*Ngrid3;
         Index_Rcv_Grid_B2D[ID][Num_Rcv_Grid_B2D[ID]] = BN; 
         Num_Rcv_Grid_B2D[ID]++;
      }
    }
  }

  /* MPI: Index_Rcv_Grid_B2D */

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    tag = 999;

    if (Num_Rcv_Grid_B2D[IDS]!=0){
      MPI_Isend( &Index_Rcv_Grid_B2D[IDS][0], Num_Rcv_Grid_B2D[IDS],
                 MPI_INT, IDS, tag, mpi_comm_level1, &request);
    }

    if (Num_Snd_Grid_B2D[IDR]!=0){ 
      MPI_Recv( &Index_Snd_Grid_B2D[IDR][0], Num_Snd_Grid_B2D[IDR], 
                MPI_INT, IDR, tag, mpi_comm_level1, &stat);
    }

    if (Num_Rcv_Grid_B2D[IDS]!=0)  MPI_Wait(&request,&stat);
  }

  /* reconstruct Index_Rcv_Grid_B2D
     Index_Rcv_Grid_B2D stores DN.
  */

  for (ID=0; ID<numprocs; ID++) Num_Rcv_Grid_B2D[ID] = 0;

  N2D = Ngrid1*Ngrid2;
  DN = 0;
  for (n1=Min_Grid_Index_D[1]; n1<=Max_Grid_Index_D[1]; n1++){
    for (n2=Min_Grid_Index_D[2]; n2<=Max_Grid_Index_D[2]; n2++){
      for (n3=Min_Grid_Index_D[3]; n3<=Max_Grid_Index_D[3]; n3++){

         Find_CGrids(0,n1,n2,n3,Cxyz,N3);
         n2D = N3[1]*Ngrid2 + N3[2];
         ID = n2D*numprocs/N2D;
         Index_Rcv_Grid_B2D[ID][Num_Rcv_Grid_B2D[ID]] = DN; 
         Num_Rcv_Grid_B2D[ID]++;

         DN++;
      }
    }
  }

  /* find NN_B2C_S and NN_B2C_R
     and set ID_NN_B2C_S
             ID_NN_B2C_R
             GP_B2C_S
             GP_B2C_R 
  */ 

  NN_B2D_S = 0;
  NN_B2D_R = 0;

  for (ID=0; ID<numprocs; ID++){
    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;
    if (Num_Snd_Grid_B2D[IDS]!=0) NN_B2D_S++;
    if (Num_Rcv_Grid_B2D[IDR]!=0) NN_B2D_R++;
  }

  ID_NN_B2D_S = (int*)malloc(sizeof(int)*NN_B2D_S);
  ID_NN_B2D_R = (int*)malloc(sizeof(int)*NN_B2D_R);
  GP_B2D_S = (int*)malloc(sizeof(int)*(NN_B2D_S+1));
  GP_B2D_R = (int*)malloc(sizeof(int)*(NN_B2D_R+1));

  NN_B2D_S = 0;
  NN_B2D_R = 0;
  GP_B2D_S[0] = 0;
  GP_B2D_R[0] = 0;

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (Num_Snd_Grid_B2D[IDS]!=0){
      ID_NN_B2D_S[NN_B2D_S] = IDS;
      NN_B2D_S++;
      GP_B2D_S[NN_B2D_S] = GP_B2D_S[NN_B2D_S-1] + Num_Snd_Grid_B2D[IDS];
    }

    if (Num_Rcv_Grid_B2D[IDR]!=0){
      ID_NN_B2D_R[NN_B2D_R] = IDR;
      NN_B2D_R++;
      GP_B2D_R[NN_B2D_R] = GP_B2D_R[NN_B2D_R-1] + Num_Rcv_Grid_B2D[IDR];
    }
  }

  /* set the number of grids allocated to the D for myid */

  My_NumGridD = DN;

  /****************************************************************
      AB to CA in the partition B for MPI communication in FFT 
  ****************************************************************/

  /* set GNs */

  N2D = Ngrid1*Ngrid2;
  GNs = ((myid*N2D+numprocs-1)/numprocs)*Ngrid3;

  /* find Num_Snd_Grid_B_AB2CA */

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_B_AB2CA[ID] = 0;

  N2D = Ngrid3*Ngrid1; /* for CA */

  for (BN_AB=0; BN_AB<My_NumGridB_AB; BN_AB++){

    GN_AB = BN_AB + GNs;

    n1 = GN_AB/(Ngrid2*Ngrid3);
    n2 = (GN_AB - n1*(Ngrid2*Ngrid3))/Ngrid3;
    n3 = GN_AB - n1*(Ngrid2*Ngrid3) - n2*Ngrid3;
    n2D = n3*Ngrid1 + n1;
    ID = n2D*numprocs/N2D;
    Num_Snd_Grid_B_AB2CA[ID]++;
  }

  /* MPI: Num_Snd_Grid_B_AB2CA */  

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    tag = 999;
    MPI_Isend(&Num_Snd_Grid_B_AB2CA[IDS], 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
    MPI_Recv(&Num_Rcv_Grid_B_AB2CA[IDR], 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
    MPI_Wait(&request,&stat);
  }

  /* allocation of arrays */  

  if (alloc_first[28]==0){

    for (ID=0; ID<numprocs; ID++){
      free(Index_Snd_Grid_B_AB2CA[ID]);
    }  
    free(Index_Snd_Grid_B_AB2CA);
  
    for (ID=0; ID<numprocs; ID++){
      free(Index_Rcv_Grid_B_AB2CA[ID]);
    }  
    free(Index_Rcv_Grid_B_AB2CA);

    free(ID_NN_B_AB2CA_S);
    free(ID_NN_B_AB2CA_R);
    free(GP_B_AB2CA_S);
    free(GP_B_AB2CA_R);
  }

  size_Index_Snd_Grid_B_AB2CA = 0;
  Index_Snd_Grid_B_AB2CA = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Index_Snd_Grid_B_AB2CA[ID] = (int*)malloc(sizeof(int)*Num_Snd_Grid_B_AB2CA[ID]);
    size_Index_Snd_Grid_B_AB2CA += Num_Snd_Grid_B_AB2CA[ID];
  }  
  
  size_Index_Rcv_Grid_B_AB2CA = 0; 
  Index_Rcv_Grid_B_AB2CA = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Index_Rcv_Grid_B_AB2CA[ID] = (int*)malloc(sizeof(int)*Num_Rcv_Grid_B_AB2CA[ID]);
    size_Index_Rcv_Grid_B_AB2CA += Num_Rcv_Grid_B_AB2CA[ID]; 
  }  

  alloc_first[28] = 0;

  /* construct Index_Snd_Grid_B_AB2CA */  

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_B_AB2CA[ID] = 0;

  N2D = Ngrid3*Ngrid1;  /* for CA */

  for (BN_AB=0; BN_AB<My_NumGridB_AB; BN_AB++){

    GN_AB = BN_AB + GNs;

    n1 = GN_AB/(Ngrid2*Ngrid3);
    n2 = (GN_AB - n1*(Ngrid2*Ngrid3))/Ngrid3;
    n3 = GN_AB - n1*(Ngrid2*Ngrid3) - n2*Ngrid3;
    n2D = n3*Ngrid1 + n1;
    ID = n2D*numprocs/N2D;
    GN_CA = n3*Ngrid1*Ngrid2 + n1*Ngrid2 + n2;
    BN_CA = GN_CA - ((ID*N2D+numprocs-1)/numprocs)*Ngrid2;
    Index_Snd_Grid_B_AB2CA[ID][Num_Snd_Grid_B_AB2CA[ID]] = BN_CA;
    Num_Snd_Grid_B_AB2CA[ID]++;
  }

  /* MPI: Index_Snd_Grid_B_AB2CA */

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    tag = 999;
   
    if (Num_Snd_Grid_B_AB2CA[IDS]!=0){
      MPI_Isend( &Index_Snd_Grid_B_AB2CA[IDS][0], Num_Snd_Grid_B_AB2CA[IDS], 
		 MPI_INT, IDS, tag, mpi_comm_level1, &request);
    }

    if (Num_Rcv_Grid_B_AB2CA[IDR]!=0){
      MPI_Recv( &Index_Rcv_Grid_B_AB2CA[IDR][0], Num_Rcv_Grid_B_AB2CA[IDR], 
		MPI_INT, IDR, tag, mpi_comm_level1, &stat);
    }

    if (Num_Snd_Grid_B_AB2CA[IDS]!=0)  MPI_Wait(&request,&stat);
  }

  /* reset: Index_Snd_Grid_B_AB2CA

  Index_Snd_Grid_B_AB2CA:  BN_AB
  Index_Rcv_Grid_B_AB2CA:  BN_CA
  */

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_B_AB2CA[ID] = 0;

  N2D = Ngrid3*Ngrid1; /* for CA */

  for (BN_AB=0; BN_AB<My_NumGridB_AB; BN_AB++){

    GN_AB = BN_AB + GNs;
    n1 = GN_AB/(Ngrid2*Ngrid3);
    n2 = (GN_AB - n1*(Ngrid2*Ngrid3))/Ngrid3;
    n3 = GN_AB - n1*(Ngrid2*Ngrid3) - n2*Ngrid3;
    n2D = n3*Ngrid1 + n1;
    ID = n2D*numprocs/N2D;
    Index_Snd_Grid_B_AB2CA[ID][Num_Snd_Grid_B_AB2CA[ID]] = BN_AB;
    Num_Snd_Grid_B_AB2CA[ID]++;
  }

  /* find the maximum Num_Snd_Grid_B_AB2CA and Num_Rcv_Grid_B_AB2CA */

  Max_Num_Snd_Grid_B_AB2CA = 0;
  Max_Num_Rcv_Grid_B_AB2CA = 0;

  for (ID=0; ID<numprocs; ID++){
    if (Max_Num_Snd_Grid_B_AB2CA<Num_Snd_Grid_B_AB2CA[ID]){ 
      Max_Num_Snd_Grid_B_AB2CA = Num_Snd_Grid_B_AB2CA[ID];
    }

    if (Max_Num_Rcv_Grid_B_AB2CA<Num_Rcv_Grid_B_AB2CA[ID]){ 
      Max_Num_Rcv_Grid_B_AB2CA = Num_Rcv_Grid_B_AB2CA[ID];
    }
  }

  /* find NN_B_AB2CA_S and NN_B_AB2CA_R 
     and set ID_NN_B_AB2CA_S, 
             ID_NN_B_AB2CA_R,
             GP_B_AB2CA_S,
             GP_B_AB2CA_R
  */

  NN_B_AB2CA_S = 0;
  NN_B_AB2CA_R = 0;
 
  for (ID=0; ID<numprocs; ID++){
    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;
    if (Num_Snd_Grid_B_AB2CA[IDS]!=0) NN_B_AB2CA_S++;
    if (Num_Rcv_Grid_B_AB2CA[IDR]!=0) NN_B_AB2CA_R++;
  }

  ID_NN_B_AB2CA_S = (int*)malloc(sizeof(int)*NN_B_AB2CA_S);
  ID_NN_B_AB2CA_R = (int*)malloc(sizeof(int)*NN_B_AB2CA_R);
  GP_B_AB2CA_S = (int*)malloc(sizeof(int)*(NN_B_AB2CA_S+1));
  GP_B_AB2CA_R = (int*)malloc(sizeof(int)*(NN_B_AB2CA_R+1));

  NN_B_AB2CA_S = 0;
  NN_B_AB2CA_R = 0;
  GP_B_AB2CA_S[0] = 0;
  GP_B_AB2CA_R[0] = 0;

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (Num_Snd_Grid_B_AB2CA[IDS]!=0){
      ID_NN_B_AB2CA_S[NN_B_AB2CA_S] = IDS;
      NN_B_AB2CA_S++;
      GP_B_AB2CA_S[NN_B_AB2CA_S] = GP_B_AB2CA_S[NN_B_AB2CA_S-1] + Num_Snd_Grid_B_AB2CA[IDS];
    }

    if (Num_Rcv_Grid_B_AB2CA[IDR]!=0){
      ID_NN_B_AB2CA_R[NN_B_AB2CA_R] = IDR;
      NN_B_AB2CA_R++;
      GP_B_AB2CA_R[NN_B_AB2CA_R] = GP_B_AB2CA_R[NN_B_AB2CA_R-1] + Num_Rcv_Grid_B_AB2CA[IDR];
    }
  }

  /* set My_NumGridB_CA */

  N2D = Ngrid3*Ngrid1;
  My_NumGridB_CA = (((myid+1)*N2D+numprocs-1)/numprocs)*Ngrid2
                  - ((myid*N2D+numprocs-1)/numprocs)*Ngrid2;

  /****************************************************************
      CA to CB in the partition B for MPI communication in FFT 
  ****************************************************************/

  /* set GNs */

  N2D = Ngrid3*Ngrid1;
  GNs = ((myid*N2D+numprocs-1)/numprocs)*Ngrid2;

  /* find Num_Snd_Grid_B_CA2CB */

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_B_CA2CB[ID] = 0;

  N2D = Ngrid3*Ngrid2; /* for CB */

  for (BN_CA=0; BN_CA<My_NumGridB_CA; BN_CA++){

    GN_CA = GNs + BN_CA;    
    n3 = GN_CA/(Ngrid1*Ngrid2);
    n1 = (GN_CA - n3*(Ngrid1*Ngrid2))/Ngrid2;
    n2 = GN_CA - n3*(Ngrid1*Ngrid2) - n1*Ngrid2;
    n2D = n3*Ngrid2 + n2;
    ID = n2D*numprocs/N2D;
    Num_Snd_Grid_B_CA2CB[ID]++;
  }

  /* MPI: Num_Snd_Grid_B_CA2CB */  

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    tag = 999;
    MPI_Isend(&Num_Snd_Grid_B_CA2CB[IDS], 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
    MPI_Recv(&Num_Rcv_Grid_B_CA2CB[IDR], 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
    MPI_Wait(&request,&stat);
  }

  /* allocation of arrays */  

  if (alloc_first[29]==0){

    for (ID=0; ID<numprocs; ID++){
      free(Index_Snd_Grid_B_CA2CB[ID]);
    }  
    free(Index_Snd_Grid_B_CA2CB);
  
    for (ID=0; ID<numprocs; ID++){
      free(Index_Rcv_Grid_B_CA2CB[ID]);
    }  
    free(Index_Rcv_Grid_B_CA2CB);

    free(ID_NN_B_CA2CB_S);
    free(ID_NN_B_CA2CB_R);
    free(GP_B_CA2CB_S);
    free(GP_B_CA2CB_R);
  }

  size_Index_Snd_Grid_B_CA2CB = 0;
  Index_Snd_Grid_B_CA2CB = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Index_Snd_Grid_B_CA2CB[ID] = (int*)malloc(sizeof(int)*Num_Snd_Grid_B_CA2CB[ID]);
    size_Index_Snd_Grid_B_CA2CB += Num_Snd_Grid_B_CA2CB[ID];
  }  
  
  size_Index_Rcv_Grid_B_CA2CB = 0; 
  Index_Rcv_Grid_B_CA2CB = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Index_Rcv_Grid_B_CA2CB[ID] = (int*)malloc(sizeof(int)*Num_Rcv_Grid_B_CA2CB[ID]);
    size_Index_Rcv_Grid_B_CA2CB += Num_Rcv_Grid_B_CA2CB[ID]; 
  }  

  alloc_first[29] = 0;

  /* construct Index_Snd_Grid_B_CA2CB */  

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_B_CA2CB[ID] = 0;

  N2D = Ngrid3*Ngrid2; /* for CB */

  for (BN_CA=0; BN_CA<My_NumGridB_CA; BN_CA++){

    GN_CA = GNs + BN_CA;    
    n3 = GN_CA/(Ngrid1*Ngrid2);
    n1 = (GN_CA - n3*(Ngrid1*Ngrid2))/Ngrid2;
    n2 = GN_CA - n3*(Ngrid1*Ngrid2) - n1*Ngrid2;
    n2D = n3*Ngrid2 + n2;
    ID = n2D*numprocs/N2D;
    GN_CB = n3*Ngrid2*Ngrid1 + n2*Ngrid1 + n1;
    BN_CB = GN_CB - ((ID*N2D+numprocs-1)/numprocs)*Ngrid1;
    Index_Snd_Grid_B_CA2CB[ID][Num_Snd_Grid_B_CA2CB[ID]] = BN_CB;
    Num_Snd_Grid_B_CA2CB[ID]++;
  }

  /* MPI: Index_Snd_Grid_B_CA2CB */

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    tag = 999;
   
    if (Num_Snd_Grid_B_CA2CB[IDS]!=0){
      MPI_Isend( &Index_Snd_Grid_B_CA2CB[IDS][0], Num_Snd_Grid_B_CA2CB[IDS], 
                 MPI_INT, IDS, tag, mpi_comm_level1, &request);
    }

    if (Num_Rcv_Grid_B_CA2CB[IDR]!=0){
      MPI_Recv( &Index_Rcv_Grid_B_CA2CB[IDR][0], Num_Rcv_Grid_B_CA2CB[IDR], 
                MPI_INT, IDR, tag, mpi_comm_level1, &stat);
    }

    if (Num_Snd_Grid_B_CA2CB[IDS]!=0)  MPI_Wait(&request,&stat);
  }

  /* reset: Index_Snd_Grid_B_CA2CB

     Index_Snd_Grid_B_CB2CA:  BN_CA
     Index_Rcv_Grid_B_CB2CA:  BN_CB
  */

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_B_CA2CB[ID] = 0;

  N2D = Ngrid3*Ngrid2; /* for CB */

  for (BN_CA=0; BN_CA<My_NumGridB_CA; BN_CA++){

    GN_CA = GNs + BN_CA;    
    n3 = GN_CA/(Ngrid1*Ngrid2);
    n1 = (GN_CA - n3*(Ngrid1*Ngrid2))/Ngrid2;
    n2 = GN_CA - n3*(Ngrid1*Ngrid2) - n1*Ngrid2;
    n2D = n3*Ngrid2 + n2;
    ID = n2D*numprocs/N2D;
    Index_Snd_Grid_B_CA2CB[ID][Num_Snd_Grid_B_CA2CB[ID]] = BN_CA;
    Num_Snd_Grid_B_CA2CB[ID]++;
  }

  /* find the maximum Num_Snd_Grid_B_CA2CB and Num_Rcv_Grid_B_CA2CB */

  Max_Num_Snd_Grid_B_CA2CB = 0;
  Max_Num_Rcv_Grid_B_CA2CB = 0;

  for (ID=0; ID<numprocs; ID++){
    if (Max_Num_Snd_Grid_B_CA2CB<Num_Snd_Grid_B_CA2CB[ID]){ 
      Max_Num_Snd_Grid_B_CA2CB = Num_Snd_Grid_B_CA2CB[ID];
    }

    if (Max_Num_Rcv_Grid_B_CA2CB<Num_Rcv_Grid_B_CA2CB[ID]){ 
      Max_Num_Rcv_Grid_B_CA2CB = Num_Rcv_Grid_B_CA2CB[ID];
    }
  }

  /* find NN_B_CA2CB_S and NN_B_CA2CB_R 
     and set ID_NN_B_CA2CB_S, 
             ID_NN_B_CA2CB_R,
             GP_B_CA2CB_S,
             GP_B_CA2CB_R
  */

  NN_B_CA2CB_S = 0;
  NN_B_CA2CB_R = 0;

  for (ID=0; ID<numprocs; ID++){
    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;
    if (Num_Snd_Grid_B_CA2CB[IDS]!=0) NN_B_CA2CB_S++;
    if (Num_Rcv_Grid_B_CA2CB[IDR]!=0) NN_B_CA2CB_R++;
  }

  ID_NN_B_CA2CB_S = (int*)malloc(sizeof(int)*NN_B_CA2CB_S);
  ID_NN_B_CA2CB_R = (int*)malloc(sizeof(int)*NN_B_CA2CB_R);
  GP_B_CA2CB_S = (int*)malloc(sizeof(int)*(NN_B_CA2CB_S+1));
  GP_B_CA2CB_R = (int*)malloc(sizeof(int)*(NN_B_CA2CB_R+1));

  NN_B_CA2CB_S = 0;
  NN_B_CA2CB_R = 0;
  GP_B_CA2CB_S[0] = 0;
  GP_B_CA2CB_R[0] = 0;

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (Num_Snd_Grid_B_CA2CB[IDS]!=0){
      ID_NN_B_CA2CB_S[NN_B_CA2CB_S] = IDS;
      NN_B_CA2CB_S++;
      GP_B_CA2CB_S[NN_B_CA2CB_S] = GP_B_CA2CB_S[NN_B_CA2CB_S-1] + Num_Snd_Grid_B_CA2CB[IDS];
    }

    if (Num_Rcv_Grid_B_CA2CB[IDR]!=0){
      ID_NN_B_CA2CB_R[NN_B_CA2CB_R] = IDR;
      NN_B_CA2CB_R++;
      GP_B_CA2CB_R[NN_B_CA2CB_R] = GP_B_CA2CB_R[NN_B_CA2CB_R-1] + Num_Rcv_Grid_B_CA2CB[IDR];
    }
  }

  /* set My_NumGridB_CB */

  N2D = Ngrid3*Ngrid2;
  My_NumGridB_CB = (((myid+1)*N2D+numprocs-1)/numprocs)*Ngrid1
                  - ((myid*N2D+numprocs-1)/numprocs)*Ngrid1;

  /* find My_Max_NumGridB */

  My_Max_NumGridB = 0;
  if (My_Max_NumGridB<My_NumGridB_AB) My_Max_NumGridB = My_NumGridB_AB;
  if (My_Max_NumGridB<My_NumGridB_CA) My_Max_NumGridB = My_NumGridB_CA;
  if (My_Max_NumGridB<My_NumGridB_CB) My_Max_NumGridB = My_NumGridB_CB;

  /* PrintMemory */
  if (firsttime){
  PrintMemory("truncation: Index_Snd_Grid_A2B",     sizeof(int)*size_Index_Snd_Grid_A2B,  NULL);
  PrintMemory("truncation: Index_Rcv_Grid_A2B",     sizeof(int)*size_Index_Rcv_Grid_A2B,  NULL);
  PrintMemory("truncation: Index_Snd_Grid_B2C",     sizeof(int)*size_Index_Snd_Grid_B2C,  NULL);
  PrintMemory("truncation: Index_Rcv_Grid_B2C",     sizeof(int)*size_Index_Rcv_Grid_B2C,  NULL);
  PrintMemory("truncation: Index_Snd_Grid_B2D",     sizeof(int)*size_Index_Snd_Grid_B2D,  NULL);
  PrintMemory("truncation: Index_Rcv_Grid_B2D",     sizeof(int)*size_Index_Rcv_Grid_B2D,  NULL);
  PrintMemory("truncation: Index_Snd_Grid_B_AB2CA", sizeof(int)*size_Index_Snd_Grid_B_AB2CA, NULL);
  PrintMemory("truncation: Index_Rcv_Grid_B_AB2CA", sizeof(int)*size_Index_Rcv_Grid_B_AB2CA, NULL);
  PrintMemory("truncation: Index_Snd_Grid_B_CA2CB", sizeof(int)*size_Index_Snd_Grid_B_CA2CB, NULL);
  PrintMemory("truncation: Index_Rcv_Grid_B_CA2CB", sizeof(int)*size_Index_Rcv_Grid_B_CA2CB, NULL);
  PrintMemory("truncation: ID_NN_B_AB2CA_S", sizeof(int)*NN_B_AB2CA_S, NULL);
  PrintMemory("truncation: ID_NN_B_AB2CA_R", sizeof(int)*NN_B_AB2CA_R, NULL);
  PrintMemory("truncation: ID_NN_B_CA2CB_S", sizeof(int)*NN_B_CA2CB_S, NULL);
  PrintMemory("truncation: ID_NN_B_CA2CB_R", sizeof(int)*NN_B_CA2CB_R, NULL);
  PrintMemory("truncation: GP_B_AB2CA_S", sizeof(int)*(NN_B_AB2CA_S+1), NULL);
  PrintMemory("truncation: GP_B_AB2CA_R", sizeof(int)*(NN_B_AB2CA_R+1), NULL);
  PrintMemory("truncation: GP_B_CA2CB_S", sizeof(int)*(NN_B_CA2CB_S+1), NULL);
  PrintMemory("truncation: GP_B_CA2CB_R", sizeof(int)*(NN_B_CA2CB_R+1), NULL);
  firsttime = 0;
  }
}

