/**********************************************************************
  Cluster_DFT_ON2.c:

     Cluster_DFT_ON2.c is a subroutine to perform cluster calculations
     with an O(N^2) method.

  Log of Cluster_DFT_ON2.c:

     11/Sep./2009  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "openmx_common.h"
#include "lapack_prototypes.h"
#include "cs.h"
#include <metis.h>

#define  measure_time       2
#define  sparse_solver      3    /* 1: CXSparse, 2: WSMP, 3: OND */

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

typedef struct {
   long int a;
   double complex b;
} clists;

static void qsort_complex(long n, long *a, double complex *b);
static int clists_cmp(const clists *x, const clists *y);
int Lapack_LU_Zinverse(int , dcomplex *);

static void Partition_System(
    int sizeflag, int *nl0, int *mpmax0,
    int Nmat, int TNZE2, 
    idxtype *xadj, idxtype *adjncy, 
    int **Index_BC,
    int **SindB,int **EindB,
    int **SindC,int **EindC);

static double Cluster_collinear_ON2(
				    char *mode,
				    int SCF_iter,
				    int SpinP_switch,
				    double ***C,
				    double **ko,
				    double *****nh, double ****CntOLP,
				    double *****CDM,
				    double *****EDM,
				    double Eele0[2], double Eele1[2]);

static void OND_Solver(
  int myid,
  int numprocs,
  int Nloop1,
  MPI_Comm MPI_Comm,
  int ChemP_point,
  int spin,
  int Epoint,
  double ChemP_trial[3],
  int Nmat,
  int *NZE,
  int *PZE,
  double **Hall,
  double *Sall,
  double ***DMall,
  int *conv_index_row,
  int *conv_index_col,
  int **conv_ind2,
  int TNZE);



static void Get_EigenStates_Near_pChemP(
  /* input */
  int myid0,
  int numprocs0, 
  int myworld1,
  int Num_Comm_World1, 
  int *NPROCS_ID1,
  int *NPROCS_WD1,
  int *Comm_World1, 
  int *Comm_World_StartID1,
  MPI_Comm *MPI_CommWD1,
  int NEV, 
  int Nmat, 
  double pChemP, 
  int *NZE,
  int *PZE,
  double **Hall,
  double *Sall,
  double **DMall,
  int *conv_index_row,
  int *conv_index_col,
  int **conv_ind2,
  int TNZE, 
  /* output */
  double **eval, 
  double **evec);



double Cluster_DFT_ON2(char *mode,
		       int SCF_iter,
		       int SpinP_switch,
		       double ***Cluster_ReCoes,
		       double **Cluster_ko,
		       double *****nh,
		       double *****ImNL,
		       double ****CntOLP,
		       double *****CDM,
		       double *****EDM,
		       double Eele0[2], double Eele1[2])
{
  static double time0;

  /****************************************************
         collinear without spin-orbit coupling
  ****************************************************/

  if ( (SpinP_switch==0 || SpinP_switch==1) && SO_switch==0 ){

    time0 = Cluster_collinear_ON2(mode,SCF_iter,SpinP_switch,Cluster_ReCoes,Cluster_ko,
                                  nh,CntOLP,CDM,EDM,Eele0,Eele1);

  }

  /****************************************************
           collinear with spin-orbit coupling
  ****************************************************/

  else if ( (SpinP_switch==0 || SpinP_switch==1) && SO_switch==1 ){
    printf("Spin-orbit coupling is not supported for collinear DFT calculations.\n");
    MPI_Finalize();
    exit(1);
  }

  /****************************************************
   non-collinear with and without spin-orbit coupling
  ****************************************************/

  else if (SpinP_switch==3){
    /*
    time0 = Cluster_non_collinear(mode,SCF_iter,SpinP_switch,nh,ImNL,CntOLP,CDM,EDM,Eele0,Eele1);
    */
  }

  return time0;
}



static double Cluster_collinear_ON2(
				    char *mode,
				    int SCF_iter,
				    int SpinP_switch,
				    double ***C,
				    double **ko,
				    double *****nh, double ****CntOLP,
				    double *****CDM,
				    double *****EDM,
				    double Eele0[2], double Eele1[2])
{
  static int firsttime=1;
  double TStime,TEtime;
  double stime,etime;
  int numprocs,myid,ID,tag=999;
  int i,j,k,k2,wan,wanA,wanB,Epoint,column;
  int tnoA,tnoB,LB_AN,i1,j1,l,ChemP_point;
  int MA_AN,GA_AN,GB_AN,Gc_AN,h_AN,ni,nj,TNZE;
  int Cwan,Hwan,Gh_AN,Anum,tnum,num;
  int loop1,Nloop1,Nmat,spin,NEV;
  int *NZE,*PZE,*MP;
  double **Hall,*Sall;
  double ***DMall;
  double **eval,**evec;
  double a,b,c,x0,x1,x2,y0,y1,y2,d;
  double xopt,weight_DM[3];
  double eachsum[3],totalsum[3],pChemP; 
  double TZ,ChemP_MIN,ChemP_MAX;
  double Dnum,spin_degeneracy;
  double f0,f1,totalsum1;
  double My_Eele1[2]; 
  double ChemP_trial[3];
  int po,loopN; 
  int *conv_index_row;
  int *conv_index_col;
  int **conv_ind2;
  int *My_NZeros;
  int *is1,*ie1;
  /* for MPI */
  int myworld1;
  int numprocs1,myid1;
  int Num_Comm_World1;
  int *NPROCS_ID1,*NPROCS_WD1;
  int *Comm_World1;
  int *Comm_World_StartID1;
  MPI_Comm *MPI_CommWD1;

  int myworld2B;
  int numprocs2B,myid2B;
  int Num_Comm_World2B;
  int *NPROCS_ID2B,*NPROCS_WD2B;
  int *Comm_World2B;
  int *Comm_World_StartID2B;
  MPI_Comm *MPI_CommWD2B;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /* for elapsed time */
  MPI_Barrier(mpi_comm_level1);
  dtime(&TStime);

  if (0<measure_time) dtime(&stime);

  /****************************************************
         set up the first communication world
  ****************************************************/

  Num_Comm_World1 = SpinP_switch + 1; 

  NPROCS_ID1 = (int*)malloc(sizeof(int)*numprocs); 
  Comm_World1 = (int*)malloc(sizeof(int)*numprocs); 
  NPROCS_WD1 = (int*)malloc(sizeof(int)*Num_Comm_World1); 
  Comm_World_StartID1 = (int*)malloc(sizeof(int)*Num_Comm_World1); 
  MPI_CommWD1 = (MPI_Comm*)malloc(sizeof(MPI_Comm)*Num_Comm_World1);

  Make_Comm_Worlds(mpi_comm_level1, myid, numprocs, Num_Comm_World1, &myworld1, MPI_CommWD1, 
                   NPROCS_ID1, Comm_World1, NPROCS_WD1, Comm_World_StartID1);

  /****************************************************
           calculation of dimension of matrix
  ****************************************************/

  Nmat = 0;
  for (i=1; i<=atomnum; i++){
    wanA  = WhatSpecies[i];
    Nmat += Spe_Total_CNO[wanA];
  }

  /****************************************************
                  total core charge
  ****************************************************/

  TZ = 0.0;
  for (i=1; i<=atomnum; i++){
    wan = WhatSpecies[i];
    TZ = TZ + Spe_Core_Charge[wan];
  }

  /****************************************************
                   calculate TNZE
  ****************************************************/

  TNZE = 0;
  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    Cwan = WhatSpecies[Gc_AN];

    nj = 0;
    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
      Gh_AN = natn[Gc_AN][h_AN];
      Hwan = WhatSpecies[Gh_AN];
      nj += Spe_Total_NO[Hwan];
    }

    TNZE += Spe_Total_NO[Cwan]*nj;
  }

  /****************************************************
                   allocate arrays
  ****************************************************/

  conv_ind2 = (int**)malloc(sizeof(int*)*Nmat);
  for (i=0; i<Nmat; i++){
    conv_ind2[i] = (int*)malloc(sizeof(int)*Nmat);
  } 

  DMall = (double***)malloc(sizeof(double**)*3);
  for (ChemP_point=0; ChemP_point<3; ChemP_point++){
    DMall[ChemP_point] = (double**)malloc(sizeof(double*)*(SpinP_switch+1));
    for (spin=0; spin<(SpinP_switch+1); spin++){
      DMall[ChemP_point][spin] = (double*)malloc(sizeof(double)*TNZE);
      for (i=0; i<TNZE; i++) DMall[ChemP_point][spin][i] = 0.0;
    }
  }

  Hall = (double**)malloc(sizeof(double*)*(SpinP_switch+1));
  for (spin=0; spin<(SpinP_switch+1); spin++){
    Hall[spin] = (double*)malloc(sizeof(double)*TNZE);
  }

  Sall = (double*)malloc(sizeof(double)*TNZE);

  MP = (int*)malloc(sizeof(int)*(atomnum+1));
    
  /* set MP */

  Anum = 1;
  for (i=1; i<=atomnum; i++){
    MP[i] = Anum;
    wanA = WhatSpecies[i];
    Anum += Spe_Total_CNO[wanA];
  }

  conv_index_row = (int*)malloc(sizeof(int)*TNZE);
  conv_index_col = (int*)malloc(sizeof(int)*TNZE);

  /******************************************
      MPI communication of nh and CntOLP
  ******************************************/

  /* allocation of arrays */

  My_NZeros = (int*)malloc(sizeof(int)*numprocs);
  is1 = (int*)malloc(sizeof(int)*numprocs);
  ie1 = (int*)malloc(sizeof(int)*numprocs);

  /* find my total number of non-zero elements in myid */

  My_NZeros[myid] = 0;
  for (MA_AN=1; MA_AN<=Matomnum; MA_AN++){
    GA_AN = M2G[MA_AN];
    wanA = WhatSpecies[GA_AN];
    tnoA = Spe_Total_CNO[wanA];

    num = 0;      
    for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
      GB_AN = natn[GA_AN][LB_AN];
      wanB = WhatSpecies[GB_AN];
      tnoB = Spe_Total_CNO[wanB];
      num += tnoB;
    }

    My_NZeros[myid] += tnoA*num;
  }

  for (ID=0; ID<numprocs; ID++){
    MPI_Bcast(&My_NZeros[ID],1,MPI_INT,ID,mpi_comm_level1);
  }

  tnum = 0;
  for (ID=0; ID<numprocs; ID++){
    tnum += My_NZeros[ID];
  }  

  is1[0] = 0;
  ie1[0] = My_NZeros[0] - 1;

  for (ID=1; ID<numprocs; ID++){
    is1[ID] = ie1[ID-1] + 1;
    ie1[ID] = is1[ID] + My_NZeros[ID] - 1;
  }  

  /******************
       for nh
  ******************/

  for (spin=0; spin<=SpinP_switch; spin++){

    /* set nh */

    k = is1[myid];
    for (MA_AN=1; MA_AN<=Matomnum; MA_AN++){
      GA_AN = M2G[MA_AN];
      wanA = WhatSpecies[GA_AN];
      tnoA = Spe_Total_CNO[wanA];
      for (i=0; i<tnoA; i++){
	for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
	  GB_AN = natn[GA_AN][LB_AN];
	  wanB = WhatSpecies[GB_AN];
	  tnoB = Spe_Total_CNO[wanB];
	  for (j=0; j<tnoB; j++){
	    Hall[spin][k] = nh[spin][MA_AN][LB_AN][i][j]; 
	    k++;
	  }
	}
      }
    }

    /* MPI Hall */

    MPI_Barrier(mpi_comm_level1);
    
    for (ID=0; ID<numprocs; ID++){
      k = is1[ID];
      MPI_Bcast(&Hall[spin][k], My_NZeros[ID], MPI_DOUBLE, ID, mpi_comm_level1);
    }

  }

  /******************
      for CntOLP
  ******************/

  k = is1[myid];
  for (MA_AN=1; MA_AN<=Matomnum; MA_AN++){
    GA_AN = M2G[MA_AN];    
    wanA = WhatSpecies[GA_AN];
    tnoA = Spe_Total_CNO[wanA];
    for (i=0; i<tnoA; i++){
      for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
        GB_AN = natn[GA_AN][LB_AN];
        wanB = WhatSpecies[GB_AN];
        tnoB = Spe_Total_CNO[wanB];
        for (j=0; j<tnoB; j++){
          Sall[k] = CntOLP[MA_AN][LB_AN][i][j]; 
          k++;
	}
      }
    }
  }

  /* MPI S1 */

  MPI_Barrier(mpi_comm_level1);

  for (ID=0; ID<numprocs; ID++){
    k = is1[ID];
    MPI_Bcast(&Sall[k], My_NZeros[ID], MPI_DOUBLE, ID, mpi_comm_level1);
  }

  /******************************************
   set up conv_index_row and conv_index_col
  ******************************************/

  k = is1[myid];

  for (MA_AN=1; MA_AN<=Matomnum; MA_AN++){
    GA_AN = M2G[MA_AN];    
    wanA = WhatSpecies[GA_AN];
    tnoA = Spe_Total_CNO[wanA];
    for (i=0; i<tnoA; i++){
      for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
        GB_AN = natn[GA_AN][LB_AN];
        wanB = WhatSpecies[GB_AN];
        tnoB = Spe_Total_CNO[wanB];
        for (j=0; j<tnoB; j++){

          conv_index_row[k] = MP[GA_AN] + i - 1;
          conv_index_col[k] = MP[GB_AN] + j - 1;
          k++;
	}
      }
    }
  }

  for (ID=0; ID<numprocs; ID++){
    k = is1[ID];
    MPI_Bcast(&conv_index_row[k], My_NZeros[ID], MPI_INT, ID, mpi_comm_level1);
  }

  for (ID=0; ID<numprocs; ID++){
    k = is1[ID];
    MPI_Bcast(&conv_index_col[k], My_NZeros[ID], MPI_INT, ID, mpi_comm_level1);
  }

  /****************************************************
                 set up NZE and PZE
  ****************************************************/

  NZE = (int*)malloc(sizeof(int)*Nmat);
  PZE = (int*)malloc(sizeof(int)*Nmat);

  for (i=0; i<Nmat; i++){
    NZE[i] = 0;
    PZE[i] = -1;
  }

  for (k=0; k<TNZE; k++){

    i = conv_index_row[k];
    
    NZE[i]++;
    if (PZE[i]==-1) PZE[i] = k;
  }

  /* freeing of arrays */

  free(My_NZeros);
  free(ie1);

  if (0<measure_time){
    dtime(&etime);
    printf("myid=%2d time1=%15.12f\n",myid,etime-stime);fflush(stdout);
  }

  /*************************************************************************
   *************************************************************************

     Part 1: contour integration with a chemical potential \mu_0

     main loop which is decomposed to 
     spin and energy point

   *************************************************************************
  *************************************************************************/

  MPI_Barrier(mpi_comm_level1);

  if (0<measure_time) dtime(&stime);

  if (SCF_iter<4){
    ChemP_trial[0] = ChemP - 0.05;
    ChemP_trial[1] = ChemP;
    ChemP_trial[2] = ChemP + 0.05;
  }
  else if (SCF_iter<10){
    ChemP_trial[0] = ChemP - 0.002;
    ChemP_trial[1] = ChemP;
    ChemP_trial[2] = ChemP + 0.002;
  }
  else{
    ChemP_trial[0] = ChemP - 0.0002;
    ChemP_trial[1] = ChemP;
    ChemP_trial[2] = ChemP + 0.0002;
  }



  Nloop1 = (SpinP_switch+1)*ON2_Npoles*3;

  if (Nloop1<=numprocs){

    Num_Comm_World2B = Nloop1;
                  
    NPROCS_ID2B = (int*)malloc(sizeof(int)*numprocs);
    Comm_World2B = (int*)malloc(sizeof(int)*numprocs);
    NPROCS_WD2B = (int*)malloc(sizeof(int)*Num_Comm_World2B);
    Comm_World_StartID2B = (int*)malloc(sizeof(int)*Num_Comm_World2B);
    MPI_CommWD2B = (MPI_Comm*)malloc(sizeof(MPI_Comm)*Num_Comm_World2B);

    Make_Comm_Worlds(mpi_comm_level1, myid, numprocs, Num_Comm_World2B, &myworld2B, MPI_CommWD2B, 
		     NPROCS_ID2B, Comm_World2B, NPROCS_WD2B, Comm_World_StartID2B);

    loop1 = myworld2B; 

    /* get ChemP_point, spin, and Epoint */

    ChemP_point = loop1/(ON2_Npoles*(SpinP_switch+1));
    spin = (loop1 - ChemP_point*ON2_Npoles*(SpinP_switch+1))/ON2_Npoles;
    Epoint = loop1 - ChemP_point*ON2_Npoles*(SpinP_switch+1) - spin*ON2_Npoles;

    OND_Solver(myid,numprocs,
	       Nloop1,
	       MPI_CommWD2B[myworld2B],
               ChemP_point,
	       spin, 
               Epoint,
               ChemP_trial, 
	       Nmat,NZE,PZE,
	       Hall,Sall,DMall,
	       conv_index_row, 
	       conv_index_col,
	       conv_ind2,
	       TNZE);
  }  

  else {

    for ( loop1=myid; loop1<Nloop1; loop1+=numprocs ){

      /* get ChemP_point, spin, and Epoint */

      ChemP_point = loop1/(ON2_Npoles*(SpinP_switch+1));
      spin = (loop1 - ChemP_point*ON2_Npoles*(SpinP_switch+1))/ON2_Npoles;
      Epoint = loop1 - ChemP_point*ON2_Npoles*(SpinP_switch+1) - spin*ON2_Npoles;

      printf("loop1=%2d ChemP_point=%2d spin=%2d Epoint=%2d\n",
              loop1,ChemP_point,spin,Epoint);fflush(stdout);


      OND_Solver( myid,
                  numprocs,
		  Nloop1,
		  mpi_comm_level1,
                  ChemP_point,
		  spin,
                  Epoint,
                  ChemP_trial, 
		  Nmat,NZE,PZE,
		  Hall,Sall,DMall,
		  conv_index_row, 
		  conv_index_col,
		  conv_ind2,
		  TNZE);

      /*
      for (i=0; i<TNZE; i++){
        printf("i=%2d DMall=%15.12f %15.12f\n",i,creal(DMall[0][i]),cimag(DMall[0][i]));
      }     

      MPI_Finalize();
      exit(0);
      */

    }
  }

  if (0<measure_time){
    dtime(&etime);
    printf("myid=%2d time2=%15.12f\n",myid,etime-stime);fflush(stdout);
  }

  /******************************************
    calculate "each" number of electrons
    on each processor, and MPI communicate 
    to get the total sum.
  ******************************************/

  MPI_Barrier(mpi_comm_level1);

  if (0<measure_time) dtime(&stime);

  for (ChemP_point=0; ChemP_point<3; ChemP_point++){
    eachsum[ChemP_point] = 0.0;
    for (spin=0; spin<=SpinP_switch; spin++){
      for (i=0; i<TNZE; i++){
	eachsum[ChemP_point] += DMall[ChemP_point][spin][i]*Sall[i];
      }
    }

    MPI_Allreduce( &eachsum[ChemP_point], &totalsum[ChemP_point],
                   1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1 );
  }

  if (SpinP_switch==0){
    for (ChemP_point=0; ChemP_point<3; ChemP_point++){
      eachsum[ChemP_point]  *= 2.0;
      totalsum[ChemP_point] *= 2.0;
    }
  }

  if (0<measure_time){
    dtime(&etime);
    printf("myid=%2d time3=%15.12f\n",myid,etime-stime);fflush(stdout);
  }

  for (ChemP_point=0; ChemP_point<3; ChemP_point++){
    printf("myid=%2d eachsum=%15.12f totalsum=%15.12f\n",
           myid,eachsum[ChemP_point],totalsum[ChemP_point]);
  }

  

  /******************************************
  ******************************************/

  x0 = ChemP_trial[0];  
  x1 = ChemP_trial[1];  
  x2 = ChemP_trial[2];  

  y0 = totalsum[0];
  y1 = totalsum[1];
  y2 = totalsum[2];

  printf("x0=%15.12f\n",x0);
  printf("x1=%15.12f\n",x1);
  printf("x2=%15.12f\n",x2);
 
  printf("y0=%15.12f\n",y0);
  printf("y1=%15.12f\n",y1);
  printf("y2=%15.12f\n",y2);
   
  d = x0*x0*x1 + x1*x1*x2 + x2*x2*x0
     -x1*x2*x2 - x0*x0*x2 - x0*x1*x1;

  a = ((x1-x2)*y0 + (x2-x0)*y1 + (x0-x1)*y2)/d;
  b = ((x2*x2-x1*x1)*y0 + (x0*x0-x2*x2)*y1 + (x1*x1-x0*x0)*y2)/d;
  c = (x1*x2*(x1-x2)*y0 + x2*x0*(x2-x0)*y1 + x0*x1*(x0-x1)*y2)/d;

  xopt = (-b+sqrt(b*b-4*a*(c-(TZ-system_charge))))/(2.0*a); 

  printf("a=%15.12f\n",a);
  printf("b=%15.12f\n",b);
  printf("c=%15.12f\n",c);
  printf("xopt=%15.12f\n",xopt);

  d = x1*y0 - x2*y0 - x0*y1 + x2*y1 + x0*y2 - x1*y2;

  weight_DM[0] = ((x2*y1-x1*y2) + (x1-x2)*(TZ-system_charge) + (y2-y1)*xopt)/d;
  weight_DM[1] = ((x0*y2-x2*y0) + (x2-x0)*(TZ-system_charge) + (y0-y2)*xopt)/d;
  weight_DM[2] = ((x1*y0-x0*y1) + (x0-x1)*(TZ-system_charge) + (y1-y0)*xopt)/d;

  printf("p0=%15.12f\n",weight_DM[0]);
  printf("p1=%15.12f\n",weight_DM[1]);
  printf("p2=%15.12f\n",weight_DM[2]);

  ChemP = xopt;

  /*************************************************************************
           store the density matrix by the part 1 into an array DM
  **************************************************************************/

  MPI_Barrier(mpi_comm_level1);

  for (spin=0; spin<=SpinP_switch; spin++){
    for (MA_AN=1; MA_AN<=Matomnum; MA_AN++){
      GA_AN = M2G[MA_AN];
      wanA = WhatSpecies[GA_AN];
      tnoA = Spe_Total_CNO[wanA];
      for (i=0; i<tnoA; i++){
	for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
	  GB_AN = natn[GA_AN][LB_AN];
	  wanB = WhatSpecies[GB_AN];
	  tnoB = Spe_Total_CNO[wanB];
	  for (j=0; j<tnoB; j++){
	    CDM[spin][MA_AN][LB_AN][i][j] = 0.0;
	  }
	}
      }
    }
  }

  for (ChemP_point=0; ChemP_point<3; ChemP_point++){
    for (spin=0; spin<=SpinP_switch; spin++){

      MPI_Allreduce( &DMall[ChemP_point][spin][0], &Sall[0], TNZE, 
                     MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

      k = is1[myid];
      for (MA_AN=1; MA_AN<=Matomnum; MA_AN++){
	GA_AN = M2G[MA_AN];
	wanA = WhatSpecies[GA_AN];
	tnoA = Spe_Total_CNO[wanA];
	for (i=0; i<tnoA; i++){
	  for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
	    GB_AN = natn[GA_AN][LB_AN];
	    wanB = WhatSpecies[GB_AN];
	    tnoB = Spe_Total_CNO[wanB];
	    for (j=0; j<tnoB; j++){
	      CDM[spin][MA_AN][LB_AN][i][j] += weight_DM[ChemP_point]*Sall[k];
	      k++;
	    }
	  }
	}
      }
    }
  }

  /*************************************************************************
   *************************************************************************

     Part 2: shifted invert iteration for correction of 
             the chemical potential

   *************************************************************************
  **************************************************************************/

  /*

  if (0<measure_time) dtime(&stime);

  NEV = abs(TZ-totalsum-system_charge) + 50;
  if (Nmat<NEV) NEV = Nmat;

  eval = (double**)malloc(sizeof(double*)*(SpinP_switch+1));
  for (spin=0; spin<(SpinP_switch+1); spin++){
    eval[spin] = (double*)malloc(sizeof(double)*NEV);
    for (i=0; i<NEV; i++) eval[spin][i] = 0.0;
  }
  
  evec = (double**)malloc(sizeof(double*)*(SpinP_switch+1));
  for (spin=0; spin<(SpinP_switch+1); spin++){
    evec[spin] = (double*)malloc(sizeof(double)*NEV*Nmat);
    for (i=0; i<NEV*Nmat; i++) evec[spin][i] = 0.0;
  }

  printf("myid=%2d NEV=%2d\n",myid,NEV);

  */

  /*********************************************************
              call Get_EigenStates_Near_pChemP 
  *********************************************************/

  /* store the chemical potential used in the part 1 */

  /*

  pChemP = ChemP;

  Get_EigenStates_Near_pChemP( myid, numprocs,
                               myworld1,
                               Num_Comm_World1, NPROCS_ID1, Comm_World1, 
                               NPROCS_WD1, Comm_World_StartID1, MPI_CommWD1, 
                               NEV, Nmat, pChemP, NZE, PZE, Hall, Sall, DMall, 
                               conv_index_row, conv_index_col, conv_ind2, TNZE, 
                               eval, evec );

  if (0<measure_time){
    dtime(&etime);
    printf("myid=%2d time4=%15.12f\n",myid,etime-stime);fflush(stdout);
  }

  */

  /*********************************************************
            find the correct chemical potential
  *********************************************************/

  /*

  if      (SpinP_switch==0) spin_degeneracy = 2.0;
  else if (SpinP_switch==1) spin_degeneracy = 1.0;

  ChemP_MIN = pChemP - 10.0;
  ChemP_MAX = pChemP + 10.0;

  po = 0;
  loopN = 0;

  do {

    ChemP = 0.5*(ChemP_MIN+ChemP_MAX);

    totalsum1 = 0.0;
    for (spin=0; spin<=SpinP_switch; spin++){
      for (i=0; i<NEV; i++){
	f0 = 1.0/(1.0+exp((eval[spin][i] - pChemP)*Beta)); 
	f1 = 1.0/(1.0+exp((eval[spin][i] -  ChemP)*Beta)); 
	totalsum1 += spin_degeneracy*(f1 - f0);
      }
    }

    Dnum = (TZ - totalsum - totalsum1) - system_charge;

    if (0.0<=Dnum) ChemP_MIN = ChemP;
    else           ChemP_MAX = ChemP;
    if (fabs(Dnum)<1.0e-14) po = 1;

    loopN++;

  }
  while (po==0 && loopN<1000); 

  printf("pChemP=%18.15f\n",pChemP);
  printf(" ChemP=%18.15f\n",ChemP);

  */

  /*********************************************************
            add the correction by the part 2 to DM
  *********************************************************/

  /*

  for (spin=0; spin<=SpinP_switch; spin++){

    for (k=0; k<NEV; k++){

      f0 = 1.0/(1.0+exp((eval[spin][k] - pChemP)*Beta)); 
      f1 = 1.0/(1.0+exp((eval[spin][k] -  ChemP)*Beta)); 

      for (MA_AN=1; MA_AN<=Matomnum; MA_AN++){

	GA_AN = M2G[MA_AN];
	wanA = WhatSpecies[GA_AN];
	tnoA = Spe_Total_CNO[wanA];

	for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){

	  GB_AN = natn[GA_AN][LB_AN];
	  wanB = WhatSpecies[GB_AN];
	  tnoB = Spe_Total_CNO[wanB];

	  for (i=0; i<tnoA; i++){

            i1 = MP[GA_AN] + i - 1;

	    for (j=0; j<tnoB; j++){
 
	      j1 = MP[GB_AN] + j - 1;
	       
  	      CDM[spin][MA_AN][LB_AN][i][j] += (f1 - f0)*evec[spin][k*Nmat+i1]*evec[spin][k*Nmat+j1];
	    }
	  }
	}
      }
    }
  }

  */

  /****************************************************
               calculation of band energy
  ****************************************************/

  My_Eele1[0] = 0.0;
  My_Eele1[1] = 0.0;

  for (spin=0; spin<=SpinP_switch; spin++){
    for (MA_AN=1; MA_AN<=Matomnum; MA_AN++){
      GA_AN = M2G[MA_AN];
      wanA = WhatSpecies[GA_AN];
      tnoA = Spe_Total_CNO[wanA];
      for (j=0; j<=FNAN[GA_AN]; j++){
	wanB = WhatSpecies[natn[GA_AN][j]];
	tnoB = Spe_Total_CNO[wanB];
	for (k=0; k<tnoA; k++){
	  for (l=0; l<tnoB; l++){
	    My_Eele1[spin] += CDM[spin][MA_AN][j][k][l]*nh[spin][MA_AN][j][k][l];
	  }
	}
      }
    }
  }

  MPI_Barrier(mpi_comm_level1);

  /* MPI, My_Eele1 */

  for (spin=0; spin<=SpinP_switch; spin++){
    MPI_Allreduce(&My_Eele1[spin], &Eele1[spin], 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
  }

  if (SpinP_switch==0) Eele1[1] = Eele1[0];

  /****************************************************
                   freeing of arrays
  ****************************************************/

  /*
  for (spin=0; spin<(SpinP_switch+1); spin++){
    free(evec[spin]);
  }
  free(evec);

  for (spin=0; spin<(SpinP_switch+1); spin++){
    free(eval[spin]);
  }
  free(eval);
  */

  free(is1);

  free(NZE);
  free(PZE);

  for (i=0; i<Nmat; i++){
    free(conv_ind2[i]);
  } 
  free(conv_ind2);

  for (ChemP_point=0; ChemP_point<3; ChemP_point++){
    for (spin=0; spin<(SpinP_switch+1); spin++){
      free(DMall[ChemP_point][spin]);
    }
    free(DMall[ChemP_point]);
  }
  free(DMall);

  for (spin=0; spin<(SpinP_switch+1); spin++){
    free(Hall[spin]);
  }
  free(Hall);

  free(Sall);

  free(MP);
  free(conv_index_row);
  free(conv_index_col);

  /* freeing of arrays for the first world */

  if (Num_Comm_World1<=numprocs){
    MPI_Comm_free(&MPI_CommWD1[myworld1]);
  }

  free(NPROCS_ID1);
  free(Comm_World1);
  free(NPROCS_WD1);
  free(Comm_World_StartID1);
  free(MPI_CommWD1);

  /* freeing of arrays for the second B world */

  if (Nloop1<=numprocs){

    MPI_Comm_free(&MPI_CommWD2B[myworld2B]);
    free(NPROCS_ID2B);
    free(Comm_World2B);
    free(NPROCS_WD2B);
    free(Comm_World_StartID2B);
    free(MPI_CommWD2B);
  }

  /* for elapsed time */
  MPI_Barrier(mpi_comm_level1);
  dtime(&TEtime);
  return (TEtime-TStime);
}








static void Get_EigenStates_Near_pChemP(
  /* input */
  int myid0,
  int numprocs0, 
  int myworld1,
  int Num_Comm_World1, 
  int *NPROCS_ID1,
  int *NPROCS_WD1,
  int *Comm_World1, 
  int *Comm_World_StartID1,
  MPI_Comm *MPI_CommWD1,
  int NEV, 
  int Nmat, 
  double pChemP, 
  int *NZE,
  int *PZE,
  double **Hall,
  double *Sall,
  double **DMall,
  int *conv_index_row,
  int *conv_index_col,
  int **conv_ind2,
  int TNZE, 
  /* output */
  double **eval, 
  double **evec)
{
  UF_long i,j,k,l,l1,spin;
  UF_long order,iev;
  cs_dl *A,*T,*Scomp,*Hcomp;   
  cs_dls *s;
  cs_dln *N;
  UF_long ok;
  int po,nloop,NEV0;
  int TNZE2,col,row;
  double tol,sum,coe;
  double *x,*b,*c,**evec_tmp;
  double *Ssub,*Hsub;
  double new_sumev,old_sumev;
  double threshold=1.0e-13; 
  INTEGER itype;
  INTEGER n,lda,ldb,lwork,info; 
  double *w,*work,*lambda;
  double *order_w;
  char jobz = 'V';
  char uplo ='U';
  int *is1,*ie1;
  int numprocs1,myid1,ID,num;
  double av_num;

  /* allocation of arrays */

  T = cs_dl_spalloc (Nmat, Nmat, TNZE, 1, 1) ; 
  x = cs_dl_malloc (Nmat, sizeof(double)) ;  /* get workspace */
  b = cs_dl_malloc (Nmat, sizeof(double)) ;  /* get workspace */
  c = cs_dl_malloc (Nmat, sizeof(double)) ;  /* get workspace */

  Ssub = (double*)malloc(sizeof(double)*NEV*NEV);
  Hsub = (double*)malloc(sizeof(double)*NEV*NEV);

  evec_tmp = (double**)malloc(sizeof(double*)*NEV);
  for (i=0; i<NEV; i++){
    evec_tmp[i] = (double*)malloc(sizeof(double)*Nmat);
  }

  w = (double*)malloc(sizeof(double)*NEV);
  lambda = (double*)malloc(sizeof(double)*(NEV+2));
  work = (double*)malloc(sizeof(double)*3*NEV);
  order_w = (double*)malloc(sizeof(double)*(NEV+2));

  itype = 1;
  n = NEV; 
  lda = NEV;
  ldb = NEV;
  lwork = 3*NEV;   

  /* find the numbers of partions for MPI */

  MPI_Comm_size(MPI_CommWD1[myworld1],&numprocs1);
  MPI_Comm_rank(MPI_CommWD1[myworld1],&myid1);
  
  is1 = (int*)malloc(sizeof(int)*numprocs1);
  ie1 = (int*)malloc(sizeof(int)*numprocs1);

  if ( numprocs1<=NEV ){

    av_num = (double)NEV/(double)numprocs1;

    for (ID=0; ID<numprocs1; ID++){
      is1[ID] = (int)(av_num*(double)ID); 
      ie1[ID] = (int)(av_num*(double)(ID+1))-1; 
    }

    is1[0] = 0;
    ie1[numprocs1-1] = NEV - 1;
  }

  else {
    for (ID=0; ID<NEV; ID++){
      is1[ID] = ID; 
      ie1[ID] = ID;
    }
    for (ID=NEV; ID<numprocs1; ID++){
      is1[ID] =  0;
      ie1[ID] = -2;
    }
  }

  /* set a matrix Scomp */

  for (i=0; i<Nmat; i++){
    for (j=0; j<Nmat; j++){
      conv_ind2[i][j] = -1;
    }  
  }  

  TNZE2 = 0;

  for (i=0; i<TNZE; i++){

    row = conv_index_row[i];                /* row index    */
    col = conv_index_col[i];                /* column index */

    if (conv_ind2[row][col]==-1){  

      conv_ind2[row][col] = TNZE2;
      k = TNZE2;
      T->x[k] = 0.0 + 0.0*I;
      TNZE2++;
    }
    else {
      k = conv_ind2[row][col];
    }

    T->p[k]  = col;             /* column index */
    T->i[k]  = row;             /* row index */
    T->x[k] += Sall[i];         /* value */
  }

  T->nzmax = TNZE2;
  T->nz = TNZE2;

  Scomp = cs_dl_compress(T) ;      /* allocation */

  /* loop for spin */

  spin = myworld1;  /* spin=myworld1 */

 spinloop:

  /******************************************************************
                               LU decomposition
  *******************************************************************/

  /* set a matrix T */

  for (i=0; i<Nmat; i++){
    for (j=0; j<Nmat; j++){
      conv_ind2[i][j] = -1;
    }  
  }  

  TNZE2 = 0;

  for (i=0; i<TNZE; i++){

    row = conv_index_row[i];                /* row index    */
    col = conv_index_col[i];                /* column index */

    if (conv_ind2[row][col]==-1){  

      conv_ind2[row][col] = TNZE2;
      k = TNZE2;
      T->x[k] = 0.0 + 0.0*I;
      TNZE2++;
    }
    else {
      k = conv_ind2[row][col];
    }

    T->p[k]  = col;                                /* column index */
    T->i[k]  = row;                                /* row index */
    T->x[k] += (pChemP*Sall[i] - Hall[spin][i]);   /* value */
  }

  T->nzmax = TNZE2;
  T->nz = TNZE2;

  /* A = compressed-column form of T */

  A = cs_dl_compress(T) ;      /* allocation */

  /* decompose A into L*U */

  tol = 1.0e-5;
  order = 1;

  /* ordering and symbolic analysis */

  s = cs_dl_sqr (order, A, 0) ; /* allocation */      

  /* numeric LU factorization */

  N = cs_dl_lu (A, s, tol) ;   /* allocation */                    

  /******************************************************************
              shift-invert method using the LU decomposition.
  *******************************************************************/

  /* set a matrix Hcomp */

  for (i=0; i<Nmat; i++){
    for (j=0; j<Nmat; j++){
      conv_ind2[i][j] = -1;
    }  
  }  

  TNZE2 = 0;

  for (i=0; i<TNZE; i++){

    row = conv_index_row[i];                /* row index    */
    col = conv_index_col[i];                /* column index */

    if (conv_ind2[row][col]==-1){  

      conv_ind2[row][col] = TNZE2;
      k = TNZE2;
      T->x[k] = 0.0 + 0.0*I;
      TNZE2++;
    }
    else {
      k = conv_ind2[row][col];
    }

    T->p[k]  = col;                 /* column index */
    T->i[k]  = row;                 /* row index */
    T->x[k] += Hall[spin][i];       /* value */
  }

  T->nzmax = TNZE2;
  T->nz = TNZE2;

  Hcomp = cs_dl_compress(T);      /* allocation */

  /* set up initial vectors */

  for (i=0; i<NEV; i++){
    for (k=0; k<Nmat; k++){
      evec[spin][i*Nmat+k] = 0.0;
    }
    evec[spin][i*Nmat+i] = 1.0;
  }

  /* iteration */

  po = 0;
  nloop = 0;
  new_sumev = 1.0e+10;

  do {

    /* (H-ep*S)^{-1}*S*c */

    for (iev=is1[myid1]; iev<=ie1[myid1]; iev++){

      /* set c */
      for (j=0; j<Nmat; j++){
	c[j] = evec[spin][iev*Nmat+j];
      }

      /* calculate S*c */

      for (k=0; k<Nmat; k++) b[k] = 0.0;
      cs_dl_gaxpy(Scomp, c, b);

      /* x = b(p) */
      cs_dl_ipvec (N->pinv, b, x, Nmat);

      /* x = L\x */
      cs_dl_lsolve (N->L, x);               
  
      /* x = U\x */
      cs_dl_usolve (N->U, x);               

      /* evec[spin][iev](q) = x */
      cs_dl_ipvec (s->q, x, &evec[spin][iev*Nmat], Nmat);

    } /* iev */

    /* MPI: evec */

    for (ID=0; ID<numprocs1; ID++){

      num = ie1[ID] - is1[ID] + 1;

      if (0<num){ 
        MPI_Bcast(&evec[spin][ is1[ID]*Nmat ], num*Nmat, MPI_DOUBLE, ID, MPI_CommWD1[myworld1]);
      }
    }

    /* solve a generalized eigenvalue problem using evec */    

    /* calculate the effective S */

    for (j=is1[myid1]; j<=ie1[myid1]; j++){

      for (k=0; k<Nmat; k++) b[k] = 0.0;
      cs_dl_gaxpy(Scomp, &evec[spin][j*Nmat], b);

      for (i=0; i<NEV; i++){

	sum = 0.0;
	for (k=0; k<Nmat; k++){
	  sum += evec[spin][i*Nmat+k]*b[k];
	}

	Ssub[NEV*j+i] = sum;
      }
    }

    /* MPI: Ssub */

    for (ID=0; ID<numprocs1; ID++){

      num = ie1[ID] - is1[ID] + 1;

      if (0<num){ 
        MPI_Bcast(&Ssub[ is1[ID]*NEV ], num*NEV, MPI_DOUBLE, ID, MPI_CommWD1[myworld1]);
      }
    }

    /* calculate the effective H */

    for (j=is1[myid1]; j<=ie1[myid1]; j++){

      for (k=0; k<Nmat; k++) b[k] = 0.0;
      cs_dl_gaxpy(Hcomp, &evec[spin][j*Nmat], b);

      for (i=0; i<NEV; i++){

	sum = 0.0;
	for (k=0; k<Nmat; k++){
	  sum += evec[spin][i*Nmat+k]*b[k];
	}

	Hsub[NEV*j+i] = sum;
      }
    }

    /* MPI: Hsub */

    for (ID=0; ID<numprocs1; ID++){

      num = ie1[ID] - is1[ID] + 1;

      if (0<num){ 
        MPI_Bcast(&Hsub[ is1[ID]*NEV ], num*NEV, MPI_DOUBLE, ID, MPI_CommWD1[myworld1]);
      }
    }

    /* diagonalize Hsub d = ep Ssub d */

    F77_NAME(dsygv,DSYGV)(&itype, &jobz, &uplo, &n, Hsub, &lda, Ssub, &ldb, w, work, &lwork, &info); 

    /* calculate new trial vectors */

    for (i=is1[myid1]; i<=ie1[myid1]; i++){

      for (k=0; k<Nmat; k++) c[k] = 0.0;

      for (j=0; j<NEV; j++){
         
	coe = Hsub[i*NEV+j];

	for (k=0; k<Nmat; k++){
	  c[k] += coe*evec[spin][j*Nmat+k];  
	}
      }

      for (k=0; k<Nmat; k++){
	evec_tmp[i][k] = c[k];
      }
    }

    /*
      printf("evec\n");
    */

    for (i=is1[myid1]; i<=ie1[myid1]; i++){
      for (k=0; k<Nmat; k++){
	evec[spin][i*Nmat+k] = evec_tmp[i][k];

	/*
          printf("i=%2d k=%2d evec=%15.12f\n",i,k,evec[spin][i*Nmat+k]);
	*/

      }
    }

    /***************************************************
       check convergence for the top NEV/2 of eigenvalues 
       being close to pChemP 
    ***************************************************/

    for (i=0; i<NEV; i++){
      lambda[i+1] = fabs(w[i]-pChemP);
      order_w[i+1] = (double)i;
    }      

    qsort_double(NEV,lambda,order_w);

    old_sumev = new_sumev;
    new_sumev = 0.0;
    NEV0 = NEV/2;

    for (i=0; i<NEV0; i++){
      j = (int)order_w[i+1]; 
      new_sumev += fabs(w[j]);
    }      


    /*
    printf("        pChemP=%18.15f\n",pChemP);fflush(stdout);
    for (i=0; i<NEV; i++){
      printf("i=%3d        w=%18.15f  order=%2d diff=%18.15f\n",
      i,w[i],(int)order_w[i+1],fabs(w[i]-pChemP) );fflush(stdout);
    }
    */


    /*
      old_sumev = new_sumev;
      new_sumev = 0.0;
      for (i=0; i<NEV; i++){
      new_sumev += fabs(w[i]);
      } 
    */
     

    if (fabs(new_sumev-old_sumev)<threshold){

      po = 1;

      /* MPI: evec */

      for (ID=0; ID<numprocs1; ID++){

	num = ie1[ID] - is1[ID] + 1;

	if (0<num){ 
	  MPI_Bcast(&evec[spin][ is1[ID]*Nmat ], num*Nmat, MPI_DOUBLE, ID, MPI_CommWD1[myworld1]);
	}
      }
    }


    printf("spin=%2d nloop=%3d diff     =%18.15f po=%2d\n",
            spin,nloop,fabs(new_sumev-old_sumev),po);fflush(stdout);


    /* increment nloop */

    nloop++;

  } while (po==0 && nloop<1000);

  if (po==0){
    printf("Error in Cluster_DFT_ON2: #2\n");
    exit(0);
  }

  printf("spin=%2d nloop=%3d threshold=%18.15f\n",spin,nloop,threshold);fflush(stdout);

  /* store eigenvalues */

  for (i=0; i<NEV; i++){
    eval[spin][i] = w[i];
  }

  /* free arrays */

  cs_dl_spfree(Hcomp);
  cs_dl_nfree (N) ;
  cs_dl_sfree (s) ;
  cs_dl_spfree(A);

  if (SpinP_switch==1 && numprocs0==1 && spin==0){
    spin++;  
    goto spinloop;
  }

  /***********************************************
              MPI communication: eval
  ***********************************************/

  if (SpinP_switch==1 && numprocs0!=1){

    MPI_Bcast(&eval[0][0], NEV, MPI_DOUBLE, Comm_World_StartID1[0], mpi_comm_level1);
    MPI_Bcast(&eval[1][0], NEV, MPI_DOUBLE, Comm_World_StartID1[1], mpi_comm_level1);

    MPI_Bcast(&evec[0][0], NEV*Nmat, MPI_DOUBLE, Comm_World_StartID1[0], mpi_comm_level1);
    MPI_Bcast(&evec[1][0], NEV*Nmat, MPI_DOUBLE, Comm_World_StartID1[1], mpi_comm_level1);
  }

  /* freeing of arrays */

  free(ie1);
  free(is1);

  cs_dl_spfree(Scomp);

  free(work);
  free(w);
  free(lambda);
  free(order_w);

  for (i=0; i<NEV; i++){
    free(evec_tmp[i]);
  }
  free(evec_tmp);

  free(Hsub);
  free(Ssub);

  cs_dl_free(c);
  cs_dl_free(b);
  cs_dl_free(x);
  cs_dl_spfree(T);
}










static void OND_Solver(
  int myid,
  int numprocs,
  int Nloop1,
  MPI_Comm MPI_Comm,
  int ChemP_point,
  int spin,
  int Epoint,
  double ChemP_trial[3],
  int Nmat,
  int *NZE,
  int *PZE,
  double **Hall,
  double *Sall,
  double ***DMall,
  int *conv_index_row,
  int *conv_index_col,
  int **conv_ind2,
  int TNZE)
{
  int po,TNZE2,TNZE3,col,row;
  UF_long i,j,k,l,l1;
  UF_long kg,jg,j1,k1,order;
  cs_cl *A,*T; 
  cs_cls *s;
  cs_cln *N;
  double tol;
  double stime,etime;
  double stime1,etime1;
  double complex alpha,weight;
  double complex ctmp,ctmp2,ctmp3;
  int numprocs2B,myid2B,ID;
  double av_num;
  cs_complex_t *x,*b;
  idxtype *xadj,*adjncy;
  int io,jo,jo_min,n;
  int m0,m1,m2,m3,m4,m5;
  int n0,n1,n2,mp0,nsr;
  int nb0,nb1,nc,k2,kk;
  int nl,mpmax,mp,mm,numa;
  UF_long ok;
  int **Index_BC;
  int **SindB,**EindB;
  int **SindC,**EindC;
  int *invp;
  int ***Bnum[2],***Cnum,**Anum;
  int ****Bi[2],****Ci,***Ai;
  double complex ****Bx[2];
  double complex ****Cx;
  double complex ****IBx[2];
  double complex ****ICx;
  double complex ***Ax;
  double complex ***IAx;
  int m,i1,i2,j2,i3,j3,nsa,nsc,nsc2,nsb;
  int max_nsc,max_nsb;
  dcomplex **IA,***IS;
  dcomplex ***Lvec[2];
  dcomplex dcsum0,dcsum1,dcsum2,dcsum3,dcsum4;
  double complex **Vvec,**Vvec2;
  double complex **Q0,**Q1;
  double complex csum0,csum1,csum2,csum3,csum4,ctmp1;
  double time1,time2,time3,time4,time5,time6,time41,time42;
  double time31,time32,time33,time51,time52,time53,time54,time55;
  long int numop;
  dcomplex al,be;

  al.r = 1.0;
  al.i = 0.0;
  be.r = 0.0;
  be.i = 0.0;

  numop = 0;

  time1 = 0.0;
  time2 = 0.0;
  time3 = 0.0;
  time31= 0.0;
  time32= 0.0;
  time33= 0.0;
  time4 = 0.0;
  time41= 0.0;
  time42= 0.0;
  time5 = 0.0;
  time51= 0.0;
  time52= 0.0;
  time53= 0.0;
  time54= 0.0;
  time55= 0.0;
  time6 = 0.0;

  dtime(&stime);

  /* set alpha and weight */
  
  if (ON2_method[Epoint]==1){ /* poles */
    alpha  = ChemP_trial[ChemP_point] + I*(ON2_zp[Epoint].i/Beta);
    weight = -2.0*ON2_Rp[Epoint].r/Beta;
  }
  else if (ON2_method[Epoint]==2){ /* zeroth moment */
    alpha  = ON2_zp[Epoint].r + I*ON2_zp[Epoint].i;
    weight = I*ON2_Rp[Epoint].i;
  }

  /* set a matrix T */

  T = cs_cl_spalloc (Nmat, Nmat, TNZE, 1, 1) ; 

  if (1<measure_time){
    printf("Epoint =%2d %10.5f %10.5f\n",Epoint,ChemP_trial[ChemP_point],cimag(alpha) );
  }

  for (i=0; i<Nmat; i++){
    for (j=0; j<Nmat; j++){
      conv_ind2[i][j] = -1;
    }  
  }  

  TNZE2 = 0;

  for (i=0; i<TNZE; i++){

    row = conv_index_row[i];                /* row index    */
    col = conv_index_col[i];                /* column index */

    if (conv_ind2[row][col]==-1){  

      conv_ind2[row][col] = TNZE2;
      k = TNZE2;
      T->x[k] = 0.0;
      TNZE2++;
    }
    else {
      k = conv_ind2[row][col];
    }

    T->p[k]  = col;                                /* column index */
    T->i[k]  = row;                                /* row index    */
    T->x[k] += (alpha*Sall[i] - Hall[spin][i]);    /* value        */
  }

  T->nzmax = TNZE2;
  T->nz = TNZE2;

  /* A = compressed-column form of T */

  A = cs_cl_compress(T) ;
  cs_cl_spfree(T);

  /* sorting of A->i and A->x */

  for (i=0; i<Nmat; i++){
    n = A->p[i+1] - A->p[i];
    qsort_complex( (long)n, &(A->i[A->p[i]]), &(A->x[A->p[i]]) );
  }  

  /**************************************************
                   nested dissection 
  **************************************************/

  /* allocation of arrays */

  xadj = (idxtype*)malloc(sizeof(idxtype)*(Nmat+3));
  adjncy = (idxtype*)malloc(sizeof(idxtype)*(TNZE2+3));

  /* eliminate diagonal terms */

  TNZE3 = 0;
  for (i=0; i<Nmat; i++){
  
    xadj[i] = (idxtype)TNZE3;

    for (j=A->p[i]; j<A->p[i+1]; j++){
      if (A->i[j]!=i){
        adjncy[TNZE3] = (idxtype)A->i[j];         
        TNZE3++;
      }
    }
  }
  xadj[Nmat] = TNZE3;

  /*
  for (i=0; i<=Nmat; i++){
    printf("A i=%2d A->p=%2d\n",i,A->p[i]);fflush(stdout);
  }
   
  for (i=0; i<TNZE2; i++){
    printf("A i=%2d A->i=%2d\n",i,A->i[i]);fflush(stdout);
  }
  
  for (i=0; i<=Nmat; i++){
    printf("A i=%2d xadj=%2d\n",i,xadj[i]);fflush(stdout);
  }
  
  for (i=0; i<TNZE3; i++){
    printf("A i=%2d adjncy=%2d\n",i,adjncy[i]);fflush(stdout);
  }
  */

  /* get nl and mpmax */

  Partition_System( 1, &nl, &mpmax, Nmat, TNZE2, xadj, adjncy,
                    Index_BC, SindB, EindB, SindC, EindC );
  
  Index_BC = (int**)malloc(sizeof(int*)*(nl+1));
  for (n=0; n<(nl+1); n++){
    Index_BC[n] = (int*)malloc(sizeof(int)*Nmat);
    for (i=0; i<Nmat; i++){
      Index_BC[n][i] = i;
    }
  }

  SindB = (int**)malloc(sizeof(int*)*(nl+1));
  mp = mpmax;
  for (n=0; n<(nl+1); n++){
    SindB[n] = (int*)malloc(sizeof(int)*mp);
    for (m=0; m<mp; m++) SindB[n][m] = 0;
    mp /= 2;
  }  

  EindB = (int**)malloc(sizeof(int*)*(nl+1));
  mp = mpmax;
  for (n=0; n<(nl+1); n++){
    EindB[n] = (int*)malloc(sizeof(int)*mp);
    for (m=0; m<mp; m++) EindB[n][m] = 0;
    mp /= 2;
  }  

  SindC = (int**)malloc(sizeof(int*)*nl);
  mp = mpmax/2;
  for (n=0; n<nl; n++){
    SindC[n] = (int*)malloc(sizeof(int)*mp);
    for (m=0; m<mp; m++) SindC[n][m] = 0;
    mp /= 2;
  }  

  EindC = (int**)malloc(sizeof(int*)*nl);
  mp = mpmax/2;
  for (n=0; n<nl; n++){
    EindC[n] = (int*)malloc(sizeof(int)*mp);
    for (m=0; m<mp; m++) EindC[n][m] = 0;
    mp /= 2;
  }  

  /* set up Index_BC, SindB, EindB, SindC, and EindC */
 
  Partition_System( 0, &nl, &mpmax, Nmat, TNZE2, xadj, adjncy, 
                    Index_BC, SindB, EindB, SindC, EindC );

  /* make invp */

  invp = (int*)malloc(sizeof(int)*Nmat);

  for (i=0; i<Nmat; i++){
    j = Index_BC[nl][i];
    invp[j] = i;
  }

  /***********************************************************
             set up block matrices, B0, B1, and C
  ***********************************************************/

  Bnum[0] = (int***)malloc(sizeof(int**)*nl);
  Bnum[1] = (int***)malloc(sizeof(int**)*nl);
  Cnum  = (int***)malloc(sizeof(int**)*nl);

  Bi[0] = (int****)malloc(sizeof(int***)*nl);
  Bi[1] = (int****)malloc(sizeof(int***)*nl);
  Ci  = (int****)malloc(sizeof(int***)*nl);

  Bx[0] = (double complex****)malloc(sizeof(double complex***)*nl);
  Bx[1] = (double complex****)malloc(sizeof(double complex***)*nl);
  Cx  = (double complex****)malloc(sizeof(double complex***)*nl);

  IBx[0] = (double complex****)malloc(sizeof(double complex***)*nl);
  IBx[1] = (double complex****)malloc(sizeof(double complex***)*nl);
  ICx  = (double complex****)malloc(sizeof(double complex***)*nl);

  mp = mpmax/2;

  for (n=0; n<nl; n++){

    Bnum[0][n] = (int**)malloc(sizeof(int*)*mp);
    Bnum[1][n] = (int**)malloc(sizeof(int*)*mp);
    Cnum[n]  = (int**)malloc(sizeof(int*)*mp);

    Bi[0][n] = (int***)malloc(sizeof(int**)*mp);
    Bi[1][n] = (int***)malloc(sizeof(int**)*mp);
    Ci[n]  = (int***)malloc(sizeof(int**)*mp);

    Bx[0][n] = (double complex***)malloc(sizeof(double complex**)*mp);
    Bx[1][n] = (double complex***)malloc(sizeof(double complex**)*mp);
    Cx[n]  = (double complex***)malloc(sizeof(double complex**)*mp);

    IBx[0][n] = (double complex***)malloc(sizeof(double complex**)*mp);
    IBx[1][n] = (double complex***)malloc(sizeof(double complex**)*mp);
    ICx[n]  = (double complex***)malloc(sizeof(double complex**)*mp);

    for (m=0; m<mp; m++){

      Bnum[0][n][m] = (int*)malloc(sizeof(int)*(EindC[n][m]-SindC[n][m]+1));
      Bnum[1][n][m] = (int*)malloc(sizeof(int)*(EindC[n][m]-SindC[n][m]+1));
      Cnum[n][m]  = (int*)malloc(sizeof(int)*(EindC[n][m]-SindC[n][m]+1));

      Bi[0][n][m] = (int**)malloc(sizeof(int*)*(EindC[n][m]-SindC[n][m]+1));
      Bi[1][n][m] = (int**)malloc(sizeof(int*)*(EindC[n][m]-SindC[n][m]+1));
      Ci[n][m]  = (int**)malloc(sizeof(int*)*(EindC[n][m]-SindC[n][m]+1));

      Bx[0][n][m] = (double complex**)malloc(sizeof(double complex*)*(EindC[n][m]-SindC[n][m]+1));
      Bx[1][n][m] = (double complex**)malloc(sizeof(double complex*)*(EindC[n][m]-SindC[n][m]+1));
      Cx[n][m]  = (double complex**)malloc(sizeof(double complex*)*(EindC[n][m]-SindC[n][m]+1));

      IBx[0][n][m] = (double complex**)malloc(sizeof(double complex*)*(EindC[n][m]-SindC[n][m]+1));
      IBx[1][n][m] = (double complex**)malloc(sizeof(double complex*)*(EindC[n][m]-SindC[n][m]+1));
      ICx[n][m]  = (double complex**)malloc(sizeof(double complex*)*(EindC[n][m]-SindC[n][m]+1));

      for (i=SindC[n][m]; i<=EindC[n][m]; i++){

        i1 = Index_BC[n][i];      

        nb0 = 0;
        nb1 = 0;
        nc  = 0;

        for (j=A->p[i1]; j<A->p[i1+1]; j++){

          j1 = invp[A->i[j]];

          if      ( SindB[n][2*m  ]<=j1 && j1<=EindB[n][2*m  ] ) nb0++;
          else if ( SindB[n][2*m+1]<=j1 && j1<=EindB[n][2*m+1] ) nb1++;
          else if ( SindC[n][m    ]<=j1 && j1<=EindC[n][m    ] ) nc++;
        }

        Bnum[0][n][m][i-SindC[n][m]] = nb0;
        Bnum[1][n][m][i-SindC[n][m]] = nb1;
        Cnum[n][m][i-SindC[n][m]] = nc;

        Bi[0][n][m][i-SindC[n][m]] = (int*)malloc(sizeof(int)*nb0);
        Bi[1][n][m][i-SindC[n][m]] = (int*)malloc(sizeof(int)*nb1);
        Ci[n][m][i-SindC[n][m]] = (int*)malloc(sizeof(int)*nc);

        Bx[0][n][m][i-SindC[n][m]] = (double complex*)malloc(sizeof(double complex)*nb0);
        Bx[1][n][m][i-SindC[n][m]] = (double complex*)malloc(sizeof(double complex)*nb1);
        Cx[n][m][i-SindC[n][m]] = (double complex*)malloc(sizeof(double complex)*nc);

        IBx[0][n][m][i-SindC[n][m]] = (double complex*)malloc(sizeof(double complex)*nb0);
        IBx[1][n][m][i-SindC[n][m]] = (double complex*)malloc(sizeof(double complex)*nb1);
        ICx[n][m][i-SindC[n][m]] = (double complex*)malloc(sizeof(double complex)*nc);

        for (k=0; k<nb0; k++) IBx[0][n][m][i-SindC[n][m]][k] = 0.0;
        for (k=0; k<nb1; k++) IBx[1][n][m][i-SindC[n][m]][k] = 0.0;
        for (k=0; k<nc; k++)  ICx[n][m][i-SindC[n][m]][k] = 0.0;

        nb0 = 0;
        nb1 = 0;
        nc  = 0;

        for (j=A->p[i1]; j<A->p[i1+1]; j++){

          j1 = invp[A->i[j]];
          i2 = i - SindC[n][m];           

          if      ( SindB[n][2*m]<=j1 && j1<=EindB[n][2*m] ){

            Bi[0][n][m][i2][nb0] = j1 - SindB[n][2*m];
            Bx[0][n][m][i2][nb0] = A->x[j];

            nb0++;
	  }

          else if ( SindB[n][2*m+1]<=j1 && j1<=EindB[n][2*m+1] ){

            Bi[1][n][m][i2][nb1] = j1 - SindB[n][2*m+1];
            Bx[1][n][m][i2][nb1] = A->x[j];

            nb1++;
	  }
          else if ( SindC[n][m]<=j1 && j1<=EindC[n][m] ){

            Ci[n][m][i2][nc] = j1 - SindC[n][m];
            Cx[n][m][i2][nc] = A->x[j];

            nc++;
	  }

        } /* j */
      } /* i */
    } /* m */

    mp /= 2;

  } /* n */

  /***********************************************************
             calculate A^-1 for the smallest blocks
  ***********************************************************/

  Anum = (int**)malloc(sizeof(int*)*mpmax);
  Ai = (int***)malloc(sizeof(int**)*mpmax);
  Ax = (double complex***)malloc(sizeof(double complex**)*mpmax);
  IAx = (double complex***)malloc(sizeof(double complex**)*mpmax);

  for (m=0; m<mpmax; m++){

    nsa = EindB[0][m] - SindB[0][m] + 1;

    Anum[m] = (int*)malloc(sizeof(int)*nsa);
    Ai[m] = (int**)malloc(sizeof(int*)*nsa);
    Ax[m] = (double complex**)malloc(sizeof(double complex*)*nsa);
    IAx[m] = (double complex**)malloc(sizeof(double complex*)*nsa);

    for (i=SindB[0][m]; i<=EindB[0][m]; i++){

      i1 = Index_BC[0][i];
      i2 = i - SindB[0][m];
      numa = 0;

      for (j=A->p[i1]; j<A->p[i1+1]; j++){
        j1 = invp[A->i[j]];
	if ( SindB[0][m]<=j1 && j1<=EindB[0][m] ) numa++;
      }

      Anum[m][i2] = numa;

      Ai[m][i2] = (int*)malloc(sizeof(int)*numa);
      Ax[m][i2] = (double complex*)malloc(sizeof(double complex)*numa);
      IAx[m][i2] = (double complex*)malloc(sizeof(double complex)*numa);

      numa = 0;

      for (j=A->p[i1]; j<A->p[i1+1]; j++){

        j1 = invp[A->i[j]];

	if ( SindB[0][m]<=j1 && j1<=EindB[0][m] ){

          Ai[m][i2][numa] = j1 - SindB[0][m];
          Ax[m][i2][numa] = A->x[j];

          numa++; 
	}
      }
    }
  }

  /* calculate IA */

  mp = mpmax;
  n = 0;
  IA = (dcomplex**)malloc(sizeof(dcomplex*)*mp);
  for (m=0; m<mp; m++){
    nsa = EindB[n][m] - SindB[n][m] + 1;
    IA[m] = (dcomplex*)malloc(sizeof(dcomplex)*nsa*nsa);
    for (i=0; i<nsa*nsa; i++) IA[m][i] = Complex(0.0,0.0);
  }

  mp = mpmax;
  n = 0;
  for (m=0; m<mp; m++){

    nsa = EindB[n][m] - SindB[n][m] + 1;

    for (i=0; i<(EindB[n][m]-SindB[n][m]+1); i++){
      for (j=0; j<Anum[m][i]; j++){

        j2 = Ai[m][i][j];
        IA[m][nsa*j2+i].r = creal(Ax[m][i][j]); 
        IA[m][nsa*j2+i].i = cimag(Ax[m][i][j]); 
      }
    }

    if (0<nsa) Lapack_LU_Zinverse(nsa,IA[m]);

    /* store IA into IAx */

    for (i=0; i<(EindB[n][m]-SindB[n][m]+1); i++){
      for (j=0; j<Anum[m][i]; j++){
        j2 = Ai[m][i][j];
        IAx[m][i][j] = IA[m][nsa*j2+i].r + I*IA[m][nsa*j2+i].i;
      }
    }

  } /* m */

  /***********************************************************
                  allocate IS, Vvec, and Lvec
  ***********************************************************/

  mp = mpmax/2;
  IS = (dcomplex***)malloc(sizeof(dcomplex**)*nl);
  for (n=0; n<nl; n++){
    IS[n] = (dcomplex**)malloc(sizeof(dcomplex*)*mp);
    for (m=0; m<mp; m++){
      nsa = EindC[n][m] - SindC[n][m] + 1;

      IS[n][m] = (dcomplex*)malloc(sizeof(dcomplex)*nsa*nsa);
      for (i=0; i<nsa*nsa; i++) IS[n][m][i] = Complex(0.0,0.0);
    }
    mp /= 2;
  }  

  mp = mpmax/2;
  max_nsc = 0;
  for (n=0; n<nl; n++){
    for (m=0; m<mp; m++){
      nsc = EindC[n][m] - SindC[n][m] + 1;
      if (max_nsc<nsc) max_nsc = nsc;
    }
    mp /= 2;
  }

  mp = mpmax/2;
  max_nsb = 0;
  for (n=0; n<nl; n++){
    for (m=0; m<mp; m++){
      nsb = EindB[n][2*m]-SindB[n][2*m] + 1;
      if (max_nsb<nsb) max_nsb = nsb;

      nsb = EindB[n][2*m+1]-SindB[n][2*m+1] + 1;
      if (max_nsb<nsb) max_nsb = nsb;
    }
    mp /= 2;
  }


  Vvec = (double complex**)malloc(sizeof(double complex*)*max_nsc);
  for (i=0; i<max_nsc; i++){
    Vvec[i] = (double complex*)malloc(sizeof(double complex)*Nmat);
  }

  Vvec2 = (double complex**)malloc(sizeof(double complex*)*Nmat);
  for (i=0; i<Nmat; i++){
    Vvec2[i] = (double complex*)malloc(sizeof(double complex)*max_nsc);
  }

  mp = mpmax/2;
  Lvec[0] = (dcomplex***)malloc(sizeof(dcomplex**)*nl);
  for (n=0; n<nl; n++){
    Lvec[0][n] = (dcomplex**)malloc(sizeof(dcomplex*)*mp);
    for (m=0; m<mp; m++){
      kk = (EindB[n][2*m]-SindB[n][2*m]+1)*(EindC[n][m]-SindC[n][m]+1);
      Lvec[0][n][m] = (dcomplex*)malloc(sizeof(dcomplex)*kk);
    }
    mp /= 2;
  }

  mp = mpmax/2;
  Lvec[1] = (dcomplex***)malloc(sizeof(dcomplex**)*nl);
  for (n=0; n<nl; n++){
    Lvec[1][n] = (dcomplex**)malloc(sizeof(dcomplex*)*mp);
    for (m=0; m<mp; m++){
      kk = (EindB[n][2*m+1]-SindB[n][2*m+1]+1)*(EindC[n][m]-SindC[n][m]+1);
      Lvec[1][n][m] = (dcomplex*)malloc(sizeof(dcomplex)*kk);
    }
    mp /= 2;
  }

  Q0 = (double complex**)malloc(sizeof(double complex*)*max_nsc);
  for (i=0; i<max_nsc; i++){
    Q0[i] = (double complex*)malloc(sizeof(double complex)*max_nsc);
  }

  Q1 = (double complex**)malloc(sizeof(double complex*)*max_nsc);
  for (i=0; i<max_nsc; i++){
    Q1[i] = (double complex*)malloc(sizeof(double complex)*max_nsc);
  }

  dtime(&etime);

  time1 = etime - stime;
  printf("OND time1 =%15.12f\n",time1); fflush(stdout);

  /***********************************************************
         main calculations by the recurrence relations
  ***********************************************************/
 
  mp = mpmax/2;
  m2 = 1;

  for (n1=0; n1<nl; n1++){

    /*******************************
        calculate the inital V0 
    *******************************/

    dtime(&stime);

    for (m=0; m<mpmax/2; m++){
      for (m0=0; m0<2; m0++){

	nsa = EindB[0][2*m+m0] - SindB[0][2*m+m0] + 1;

	m3 = (2*m+m0)/(m2*2);
	m4 = ((2*m+m0)/m2)%2;

        nsc = EindC[n1][m3] - SindC[n1][m3] + 1;

	for (i=0; i<nsa; i++){
	  for (j=0; j<nsc; j++){

  	    csum1 = 0.0;
            
	    for (k=0; k<Bnum[m4][n1][m3][j]; k++){

  	      k1 = Bi[m4][n1][m3][j][k] + SindB[n1][2*m3+m4];
              k2 = k1 - SindB[0][2*m+m0];

              if (SindB[0][2*m+m0]<=k1 && k1<=EindB[0][2*m+m0]){

      	        ctmp1 = IA[2*m+m0][nsa*k2+i].r + I*IA[2*m+m0][nsa*k2+i].i;
  	        csum1 += ctmp1*Bx[m4][n1][m3][j][k];
	      }
	    }

 	    Vvec[j][SindB[0][2*m+m0]+i] = csum1;
	  }
	}     
      } /* m0 */
    } /* m */

    dtime(&etime);
    time2 += etime - stime;

    /*******************************
          recurrence relations
    *******************************/

    dtime(&stime);

    for (m=0; m<mp; m++){

      mp0 = m2; 

      for (n2=0; n2<n1; n2++){
        for (m1=0; m1<mp0; m1++){

          mm = mp0*m + m1;

          dtime(&stime1);

          nsc = EindC[n1][m] - SindC[n1][m] + 1;
          nsr = EindC[n2][mm] - SindC[n2][mm] + 1;

	  for (i=0; i<nsr; i++){
	    for (j=0; j<nsc; j++){

              csum1 = 0.0;

              /* B_{even} V_{even} */

  	      for (k=0; k<Bnum[0][n2][mm][i]; k++){
                k1 = SindB[n2][2*mm] + Bi[0][n2][mm][i][k];   
                csum1 += Bx[0][n2][mm][i][k]*Vvec[j][k1]; 
	      }

              numop += Bnum[0][n2][mm][i]; 

              /* B_{odd} V_{odd} */

  	      for (k=0; k<Bnum[1][n2][mm][i]; k++){
                k1 = SindB[n2][2*mm+1] + Bi[1][n2][mm][i][k];
                csum1 += Bx[1][n2][mm][i][k]*Vvec[j][k1]; 
	      }

              numop += Bnum[1][n2][mm][i]; 

              /* B[C] */

              m5 = 1;
              for (l=0; l<(n1-n2); l++) m5 *= 2;

              m3 = mm/m5;
              m4 = (mm/(m5/2))%2;

  	      for (k=0; k<Bnum[m4][n1][m3][j]; k++){
    	        k1 = Bi[m4][n1][m3][j][k] + SindB[n1][2*m3+m4] - SindC[n2][mm];
                if (k1==i) csum1 -= Bx[m4][n1][m3][j][k]; 
	      }

              /* store the temporal result into Q0 */

              Q0[j][i] = csum1;

	    } /* j */
	  } /* i */

          dtime(&etime1);
          time31 += etime1 - stime1;

          /* calculate Q */

          dtime(&stime1);

	  for (i=0; i<nsr; i++){
	    for (j=0; j<nsc; j++){

              csum1 = 0.0;
              for (k=0; k<nsr; k++){
                csum1 += (IS[n2][mm][i*nsr+k].r+I*IS[n2][mm][i*nsr+k].i)*Q0[j][k];  
              }           

              Q1[j][i] = csum1;
	    }
	  }

          dtime(&etime1);
          time32 += etime1 - stime1;

          /* update V */

          dtime(&stime1);

          for (j=0; j<nsc; j++){

	    /* V_{even} */

            for (i=0; i<(EindB[n2][2*mm]-SindB[n2][2*mm]+1); i++){

              i2 = i*nsr; 
  	      csum0 = 0.0;

              for (k=0; k<nsr; k++){
                csum0 += (Lvec[0][n2][mm][i2+k+0].r+I*Lvec[0][n2][mm][i2+k+0].i)*Q1[j][k+0];
	      }

              Vvec[j][SindB[n2][2*mm]+i] += csum0;
	    }

	    /* V_{odd} */

            for (i=0; i<(EindB[n2][2*mm+1]-SindB[n2][2*mm+1]+1); i++){

              i2 = i*nsr; 
  	      csum0 = 0.0;

              for (k=0; k<nsr; k++){
                csum0 += (Lvec[1][n2][mm][i2+k+0].r+I*Lvec[1][n2][mm][i2+k+0].i)*Q1[j][k+0];
	      }

              Vvec[j][SindB[n2][2*mm+1]+i] += csum0;
	    }

	    /* add the contribution of -Q */

            for (i=0; i<nsr; i++){

              i1 = SindC[n2][mm] + i;
              Vvec[j][i1] = -Q1[j][i];
	    }

	  } /* j */   

          dtime(&etime1);
          time33 += etime1 - stime1;

	} /* m1 */

        mp0 /= 2;

      } /* n2 */

      /*******************************
            store Vvec as Lvec
      *******************************/

      nsc = EindC[n1][m] - SindC[n1][m] + 1;

      for (i=SindB[n1][2*m]; i<=EindB[n1][2*m]; i++){
        for (j=0; j<nsc; j++){
	  Lvec[0][n1][m][ (i-SindB[n1][2*m])*nsc+j ].r = creal(Vvec[j][i]);
	  Lvec[0][n1][m][ (i-SindB[n1][2*m])*nsc+j ].i = cimag(Vvec[j][i]);
	}        
      }    

      for (i=SindB[n1][2*m+1]; i<=EindB[n1][2*m+1]; i++){
        for (j=0; j<nsc; j++){
	  Lvec[1][n1][m][ (i-SindB[n1][2*m+1])*nsc+j ].r = creal(Vvec[j][i]);
	  Lvec[1][n1][m][ (i-SindB[n1][2*m+1])*nsc+j ].i = cimag(Vvec[j][i]);
	}        
      }    
    } /* m */

    dtime(&etime);
    time3 += etime - stime;

    /*******************************
      calculate the inverse of S
    *******************************/

    dtime(&stime);

    for (m=0; m<mp; m++){

      nsc = EindC[n1][m] - SindC[n1][m] + 1;

      /* initialize IS */

      for (i=0; i<nsc; i++){
	for (j=0; j<nsc; j++){
	  IS[n1][m][j*nsc+i] = Complex(0.0,0.0);
	}
      }

      /* put C to IS */

      for (i=0; i<nsc; i++){
	for (j=0; j<Cnum[n1][m][i]; j++){
	  j1 = Ci[n1][m][i][j];
	  IS[n1][m][j1*nsc+i].r = creal(Cx[n1][m][i][j]);
	  IS[n1][m][j1*nsc+i].i = cimag(Cx[n1][m][i][j]);
	}
      }

      /* calculate the inner products */

      dtime(&stime1);

      for (i=0; i<nsc; i++){
	for (j=0; j<nsc; j++){

	  csum1 = 0.0;

          /* contribution of even m */

	  for (k=0; k<Bnum[0][n1][m][i]; k++){

	    k1 = SindB[n1][2*m] + Bi[0][n1][m][i][k];
	    csum1 += Bx[0][n1][m][i][k]*Vvec[j][k1];
	  }

          /* contribution of odd m */

	  for (k=0; k<Bnum[1][n1][m][i]; k++){

	    k1 = SindB[n1][2*m+1] + Bi[1][n1][m][i][k];
	    csum1 += Bx[1][n1][m][i][k]*Vvec[j][k1];
	  }

	  IS[n1][m][j*nsc+i].r -= creal(csum1);
	  IS[n1][m][j*nsc+i].i -= cimag(csum1);
	}
      }
 
     dtime(&etime1);
      time41 += etime1 - stime1;

      dtime(&stime1);

      if (0<nsc) Lapack_LU_Zinverse(nsc,IS[n1][m]);

      dtime(&etime1);
      time42 += etime1 - stime1;

    } /* m */

    dtime(&etime);
    time4 += etime - stime;

    /*****************************************
        update entries in the inverse of A
    *****************************************/

    dtime(&stime);

    for (m=0; m<mp; m++){

      nsc = EindC[n1][m] - SindC[n1][m] + 1;

      printf("n1=%2d m=%2d nsc=%2d\n",n1,m,nsc);

      if (0<nsc){

	/* calculate L_{even}^t S^{-1} */

        dtime(&stime1);

        for (i=0; i<(EindB[n1][2*m]-SindB[n1][2*m]+1); i++){

          i2 = i*nsc;

  	  for (j=0; j<nsc; j++){
             
            j2 = j*nsc;

	    dcsum0.r = 0.0; dcsum0.i = 0.0;
	    dcsum1.r = 0.0; dcsum1.i = 0.0;
	    dcsum2.r = 0.0; dcsum2.i = 0.0;
	    dcsum3.r = 0.0; dcsum3.i = 0.0;
	    dcsum4.r = 0.0; dcsum4.i = 0.0;
    
	    for (k=0; k<nsc-3; k+=4){

              i3 = i2 + k;
              j3 = j2 + k;   

	      dcsum0.r += IS[n1][m][j3+0].r*Lvec[0][n1][m][i3+0].r - IS[n1][m][j3+0].i*Lvec[0][n1][m][i3+0].i;
	      dcsum0.i += IS[n1][m][j3+0].i*Lvec[0][n1][m][i3+0].r + IS[n1][m][j3+0].r*Lvec[0][n1][m][i3+0].i;

	      dcsum1.r += IS[n1][m][j3+1].r*Lvec[0][n1][m][i3+1].r - IS[n1][m][j3+1].i*Lvec[0][n1][m][i3+1].i;
	      dcsum1.i += IS[n1][m][j3+1].i*Lvec[0][n1][m][i3+1].r + IS[n1][m][j3+1].r*Lvec[0][n1][m][i3+1].i;

	      dcsum2.r += IS[n1][m][j3+2].r*Lvec[0][n1][m][i3+2].r - IS[n1][m][j3+2].i*Lvec[0][n1][m][i3+2].i;
	      dcsum2.i += IS[n1][m][j3+2].i*Lvec[0][n1][m][i3+2].r + IS[n1][m][j3+2].r*Lvec[0][n1][m][i3+2].i;

	      dcsum3.r += IS[n1][m][j3+3].r*Lvec[0][n1][m][i3+3].r - IS[n1][m][j3+3].i*Lvec[0][n1][m][i3+3].i;
	      dcsum3.i += IS[n1][m][j3+3].i*Lvec[0][n1][m][i3+3].r + IS[n1][m][j3+3].r*Lvec[0][n1][m][i3+3].i;
	    }
        
	    for ( ; k<nsc; k++){

              i3 = i2 + k;
              j3 = j2 + k;   

	      dcsum4.r += IS[n1][m][j3].r*Lvec[0][n1][m][i3].r - IS[n1][m][j3].i*Lvec[0][n1][m][i3].i;
	      dcsum4.i += IS[n1][m][j3].i*Lvec[0][n1][m][i3].r + IS[n1][m][j3].r*Lvec[0][n1][m][i3].i;
	    }

	    Vvec2[i+SindB[n1][2*m]][j] = (dcsum0.r+dcsum1.r+dcsum2.r+dcsum3.r+dcsum4.r)
                                      +I*(dcsum0.i+dcsum1.i+dcsum2.i+dcsum3.i+dcsum4.i);
	  }        

          numop += nsc*nsc; 
	}    

	/* calculate L_{odd}^t S^{-1} */

        for (i=0; i<(EindB[n1][2*m+1]-SindB[n1][2*m+1]+1); i++){

          i2 = i*nsc;

	  for (j=0; j<nsc; j++){

            j2 = j*nsc;

	    dcsum0.r = 0.0; dcsum0.i = 0.0;
	    dcsum1.r = 0.0; dcsum1.i = 0.0;
	    dcsum2.r = 0.0; dcsum2.i = 0.0;
	    dcsum3.r = 0.0; dcsum3.i = 0.0;
	    dcsum4.r = 0.0; dcsum4.i = 0.0;

	    for (k=0; k<nsc-3; k+=4){

              j3 = j2 + k;   
              i3 = i2 + k;

	      dcsum0.r += IS[n1][m][j3+0].r*Lvec[1][n1][m][i3+0].r - IS[n1][m][j3+0].i*Lvec[1][n1][m][i3+0].i;
	      dcsum0.i += IS[n1][m][j3+0].i*Lvec[1][n1][m][i3+0].r + IS[n1][m][j3+0].r*Lvec[1][n1][m][i3+0].i;

	      dcsum1.r += IS[n1][m][j3+1].r*Lvec[1][n1][m][i3+1].r - IS[n1][m][j3+1].i*Lvec[1][n1][m][i3+1].i;
	      dcsum1.i += IS[n1][m][j3+1].i*Lvec[1][n1][m][i3+1].r + IS[n1][m][j3+1].r*Lvec[1][n1][m][i3+1].i;

	      dcsum2.r += IS[n1][m][j3+2].r*Lvec[1][n1][m][i3+2].r - IS[n1][m][j3+2].i*Lvec[1][n1][m][i3+2].i;
	      dcsum2.i += IS[n1][m][j3+2].i*Lvec[1][n1][m][i3+2].r + IS[n1][m][j3+2].r*Lvec[1][n1][m][i3+2].i;

	      dcsum3.r += IS[n1][m][j3+3].r*Lvec[1][n1][m][i3+3].r - IS[n1][m][j3+3].i*Lvec[1][n1][m][i3+3].i;
	      dcsum3.i += IS[n1][m][j3+3].i*Lvec[1][n1][m][i3+3].r + IS[n1][m][j3+3].r*Lvec[1][n1][m][i3+3].i;
	    }

	    for ( ; k<nsc; k++){

              j3 = j2 + k;   
              i3 = i2 + k;

	      dcsum4.r += IS[n1][m][j3].r*Lvec[1][n1][m][i3].r - IS[n1][m][j3].i*Lvec[1][n1][m][i3].i;
	      dcsum4.i += IS[n1][m][j3].i*Lvec[1][n1][m][i3].r + IS[n1][m][j3].r*Lvec[1][n1][m][i3].i;
	    }

	    Vvec2[i+SindB[n1][2*m+1]][j] = (dcsum0.r+dcsum1.r+dcsum2.r+dcsum3.r+dcsum4.r)
                                        +I*(dcsum0.i+dcsum1.i+dcsum2.i+dcsum3.i+dcsum4.i);
	  }        

          numop += nsc*nsc; 
	}    

        dtime(&etime1);
        time51 += etime1 - stime1;

	/* update IAx */

        dtime(&stime1);

	n2 = 0;

	for (m1=0; m1<m2; m1++){

	  m3 = m2*m + m1;

	  /* for even */  

	  if (m2==1) k2 = 0;
	  else       k2 = ((m2/2)<=m1);

	  for (i=SindB[n2][2*m3]; i<=EindB[n2][2*m3]; i++){

	    i2 = i - SindB[n2][2*m3];

	    for (j=0; j<Anum[2*m3][i2]; j++){

	      j2 = Ai[2*m3][i2][j] + SindB[0][2*m3] - SindB[n1][2*m+k2];
	      csum1 = 0.0;

              nsc2 = EindC[n1][m]-SindC[n1][m]+1;

	      for (k=0; k<nsc2; k++){
                ctmp2 = Lvec[k2][n1][m][j2*nsc2+k].r + I*Lvec[k2][n1][m][j2*nsc2+k].i;
		csum1 += Vvec2[i][k]*ctmp2;
	      }

	      IAx[2*m3][i2][j] += csum1; 
	    }      
	  }

	  /* for odd */  

	  if (m2==1) k2 = 1;
	  else       k2 = ((m2/2)<=m1);

	  for (i=SindB[n2][2*m3+1]; i<=EindB[n2][2*m3+1]; i++){

	    i2 = i - SindB[n2][2*m3+1];

	    for (j=0; j<Anum[2*m3+1][i2]; j++){

	      j2 = Ai[2*m3+1][i2][j] + SindB[0][2*m3+1] - SindB[n1][2*m+k2];
	      csum1 = 0.0;

              nsc2 = EindC[n1][m] - SindC[n1][m] + 1;

	      for (k=0; k<nsc2; k++){
                ctmp2 = Lvec[k2][n1][m][j2*nsc2+k].r + I*Lvec[k2][n1][m][j2*nsc2+k].i;
		csum1 += Vvec2[i][k]*ctmp2;
	      }

	      IAx[2*m3+1][i2][j] += csum1; 
	    }      
	  }

	} /* m1 */

        dtime(&etime1);
        time52 += etime1 - stime1;

	/* update the outer IBx */

        dtime(&stime1);

	for (i=0; i<(EindC[n1][m]-SindC[n1][m]+1); i++){

	  /* for even */  

	  for (j=0; j<Bnum[0][n1][m][i]; j++){
             
	    j2 = SindB[n1][2*m] + Bi[0][n1][m][i][j];
	    IBx[0][n1][m][i][j] = -Vvec2[j2][i];
	  }

	  /* for odd */  

	  for (j=0; j<Bnum[1][n1][m][i]; j++){
             
	    j2 = SindB[n1][2*m+1] + Bi[1][n1][m][i][j];
	    IBx[1][n1][m][i][j] = -Vvec2[j2][i];

	  } /* j */
	} /* i */             

        dtime(&etime1);
        time53 += etime1 - stime1;

	/* update the inner IBx */

        dtime(&stime1);

	mp0 = m2;
      
	for (n2=0; n2<n1; n2++){
	  for (m1=0; m1<mp0; m1++){

	    m3 = mp0*m + m1;
	    k2 = ((mp0/2)<=m1);

	    for (i=SindC[n2][m3]; i<=EindC[n2][m3]; i++){

	      i2 = i - SindC[n2][m3];

	      /* for even */  

	      for (j=0; j<Bnum[0][n2][m3][i2]; j++){

		j2 = Bi[0][n2][m3][i2][j] + SindB[n2][2*m3] - SindB[n1][2*m+k2]; 
		csum1 = 0.0;

                nsc2 = EindC[n1][m]-SindC[n1][m]+1;

		for (k=0; k<nsc2; k++){

                  ctmp1 = Lvec[k2][n1][m][j2*nsc2+k].r + I*Lvec[k2][n1][m][j2*nsc2+k].i;
		  csum1 += Vvec2[i][k]*ctmp1;
		}

		IBx[0][n2][m3][i2][j] += csum1;
	      }

	      /* for odd */  

	      for (j=0; j<Bnum[1][n2][m3][i2]; j++){

		j2 = Bi[1][n2][m3][i2][j] + SindB[n2][2*m3+1] - SindB[n1][2*m+k2]; 
		csum1 = 0.0;

                nsc2 = EindC[n1][m]-SindC[n1][m]+1;  

		for (k=0; k<nsc2; k++){
                  ctmp1 = Lvec[k2][n1][m][j2*nsc2+k].r + I*Lvec[k2][n1][m][j2*nsc2+k].i;
		  csum1 += Vvec2[i][k]*ctmp1;
		}

		IBx[1][n2][m3][i2][j] += csum1;
	      }
	    }

	  } /* m1 */             

	  mp0 /= 2;

	} /* n2 */

        dtime(&etime1);
        time54 += etime1 - stime1;

	/* update ICx */

        dtime(&stime1);

	mp0 = m2;

	for (n2=0; n2<n1; n2++){
	  for (m1=0; m1<mp0; m1++){

	    m3 = mp0*m + m1;
	    k2 = ((mp0/2)<=m1);

	    for (i=SindC[n2][m3]; i<=EindC[n2][m3]; i++){

	      i2 = i - SindC[n2][m3];

	      for (j=0; j<Cnum[n2][m3][i2]; j++){

		j2 = Ci[n2][m3][i2][j] + SindC[n2][m3] - SindB[n1][2*m+k2]; 
		csum1 = 0.0;

                nsc2 = EindC[n1][m]-SindC[n1][m]+1;

		for (k=0; k<nsc2; k++){
                  ctmp1 = Lvec[k2][n1][m][j2*nsc2+k].r + I*Lvec[k2][n1][m][j2*nsc2+k].i;
		  csum1 += Vvec2[i][k]*ctmp1;
		}

		ICx[n2][m3][i2][j] += csum1; 

	      } /* j */
	    } /* i */          
	  } /* m1 */

	  mp0 /= 2;

	} /* n2 */


        dtime(&etime1);
        time55 += etime1 - stime1;

        nsc = EindC[n1][m] - SindC[n1][m] + 1;

	for (i=0; i<nsc; i++){
	  for (j=0; j<Cnum[n1][m][i]; j++){

	    j2 = Ci[n1][m][i][j];
	    ICx[n1][m][i][j] = IS[n1][m][j2*nsc+i].r + I*IS[n1][m][j2*nsc+i].i;
	  }
	}

      } /* if (0<nsc) */

    } /* m */                    

    dtime(&etime);
    time5 += etime - stime;
  
    /* update m2 and mp */

    m2 *= 2;
    mp /= 2;

  } /* n1 */

  /***********************************************************
              add the part of inverse to DMall
  ***********************************************************/

  dtime(&stime);

  mp = mpmax/2;
  for (n=0; n<nl; n++){
    for (m=0; m<mp; m++){
      for (i=SindC[n][m]; i<=EindC[n][m]; i++){
        for (j=0; j<Bnum[0][n][m][i-SindC[n][m]]; j++){

          j1 = Bi[0][n][m][i-SindC[n][m]][j];
          j2 = j1 + SindB[n][2*m];
          i2 = i - SindC[n][m];

          i3 = Index_BC[nl][i];          
          j3 = Index_BC[nl][j2]; 
          k = conv_ind2[i3][j3];   
  	  A->x[k] = creal(weight*IBx[0][n][m][i2][j]);

          k = conv_ind2[j3][i3];   
  	  A->x[k] = creal(weight*IBx[0][n][m][i2][j]);
	}
      }
    }

    mp /= 2;
  }

  mp = mpmax/2;
  for (n=0; n<nl; n++){
    for (m=0; m<mp; m++){
      for (i=SindC[n][m]; i<=EindC[n][m]; i++){
        for (j=0; j<Bnum[1][n][m][i-SindC[n][m]]; j++){

          j1 = Bi[1][n][m][i-SindC[n][m]][j];
          j2 = j1 + SindB[n][2*m+1];
          i2 = i - SindC[n][m];

          i3 = Index_BC[nl][i];          
          j3 = Index_BC[nl][j2]; 
          k = conv_ind2[i3][j3];   
  	  A->x[k] = creal(weight*IBx[1][n][m][i2][j]);

          k = conv_ind2[j3][i3];   
  	  A->x[k] = creal(weight*IBx[1][n][m][i2][j]);
	}
      }
    }

    mp /= 2;
  }


  mp = mpmax/2;
  for (n=0; n<nl; n++){
    for (m=0; m<mp; m++){
      for (i=SindC[n][m]; i<=EindC[n][m]; i++){
        for (j=0; j<Cnum[n][m][i-SindC[n][m]]; j++){

          j1 = Ci[n][m][i-SindC[n][m]][j];
          j2 = j1 + SindC[n][m];
          i2 = i - SindC[n][m];

          i3 = Index_BC[nl][i];          
          j3 = Index_BC[nl][j2]; 
          k = conv_ind2[i3][j3];   

  	  A->x[k] = creal(weight*ICx[n][m][i2][j]);
	}
      }
    }

    mp /= 2;
  }

  mp = mpmax;
  for (m=0; m<mp; m++){
    for (i=SindB[0][m]; i<=EindB[0][m]; i++){
      for (j=0; j<Anum[m][i-SindB[0][m]]; j++){

	j1 = Ai[m][i-SindB[0][m]][j];
	j2 = j1 + SindB[0][m];
	i2 = i - SindB[0][m];

	i3 = Index_BC[nl][i];          
	j3 = Index_BC[nl][j2]; 

	k = conv_ind2[i3][j3]; 

	A->x[k] = creal(weight*IAx[m][i2][j]);
      }
    }
  }


  for (i=0; i<TNZE; i++){

    row = conv_index_row[i];                /* row index    */
    col = conv_index_col[i];                /* column index */
    k = conv_ind2[row][col];

    DMall[ChemP_point][spin][i] += A->x[k];
  }  

  dtime(&etime);
  time6 = etime - stime;

  printf("OND numop=%2d\n",numop); fflush(stdout);
  printf("OND time2 =%15.12f\n",time2); fflush(stdout);
  printf("OND time3 =%15.12f\n",time3); fflush(stdout);
  printf("OND time31=%15.12f\n",time31); fflush(stdout);
  printf("OND time32=%15.12f\n",time32); fflush(stdout);
  printf("OND time33=%15.12f\n",time33); fflush(stdout);
  printf("OND time4 =%15.12f\n",time4); fflush(stdout);
  printf("OND time41=%15.12f\n",time41); fflush(stdout);
  printf("OND time42=%15.12f\n",time42); fflush(stdout);
  printf("OND time5 =%15.12f\n",time5); fflush(stdout);
  printf("OND time51=%15.12f\n",time51); fflush(stdout);
  printf("OND time52=%15.12f\n",time52); fflush(stdout);
  printf("OND time53=%15.12f\n",time53); fflush(stdout);
  printf("OND time54=%15.12f\n",time54); fflush(stdout);
  printf("OND time55=%15.12f\n",time55); fflush(stdout);
  printf("OND time6 =%15.12f\n",time6); fflush(stdout);

  /* free arrays */

  for (i=0; i<max_nsc; i++){
    free(Q0[i]);
  }
  free(Q0);

  for (i=0; i<max_nsc; i++){
    free(Q1[i]);
  }
  free(Q1);

  for (i=0; i<max_nsc; i++){
    free(Vvec[i]);
  }
  free(Vvec);

  for (i=0; i<Nmat; i++){
    free(Vvec2[i]);
  }
  free(Vvec2);

  mp = mpmax/2;
  for (n=0; n<nl; n++){
    for (m=0; m<mp; m++){
      free(Lvec[0][n][m]);
    }
    free(Lvec[0][n]);
    mp /= 2;
  }
  free(Lvec[0]);

  mp = mpmax/2;
  for (n=0; n<nl; n++){
    for (m=0; m<mp; m++){
      free(Lvec[1][n][m]);
    }
    free(Lvec[1][n]);
    mp /= 2;
  }
  free(Lvec[1]);

  for (m=0; m<mpmax; m++){

    nsa = EindB[0][m] - SindB[0][m] + 1;

    Anum[m] = (int*)malloc(sizeof(int)*nsa);

    for (i=SindB[0][m]; i<=EindB[0][m]; i++){

      i2 = i - SindB[0][m];

      free(Ai[m][i2]);
      free(Ax[m][i2]);
      free(IAx[m][i2]);
    }

    free(Ai[m]);
    free(Ax[m]);
    free(IAx[m]);
    free(Anum[m]);
  }
  free(Ai);
  free(Ax);
  free(IAx);
  free(Anum);

  mp = mpmax;
  for (m=0; m<mp; m++){
    free(IA[m]);
  }
  free(IA);

  mp = mpmax/2;
  for (n=0; n<nl; n++){
    for (m=0; m<mp; m++){
      free(IS[n][m]);
    }
    free(IS[n]);
    mp /= 2;
  }  
  free(IS);

  mp = mpmax/2;

  for (n=0; n<nl; n++){
    for (m=0; m<mp; m++){
      for (i=SindC[n][m]; i<=EindC[n][m]; i++){

        free(IBx[0][n][m][i-SindC[n][m]]);
        free(IBx[1][n][m][i-SindC[n][m]]);
        free(ICx[n][m][i-SindC[n][m]]);

        free(Bx[0][n][m][i-SindC[n][m]]);
        free(Bx[1][n][m][i-SindC[n][m]]);
        free(Cx[n][m][i-SindC[n][m]]);

        free(Bi[0][n][m][i-SindC[n][m]]);
        free(Bi[1][n][m][i-SindC[n][m]]);
        free(Ci[n][m][i-SindC[n][m]]);

      } /* i */

      free(IBx[0][n][m]);
      free(IBx[1][n][m]);
      free(ICx[n][m]);

      free(Bx[0][n][m]);
      free(Bx[1][n][m]);
      free(Cx[n][m]);

      free(Bi[0][n][m]);
      free(Bi[1][n][m]);
      free(Ci[n][m]);

      free(Bnum[0][n][m]);
      free(Bnum[1][n][m]);
      free(Cnum[n][m]);

    } /* m */

    free(IBx[0][n]);
    free(IBx[1][n]);
    free(ICx[n]);

    free(Bx[0][n]);
    free(Bx[1][n]);
    free(Cx[n]);

    free(Bi[0][n]);
    free(Bi[1][n]);
    free(Ci[n]);

    free(Bnum[0][n]);
    free(Bnum[1][n]);
    free(Cnum[n]);

    mp /= 2;

  } /* n */

  free(IBx[0]);
  free(IBx[1]);
  free(ICx);

  free(Bx[0]);
  free(Bx[1]);
  free(Cx);

  free(Bi[0]);
  free(Bi[1]);
  free(Ci);

  free(Bnum[0]);
  free(Bnum[1]);
  free(Cnum);

  free(invp);

  for (n=0; n<(nl+1); n++){
    free(Index_BC[n]);
  }
  free(Index_BC);

  for (n=0; n<(nl+1); n++){
    free(SindB[n]);
  }  
  free(SindB);

  for (n=0; n<(nl+1); n++){
    free(EindB[n]);
  }  
  free(EindB);

  for (n=0; n<nl; n++){
    free(SindC[n]);
  }  
  free(SindC);

  for (n=0; n<nl; n++){
    free(EindC[n]);
  }  
  free(EindC);

  free(adjncy);
  free(xadj);

  cs_cl_free(x);
  cs_cl_free(b);

  cs_cl_nfree (N) ;
  cs_cl_sfree (s) ;
  cs_cl_spfree(A);
}








static void qsort_complex(long n, long *a, double complex *b)
{
  long int i;
  clists *AB;

  AB = (clists*)malloc(sizeof(clists)*n);

  for (i=0; i<n; i++){
    AB[i].a = a[i];     
    AB[i].b = b[i];
  }

  qsort(AB, n, sizeof(clists), (int(*)(const void*, const void*))clists_cmp);

  for (i=0; i<n; i++){
    a[i] = AB[i].a;
    b[i] = AB[i].b;
  }

  free(AB);
}


static int clists_cmp(const clists *x, const clists *y)
{
  return (x->a < y->a ? -1 :
          y->a < x->a ?  1 : 0);
}






static void Partition_System(
    int sizeflag, int *nl0, int *mpmax0,
    int Nmat, int TNZE2, 
    idxtype *xadj, idxtype *adjncy, 
    int **Index_BC,
    int **SindB, int **EindB,
    int **SindC, int **EindC)
{
  int i,j,j1,k1,n0,n1,k;
  int mpmax,Nmat_tmp,i0,p;
  int wgtflag,npart,edgecut;
  int numflag,options[10];
  int avsize_cluster;
  int n,np,nl,m,mp,po;
  int *xadj_tmp,*adjncy_tmp;
  idxtype *perm, *iperm;

  avsize_cluster = 4;
  np = Nmat/avsize_cluster;
  nl = log((double)np)/log(2.0);
  *nl0 = nl; 

  mpmax = 1;
  for (i=0; i<nl; i++) mpmax = 2*mpmax;
  *mpmax0 = mpmax;

  if (sizeflag==1) goto Last_Partition;

  /* allocation of arrays */

  xadj_tmp = (int*)malloc(sizeof(int)*(Nmat+3));
  adjncy_tmp = (int*)malloc(sizeof(int)*(TNZE2+3));

  perm = (idxtype*)malloc(sizeof(idxtype)*(Nmat+3));
  iperm = (idxtype*)malloc(sizeof(idxtype)*(Nmat+3));

  for (i=0; i<Nmat; i++){
    perm[i] = 0;
    iperm[i] = 0;
  }

  /**********************************************
   perform a nested dissection by the bisection 
  **********************************************/

  numflag = 0;
  wgtflag = 0;
  npart = 2;
  options[0] = 0;

  SindB[nl][0] = 0;
  EindB[nl][0] = Nmat-1;

  mp = 1;

  /* n: level of hierarchy */

  for (n=(nl-1); 0<=n; n--){

    /* m: child number */

    for (m=0; m<mp; m++){

      /* construct the adjacency vector */

      xadj_tmp[0] = 0;
      i0 = 0; 

      for (i=SindB[n+1][m]; i<=EindB[n+1][m]; i++){

        j = Index_BC[n+1][i];

        for (k=xadj[j]; k<xadj[j+1]; k++){

          po = 0;
          for (p=SindB[n+1][m]; p<=EindB[n+1][m]; p++){
            if (adjncy[k]==Index_BC[n+1][p]) { k1=p-SindB[n+1][m]; po=1; break; }
	  }

          if (po==1){
            adjncy_tmp[i0] = k1;
            i0++;
	  }
        }

        xadj_tmp[i-SindB[n+1][m]+1] = i0;
      }

      Nmat_tmp = EindB[n+1][m] - SindB[n+1][m] + 1;

      /*
      printf("\n\n");      
      printf("Z1 n=%2d m=%2d Nmat_tmp=%2d\n",n,m,Nmat_tmp);

      for (i=0; i<=Nmat_tmp; i++){
        printf("Z2 n=%2d m=%2d i=%2d xadj_tmp=%2d\n",n,m,i,xadj_tmp[i]);
      }

      for (i=0; i<Nmat_tmp; i++){
        for (j=xadj_tmp[i]; j<xadj_tmp[i+1]; j++){
          printf("Z3 n=%2d m=%2d i=%2d j=%2d adjncy_tmp=%2d\n",n,m,i,j,adjncy_tmp[j]);
	}
      }
      */

      /* partitioning */

      METIS_PartGraphRecursive( &Nmat_tmp, xadj_tmp, adjncy_tmp, NULL, NULL, 
				&wgtflag, &numflag, &npart, 
				options, &edgecut, perm );

      for (i=0; i<Nmat_tmp; i++){
        iperm[i] = perm[i]; 
      }

      /*
 
      for (i=0; i<Nmat_tmp; i++){
	printf("Z4 n=%2d m=%2d i=%2d part=%2d\n",n,m,i,perm[i]);fflush(stdout);
      }
      */

      /* construct Index_BC for the C part */

      EindC[n][m] = EindB[n+1][m]; 
      SindC[n][m] = EindC[n][m] + 1;

      for (i=(Nmat_tmp-1); 0<=i; i--){

        po = 0; 
	for (j=xadj_tmp[i]; j<xadj_tmp[i+1]; j++){
	  j1 = adjncy_tmp[j]; 
	  if (perm[i]!=perm[j1]) { po=1; break; }
	}

        if (po==1){

	  SindC[n][m]--;
	  Index_BC[n][SindC[n][m]] = Index_BC[n+1][SindB[n+1][m]+i];

          iperm[i] = -1;
	}
      }

      /*
      for (i=0; i<Nmat_tmp; i++){
	printf("Z5 n=%2d m=%2d i=%2d iperm=%2d\n",n,m,i,iperm[i]);fflush(stdout);
      }
      */

      /* construct Index_BC for the B parts */

      n0 = 0;
      n1 = 0;
      for (i=0; i<Nmat_tmp; i++){
        if      (iperm[i]==0) n0++;
        else if (iperm[i]==1) n1++;
      }

      SindB[n][2*m] = SindB[n+1][m];
      EindB[n][2*m] = SindB[n][2*m] + n0 - 1;

      SindB[n][2*m+1] = EindB[n][2*m] + 1;
      EindB[n][2*m+1] = SindB[n][2*m+1] + n1 - 1;

      /*
      printf("n=%2d m=%2d SindB0=%2d EindB0=%2d\n",n,m,SindB[n][2*m],EindB[n][2*m]);
      printf("n=%2d m=%2d SindB1=%2d EindB1=%2d\n",n,m,SindB[n][2*m+1],EindB[n][2*m+1]);
      */

      n0 = 0;
      n1 = 0;
      for (i=0; i<Nmat_tmp; i++){

        if (iperm[i]==0){
          Index_BC[n][SindB[n][2*m]+n0] = Index_BC[n+1][SindB[n+1][m]+i];
          n0++;
	}
        else if (iperm[i]==1){
          Index_BC[n][SindB[n][2*m+1]+n1] = Index_BC[n+1][SindB[n+1][m]+i];
          n1++;
	}
      }

    } /* m */

    mp *= 2; 

  } /* n */

  /********************************************************
   reorder Index_BC so that the order of the grobal index 
   can be equivalent in all the Index_BCs 
  ********************************************************/

  mp = mpmax;

  for (n=0; n<nl; n++){

    for (m=0; m<mp; m++){
      for (i=SindB[n][m]; i<=EindB[n][m]; i++){
	Index_BC[n+1][i] = Index_BC[n][i];
      }
    }

    for (m=0; m<(mp/2); m++){
      for (i=SindC[n][m]; i<=EindC[n][m]; i++){
	Index_BC[n+1][i] = Index_BC[n][i];
      }
    }

    mp /= 2;
  }  


  printf("\n\n");
  printf("nl=%2d mpmax=%2d\n",nl,mpmax);

  /*
  mp = mpmax;
  for (n=0; n<(nl+1); n++){

    for (m=0; m<mp; m++){
      for (i=SindB[n][m]; i<=EindB[n][m]; i++){
        printf("B n=%2d m=%2d i=%2d Index_BC=%2d\n",n,m,i,Index_BC[n][i]);
      }
    }

    mp /= 2;
  }

  mp = mpmax/2;
  for (n=0; n<nl; n++){

    for (m=0; m<mp; m++){
      for (i=SindC[n][m]; i<=EindC[n][m]; i++){
        printf("C n=%2d m=%2d i=%2d Index_BC=%2d\n",n,m,i,Index_BC[n][i]);
      }
    }

    mp /= 2;
  }
  */

  /* freeing of arrays */

  free(iperm);
  free(perm);
  free(xadj_tmp);
  free(adjncy_tmp);

Last_Partition:

}
