/**********************************************************************
  Total_Energy.c:

     Total_Energy.c is a subrutine to calculate the total energy

  Log of Total_Energy.c:

     22/Nov/2001  Released by T.Ozaki
     19/Feb/2006  The subroutine name 'Correction_Energy' was changed 
                  to 'Total_Energy'

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "openmx_common.h"

#define  measure_time   0

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


static double Calc_Ecore();
static double Calc_EH0(int MD_iter);
static double Calc_Ekin();
static double Calc_Ena();
static double Calc_Enl();
static void Calc_EXC_EH1(double ECE[]);
static double Calc_Ehub();  /* --- added by MJ  */
static double Calc_EdftD(); /* added by okuno */
static void EH0_TwoCenter(int Gc_AN, int h_AN, double VH0ij[4]);
static void EH0_TwoCenter_at_Cutoff(int wan1, int wan2, double VH0ij[4]);

/* for OpenMP */
int OneD_Nloop,*OneD2Mc_AN,*OneD2h_AN;


double Total_Energy(int MD_iter, double *****CDM, double ECE[])
{ 
  double time0;
  double TStime,TEtime;
  int numprocs,myid;
  int Mc_AN,Gc_AN,h_AN;
  double stime,etime;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  dtime(&TStime);

  /****************************************************
   For OpenMP:
   making of arrays of the one-dimensionalized loop
  ****************************************************/

  OneD_Nloop = 0;
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];    
    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
      OneD_Nloop++;
    }
  }  

  OneD2Mc_AN = (int*)malloc(sizeof(int)*(OneD_Nloop+1));
  OneD2h_AN = (int*)malloc(sizeof(int)*(OneD_Nloop+1));

  OneD_Nloop = 0;
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];    
    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
      OneD2Mc_AN[OneD_Nloop] = Mc_AN; 
      OneD2h_AN[OneD_Nloop] = h_AN; 
      OneD_Nloop++;
    }
  }

  /****************************************************
               core-core repulsion energy
  ****************************************************/

  dtime(&stime);

  ECE[0] = Calc_Ecore();

  dtime(&etime);
  if(myid==0 && measure_time){
    printf("Time for Ecore=%18.5f\n",etime-stime);fflush(stdout);
  } 
  
  /****************************************************
              EH0 = -1/2\int n^a(r) V^a_H dr
  ****************************************************/

  dtime(&stime);

  ECE[1] = Calc_EH0(MD_iter);

  dtime(&etime);
  if(myid==0 && measure_time){
    printf("Time for EH0=%18.5f\n",etime-stime);fflush(stdout);
  } 

  /****************************************************
                    kinetic energy
  ****************************************************/

  dtime(&stime);

  if (F_Kin_flag==1)  ECE[2] = Calc_Ekin();

  dtime(&etime);
  if(myid==0 && measure_time){
    printf("Time for Ekin=%18.5f\n",etime-stime);fflush(stdout);
  } 

  /****************************************************
              neutral atom potential energy
  ****************************************************/

  dtime(&stime);

  if (F_VNA_flag==1 && ProExpn_VNA==1)  ECE[3] = Calc_Ena();

  dtime(&etime);
  if(myid==0 && measure_time){
    printf("Time for Ena=%18.5f\n",etime-stime);fflush(stdout);
  } 

  /****************************************************
           non-local pseudo potential energy
  ****************************************************/

  dtime(&stime);

  if (F_NL_flag==1)  ECE[4] = Calc_Enl();

  dtime(&etime);
  if(myid==0 && measure_time){
    printf("Time for Enl=%18.5f\n",etime-stime);fflush(stdout);
  } 

  /****************************************************
     EXC = \sum_{\sigma}
           n_{\sigma}(\epsilon_{xc}-\mu_{xc,\sigma})

     EH1 = -1/2\int {n(r)+n^a(r)} \delta V_H dr
  ****************************************************/

  dtime(&stime);

  Calc_EXC_EH1(ECE);

  dtime(&etime);
  if(myid==0 && measure_time){
    printf("Time for EXC_EH1=%18.5f\n",etime-stime);fflush(stdout);
  } 

  /****************************************************
    LDA+U energy   --- added by MJ
  ****************************************************/

  if (F_U_flag==1){
    if (Hub_U_switch==1)  ECE[8] = Calc_Ehub();
    else                  ECE[8] = 0.0;
  }

  /****************************************************

  ****************************************************/

  if(F_dftD_flag==1){
    if(dftD_switch==1) ECE[13] = Calc_EdftD();
    else               ECE[13] = 0.0;
  }

  /****************************************************
   freeing of arrays
  ****************************************************/

  free(OneD2Mc_AN);
  free(OneD2h_AN);

  /* computational time */

  dtime(&TEtime);
  time0 = TEtime - TStime;
  return time0;
}








double Calc_Ekin()
{
  /****************************************************
          calculate the kinetic energy, Ekin
  ****************************************************/

  int i,j,spin;
  int Mc_AN,Gc_AN,Cwan,Rn,h_AN,Gh_AN,Hwan;
  double My_Ekin,Ekin,Zc,Zh,dum,dum2;
  int numprocs,myid;
  double Stime_atom, Etime_atom;
  double *My_Ekin_threads;
  /* for OpenMP */
  int OMPID,Nthrds,Nthrds0,Nprocs,Nloop;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /* get Nthrds0 */  
#pragma omp parallel shared(Nthrds0)
  {
    Nthrds0 = omp_get_num_threads();
  }

  /* allocation of array */
  My_Ekin_threads = (double*)malloc(sizeof(double)*Nthrds0);
  for (Nloop=0; Nloop<Nthrds0; Nloop++) My_Ekin_threads[Nloop] = 0.0;

  if (SpinP_switch==0 || SpinP_switch==1){

    for (spin=0; spin<=SpinP_switch; spin++){

#pragma omp parallel shared(SpinP_switch,time_per_atom,spin,CntH0,H0,DM,My_Ekin_threads,Spe_Total_CNO,Cnt_switch,natn,WhatSpecies,M2G,OneD2h_AN,OneD2Mc_AN,OneD_Nloop) private(OMPID,Nthrds,Nprocs,Nloop,Stime_atom,Mc_AN,h_AN,Gc_AN,Cwan,Gh_AN,Hwan,i,j,Etime_atom)
      {

	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	/* one-dimensionalized loop */

	for (Nloop=OMPID*OneD_Nloop/Nthrds; Nloop<(OMPID+1)*OneD_Nloop/Nthrds; Nloop++){

	  dtime(&Stime_atom);

	  /* get Mc_AN and h_AN */

	  Mc_AN = OneD2Mc_AN[Nloop];
	  h_AN  = OneD2h_AN[Nloop];

	  /* set data on Mc_AN */

	  Gc_AN = M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];

	  /* set data on h_AN */

	  Gh_AN = natn[Gc_AN][h_AN];
	  Hwan = WhatSpecies[Gh_AN];

	  if (Cnt_switch==0){
	    for (i=0; i<Spe_Total_CNO[Cwan]; i++){
	      for (j=0; j<Spe_Total_CNO[Hwan]; j++){
		My_Ekin_threads[OMPID] += DM[0][spin][Mc_AN][h_AN][i][j]*H0[0][Mc_AN][h_AN][i][j];
	      }
	    }
	  }

	  else{

	    for (i=0; i<Spe_Total_CNO[Cwan]; i++){
	      for (j=0; j<Spe_Total_CNO[Hwan]; j++){
		My_Ekin_threads[OMPID] += DM[0][spin][Mc_AN][h_AN][i][j]*CntH0[0][Mc_AN][h_AN][i][j];
	      }
	    }
	  }

	  dtime(&Etime_atom);
	  time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
	}

        if (SpinP_switch==0) My_Ekin_threads[OMPID] = 2.0*My_Ekin_threads[OMPID];

      } /* #pragma omp parallel */
    } /* spin */ 

  }
  else if (SpinP_switch==3){

#pragma omp parallel shared(time_per_atom,H0,DM,My_Ekin_threads,Spe_Total_CNO,natn,WhatSpecies,M2G,OneD2h_AN,OneD2Mc_AN,OneD_Nloop) private(OMPID,Nthrds,Nprocs,Nloop,Stime_atom,Mc_AN,Gc_AN,Cwan,h_AN,Gh_AN,Hwan,i,j,Etime_atom)
    {

      /* get info. on OpenMP */ 

      OMPID = omp_get_thread_num();
      Nthrds = omp_get_num_threads();
      Nprocs = omp_get_num_procs();

      /* one-dimensionalized loop */

      for (Nloop=OMPID*OneD_Nloop/Nthrds; Nloop<(OMPID+1)*OneD_Nloop/Nthrds; Nloop++){

	dtime(&Stime_atom);

	/* get Mc_AN and h_AN */

	Mc_AN = OneD2Mc_AN[Nloop];
	h_AN  = OneD2h_AN[Nloop];

	/* set data on Mc_AN */

	Gc_AN = M2G[Mc_AN];
	Cwan = WhatSpecies[Gc_AN];

	/* set data on h_AN */

	Gh_AN = natn[Gc_AN][h_AN];
	Hwan = WhatSpecies[Gh_AN];

	for (i=0; i<Spe_Total_CNO[Cwan]; i++){
	  for (j=0; j<Spe_Total_CNO[Hwan]; j++){
	    My_Ekin_threads[OMPID] += (DM[0][0][Mc_AN][h_AN][i][j] + DM[0][1][Mc_AN][h_AN][i][j])*H0[0][Mc_AN][h_AN][i][j];
	  }
	}

	dtime(&Etime_atom);
	time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
      }

    } /* #pragma omp parallel */
  }

  /* sum of My_Ekin_threads */
  My_Ekin = 0.0;
  for (Nloop=0; Nloop<Nthrds0; Nloop++){
    My_Ekin += My_Ekin_threads[Nloop];
  }

  /* sum of My_Ekin */
  MPI_Allreduce(&My_Ekin, &Ekin, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

  /* freeing of array */
  free(My_Ekin_threads);

  return Ekin;  
}





double Calc_Ena()
{
  /****************************************************
     calculate the neutral atom potential energy, Ena
  ****************************************************/

  int i,j,spin;
  int Mc_AN,Gc_AN,Cwan,Rn,h_AN,Gh_AN,Hwan;
  double My_Ena,Ena,Zc,Zh,dum,dum2;
  int numprocs,myid;
  double Stime_atom, Etime_atom;
  double *My_Ena_threads;
  /* for OpenMP */
  int OMPID,Nthrds,Nthrds0,Nprocs,Nloop;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /* get Nthrds0 */  
#pragma omp parallel shared(Nthrds0)
  {
    Nthrds0 = omp_get_num_threads();
  }

  /* allocation of array */
  My_Ena_threads = (double*)malloc(sizeof(double)*Nthrds0);
  for (Nloop=0; Nloop<Nthrds0; Nloop++) My_Ena_threads[Nloop] = 0.0;

  if (SpinP_switch==0 || SpinP_switch==1){

    if (Cnt_switch==1){
      /* temporaly, we borrow the CntH0 matrix */
      Cont_Matrix0(HVNA,CntH0[0]);
    }

    for (spin=0; spin<=SpinP_switch; spin++){

#pragma omp parallel shared(SpinP_switch,time_per_atom,CntH0,HVNA,DM,My_Ena_threads,Spe_Total_CNO,Cnt_switch,natn,WhatSpecies,M2G,OneD2h_AN,OneD2Mc_AN,OneD_Nloop) private(OMPID,Nthrds,Nprocs,Nloop,Stime_atom,Mc_AN,h_AN,Gc_AN,Cwan,Gh_AN,Hwan,i,j,Etime_atom)
      {

	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	/* one-dimensionalized loop */

	for (Nloop=OMPID*OneD_Nloop/Nthrds; Nloop<(OMPID+1)*OneD_Nloop/Nthrds; Nloop++){

	  dtime(&Stime_atom);

	  /* get Mc_AN and h_AN */

	  Mc_AN = OneD2Mc_AN[Nloop];
	  h_AN  = OneD2h_AN[Nloop];

	  /* set data on Mc_AN */

	  Gc_AN = M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];

	  /* set data on h_AN */

	  Gh_AN = natn[Gc_AN][h_AN];
	  Hwan = WhatSpecies[Gh_AN];

	  if (Cnt_switch==0){
	    for (i=0; i<Spe_Total_CNO[Cwan]; i++){
	      for (j=0; j<Spe_Total_CNO[Hwan]; j++){
		My_Ena_threads[OMPID] += DM[0][spin][Mc_AN][h_AN][i][j]*HVNA[Mc_AN][h_AN][i][j];
	      }
	    }
	  }

	  else {

	    for (i=0; i<Spe_Total_CNO[Cwan]; i++){
	      for (j=0; j<Spe_Total_CNO[Hwan]; j++){
		My_Ena_threads[OMPID] += DM[0][spin][Mc_AN][h_AN][i][j]*CntH0[0][Mc_AN][h_AN][i][j];
	      }
	    }
	  }

	  dtime(&Etime_atom);
	  time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
	}

        if (SpinP_switch==0) My_Ena_threads[OMPID] = 2.0*My_Ena_threads[OMPID];

      } /* #pragma omp parallel */
    } /* spin */

  }
  else if (SpinP_switch==3){

#pragma omp parallel shared(time_per_atom,HVNA,DM,Spe_Total_CNO,natn,WhatSpecies,M2G,OneD2h_AN,OneD2Mc_AN,OneD_Nloop) private(OMPID,Nthrds,Nprocs,Nloop,Stime_atom,Mc_AN,Gc_AN,h_AN,Cwan,Gh_AN,Hwan,i,j,Etime_atom)
    {

      /* get info. on OpenMP */ 

      OMPID  = omp_get_thread_num();
      Nthrds = omp_get_num_threads();
      Nprocs = omp_get_num_procs();

      /* one-dimensionalized loop */

      for (Nloop=OMPID*OneD_Nloop/Nthrds; Nloop<(OMPID+1)*OneD_Nloop/Nthrds; Nloop++){

	dtime(&Stime_atom);

	/* get Mc_AN and h_AN */

	Mc_AN = OneD2Mc_AN[Nloop];
	h_AN  = OneD2h_AN[Nloop];

	/* set data on Mc_AN */

	Gc_AN = M2G[Mc_AN];
	Cwan = WhatSpecies[Gc_AN];

	/* set data on h_AN */

	Gh_AN = natn[Gc_AN][h_AN];
	Hwan = WhatSpecies[Gh_AN];

	for (i=0; i<Spe_Total_CNO[Cwan]; i++){
	  for (j=0; j<Spe_Total_CNO[Hwan]; j++){
	    My_Ena_threads[OMPID] += (DM[0][0][Mc_AN][h_AN][i][j]+DM[0][1][Mc_AN][h_AN][i][j])*HVNA[Mc_AN][h_AN][i][j];
	  }
	}

	dtime(&Etime_atom);
	time_per_atom[Gc_AN] += Etime_atom - Stime_atom;

      } /* Nloop */
    } /* #pragma omp parallel */ 
  } /* else if (SpinP_switch==3) */


  /* sum of My_Ena_threads */
  My_Ena = 0.0;
  for (Nloop=0; Nloop<Nthrds0; Nloop++){
    My_Ena += My_Ena_threads[Nloop];
  }

  /* sum of My_Ena */
  MPI_Allreduce(&My_Ena, &Ena, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

  /* freeing of array */
  free(My_Ena_threads);

  return Ena;  
}





double Calc_Enl()
{
  /****************************************************
     calculate the non-local pseudo potential energy
  ****************************************************/

  int i,j,spin;
  int Mc_AN,Gc_AN,Cwan,Rn,h_AN,Gh_AN,Hwan;
  double My_Enl,Enl,Zc,Zh,dum,dum2;
  int numprocs,myid;
  double Stime_atom, Etime_atom;
  double *My_Enl_threads;
  /* for OpenMP */
  int OMPID,Nthrds,Nthrds0,Nprocs,Nloop;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /* get Nthrds0 */  
#pragma omp parallel shared(Nthrds0)
  {
    Nthrds0 = omp_get_num_threads();
  }

  /* allocation of array */
  My_Enl_threads = (double*)malloc(sizeof(double)*Nthrds0);
  for (Nloop=0; Nloop<Nthrds0; Nloop++) My_Enl_threads[Nloop] = 0.0;

  if (SpinP_switch==0 || SpinP_switch==1){

    for (spin=0; spin<=SpinP_switch; spin++){

      if (Cnt_switch==1){
        /* temporaly, borrow the CntH0 matrix */
        Cont_Matrix0(HNL[spin],CntH0[0]);
      }

#pragma omp parallel shared(spin,SpinP_switch,time_per_atom,CntH0,HNL,DM,My_Enl_threads,Spe_Total_CNO,Cnt_switch,natn,WhatSpecies,M2G,OneD2h_AN,OneD2Mc_AN,OneD_Nloop) private(OMPID,Nthrds,Nprocs,Nloop,Stime_atom,Mc_AN,h_AN,Gc_AN,Cwan,Gh_AN,Hwan,Etime_atom,i,j)
      {

	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	/* one-dimensionalized loop */

	for (Nloop=OMPID*OneD_Nloop/Nthrds; Nloop<(OMPID+1)*OneD_Nloop/Nthrds; Nloop++){

	  dtime(&Stime_atom);

	  /* get Mc_AN and h_AN */

	  Mc_AN = OneD2Mc_AN[Nloop];
	  h_AN  = OneD2h_AN[Nloop];

	  /* set data on Mc_AN */

	  Gc_AN = M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];

	  /* set data on h_AN */

	  Gh_AN = natn[Gc_AN][h_AN];
	  Hwan = WhatSpecies[Gh_AN];

	  if (Cnt_switch==0){
	    for (i=0; i<Spe_Total_CNO[Cwan]; i++){
	      for (j=0; j<Spe_Total_CNO[Hwan]; j++){
		My_Enl_threads[OMPID] += DM[0][spin][Mc_AN][h_AN][i][j]*HNL[spin][Mc_AN][h_AN][i][j];
	      }
	    }
	  }

	  else {
	    for (i=0; i<Spe_Total_CNO[Cwan]; i++){
	      for (j=0; j<Spe_Total_CNO[Hwan]; j++){
		My_Enl_threads[OMPID] += DM[0][spin][Mc_AN][h_AN][i][j]*CntH0[0][Mc_AN][h_AN][i][j];
	      }
	    }
	  }

	  dtime(&Etime_atom);
	  time_per_atom[Gc_AN] += Etime_atom - Stime_atom;

	} /* Nloop */

        if (SpinP_switch==0) My_Enl_threads[OMPID] = 2.0*My_Enl_threads[OMPID];

      } /* #pragma omp parallel */
    } /* spin */

  }
  else if (SpinP_switch==3){

#pragma omp parallel shared(time_per_atom,iHNL0,HNL,iDM,DM,My_Enl_threads,Spe_Total_CNO,natn,WhatSpecies,M2G,OneD2h_AN,OneD2Mc_AN,OneD_Nloop) private(OMPID,Nthrds,Nprocs,Nloop,Stime_atom,Etime_atom,Mc_AN,h_AN,Gc_AN,Cwan,Gh_AN,Hwan,i,j)
    {

      /* get info. on OpenMP */ 

      OMPID  = omp_get_thread_num();
      Nthrds = omp_get_num_threads();
      Nprocs = omp_get_num_procs();

      /* one-dimensionalized loop */

      for (Nloop=OMPID*OneD_Nloop/Nthrds; Nloop<(OMPID+1)*OneD_Nloop/Nthrds; Nloop++){

	dtime(&Stime_atom);

	/* get Mc_AN and h_AN */

	Mc_AN = OneD2Mc_AN[Nloop];
	h_AN  = OneD2h_AN[Nloop];

	/* set data on Mc_AN */

	Gc_AN = M2G[Mc_AN];
	Cwan = WhatSpecies[Gc_AN];

	/* set data on h_AN */

	Gh_AN = natn[Gc_AN][h_AN];
	Hwan = WhatSpecies[Gh_AN];

	for (i=0; i<Spe_Total_CNO[Cwan]; i++){
	  for (j=0; j<Spe_Total_CNO[Hwan]; j++){

	    My_Enl_threads[OMPID] += 
	         DM[0][0][Mc_AN][h_AN][i][j]*  HNL[0][Mc_AN][h_AN][i][j]
	      - iDM[0][0][Mc_AN][h_AN][i][j]*iHNL0[0][Mc_AN][h_AN][i][j]
	      +  DM[0][1][Mc_AN][h_AN][i][j]*  HNL[1][Mc_AN][h_AN][i][j]
	      - iDM[0][1][Mc_AN][h_AN][i][j]*iHNL0[1][Mc_AN][h_AN][i][j]
	   + 2.0*DM[0][2][Mc_AN][h_AN][i][j]*  HNL[2][Mc_AN][h_AN][i][j] 
	   - 2.0*DM[0][3][Mc_AN][h_AN][i][j]*iHNL0[2][Mc_AN][h_AN][i][j];
 
	  }
	}

	dtime(&Etime_atom);
	time_per_atom[Gc_AN] += Etime_atom - Stime_atom;

      } /* Nloop */
    } /* #pragma omp parallel */
  }

  /* sum of My_Enl_threads */
  My_Enl = 0.0;
  for (Nloop=0; Nloop<Nthrds0; Nloop++){
    My_Enl += My_Enl_threads[Nloop];
  }

  /* sum of My_Enl */
  MPI_Allreduce(&My_Enl, &Enl, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

  /* freeing of array */
  free(My_Enl_threads);

  return Enl;  
}






double Calc_Ecore()
{
  /****************************************************
                         Ecore
  ****************************************************/

  int Mc_AN,Gc_AN,Cwan,Rn,h_AN,Gh_AN,Hwan;
  double My_Ecore,Ecore,Zc,Zh,dum,dum2;
  double *My_Ecore_threads;
  double dEx,dEy,dEz,r,lx,ly,lz;
  int numprocs,myid;
  double Stime_atom,Etime_atom;
  /* for OpenMP */
  int OMPID,Nthrds,Nthrds0,Nprocs,Nloop;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (myid==Host_ID && 0<level_stdout){
    printf("  Force calculation #6\n");fflush(stdout);
  }

  /* get Nthrds0 */  
#pragma omp parallel shared(Nthrds0)
  {
    Nthrds0 = omp_get_num_threads();
  }

  /* allocation of array */
  My_Ecore_threads = (double*)malloc(sizeof(double)*Nthrds0);
  for (Nloop=0; Nloop<Nthrds0; Nloop++) My_Ecore_threads[Nloop] = 0.0;

#pragma omp parallel shared(level_stdout,time_per_atom,atv,Gxyz,Dis,ncn,natn,FNAN,Spe_Core_Charge,WhatSpecies,M2G,Matomnum,My_Ecore_threads) private(OMPID,Nthrds,Nprocs,Mc_AN,Stime_atom,Gc_AN,Cwan,Zc,dEx,dEy,dEz,h_AN,Gh_AN,Rn,Hwan,Zh,r,lx,ly,lz,dum,dum2,Etime_atom)
  {

    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    for (Mc_AN=(OMPID*Matomnum/Nthrds+1); Mc_AN<((OMPID+1)*Matomnum/Nthrds+1); Mc_AN++){

      dtime(&Stime_atom);

      Gc_AN = M2G[Mc_AN];
      Cwan = WhatSpecies[Gc_AN];
      Zc = Spe_Core_Charge[Cwan];
      dEx = 0.0;
      dEy = 0.0;
      dEz = 0.0;

      for (h_AN=1; h_AN<=FNAN[Gc_AN]; h_AN++){

	Gh_AN = natn[Gc_AN][h_AN];
	Rn = ncn[Gc_AN][h_AN];
	Hwan = WhatSpecies[Gh_AN];
	Zh = Spe_Core_Charge[Hwan];
	r = Dis[Gc_AN][h_AN];

	/* for empty atoms or finite elemens basis */
	if (r<1.0e-10) r = 1.0e-10;

	lx = (Gxyz[Gc_AN][1] - Gxyz[Gh_AN][1] - atv[Rn][1])/r;
	ly = (Gxyz[Gc_AN][2] - Gxyz[Gh_AN][2] - atv[Rn][2])/r;
	lz = (Gxyz[Gc_AN][3] - Gxyz[Gh_AN][3] - atv[Rn][3])/r;
	dum = Zc*Zh/r;
	dum2 = dum/r;
	My_Ecore_threads[OMPID] += dum;
	dEx = dEx - lx*dum2;
	dEy = dEy - ly*dum2;
	dEz = dEz - lz*dum2;
      }

      /****************************************************
                        #6 of force
         Contribution from the core-core repulsions
      ****************************************************/

      Gxyz[Gc_AN][17] += dEx;
      Gxyz[Gc_AN][18] += dEy;
      Gxyz[Gc_AN][19] += dEz;

      if (2<=level_stdout){
	printf("<Total_Ene>  force(6) myid=%2d  Mc_AN=%2d Gc_AN=%2d  %15.12f %15.12f %15.12f\n",
	       myid,Mc_AN,Gc_AN,dEx,dEy,dEz);fflush(stdout);
      }

      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
    }

    My_Ecore_threads[OMPID] *= 0.50;

  } /* #pragma omp parallel */

  My_Ecore = 0.0;
  for (Nloop=0; Nloop<Nthrds0; Nloop++){
    My_Ecore += My_Ecore_threads[Nloop];
  }

  MPI_Allreduce(&My_Ecore, &Ecore, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

  /* freeing of array */
  free(My_Ecore_threads);

  return Ecore;  
}








double Calc_EH0(int MD_iter)
{
  /****************************************************
              EH0 = -1/2\int n^a(r) V^a_H dr
  ****************************************************/

  int Mc_AN,Gc_AN,h_AN,Gh_AN,num,gnum;
  int wan,wan1,wan2,Nd,n1,n2,n3;
  double bc,dx,x,y,z,r1,r2,rho0;
  double Scale_Grid_Ecut;
  double EH0ij[4],My_EH0,EH0,tmp0;
  double *Fx,*Fy,*Fz,*g0;
  double dEx,dEy,dEz,Dx,Sx;
  double Z1,Z2,factor;
  double My_dEx,My_dEy,My_dEz;
  int numprocs,myid,ID;
  double stime,etime;
  double Stime_atom, Etime_atom;
  /* for OpenMP */
  int OMPID,Nthrds,Nthrds0,Nprocs,Nloop;
  double *My_EH0_threads;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
     allocation of arrays:

    double Fx[Matomnum+1];
    double Fy[Matomnum+1];
    doubel Fz[Matomnum+1];
  ****************************************************/

  Fx = (double*)malloc(sizeof(double)*(Matomnum+1));
  Fy = (double*)malloc(sizeof(double)*(Matomnum+1));
  Fz = (double*)malloc(sizeof(double)*(Matomnum+1));

  /****************************************************
             Set of atomic density on grids
  ****************************************************/

  if (MD_iter==1){

    dtime(&stime);

    Scale_Grid_Ecut = 16.0;

    /* estimate the size of an array g0 */

    Max_Nd = 0;
    for (wan=0; wan<SpeciesNum; wan++){
      Spe_Spe2Ban[wan] = wan;
      bc = Spe_Atom_Cut1[wan];
      dx = PI/sqrt(Grid_Ecut*Scale_Grid_Ecut);
      Nd = 2*(int)(bc/dx) + 1;
      if (Max_Nd<Nd) Max_Nd = Nd;
    }

    /* estimate sizes of arrays GridX,Y,Z_EH0, Arho_EH0, and Wt_EH0 */

    Max_TGN_EH0 = 0;
    for (wan=0; wan<SpeciesNum; wan++){

      Spe_Spe2Ban[wan] = wan;
      bc = Spe_Atom_Cut1[wan];
      dx = PI/sqrt(Grid_Ecut*Scale_Grid_Ecut);
      Nd = 2*(int)(bc/dx) + 1;
      dx = 2.0*bc/(double)(Nd-1);
      gnum = Nd*CoarseGL_Mesh;

      if (Max_TGN_EH0<gnum) Max_TGN_EH0 = gnum;

      if (2<=level_stdout){
        printf("<Calc_EH0> A spe=%2d 1D-grids=%2d 3D-grids=%2d\n",wan,Nd,gnum);fflush(stdout);
      }
    }
    
    /* allocation of arrays GridX,Y,Z_EH0, Arho_EH0, and Wt_EH0 */

    Max_TGN_EH0 += 10; 
    Allocate_Arrays(4);

    /* calculate GridX,Y,Z_EH0 and Wt_EH0 */

#pragma omp parallel shared(Max_Nd,level_stdout,TGN_EH0,Wt_EH0,Arho_EH0,GridZ_EH0,GridY_EH0,GridX_EH0,Scale_Grid_Ecut,Grid_Ecut,Spe_Atom_Cut1,dv_EH0,Spe_Spe2Ban,SpeciesNum,CoarseGL_Abscissae,CoarseGL_Weight) private(OMPID,Nthrds,Nprocs,wan,bc,dx,Nd,gnum,n1,n2,n3,x,y,z,tmp0,r1,rho0,g0,Sx,Dx)
    {

      /* allocation of arrays g0 */

      g0 = (double*)malloc(sizeof(double)*Max_Nd);

      /* get info. on OpenMP */ 

      OMPID = omp_get_thread_num();
      Nthrds = omp_get_num_threads();
      Nprocs = omp_get_num_procs();
    
      for (wan=OMPID*SpeciesNum/Nthrds; wan<(OMPID+1)*SpeciesNum/Nthrds; wan++){

	Spe_Spe2Ban[wan] = wan;
	bc = Spe_Atom_Cut1[wan];
	dx = PI/sqrt(Grid_Ecut*Scale_Grid_Ecut);
	Nd = 2*(int)(bc/dx) + 1;
	dx = 2.0*bc/(double)(Nd-1);
	dv_EH0[wan] = dx;

	for (n1=0; n1<Nd; n1++){
	  g0[n1] = dx*(double)n1 - bc;
	}

	gnum = 0; 
        y = 0.0;

        Sx = Spe_Atom_Cut1[wan] + 0.0;
        Dx = Spe_Atom_Cut1[wan] - 0.0;

        for (n3=0; n3<Nd; n3++){

          z = g0[n3];
          tmp0 = z*z;

  	  for (n1=0; n1<CoarseGL_Mesh; n1++){

            x = 0.50*(Dx*CoarseGL_Abscissae[n1] + Sx);
	    r1 = sqrt(x*x + tmp0);

	    GridX_EH0[wan][gnum] = x;
	    GridY_EH0[wan][gnum] = y;
	    GridZ_EH0[wan][gnum] = z;

      	    rho0 = AtomicDenF(wan,r1);
	    Arho_EH0[wan][gnum] = rho0;
   	    Wt_EH0[wan][gnum] = PI*x*CoarseGL_Weight[n1]*Dx;

	    gnum++;  
	  }

	} /* n3 */

	TGN_EH0[wan] = gnum;

	if (2<=level_stdout){
	  printf("<Calc_EH0> B spe=%2d 1D-grids=%2d 3D-grids=%2d\n",wan,Nd,gnum);fflush(stdout);
	}

      } /* wan */

      /* free */
      free(g0);

    } /* #pragma omp parallel */

    dtime(&etime);
    if(myid==0 && measure_time){
      printf("Time for part1 of EH0=%18.5f\n",etime-stime);fflush(stdout);
    } 

  } /* if (MD_iter==1) */

  /****************************************************
    calculation of scaling factors:
  ****************************************************/

  if (MD_iter==1){

    for (wan1=0; wan1<SpeciesNum; wan1++){

      r1 = Spe_Atom_Cut1[wan1];
      Z1 = Spe_Core_Charge[wan1];

      for (wan2=0; wan2<SpeciesNum; wan2++){

	/* EH0_TwoCenter_at_Cutoff is parallelized by OpenMP */      
	EH0_TwoCenter_at_Cutoff(wan1, wan2, EH0ij);

	r2 = Spe_Atom_Cut1[wan2];
	Z2 = Spe_Core_Charge[wan2];
	tmp0 = Z1*Z2/(r1+r2);


	if (1.0e-20<fabs(EH0ij[0])){ 
	  EH0_scaling[wan1][wan2] = tmp0/EH0ij[0];
	}
	else{
	  EH0_scaling[wan1][wan2] = 0.0;
	}

      }
    }
  }

  /****************************************************
                -1/2\int n^a(r) V^a_H dr
  ****************************************************/

  dtime(&stime);

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Fx[Mc_AN] = 0.0;
    Fy[Mc_AN] = 0.0;
    Fz[Mc_AN] = 0.0;
  }

  /* get Nthrds0 */  
#pragma omp parallel shared(Nthrds0)
  {
    Nthrds0 = omp_get_num_threads();
  }

  /* allocation of array */
  My_EH0_threads = (double*)malloc(sizeof(double)*Nthrds0);
  for (Nloop=0; Nloop<Nthrds0; Nloop++) My_EH0_threads[Nloop] = 0.0;

#pragma omp parallel shared(time_per_atom,RMI1,EH0_scaling,natn,FNAN,WhatSpecies,M2G,Matomnum,My_EH0_threads) private(OMPID,Nthrds,Nprocs,Mc_AN,Stime_atom,Gc_AN,wan1,h_AN,Gh_AN,wan2,factor,EH0ij,Etime_atom)
  {

    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();
  
    for (Mc_AN=(OMPID*Matomnum/Nthrds+1); Mc_AN<((OMPID+1)*Matomnum/Nthrds+1); Mc_AN++){

      dtime(&Stime_atom);

      Gc_AN = M2G[Mc_AN];
      wan1 = WhatSpecies[Gc_AN];

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	Gh_AN = natn[Gc_AN][h_AN];
	wan2 = WhatSpecies[Gh_AN];

	if (h_AN==0) factor = 1.0;
	else         factor = EH0_scaling[wan1][wan2];

	EH0_TwoCenter(Gc_AN, h_AN, EH0ij);

	My_EH0_threads[OMPID] -= 0.250*factor*EH0ij[0];
	Fx[Mc_AN] = Fx[Mc_AN] - 0.5*factor*EH0ij[1];
	Fy[Mc_AN] = Fy[Mc_AN] - 0.5*factor*EH0ij[2];
	Fz[Mc_AN] = Fz[Mc_AN] - 0.5*factor*EH0ij[3];

	if (h_AN==0) factor = 1.0;
	else         factor = EH0_scaling[wan2][wan1];

	EH0_TwoCenter(Gh_AN, RMI1[Mc_AN][h_AN][0], EH0ij);

	My_EH0_threads[OMPID] -= 0.250*factor*EH0ij[0];
	Fx[Mc_AN] = Fx[Mc_AN] + 0.5*factor*EH0ij[1];
	Fy[Mc_AN] = Fy[Mc_AN] + 0.5*factor*EH0ij[2];
	Fz[Mc_AN] = Fz[Mc_AN] + 0.5*factor*EH0ij[3];

      } /* h_AN */

      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;

    } /* Mc_AN */
  } /* #pragma omp parallel */

  /* sum of My_EH0_threads */
  My_EH0 = 0.0;
  for (Nloop=0; Nloop<Nthrds0; Nloop++){
    My_EH0 += My_EH0_threads[Nloop];
  }

  /* sum of My_EH0 */
  MPI_Allreduce(&My_EH0, &EH0, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

  /* freeing of array */
  free(My_EH0_threads);

  dtime(&etime);
  if(myid==0 && measure_time){
    printf("Time for part2 of EH0=%18.5f\n",etime-stime);fflush(stdout);
  } 

  /****************************************************
                      #7 of force
     contribution from the Hartree energy between
   the neutral atomic charge and the neutral potential 
  ****************************************************/

  if (myid==Host_ID && 0<level_stdout){
    printf("  Force calculation #7\n");fflush(stdout);
  }

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];

    if (2<=level_stdout){
      printf("<Total_Ene>  force(7) myid=%2d  Mc_AN=%2d Gc_AN=%2d  %15.12f %15.12f %15.12f\n",
              myid,Mc_AN,Gc_AN,Fx[Mc_AN],Fy[Mc_AN],Fz[Mc_AN]);fflush(stdout);
    }

    Gxyz[Gc_AN][17] += Fx[Mc_AN];
    Gxyz[Gc_AN][18] += Fy[Mc_AN];
    Gxyz[Gc_AN][19] += Fz[Mc_AN];
  }

  /****************************************************
   MPI, Gxyz[Gc_AN][17-19]
  ****************************************************/

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    ID = G2ID[Gc_AN];
    MPI_Bcast(&Gxyz[Gc_AN][17], 3, MPI_DOUBLE, ID, mpi_comm_level1);
  }

  /****************************************************
    freeing of arrays:

    double Fx[Matomnum+1];
    double Fy[Matomnum+1];
    doubel Fz[Matomnum+1];
  ****************************************************/

  free(Fx);
  free(Fy);
  free(Fz);

  /* return */

  return EH0;  
}


void Calc_EXC_EH1(double ECE[])
{
  /****************************************************
     EXC = \sum_{\sigma}
           n_{\sigma}(\epsilon_{xc}-\mu_{xc,\sigma})

     EH1 = -1/2\int {n(r)+n^a(r)} \delta V_H dr
  ****************************************************/

  int i,spin,XC_P_switch,max_spin;
  int numS,numR,My_GNum;
  int MN,MN0,MN1,MN2,n1,n2,n3;
  int MNN0,MNN1,nn0,nn1;
  double EXC[2],EH1,sum;
  double My_EXC[2],My_EH1;
  double My_Eef,Eef;
  double My_Eva,Eva;
  double sum_charge,My_charge;
  double *tmp_array0,*tmp_array1;
  double **Den,*ADen,*PDen;
  double *dVHart,*VNA,**Vxc;
  double *Vef,*RefVxc;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  double Stime_atom, Etime_atom;

  /* dipole moment */
  int Gc_AN,Mc_AN,spe;
  double x,y,z,den,charge,cden_BG;
  double E_dpx,E_dpy,E_dpz; 
  double E_dpx_BG,E_dpy_BG,E_dpz_BG; 
  double C_dpx,C_dpy,C_dpz;
  double My_E_dpx_BG,My_E_dpy_BG,My_E_dpz_BG; 
  double My_E_dpx,My_E_dpy,My_E_dpz; 
  double My_C_dpx,My_C_dpy,My_C_dpz;
  double AU2Debye,AbsD;
  double cot,sit,x0,y0,z0,r;
  double rs,re,ts,te,ps,pe;
  double Sp,Dp,St,Dt,Sr,Dr;
  double r1,dx,dy,dz,dx1,dy1,dz1;
  double x1,y1,z1,den0,exc0;
  double theta,phi,sump,sumr,sumt;
  double sumpx,sumrx,sumtx;
  double sumpy,sumry,sumty;
  double sumpz,sumrz,sumtz;
  double gden0,vxc0;
  double *My_sumr,*My_sumrx,*My_sumry,*My_sumrz;
  int ip,it,ir,Cwan,Rn,Hwan,Gh_AN,h_AN;
  char file_DPM[YOUSO10] = ".dpm";
  FILE *fp_DPM;
  char buf[fp_bsize];          /* setvbuf */
  MPI_Status stat;
  MPI_Request request;
  /* for OpenMP */
  int OMPID,Nthrds,Nthrds0,Nprocs,Nloop;
  double *My_Eef_threads;
  double *My_EH1_threads;
  double *My_EXC_threads;
  double *My_E_dpx_threads;
  double *My_E_dpy_threads;
  double *My_E_dpz_threads;
  double *My_E_dpx_BG_threads;
  double *My_E_dpy_BG_threads;
  double *My_E_dpz_BG_threads;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
   set Vxc_Grid
  ****************************************************/

  XC_P_switch = 0;
  Set_XC_Grid( XC_P_switch,XC_switch,
               Density_Grid[0],Density_Grid[1],
               Density_Grid[2],Density_Grid[3],
               Vxc_Grid[0], Vxc_Grid[1],
               Vxc_Grid[2], Vxc_Grid[3] );

  /****************************************************
   set RefVxc_Grid, where the CA-LDA exchange-correlation 
   functional is alway used.
  ****************************************************/

  XC_P_switch = 0;
  Set_XC_Grid( XC_P_switch, 1,
               ADensity_Grid,ADensity_Grid,
               ADensity_Grid,ADensity_Grid,
               RefVxc_Grid, RefVxc_Grid,
               RefVxc_Grid, RefVxc_Grid );

  /****************************************************
   MPI:

   Density_Grid[0]
   Density_Grid[1]
   ADensity_Grid
   PCCDensity_Grid
   dVHart_Grid
   Vxc_Grid[0]
   Vxc_Grid[1]
   Vef
   RefVxc
  ****************************************************/
  
  My_GNum = My_NGrid1_Poisson*Ngrid2*Ngrid3;  

  Den = (double**)malloc(sizeof(double*)*2); 
  Den[0] = (double*)malloc(sizeof(double)*My_GNum); 
  Den[1] = (double*)malloc(sizeof(double)*My_GNum); 
  ADen = (double*)malloc(sizeof(double)*My_GNum); 
  PDen = (double*)malloc(sizeof(double)*My_GNum); 
  dVHart = (double*)malloc(sizeof(double)*My_GNum); 
  VNA = (double*)malloc(sizeof(double)*My_GNum);   
  Vxc    = (double**)malloc(sizeof(double*)*2);
  Vxc[0] = (double*)malloc(sizeof(double)*My_GNum); 
  Vxc[1] = (double*)malloc(sizeof(double)*My_GNum);
  Vef = (double*)malloc(sizeof(double)*My_GNum);
  RefVxc = (double*)malloc(sizeof(double)*My_GNum);

  /* initialize */
  for (MN1=0; MN1<My_GNum; MN1++){
    Den[0][MN1] = 0.0;
    Den[1][MN1] = 0.0;
    ADen[MN1]   = 0.0;
    PDen[MN1]   = 0.0;
    dVHart[MN1] = 0.0;
    VNA[MN1]    = 0.0;
    Vxc[0][MN1] = 0.0;
    Vxc[1][MN1] = 0.0;
    Vef[MN1]    = 0.0;
    RefVxc[MN1] = 0.0;
  }

  /* use their densities using MPI  */
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    /* Isend */
    if (Num_Snd_Grid1[IDS]!=0){

      numS = Num_Snd_Grid1[IDS]*Ngrid2*Ngrid3;        
      tmp_array0 = (double*)malloc(sizeof(double)*numS*10);

      for (i=0; i<Num_Snd_Grid1[IDS]; i++){ 
	n1 = Snd_Grid1[IDS][i];
        nn1 = My_Cell0[n1];

        MN1 = nn1*Ngrid2*Ngrid3;
        MN0 = i*Ngrid2*Ngrid3;
        for (n2=0; n2<Ngrid2; n2++){
          MN2 = n2*Ngrid3;
          for (n3=0; n3<Ngrid3; n3++){
            MN = MN1 + MN2 + n3;
            MNN0 = MN0 + MN2 + n3; 
            tmp_array0[         MNN0] = Density_Grid[0][MN];
            tmp_array0[  numS + MNN0] = Density_Grid[1][MN];
            tmp_array0[2*numS + MNN0] = ADensity_Grid[MN];
            tmp_array0[3*numS + MNN0] = PCCDensity_Grid[MN];
            tmp_array0[4*numS + MNN0] = dVHart_Grid[MN];
            tmp_array0[5*numS + MNN0] = VNA_Grid[MN];
            tmp_array0[6*numS + MNN0] = Vxc_Grid[0][MN];
            tmp_array0[7*numS + MNN0] = Vxc_Grid[1][MN];
            tmp_array0[8*numS + MNN0] = VEF_Grid[MN];
            tmp_array0[9*numS + MNN0] = RefVxc_Grid[MN];
	  }
	}
      }

      MPI_Isend(&tmp_array0[0], numS*10, MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);
    }

    /* Recv */
    if (Num_Rcv_Grid1[IDR]!=0){

      numR = Num_Rcv_Grid1[IDR]*Ngrid2*Ngrid3;        

      tmp_array1 = (double*)malloc(sizeof(double)*numR*10);

      MPI_Recv(&tmp_array1[0], numR*10, MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

      for (i=0; i<Num_Rcv_Grid1[IDR]; i++){ 
	n1 = Rcv_Grid1[IDR][i];
        nn1 = My_Cell0[n1];
        nn0 = n1 - Start_Grid1[myid];
        MN0 = i*Ngrid2*Ngrid3;
        MN1 = nn0*Ngrid2*Ngrid3;
        for (n2=0; n2<Ngrid2; n2++){
          MN2 = n2*Ngrid3;
          for (n3=0; n3<Ngrid3; n3++){
            MNN1 = MN1 + MN2 + n3;
            MNN0 = MN0 + MN2 + n3;
            Den[0][MNN1] = tmp_array1[         MNN0];
            Den[1][MNN1] = tmp_array1[  numR + MNN0];
            ADen[MNN1]   = tmp_array1[2*numR + MNN0];
            PDen[MNN1]   = tmp_array1[3*numR + MNN0];
            dVHart[MNN1] = tmp_array1[4*numR + MNN0];
            VNA[MNN1]    = tmp_array1[5*numR + MNN0];
            Vxc[0][MNN1] = tmp_array1[6*numR + MNN0];
            Vxc[1][MNN1] = tmp_array1[7*numR + MNN0];
            Vef[MNN1]    = tmp_array1[8*numR + MNN0];
            RefVxc[MNN1] = tmp_array1[9*numR + MNN0];
	  }
	}
      }

      free(tmp_array1);
    }

    if (Num_Snd_Grid1[IDS]!=0){
      MPI_Wait(&request,&stat);
      free(tmp_array0);
    }
  }

  /* use own densities */

  for (n1=0; n1<Ngrid1; n1++){

    if (Cell_ID0[n1]==myid){

      nn0 = My_Cell0[n1];
      nn1 = n1 - Start_Grid1[myid];

      if (nn0!=-1 && Start_Grid1[myid]<=n1 && n1<=End_Grid1[myid]){

        MN0 = nn0*Ngrid2*Ngrid3;  
        MN1 = nn1*Ngrid2*Ngrid3;

        for (n2=0; n2<Ngrid2; n2++){
          MN2 = n2*Ngrid3;
          for (n3=0; n3<Ngrid3; n3++){
            MNN0 = MN0 + MN2 + n3;
            MNN1 = MN1 + MN2 + n3;

            Den[0][MNN1] = Density_Grid[0][MNN0];
            Den[1][MNN1] = Density_Grid[1][MNN0];
            ADen[MNN1]   = ADensity_Grid[MNN0];
            PDen[MNN1]   = PCCDensity_Grid[MNN0];
            dVHart[MNN1] = dVHart_Grid[MNN0];
            VNA[MNN1]    = VNA_Grid[MNN0];
            Vxc[0][MNN1] = Vxc_Grid[0][MNN0];
            Vxc[1][MNN1] = Vxc_Grid[1][MNN0];
            Vef[MNN1]    = VEF_Grid[MNN0];
            RefVxc[MNN1] = RefVxc_Grid[MNN0];
          }    
        }
      }
    }
  }

  /****************************************************
    if (ProExpn_VNA==off), Ena is calculated here.
  ****************************************************/

  if (ProExpn_VNA==0){

    sum = 0.0;
    for (MN=0; MN<My_GNum; MN++){
      sum += (Den[0][MN] + Den[1][MN])*VNA[MN];
    }

    My_Eva = sum*GridVol;

    MPI_Barrier(mpi_comm_level1);
    MPI_Allreduce(&My_Eva, &Eva, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

    ECE[3] = Eva;
  }

  /****************************************************
           electric energy by electric field 
  ****************************************************/

  if (E_Field_switch==1){

  /* get Nthrds0 */  
#pragma omp parallel shared(Nthrds0)
    {
      Nthrds0 = omp_get_num_threads();
    }

    /* allocation of array */
    My_Eef_threads = (double*)malloc(sizeof(double)*Nthrds0);
    for (Nloop=0; Nloop<Nthrds0; Nloop++) My_Eef_threads[Nloop] = 0.0;

#pragma omp parallel shared(My_GNum,My_Eef_threads,Den,Vef) private(OMPID,Nthrds,Nprocs,MN)
    {

      /* get info. on OpenMP */ 

      OMPID = omp_get_thread_num();
      Nthrds = omp_get_num_threads();
      Nprocs = omp_get_num_procs();

      for (MN=OMPID*My_GNum/Nthrds; MN<(OMPID+1)*My_GNum/Nthrds; MN++){
	My_Eef_threads[OMPID] += (Den[0][MN] + Den[1][MN])*Vef[MN];
      }
    } /* #pragma omp parallel */ 

    /* sum of My_Eef_threads */
    My_Eef = 0.0;
    for (Nloop=0; Nloop<Nthrds0; Nloop++){
      My_Eef += My_Eef_threads[Nloop];
    }
    My_Eef *= GridVol;

    /* sum of My_Eef */
    MPI_Barrier(mpi_comm_level1);
    MPI_Allreduce(&My_Eef, &Eef, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

    ECE[12] = Eef;

    /* freeing of array */
    free(My_Eef_threads);
  }
  else {
    ECE[12] = 0.0;
  }

  if (F_VEF_flag==0){
    ECE[12] = 0.0;
  }

  /****************************************************
     EH1 = 1/2\int \delta n(r) \delta V_H dr
  ****************************************************/

  /* get Nthrds0 */  
#pragma omp parallel shared(Nthrds0)
  {
    Nthrds0 = omp_get_num_threads();
  }

  /* allocation of array */
  My_EH1_threads = (double*)malloc(sizeof(double)*Nthrds0);
  for (Nloop=0; Nloop<Nthrds0; Nloop++) My_EH1_threads[Nloop] = 0.0;

#pragma omp parallel shared(My_GNum,My_EH1_threads,Den,ADen,dVHart) private(OMPID,Nthrds,Nprocs,MN)
  {

    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    for (MN=OMPID*My_GNum/Nthrds; MN<(OMPID+1)*My_GNum/Nthrds; MN++){
      My_EH1_threads[OMPID] += (Den[0][MN] + Den[1][MN] - 2.0*ADen[MN])*dVHart[MN];
    }

  } /* #pragma omp parallel */ 

  /* sum of My_Eef_threads */
  My_EH1 = 0.0;
  for (Nloop=0; Nloop<Nthrds0; Nloop++){
    My_EH1 += My_EH1_threads[Nloop];
  }
  My_EH1 *= (0.5*GridVol);

  /* freeing of array */
  free(My_EH1_threads);

  /*
  My_charge = 0.0;
  for (MN=0; MN<My_GNum; MN++){
    My_charge += Den[0][MN] + Den[1][MN];
  }
  My_charge = My_charge*GridVol;
  MPI_Allreduce(&My_charge, &sum_charge, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
  printf("sum_charge=%15.12f\n",sum_charge);
  */

  /****************************************************
   EXC = \sum_{\sigma} n_{\sigma}\epsilon_{xc}
         - n_{atom}\epsilon_{xc}(n_{atom})

   calculation of the difference between the xc energies 
   calculated by wave-function-charge and atomic charge.
   on the coarse grid.    
  ****************************************************/

  if      (SpinP_switch==0) max_spin = 0;
  else if (SpinP_switch==1) max_spin = 1;
  else if (SpinP_switch==3) max_spin = 1;

  if (PCC_switch==0){

    /* get Nthrds0 */  
#pragma omp parallel shared(Nthrds0)
    {
      Nthrds0 = omp_get_num_threads();
    }

    /* allocation of array */
    My_EXC_threads = (double*)malloc(sizeof(double)*Nthrds0);

    for (spin=0; spin<=max_spin; spin++){

      for (Nloop=0; Nloop<Nthrds0; Nloop++) My_EXC_threads[Nloop] = 0.0;

#pragma omp parallel shared(spin,My_GNum,My_EXC_threads,Den,Vxc,ADen,RefVxc) private(OMPID,Nthrds,Nprocs,MN)
      {

	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	for (MN=OMPID*My_GNum/Nthrds; MN<(OMPID+1)*My_GNum/Nthrds; MN++){
	    My_EXC_threads[OMPID] += Den[spin][MN]*Vxc[spin][MN]
                                    -ADen[MN]*RefVxc[MN];
	}
      } /* #pragma omp parallel */

      /* sum of My_EXC_threads */
      My_EXC[spin] = 0.0;
      for (Nloop=0; Nloop<Nthrds0; Nloop++){
	My_EXC[spin] += My_EXC_threads[Nloop];
      }
      My_EXC[spin] *= GridVol;

    } /* spin */

    /* allocation of array */
    free(My_EXC_threads);

  } /* if (PCC_switch==0) */

  else if (PCC_switch==1){

    /* get Nthrds0 */  
#pragma omp parallel shared(Nthrds0)
    {
      Nthrds0 = omp_get_num_threads();
    }

    /* allocation of array */
    My_EXC_threads = (double*)malloc(sizeof(double)*Nthrds0);

    for (spin=0; spin<=max_spin; spin++){

      for (Nloop=0; Nloop<Nthrds0; Nloop++) My_EXC_threads[Nloop] = 0.0;

#pragma omp parallel shared(spin,My_GNum,My_EXC_threads,Den,PDen,Vxc,ADen,RefVxc) private(OMPID,Nthrds,Nprocs,MN)
      {

	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	for (MN=OMPID*My_GNum/Nthrds; MN<(OMPID+1)*My_GNum/Nthrds; MN++){
            My_EXC_threads[OMPID] += (Den[spin][MN] + PDen[MN])*Vxc[spin][MN]
                                    -(ADen[MN] + PDen[MN])*RefVxc[MN];
        }
      } /* #pragma omp parallel */

      /* sum of My_EXC_threads */
      My_EXC[spin] = 0.0;
      for (Nloop=0; Nloop<Nthrds0; Nloop++){
	My_EXC[spin] += My_EXC_threads[Nloop];
      }
      My_EXC[spin] *= GridVol;

    } /* spin */

    /* freeing of array */
    free(My_EXC_threads);
  } /* else if (PCC_switch==1) */

  /****************************************************
    calculation of Exc^(0) and its contribution 
    to forces on the fine mesh
  ****************************************************/

  /* get Nthrds0 */  
#pragma omp parallel shared(Nthrds0)
  {
    Nthrds0 = omp_get_num_threads();
  }

  /* allocation of arrays */
  My_sumr  = (double*)malloc(sizeof(double)*Nthrds0);
  My_sumrx = (double*)malloc(sizeof(double)*Nthrds0);
  My_sumry = (double*)malloc(sizeof(double)*Nthrds0);
  My_sumrz = (double*)malloc(sizeof(double)*Nthrds0);

  /* start calc. */

  rs = 0.0;

  ts = 0.0;
  te = PI;
  St = te + ts;
  Dt = te - ts;

  ps = 0.0;
  pe = 2.0*PI;
  Sp = pe + ps;
  Dp = pe - ps;

  sum  = 0.0;

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];
    Cwan = WhatSpecies[Gc_AN];
    re = Spe_Atom_Cut1[Cwan];
    Sr = re + rs;
    Dr = re - rs;

    for (Nloop=0; Nloop<Nthrds0; Nloop++){
      My_sumr[Nloop]  = 0.0;
      My_sumrx[Nloop] = 0.0;
      My_sumry[Nloop] = 0.0;
      My_sumrz[Nloop] = 0.0;
    }

#pragma omp parallel shared(Dr,Sr,CoarseGL_Abscissae,CoarseGL_Weight,Dt,St,Exc0_GL_Abscissae1,Dp,Sp,Exc0_GL_Abscissae2,Gxyz,Gc_AN,FNAN,natn,ncn,WhatSpecies,atv,F_Vxc_flag,Cwan,Exc0_GL_Weight2,Exc0_GL_Weight1) private(OMPID,Nthrds,Nprocs,ir,r,sumt,sumtx,sumty,sumtz,it,theta,sit,cot,sump,sumpx,sumpy,sumpz,ip,phi,x0,y0,z0,h_AN,Gh_AN,Hwan,x1,y1,z1,dx,dy,dz,r1,den,den0,gden0,dx1,dy1,dz1,exc0,vxc0)
    {

      /* get info. on OpenMP */ 

      OMPID = omp_get_thread_num();
      Nthrds = omp_get_num_threads();
      Nprocs = omp_get_num_procs();

      for (ir=(OMPID*CoarseGL_Mesh/Nthrds); ir<((OMPID+1)*CoarseGL_Mesh/Nthrds); ir++){

	r = 0.50*(Dr*CoarseGL_Abscissae[ir] + Sr);
	sumt  = 0.0; 
	sumtx = 0.0; 
	sumty = 0.0; 
	sumtz = 0.0; 

	for (it=0; it<Exc0_GL_Mesh1; it++){

	  theta = 0.50*(Dt*Exc0_GL_Abscissae1[it] + St);
	  sit = sin(theta); 
	  cot = cos(theta);
	  sump  = 0.0;
	  sumpx = 0.0;
	  sumpy = 0.0;
	  sumpz = 0.0;

	  for (ip=0; ip<Exc0_GL_Mesh2; ip++){

	    phi = 0.50*(Dp*Exc0_GL_Abscissae2[ip] + Sp);
	    x0 = r*sit*cos(phi) + Gxyz[Gc_AN][1];
	    y0 = r*sit*sin(phi) + Gxyz[Gc_AN][2];
	    z0 = r*cot          + Gxyz[Gc_AN][3];

	    /* calculate rho_atom + rho_pcc */ 

	    den = 0.0;
	    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	      Gh_AN = natn[Gc_AN][h_AN];
	      Rn = ncn[Gc_AN][h_AN]; 
	      Hwan = WhatSpecies[Gh_AN];

	      x1 = Gxyz[Gh_AN][1] + atv[Rn][1];
	      y1 = Gxyz[Gh_AN][2] + atv[Rn][2];
	      z1 = Gxyz[Gh_AN][3] + atv[Rn][3];
            
	      dx = x1 - x0;
	      dy = y1 - y0;
	      dz = z1 - z0;
	      r1 = sqrt(dx*dx + dy*dy + dz*dz);

	      /* calculate density */

	      den += ( AtomicDenF(Hwan,r1) + (double)PCC_switch*AtomicPCCF(Hwan,r1) )*F_Vxc_flag;
	      if (h_AN==0) den0 = den;

	      /* calculate gradient of density */

	      if (h_AN==0){
		gden0 = ( Dr_AtomicDenF(Cwan,r1) + (double)PCC_switch*Dr_AtomicPCCF(Cwan,r1) )*F_Vxc_flag;

		dx1 = gden0/r1*dx;
		dy1 = gden0/r1*dy;
		dz1 = gden0/r1*dz;
	      }
	    }

	    /* calculate the CA-LDA exchange-correlation energy density */
	    exc0 = XC_Ceperly_Alder(den,0);

	    /* calculate the CA-LDA exchange-correlation potential */
	    vxc0 = XC_Ceperly_Alder(den,1);

	    /* phi for Gauss-Legendre quadrature */
	    sump  += Exc0_GL_Weight2[ip]*den0*exc0;
	    sumpx += Exc0_GL_Weight2[ip]*vxc0*dx1;
	    sumpy += Exc0_GL_Weight2[ip]*vxc0*dy1;
	    sumpz += Exc0_GL_Weight2[ip]*vxc0*dz1;

	  } /* ip */

	  /* theta for Gauss-Legendre quadrature */
	  sumt  += 0.5*Dp*sump*sit*Exc0_GL_Weight1[it];
	  sumtx += 0.5*Dp*sumpx*sit*Exc0_GL_Weight1[it];
	  sumty += 0.5*Dp*sumpy*sit*Exc0_GL_Weight1[it];
	  sumtz += 0.5*Dp*sumpz*sit*Exc0_GL_Weight1[it];

	} /* it */

	/* r for Gauss-Legendre quadrature */
	My_sumr[OMPID]  += 0.5*Dt*sumt*r*r*CoarseGL_Weight[ir];
	My_sumrx[OMPID] += 0.5*Dt*sumtx*r*r*CoarseGL_Weight[ir];
	My_sumry[OMPID] += 0.5*Dt*sumty*r*r*CoarseGL_Weight[ir];
	My_sumrz[OMPID] += 0.5*Dt*sumtz*r*r*CoarseGL_Weight[ir];

      } /* ir */
    } /* #pragma omp */

    sumr = 0.0;
    for (Nloop=0; Nloop<Nthrds0; Nloop++){
      sumr += My_sumr[Nloop];
    }

    sum += 0.5*Dr*sumr;

    /* add force */

    sumrx = 0.0;
    sumry = 0.0;
    sumrz = 0.0;
    for (Nloop=0; Nloop<Nthrds0; Nloop++){
      sumrx += My_sumrx[Nloop];
      sumry += My_sumry[Nloop];
      sumrz += My_sumrz[Nloop];
    }

    if (2<=level_stdout){
      printf("<Total_Ene>  force(8) myid=%2d  Mc_AN=%2d Gc_AN=%2d  %15.12f %15.12f %15.12f\n",
              myid,Mc_AN,Gc_AN,
              0.5*Dr*sumrx, 0.5*Dr*sumry, 0.5*Dr*sumrz);fflush(stdout);
    }

    Gxyz[Gc_AN][17] += 0.5*Dr*sumrx;
    Gxyz[Gc_AN][18] += 0.5*Dr*sumry;
    Gxyz[Gc_AN][19] += 0.5*Dr*sumrz;

  } /* Mc_AN */

  /* add Exc^0 calculated on the fine mesh to My_EXC */

  My_EXC[0] += 0.5*sum;
  My_EXC[1] += 0.5*sum;

  /* freeing of arrays */
  free(My_sumr);
  free(My_sumrx);
  free(My_sumry);
  free(My_sumrz);

  /****************************************************
   MPI, Gxyz[Gc_AN][17-19]
  ****************************************************/

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    ID = G2ID[Gc_AN];
    MPI_Bcast(&Gxyz[Gc_AN][17], 3, MPI_DOUBLE, ID, mpi_comm_level1);

    if (2<=level_stdout && myid==Host_ID){
      printf("<Total_Ene>  force(t) myid=%2d Gc_AN=%2d  %15.12f %15.12f %15.12f\n",
              myid,Gc_AN,Gxyz[Gc_AN][17],Gxyz[Gc_AN][18],Gxyz[Gc_AN][19]);fflush(stdout);
    }
  }

  /****************************************************
   MPI:

   EH1, EXC
  ****************************************************/

  MPI_Barrier(mpi_comm_level1);
  MPI_Allreduce(&My_EH1, &EH1, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
  for (spin=0; spin<=max_spin; spin++){
    MPI_Allreduce(&My_EXC[spin], &EXC[spin], 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
  }

  if (SpinP_switch==0){
    ECE[5] = EH1;
    ECE[6] = EXC[0];
    ECE[7] = EXC[0];
  }
  else if (SpinP_switch==1 || SpinP_switch==3) {
    ECE[5] = EH1;
    ECE[6] = EXC[0];
    ECE[7] = EXC[1];
  }

  if (F_dVHart_flag==0){
    ECE[5] = 0.0;
  }

  if (F_Vxc_flag==0){
    ECE[6] = 0.0;
    ECE[7] = 0.0;
  }

  /****************************************************
             calculation of dipole moment
  ****************************************************/

  /* contribution from electron density */

  /* get Nthrds0 */  
#pragma omp parallel shared(Nthrds0)
  {
    Nthrds0 = omp_get_num_threads();
  }

  /* allocation of array */
  My_E_dpx_threads = (double*)malloc(sizeof(double)*Nthrds0);
  My_E_dpy_threads = (double*)malloc(sizeof(double)*Nthrds0);
  My_E_dpz_threads = (double*)malloc(sizeof(double)*Nthrds0);
  My_E_dpx_BG_threads = (double*)malloc(sizeof(double)*Nthrds0);
  My_E_dpy_BG_threads = (double*)malloc(sizeof(double)*Nthrds0);
  My_E_dpz_BG_threads = (double*)malloc(sizeof(double)*Nthrds0);

  for (Nloop=0; Nloop<Nthrds0; Nloop++){
    My_E_dpx_threads[Nloop] = 0.0;
    My_E_dpy_threads[Nloop] = 0.0; 
    My_E_dpz_threads[Nloop] = 0.0; 
    My_E_dpx_BG_threads[Nloop] = 0.0; 
    My_E_dpy_BG_threads[Nloop] = 0.0; 
    My_E_dpz_BG_threads[Nloop] = 0.0;
  }

#pragma omp parallel shared(My_E_dpz_BG_threads,My_E_dpy_BG_threads,My_E_dpx_BG_threads,My_E_dpz_threads,My_E_dpy_threads,My_E_dpx_threads,Den,Grid_Origin,gtv,myid,Start_Grid1,Ngrid3,Ngrid2,My_GNum) private(OMPID,Nthrds,Nprocs,MN,n1,n2,n3,x,y,z,den)
  {

    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    for (MN=OMPID*My_GNum/Nthrds; MN<(OMPID+1)*My_GNum/Nthrds; MN++){
    
      n1 = MN/(Ngrid2*Ngrid3);
      n2 = (MN - n1*(Ngrid2*Ngrid3))/Ngrid3;
      n3 = MN - n1*(Ngrid2*Ngrid3) - n2*Ngrid3;
      n1 = n1 + Start_Grid1[myid];
    
      x = (double)n1*gtv[1][1] + (double)n2*gtv[2][1]
	+ (double)n3*gtv[3][1] + Grid_Origin[1];
      y = (double)n1*gtv[1][2] + (double)n2*gtv[2][2]
	+ (double)n3*gtv[3][2] + Grid_Origin[2];
      z = (double)n1*gtv[1][3] + (double)n2*gtv[2][3]
	+ (double)n3*gtv[3][3] + Grid_Origin[3];
    
      den = Den[0][MN] + Den[1][MN];

      My_E_dpx_threads[OMPID] += den*x;
      My_E_dpy_threads[OMPID] += den*y;
      My_E_dpz_threads[OMPID] += den*z; 

      My_E_dpx_BG_threads[OMPID] += x;
      My_E_dpy_BG_threads[OMPID] += y;
      My_E_dpz_BG_threads[OMPID] += z; 

    } /* MN */
  } /* #pragma omp parallel */

  /* sum of My_EXC_threads */

  My_E_dpx = 0.0;
  My_E_dpy = 0.0;
  My_E_dpz = 0.0;
  My_E_dpx_BG = 0.0;
  My_E_dpy_BG = 0.0;
  My_E_dpz_BG = 0.0;

  for (Nloop=0; Nloop<Nthrds0; Nloop++){
    My_E_dpx    += My_E_dpx_threads[Nloop];
    My_E_dpy    += My_E_dpy_threads[Nloop];
    My_E_dpz    += My_E_dpz_threads[Nloop];
    My_E_dpx_BG += My_E_dpx_BG_threads[Nloop];
    My_E_dpy_BG += My_E_dpy_BG_threads[Nloop];
    My_E_dpz_BG += My_E_dpz_BG_threads[Nloop];
  }

  /* freeing of array */
  free(My_E_dpz_BG_threads);
  free(My_E_dpy_BG_threads);
  free(My_E_dpx_BG_threads);
  free(My_E_dpz_threads);
  free(My_E_dpy_threads);
  free(My_E_dpx_threads);

  MPI_Allreduce(&My_E_dpx, &E_dpx, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
  MPI_Allreduce(&My_E_dpy, &E_dpy, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
  MPI_Allreduce(&My_E_dpz, &E_dpz, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

  MPI_Allreduce(&My_E_dpx_BG, &E_dpx_BG, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
  MPI_Allreduce(&My_E_dpy_BG, &E_dpy_BG, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
  MPI_Allreduce(&My_E_dpz_BG, &E_dpz_BG, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

  E_dpx = E_dpx*GridVol;
  E_dpy = E_dpy*GridVol;
  E_dpz = E_dpz*GridVol;

  cden_BG = system_charge/Cell_Volume; 

  E_dpx_BG = E_dpx_BG*GridVol*cden_BG;
  E_dpy_BG = E_dpy_BG*GridVol*cden_BG;
  E_dpz_BG = E_dpz_BG*GridVol*cden_BG;

  /* contribution from core charge */

  My_C_dpx = 0.0;
  My_C_dpy = 0.0;
  My_C_dpz = 0.0;

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    x = Gxyz[Gc_AN][1];
    y = Gxyz[Gc_AN][2];
    z = Gxyz[Gc_AN][3];

    spe = WhatSpecies[Gc_AN];
    charge = Spe_Core_Charge[spe];
    My_C_dpx += charge*x;
    My_C_dpy += charge*y;
    My_C_dpz += charge*z;
  }

  MPI_Allreduce(&My_C_dpx, &C_dpx, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
  MPI_Allreduce(&My_C_dpy, &C_dpy, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
  MPI_Allreduce(&My_C_dpz, &C_dpz, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

  AU2Debye = 2.54174776;

  dipole_moment[0][1] = AU2Debye*(C_dpx - E_dpx - E_dpx_BG);
  dipole_moment[0][2] = AU2Debye*(C_dpy - E_dpy - E_dpy_BG);
  dipole_moment[0][3] = AU2Debye*(C_dpz - E_dpz - E_dpz_BG);

  dipole_moment[1][1] = AU2Debye*C_dpx;
  dipole_moment[1][2] = AU2Debye*C_dpy;
  dipole_moment[1][3] = AU2Debye*C_dpz;

  dipole_moment[2][1] = -AU2Debye*E_dpx;
  dipole_moment[2][2] = -AU2Debye*E_dpy;
  dipole_moment[2][3] = -AU2Debye*E_dpz;

  dipole_moment[3][1] = -AU2Debye*E_dpx_BG;
  dipole_moment[3][2] = -AU2Debye*E_dpy_BG;
  dipole_moment[3][3] = -AU2Debye*E_dpz_BG;

  AbsD = sqrt( dipole_moment[0][1]*dipole_moment[0][1]
             + dipole_moment[0][2]*dipole_moment[0][2]
             + dipole_moment[0][3]*dipole_moment[0][3] );

  if (myid==Host_ID){

    if (0<level_stdout){
      printf("\n*******************************************************\n"); fflush(stdout);
      printf("                  Dipole moment (Debye)                 \n");  fflush(stdout);
      printf("*******************************************************\n\n"); fflush(stdout);

      printf(" Absolute D %17.8f\n\n",AbsD);
      printf("                      Dx                Dy                Dz\n"); fflush(stdout);
      printf(" Total       %17.8f %17.8f %17.8f\n",
	     dipole_moment[0][1],dipole_moment[0][2],dipole_moment[0][3]);fflush(stdout);
      printf(" Core        %17.8f %17.8f %17.8f\n",
	     dipole_moment[1][1],dipole_moment[1][2],dipole_moment[1][3]);fflush(stdout);
      printf(" Electron    %17.8f %17.8f %17.8f\n",
	     dipole_moment[2][1],dipole_moment[2][2],dipole_moment[2][3]);fflush(stdout);
      printf(" Back ground %17.8f %17.8f %17.8f\n",
	     dipole_moment[3][1],dipole_moment[3][2],dipole_moment[3][3]);fflush(stdout);
    }

    /********************************************************
             write the dipole moments to a file
    ********************************************************/

    fnjoint(filepath,filename,file_DPM);

    if ((fp_DPM = fopen(file_DPM,"w")) != NULL){

#ifdef xt3
      setvbuf(fp_DPM,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      fprintf(fp_DPM,"\n");
      fprintf(fp_DPM,"***********************************************************\n");
      fprintf(fp_DPM,"***********************************************************\n");
      fprintf(fp_DPM,"                    Dipole moment (Debye)                  \n");
      fprintf(fp_DPM,"***********************************************************\n");
      fprintf(fp_DPM,"***********************************************************\n\n");

      fprintf(fp_DPM," Absolute D %17.8f\n\n",AbsD);
      fprintf(fp_DPM,"                      Dx                Dy                Dz\n");
      fprintf(fp_DPM," Total       %17.8f %17.8f %17.8f\n",
                dipole_moment[0][1],dipole_moment[0][2],dipole_moment[0][3]);
      fprintf(fp_DPM," Core        %17.8f %17.8f %17.8f\n",
                dipole_moment[1][1],dipole_moment[1][2],dipole_moment[1][3]);
      fprintf(fp_DPM," Electron    %17.8f %17.8f %17.8f\n",
                dipole_moment[2][1],dipole_moment[2][2],dipole_moment[2][3]);
      fprintf(fp_DPM," Back ground %17.8f %17.8f %17.8f\n",
                dipole_moment[3][1],dipole_moment[3][2],dipole_moment[3][3]);

      fclose(fp_DPM);
    }
    else{
      printf("Failure of saving the DPM file.\n");fflush(stdout);
    }
  }

  /****************************************************
   freeing of arrays:

    free(Den[0]);
    free(Den[1]);
    free(Den);
    free(ADen);
    free(PDen);
    free(dVHart);
    free(VNA);
    free(Vxc[0]);
    free(Vxc[1]);
    free(Vxc);
    free(Vef);
  ****************************************************/

  free(RefVxc);
  free(Vef);
  free(Vxc[1]);
  free(Vxc[0]);
  free(Vxc);
  free(VNA);
  free(dVHart);
  free(PDen);
  free(ADen);
  free(Den[1]);
  free(Den[0]);
  free(Den);

}



void EH0_TwoCenter(int Gc_AN, int h_AN, double VH0ij[4])
{ 
  int n1,ban;
  int Gh_AN,Rn,wan1,wan2;
  double dv,x,y,z,r,r2,va0,rho0,dr_va0;
  double z2,sum,sumr,sumx,sumy,sumz,wt;

  Gh_AN = natn[Gc_AN][h_AN];
  Rn = ncn[Gc_AN][h_AN];
  wan1 = WhatSpecies[Gc_AN];
  ban = Spe_Spe2Ban[wan1];
  wan2 = WhatSpecies[Gh_AN];
  dv = dv_EH0[ban];
  
  sum = 0.0;
  sumr = 0.0;

  for (n1=0; n1<TGN_EH0[ban]; n1++){
    x = GridX_EH0[ban][n1];
    y = GridY_EH0[ban][n1];
    z = GridZ_EH0[ban][n1];
    rho0 = Arho_EH0[ban][n1];
    wt = Wt_EH0[ban][n1];
    z2 = z - Dis[Gc_AN][h_AN];
    r2 = sqrt(x*x + y*y + z2*z2);

    /* for empty atoms or finite elemens basis */
    if (r2<1.0e-10) r2 = 1.0e-10;

    va0 = VH_AtomF(wan2,r2);
    sum = sum + wt*va0*rho0;

    if (h_AN!=0 && 1.0e-14<r2){
      dr_va0 = Dr_VH_AtomF(wan2,r2);
      sumr = sumr - wt*rho0*dr_va0*z2/r2;
    }
  }

  sum  = sum*dv;

  if (h_AN!=0){

    /* for empty atoms or finite elemens basis */
    r = Dis[Gc_AN][h_AN];
    if (r<1.0e-10) r = 1.0e-10;

    x = Gxyz[Gc_AN][1] - (Gxyz[Gh_AN][1] + atv[Rn][1]);
    y = Gxyz[Gc_AN][2] - (Gxyz[Gh_AN][2] + atv[Rn][2]);
    z = Gxyz[Gc_AN][3] - (Gxyz[Gh_AN][3] + atv[Rn][3]);
    sumr = sumr*dv;
    sumx = sumr*x/r;
    sumy = sumr*y/r;
    sumz = sumr*z/r;
  }
  else{
    sumx = 0.0;
    sumy = 0.0;
    sumz = 0.0;
  }

  VH0ij[0] = sum;
  VH0ij[1] = sumx;
  VH0ij[2] = sumy;
  VH0ij[3] = sumz;
}














void EH0_TwoCenter_at_Cutoff(int wan1, int wan2, double VH0ij[4])
{ 
  int n1,ban;
  double dv,x,y,z,r1,r2,va0,rho0,dr_va0,rcut;
  double z2,sum,sumr,sumx,sumy,sumz,wt;
  /* for OpenMP */
  int OMPID,Nthrds,Nthrds0,Nprocs,Nloop;
  double *my_sum_threads;

  ban = Spe_Spe2Ban[wan1];
  dv  = dv_EH0[ban];

  rcut = Spe_Atom_Cut1[wan1] + Spe_Atom_Cut1[wan2];

  /* get Nthrds0 */  
#pragma omp parallel shared(Nthrds0)
  {
    Nthrds0 = omp_get_num_threads();
  }

  /* allocation of array */
  my_sum_threads = (double*)malloc(sizeof(double)*Nthrds0);

  for (Nloop=0; Nloop<Nthrds0; Nloop++){
    my_sum_threads[Nloop] = 0.0;
  }

#pragma omp parallel shared(wan2,Wt_EH0,my_sum_threads,rcut,Arho_EH0,GridZ_EH0,GridY_EH0,GridX_EH0,TGN_EH0,ban) private(n1,OMPID,Nthrds,Nprocs,x,y,z,rho0,wt,z2,r2,va0)
  {
    /* get info. on OpenMP */

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    for (n1=OMPID*TGN_EH0[ban]/Nthrds; n1<(OMPID+1)*TGN_EH0[ban]/Nthrds; n1++){

      x = GridX_EH0[ban][n1];
      y = GridY_EH0[ban][n1];
      z = GridZ_EH0[ban][n1];
      rho0 = Arho_EH0[ban][n1];
      wt = Wt_EH0[ban][n1];
      z2 = z - rcut;
      r2 = sqrt(x*x + y*y + z2*z2);
      va0 = VH_AtomF(wan2,r2);

      my_sum_threads[OMPID] += wt*va0*rho0;
    }

  } /* #pragma omp parallel */

  sum  = 0.0;
  for (Nloop=0; Nloop<Nthrds0; Nloop++){
    sum += my_sum_threads[Nloop];
  }

  sum  = sum*dv;
  sumx = 0.0;
  sumy = 0.0;
  sumz = 0.0;

  VH0ij[0] = sum;
  VH0ij[1] = sumx;
  VH0ij[2] = sumy;
  VH0ij[3] = sumz;

  /* freeing of array */
  free(my_sum_threads);
}







double Calc_Ehub()
{   
 /****************************************************
         LDA+U energy correction added by MJ
  ****************************************************/

  int Mc_AN,Gc_AN,wan1;
  int cnt1,cnt2,l1,mul1,m1,l2,mul2,m2;
  int spin,max_spin;
  double My_Ehub,Ehub,Uvalue;
  int numprocs,myid,ID;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

 /****************************************************
                 caculation of My_Ehub
  ****************************************************/

  if      (SpinP_switch==0) max_spin = 0;
  else if (SpinP_switch==1) max_spin = 1;
  else if (SpinP_switch==3) max_spin = 1;

  My_Ehub = 0.0;

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    wan1 = WhatSpecies[Gc_AN];

   /****************************************************
                     collinear case
    ****************************************************/

    if (SpinP_switch!=3){

      for (spin=0; spin<=max_spin; spin++){


        /* Hubbard term, 0.5*Tr(N) */

	cnt1 = 0;
	for(l1=0; l1<=Spe_MaxL_Basis[wan1]; l1++ ){
	  for(mul1=0; mul1<Spe_Num_Basis[wan1][l1]; mul1++){

	    Uvalue = Hub_U_Basis[wan1][l1][mul1];
	    for(m1=0; m1<(2*l1+1); m1++){

	      My_Ehub += 0.5*Uvalue*DM_onsite[0][spin][Mc_AN][cnt1][cnt1];

	      cnt1++;
	    }
	  }
	}


        /* Hubbard term, -0.5*Tr(N*N) */

	cnt1 = 0;
	for(l1=0; l1<=Spe_MaxL_Basis[wan1]; l1++ ){
	  for(mul1=0; mul1<Spe_Num_Basis[wan1][l1]; mul1++){
	    for(m1=0; m1<(2*l1+1); m1++){

	      cnt2 = 0;
	      for(l2=0; l2<=Spe_MaxL_Basis[wan1]; l2++ ){
		for(mul2=0; mul2<Spe_Num_Basis[wan1][l2]; mul2++){
		  for(m2=0; m2<(2*l2+1); m2++){

		    if (l1==l2 && mul1==mul2){
		      Uvalue = Hub_U_Basis[wan1][l1][mul1];
		      My_Ehub -= 0.5*Uvalue*DM_onsite[0][spin][Mc_AN][cnt1][cnt2]*
			                    DM_onsite[0][spin][Mc_AN][cnt2][cnt1];
		    }

		    cnt2++;
		  }
		}
	      }

	      cnt1++;
	    }
	  }
	}
      }
    }

   /****************************************************
                     non-collinear case
    ****************************************************/

    else {

      /* Hubbard term, 0.5*Tr(N) */

      cnt1 = 0;
      for(l1=0; l1<=Spe_MaxL_Basis[wan1]; l1++ ){
	for(mul1=0; mul1<Spe_Num_Basis[wan1][l1]; mul1++){

	  Uvalue = Hub_U_Basis[wan1][l1][mul1];
	  for(m1=0; m1<(2*l1+1); m1++){

	    My_Ehub += 0.5*Uvalue*( NC_OcpN[0][0][0][Mc_AN][cnt1][cnt1].r
				  + NC_OcpN[0][1][1][Mc_AN][cnt1][cnt1].r);

	    cnt1++;
	  }
	}
      }

      /* Hubbard term, -0.5*Tr(N*N) */

      cnt1 = 0;
      for(l1=0; l1<=Spe_MaxL_Basis[wan1]; l1++ ){
	for(mul1=0; mul1<Spe_Num_Basis[wan1][l1]; mul1++){
	  for(m1=0; m1<(2*l1+1); m1++){

	    cnt2 = 0;
	    for(l2=0; l2<=Spe_MaxL_Basis[wan1]; l2++ ){
	      for(mul2=0; mul2<Spe_Num_Basis[wan1][l2]; mul2++){
		for(m2=0; m2<(2*l2+1); m2++){

		  if (l1==l2 && mul1==mul2){

		    Uvalue = Hub_U_Basis[wan1][l1][mul1];

		    My_Ehub -= 0.5*Uvalue*( NC_OcpN[0][0][0][Mc_AN][cnt1][cnt2].r*
					    NC_OcpN[0][0][0][Mc_AN][cnt1][cnt2].r
					    +
					    NC_OcpN[0][0][0][Mc_AN][cnt1][cnt2].i*
					    NC_OcpN[0][0][0][Mc_AN][cnt1][cnt2].i
					    +
					    NC_OcpN[0][0][1][Mc_AN][cnt1][cnt2].r*
					    NC_OcpN[0][0][1][Mc_AN][cnt1][cnt2].r
					    +
					    NC_OcpN[0][0][1][Mc_AN][cnt1][cnt2].i*
					    NC_OcpN[0][0][1][Mc_AN][cnt1][cnt2].i
					    +
					    NC_OcpN[0][1][0][Mc_AN][cnt1][cnt2].r*
					    NC_OcpN[0][1][0][Mc_AN][cnt1][cnt2].r
					    +
					    NC_OcpN[0][1][0][Mc_AN][cnt1][cnt2].i*
					    NC_OcpN[0][1][0][Mc_AN][cnt1][cnt2].i
					    +
					    NC_OcpN[0][1][1][Mc_AN][cnt1][cnt2].r*
					    NC_OcpN[0][1][1][Mc_AN][cnt1][cnt2].r
					    +
					    NC_OcpN[0][1][1][Mc_AN][cnt1][cnt2].i*
					    NC_OcpN[0][1][1][Mc_AN][cnt1][cnt2].i );

		  }

		  cnt2++;
		}
	      }
	    }

	    cnt1++;
	  }
	}
      }
    }
  }

  if (SpinP_switch==0) My_Ehub = 2.0*My_Ehub;

 /****************************************************
                      MPI My_Ehub
  ****************************************************/

  MPI_Allreduce(&My_Ehub, &Ehub, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

  /* if (F_U_flag==0) */
  if (F_U_flag==0) Ehub = 0.0;
  
  return Ehub;  
}



/* okuno */
double Calc_EdftD()
{
  /*********************************************
  The subroutine calculates the semiemprical 
  vdW correction to DFT-GGA proposed by 
  S. Grimme, J. Comput. Chem. 27, 1787 (2006).
  *********************************************/

  double My_EdftD,EdftD;
  double rij[4],fdamp,fdamp2;
  double rij0[4],par;
  double dist,dist6,dist2;
  double exparg,expval;
  double rcut_dftD2;
  int numprocs,myid,ID;
  int Mc_AN,Gc_AN,wanA,wanB;
  int Gc_BN;
  int nrm,nr;
  int i,j;
  int n1,n2,n3;
  int per_flag1,per_flag2;
  int n1_max,n2_max,n3_max; 
  double test_ene;
  double dblcnt_factor;
  double dEx,dEy,dEz,dist7;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  My_EdftD = 0.0;
  EdftD    = 0.0;
  rcut_dftD2 = rcut_dftD*rcut_dftD;

  dblcnt_factor=0.5;

  /* here we calculate DFT-D diserpesion energy */
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    dEx = 0.0;
    dEy = 0.0;
    dEz = 0.0;

    Gc_AN = M2G[Mc_AN];
    wanA = WhatSpecies[Gc_AN];
    per_flag1 = (int)Gxyz[Gc_AN][60]; 

    for(Gc_BN=1; Gc_BN<=atomnum; Gc_BN++){

      wanB = WhatSpecies[Gc_BN];
      per_flag2 = (int)Gxyz[Gc_BN][60]; 

      rij0[1] = Gxyz[Gc_AN][1] - Gxyz[Gc_BN][1];
      rij0[2] = Gxyz[Gc_AN][2] - Gxyz[Gc_BN][2];
      rij0[3] = Gxyz[Gc_AN][3] - Gxyz[Gc_BN][3];

      par = beta_dftD/(Rsum_dftD[wanA][wanB]);

      if (per_flag1==0 && per_flag2==0){
        n1_max = 0;
        n2_max = 0;
        n3_max = 0;
      }
      else if (per_flag1==0 && per_flag2==1){
        n1_max = n1_DFT_D;
        n2_max = n2_DFT_D;
        n3_max = n3_DFT_D;
      }
      else if (per_flag1==1 && per_flag2==0){
        n1_max = 0;
        n2_max = 0;
        n3_max = 0;
      }
      else if (per_flag1==1 && per_flag2==1){
        n1_max = n1_DFT_D;
        n2_max = n2_DFT_D;
        n3_max = n3_DFT_D;
      }
/*
      printf("Gc_AN=%2d Gc_BN=%2d %2d %2d %2d %2d %2d\n",Gc_AN,Gc_BN,per_flag1,per_flag2,n1_max,n2_max,n3_max);
*/
      for (n1=-n1_max; n1<=n1_max; n1++){
	for (n2=-n2_max; n2<=n2_max; n2++){
	  for (n3=-n3_max; n3<=n3_max; n3++){
            
            /* for double counting */
            if((!(abs(n1)+abs(n2)+abs(n3))==0) && (per_flag1==0 && per_flag2==1) ){
	      dblcnt_factor = 1.0;
	    }
            else{
	      dblcnt_factor = 0.5;
            }

	    rij[1] = rij0[1] - ( (double)n1*tv[1][1]
	                       + (double)n2*tv[2][1] 
	                       + (double)n3*tv[3][1] ); 

	    rij[2] = rij0[2] - ( (double)n1*tv[1][2]
			       + (double)n2*tv[2][2] 
			       + (double)n3*tv[3][2] ); 

	    rij[3] = rij0[3] - ( (double)n1*tv[1][3]
			       + (double)n2*tv[2][3] 
			       + (double)n3*tv[3][3] ); 

            dist2 = rij[1]*rij[1] + rij[2]*rij[2] + rij[3]*rij[3];

            if (0.1<dist2 && dist2<=rcut_dftD2){

	      dist  = sqrt(dist2); 
	      dist6 = dist2*dist2*dist2;            

	      /* calculate the vdW energy */
	      exparg = -beta_dftD*((dist/Rsum_dftD[wanA][wanB])-1.0);
	      expval = exp(exparg);
	      fdamp = scal6_dftD/(1.0+expval);
	      My_EdftD -= dblcnt_factor*C6ij_dftD[wanA][wanB]/dist6*fdamp;

	      /* calculate the gradient of the vdW energy */

              dist7 = dist6 * dist;
	      fdamp2 = C6ij_dftD[wanA][wanB]*fdamp/dist6*(expval*par/(1.0+expval) - 6.0/dist);
              dEx -= fdamp2*rij[1]/dist;
              dEy -= fdamp2*rij[2]/dist;
              dEz -= fdamp2*rij[3]/dist;
	    }

	  } /* n3 */
	} /* n2 */
      } /* n1 */
    } /* Gc_BN */

    Gxyz[Gc_AN][17] += dEx;
    Gxyz[Gc_AN][18] += dEy;
    Gxyz[Gc_AN][19] += dEz;

    /*
    printf("Gc_AN=%2d dEx=%15.12f dEy=%15.12f dEz=%15.12f\n",Gc_AN,dEx,dEy,dEz);
    */

  } /* Mc_AN */

  MPI_Allreduce(&My_EdftD, &EdftD, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

  return EdftD;
}
/* okuno */
