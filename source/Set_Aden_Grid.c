/**********************************************************************
  Set_Aden_Grid.c:

     Set_Aden_Grid.c is a subroutine to calculate a charge density 
     superposed atomic densities on grid.

  Log of Set_Aden_Grid.c:

     22/Nov/2001  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
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
#include <omp.h>
#endif


double Set_Aden_Grid(int init_density)
{
  /****************************************************
          Densities by the atomic superposition
                   densities on grids
  ****************************************************/

  static int firsttime=1;
  int i,k,MN,ct_AN,Gc_AN,Mc_AN,top_num;
  int Rn,Cwan,Nc,GNc,GRc,n,size1,size2;
  int My_Max,Max_Size;
  int size_AtomDen_Grid;
  double time0;
  double x,y,z,DenA,DenPCC,dDenA,dDenPCC;
  double tmp0,dx,dy,dz,r,rmin=10e-14;
  double Nele,Nu,Nd,M,ocupcy_u,ocupcy_d;
  double rho,mag,magx,magy,magz,theta,phi;
  double TStime,TEtime;
  double Cxyz[4];
  double S_coordinate[3];
  double *tmp_array;
  double *tmp_array2;
  double **AtomDen_Grid;
  double **PCCDen_Grid;
  double *AtomDen2_Grid;
  double *PCCDen2_Grid;
  int *Snd_Size,*Rcv_Size;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  double Stime_atom, Etime_atom;
  double Nup,Ndown,sit,cot,sip,cop;
  dcomplex U[2][2];

  MPI_Status stat;
  MPI_Request request;

  /* for OpenMP */
  int OMPID,Nthrds,Nprocs;

  /* MPI */
  if (atomnum<=MYID_MPI_COMM_WORLD) return 0.0;
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  
  dtime(&TStime);

  /****************************************************
    allocation of arrays:

    int Snd_Size[numprocs]
    int Rcv_Size[numprocs]

    double AtomDen_Grid[Matomnum+MatomnumF+1]
                       [GridN_Atom[Gc_AN]]

    double PCCDen_Grid[Matomnum+MatomnumF+1]
                       [GridN_Atom[Gc_AN]]
  ****************************************************/

  Snd_Size = (int*)malloc(sizeof(int)*numprocs); 
  Rcv_Size = (int*)malloc(sizeof(int)*numprocs); 

  size_AtomDen_Grid = Matomnum+MatomnumF+1;
  AtomDen_Grid = (double**)malloc(sizeof(double*)*(Matomnum+MatomnumF+1)); 
  AtomDen_Grid[0] = (double*)malloc(sizeof(double)*1); 
  for (Mc_AN=1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
    Gc_AN = F_M2G[Mc_AN];
    AtomDen_Grid[Mc_AN] = (double*)malloc(sizeof(double)*GridN_Atom[Gc_AN]);
    size_AtomDen_Grid += GridN_Atom[Gc_AN];
  }

  PCCDen_Grid = (double**)malloc(sizeof(double*)*(Matomnum+MatomnumF+1)); 
  PCCDen_Grid[0] = (double*)malloc(sizeof(double)*1); 
  for (Mc_AN=1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
    Gc_AN = F_M2G[Mc_AN];
    PCCDen_Grid[Mc_AN] = (double*)malloc(sizeof(double)*GridN_Atom[Gc_AN]); 
  }

  /* PrintMemory */

  if (firsttime==1){
    PrintMemory("Set_Aden_Grid: AtomDen_Grid", sizeof(double)*size_AtomDen_Grid, NULL);
    PrintMemory("Set_Aden_Grid: PCCDen_Grid",  sizeof(double)*size_AtomDen_Grid, NULL);
    firsttime = 0;
  }

  /******************************************************
                 setting of AtomDen_Grid
  ******************************************************/
 
  /**************************************
   for spin non-collinear
   1. set rho, mx, my, mz
   2. calculate theta and phi 
   3. n_up = (rho+m)/2 
      n_dn = (rho-m)/2 
  **************************************/

  if (SpinP_switch==3){
    for (MN=0; MN<My_NumGrid1; MN++){
      Density_Grid[0][MN] = 0.0;
      Density_Grid[1][MN] = 0.0;
      Density_Grid[2][MN] = 0.0;
      Density_Grid[3][MN] = 0.0;
      ADensity_Grid[MN]   = 0.0;
    }
  }

  /* spin collinear */
  else{ 
    for (MN=0; MN<My_NumGrid1; MN++){
      Density_Grid[0][MN] = 0.0;
      Density_Grid[1][MN] = 0.0;
      ADensity_Grid[MN]   = 0.0;
    }
  }  

  /* PCC */
  if (PCC_switch==1) {
    for (MN=0; MN<My_NumGrid1; MN++){
      PCCDensity_Grid[MN] = 0.0;
    }
  }

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    dtime(&Stime_atom);

    Gc_AN = M2G[Mc_AN];    
    Cwan = WhatSpecies[Gc_AN];
 
#pragma omp parallel shared(PCCDen_Grid,PCC_switch,AtomDen_Grid,Cwan,Gxyz,atv,Mc_AN,CellListAtom,GridListAtom,GridN_Atom,Gc_AN) private(OMPID,Nthrds,Nprocs,Nc,GNc,GRc,Cxyz,dx,dy,dz,r,DenA)
    {

      /* get info. on OpenMP */ 

      OMPID = omp_get_thread_num();
      Nthrds = omp_get_num_threads();
      Nprocs = omp_get_num_procs();

      for (Nc=OMPID*GridN_Atom[Gc_AN]/Nthrds; Nc<(OMPID+1)*GridN_Atom[Gc_AN]/Nthrds; Nc++){
      
	GNc = GridListAtom[Mc_AN][Nc];
	GRc = CellListAtom[Mc_AN][Nc];
      
	Get_Grid_XYZ(GNc,Cxyz);
	dx = Cxyz[1] + atv[GRc][1] - Gxyz[Gc_AN][1];
	dy = Cxyz[2] + atv[GRc][2] - Gxyz[Gc_AN][2];
	dz = Cxyz[3] + atv[GRc][3] - Gxyz[Gc_AN][3];
      
	r = sqrt(dx*dx + dy*dy + dz*dz); 

	/* AtomicDenF */

	DenA = AtomicDenF(Cwan,r);
	AtomDen_Grid[Mc_AN][Nc] = DenA;
      
	/*  partial core correction */
	if (PCC_switch==1) {
	  PCCDen_Grid[Mc_AN][Nc] = AtomicPCCF(Cwan,r);
	}
      }

    } /* #pragma omp parallel */

    dtime(&Etime_atom);
    time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
  }

  /******************************************************
   MPI:
        AtomDen_Grid
  ******************************************************/

  /* find data size for sending and recieving */

  tag = 999;
  My_Max = -10000;
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){
      /*  sending size */
      if (F_Snd_Num[IDS]!=0){
        /* find data size  */
        size1 = 0; 
        for (n=0; n<F_Snd_Num[IDS]; n++){
          Gc_AN = Snd_GAN[IDS][n];
          size1 += GridN_Atom[Gc_AN];
	}

        Snd_Size[IDS] = size1;
        MPI_Isend(&size1, 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
      }
      else{
        Snd_Size[IDS] = 0;
      }

      /*  receiving size */
      if (F_Rcv_Num[IDR]!=0){
        MPI_Recv(&size2, 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
        Rcv_Size[IDR] = size2;
      }
      else{
        Rcv_Size[IDR] = 0;
      }

      if (F_Snd_Num[IDS]!=0) MPI_Wait(&request,&stat);

    } 
    else{
      Snd_Size[IDS] = 0;
      Rcv_Size[IDR] = 0;
    }

    if (My_Max<Snd_Size[IDS]) My_Max = Snd_Size[IDS];
    if (My_Max<Rcv_Size[IDR]) My_Max = Rcv_Size[IDR];
  }  

  MPI_Allreduce(&My_Max, &Max_Size, 1, MPI_INT, MPI_MAX, mpi_comm_level1);
  tmp_array  = (double*)malloc(sizeof(double)*Max_Size);
  tmp_array2 = (double*)malloc(sizeof(double)*Max_Size);

  /* send and recieve AtomDen_Grid */

  tag = 999;
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){

      /*****************************
              sending of data 
      *****************************/

      if (F_Snd_Num[IDS]!=0){

        /* find data size  */
        size1 = Snd_Size[IDS];

        /* multidimentional array to vector array */
        k = 0; 
        for (n=0; n<F_Snd_Num[IDS]; n++){
          Mc_AN = Snd_MAN[IDS][n];
          Gc_AN = Snd_GAN[IDS][n];

          for (i=0; i<GridN_Atom[Gc_AN]; i++){
            tmp_array[k] = AtomDen_Grid[Mc_AN][i];
            k++;
          }          
	} 

        /* MPI_Isend */
        MPI_Isend(&tmp_array[0], size1, MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);
      }

      /*****************************
         receiving of block data
      *****************************/

      if (F_Rcv_Num[IDR]!=0){

        /* find data size */
        size2 = Rcv_Size[IDR]; 

        /* MPI_Recv */
        MPI_Recv(&tmp_array2[0], size2, MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

        k = 0;
        Mc_AN = F_TopMAN[IDR] - 1;
        for (n=0; n<F_Rcv_Num[IDR]; n++){
          Mc_AN++;
          Gc_AN = Rcv_GAN[IDR][n];

          for (i=0; i<GridN_Atom[Gc_AN]; i++){
            AtomDen_Grid[Mc_AN][i] = tmp_array2[k];
            k++;
          }
        }
      }
      if (F_Snd_Num[IDS]!=0) MPI_Wait(&request,&stat);
    } 
  }  

  /******************************************************
   MPI:
         PCCDen_Grid
  ******************************************************/

  /* send and receive PCCDen_Grid */

  tag = 999;
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){

      /*****************************
              sending of data 
      *****************************/

      if (F_Snd_Num[IDS]!=0){

        /* find data size  */
        size1 = Snd_Size[IDS];

        /* multidimentional array to vector array */
        k = 0; 
        for (n=0; n<F_Snd_Num[IDS]; n++){
          Mc_AN = Snd_MAN[IDS][n];
          Gc_AN = Snd_GAN[IDS][n];

          for (i=0; i<GridN_Atom[Gc_AN]; i++){
            tmp_array[k] = PCCDen_Grid[Mc_AN][i];
            k++;
          }
	} 

        /* MPI_Isend */
        MPI_Isend(&tmp_array[0], size1, MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);
      }

      /*****************************
         receiving of block data
      *****************************/

      if (F_Rcv_Num[IDR]!=0){

        /* find data size */
        size2 = Rcv_Size[IDR]; 

        /* MPI_Recv */
        MPI_Recv(&tmp_array2[0], size2, MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

        k = 0;
        Mc_AN = F_TopMAN[IDR] - 1;
        for (n=0; n<F_Rcv_Num[IDR]; n++){
          Mc_AN++;
          Gc_AN = Rcv_GAN[IDR][n];

          for (i=0; i<GridN_Atom[Gc_AN]; i++){
            PCCDen_Grid[Mc_AN][i] = tmp_array2[k];
            k++;
          }          
        }
      }
      if (F_Snd_Num[IDS]!=0) MPI_Wait(&request,&stat);
    } 
  }  

  /****************************************************
    freeing of arrays:

     tmp_array
     tmp_array2
  ****************************************************/

  free(tmp_array);
  free(tmp_array2);

  /******************************************************
            superposition of atomic densities
  ******************************************************/

  for (Mc_AN=1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){

    dtime(&Stime_atom);

    Gc_AN = F_M2G[Mc_AN];
    Cwan = WhatSpecies[Gc_AN];
    Nele = InitN_USpin[Gc_AN] + InitN_DSpin[Gc_AN];
    Nu = InitN_USpin[Gc_AN];
    Nd = InitN_DSpin[Gc_AN];

    if (1.0e-15<Nele){

      ocupcy_u = Nu/Nele;
      ocupcy_d = Nd/Nele;
    }
    else{
      ocupcy_u = 0.0;
      ocupcy_d = 0.0;
    }
 
    if (2<=level_stdout){
      printf("  Mc_AN=%3d Gc_AN=%3d ocupcy_u=%15.12f ocupcy_d=%15.12f\n",
                Mc_AN,Gc_AN,ocupcy_u,ocupcy_d);
    }

    for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){

      MN = MGridListAtom[Mc_AN][Nc];
      GRc = CellListAtom[Mc_AN][Nc];

      if ( (0<=MN && Solver!=4) || (0<=MN && Solver==4 && atv_ijk[GRc][1]==0) ){

	DenA = AtomDen_Grid[Mc_AN][Nc];

	/* spin collinear */
	if ( init_density==1 && SpinP_switch==1 ){
	  Density_Grid[0][MN] += ocupcy_u*DenA;
	  Density_Grid[1][MN] += ocupcy_d*DenA;
	} 

	/* spin non-collinear */
	else if ( init_density==1 && SpinP_switch==3 ){

	  theta = Angle0_Spin[Gc_AN];
	  phi   = Angle1_Spin[Gc_AN];
          
	  rho = DenA;
	  mag = (ocupcy_u - ocupcy_d)*DenA;             
	  magx = mag*sin(theta)*cos(phi);
	  magy = mag*sin(theta)*sin(phi);
	  magz = mag*cos(theta);

	  Density_Grid[0][MN] += rho;
	  Density_Grid[1][MN] += magx;
	  Density_Grid[2][MN] += magy;
	  Density_Grid[3][MN] += magz;
	} 

	else if (init_density==1){
	  Density_Grid[0][MN] += 0.5*DenA;
	}
	else{
	  ADensity_Grid[MN] += 0.5*DenA;
	}

	/*  partial core correction  */
	if (PCC_switch==1) {
	  DenPCC = PCCDen_Grid[Mc_AN][Nc];
	  /* later add this in Set_XC_Grid */
	  PCCDensity_Grid[MN] += 0.5*DenPCC;
	}

      } 
    } /* Nc */

    dtime(&Etime_atom);
    time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
  }

  /******************************************************
    atomic densities (AtomDen2_Grid) in terms of FNAN2 
  ******************************************************/

  AtomDen2_Grid = (double*)malloc(sizeof(double)*FNAN2_Grid);

  /* MPI */
  tag = 999;
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){

      /*****************************
              sending of data 
      *****************************/

      if (Num_Snd_FNAN2_Grid[IDS]!=0){

        tmp_array = (double*)malloc(sizeof(double)*Num_Snd_FNAN2_Grid[IDS]);

        /* vector array */
        for (i=0; i<Num_Snd_FNAN2_Grid[IDS]; i++){
          Gc_AN = Snd_FNAN2_At[IDS][i];
          Mc_AN = F_G2M[Gc_AN];
          Nc    = Snd_FNAN2_Nc[IDS][i];
          tmp_array[i] = AtomDen_Grid[Mc_AN][Nc];
        }

        /* MPI_Isend */
        MPI_Isend(&tmp_array[0], Num_Snd_FNAN2_Grid[IDS], MPI_DOUBLE,
                  IDS, tag, mpi_comm_level1, &request);
      }

      /*****************************
            receiving of data
      *****************************/

      if (Num_Rcv_FNAN2_Grid[IDR]!=0){
        top_num = TopMAN2_Grid[IDR];
        /* MPI_Recv */
        MPI_Recv(&AtomDen2_Grid[top_num], Num_Rcv_FNAN2_Grid[IDR], MPI_DOUBLE,
                  IDR, tag, mpi_comm_level1, &stat);
      }

      if (Num_Snd_FNAN2_Grid[IDS]!=0){
         MPI_Wait(&request,&stat);
         free(tmp_array);
      }

    }
  }

  /* Density_Grid += AtomDen2_Grid */

  for (i=0; i<FNAN2_Grid; i++){

    DenA = AtomDen2_Grid[i];
    MN    = Rcv_FNAN2_MN[i];
    Gc_AN = Rcv_FNAN2_GA[i];
    GRc   = Rcv_FNAN2_GRc[i];
    Nele = InitN_USpin[Gc_AN] + InitN_DSpin[Gc_AN];
    Nu = InitN_USpin[Gc_AN];
    Nd = InitN_DSpin[Gc_AN];

    if (1.0e-15<Nele){
      ocupcy_u = Nu/Nele;
      ocupcy_d = Nd/Nele;
    }
    else{
      ocupcy_u = 0.0;
      ocupcy_d = 0.0;
    }

    if (Solver!=4 || (Solver==4 && atv_ijk[GRc][1]==0 )){

      /* spin collinear */
      if ( init_density==1 && SpinP_switch==1 ){
	Density_Grid[0][MN] += ocupcy_u*DenA;
	Density_Grid[1][MN] += ocupcy_d*DenA;
      } 

      /* spin non-collinear */
      else if ( init_density==1 && SpinP_switch==3 ){

	theta = Angle0_Spin[Gc_AN];
	phi   = Angle1_Spin[Gc_AN];

	rho = DenA;
	mag = (ocupcy_u - ocupcy_d)*DenA;
	magx = mag*sin(theta)*cos(phi);
	magy = mag*sin(theta)*sin(phi);
	magz = mag*cos(theta);

	Density_Grid[0][MN] += rho;
	Density_Grid[1][MN] += magx;
	Density_Grid[2][MN] += magy;
	Density_Grid[3][MN] += magz;
      } 

      else if (init_density==1){
	Density_Grid[0][MN] += 0.5*DenA;
      }

    }
  }

  /******************************************************
    pcc densities (PCCDen2_Grid) in terms of FNAN2 
  ******************************************************/

  PCCDen2_Grid = (double*)malloc(sizeof(double)*FNAN2_Grid);

  if (PCC_switch==1) {

    /* MPI */
    tag = 999;
    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      if (ID!=0){

        /*****************************
                sending of data 
        *****************************/

        if (Num_Snd_FNAN2_Grid[IDS]!=0){

          tmp_array = (double*)malloc(sizeof(double)*Num_Snd_FNAN2_Grid[IDS]);

          /* vector array */
          for (i=0; i<Num_Snd_FNAN2_Grid[IDS]; i++){
            Gc_AN = Snd_FNAN2_At[IDS][i];
            Mc_AN = F_G2M[Gc_AN];
            Nc    = Snd_FNAN2_Nc[IDS][i];
            tmp_array[i] = PCCDen_Grid[Mc_AN][Nc];
         }

          /* MPI_Isend */
          MPI_Isend(&tmp_array[0], Num_Snd_FNAN2_Grid[IDS], MPI_DOUBLE,
                    IDS, tag, mpi_comm_level1, &request);
        }

        /*****************************
              receiving of data
        *****************************/

        if (Num_Rcv_FNAN2_Grid[IDR]!=0){
          top_num = TopMAN2_Grid[IDR];
          /* MPI_Recv */
          MPI_Recv(&PCCDen2_Grid[top_num], Num_Rcv_FNAN2_Grid[IDR], MPI_DOUBLE,
                    IDR, tag, mpi_comm_level1, &stat);
        }

        if (Num_Snd_FNAN2_Grid[IDS]!=0){
           MPI_Wait(&request,&stat);
           free(tmp_array);
        }

      }
    }

    /* PCCDensity_Grid += PCCDen2_Grid */
    for (i=0; i<FNAN2_Grid; i++){
      MN = Rcv_FNAN2_MN[i];
      PCCDensity_Grid[MN] += 0.5*PCCDen2_Grid[i];
    }
  }

  /****************************************************
     initialize diagonal and off-diagonal densities
           in case of spin non-collinear DFT
  ****************************************************/

  if (init_density==1 && SpinP_switch==3){
    for (MN=0; MN<My_NumGrid1; MN++){

      rho  = Density_Grid[0][MN];
      magx = Density_Grid[1][MN];
      magy = Density_Grid[2][MN];
      magz = Density_Grid[3][MN];

      Density_Grid[0][MN] = 0.5*(rho + magz);
      Density_Grid[1][MN] = 0.5*(rho - magz);
      Density_Grid[2][MN] = 0.5*magx;
      Density_Grid[3][MN] =-0.5*magy;
    }
  }

  /******************************************************
               Density_Grid to ADensity_Grid
  ******************************************************/

  if ( init_density==1 && (SpinP_switch==1 || SpinP_switch==3) ){
    for (MN=0; MN<My_NumGrid1; MN++){
      ADensity_Grid[MN] = 0.5*(Density_Grid[0][MN] + Density_Grid[1][MN]);
    }
  }
  else if (init_density==1){
    for (MN=0; MN<My_NumGrid1; MN++){
      ADensity_Grid[MN] = Density_Grid[0][MN];
    }
  } 

  /****************************************************
            in case of non-spin polarization
                up-spin to down-spin
  ****************************************************/

  if (init_density==1 && SpinP_switch==0){
  
    for (MN=0; MN<My_NumGrid1; MN++){
      Density_Grid[1][MN] = Density_Grid[0][MN];
    }
  }

  /****************************************************
    freeing of arrays:

     Snd_Size
     Rcv_Size
     AtomDen_Grid
     PCCDen_Grid
  ****************************************************/

  free(Snd_Size);
  free(Rcv_Size);

  for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
    free(AtomDen_Grid[Mc_AN]);
  }
  free(AtomDen_Grid);

  for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
    free(PCCDen_Grid[Mc_AN]);
  }
  free(PCCDen_Grid);

  free(AtomDen2_Grid);
  free(PCCDen2_Grid);

  /* elapsed time */
  dtime(&TEtime);
  time0 = TEtime - TStime;
  return time0;
}
