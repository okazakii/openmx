/**********************************************************************
  Set_Orbitals_Grid.c:

   Set_Orbitals_Grid.c is a subroutine to calculate the value of basis
   functions on each grid point.

  Log of Set_Orbitals_Grid.c:

     22/Nov/2001  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
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



double Set_Orbitals_Grid(int Cnt_kind)
{
  int i,n,Mc_AN,Gc_AN,Cwan,NO0,GNc,GRc;
  long int size1,size2; 
  long int k,Nc;
  long int My_Max,Max_Size;
  double time0;
  double x,y,z,dx,dy,dz;
  double *Chi0;
  double TStime,TEtime;
  double Cxyz[4];
  Type_Orbs_Grid *tmp_array;
  Type_Orbs_Grid *tmp_array2;
  int *Snd_Size,*Rcv_Size;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  double Stime_atom, Etime_atom;

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
    allocation of array:

    int Snd_Size[numprocs]
    int Rcv_Size[numprocs]
  ****************************************************/

  Snd_Size = (int*)malloc(sizeof(int)*numprocs); 
  Rcv_Size = (int*)malloc(sizeof(int)*numprocs); 

  /*****************************************************
                Calculate orbitals on grids
  *****************************************************/

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    dtime(&Stime_atom);

    Gc_AN = M2G[Mc_AN];    
    Cwan = WhatSpecies[Gc_AN];

    if (Cnt_kind==0)
      NO0 = Spe_Total_NO[Cwan];
    else
      NO0 = Spe_Total_CNO[Cwan]; 

#pragma omp parallel shared(Orbs_Grid,Cnt_kind,Gxyz,atv,CellListAtom,GridListAtom,GridN_Atom,Gc_AN,Cwan,Mc_AN,NO0) private(OMPID,Nthrds,Nprocs,Nc,GNc,GRc,Cxyz,x,y,z,dx,dy,dz,Chi0,i)
    {
 
      /* allocation of array */

      Chi0 = (double*)malloc(sizeof(double)*List_YOUSO[7]);

      /* get info. on OpenMP */ 

      OMPID = omp_get_thread_num();
      Nthrds = omp_get_num_threads();
      Nprocs = omp_get_num_procs();

      for (Nc=OMPID*GridN_Atom[Gc_AN]/Nthrds; Nc<(OMPID+1)*GridN_Atom[Gc_AN]/Nthrds; Nc++){

	GNc = GridListAtom[Mc_AN][Nc]; 
	GRc = CellListAtom[Mc_AN][Nc];

	Get_Grid_XYZ(GNc,Cxyz);
	x = Cxyz[1] + atv[GRc][1]; 
	y = Cxyz[2] + atv[GRc][2]; 
	z = Cxyz[3] + atv[GRc][3];

	dx = x - Gxyz[Gc_AN][1];
	dy = y - Gxyz[Gc_AN][2];
	dz = z - Gxyz[Gc_AN][3];

	/* no contraction */
	if (Cnt_kind==0){
	  Get_Orbitals(Cwan,dx,dy,dz,Chi0);
	} 

	/* contraction */
	else{
	  Get_Cnt_Orbitals(Mc_AN,dx,dy,dz,Chi0);
	}

	for (i=0; i<NO0; i++){
	  Orbs_Grid[Mc_AN][i][Nc] = (Type_Orbs_Grid)Chi0[i];
	}

      } /* Nc */

      /* freeing of array */

      free(Chi0);

    } /* #pragma omp parallel */

    dtime(&Etime_atom);
    time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
  }

  /****************************************************
     MPI:

     Orbs_Grid
  ****************************************************/

  /* find data size for sending and receiving */

  tag = 999;
  My_Max = -10000;
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){
      /*  sending size */
      if (F_Snd_Num[IDS]!=0){
        /* find data size */ 
        size1 = 0; 
        for (n=0; n<F_Snd_Num[IDS]; n++){
          Gc_AN = Snd_GAN[IDS][n];
          Cwan = WhatSpecies[Gc_AN];

          if (Cnt_kind==0)
            size1 += GridN_Atom[Gc_AN]*Spe_Total_NO[Cwan];
          else 
            size1 += GridN_Atom[Gc_AN]*Spe_Total_CNO[Cwan];
	}

        Snd_Size[IDS] = size1;
        MPI_Isend(&size1, 1, MPI_LONG, IDS, tag, mpi_comm_level1, &request);
      }
      else{
        Snd_Size[IDS] = 0;
      }

      /*  receiving size */
      if (F_Rcv_Num[IDR]!=0){
        MPI_Recv(&size2, 1, MPI_LONG, IDR, tag, mpi_comm_level1, &stat);
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

  MPI_Allreduce(&My_Max, &Max_Size, 1, MPI_LONG, MPI_MAX, mpi_comm_level1);
  tmp_array  = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*Max_Size);
  tmp_array2 = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*Max_Size);

  /* send and receive Orbs_Grid */

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
          Cwan = WhatSpecies[Gc_AN];

          if (Cnt_kind==0)
            NO0 = Spe_Total_NO[Cwan];
          else
            NO0 = Spe_Total_CNO[Cwan]; 

          for (i=0; i<NO0; i++){
            for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
              tmp_array[k] = Orbs_Grid[Mc_AN][i][Nc];
              k++;
            }          
          }
	} 

        /* MPI_Isend */

        MPI_Isend(&tmp_array[0], size1, MPI_Type_Orbs_Grid, IDS, tag, mpi_comm_level1, &request);
      }

      /*****************************
         receiving of block data
      *****************************/

      if (F_Rcv_Num[IDR]!=0){

        /* find data size */
        size2 = Rcv_Size[IDR]; 

        /* MPI_Recv */

        MPI_Recv(&tmp_array2[0], size2, MPI_Type_Orbs_Grid, IDR, tag, mpi_comm_level1, &stat);

        k = 0;
        Mc_AN = F_TopMAN[IDR] - 1;
        for (n=0; n<F_Rcv_Num[IDR]; n++){
          Mc_AN++;
          Gc_AN = Rcv_GAN[IDR][n];
          Cwan = WhatSpecies[Gc_AN];

          if (Cnt_kind==0)
            NO0 = Spe_Total_NO[Cwan];
          else
            NO0 = Spe_Total_CNO[Cwan]; 


          for (i=0; i<NO0; i++){
            for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
              Orbs_Grid[Mc_AN][i][Nc] = tmp_array2[k];
              k++;
            }          
          }

        }
      }

      if (F_Snd_Num[IDS]!=0) MPI_Wait(&request,&stat);
    } 
  }  

  /****************************************************
     freeing of array:

     tmp_array
     tmp_array2
     Snd_Size
     Rcv_Size
  ****************************************************/

  free(tmp_array);
  free(tmp_array2);

  free(Snd_Size);
  free(Rcv_Size);

  /* time */

  dtime(&TEtime);
  time0 = TEtime - TStime;

  return time0;
}
