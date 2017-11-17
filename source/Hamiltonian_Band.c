/**********************************************************************
  Hamiltonian_Band.c:

     Hamiltonian_Band.c is a subroutine to make a Hamiltonian matrix
     for a periodic boundary system using Bloch theorem.

  Log of Hamiltonian_Band.c:

     22/Nov/2001  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "openmx_common.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

void Hamiltonian_Band(int Host_ID1,
                      double ****RH, 
                      dcomplex **H, int *MP,
                      double k1, double k2, double k3)
{
  static int firsttime=1;
  int i,j,k,wanA,wanB,tnoA,tnoB,Anum,Bnum;
  int NUM,MA_AN,GA_AN,LB_AN,GB_AN;
  int l1,l2,l3,Rn,n2;
  double **H1,**H2;
  double *tmp_array1,*tmp_array2;
  double kRn,si,co,rh,ih;
  int ID,myid,numprocs,tag=999;

  MPI_Status stat;
  MPI_Request request;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  /* set MP */
  Anum = 1;
  for (i=1; i<=atomnum; i++){
    MP[i] = Anum;
    wanA = WhatSpecies[i];
    tnoA = Spe_Total_CNO[wanA];
    Anum = Anum + tnoA;
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

  if (firsttime) 
  PrintMemory("Hamiltonian_Band: H1",sizeof(double)*n2*n2,NULL);

  H2 = (double**)malloc(sizeof(double*)*n2);
  for (i=0; i<n2; i++){
    H2[i] = (double*)malloc(sizeof(double)*n2);
  }
  if (firsttime)
  PrintMemory("Hamiltonian_Band: H2",sizeof(double)*n2*n2,NULL);

  /* for PrintMemory */
  firsttime=0;

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

  for (MA_AN=1; MA_AN<=Matomnum; MA_AN++){
    GA_AN = M2G[MA_AN];    
    wanA = WhatSpecies[GA_AN];
    tnoA = Spe_Total_CNO[wanA];
    Anum = MP[GA_AN];

    for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
      GB_AN = natn[GA_AN][LB_AN];
      Rn = ncn[GA_AN][LB_AN];
      wanB = WhatSpecies[GB_AN];
      tnoB = Spe_Total_CNO[wanB];

/*
      kRn = k1*( (double)atv_ijk[Rn][1] + Cell_Gxyz[GB_AN][1] - Cell_Gxyz[GA_AN][1] )
          + k2*( (double)atv_ijk[Rn][2] + Cell_Gxyz[GB_AN][2] - Cell_Gxyz[GA_AN][2] )
          + k3*( (double)atv_ijk[Rn][3] + Cell_Gxyz[GB_AN][3] - Cell_Gxyz[GA_AN][3] );
*/

      l1 = atv_ijk[Rn][1];
      l2 = atv_ijk[Rn][2];
      l3 = atv_ijk[Rn][3];
      kRn = k1*(double)l1 + k2*(double)l2 + k3*(double)l3;

      si = sin(2.0*PI*kRn);
      co = cos(2.0*PI*kRn);
      Bnum = MP[GB_AN];

      for (i=0; i<=(tnoA-1); i++){
	for (j=0; j<=(tnoB-1); j++){
	  rh = RH[MA_AN][LB_AN][i][j];
	  H1[Anum+i][Bnum+j] += rh*co;
	  H2[Anum+i][Bnum+j] += rh*si;
	}
      }
    }
  }

  /******************************************************
    MPI: H1 and H2
  ******************************************************/

  tmp_array1 = (double*)malloc(sizeof(double)*2*n2*List_YOUSO[7]);
  tmp_array2 = (double*)malloc(sizeof(double)*2*n2*List_YOUSO[7]);

  /* H1 and H2 */

  if (myid!=Host_ID1){
    for (MA_AN=1; MA_AN<=Matomnum; MA_AN++){
      GA_AN = M2G[MA_AN];    
      wanA = WhatSpecies[GA_AN];
      tnoA = Spe_Total_CNO[wanA];
      Anum = MP[GA_AN];

      k = 0;  

      for (i=0; i<tnoA; i++){
        for (j=0; j<n2; j++){
          tmp_array1[k] = H1[Anum+i][j];          
	  k++; 
	}
      }
      for (i=0; i<tnoA; i++){
        for (j=0; j<n2; j++){
          tmp_array1[k] = H2[Anum+i][j];          
	  k++; 
	}
      }

      tag = 999;
      MPI_Isend(&tmp_array1[0], 2*tnoA*n2, MPI_DOUBLE, Host_ID1,
                tag, mpi_comm_level1, &request);
      MPI_Wait(&request,&stat);

    }
  }
  else{
    for (GA_AN=1; GA_AN<=atomnum; GA_AN++){
      wanA = WhatSpecies[GA_AN];
      tnoA = Spe_Total_CNO[wanA];
      Anum = MP[GA_AN];
      ID = G2ID[GA_AN];
      if (ID!=Host_ID1){

        tag = 999;
        MPI_Recv(&tmp_array2[0], 2*tnoA*n2, MPI_DOUBLE, ID, tag, mpi_comm_level1, &stat);

	k = 0;
        for (i=0; i<tnoA; i++){
          for (j=0; j<n2; j++){
            H1[Anum+i][j] = tmp_array2[k];
	    k++;
	  }
	}
        for (i=0; i<tnoA; i++){
          for (j=0; j<n2; j++){
            H2[Anum+i][j] = tmp_array2[k];
	    k++;
	  }
	}

      }    
    }
  }

  free(tmp_array1);
  free(tmp_array2);

  /******************************************************
    the full complex matrix of H
  ******************************************************/

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
  }
  free(H1);

  for (i=0; i<n2; i++){
    free(H2[i]);
  }
  free(H2);

}

