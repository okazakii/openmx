/**********************************************************************
  BroadCast_ComplexMatrix.c:

     BroadCast_ComplexMatrix.c is a subroutine to broadcast a matrix "Mat"
     which is distributed by row in each processor.

  Log of BroadCast_ComplexMatrix.c:

     14/Dec/2004  Released by T.Ozaki

***********************************************************************/

#include <stdio.h> 
#include <math.h>
#include <stdlib.h>
#include "openmx_common.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

void BroadCast_ComplexMatrix(MPI_Comm MPI_Curret_Comm_WD, 
                             dcomplex **Mat, int n, int *is1, int *ie1, int myid, int numprocs, 
                             MPI_Status *stat_send,
                             MPI_Request *request_send,
                             MPI_Request *request_recv)
{
  int num0,num1,k0,k1;
  int i,j,k,num,ID,tag=999;
  double *Mat1;

  MPI_Status stat;
  MPI_Request request;

  /*********************************************
     Elemements are stored from 1 to n in Mat. 
  **********************************************/

  if (numprocs!=1){

    Mat1 = (double*)malloc(sizeof(double)*(n+1)*(n+1));

    /********************************
           Real part of Mat 
    ********************************/

    for (i=is1[myid]; i<=ie1[myid]; i++){
      for (j=1; j<=n; j++){
	k = (i-1)*n + j - 1;
	Mat1[k] = Mat[i][j].r;
      }
    }

    /* sending */

    k0 = (is1[myid]-1)*n;
    if (k0<0) k0 = 0;  
    num0 = (ie1[myid] - is1[myid] + 1)*n;
    if (num0<0) num0 = 0;

    for (ID=0; ID<numprocs; ID++){
      MPI_Isend(&Mat1[k0], num0, MPI_DOUBLE, ID, tag, MPI_Curret_Comm_WD, &request_send[ID]);
    }

    /* receiving */

    for (ID=0; ID<numprocs; ID++){
      k1 = (is1[ID]-1)*n;
      if (k1<0) k1 = 0;  
      num1 = (ie1[ID] - is1[ID] + 1)*n;
      if (num1<0) num1 = 0;
      MPI_Irecv(&Mat1[k1], num1, MPI_DOUBLE, ID, tag, MPI_Curret_Comm_WD, &request_recv[ID]);
    }

    /* waitall */

    MPI_Waitall(numprocs,request_recv,stat_send);
    MPI_Waitall(numprocs,request_send,stat_send);

    for (ID=0; ID<numprocs; ID++){
      for (i=is1[ID]; i<=ie1[ID]; i++){
	for (j=1; j<=n; j++){
	  k = (i-1)*n + j - 1;
	  Mat[i][j].r = Mat1[k];
	}
      }
    }

    /********************************
          Imaginary part of Mat 
    ********************************/

    for (i=is1[myid]; i<=ie1[myid]; i++){
      for (j=1; j<=n; j++){
	k = (i-1)*n + j - 1;
	Mat1[k] = Mat[i][j].i;
      }
    }

    /* sending */

    k0 = (is1[myid]-1)*n;
    if (k0<0) k0 = 0;  
    num0 = (ie1[myid] - is1[myid] + 1)*n;
    if (num0<0) num0 = 0;

    for (ID=0; ID<numprocs; ID++){
      MPI_Isend(&Mat1[k0], num0, MPI_DOUBLE, ID, tag, MPI_Curret_Comm_WD, &request_send[ID]);
    }

    /* receiving */

    for (ID=0; ID<numprocs; ID++){
      k1 = (is1[ID]-1)*n;
      if (k1<0) k1 = 0;  
      num1 = (ie1[ID] - is1[ID] + 1)*n;
      if (num1<0) num1 = 0;
      MPI_Irecv(&Mat1[k1], num1, MPI_DOUBLE, ID, tag, MPI_Curret_Comm_WD, &request_recv[ID]);
    }

    /* waitall */

    MPI_Waitall(numprocs,request_recv,stat_send);
    MPI_Waitall(numprocs,request_send,stat_send);

    for (ID=0; ID<numprocs; ID++){
      for (i=is1[ID]; i<=ie1[ID]; i++){
	for (j=1; j<=n; j++){
	  k = (i-1)*n + j - 1;
	  Mat[i][j].i = Mat1[k];
	}
      }
    }

    free(Mat1);
  }
}


