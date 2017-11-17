/**********************************************************************
    lapack_dstedc1.c:

    lapack_dstedc1.c is a subroutine to find eigenvalues and eigenvectors
    of tridiagonlized real matrix using lapack's routines dstedc.

    Log of lapack_dstedc1.c:

       Dec/24/2004  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "openmx_common.h"
#include "lapack_prototypes.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

void lapack_dstedc1(INTEGER N, double *D, double *E, double *W, double **ev)
{
  int i,j;
  char  *COMPZ="I";
  INTEGER IL,IU; 
  double ABSTOL_GR=1.0e-14;
  double *Z;
  INTEGER LDZ;
  INTEGER LWORK;
  double *WORK;
  INTEGER *IWORK;
  INTEGER LIWORK;
  INTEGER INFO;

  LDZ = N;

  Z = (double*)malloc(sizeof(double)*LDZ*N);

  LWORK = 1 + 4*N + N*N;
  LIWORK = 3 + 5*N;
  WORK = (double*)malloc(sizeof(double)*LWORK);
  IWORK = (INTEGER*)malloc(sizeof(INTEGER)*LIWORK);

  F77_NAME(dstedc,DSTEDC)( COMPZ, &N, D, E, Z, &LDZ, WORK, &LWORK, IWORK, &LIWORK, &INFO );

  /* store eigenvectors */

  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      ev[i+1][j+1]= Z[i*N+j];
    }
  }

  /* shift ko by 1 */
  for (i=N; i>=1; i--){
    W[i]= D[i-1];
  }

  if (INFO>0) {
    /*
    printf("\n error in dstedc_, info=%d\n\n",INFO);
    */
  }
  if (INFO<0) {
    printf("info=%d in dstedc_\n",INFO);
    MPI_Finalize();
    exit(0);
  }

  free(Z);
  free(WORK);
  free(IWORK);
}
