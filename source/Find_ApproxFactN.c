/**********************************************************************
  Find_ApproxFactN.c

   Find_ApproxFactN.c is a subrutine to find the number of grids along 
   the a-, b-, and c-axes which satisies the required cutoff energy
   approximately.

  Log of Find_ApproxFactN.c:

     9/Dec./2006  Released by H. Kino

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

#include "tran_prototypes.h"

int factorize(int num0, int N, int *fund, int *pow, int must );


void Find_ApproxFactN(double tv[4][4],double *GEcut,
                      int *Ng1,int *Ng2,int *Ng3,
                      double *A2, double *B2, double *C2)
{
  int i,iaxis;
  double fac,Ecut;
  int imin,imax;
  double tolfac=1.50;

  /* radix of FFT */
  int mustfact[3]; 
  int pow[4];
  int num,must[4];

  int  numprocs,myid;

  int numGrid0[3];
  double Ecut0[3];

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  numGrid0[0]= *Ng1;
  numGrid0[1]= *Ng2;
  numGrid0[2]= *Ng3;
  Ecut0[0] = *A2;
  Ecut0[1] = *B2;
  Ecut0[2] = *C2;

  /* set a divisible number by numprocs for Ng1 */
  if (MPI_tunedgrid_flag){
    must[0]=numprocs;
    must[1]=1;
    must[2]=1;
  }
  /* set an optimized number in one cpu for Ng1 */
  else{
    must[0]=1;
    must[1]=1;
    must[2]=1;
  }

  for (iaxis=0; iaxis<3; iaxis++) {

    if (2<=numGrid0[iaxis]) {

      for (num=numGrid0[iaxis]/2;;num++) {

	if ( factorize(num,NfundamentalNum,fundamentalNum,pow,must[iaxis]) ) {

	  fac = (double)num/(double)numGrid0[iaxis];
	  Ecut = Ecut0[iaxis] * fac*fac;

	  if ( Ecut > *GEcut ) {
            numGrid0[iaxis]=num;
            Ecut0[iaxis] = Ecut;
            break;
	  }
	}
      }
    }
  }

  /* tune GridX */

  imin =0;
  for (i=1;i<3;i++) {
    if (Ecut0[imin] > Ecut0[i] ) imin=i;
  }

  imax =0;
  for (i=1;i<3;i++) {
    if (Ecut0[imax] < Ecut0[i] ) imax=i;
  }

  if (Ecut0[imax]/Ecut0[imin] > tolfac) { tolfac = Ecut0[imax]/Ecut0[imin]; }

  /* fine tuning */
   
  if ( (Ecut0[imax]/Ecut0[imin]) > 1.10 ){

    int Fnum=10;
    int FGrid[3][10];
    double FEcut[3][10];

    int maxnum[3];
    int inum;
    double score[10][10][10],minval;
    int i,j,k;
    int myid,numprocs;
    int must[3]; /* for factorize */

    MPI_Comm_size(mpi_comm_level1,&numprocs);
    MPI_Comm_rank(mpi_comm_level1,&myid);

    if (myid==Host_ID && 0<level_stdout){    
      printf("<Find_ApproxFactN> diff. of Cutoff is too large. tune Cutoff(Ng)\n");
    }

    /* set a divisible number by numprocs for Ng1 */
    if (MPI_tunedgrid_flag){
      must[0]=numprocs;
      must[1]=1;
      must[2]=1;
    }
    /* set an optimized number in one cpu for Ng1 */
    else{
      must[0]=1;
      must[1]=1;
      must[2]=1;
    }

    /* find 1st- Fnum-th  Grid larger than numGrid */

    for (iaxis=0;iaxis<3;iaxis++) {

      inum=0;
      FGrid[iaxis][inum] = numGrid0[iaxis];
      FEcut[iaxis][inum] = Ecut0[iaxis];
      inum++; /* num. of grids found so far */

      for (num=numGrid0[iaxis]+1; ; num++) {

	if ( factorize(num,NfundamentalNum,fundamentalNum,pow,must[iaxis]) ) {
          fac = (double)num/(double)numGrid0[iaxis];
          Ecut = Ecut0[iaxis] * fac*fac;
          FGrid[iaxis][inum] = num;
	  FEcut[iaxis][inum] = Ecut;

	  if ( (inum+1)>=Fnum|| FEcut[iaxis][inum] > (*GEcut)*tolfac) {
	    maxnum[iaxis]=inum + 1;
	    break;
	  }
	  inum++;
        }
      }

      if (myid==Host_ID && 0<level_stdout){    
        for(i=0;i<maxnum[iaxis];i++) {
	  printf("%f(%i) ",FEcut[iaxis][i],FGrid[iaxis][i]);
        }
        printf("\n");
      }
    }

    for (i=0;i<maxnum[0];i++) {
      for (j=0;j<maxnum[1];j++) {
        for (k=0;k<maxnum[2];k++) {
           score[i][j][k] = ( sqr( FEcut[0][i]-FEcut[1][j] ) 
                            + sqr( FEcut[0][i]-FEcut[2][k] ) 
                            + sqr( FEcut[1][j]-FEcut[2][k] )
                            + 0.0001*sqr( FEcut[0][i]-*GEcut )
                            + 0.0001*sqr( FEcut[1][j]-*GEcut )
                            + 0.0001*sqr( FEcut[2][k]-*GEcut )    );
        }
      }
    }
    minval=score[0][0][0];
    
    for (i=0;i<maxnum[0];i++) {
      for (j=0;j<maxnum[1];j++) {
        for (k=0;k<maxnum[2];k++) {
         /* printf("%i %i %i score=%lf\n",i,j,k,score[i][j][k]); */
          if (minval > score[i][j][k] ) { minval=score[i][j][k] ; }
        }
      }
    }

    for (i=0;i<maxnum[0];i++) {
      for (j=0;j<maxnum[1];j++) {
        for (k=0;k<maxnum[2];k++) {
          if ( fabs(minval-score[i][j][k])<1.0e-50 ) goto foundmin;
        }
      }
    }
     
foundmin:

    numGrid0[0]= FGrid[0][i];
    numGrid0[1]= FGrid[1][j];
    numGrid0[2]= FGrid[2][k];
    Ecut0[0] = FEcut[0][i];
    Ecut0[1] = FEcut[1][j];
    Ecut0[2] = FEcut[2][k];

  }

  *Ng1 = numGrid0[0];
  *Ng2 = numGrid0[1];
  *Ng3 = numGrid0[2];
  *A2 = Ecut0[0];
  *B2 = Ecut0[1];
  *C2 = Ecut0[2];

}



int factorize(int num0, int N, int *fund, int *pow, int must )
{
  int i;
  int a,b;
  int num;
  int ret;

  /* must exclude  division 0 */
  if (must==0)  return 0;

  if (must==1) {

    for (i=0;i<N;i++) {
      pow[i] = 0;
    }

  }
  else {
    num=num0%must;
    if ( num==0 ) {
      ret=factorize(  must, N, fund, pow, 1);
      if (ret==0) {
	return ret;
      }
    }
    else {
      return 0;
    }
  }

  num=num0/must ;
  for (i=0; i<N; i++) {
    while (1) {
      b = num%fund[i];

      if (b==0) {
	num /= fund[i];
	pow[i]++;
      }
      else {
	break;
      }
    }

  }

  if (num==1) {
    return num0;
  }
  else {
    return 0;
  }

}
 
