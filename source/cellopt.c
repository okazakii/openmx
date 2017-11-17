/**********************************************************************
  cellopt.c:

     cellopt.c is a subroutine to perform cell optimization.

  Log of cellopt.c:

    07/July/2014  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "openmx_common.h"
#include "Inputtools.h"
#include "lapack_prototypes.h"
#include "mpi.h"
#include <omp.h>

static void CellTet(char *argv[], double **CompTime);
static void CellCub(char *argv[], double **CompTime);
  

void cellopt(char *argv[], double **CompTime) 
{ 
  int i,j,flag;
  char *s_vec[40];
  int i_vec[40];

  if (CellOpt_switch==1) CellCub(argv,CompTime);
  if (CellOpt_switch==2) CellTet(argv,CompTime);

  MPI_Finalize();
  exit(0);
} 


void CellCub(char *argv[], double **CompTime)
{
  int i,j,I,INFO,Gc_AN,po;
  int *IPIV,LDA,LDB,N,NRHS;
  int MD_iter,CellOpt_iter,CellOpt_Maxiter;
  int myid,numprocs,ID;
  double His_Utot[300],His_Vol[300];
  double *a,*b,len,len0,gf,gf0,gf1,v_step,v,v0,v1;
  double diff0,diff1;
  static double coef = 0.01;
  double tv0[4][4],Cell_Volume0;
  static int firsttime=1;
  static char fileE[YOUSO10] = ".ene"; 
  static char fileDRC[YOUSO10] = ".md";

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /***********************************
    allocation of arrays
  ***********************************/

  a = (double*)malloc(sizeof(double)*20);
  b = (double*)malloc(sizeof(double)*20);
  IPIV = (INTEGER*)malloc(sizeof(INTEGER)*20);

  /***********************************
    assumption: 
    a = b = c
  ***********************************/

  CellOpt_iter = 1;
  CellOpt_Maxiter = 200;

  do {
    
    MD_iter = 1;
    
    do {
    
      CompTime[myid][2] += truncation(MD_iter,1);
      CompTime[myid][3] += DFT(MD_iter,(MD_iter-1)%orbitalOpt_per_MDIter+1);
      iterout(MD_iter+MD_Current_Iter,MD_TimeStep*(MD_iter+MD_Current_Iter-1),filepath,filename);
      CompTime[myid][4] += MD_pac(MD_iter,argv[1]);
      MD_iter++;
    
    } while(MD_Opt_OK==0 && MD_iter<=MD_IterNumber);
    
    His_Utot[CellOpt_iter] = Utot;
    His_Vol[CellOpt_iter]  = Cell_Volume;

    if (myid==Host_ID){
      printf("CellCub %2d %15.12f %15.12f\n",
	     CellOpt_iter,
	     His_Vol[CellOpt_iter],His_Utot[CellOpt_iter]);
    }

    if (CellOpt_iter<=3){

      if (firsttime){
      
	/* store tv into tv0 */ 
	for (i=1; i<=3; i++){
	  for (j=1; j<=3; j++){
	    tv0[i][j] = tv[i][j];
	  }
	}
	Cell_Volume0 = Cell_Volume;
	firsttime = 0;
      }

      /* update tv */ 
      for (i=1; i<=3; i++){
	for (j=1; j<=3; j++){
	  tv[i][j] += 0.01*tv0[i][j];
	}
      }

    } /* if (CellOpt_iter<=3) */

    /* if (3<CellOpt_iter) */

    else{

      /***********************************************************
         least square fitting using a polynomial function 
         including up to third-order terms 
      ***********************************************************/

      /* making the matrix 'a' */

      for (j=0; j<4; j++){
        for (i=0; i<4; i++){

          a[4*j+i] = 0.0; 
          for (I=1; I<=CellOpt_iter; I++){
            a[4*j+i] += pow(His_Vol[I],(double)i)*pow(His_Vol[I],(double)j);
	  }      
        }
      }

      /* making the vector 'b' */

      for (i=0; i<4; i++){

	b[i] = 0.0; 
	for (I=1; I<=CellOpt_iter; I++){
	  b[i] += His_Utot[I]*pow(His_Vol[I],(double)i);
	}      
      }

      /* solving the linear equation */

      N = 4;
      NRHS = 1;
      LDA = N;
      LDB = N;

      F77_NAME(dgesv,DGESV)(&N, &NRHS, a, &LDA, IPIV, b, &LDB, &INFO);


      /*
      for (i=0; i<4; i++){
        printf("i=%2d b=%10.7f\n",i,b[i]);
      }      
      MPI_Finalize();
      exit(0);      
      */

      /************************************************
          update the volume and lattice constant 
      ************************************************/

      diff0 = His_Vol[CellOpt_iter  ] - His_Vol[CellOpt_iter-1];
      diff1 = His_Vol[CellOpt_iter-1] - His_Vol[CellOpt_iter-2];

      if (diff0*diff1<0.0) coef *= 0.5;

      if (myid==Host_ID){
        printf("CellCub %7.3f %7.3f %7.3f  %7.3f %7.3f coef=%10.5f\n",
               His_Vol[CellOpt_iter  ],His_Vol[CellOpt_iter-1],His_Vol[CellOpt_iter-2],diff0,diff1,coef);
      }

      v = Cell_Volume;
      gf = b[1] + 2.0*b[2]*v + 3.0*b[3]*v*v;
      v = v - coef*v*sgn(gf);

      if (myid==Host_ID){
        printf("CellCub gf=%15.12f v=%15.12f\n",gf,v);
      }
   
      /* update the cell vectors */

      len = pow(v,0.33333333333333333333);

      len0 = sqrt(tv[1][1]*tv[1][1] + tv[1][2]*tv[1][2] + tv[1][3]*tv[1][3]); 
      tv[1][1] = len*tv[1][1]/len0; tv[1][2] = len*tv[1][2]/len0; tv[1][3] = len*tv[1][3]/len0;

      len0 = sqrt(tv[2][1]*tv[2][1] + tv[2][2]*tv[2][2] + tv[2][3]*tv[2][3]); 
      tv[2][1] = len*tv[2][1]/len0; tv[2][2] = len*tv[2][2]/len0; tv[2][3] = len*tv[2][3]/len0;  

      len0 = sqrt(tv[3][1]*tv[3][1] + tv[3][2]*tv[3][2] + tv[3][3]*tv[3][3]); 
      tv[3][1] = len*tv[3][1]/len0; tv[3][2] = len*tv[3][2]/len0; tv[3][3] = len*tv[3][3]/len0;  
    }

    /***********************************************************
                      update cartesian coordinates
    ***********************************************************/

    for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){

      Gxyz[Gc_AN][1] =  Cell_Gxyz[Gc_AN][1]*tv[1][1]
                      + Cell_Gxyz[Gc_AN][2]*tv[2][1]
                      + Cell_Gxyz[Gc_AN][3]*tv[3][1] + Grid_Origin[1];

      Gxyz[Gc_AN][2] =  Cell_Gxyz[Gc_AN][1]*tv[1][2]
                      + Cell_Gxyz[Gc_AN][2]*tv[2][2]
                      + Cell_Gxyz[Gc_AN][3]*tv[3][2] + Grid_Origin[2];

      Gxyz[Gc_AN][3] =  Cell_Gxyz[Gc_AN][1]*tv[1][3]
                      + Cell_Gxyz[Gc_AN][2]*tv[2][3]
                      + Cell_Gxyz[Gc_AN][3]*tv[3][3] + Grid_Origin[3];
    }

    /* increment CellOpt_iter */

    CellOpt_iter++;

  } while (CellOpt_iter<CellOpt_Maxiter); 

  /***********************************
    freeing of arrays
  ***********************************/

  free(a);
  free(b);
  free(IPIV);
}




void CellTet(char *argv[], double **CompTime)
{
  int i,j,i1,i2,j1,j2,k1,k2,I,INFO,Gc_AN,po,Imin;
  int *IPIV,LDA,LDB,N,NRHS;
  int MD_iter,CellOpt_iter,CellOpt_Maxiter;
  int myid,numprocs,ID;
  double His_Utot[300],His_Vol[300];
  double His_a[300],His_c[300];
  double *a,*b,len,len0,gf,gf0,gf1,v_step,v,v0,v1;
  double w1,w2,w3,w4,A,C,ga,gc,Emin;
  double tv0[4][4],Cell_Volume0;
  static int firsttime=1;
  static char fileE[YOUSO10] = ".ene"; 
  static char fileDRC[YOUSO10] = ".md";

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /***********************************
    allocation of arrays
  ***********************************/

  a = (double*)malloc(sizeof(double)*100);
  b = (double*)malloc(sizeof(double)*100);
  IPIV = (INTEGER*)malloc(sizeof(INTEGER)*100);

  /***********************************
    assumption: 
    a = b != c
  ***********************************/

  CellOpt_iter = 1;
  CellOpt_Maxiter = 200;

  do {
    
    MD_iter = 1;
    
    do {
    
      CompTime[myid][2] += truncation(MD_iter,1);
      CompTime[myid][3] += DFT(MD_iter,(MD_iter-1)%orbitalOpt_per_MDIter+1);
      iterout(MD_iter+MD_Current_Iter,MD_TimeStep*(MD_iter+MD_Current_Iter-1),filepath,filename);
      CompTime[myid][4] += MD_pac(MD_iter,argv[1]);
      MD_iter++;
    
    } while(MD_Opt_OK==0 && MD_iter<=MD_IterNumber);
    
    His_Utot[CellOpt_iter] = Utot;
    His_Vol[CellOpt_iter]  = Cell_Volume;
    His_a[CellOpt_iter] = sqrt(tv[1][1]*tv[1][1] + tv[1][2]*tv[1][2] + tv[1][3]*tv[1][3]);
    His_c[CellOpt_iter] = sqrt(tv[3][1]*tv[3][1] + tv[3][2]*tv[3][2] + tv[3][3]*tv[3][3]);

    if (CellOpt_iter<=8){
      
      if (firsttime){

	/* store tv into tv0 */ 
	for (i=1; i<=3; i++){
	  for (j=1; j<=3; j++){
	    tv0[i][j] = tv[i][j];
	  }
	}
	Cell_Volume0 = Cell_Volume;
	firsttime = 0;
      }

      /* update tv */ 

      switch (CellOpt_iter){

      case 1:

        /* change a&b (a=b) */

        for (i=1; i<=2; i++){
	  for (j=1; j<=3; j++){
	    tv[i][j] = 1.02*tv0[i][j];
	  }
        }

        i = 3;
	for (j=1; j<=3; j++){
	  tv[i][j] = tv0[i][j];
	}

        break;

      case 2: 

        /* change c */

        for (i=1; i<=2; i++){
	  for (j=1; j<=3; j++){
	    tv[i][j] = tv0[i][j];
	  }
        }

        i = 3;
	for (j=1; j<=3; j++){
	  tv[i][j] = 1.02*tv0[i][j];
	}

        break;

      case 3:

        /* change a, b, and c (a=b!=c) */

        for (i=1; i<=3; i++){
	  for (j=1; j<=3; j++){
	    tv[i][j] = 1.02*tv0[i][j];
	  }
        }

        break;

      case 4:

        /* change a&b (a=b) */

        for (i=1; i<=2; i++){
	  for (j=1; j<=3; j++){
	    tv[i][j] = 0.98*tv0[i][j];
	  }
        }

        i = 3;
	for (j=1; j<=3; j++){
	  tv[i][j] = tv0[i][j];
	}

        break;

      case 5: 

        /* change c */

        for (i=1; i<=2; i++){
	  for (j=1; j<=3; j++){
	    tv[i][j] = tv0[i][j];
	  }
        }

        i = 3;
	for (j=1; j<=3; j++){
	  tv[i][j] = 0.98*tv0[i][j];
	}

        break;

      case 6:

        /* change a, b, and c (a=b!=c) */

        for (i=1; i<=3; i++){
	  for (j=1; j<=3; j++){
	    tv[i][j] = 0.98*tv0[i][j];
	  }
        }

        break;

      case 7:

        /* change a, b, and c (a=b!=c) */

        for (i=1; i<=2; i++){
	  for (j=1; j<=3; j++){
	    tv[i][j] = 1.02*tv0[i][j];
	  }
        }

        for (i=3; i<=3; i++){
	  for (j=1; j<=3; j++){
	    tv[i][j] = 0.98*tv0[i][j];
	  }
        }

        break;

      case 8:

        /* change a, b, and c (a=b!=c) */

        for (i=1; i<=2; i++){
	  for (j=1; j<=3; j++){
	    tv[i][j] = 0.98*tv0[i][j];
	  }
        }

        for (i=3; i<=3; i++){
	  for (j=1; j<=3; j++){
	    tv[i][j] = 1.02*tv0[i][j];
	  }
        }

        break;

      }

      if (myid==Host_ID){
	printf("CellTet %2d %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f\n",
	       CellOpt_iter,
	       His_a[CellOpt_iter],His_c[CellOpt_iter],
               0.0,0.0,
	       His_Vol[CellOpt_iter],His_Utot[CellOpt_iter]);
      }

    } /* if (CellOpt_iter<=8) */

    /* if (8<CellOpt_iter) */

    else{

      /***********************************************************
         update 'a' and 'c' using a polynomial function 
         including up to second-order terms
      ***********************************************************/

      /* making the matrix 'a' */

      for (i1=0; i1<3; i1++){
        for (j1=0; j1<3; j1++){

          k1 = 3*i1 + j1;

	  for (i2=0; i2<3; i2++){
	    for (j2=0; j2<3; j2++){

              k2 = 3*i2 + j2;

	      a[9*k2+k1] = 0.0; 

	      for (I=CellOpt_iter-7; I<=CellOpt_iter; I++){

                w1 = pow(His_a[I],(double)i1);
                w2 = pow(His_c[I],(double)j1);

                w3 = pow(His_a[I],(double)i2);
                w4 = pow(His_c[I],(double)j2);

		a[9*k2+k1] += w1*w2*w3*w4;
	      }      
	    }
	  }
        }
      }

      /* making the vector 'b' */

      for (i1=0; i1<3; i1++){
        for (j1=0; j1<3; j1++){

          k1 = 3*i1 + j1;

	  b[k1] = 0.0; 
	  for (I=CellOpt_iter-7; I<=CellOpt_iter; I++){
	    b[k1] += His_Utot[I]*pow(His_a[I],(double)i1)*pow(His_c[I],(double)j1);
	  }      
	}
      }

      /* solving the linear equation */

      N = 9;
      NRHS = 1;
      LDA = N;
      LDB = N;
      IPIV = (INTEGER*)malloc(sizeof(INTEGER)*N);

      F77_NAME(dgesv,DGESV)(&N, &NRHS, a, &LDA, IPIV, b, &LDB, &INFO);

      /*
      for (i=0; i<9; i++){
        printf("i=%2d b=%10.7f\n",i,b[i]);
      }      
      MPI_Finalize();
      exit(0);      
      */

      /************************************************
           find the point having the lowest energy 
      ************************************************/

      Emin = 10000.0;
      for (I=CellOpt_iter-7; I<=CellOpt_iter; I++){
        if (His_Utot[I]<Emin){
          Emin = His_Utot[I];
          Imin = I;
	}
      }

      A = His_a[Imin];
      C = His_c[Imin];

      /********************************************************
        calulate the gradient of the energy at the point (a,c)
      ********************************************************/

      ga = b[3] + b[4]*C + b[5]*C*C + 2.0*b[6]*A + 2.0*b[7]*A*C + 2.0*b[8]*A*C*C;
      gc = b[1] + 2.0*b[2]*C + b[4]*A + 2.0*b[5]*A*C + b[7]*A*A + 2.0*b[8]*A*A*C;

      /********************************************************
                          update a, c, and tv 
      ********************************************************/

      if (myid==Host_ID){
	printf("CellTet %2d %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f\n",
	       CellOpt_iter,
	       His_a[CellOpt_iter],His_c[CellOpt_iter],
               ga,gc,
	       His_Vol[CellOpt_iter],His_Utot[CellOpt_iter]);
      }

      A = A - 0.2*ga;
      C = C - 0.2*gc;

      len0 = sqrt(tv[1][1]*tv[1][1] + tv[1][2]*tv[1][2] + tv[1][3]*tv[1][3]); 
      tv[1][1] = A*tv[1][1]/len0; tv[1][2] = A*tv[1][2]/len0; tv[1][3] = A*tv[1][3]/len0;

      len0 = sqrt(tv[2][1]*tv[2][1] + tv[2][2]*tv[2][2] + tv[2][3]*tv[2][3]); 
      tv[2][1] = A*tv[2][1]/len0; tv[2][2] = A*tv[2][2]/len0; tv[2][3] = A*tv[2][3]/len0; 

      len0 = sqrt(tv[3][1]*tv[3][1] + tv[3][2]*tv[3][2] + tv[3][3]*tv[3][3]);
      tv[3][1] = C*tv[3][1]/len0; tv[3][2] = C*tv[3][2]/len0; tv[3][3] = C*tv[3][3]/len0;  
    }

    /***********************************************************
                      update cartesian coordinates
    ***********************************************************/

    for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){

      Gxyz[Gc_AN][1] =  Cell_Gxyz[Gc_AN][1]*tv[1][1]
                      + Cell_Gxyz[Gc_AN][2]*tv[2][1]
                      + Cell_Gxyz[Gc_AN][3]*tv[3][1] + Grid_Origin[1];

      Gxyz[Gc_AN][2] =  Cell_Gxyz[Gc_AN][1]*tv[1][2]
                      + Cell_Gxyz[Gc_AN][2]*tv[2][2]
                      + Cell_Gxyz[Gc_AN][3]*tv[3][2] + Grid_Origin[2];

      Gxyz[Gc_AN][3] =  Cell_Gxyz[Gc_AN][1]*tv[1][3]
                      + Cell_Gxyz[Gc_AN][2]*tv[2][3]
                      + Cell_Gxyz[Gc_AN][3]*tv[3][3] + Grid_Origin[3];
    }

    /* increment CellOpt_iter */

    CellOpt_iter++;

  } while (CellOpt_iter<=CellOpt_Maxiter); 

  /***********************************
    freeing of arrays
  ***********************************/

  free(a);
  free(b);
  free(IPIV);

}
