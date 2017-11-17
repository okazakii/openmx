/**********************************************************************
    Eigen_PReHH.c:

    Eigen_PReHH.c is a MPI parallerized subroutine to solve a standard
    eigenvalue problem with a real symmetric marix using Householder
    method and lapack's dstevx_(), dstegr_(), or dstedc_().

    Log of Eigen_PReHH.c:

       Nov/26/2004  Released by T.Ozaki

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

#ifdef noomp
#include "mimic_omp.h"
#else
#include <omp.h>
#endif

 
#define  measure_time   0


static void Eigen_Original_PReHH(MPI_Comm MPI_Curret_Comm_WD, 
			 double **ac, double *ko, int n, int EVmax, int bcast_flag);

static void Eigen_Improved_PReHH(MPI_Comm MPI_Curret_Comm_WD, 
			double **ac, double *ko, int n, int EVmax, int bcast_flag);

static void myHH( MPI_Comm MPI_Curret_Comm_WD, int numprocs, int myid, int n, double **ac, double **ad, 
		  double *uu, double *b1, 
		  MPI_Request *request_send, 
		  MPI_Request *request_recv, 
		  MPI_Status *stat_send);

static void call_dsytrd( MPI_Comm MPI_Curret_Comm_WD, int numprocs, int myid, int n, double **ac, double **ad, 
			 double *uu, double *p, double *q, double *b1, 
			 MPI_Request *request_send, 
			 MPI_Request *request_recv, 
			 MPI_Status *stat_send);


void Eigen_PReHH(MPI_Comm MPI_Curret_Comm_WD, 
                 double **ac, double *ko, int n, int EVmax, int bcast_flag)
{
  Eigen_Improved_PReHH(MPI_Curret_Comm_WD, ac, ko, n, EVmax, bcast_flag); 
}








#pragma optimization_level 1
void Eigen_Improved_PReHH(MPI_Comm MPI_Curret_Comm_WD, 
                          double **ac, double *ko, int n, int EVmax, int bcast_flag)
{
  double ABSTOL=1.0e-13;
  double *q,*p;
  double **ad,*b1,*uu,
    s1,s2,s3,ss,u1,u2,r,p1,my_r,my_s2, 
    xsum,bunbo,si,co,sum,
    a1,a2,a3,a4,a5,a6,b7,r1,r2,
    r3,x1,x2,xap,tmp1,tmp2,
    bb,bb1,ui,uj,uij,
    ss0,ss1,ss2,ss3,
    s20,s21,s22,s23;

  int jj,jj1,jj2,k,k0,k1,ii,ll,i3,i2,j2,
    i,j,i1,j1,n1,n2,ik,ks,i1s,
    jk,po1,nn,count,num,num0,num1;

  int is0;
  int *is1,*ie1;
  double av_num;
  double Stime, Etime;
  int numprocs,nump,myid,tag=999,ID,IDS,IDR,ID0,ID1;
  double Stime1, Etime1,time1,time2;
  double time1a,time1b,time1c,time1d,time1e;

  MPI_Status stat;
  MPI_Request request;

  MPI_Status *stat_send;
  MPI_Request *request_send;
  MPI_Request *request_recv;

  /* MPI */
  MPI_Comm_size(MPI_Curret_Comm_WD,&numprocs);
  MPI_Comm_rank(MPI_Curret_Comm_WD,&myid);

  stat_send = malloc(sizeof(MPI_Status)*numprocs);
  request_send = malloc(sizeof(MPI_Request)*numprocs);
  request_recv = malloc(sizeof(MPI_Request)*numprocs);
  
  /****************************************************
         find the numbers of partions for MPI
  ****************************************************/

  is1 = (int*)malloc(sizeof(int)*numprocs);
  ie1 = (int*)malloc(sizeof(int)*numprocs);

  if ( numprocs<=EVmax ){

    av_num = (double)EVmax/(double)numprocs;

    for (ID=0; ID<numprocs; ID++){
      is1[ID] = (int)(av_num*(double)ID) + 1; 
      ie1[ID] = (int)(av_num*(double)(ID+1)); 
    }

    is1[0] = 1;
    ie1[numprocs-1] = EVmax; 
  }

  else {
    for (ID=0; ID<EVmax; ID++){
      is1[ID] = ID + 1; 
      ie1[ID] = ID + 1;
    }
    for (ID=EVmax; ID<numprocs; ID++){
      is1[ID] =  1;
      ie1[ID] = -2;
    }
  }

  /****************************************************
    allocation of arrays:
  ****************************************************/

  n2 = n + 5;

  ad = (double**)malloc(sizeof(double*)*n2);
  for (i=0; i<n2; i++){
    ad[i] = (double*)malloc(sizeof(double)*n2);
  }

  b1 = (double*)malloc(sizeof(double)*n2);
  uu = (double*)malloc(sizeof(double)*n2);
  q = (double*)malloc(sizeof(double)*n2);
  p = (double*)malloc(sizeof(double)*2*n2);

  for (i=1; i<=(n+2); i++){
    uu[i] = 0.0;
    b1[i] = 0.0;
  }

  if (measure_time==1) printf("size n=%3d EVmax=%2d\n",n,EVmax);

  time1 = 0.0;
  time2 = 0.0;

  /****************************************************
                   Householder method
  ****************************************************/

  myHH(MPI_Curret_Comm_WD, numprocs, myid, n, ac, ad, uu, b1, request_send, request_recv, stat_send);

  /*
  call_dsytrd(MPI_Curret_Comm_WD, numprocs, myid, n, ac, ad, uu, p, q, b1, request_send, request_recv, stat_send);
  */

  /****************************************************
                   call a lapack routine
  ****************************************************/

  if (measure_time==1) dtime(&Stime);

  for (i=1; i<=n; i++){
    p[i-1] = ad[i][i];
    q[i-1] = ad[i][i+1];
  }

  if      (dste_flag==0) lapack_dstegr1(n,EVmax,p,q,ko,ac);
  else if (dste_flag==1) lapack_dstedc1(n,p,q,ko,ac);
  else if (dste_flag==2) lapack_dstevx1(n,EVmax,p,q,ko,ac);
  else if (dste_flag==3) lapack_dsteqr1(n,p,q,ko,ac);

  if (measure_time==1){
    dtime(&Etime);
    printf("T2 myid=%2d   %15.12f\n",myid,Etime-Stime); 
  }
 
  /****************************************************
    transformation of eigenvectors to original space
  ****************************************************/

  if (measure_time==1) dtime(&Stime);

  for (i=2; i<=(n-1); i++){
    ad[i-1][i] = b1[i-1];
  }

  for (k=is1[myid]; k<=(ie1[myid]-3); k+=4){
    for (nn=2; nn<=(n-1); nn++){

      if ( (1.0e-3*ABSTOL)<fabs(uu[n-nn])){

	ss0 = 0.0;
	ss1 = 0.0;
	ss2 = 0.0;
	ss3 = 0.0;

	for (i=n-nn+1; i<=n; i++){
	  ss0 += ad[n-nn][i] * ac[k+0][i];
	  ss1 += ad[n-nn][i] * ac[k+1][i];
	  ss2 += ad[n-nn][i] * ac[k+2][i];
	  ss3 += ad[n-nn][i] * ac[k+3][i];
	}

	ss0 = 2.0*ss0/uu[n-nn];
	ss1 = 2.0*ss1/uu[n-nn];
	ss2 = 2.0*ss2/uu[n-nn];
	ss3 = 2.0*ss3/uu[n-nn];

	for (i=n-nn+1; i<=n; i++){
	  ac[k+0][i] -= ss0 * ad[n-nn][i];
	  ac[k+1][i] -= ss1 * ad[n-nn][i];
	  ac[k+2][i] -= ss2 * ad[n-nn][i];
	  ac[k+3][i] -= ss3 * ad[n-nn][i];
	}

      }
    }
  }

  ks = ie1[myid] + 1 - (ie1[myid]+1-is1[myid])%4;

  for (k=ks; (is1[myid]<=k && k<=ie1[myid]); k++){
    for (nn=2; nn<=(n-1); nn++){

      if ( (1.0e-3*ABSTOL)<fabs(uu[n-nn])){

	ss = 0.0;
	for (i=n-nn+1; i<=n; i++){
	  ss += ad[n-nn][i] * ac[k][i];
	}

	ss = 2.0*ss/uu[n-nn];
	for (i=n-nn+1; i<=n; i++){
	  ac[k][i] -= ss * ad[n-nn][i];
	}

      }
    }
  }

  if (measure_time==1){
    dtime(&Etime);
    printf("T3 myid=%2d   %15.12f\n",myid,Etime-Stime);
  }

  /****************************************************
                     normalization
  ****************************************************/

  if (measure_time==1) dtime(&Stime);

  for (j=is1[myid]; j<=ie1[myid]; j++){
    sum = 0.0;
    for (i=1; i<=n; i++){
      sum = sum + ac[j][i] * ac[j][i];
    }
    sum = 1.0/sqrt(sum);
    for (i=1; i<=n; i++){
      ac[j][i] = ac[j][i] * sum;
    }
  }

  if (measure_time==1){
    dtime(&Etime);
    printf("T4 myid=%2d   %15.12f\n",myid,Etime-Stime);
  }

  /******************************************************
   MPI:

    broadcast the full ac matrix in Host_ID
  ******************************************************/

  if (bcast_flag==1){

    if (measure_time==1) dtime(&Stime);

    BroadCast_ReMatrix(MPI_Curret_Comm_WD,ac,n,is1,ie1,myid,numprocs,
                       stat_send,request_send,request_recv);

    if (measure_time==1){
      dtime(&Etime);
      printf("T5 myid=%2d   %15.12f\n",myid,Etime-Stime);
    }
  }

  /****************************************************
            Eigenvectors to the "ac" array

             ac is distributed by column
             row:     1 to n
             column:  1 to MaxN
  ****************************************************/

  for (i=1; i<=n; i++){
    for (j=(i+1); j<=n; j++){
      tmp1 = ac[i][j];
      tmp2 = ac[j][i];
      ac[i][j] = tmp2;
      ac[j][i] = tmp1;
    }
  }

  /****************************************************
                  freeing of arrays:
  ****************************************************/

  free(stat_send);
  free(request_send);
  free(request_recv);

  free(is1);
  free(ie1);

  for (i=0; i<n2; i++){
    free(ad[i]);
  }
  free(ad);

  free(b1);
  free(uu);
  free(q);
  free(p);
}



void call_dsytrd( MPI_Comm MPI_Curret_Comm_WD, int numprocs, int myid, int n, double **ac, double **ad, 
		  double *uu, double *p, double *q, double *b1, 
		  MPI_Request *request_send, 
		  MPI_Request *request_recv, 
		  MPI_Status *stat_send)
{
  char  *UPLO="U";
  int i,j,lwork,info;
  double *A,*work;

  lwork = 2*n;

  A = (double*)malloc(sizeof(double)*n*n) ;
  work = (double*)malloc(sizeof(double)*lwork) ;

  for (i=1; i<=n; i++){
    for (j=1; j<=n; j++){
       A[(j-1)*n+(i-1)]= ac[i][j];
    }
  }
       
  dsytrd_( UPLO, &n, A, &n, p, q, b1, work, &lwork, &info);

  for (j=0; j<n; j++) {
    for (i=0; i<=j; i++) {
      ad[i+1][j+1] = A[j*n+i];
      ad[j+1][i+1] = A[j*n+i];
    }
  }



}




void myHH( MPI_Comm MPI_Curret_Comm_WD, int numprocs, int myid, int n, double **ac, double **ad, 
           double *uu, double *b1, 
           MPI_Request *request_send, 
           MPI_Request *request_recv, 
           MPI_Status *stat_send)
{
  double ABSTOL=1.0e-13; 
  int jj1,i,j,i1,j1,i1s,is0,ID1,ID0,k1,i2,n2;
  int nump,num0,k0,num,count,ID;
  int IDS,IDR,num1,tag=999;
  int *is2,*ie2,*row_flag;
  double p0,p1,p2,p3;
  double *p,*q,*u;
  double u1,u2,ss,my_r,r;
  double s1,s2,s3,s20,s21,s22,s23,tmp1,tmp2;
  double time1a,time1b,time1c,time1d,time1e;
  double Stime1, Etime1,av_num;
  double Stime, Etime;
  int OMPID,Nthrds; 

  MPI_Status stat;
  MPI_Request request;

  n2 = n + 5;

  /* allocation of arrays */
  
  u = (double*)malloc(sizeof(double)*n2);
  is2 = (int*)malloc(sizeof(int)*numprocs);
  ie2 = (int*)malloc(sizeof(int)*numprocs);
  row_flag = (int*)malloc(sizeof(int)*n2);
  q = (double*)malloc(sizeof(double)*n2);
  p = (double*)malloc(sizeof(double)*2*n2);

  time1a = 0.0;
  time1b = 0.0;
  time1c = 0.0;
  time1d = 0.0;
  time1e = 0.0;

  if (measure_time==1) dtime(&Stime);

  /****************************************************
                   Householder method
  ****************************************************/

  /* make is2 and ie2 */

  if ( numprocs<=(n-1) ){

    av_num = (double)(n-1)/(double)numprocs;

    for (ID=0; ID<numprocs; ID++){
      is2[ID] = (int)(av_num*(double)ID) + 2;
      ie2[ID] = (int)(av_num*(double)(ID+1)) + 1;
    }

    is2[0] = 2;
    ie2[numprocs-1] = n;
  }

  else{
    for (ID=0; ID<(n-1); ID++){
      is2[ID] = ID + 2;
      ie2[ID] = ID + 2;
    }
    for (ID=(n-1); ID<numprocs; ID++){
      is2[ID] = 2;
      ie2[ID] = 0;
    }
  }

  row_flag[1] = 0;
  for (ID=0; ID<numprocs; ID++){
    for (i=is2[ID]; i<=ie2[ID]; i++){
      row_flag[i] = ID;
    }
  }

  /* loop for i */

  for (i=1; i<=(n-2); i++){

    if (measure_time==1) dtime(&Stime1);

    /* make u vector */

    u[i+1] = ac[i][i+1];

    /* original version */    
    /*
    s2 = 0.0;
    for (i1=i+2; i1<=n; i1++){
      u[i1] = ac[i][i1];
      s2 += u[i1]*u[i1];
    }
    */

    /* unrolling version */    

    s20 = 0.0;
    s21 = 0.0;
    s22 = 0.0;
    s23 = 0.0;

    for (i1=i+2; i1<=(n-3); i1+=4){

      s20 += ac[i][i1+0]*ac[i][i1+0];
      s21 += ac[i][i1+1]*ac[i][i1+1];
      s22 += ac[i][i1+2]*ac[i][i1+2];
      s23 += ac[i][i1+3]*ac[i][i1+3];

      u[i1+0] = ac[i][i1+0];
      u[i1+1] = ac[i][i1+1];
      u[i1+2] = ac[i][i1+2];
      u[i1+3] = ac[i][i1+3];
    }

    i1s = n + 1 - (n+1-(i+2))%4;
    for (i1=i1s; ((i+2)<=i1 && i1<=n); i1++){
      s20 += ac[i][i1]*ac[i][i1];
      u[i1] = ac[i][i1];
    }

    s2 = s20 + s21 + s22 + s23;

    u[i] = s2;
      
    ID = row_flag[i];
    count = n - i + 1;

    if (myid==ID){
      for (IDS=ID+1; IDS<numprocs; IDS++){
        MPI_Isend(&u[i], count, MPI_DOUBLE, IDS, tag, MPI_Curret_Comm_WD, &request_send[IDS-(ID+1)]);
      }

      num = (numprocs-1) - (ID+1) + 1;  
      if (1<=num){
        MPI_Waitall(num,request_send,stat_send);
      }
    }
    else if ( (ID+1)<=myid ) {
      MPI_Recv(&u[i], count, MPI_DOUBLE, ID, tag, MPI_Curret_Comm_WD, &stat);
    }

    s2 = u[i];

    if (measure_time==1){ 
      dtime(&Etime1);
      time1a += Etime1 - Stime1;      
    }

    s1 = u[i+1]*u[i+1];
    s3 = fabs(s1 + s2);

    if (ABSTOL<fabs(u[i+1])){
      if (u[i+1]<0.0)    s3 =  sqrt(s3);
      else               s3 = -sqrt(s3);
    }
    else{
      s3 = sqrt(s3);
    }

    if (ABSTOL<fabs(s2)){

      if (measure_time==1) dtime(&Stime1);

      ss = u[i+1];
      ac[i+1][i] = s3;
      ac[i][i+1] = s3;
      u[i+1] = u[i+1] - s3;
      u1 = s3 * s3 - ss * s3;
      u2 = 2.0 * u1;

      uu[i] = u2;
      b1[i] = ss - s3;

      my_r = 0.0;
      ID0 = row_flag[i+1];

      if ( (i+1)<=is2[myid] || (i+1)<=ie2[myid] ){

	if (is2[myid]<(i+1)) is0 = i + 1;
	else                 is0 = is2[myid];

	for (i1=is0; i1<=ie2[myid]; i1++){
 
          i2 = i1 + (myid - ID0);

          p0 = 0.0;
          p1 = 0.0;
          p2 = 0.0;
          p3 = 0.0;

	  for (j=i+1; j<=(n-3); j+=4){
	    p0 += ac[i1][j  ] * u[j  ];
	    p1 += ac[i1][j+1] * u[j+1];
	    p2 += ac[i1][j+2] * u[j+2];
	    p3 += ac[i1][j+3] * u[j+3];
	  }

	  for ( ; j<=n; j++) p0 += ac[i1][j] * u[j];
	  p[i2] = (p0+p1+p2+p3) / u1;

	  /*
          p[i2] = 0.0;
	  for (j=i+1; j<=n; j++){
	    p[i2] += ac[i1][j] * u[j];
	  }
	  p[i2] = p[i2] / u1;
	  */

	  my_r += u[i1] * p[i2];
	}

        p[ie2[myid]+(myid-ID0)+1] = my_r; 
      }

      if (measure_time==1){
	dtime(&Etime1);
	time1b += Etime1 - Stime1;      
      }
      
      if (measure_time==1) dtime(&Stime1);

      /* MPI: p */

      if (ID0<=myid){ 

        /* the number of working processors */

	nump = numprocs - ID0;

        /* sending */

	num0 = ie2[myid] - is2[myid] + 2;
	k0 = is2[myid] + (myid - ID0);

	for (ID=0; ID<nump; ID++){
	  IDS = (myid + ID) % nump + ID0; 
	  MPI_Isend(&p[k0], num0, MPI_DOUBLE, IDS, tag, MPI_Curret_Comm_WD, &request_send[IDS-ID0]);
	}

        /* receiving */

	for (ID=0; ID<nump; ID++){

	  IDR = (myid - ID + nump) % nump + ID0;
	  num1 = ie2[IDR] - is2[IDR] + 2;
	  k1 = is2[IDR] + (IDR - ID0);
	  MPI_Irecv(&p[k1], num1, MPI_DOUBLE, IDR, tag, MPI_Curret_Comm_WD, &request_recv[IDR-ID0]);
	}

        /* waitall */

	MPI_Waitall(nump,request_recv,stat_send);
	MPI_Waitall(nump,request_send,stat_send);
      }

      if (measure_time==1){
	dtime(&Etime1);
	time1c += Etime1-Stime1;      
      }

      if (measure_time==1) dtime(&Stime1);

      /* calculate r */

      r = 0.0;
      for (ID=ID0; ID<numprocs; ID++){
        r += p[ie2[ID]+(ID-ID0)+1];
      }

      /* shift p */
      
      for (ID=(ID0+1); ID<numprocs; ID++){
        ID1 = ID - ID0;
	for (i1=is2[ID]; i1<=ie2[ID]; i1++){
          p[i1] = p[i1+ID1];
	}
      }

      /* calculate q */

      if ( (i+1)<=is2[myid] || (i+1)<=ie2[myid] ){
	r = r / u2;
	for (i1=i+1; i1<=n; i1++){
	  q[i1] = p[i1] - r * u[i1];
	}
      }
      
      if ( (i+1)<=is2[myid] || (i+1)<=ie2[myid] ){
        
	if (is2[myid]<(i+1)) is0 = i + 1;
	else                 is0 = is2[myid];
        
	/*
	for (i1=is0; i1<=ie2[myid]; i1++){
          
	  tmp1 = u[i1];
	  tmp2 = q[i1];

	  for (j1=i+1; j1<=n; j1++){
	    ac[i1][j1] -= tmp1 * q[j1] + tmp2 * u[j1];
	  }
	}
	*/

	for (i1=is0; i1<=ie2[myid]; i1++){
          
	  tmp1 = u[i1];
	  tmp2 = q[i1];

	  for (j1=i+1; j1<=(n-3); j1+=4){
	    ac[i1][j1  ] -= tmp1 * q[j1  ] + tmp2 * u[j1  ];
	    ac[i1][j1+1] -= tmp1 * q[j1+1] + tmp2 * u[j1+1];
	    ac[i1][j1+2] -= tmp1 * q[j1+2] + tmp2 * u[j1+2];
	    ac[i1][j1+3] -= tmp1 * q[j1+3] + tmp2 * u[j1+3];
	  }

	  for ( ; j1<=n; j1++){
	    ac[i1][j1] -= tmp1 * q[j1] + tmp2 * u[j1];
	  }
	}


      }

      if (measure_time==1){
	dtime(&Etime1);
	time1d += Etime1 - Stime1;      
      }

    }

    if (myid<row_flag[i+1]) break;
  }

  MPI_Barrier(MPI_Curret_Comm_WD);

  if (measure_time==1) dtime(&Stime1);

  /* broadcast uu, b1, and ac */

  MPI_Bcast(&uu[0],n+1,MPI_DOUBLE,numprocs-1,MPI_Curret_Comm_WD);
  MPI_Bcast(&b1[0],n+1,MPI_DOUBLE,numprocs-1,MPI_Curret_Comm_WD);

  BroadCast_ReMatrix(MPI_Curret_Comm_WD,ac,n,is2,ie2,myid,numprocs,
                     stat_send,request_send,request_recv);

  if (measure_time==1){
    dtime(&Etime1);
    time1e += Etime1-Stime1;      
  }

  for (i=1; i<=n; i++){
    for (j=1; j<=n; j++){
      ad[i][j] = ac[i][j];
    }
  }

  if (measure_time==1){
    printf("T0  myid=%2d time1a=%15.12f\n",myid,time1a);
    printf("T0  myid=%2d time1b=%15.12f\n",myid,time1b);
    printf("T0  myid=%2d time1c=%15.12f\n",myid,time1c);
    printf("T0  myid=%2d time1d=%15.12f\n",myid,time1d);
    printf("T0  myid=%2d time1e=%15.12f\n",myid,time1e);
  }
  
  if (measure_time==1){
    dtime(&Etime);
    printf("T1 myid=%2d   %15.12f\n",myid,Etime-Stime);
  }

  /* freeing of arrays */

  free(row_flag);
  free(is2);
  free(ie2);
  free(u);
  free(q);
  free(p);

}








void Eigen_Original_PReHH(MPI_Comm MPI_Curret_Comm_WD, 
                  double **ac, double *ko, int n, int EVmax, int bcast_flag)
{
  double ABSTOL=1.0e-13;
  double **ad,*D,*E,
    *b1,*u,*uu,
    *p,*q,
    s1,s2,s3,ss,u1,u2,r,p1,my_r, 
    xsum,bunbo,si,co,sum,
    a1,a2,a3,a4,a5,a6,b7,r1,r2,
    r3,x1,x2,xap,tmp1,tmp2,
    bb,bb1,ui,uj,uij;

  int jj,jj1,jj2,k,ii,ll,i3,i2,j2,
    i,j,i1,j1,n1,n2,ik,
    jk,po1,nn,count,num;

  int *is1,*ie1,*is2,*ie2;
  double av_num,*ac1;
  double Stime, Etime;
  int numprocs,myid,tag=0,ID;
  double Stime1, Etime1,time1,time2,time1a,time1b,time1c,time1d;

  MPI_Status stat;
  MPI_Request request;

  MPI_Status *stat_send;
  MPI_Request *request_send;
  MPI_Request *request_recv;

  /* MPI */
  MPI_Comm_size(MPI_Curret_Comm_WD,&numprocs);
  MPI_Comm_rank(MPI_Curret_Comm_WD,&myid);

  stat_send = malloc(sizeof(MPI_Status)*numprocs);
  request_send = malloc(sizeof(MPI_Request)*numprocs);
  request_recv = malloc(sizeof(MPI_Request)*numprocs);

  /****************************************************
         find the numbers of partions for MPI
  ****************************************************/

  is1 = (int*)malloc(sizeof(int)*numprocs);
  ie1 = (int*)malloc(sizeof(int)*numprocs);

  if ( numprocs<=EVmax ){

    av_num = (double)EVmax/(double)numprocs;

    for (ID=0; ID<numprocs; ID++){
      is1[ID] = (int)(av_num*(double)ID) + 1; 
      ie1[ID] = (int)(av_num*(double)(ID+1)); 
    }

    is1[0] = 1;
    ie1[numprocs-1] = EVmax; 
  }

  else {
    for (ID=0; ID<EVmax; ID++){
      is1[ID] = ID + 1; 
      ie1[ID] = ID + 1;
    }
    for (ID=EVmax; ID<numprocs; ID++){
      is1[ID] =  1;
      ie1[ID] = -2;
    }
  }

  /*
  for (ID=0; ID<numprocs; ID++){
    printf("myid=%2d ID=%2d is1=%3d ie1=%3d\n",myid,ID,is1[ID],ie1[ID]);
  }
  */

  /****************************************************
    allocation of arrays:
  ****************************************************/

  is2 = (int*)malloc(sizeof(int)*numprocs);
  ie2 = (int*)malloc(sizeof(int)*numprocs);

  n2 = n + 5;

  ad = (double**)malloc(sizeof(double*)*n2);
  for (i=0; i<n2; i++){
    ad[i] = (double*)malloc(sizeof(double)*n2);
  }

  b1 = (double*)malloc(sizeof(double)*n2);
  u = (double*)malloc(sizeof(double)*n2);
  uu = (double*)malloc(sizeof(double)*n2);
  p = (double*)malloc(sizeof(double)*n2);
  q = (double*)malloc(sizeof(double)*n2);

  D = (double*)malloc(sizeof(double)*n2);
  E = (double*)malloc(sizeof(double)*n2);

  for (i=1; i<=(n+2); i++){
    uu[i] = 0.0;
  }

  if (measure_time==1) printf("size n=%3d EVmax=%2d\n",n,EVmax);
  if (measure_time==1) dtime(&Stime);


  time1 = 0.0;
  time2 = 0.0;
  time1a = 0.0;
  time1b = 0.0;
  time1c = 0.0;
  time1d = 0.0;

  /****************************************************
                   Householder method
  ****************************************************/

  for (i=1; i<=(n-2); i++){

    /* make u vector */

    s1 = ac[i+1][i] * ac[i+1][i];
    s2 = 0.0;
    u[i+1] = ac[i+1][i];
    for (i1=i+2; i1<=n; i1++){
      tmp1 = ac[i1][i]; 
      s2 += tmp1*tmp1;
      u[i1] = tmp1;
    }
    s3 = fabs(s1 + s2);

    if (ABSTOL<fabs(ac[i+1][i])){
      if (ac[i+1][i]<0.0)    s3 =  sqrt(s3);
      else                   s3 = -sqrt(s3);
    }
    else{
      s3 = sqrt(s3);
    }

    if (ABSTOL<fabs(s2)){

      ss = ac[i+1][i];
      ac[i+1][i] = s3;
      ac[i][i+1] = s3;
      u[i+1] = u[i+1] - s3;
      u1 = s3 * s3 - ss * s3;
      u2 = 2.0 * u1;
      uu[i] = u2;
      b1[i] = ss - s3;

      if (measure_time==1) dtime(&Stime1);

      r = 0.0;
      for (i1=i+1; i1<=n; i1++){
	p1 = 0.0;
	for (j=i+1; j<=n; j++){
	  p1 += ac[i1][j] * u[j];
	}
	p[i1] = p1 / u1;
	r += u[i1] * p[i1];
      }

      if (measure_time==1){
        dtime(&Etime1);
        time1a += Etime1 - Stime1;
      }

      if (measure_time==1) dtime(&Stime1);

      printf("i=%2d r=%18.15f\n",i,r);

      r = r / u2;
      for (i1=i+1; i1<=n; i1++){
	q[i1] = p[i1] - r * u[i1];
      }

      if (measure_time==1){
        dtime(&Etime1);
        time1c += Etime1-Stime1;      
      }

      if (measure_time==1) dtime(&Stime1);

      for (i1=i+1; i1<=n; i1++){

        tmp1 = u[i1];
        tmp2 = q[i1];

	for (j1=i+1; j1<=n; j1++){
	  ac[i1][j1] -= tmp1 * q[j1] + tmp2 * u[j1];
	}
      }

      if (measure_time==1){
        dtime(&Etime1);
        time2 += Etime1-Stime1;
      }

    }
  }

  for (i=1; i<=n; i++){
    for (j=1; j<=n; j++){
      ad[i][j] = ac[i][j];
    }
  }

  if (measure_time==1){
    printf("T0  myid=%2d time1 =%15.12f time2 =%15.12f\n",myid,time1,time2);
    printf("T0a myid=%2d time1a=%15.12f time1b=%15.12f time1c=%15.12f time1d=%15.12f\n",
            myid,time1a,time1b,time1c,time1d);
  }

  if (measure_time==1){
    dtime(&Etime);
    printf("T1 myid=%2d   %15.12f\n",myid,Etime-Stime);
  }

  /****************************************************
                   call a lapack routine
  ****************************************************/

  if (measure_time==1) dtime(&Stime);

  for (i=1; i<=n; i++){
    D[i-1] = ad[i][i];
    E[i-1] = ad[i][i+1];
  }

  if      (dste_flag==0) lapack_dstegr1(n,EVmax,D,E,ko,ac);
  else if (dste_flag==1) lapack_dstedc1(n,D,E,ko,ac);
  else if (dste_flag==2) lapack_dstevx1(n,EVmax,D,E,ko,ac);
  else if (dste_flag==3) lapack_dsteqr1(n,D,E,ko,ac);

  if (measure_time==1){
    dtime(&Etime);
    printf("T2 myid=%2d   %15.12f\n",myid,Etime-Stime); 
  }

  /****************************************************
    transformation of eigenvectors to original space
  ****************************************************/

  if (measure_time==1) dtime(&Stime);

  for (i=2; i<=n-1; i++){
    ad[i-1][i] = b1[i-1];
  }

  for (k=is1[myid]; k<=ie1[myid]; k++){
    for (nn=2; nn<=n-1; nn++){

      if ( (1.0e-3*ABSTOL)<fabs(uu[n-nn])){
	ss = 0.0;
	for (i=n-nn+1; i<=n; i++){
	  ss += ad[n-nn][i] * ac[k][i];
	}
	ss = 2.0*ss/uu[n-nn];
	for (i=n-nn+1; i<=n; i++){
	  ac[k][i] -= ss * ad[n-nn][i];
	}
      }
    }
  }

  if (measure_time==1){
    dtime(&Etime);
    printf("T3 myid=%2d   %15.12f\n",myid,Etime-Stime);
  }

  /****************************************************
                     normalization
  ****************************************************/

  if (measure_time==1) dtime(&Stime);

  for (j=is1[myid]; j<=ie1[myid]; j++){
    sum = 0.0;
    for (i=1; i<=n; i++){
      sum = sum + ac[j][i] * ac[j][i];
    }
    sum = 1.0/sqrt(sum);
    for (i=1; i<=n; i++){
      ac[j][i] = ac[j][i] * sum;
    }
  }

  if (measure_time==1){
    dtime(&Etime);
    printf("T4 myid=%2d   %15.12f\n",myid,Etime-Stime);
  }

  /******************************************************
   MPI:

    broadcast the full ac matrix in Host_ID
  ******************************************************/

  if (bcast_flag==1){

    if (measure_time==1) dtime(&Stime);

    BroadCast_ReMatrix(MPI_Curret_Comm_WD,ac,n,is1,ie1,myid,numprocs,
                       stat_send,request_send,request_recv);

    if (measure_time==1){
      dtime(&Etime);
      printf("T5 myid=%2d   %15.12f\n",myid,Etime-Stime);
    }
  }
  
  /****************************************************
            Eigenvectors to the "ac" array

             ac is distributed by column
             row:     1 to n
             column:  1 to MaxN
  ****************************************************/

  for (i=1; i<=n; i++){
    for (j=(i+1); j<=n; j++){
      tmp1 = ac[i][j];
      tmp2 = ac[j][i];
      ac[i][j] = tmp2;
      ac[j][i] = tmp1;
    }
  }

  /****************************************************
                  freeing of arrays:
  ****************************************************/

  free(stat_send);
  free(request_send);
  free(request_recv);

  free(is1);
  free(ie1);
  free(is2);
  free(ie2);

  for (i=0; i<n2; i++){
    free(ad[i]);
  }
  free(ad);

  free(b1);
  free(u);
  free(uu);
  free(p);
  free(q);
  free(D);
  free(E);
}
