#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>

#define asize1   40
#define asize2   40
#define asize3   40

#define asize4   100000000

typedef struct { double r,i; } dcomplex;

void myfunc1();
void myfunc2();
void myfunc3();

void dataset();

double *array1;

int main(int argc, char *argv[]) 
{
  int i,j,k;
  int myid,numprocs; 
  int knum_i,knum_j,knum_k;
  int T_knum,num_kloop0,S_knum,E_knum;
  double *T_KGrids1,*T_KGrids2,*T_KGrids3;
  int *Ti_KGrids1,*Tj_KGrids2,*Tk_KGrids3,*arpo;
  double *KGrids1, *KGrids2, *KGrids3;
  double tmp,k1,k2,k3,snum_i,snum_j,snum_k;


  

  knum_i = 16;
  knum_j = 16;
  knum_k = 16;
  
  numprocs = 16;



  T_knum = 0;
  for (i=0; i<=(knum_i-1); i++){
    for (j=0; j<=(knum_j-1); j++){
      for (k=0; k<=(knum_k-1); k++){
        T_knum++;
      }
    }
  }

  T_KGrids1 = (double*)malloc(sizeof(double)*T_knum);
  T_KGrids2 = (double*)malloc(sizeof(double)*T_knum);
  T_KGrids3 = (double*)malloc(sizeof(double)*T_knum);

  Ti_KGrids1 = (int*)malloc(sizeof(int)*T_knum);
  Tj_KGrids2 = (int*)malloc(sizeof(int)*T_knum);
  Tk_KGrids3 = (int*)malloc(sizeof(int)*T_knum);

  KGrids1 = (double*)malloc(sizeof(double)*knum_i);
  KGrids2 = (double*)malloc(sizeof(double)*knum_j);
  KGrids3 = (double*)malloc(sizeof(double)*knum_k);



  for (myid=0; myid<numprocs; myid++){

    snum_i = knum_i;
    snum_j = knum_j;
    snum_k = knum_k;

    /* set up  k-grids */
    for (i=0; i<=(knum_i-1); i++){
      if (knum_i==1){
	k1 = 0.0;
      }
      else {
	k1 = -0.49999999999999 + (2.0*(double)i+1.0)/(2.0*snum_i);
      }
      KGrids1[i]=k1;
    }
    for (i=0; i<=(knum_j-1); i++){
      if (knum_j==1){
	k1 = 0.0;
      }
      else {
	k1 = -0.49999999999998 + (2.0*(double)i+1.0)/(2.0*snum_j);
      }
      KGrids2[i]=k1;
    }
    for (i=0; i<=(knum_k-1); i++){
      if (knum_k==1){
	k1 = 0.0;
      }
      else {
	k1 = -0.49999999999997 + (2.0*(double)i+1.0)/(2.0*snum_k);
      }
      KGrids3[i]=k1;
    }



    T_knum = 0;
    for (i=0; i<=(knum_i-1); i++){
      for (j=0; j<=(knum_j-1); j++){
	for (k=0; k<=(knum_k-1); k++){

	  T_KGrids1[T_knum] = KGrids1[i];
	  T_KGrids2[T_knum] = KGrids2[j];
	  T_KGrids3[T_knum] = KGrids3[k];

	  Ti_KGrids1[T_knum] = i;
	  Tj_KGrids2[T_knum] = j;
	  Tk_KGrids3[T_knum] = k;

	  T_knum++;
	}
      }
    }



    if (T_knum<=myid){
      S_knum = -10;
      E_knum = -100;
      num_kloop0 = 1;
    }
    else if (T_knum<numprocs) {
      S_knum = myid;
      E_knum = myid;
      num_kloop0 = 1;
    }
    else {
      tmp = (double)T_knum/(double)numprocs; 
      num_kloop0 = (int)tmp + 1;

      S_knum = (int)((double)myid*(tmp+0.0001)); 
      E_knum = (int)((double)(myid+1)*(tmp+0.0001)) - 1;
      if (myid==(numprocs-1)) E_knum = T_knum - 1;
      if (E_knum<0)           E_knum = 0;
    }


    printf("myid=%2d num_kloop0=%2d S_knum=%2d E_knum=%2d  %2d\n",myid,num_kloop0,S_knum,E_knum,E_knum-S_knum+1);

  }

}

void dataset()
{
  static int n;
  double *array0 = array1;

  for (n=0; n<asize4; n++){
    array0[n] = 0.0;
  }

}



void myfunc1()
{
  static int i,j,k,m,n;
  static double ***a;
  static double ***b;
  static double ***c;

  /* allocation */

  a = (double***)malloc(sizeof(double**)*asize1);
  for (i=0; i<asize1; i++){
    a[i] = (double**)malloc(sizeof(double*)*asize2);
    for (j=0; j<asize2; j++){
      a[i][j] = (double*)malloc(sizeof(double)*asize3);
    }
  }

  b = (double***)malloc(sizeof(double**)*asize2);
  for (i=0; i<asize2; i++){
    b[i] = (double**)malloc(sizeof(double*)*asize2);
    for (j=0; j<asize2; j++){
      b[i][j] = (double*)malloc(sizeof(double)*asize2);
    }
  }

  c = (double***)malloc(sizeof(double**)*asize3);
  for (i=0; i<asize3; i++){
    c[i] = (double**)malloc(sizeof(double*)*asize3);
    for (j=0; j<asize3; j++){
      c[i][j] = (double*)malloc(sizeof(double)*asize3);
    }
  }

  /* someting */

  for (i=0; i<asize1; i++){
    for (j=0; j<asize2; j++){
      for (k=0; k<asize3; k++){
	a[i][j][k] = cos((double)(i+j+k));
	a[i][j][k] = cos((double)(2*i+j+k));
	a[i][j][k] = cos((double)(i+2*j+k));
	a[i][j][k] = cos((double)(i+j+2*k));
	a[i][j][k] = cos((double)(i+j+3*k));
	a[i][j][k] = cos((double)(i+j+4*k));
	a[i][j][k] = cos((double)(i+j+5*k));
      }
    }
  }

  for (i=0; i<asize2; i++){
    for (j=0; j<asize2; j++){
      for (k=0; k<asize2; k++){
	b[i][j][k] = cos((double)(i+j+k))+a[i][j][k];
      }
    }
  }

  for (i=0; i<asize3; i++){
    for (j=0; j<asize3; j++){
      for (k=0; k<asize3; k++){
	c[i][j][k] = cos((double)(i+j+k))+b[i][j][k];
      }
    }
  }
 
  /* freeing */

  for (i=0; i<asize1; i++){
    for (j=0; j<asize2; j++){
      free(a[i][j]);
    }
    free(a[i]);
  }
  free(a);

  for (i=0; i<asize2; i++){
    for (j=0; j<asize2; j++){
      free(b[i][j]);
    }
    free(b[i]);
  }
  free(b);

  for (i=0; i<asize3; i++){
    for (j=0; j<asize3; j++){
      free(c[i][j]);
    }
    free(c[i]);
  }
  free(c);

}






void myfunc2()
{
  static int i,j,k,m,n;
  static double ***a;
  static double ***b;
  static double ***c;

  /* allocation */

  a = (double***)malloc(sizeof(double**)*asize1);
  for (i=0; i<asize1; i++){
    a[i] = (double**)malloc(sizeof(double*)*asize2);
    for (j=0; j<asize2; j++){
      a[i][j] = (double*)malloc(sizeof(double)*asize3);
    }
  }

  b = (double***)malloc(sizeof(double**)*asize2);
  for (i=0; i<asize2; i++){
    b[i] = (double**)malloc(sizeof(double*)*asize2);
    for (j=0; j<asize2; j++){
      b[i][j] = (double*)malloc(sizeof(double)*asize2);
    }
  }

  c = (double***)malloc(sizeof(double**)*asize3);
  for (i=0; i<asize3; i++){
    c[i] = (double**)malloc(sizeof(double*)*asize3);
    for (j=0; j<asize3; j++){
      c[i][j] = (double*)malloc(sizeof(double)*asize3);
    }
  }

  /* someting */

  for (i=0; i<asize1; i++){
    for (j=0; j<asize2; j++){
      for (k=0; k<asize3; k++){
	a[i][j][k] = cos((double)(i+j+k));
	a[i][j][k] = cos((double)(2*i+j+k));
	a[i][j][k] = cos((double)(i+2*j+k));
	a[i][j][k] = cos((double)(i+j+2*k));
	a[i][j][k] = cos((double)(i+j+3*k));
	a[i][j][k] = cos((double)(i+j+4*k));
	a[i][j][k] = cos((double)(i+j+5*k));
      }
    }
  }

  for (i=0; i<asize2; i++){
    for (j=0; j<asize2; j++){
      for (k=0; k<asize2; k++){
	b[i][j][k] = cos((double)(i+j+k))+a[i][j][k];
      }
    }
  }

  for (i=0; i<asize3; i++){
    for (j=0; j<asize3; j++){
      for (k=0; k<asize3; k++){
	c[i][j][k] = cos((double)(i+j+k))+b[i][j][k];
      }
    }
  }
 
  /* freeing */

  for (i=0; i<asize1; i++){
    for (j=0; j<asize2; j++){
      free(a[i][j]);
    }
    free(a[i]);
  }
  free(a);

  for (i=0; i<asize2; i++){
    for (j=0; j<asize2; j++){
      free(b[i][j]);
    }
    free(b[i]);
  }
  free(b);

  for (i=0; i<asize3; i++){
    for (j=0; j<asize3; j++){
      free(c[i][j]);
    }
    free(c[i]);
  }
  free(c);

}




void myfunc3()
{
  register int i,j,k,m,n;
  register double ***a;
  register double ***b;
  register double ***c;

  /* allocation */

  a = (double***)malloc(sizeof(double**)*asize1);
  for (i=0; i<asize1; i++){
    a[i] = (double**)malloc(sizeof(double*)*asize2);
    for (j=0; j<asize2; j++){
      a[i][j] = (double*)malloc(sizeof(double)*asize3);
    }
  }

  b = (double***)malloc(sizeof(double**)*asize2);
  for (i=0; i<asize2; i++){
    b[i] = (double**)malloc(sizeof(double*)*asize2);
    for (j=0; j<asize2; j++){
      b[i][j] = (double*)malloc(sizeof(double)*asize2);
    }
  }

  c = (double***)malloc(sizeof(double**)*asize3);
  for (i=0; i<asize3; i++){
    c[i] = (double**)malloc(sizeof(double*)*asize3);
    for (j=0; j<asize3; j++){
      c[i][j] = (double*)malloc(sizeof(double)*asize3);
    }
  }

  /* someting */

  for (i=0; i<asize1; i++){
    for (j=0; j<asize2; j++){
      for (k=0; k<asize3; k++){
	a[i][j][k] = cos((double)(i+j+k));
	a[i][j][k] = cos((double)(2*i+j+k));
	a[i][j][k] = cos((double)(i+2*j+k));
	a[i][j][k] = cos((double)(i+j+2*k));
	a[i][j][k] = cos((double)(i+j+3*k));
	a[i][j][k] = cos((double)(i+j+4*k));
	a[i][j][k] = cos((double)(i+j+5*k));
      }
    }
  }

  for (i=0; i<asize2; i++){
    for (j=0; j<asize2; j++){
      for (k=0; k<asize2; k++){
	b[i][j][k] = cos((double)(i+j+k))+a[i][j][k];
      }
    }
  }

  for (i=0; i<asize3; i++){
    for (j=0; j<asize3; j++){
      for (k=0; k<asize3; k++){
	c[i][j][k] = cos((double)(i+j+k))+b[i][j][k];
      }
    }
  }
 
  /* freeing */

  for (i=0; i<asize1; i++){
    for (j=0; j<asize2; j++){
      free(a[i][j]);
    }
    free(a[i]);
  }
  free(a);

  for (i=0; i<asize2; i++){
    for (j=0; j<asize2; j++){
      free(b[i][j]);
    }
    free(b[i]);
  }
  free(b);

  for (i=0; i<asize3; i++){
    for (j=0; j<asize3; j++){
      free(c[i][j]);
    }
    free(c[i]);
  }
  free(c);

}




