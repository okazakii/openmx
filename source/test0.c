#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define   asize1   1000

double rnd(double width);

int main(int argc, char *argv[]) 
{
  int i,j,k;
  double sum;
  double **a;
  double **b;
  double **c;

  a = (double**)malloc(sizeof(double*)*asize1);
  for (i=0; i<asize1; i++){
    a[i] = (double*)malloc(sizeof(double)*asize1);
  }

  b = (double**)malloc(sizeof(double*)*asize1);
  for (i=0; i<asize1; i++){
    b[i] = (double*)malloc(sizeof(double)*asize1);
  }

  c = (double**)malloc(sizeof(double*)*asize1);
  for (i=0; i<asize1; i++){
    c[i] = (double*)malloc(sizeof(double)*asize1);
  }

  for (i=0; i<asize1; i++){
    for (j=0; j<asize1; j++){
      a[i][j] = 1.0e-10*(double)(i + j);
      b[i][j] = 1.0e-10*(double)(i - j);
    }
  }



  for (i=0; i<asize1; i++){
    for (j=0; j<asize1; j++){
      sum = 0.0;
      for (k=0; k<asize1; k++){
        sum += a[i][k]*b[j][k]; 
      }
      c[i][j] = sum;
    }
  }


}


double rnd(double width)
{
  /****************************************************
       This rnd() function generates random number
                -width/2 to width/2
  ****************************************************/

  double result;

  result = rand();

  while (width<result){
    result = result/2.0;
  }
  result = result - width*0.75;
  return result;
}

