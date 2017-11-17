/**********************************************************************
  Spherical_Bessel.c:

     Spherical_Bessel.c is a subroutine to calculate the spherical 
     Bessel functions and its derivative from 0 to lmax

  Log of Spherical_Bessel.c:

     08/Nov/2005  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "openmx_common.h"

#define xmin  0.0

void Spherical_Bessel( double x, int lmax, double *sb, double *dsb ) 
{
  int m,n,nmax;
  double *tsb;
  double j0,j1,j0p,j1p,sf,tmp,si,co,ix,ix2;

  if (x<0.0){
    printf("minus x is invalid for Spherical_Bessel\n");
    exit(0);    
  } 

  /* find an appropriate nmax */

  nmax = lmax + 3*x + 20;
  if (nmax<100) nmax = 100; 

  /* allocate tsb */

  tsb  = (double*)malloc(sizeof(double)*(nmax+1)); 
  
  /* if x is larger than xmin */

  if ( xmin < x ){

    /* initial values*/

    tsb[nmax]   = 0.0;
    tsb[nmax-1] = 1.0e-14;

    /* downward recurrence from nmax-2 to lmax+2 */

    for ( n=nmax-1; (lmax+2)<n; n-- ){

      tsb[n-1] = (2.0*n + 1.0)/x*tsb[n] - tsb[n+1];

      if (1.0e+250<tsb[n-1]){
        tmp = tsb[n-1];        
        tsb[n-1] /= tmp;
        tsb[n  ] /= tmp;
      }

      /*
      printf("n=%3d tsb[n-1]=%18.15f\n",n,tsb[n-1]);
      */
    }

    /* downward recurrence from lmax+1 to 0 */

    n = lmax + 3;
    tmp = tsb[n-1];        
    tsb[n-1] /= tmp;
    tsb[n  ] /= tmp;

    for ( n=lmax+2; 0<n; n-- ){
      tsb[n-1] = (2.0*n + 1.0)/x*tsb[n] - tsb[n+1];

      if (1.0e+250<tsb[n-1]){
        tmp = tsb[n-1];
        for (m=n-1; m<=lmax+1; m++){
          tsb[m] /= tmp;
        }
      }
    }

    /* normalization */

    si = sin(x);
    co = cos(x);
    ix = 1.0/x;
    ix2 = ix*ix;
    j0 = si*ix;
    j1 = si*ix*ix - co*ix;

    if (fabs(tsb[1])<fabs(tsb[0])) sf = j0/tsb[0];
    else                           sf = j1/tsb[1];

    /* tsb to sb */

    for ( n=0; n<=lmax+1; n++ ){
      sb[n] = tsb[n]*sf;
    }

    /* derivative of sb */

    dsb[0] = co*ix - si*ix*ix;
    for ( n=1; n<=lmax; n++ ){
      dsb[n] = ( (double)n*sb[n-1] - (double)(n+1.0)*sb[n+1] )/(2.0*(double)n + 1.0);
    }

  } 

  /* if x is smaller than xmin */

  else {

    /* sb */

    for ( n=0; n<=lmax; n++ ){
      sb[n] = 0.0;
    }
    sb[0] = 1.0;

    /* derivative of sb */

    dsb[0] = 0.0;
    for ( n=1; n<=lmax; n++ ){
      dsb[n] = ( (double)n*sb[n-1] - (double)(n+1.0)*sb[n+1] )/(2.0*(double)n + 1.0);
    }
  }

  /* free tsb */

  free(tsb);
}
