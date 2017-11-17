/**********************************************************************
  TRAN_Calc_CentGreenLesser.c:

  TRAN_Calc_CentGreenLesser.c is a subroutine to calculate the lesser 
  Green's function of the central part. 

  Log of TRAN_Calc_CentGreenLesser.c:

     24/July/2008  Released by H.Kino and T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>
#include "tran_prototypes.h"
#include "lapack_prototypes.h"
#define PI              3.1415926535897932384626



void TRAN_Calc_CentGreenLesser(
                      /* input */
                      dcomplex w,
                      double ChemP_e[2],
                      int nc, 
                      int Order_Lead_Side[2],
                      dcomplex *SigmaL,
                      dcomplex *SigmaL_Ad,
                      dcomplex *SigmaR, 
                      dcomplex *SigmaR_Ad, 
                      dcomplex *GC, 
                      dcomplex *GC_Ad, 
                      dcomplex *HCCk, 
                      dcomplex *SCC, 

                      /* work, nc*nc */
                      dcomplex *v1, 
                      dcomplex *v2,
 
                      /*  output */ 
                      dcomplex *Gless 
                      )

#define GC_ref(i,j)        GC[nc*((j)-1)+(i)-1]
#define GC_Ad_ref(i,j)     GC_Ad[nc*((j)-1)+(i)-1]
#define SigmaL_ref(i,j)    SigmaL[nc*((j)-1)+(i)-1]
#define SigmaL_Ad_ref(i,j) SigmaL_Ad[nc*((j)-1)+(i)-1]
#define SigmaR_ref(i,j)    SigmaR[nc*((j)-1)+(i)-1]
#define SigmaR_Ad_ref(i,j) SigmaR_Ad[nc*((j)-1)+(i)-1]
#define SCC_ref(i,j)       SCC[nc*((j)-1)+(i)-1]
#define HCCk_ref(i,j)      HCCk[nc*((j)-1)+(i)-1]
#define v1_ref(i,j)        v1[nc*((j)-1)+(i)-1] 
#define v2_ref(i,j)        v2[nc*((j)-1)+(i)-1] 
#define Gless_ref(i,j)     Gless[nc*((j)-1)+(i)-1]

{
  int i,j;
  int side;
  dcomplex alpha,beta;
  dcomplex ctmp;

  alpha.r = 1.0;
  alpha.i = 0.0;
  beta.r  = 0.0;
  beta.i  = 0.0;

  /******************************************************
    lesser Green's function
  ******************************************************/

  /* v1 = -\sigama_{L or R}(z^*) */

  if (Order_Lead_Side[1]==0){

    for (i=1; i<=nc; i++) {
      for (j=1; j<=nc; j++) {

	v1_ref(i,j).r = SigmaL_ref(i,j).r - SigmaL_Ad_ref(i,j).r; 
	v1_ref(i,j).i = SigmaL_ref(i,j).i - SigmaL_Ad_ref(i,j).i;
      }
    }
  }
  else{

    for (i=1; i<=nc; i++) {
      for (j=1; j<=nc; j++) {

	v1_ref(i,j).r = SigmaR_ref(i,j).r - SigmaR_Ad_ref(i,j).r; 
	v1_ref(i,j).i = SigmaR_ref(i,j).i - SigmaR_Ad_ref(i,j).i;
      }
    }
  }

  /* v2 = G(z) * v1 */

  F77_NAME(zgemm,ZGEMM)("N","N", &nc, &nc, &nc, &alpha, GC, &nc, v1, &nc, &beta, v2, &nc);

  /* Gless = G(z) * v1 * G(z^*)  */

  F77_NAME(zgemm,ZGEMM)("N","N", &nc, &nc, &nc, &alpha, v2, &nc, GC_Ad, &nc, &beta, Gless, &nc);

  /******************************************************
    -1/(i 2Pi) * Gless
  ******************************************************/

  for (i=1; i<=nc; i++) {
    for (j=1; j<=nc; j++) {
      ctmp.r = Gless_ref(i,j).r/(2.0*PI);
      ctmp.i = Gless_ref(i,j).i/(2.0*PI);
      Gless_ref(i,j).r =-ctmp.i;
      Gless_ref(i,j).i = ctmp.r;
    }
  }

}










void TRAN_Calc_CentGreenLesser_old(
                      /* input */
                      dcomplex w,
                      double ChemP_e[2],
                      int nc, 
                      int Order_Lead_Side[2],
                      dcomplex *SigmaL,
                      dcomplex *SigmaL_Ad,
                      dcomplex *SigmaR, 
                      dcomplex *SigmaR_Ad, 
                      dcomplex *GC, 
                      dcomplex *GC_Ad, 
                      dcomplex *HCCk, 
                      dcomplex *SCC, 

                      /* work, nc*nc */
                      dcomplex *v1, 
                      dcomplex *v2,
 
                      /*  output */ 
                      dcomplex *Gless 
                      )

#define GC_ref(i,j)        GC[nc*((j)-1)+(i)-1]
#define GC_Ad_ref(i,j)     GC_Ad[nc*((j)-1)+(i)-1]
#define SigmaL_ref(i,j)    SigmaL[nc*((j)-1)+(i)-1]
#define SigmaL_Ad_ref(i,j) SigmaL_Ad[nc*((j)-1)+(i)-1]
#define SigmaR_ref(i,j)    SigmaR[nc*((j)-1)+(i)-1]
#define SigmaR_Ad_ref(i,j) SigmaR_Ad[nc*((j)-1)+(i)-1]
#define SCC_ref(i,j)       SCC[nc*((j)-1)+(i)-1]
#define HCCk_ref(i,j)      HCCk[nc*((j)-1)+(i)-1]
#define v1_ref(i,j)        v1[nc*((j)-1)+(i)-1] 
#define v2_ref(i,j)        v2[nc*((j)-1)+(i)-1] 
#define Gless_ref(i,j)     Gless[nc*((j)-1)+(i)-1]

{
  int i,j;
  int side;
  dcomplex alpha,beta;
  dcomplex ctmp;

  alpha.r = 1.0;
  alpha.i = 0.0;
  beta.r  = 0.0;
  beta.i  = 0.0;

  /******************************************************
    retarded Green's function of the left or right part
  ******************************************************/

  /* v1 = 1/2z^* S - 1/2H - \sigama_{L or R}(z^*)  */

  if (Order_Lead_Side[1]==0){

    for (i=1; i<=nc; i++) {
      for (j=1; j<=nc; j++) {
	v1_ref(i,j).r = 0.0*( 0.5*w.r*SCC_ref(i,j).r + 0.5*w.i*SCC_ref(i,j).i) - 0.5*HCCk_ref(i,j).r - SigmaL_Ad_ref(i,j).r; 
	v1_ref(i,j).i = 0.0*(-0.5*w.i*SCC_ref(i,j).r + 0.5*w.r*SCC_ref(i,j).i) - 0.5*HCCk_ref(i,j).i - SigmaL_Ad_ref(i,j).i;
      }
    }
  }
  else{

    for (i=1; i<=nc; i++) {
      for (j=1; j<=nc; j++) {
	v1_ref(i,j).r = 0.0*( 0.5*w.r*SCC_ref(i,j).r + 0.5*w.i*SCC_ref(i,j).i) - 0.5*HCCk_ref(i,j).r - SigmaR_Ad_ref(i,j).r; 
	v1_ref(i,j).i = 0.0*(-0.5*w.i*SCC_ref(i,j).r + 0.5*w.r*SCC_ref(i,j).i) - 0.5*HCCk_ref(i,j).i - SigmaR_Ad_ref(i,j).i;
      }
    }
  }

  /* v2 = G(z) [1/2z^* S - 1/2 H - \sigama_{L or R}(z^*)]  */

  F77_NAME(zgemm,ZGEMM)("N","N", &nc, &nc, &nc, &alpha, GC, &nc, v1, &nc, &beta, v2, &nc);

  /* Gless = G(z) [1/2z^* S - 1/2H - \sigama_{L or R}(z^*)] G(z^*)  */

  F77_NAME(zgemm,ZGEMM)("N","N", &nc, &nc, &nc, &alpha, v2, &nc, GC_Ad, &nc, &beta, Gless, &nc);

  /******************************************************
    advanced Green's function of the left or right part
  ******************************************************/

  /* v1 = 1/2z S - 1/2 H - \sigama_{L or R}(z)  */

  if (Order_Lead_Side[1]==0){
    for (i=1; i<=nc; i++) {
      for (j=1; j<=nc; j++) {
	v1_ref(i,j).r = 0.0*(0.5*w.r*SCC_ref(i,j).r - 0.5*w.i*SCC_ref(i,j).i) - 0.5*HCCk_ref(i,j).r - SigmaL_ref(i,j).r; 
	v1_ref(i,j).i = 0.0*(0.5*w.i*SCC_ref(i,j).r + 0.5*w.r*SCC_ref(i,j).i) - 0.5*HCCk_ref(i,j).i - SigmaL_ref(i,j).i;
      }
    }
  }
  else{
    for (i=1; i<=nc; i++) {
      for (j=1; j<=nc; j++) {
	v1_ref(i,j).r = 0.0*(0.5*w.r*SCC_ref(i,j).r - 0.5*w.i*SCC_ref(i,j).i) - 0.5*HCCk_ref(i,j).r - SigmaR_ref(i,j).r; 
	v1_ref(i,j).i = 0.0*(0.5*w.i*SCC_ref(i,j).r + 0.5*w.r*SCC_ref(i,j).i) - 0.5*HCCk_ref(i,j).i - SigmaR_ref(i,j).i;
      }
    }
  }

  /* v2 = G(z) [1/2z S - 1/2 H - \sigama_{L or R}(z^*)]  */

  F77_NAME(zgemm,ZGEMM)("N","N", &nc, &nc, &nc, &alpha, GC, &nc, v1, &nc, &beta, v2, &nc);

  /* v1 = G(z) [1/2z S - 1/2 H - \sigama_{L or R}(z)] G(z^*)  */

  F77_NAME(zgemm,ZGEMM)("N","N", &nc, &nc, &nc, &alpha, v2, &nc, GC_Ad, &nc, &beta, v1, &nc);

  /******************************************************

    -1/(i 2Pi) times
    (retarded Green's function
    minus
    advanced Green's function of the left or right part)

  ******************************************************/

  for (i=1; i<=nc; i++) {
    for (j=1; j<=nc; j++) {

      ctmp.r = (Gless_ref(i,j).r - v1_ref(i,j).r)/(2.0*PI);
      ctmp.i = (Gless_ref(i,j).i - v1_ref(i,j).i)/(2.0*PI);

      Gless_ref(i,j).r =-ctmp.i;
      Gless_ref(i,j).i = ctmp.r;
    }
  }

}



