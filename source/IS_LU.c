/**********************************************************************
  IS_LU.c:

     IS_LU.c is a subroutine to calculate the inverse of overlap
     matrix using the divide (LU) method.

  Log of IS_LU.c:

     22/Nov/2001  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "openmx_common.h"
#include "lapack_prototypes.h"

static void LU(int ct_AN, double ****S0, double ****IOLP);

void IS_LU(double ****S0, double ****IOLP)
{
  static int ct_AN;
  static int wan,tno1,i,ig,ino1,k,p,q;
  static double dum1,dum2,dum;

  printf("LU\n"); 


  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    LU(ct_AN,S0,IOLP);
  }

  /****************************************************
                  Averaging of IOLP
  ****************************************************/

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    wan = WhatSpecies[ct_AN];
    tno1 = Spe_Total_CNO[wan] - 1;
    for (i=0; i<=(FNAN[ct_AN]+SNAN[ct_AN]); i++){
      ig = natn[ct_AN][i];
      ino1 = Spe_Total_CNO[WhatSpecies[ig]] - 1;
      k = RMI2[ct_AN][i][0];
      for (p=0; p<=tno1; p++){
	for (q=0; q<=ino1; q++){
	  dum1 = IOLP[ct_AN][i][p][q];
	  dum2 = IOLP[ig][k][q][p];
	  dum = 0.50*(dum1 + dum2);
	  IOLP[ct_AN][i][p][q] = dum;
	  IOLP[ig][k][q][p] = dum;
	}
      }
    } 
  }
}

void LU(int ct_AN, double ****S0, double ****IOLP)
{
  INTEGER N,lda,info,lwork;
  INTEGER *piv;
  int **jun;
  int i0,j0; 
  int i,j,k,l;
  int ig,ian,jg,jan,kg,kan,qb;
  int po,fan,san,can,wan,m,n,kl;
  double *LoS,*work;

  /****************************************************
    allocation of arrays:

     static INTEGER piv[List_YOUSO[2]*List_YOUSO[7]];
     static int jun[List_YOUSO[2]][List_YOUSO[7]];
  ****************************************************/

  piv = (INTEGER*)malloc(sizeof(INTEGER)*List_YOUSO[2]*List_YOUSO[7]);
  jun = (int**)malloc(sizeof(int*)*List_YOUSO[2]); 
  for (i=0; i<List_YOUSO[2]; i++){
    jun[i] = (int*)malloc(sizeof(int)*List_YOUSO[7]); 
  }

  /****************************************************
                  Setting of S-matrix           
  ****************************************************/

  fan = FNAN[ct_AN];
  san = SNAN[ct_AN];
  can = fan + san;

  N = 0;
  for (i=0; i<=can; i++){
    ig = natn[ct_AN][i];
    ian = Spe_Total_CNO[WhatSpecies[ig]];
    for (k=0; k<ian; k++){
      N++;
      jun[i][k] = N-1; 
    }
  }  

  LoS = (double*)malloc(sizeof(double)*N*N);
  work = (double*)malloc(sizeof(double)*N);

  for (i=0; i<=can; i++){
    ig = natn[ct_AN][i];
    ian = Spe_Total_CNO[WhatSpecies[ig]];
    for (j=0; j<=can; j++){
      kl = RMI1[ct_AN][i][j];
      jg = natn[ct_AN][j];
      jan = Spe_Total_CNO[WhatSpecies[jg]];

      if (0<=kl){
	for (m=0; m<=(ian-1); m++){
          i0 = jun[i][m];
	  for (n=0; n<=(jan-1); n++){
            j0 = jun[j][n];
            k = N*j0 + i0;
            LoS[k] = S0[ig][kl][m][n];
	  }
	}
      }
      else{
	for (m=0; m<=(ian-1); m++){
          i0 = jun[i][m];
	  for (n=0; n<=(jan-1); n++){
            j0 = jun[j][n];
            k = N*j0 + i0; 
            LoS[k] = 0.0;
	  }
	}
      }
    }
  }

  /****************************************************
                Call dgetrf_() in clapack
  ****************************************************/

  lda = N;
  F77_NAME(dgetrf,DGETRF)(&N, &N, LoS, &lda, piv, &info);

  if (info!=0){
    printf("Error in dgetrf_() which is called from IS_LU  info=%2d\n",info);
  }

  /****************************************************
                Call dgetri_() in clapack
  ****************************************************/

  lwork = N;
  F77_NAME(dgetri,DGETRI)(&N, LoS, &lda, piv, work, &lwork, &info);
  if (info!=0){
    printf("Error in dgetri_() which is called from IS_LU  info=%2d\n",info);
  }

  /****************************************************
                       LoS -> IOLP
  ****************************************************/

  i = 0;
  ig = natn[ct_AN][i];
  ian = Spe_Total_CNO[WhatSpecies[ig]];
  for (j=0; j<=can; j++){
    jg = natn[ct_AN][j];
    jan = Spe_Total_CNO[WhatSpecies[jg]];
    for (m=0; m<=(ian-1); m++){
      i0 = jun[i][m];
      for (n=0; n<=(jan-1); n++){
        j0 = jun[j][n];
        k = N*j0 + i0;
        IOLP[ct_AN][j][m][n] = LoS[k];
      }
    }
  }

  free(LoS);
  free(work);

  /****************************************************
    freeing of arrays:

     static INTEGER piv[List_YOUSO[2]*List_YOUSO[7]];
     static int jun[List_YOUSO[2]][List_YOUSO[7]];
  ****************************************************/

  free(piv);
  for (i=0; i<List_YOUSO[2]; i++){
    free(jun[i]);
  }
  free(jun);

}




