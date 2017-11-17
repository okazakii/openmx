/**********************************************************************
  IS_Hotelling.c:

     IS_Hotelling.c is a subroutine to calculate the inverse of overlap
     matrix using Hotelling's method.

  Log of IS_Hotelling.c:

     22/Nov/2001  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "openmx_common.h"

static void Hotelling(int Hotelling_switch, double ****OLP, double ****IOLP);

void IS_Hotelling(double ****OLP, double ****IOLP, int rlmax_IS)
{
  static int ct_AN,Hotelling_switch;
  static int wan,tno1,i,ig,ino1,k,p,q;
  static double dum1,dum2,dum;

  Hotelling_switch = rlmax_IS;
  Hotelling(Hotelling_switch,OLP,IOLP);

  /****************************************************
                  Averaging of IOLP
  ****************************************************/

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    wan = WhatSpecies[ct_AN];
    tno1 = Spe_Total_NO[wan] - 1;
    for (i=1; i<=(FNAN[ct_AN]+SNAN[ct_AN]); i++){
      ig = natn[ct_AN][i];
      ino1 = Spe_Total_NO[WhatSpecies[ig]] - 1;
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


void Hotelling(int Hotelling_switch, double ****OLP, double ****IOLP)
{
  static int firsttime=1;
  static int ct_AN,i,j,k,fan,san,can,order;
  static int ig,ian,jg,jan,m,n,kg,kan,kl,l,jl;
  static int wan,tno0,tno1,ino1,p,q;
  static int Cwan,Hwan,h_AN,Gh_AN,size_XP;
  static double sum,dum1,dum2,dum;
  static double ****XP;
  static double ****SXP;

  /****************************************************
    allocation of arrays:

    static double XP[atomnum+1]
                    [FNAN+SNAN+1]
                    [Spe_Total_NO[Cwan]]
                    [Spe_Total_NO[Hwan]] 

    static double SXP[atomnum+1]
                     [FNAN+SNAN+1]
                     [Spe_Total_NO[Cwan]]
                     [Spe_Total_NO[Hwan]] 
  ****************************************************/

  /* XP */

  size_XP = 0;
  FNAN[0] = 0;
  SNAN[0] = 0;

  XP = (double****)malloc(sizeof(double***)*(atomnum+1)); 
  for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

    if (ct_AN==0){
      tno0 = 1;
    }
    else{
      Cwan = WhatSpecies[ct_AN];
      tno0 = Spe_Total_NO[Cwan];  
    }    

    XP[ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN]+SNAN[ct_AN]+1)); 
    for (h_AN=0; h_AN<=(FNAN[ct_AN]+SNAN[ct_AN]); h_AN++){

      if (ct_AN==0){
	tno1 = 1;  
      }
      else{
	Gh_AN = natn[ct_AN][h_AN];
	Hwan = WhatSpecies[Gh_AN];
	tno1 = Spe_Total_NO[Hwan];
      } 

      XP[ct_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
      for (i=0; i<tno0; i++){
	XP[ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
      }
      size_XP += tno0*tno1;
    }
  }

  /* SXP */

  SXP = (double****)malloc(sizeof(double***)*(atomnum+1)); 
  FNAN[0] = 0;
  SNAN[0] = 0;
  for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

    if (ct_AN==0){
      tno0 = 1;
    }
    else{
      Cwan = WhatSpecies[ct_AN];
      tno0 = Spe_Total_NO[Cwan];  
    }    

    SXP[ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN]+SNAN[ct_AN]+1)); 
    for (h_AN=0; h_AN<=(FNAN[ct_AN]+SNAN[ct_AN]); h_AN++){

      if (ct_AN==0){
	tno1 = 1;  
      }
      else{
	Gh_AN = natn[ct_AN][h_AN];
	Hwan = WhatSpecies[Gh_AN];
	tno1 = Spe_Total_NO[Hwan];
      } 

      SXP[ct_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
      for (i=0; i<tno0; i++){
	SXP[ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
      }
    }
  }

  /* PrintMemory */
  if (firsttime){
    PrintMemory("IS_Hotelling: XP",      sizeof(double)*size_XP, NULL);
    PrintMemory("IS_Hotelling: SXP",     sizeof(double)*size_XP, NULL);
    firsttime=0;
  }

  /****************************************************
                   Initializing of XP
  ****************************************************/

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    fan = FNAN[ct_AN];
    san = SNAN[ct_AN];
    can = fan + san;
    ig = natn[ct_AN][0];
    ian = Spe_Total_NO[WhatSpecies[ig]];
    for (j=0; j<=can; j++){
      jg = natn[ct_AN][j];
      jan = Spe_Total_NO[WhatSpecies[jg]];

      if (j<=fan){
        for (m=0; m<=(ian-1); m++){
          for (n=0; n<=(jan-1); n++){
            if (j==0 && m==n){   
              XP[ct_AN][j][m][n] = 2.0 - OLP[ct_AN][j][m][n];
              IOLP[ct_AN][j][m][n] = 2.0 - OLP[ct_AN][j][m][n];
            }
            else {
              XP[ct_AN][j][m][n] = -OLP[ct_AN][j][m][n];
              IOLP[ct_AN][j][m][n] = -OLP[ct_AN][j][m][n];
            }
          }
        }
      }
      else {
        for (m=0; m<=(ian-1); m++){
          for (n=0; n<=(jan-1); n++){
            XP[ct_AN][j][m][n] = 0.0;
            IOLP[ct_AN][j][m][n] = 0.0;
          }
        }
      }

    }
  }

  /****************************************************
                       Iteration
  ****************************************************/

  for (order=1; order<=Hotelling_switch; order++){

    /****************************************************
                         2*XP -> IOLP
    ****************************************************/

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      fan = FNAN[ct_AN];
      san = SNAN[ct_AN];
      can = fan + san;
      ig = natn[ct_AN][0];
      ian = Spe_Total_NO[WhatSpecies[ig]];
      for (j=0; j<=can; j++){
        jg = natn[ct_AN][j];
        jan = Spe_Total_NO[WhatSpecies[jg]];
        for (m=0; m<=(ian-1); m++){
          for (n=0; n<=(jan-1); n++){
            IOLP[ct_AN][j][m][n] = 2.0*XP[ct_AN][j][m][n];
          }
        }
      }
    }

    /****************************************************
                         S*XP -> SXP
    ****************************************************/

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      fan = FNAN[ct_AN];
      san = SNAN[ct_AN];
      can = fan + san;
      ig = natn[ct_AN][0];
      ian = Spe_Total_NO[WhatSpecies[ig]];
      for (j=0; j<=can; j++){
        jg = natn[ct_AN][j];
        jan = Spe_Total_NO[WhatSpecies[jg]];
        for (k=0; k<=fan; k++){
          kg = natn[ct_AN][k];
          kan = Spe_Total_NO[WhatSpecies[kg]];
          kl = RMI2[ct_AN][j][k];
          for (m=0; m<=(ian-1); m++){
            for (n=0; n<=(jan-1); n++){
              sum = 0.0;
              for (l=0; l<=(kan-1); l++){
                sum = sum + OLP[ct_AN][k][m][l]*XP[jg][kl][n][l];
              }
              if (k==0){
                SXP[ct_AN][j][m][n] = sum;
              }
              else {
                SXP[ct_AN][j][m][n] = SXP[ct_AN][j][m][n] + sum;
              } 
            }
          }          
        }
      }
    }

    /****************************************************
                      IOLP = IOLP - XP*SXP
    ****************************************************/

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      fan = FNAN[ct_AN];
      san = SNAN[ct_AN];
      can = fan + san;
      ig = natn[ct_AN][0];
      ian = Spe_Total_NO[WhatSpecies[ig]];
      for (j=0; j<=can; j++){
        jg = natn[ct_AN][j];
        jan = Spe_Total_NO[WhatSpecies[jg]];
        for (k=0; k<=can; k++){
          kg = natn[ct_AN][k];
          jl = RMI2[ct_AN][k][j];
          kan = Spe_Total_NO[WhatSpecies[kg]];
          kl = RMI2[ct_AN][j][k];
          for (m=0; m<=(ian-1); m++){
            for (n=0; n<=(jan-1); n++){
              sum = 0.0;
              for (l=0; l<=(kan-1); l++){
                sum = sum + XP[ct_AN][k][m][l]*SXP[kg][jl][l][n];
              }
              IOLP[ct_AN][j][m][n] = IOLP[ct_AN][j][m][n] - sum;
            }
          }          
        }
      }
    }

    /****************************************************
                    Averaging of IOLP
    ****************************************************/
 
    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      wan = WhatSpecies[ct_AN];
      tno1 = Spe_Total_NO[wan] - 1;
      for (i=1; i<=(FNAN[ct_AN]+SNAN[ct_AN]); i++){
        ig = natn[ct_AN][i];
        ino1 = Spe_Total_NO[WhatSpecies[ig]] - 1;
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

    /****************************************************
                        IOLP -> XP
    ****************************************************/

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      fan = FNAN[ct_AN];
      san = SNAN[ct_AN];
      can = fan + san;
      ig = natn[ct_AN][0];
      ian = Spe_Total_NO[WhatSpecies[ig]];
      for (j=0; j<=can; j++){
        jg = natn[ct_AN][j];
        jan = Spe_Total_NO[WhatSpecies[jg]];
        for (m=0; m<=(ian-1); m++){
          for (n=0; n<=(jan-1); n++){
            XP[ct_AN][j][m][n] = IOLP[ct_AN][j][m][n];
          }
        }
      }
    }
  }

  /****************************************************
    freeing of arrays:

    static double XP[atomnum+1]
                    [FNAN+SNAN+1]
                    [Spe_Total_NO[Cwan]]
                    [Spe_Total_NO[Hwan]] 

    static double SXP[atomnum+1]
                     [FNAN+SNAN+1]
                     [Spe_Total_NO[Cwan]]
                     [Spe_Total_NO[Hwan]] 
  ****************************************************/

  /* XP */  

  FNAN[0] = 0;
  SNAN[0] = 0;
  for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

    if (ct_AN==0){
      tno0 = 1;
    }
    else{
      Cwan = WhatSpecies[ct_AN];
      tno0 = Spe_Total_NO[Cwan];  
    }    

    for (h_AN=0; h_AN<=(FNAN[ct_AN]+SNAN[ct_AN]); h_AN++){

      if (ct_AN==0){
	tno1 = 1;  
      }
      else{
	Gh_AN = natn[ct_AN][h_AN];
	Hwan = WhatSpecies[Gh_AN];
	tno1 = Spe_Total_NO[Hwan];
      } 

      for (i=0; i<tno0; i++){
	free(XP[ct_AN][h_AN][i]);
      }
      free(XP[ct_AN][h_AN]);
    }
    free(XP[ct_AN]);
  }
  free(XP);

  /* SXP */  

  FNAN[0] = 0;
  SNAN[0] = 0;
  for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

    if (ct_AN==0){
      tno0 = 1;
    }
    else{
      Cwan = WhatSpecies[ct_AN];
      tno0 = Spe_Total_NO[Cwan];  
    }    

    for (h_AN=0; h_AN<=(FNAN[ct_AN]+SNAN[ct_AN]); h_AN++){

      if (ct_AN==0){
	tno1 = 1;  
      }
      else{
	Gh_AN = natn[ct_AN][h_AN];
	Hwan = WhatSpecies[Gh_AN];
	tno1 = Spe_Total_NO[Hwan];
      } 

      for (i=0; i<tno0; i++){
	free(SXP[ct_AN][h_AN][i]);
      }
      free(SXP[ct_AN][h_AN]);
    }
    free(SXP[ct_AN]);
  }
  free(SXP);

}





