/**********************************************************************
  IS_Taylor.c:

     IS_Taylor.c is a subroutine to calculate the inverse of overlap
     matrix using the Taylor expansion method.

  Log of IS_Taylor.c:

     22/Nov/2001  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "openmx_common.h"

static void Taylor(int Taylor_switch, double ****OLP, double ****IOLP);

void IS_Taylor(double ****OLP, double ****IOLP, int rlmax_IS)
{
  static int ct_AN,Taylor_switch;
  static int wan,tno1,i,ig,ino1,k,p,q;
  static double dum1,dum2,dum;

  Taylor_switch = rlmax_IS;
  Taylor(Taylor_switch,OLP,IOLP);

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


void Taylor(int Taylor_switch, double ****OLP, double ****IOLP)
{
  static int firsttime=1;
  static int ct_AN,fan,san,can,ig,ian,j,jan,m,n,order;
  static int k,kg,kan,kl,l,jg,wan,tno0,tno1,i,ino1,p,q,jl;
  static int h_AN,Gh_AN,Hwan,Cwan,size_O,size_OP;

  static double dum,sum,dum1,dum2;
  static double ****O;
  static double ****OP;
  static double ****ON;

  /****************************************************
    allocation of arrays:

    static double O[atomnum+1]
                   [FNAN+1]
                   [Spe_Total_NO[Cwan]]
                   [Spe_Total_NO[Hwan]] 

    static double OP[atomnum+1]
                    [FNAN+SNAN+1]
                    [Spe_Total_NO[Cwan]]
                    [Spe_Total_NO[Hwan]] 

    static double ON[atomnum+1]
                    [FNAN+SNAN+1]
                    [Spe_Total_NO[Cwan]]
                    [Spe_Total_NO[Hwan]] 
  ****************************************************/

  /* O */

  FNAN[0] = 0;
  SNAN[0] = 0;
  size_O = 0;

  O = (double****)malloc(sizeof(double***)*(atomnum+1)); 
  for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

    if (ct_AN==0){
      tno0 = 1;
    }
    else{
      Cwan = WhatSpecies[ct_AN];
      tno0 = Spe_Total_NO[Cwan];  
    }    

    O[ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN]+1)); 
    for (h_AN=0; h_AN<=(FNAN[ct_AN]+SNAN[ct_AN]); h_AN++){

      if (ct_AN==0){
	tno1 = 1;  
      }
      else{
	Gh_AN = natn[ct_AN][h_AN];
	Hwan = WhatSpecies[Gh_AN];
	tno1 = Spe_Total_NO[Hwan];
      } 

      O[ct_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
      for (i=0; i<tno0; i++){
	O[ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
      }
      size_O += tno0*tno1;
    }
  }

  /* OP */

  size_OP = 0;
  OP = (double****)malloc(sizeof(double***)*(atomnum+1)); 
  for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

    if (ct_AN==0){
      tno0 = 1;
    }
    else{
      Cwan = WhatSpecies[ct_AN];
      tno0 = Spe_Total_NO[Cwan];  
    }    

    OP[ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN]+SNAN[ct_AN]+1)); 
    for (h_AN=0; h_AN<=(FNAN[ct_AN]+SNAN[ct_AN]); h_AN++){

      if (ct_AN==0){
	tno1 = 1;  
      }
      else{
	Gh_AN = natn[ct_AN][h_AN];
	Hwan = WhatSpecies[Gh_AN];
	tno1 = Spe_Total_NO[Hwan];
      } 

      OP[ct_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
      for (i=0; i<tno0; i++){
	OP[ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
      }
      size_OP += tno0*tno1;
    }
  }

  /* ON */

  ON = (double****)malloc(sizeof(double***)*(atomnum+1)); 
  for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

    if (ct_AN==0){
      tno0 = 1;
    }
    else{
      Cwan = WhatSpecies[ct_AN];
      tno0 = Spe_Total_NO[Cwan];  
    }    

    ON[ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN]+SNAN[ct_AN]+1)); 
    for (h_AN=0; h_AN<=(FNAN[ct_AN]+SNAN[ct_AN]); h_AN++){

      if (ct_AN==0){
	tno1 = 1;  
      }
      else{
	Gh_AN = natn[ct_AN][h_AN];
	Hwan = WhatSpecies[Gh_AN];
	tno1 = Spe_Total_NO[Hwan];
      } 

      ON[ct_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
      for (i=0; i<tno0; i++){
	ON[ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
      }
    }
  }


  /* PrintMemory */
  if (firsttime) {
    PrintMemory("Taylor: O",sizeof(O),NULL);
    PrintMemory("Taylor: OP",sizeof(OP),NULL);
    PrintMemory("Taylor: ON",sizeof(ON),NULL);
    firsttime=0;
  }

  /****************************************************
        Initializing of IOLP, OP and O matrices
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
              IOLP[ct_AN][j][m][n] = 1.0;
              O[ct_AN][j][m][n] = 0.0;
              OP[ct_AN][j][m][n] = 0.0;
            }
            else {
              IOLP[ct_AN][j][m][n] = -OLP[ct_AN][j][m][n];
              O[ct_AN][j][m][n] = OLP[ct_AN][j][m][n];
              OP[ct_AN][j][m][n] = OLP[ct_AN][j][m][n];
            }
          }
        }
      }
      else {
        for (m=0; m<=(ian-1); m++){
          for (n=0; n<=(jan-1); n++){
            IOLP[ct_AN][j][m][n] = 0.0;
            O[ct_AN][j][m][n] = 0.0;
            OP[ct_AN][j][m][n] = 0.0;
          }
        }
      }
    }
  }

  /****************************************************
                      Iteration
  ****************************************************/

  for (order=2; order<=Taylor_switch; order++){

    /****************************************************
                         O*OP -> ON
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
          jl = RMI2[ct_AN][k][j];
          kan = Spe_Total_NO[WhatSpecies[kg]];
          kl = RMI2[ct_AN][j][k];
          for (m=0; m<=(ian-1); m++){
            for (n=0; n<=(jan-1); n++){
              sum = 0.0;
              for (l=0; l<=(kan-1); l++){
                sum = sum + O[ct_AN][k][m][l]*OP[kg][jl][l][n];
              }
              if (k==0){
                ON[ct_AN][j][m][n] = sum;
              }
              else {
                ON[ct_AN][j][m][n] = ON[ct_AN][j][m][n] + sum;
              } 
            }
          }          
        }
      }
    }

    /****************************************************
                       Averaging of ON
    ****************************************************/

    if (order%2==0){
      for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
        wan = WhatSpecies[ct_AN];
        tno1 = Spe_Total_NO[wan] - 1;
        for (i=1; i<=(FNAN[ct_AN]+SNAN[ct_AN]); i++){
          ig = natn[ct_AN][i];
          ino1 = Spe_Total_NO[WhatSpecies[ig]] - 1;
          k = RMI2[ct_AN][i][0];
          for (p=0; p<=tno1; p++){
	    for (q=0; q<=ino1; q++){
	      dum1 = ON[ct_AN][i][p][q];
	      dum2 = ON[ig][k][q][p];
	      dum = 0.50*(dum1 + dum2);
	      ON[ct_AN][i][p][q] = dum;
	      ON[ig][k][q][p] = dum;
	    }
          }
        } 
      }
    }

    /****************************************************
                         ON -> OP
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
            OP[ct_AN][j][m][n] = ON[ct_AN][j][m][n];
          }
        }
      }
    }

    /****************************************************
                  IOLP = IOLP +- ON
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
            if ((order%2)==0) {
              IOLP[ct_AN][j][m][n] = IOLP[ct_AN][j][m][n]
                                    + ON[ct_AN][j][m][n];
            }
            else {
              IOLP[ct_AN][j][m][n] = IOLP[ct_AN][j][m][n]
                                    - ON[ct_AN][j][m][n];
            } 
          }
        }
      }
    }
  }


  /****************************************************
    freeing of arrays:

    static double O[atomnum+1]
                   [FNAN+1]
                   [Spe_Total_NO[Cwan]]
                   [Spe_Total_NO[Hwan]] 

    static double OP[atomnum+1]
                    [FNAN+SNAN+1]
                    [Spe_Total_NO[Cwan]]
                    [Spe_Total_NO[Hwan]] 

    static double ON[atomnum+1]
                    [FNAN+SNAN+1]
                    [Spe_Total_NO[Cwan]]
                    [Spe_Total_NO[Hwan]] 
  ****************************************************/

  /* O */

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
	free(O[ct_AN][h_AN][i]);
      }
      free(O[ct_AN][h_AN]);
    }
    free(O[ct_AN]);
  }
  free(O);

  /* OP */

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
	free(OP[ct_AN][h_AN][i]);
      }
      free(OP[ct_AN][h_AN]);
    }
    free(OP[ct_AN]);
  }
  free(OP);

  /* ON */

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
	free(ON[ct_AN][h_AN][i]);
      }
      free(ON[ct_AN][h_AN]);
    }
    free(ON[ct_AN]);
  }
  free(ON);

}




