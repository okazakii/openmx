/**********************************************************************
  RecursionS.c:

     RecursionS.c is a subroutine to solve a generalized eigenvalue
     problem with an overlap matrix for both spin and non-spin
     polarized systems using Ozaki's recursion method in O(N) operation.

  Log of RecursionS.c:

     22/Nov/2001  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "openmx_common.h"
#include "lapack_prototypes.h"

static void B_Lanczos(int ct_AN, int spin, 
                      double ****hopH, double ****OLP0, double ****HP);
static void Solve_EigenValue(int ct_AN, int spin);
static void Solve_EigenValue2(int ct_AN, int spin);
static void Search_ChemP(double ****OLP0);
static void Calc_DM(double *****CDM, double *****EDM, double ****IOLP, double *****Ham);
static int SVD_fact_inv(int rl, int n, double **a,
                        double **B, double **IB,
                        double **C, double **IC);
static void ISH(double ****IOLP, double ****nh, double ****HP);
static void SbyISH(double ****OLP0, double ****HP,
		   double ****nh, double ****DH);
static void Save_Recursion();
static void Output_RcnCof(FILE *fp);

static double *****al;
static double *****be;
static double *****ibe;
static double *****ga;
static double *****iga;
static double *****REVec;
static double ***EigenValue;
static double ******LanU;
static double Uele_OS[2],Uele_IS[2];
static int **RCN_Rank;

double RecursionS_C(int SCF_iter,
                    double *****Ham, double ****OLP0,
                    double *****CDM,
                    double *****EDM,
                    double Eele0[2], double Eele1[2])
{
  static int firsttime=1;
  static int i,ct_AN,spin,h_AN,Gh_AN,Hwan;
  static int m,n,j,Rn,k,l,rl,tno0,tno1,Cwan;
  static int size_al,size_REVec,size_EigenValue;
  static int size_LanU,size_HPS;
  static double time0;
  static double sum;
  static double *****HPS;
  static double *****DH;
  static double TStime,TEtime;

  dtime(&TStime);
  Timetool_times("RecursionS_C","start");

  /****************************************************
    allocation of arrays:

   static double al[YOUSO23][YOUSO1][YOUSO3][YOUSO7][YOUSO7];
   static double be[YOUSO23][YOUSO1][YOUSO3][YOUSO7][YOUSO7];
   static double ibe[YOUSO23][YOUSO1][YOUSO3][YOUSO7][YOUSO7];
   static double ga[YOUSO23][YOUSO1][YOUSO3][YOUSO7][YOUSO7];
   static double iga[YOUSO23][YOUSO1][YOUSO3][YOUSO7][YOUSO7];
   static double REVec[YOUSO23][YOUSO1][(YOUSO3+1)*YOUSO7][YOUSO3+1][YOUSO7];
   static double EigenValue[YOUSO23][YOUSO1][(YOUSO3+1)*YOUSO7];
   static int RCN_Rank[YOUSO1][YOUSO3];
   static double LanU[YOUSO23][YOUSO1][YOUSO3+1][YOUSO8][YOUSO7][YOUSO7];

   static double HPS[YOUSO23][YOUSO1][YOUSO2][YOUSO7][YOUSO7];
   static double DH[YOUSO23][YOUSO1][YOUSO2][YOUSO7][YOUSO7];
  ****************************************************/
   
  /* al */  

  size_al = 0;
  al = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
  for (k=0; k<=SpinP_switch; k++){
    al[k] = (double****)malloc(sizeof(double***)*(atomnum+1)); 

    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

      if (ct_AN==0){
        tno0 = 1;
      }
      else{
        Cwan = WhatSpecies[ct_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    

      al[k][ct_AN] = (double***)malloc(sizeof(double**)*List_YOUSO[3]); 
      for (rl=0; rl<List_YOUSO[3]; rl++){

        al[k][ct_AN][rl] = (double**)malloc(sizeof(double*)*tno0); 
        for (i=0; i<tno0; i++){
          al[k][ct_AN][rl][i] = (double*)malloc(sizeof(double)*tno0); 
        }
        size_al += tno0*tno0;  

      }
    }
  }

  /* be */

  be = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
  for (k=0; k<=SpinP_switch; k++){
    be[k] = (double****)malloc(sizeof(double***)*(atomnum+1)); 

    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

      if (ct_AN==0){
        tno0 = 1;
      }
      else{
        Cwan = WhatSpecies[ct_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    

      be[k][ct_AN] = (double***)malloc(sizeof(double**)*List_YOUSO[3]); 
      for (rl=0; rl<List_YOUSO[3]; rl++){

        be[k][ct_AN][rl] = (double**)malloc(sizeof(double*)*tno0); 
        for (i=0; i<tno0; i++){
          be[k][ct_AN][rl][i] = (double*)malloc(sizeof(double)*tno0); 
        }
      }
    }
  }

  /* ibe */

  ibe = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
  for (k=0; k<=SpinP_switch; k++){
    ibe[k] = (double****)malloc(sizeof(double***)*(atomnum+1)); 

    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

      if (ct_AN==0){
        tno0 = 1;
      }
      else{
        Cwan = WhatSpecies[ct_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    

      ibe[k][ct_AN] = (double***)malloc(sizeof(double**)*List_YOUSO[3]); 
      for (rl=0; rl<List_YOUSO[3]; rl++){

        ibe[k][ct_AN][rl] = (double**)malloc(sizeof(double*)*tno0); 
        for (i=0; i<tno0; i++){
          ibe[k][ct_AN][rl][i] = (double*)malloc(sizeof(double)*tno0); 
        }
      }
    }
  }

  /* ga */

  ga = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
  for (k=0; k<=SpinP_switch; k++){
    ga[k] = (double****)malloc(sizeof(double***)*(atomnum+1)); 

    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

      if (ct_AN==0){
        tno0 = 1;
      }
      else{
        Cwan = WhatSpecies[ct_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    

      ga[k][ct_AN] = (double***)malloc(sizeof(double**)*List_YOUSO[3]); 
      for (rl=0; rl<List_YOUSO[3]; rl++){

        ga[k][ct_AN][rl] = (double**)malloc(sizeof(double*)*tno0); 
        for (i=0; i<tno0; i++){
          ga[k][ct_AN][rl][i] = (double*)malloc(sizeof(double)*tno0); 
        }
      }
    }
  }

  /* iga */

  iga = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
  for (k=0; k<=SpinP_switch; k++){
    iga[k] = (double****)malloc(sizeof(double***)*(atomnum+1)); 

    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

      if (ct_AN==0){
        tno0 = 1;
      }
      else{
        Cwan = WhatSpecies[ct_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    

      iga[k][ct_AN] = (double***)malloc(sizeof(double**)*List_YOUSO[3]); 
      for (rl=0; rl<List_YOUSO[3]; rl++){

        iga[k][ct_AN][rl] = (double**)malloc(sizeof(double*)*tno0); 
        for (i=0; i<tno0; i++){
          iga[k][ct_AN][rl][i] = (double*)malloc(sizeof(double)*tno0); 
        }
      }
    }
  }

  /* REVec */

  size_REVec = 0;
  REVec = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
  for (spin=0; spin<=SpinP_switch; spin++){
    REVec[spin] = (double****)malloc(sizeof(double***)*(atomnum+1)); 
    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

      if (ct_AN==0){
        tno0 = 1;
      }
      else{
        Cwan = WhatSpecies[ct_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    

      REVec[spin][ct_AN] = (double***)malloc(sizeof(double**)*(List_YOUSO[3]+1)*tno0); 

      for (i=0; i<(List_YOUSO[3]+1)*tno0; i++){
        REVec[spin][ct_AN][i] = (double**)malloc(sizeof(double*)*(List_YOUSO[3]+1)); 
        for (j=0; j<(List_YOUSO[3]+1); j++){
          REVec[spin][ct_AN][i][j] = (double*)malloc(sizeof(double)*tno0); 
	}
        size_REVec += (List_YOUSO[3]+1)*tno0;
      }             
    }
  }


  /* EigenValue */

  size_EigenValue = 0;
  EigenValue = (double***)malloc(sizeof(double**)*(SpinP_switch+1)); 
  for (spin=0; spin<=SpinP_switch; spin++){
    EigenValue[spin] = (double**)malloc(sizeof(double*)*(atomnum+1)); 
    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

      if (ct_AN==0){
        tno0 = 1;
      }
      else{
        Cwan = WhatSpecies[ct_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    

      EigenValue[spin][ct_AN] = (double*)malloc(sizeof(double)*(List_YOUSO[3]+1)*tno0); 
      size_EigenValue += tno0;
    }
  }

  /* RCN_Rank */

  RCN_Rank = (int**)malloc(sizeof(int*)*List_YOUSO[1]); 
  for (i=0; i<List_YOUSO[1]; i++){
    RCN_Rank[i] = (int*)malloc(sizeof(int)*List_YOUSO[3]); 
  }

  /* LanU */

  size_LanU = 0;
  FNAN[0] = 0;
  LanU = (double******)malloc(sizeof(double*****)*(SpinP_switch+1)); 
  for (spin=0; spin<=SpinP_switch; spin++){
    LanU[spin] = (double*****)malloc(sizeof(double****)*(atomnum+1)); 

    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

      if (ct_AN==0){
        tno0 = 1;
      }
      else{
        Cwan = WhatSpecies[ct_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    

      LanU[spin][ct_AN] = (double****)malloc(sizeof(double***)*(List_YOUSO[3]+1)); 

      for (rl=0; rl<(List_YOUSO[3]+1); rl++){

        LanU[spin][ct_AN][rl] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN]+1)); 

        for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){

          if (ct_AN==0){
            tno1 = 1;  
          }
          else{
            Gh_AN = natn[ct_AN][h_AN];
            Hwan = WhatSpecies[Gh_AN];
            tno1 = Spe_Total_NO[Hwan];
          } 

          LanU[spin][ct_AN][rl][h_AN] = (double**)malloc(sizeof(double*)*tno1); 
          for (i=0; i<tno1; i++){
            LanU[spin][ct_AN][rl][h_AN][i] = (double*)malloc(sizeof(double)*tno0); 
            for (j=0; j<tno0; j++){
              LanU[spin][ct_AN][rl][h_AN][i][j] = 0.0;
            }
          }
          size_LanU += tno1*tno0;
        }
      }
    }
  }

  /* HPS */  

  size_HPS = 0;
  FNAN[0] = 0;
  SNAN[0] = 0;

  HPS = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
  for (k=0; k<=SpinP_switch; k++){
    HPS[k] = (double****)malloc(sizeof(double***)*(atomnum+1)); 

    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

      if (ct_AN==0){
        tno0 = 1;
      }
      else{
        Cwan = WhatSpecies[ct_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    

      HPS[k][ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN]+SNAN[ct_AN]+1));
      for (h_AN=0; h_AN<=(FNAN[ct_AN]+SNAN[ct_AN]); h_AN++){

        if (ct_AN==0){
          tno1 = 1;  
        }
        else{
          Gh_AN = natn[ct_AN][h_AN];
          Hwan = WhatSpecies[Gh_AN];
          tno1 = Spe_Total_NO[Hwan];
        } 

        HPS[k][ct_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
        for (i=0; i<tno0; i++){
          HPS[k][ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
        }
        size_HPS += tno0*tno1;
      }
    }
  }


  /* DH */  

  FNAN[0] = 0;
  SNAN[0] = 0;

  DH = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
  for (k=0; k<=SpinP_switch; k++){
    DH[k] = (double****)malloc(sizeof(double***)*(atomnum+1)); 
    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

      if (ct_AN==0){
        tno0 = 1;
      }
      else{
        Cwan = WhatSpecies[ct_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    

      DH[k][ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN]+SNAN[ct_AN]+1));
      for (h_AN=0; h_AN<=(FNAN[ct_AN]+SNAN[ct_AN]); h_AN++){

        if (ct_AN==0){
          tno1 = 1;  
        }
        else{
          Gh_AN = natn[ct_AN][h_AN];
          Hwan = WhatSpecies[Gh_AN];
          tno1 = Spe_Total_NO[Hwan];
        } 

        DH[k][ct_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
        for (i=0; i<tno0; i++){
          DH[k][ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
        }
      }
    }
  }

  /****************************************************
                      PrintMermoy
  ****************************************************/

  if (firsttime) {
    PrintMemory("RecursionS_C: al",sizeof(double)*size_al,NULL);
    PrintMemory("RecursionS_C: be",sizeof(double)*size_al,NULL);
    PrintMemory("RecursionS_C: ibe",sizeof(double)*size_al,NULL);
    PrintMemory("RecursionS_C: ga",sizeof(double)*size_al,NULL);
    PrintMemory("RecursionS_C: iga",sizeof(double)*size_al,NULL);
    PrintMemory("RecursionS_C: REVec",sizeof(double)*size_REVec,NULL);
    PrintMemory("RecursionS_C: EigenValue",sizeof(double)*size_EigenValue,NULL);
    PrintMemory("RecursionS_C: LanU",sizeof(double)*size_LanU,NULL);
    PrintMemory("RecursionS_C: HPS",sizeof(double)*size_HPS,NULL);
    PrintMemory("RecursionS_C: DH",sizeof(double)*size_HPS,NULL);
    firsttime=0;
  }

  /****************************************************
     Calculate the inverse matrix of overlap matrix
  ****************************************************/

  printf("A IS_switch=%2d\n",IS_switch);

  if (SCF_iter==1){
    if (IS_switch==1) IS_Lanczos(OLP0,IOLP,rlmax_IS);
    if (IS_switch==2) IS_LU(OLP0,IOLP);
    if (IS_switch==3) IS_Hotelling(OLP0,IOLP,rlmax_IS);
    if (IS_switch==4) IS_Taylor(OLP0,IOLP,rlmax_IS);
  }  

  /****************************************************
                   Calculate S^{-1}H
  ****************************************************/

  for (spin=0; spin<=SpinP_switch; spin++){
    ISH(IOLP,Ham[spin],HPS[spin]);
    /*
    SbyISH(OLP0,HPS[spin],Ham[spin],DH[spin]);
    */
  }

  /****************************************************
  Applying the block Lanczos transform. to each cluster
  ****************************************************/

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    for (spin=0; spin<=SpinP_switch; spin++){
      B_Lanczos(ct_AN,spin,Ham[spin],OLP0,HPS[spin]);
    }
  }

  /****************************************************
     Solve the eigenvalue problem of the matrix block
              tri-diagonalized by the TSB
  ****************************************************/

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    for (spin=0; spin<=SpinP_switch; spin++){
      Solve_EigenValue2(ct_AN,spin);
    }
  }

  /****************************************************
          Search the chemical potential to maintain
                the number of electrons
  ****************************************************/

  Search_ChemP(OLP0);

  /****************************************************
             Calculate the density matrix
  ****************************************************/

  Calc_DM(CDM,EDM,IOLP,Ham);

  printf("Uele_IS[0]=%15.12f\n",Uele_IS[0]);
  printf("Uele_IS[1]=%15.12f\n",Uele_IS[1]);

  printf("Uele_OS[0]=%15.12f\n",Uele_OS[0]);
  printf("Uele_OS[1]=%15.12f\n",Uele_OS[1]);




  if (SpinP_switch==0){
    Eele0[0] = Uele_OS[0];
    Eele0[1] = Uele_OS[0];
    Eele1[0] = Uele_IS[0];   
    Eele1[1] = Uele_IS[0];
  }
  else if (SpinP_switch==1){
    Eele0[0] = Uele_OS[0];
    Eele0[1] = Uele_OS[1];
    Eele1[0] = Uele_IS[0];   
    Eele1[1] = Uele_IS[1];
  }
  
  /****************************************************
   Output several informations of recursion algorithm
  ****************************************************/

  if (2<=level_fileout) Save_Recursion();


  /****************************************************
   freeing of arrays:

   static double al[YOUSO23][YOUSO1][YOUSO3][YOUSO7][YOUSO7];
   static double be[YOUSO23][YOUSO1][YOUSO3][YOUSO7][YOUSO7];
   static double ibe[YOUSO23][YOUSO1][YOUSO3][YOUSO7][YOUSO7];
   static double ga[YOUSO23][YOUSO1][YOUSO3][YOUSO7][YOUSO7];
   static double iga[YOUSO23][YOUSO1][YOUSO3][YOUSO7][YOUSO7];
   static double REVec[YOUSO23][YOUSO1][(YOUSO3+1)*YOUSO7][YOUSO3+1][YOUSO7];
   static double EigenValue[YOUSO23][YOUSO1][(YOUSO3+1)*YOUSO7];
   static int RCN_Rank[YOUSO1][YOUSO3];
   static double LanU[YOUSO23][YOUSO1][YOUSO3+1][YOUSO8][YOUSO7][YOUSO7];
  ****************************************************/
   
  /* al */  

  for (k=0; k<=SpinP_switch; k++){
    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

      if (ct_AN==0){
        tno0 = 1;
      }
      else{
        Cwan = WhatSpecies[ct_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    

      for (rl=0; rl<List_YOUSO[3]; rl++){
        for (i=0; i<tno0; i++){
          free(al[k][ct_AN][rl][i]);
        }
        free(al[k][ct_AN][rl]);
      }
      free(al[k][ct_AN]);
    }
    free(al[k]);
  }
  free(al);


  /* be */

  for (k=0; k<=SpinP_switch; k++){
    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

      if (ct_AN==0){
        tno0 = 1;
      }
      else{
        Cwan = WhatSpecies[ct_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    

      for (rl=0; rl<List_YOUSO[3]; rl++){
        for (i=0; i<tno0; i++){
          free(be[k][ct_AN][rl][i]);
        }
        free(be[k][ct_AN][rl]);
      }
      free(be[k][ct_AN]);
    }
    free(be[k]);
  }
  free(be);

  /* ibe */

  for (k=0; k<=SpinP_switch; k++){
    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

      if (ct_AN==0){
        tno0 = 1;
      }
      else{
        Cwan = WhatSpecies[ct_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    

      for (rl=0; rl<List_YOUSO[3]; rl++){
        for (i=0; i<tno0; i++){
          free(ibe[k][ct_AN][rl][i]);
        }
        free(ibe[k][ct_AN][rl]);
      }
      free(ibe[k][ct_AN]);
    }
    free(ibe[k]);
  }
  free(ibe);

  /* ga */

  for (k=0; k<=SpinP_switch; k++){
    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

      if (ct_AN==0){
        tno0 = 1;
      }
      else{
        Cwan = WhatSpecies[ct_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    

      for (rl=0; rl<List_YOUSO[3]; rl++){
        for (i=0; i<tno0; i++){
          free(ga[k][ct_AN][rl][i]);
        }
        free(ga[k][ct_AN][rl]);
      }
      free(ga[k][ct_AN]);
    }
    free(ga[k]);
  }
  free(ga);

  /* iga */

  for (k=0; k<=SpinP_switch; k++){
    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

      if (ct_AN==0){
        tno0 = 1;
      }
      else{
        Cwan = WhatSpecies[ct_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    

      for (rl=0; rl<List_YOUSO[3]; rl++){
        for (i=0; i<tno0; i++){
          free(iga[k][ct_AN][rl][i]);
        }
        free(iga[k][ct_AN][rl]);
      }
      free(iga[k][ct_AN]);
    }
    free(iga[k]);
  }
  free(iga);

  /* REVec */

  for (spin=0; spin<=SpinP_switch; spin++){
    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

      if (ct_AN==0){
        tno0 = 1;
      }
      else{
        Cwan = WhatSpecies[ct_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    


      for (i=0; i<(List_YOUSO[3]+1)*tno0; i++){
        for (j=0; j<(List_YOUSO[3]+1); j++){
          free(REVec[spin][ct_AN][i][j]);
	}
        free(REVec[spin][ct_AN][i]);
      }             
      free(REVec[spin][ct_AN]);
    }
    free(REVec[spin]);
  }
  free(REVec);

  /* EigenValue */

  for (spin=0; spin<=SpinP_switch; spin++){
    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

      if (ct_AN==0){
        tno0 = 1;
      }
      else{
        Cwan = WhatSpecies[ct_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    

      free(EigenValue[spin][ct_AN]);
    }
    free(EigenValue[spin]);
  }
  free(EigenValue);

  /* RCN_Rank */

  for (i=0; i<List_YOUSO[1]; i++){
    free(RCN_Rank[i]);
  }
  free(RCN_Rank);

  /* LanU */

  FNAN[0] = 0;
  for (spin=0; spin<=SpinP_switch; spin++){
    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

      if (ct_AN==0){
        tno0 = 1;
      }
      else{
        Cwan = WhatSpecies[ct_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    

      for (rl=0; rl<(List_YOUSO[3]+1); rl++){
        for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){

          if (ct_AN==0){
            tno1 = 1;  
          }
          else{
            Gh_AN = natn[ct_AN][h_AN];
            Hwan = WhatSpecies[Gh_AN];
            tno1 = Spe_Total_NO[Hwan];
          } 

          for (i=0; i<tno1; i++){
            free(LanU[spin][ct_AN][rl][h_AN][i]);
          }          
          free(LanU[spin][ct_AN][rl][h_AN]);
        }
        free(LanU[spin][ct_AN][rl]);
      }
      free(LanU[spin][ct_AN]);
    }
    free(LanU[spin]);
  }
  free(LanU);

  /* HPS */  

  FNAN[0] = 0;
  SNAN[0] = 0;

  for (k=0; k<=SpinP_switch; k++){
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
          free(HPS[k][ct_AN][h_AN][i]);
        }
        free(HPS[k][ct_AN][h_AN]);
      }
      free(HPS[k][ct_AN]);
    }
    free(HPS[k]);
  }
  free(HPS);

  /* DH */  

  FNAN[0] = 0;
  SNAN[0] = 0;

  for (k=0; k<=SpinP_switch; k++){
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
          free(DH[k][ct_AN][h_AN][i]);
        }
        free(DH[k][ct_AN][h_AN]);
      }
      free(DH[k][ct_AN]);
    }
    free(DH[k]);
  }
  free(DH);

  /* for time */

  dtime(&TEtime);
  Timetool_times("RecursionS_C","end");
  time0 = TEtime - TStime;
  return time0;
}

void Calc_DM(double *****CDM, double *****EDM, double ****IOLP, double *****Ham)
{
  static int firsttime=1;
  static int ct_AN,spin,N,wanA,wanB,tnoA,tnoB;
  static int wan,tno,tno1,itnum,ig,ino1,jg,jno1,kg,kj,kno1;
  static int i,j,k,p,l,m,n,q,rl,i1,j1,h_AN,Gh_AN,hwan,htno,rl1,rl2;
  static double tmp,sum,sum2,sumx,sumy,sumz,dtnum;
  static double **Uele_temp;
  static double ***UG;
  static double ***UGUt;
  static double x,FermiF,dum,dum1,dum2;

  /****************************************************
    allocation of arrays:

    static double Uele_temp[YOUSO23][YOUSO1];
    static double UG[YOUSO3][YOUSO7][YOUSO7];
    static double UGUt[YOUSO8][YOUSO7][YOUSO7];
  ****************************************************/

  Uele_temp = (double**)malloc(sizeof(double*)*(SpinP_switch+1)); 
  for (spin=0; spin<=SpinP_switch; spin++){
    Uele_temp[spin] = (double*)malloc(sizeof(double)*List_YOUSO[1]); 
  }

  UG = (double***)malloc(sizeof(double**)*List_YOUSO[3]); 
  for (i=0; i<List_YOUSO[3]; i++){
    UG[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
    for (j=0; j<List_YOUSO[7]; j++){
      UG[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
    }
  }

  UGUt = (double***)malloc(sizeof(double**)*List_YOUSO[8]); 
  for (i=0; i<List_YOUSO[8]; i++){
    UGUt[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
    for (j=0; j<List_YOUSO[7]; j++){
      UGUt[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
    }
  }

  /****************************************************
                      PrintMemory
  ****************************************************/

  if (firsttime) {
    PrintMemory("Calc_DM: UG",
                 sizeof(double)*List_YOUSO[3]*List_YOUSO[7]*List_YOUSO[7],NULL);
    PrintMemory("Calc_DM: UGUt",
                 sizeof(double)*List_YOUSO[8]*List_YOUSO[7]*List_YOUSO[7],NULL);

    firsttime=0;
  }

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    wan = WhatSpecies[ct_AN];
    tno = Spe_Total_CNO[wan];
    N = tno*(RNUM[ct_AN] + 1);

    for (spin=0; spin<=SpinP_switch; spin++){

      for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
        Gh_AN = natn[ct_AN][h_AN];
        hwan = WhatSpecies[Gh_AN];
        htno = Spe_Total_CNO[hwan];
        for (m=0; m<tno; m++){  
	  for (n=0; n<htno; n++){
	    CDM[spin][ct_AN][h_AN][m][n] = 0.0;
	    EDM[spin][ct_AN][h_AN][m][n] = 0.0;
	  }
        }
      }

      for (i=0; i<N; i++){

	/* U*G */

	for (rl1=0; rl1<=RNUM[ct_AN]; rl1++){
	  for (rl2=0; rl2<=RNUM[ct_AN]; rl2++){
	    for (m=0; m<tno; m++){
	      for (n=0; n<RCN_Rank[ct_AN][rl1]; n++){
		sum = 0.0;
		for (k=0; k<RCN_Rank[ct_AN][rl2]; k++){
		  sum = sum + LanU[spin][ct_AN][rl2][0][m][k]*
		                REVec[spin][ct_AN][i][rl2][k]*
		                REVec[spin][ct_AN][i][rl1][n];
		}

		if (rl2==0) UG[rl1][m][n] = sum;
		else        UG[rl1][m][n] = UG[rl1][m][n] + sum;
	      }
	    }
	  }                    
	}            

	/* U*G*U^t */

	for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
	  Gh_AN = natn[ct_AN][h_AN];
	  hwan = WhatSpecies[Gh_AN];
	  htno = Spe_Total_CNO[hwan];
	  for (rl2=0; rl2<=RNUM[ct_AN]; rl2++){
	    for (m=0; m<tno; m++){  
	      for (n=0; n<htno; n++){
		sum = 0.0;
		for (k=0; k<RCN_Rank[ct_AN][rl2]; k++){
		  sum = sum + UG[rl2][m][k]*LanU[spin][ct_AN][rl2][h_AN][n][k];
		}
		if (rl2==0) UGUt[h_AN][m][n] = sum;
		else        UGUt[h_AN][m][n] = UGUt[h_AN][m][n] + sum;
	      }
	    }            
	  }
	}

	/* Sum */

        x = (EigenValue[spin][ct_AN][i] - ChemP)*Beta;
        if (x<=-30.0) x = -30.0;
        if (30.0<=x)  x =  30.0;
        FermiF = 1.0/(1.0 + exp(x));

	for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
	  Gh_AN = natn[ct_AN][h_AN];
	  hwan = WhatSpecies[Gh_AN];
	  htno = Spe_Total_CNO[hwan];
          for (m=0; m<tno; m++){  
	    for (n=0; n<htno; n++){
              CDM[spin][ct_AN][h_AN][m][n] = CDM[spin][ct_AN][h_AN][m][n]
                                                + FermiF*UGUt[h_AN][m][n];
              EDM[spin][ct_AN][h_AN][m][n] = EDM[spin][ct_AN][h_AN][m][n]
                                                + FermiF*UGUt[h_AN][m][n]*
                                               EigenValue[spin][ct_AN][i];
	    }
	  }
	}
      }
    }
  }

  /****************************************************
               Averaging of CDM and EDM
  ****************************************************/

  for (spin=0; spin<=SpinP_switch; spin++){
    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      wan = WhatSpecies[ct_AN];
      tno1 = Spe_Total_CNO[wan] - 1;
      for (i=0; i<=FNAN[ct_AN]; i++){
        ig = natn[ct_AN][i];
        ino1 = Spe_Total_CNO[WhatSpecies[ig]] - 1;
        k = RMI1[ct_AN][i][0];
        for (p=0; p<=tno1; p++){
	  for (q=0; q<=ino1; q++){
	    dum1 = CDM[spin][ct_AN][i][p][q];
	    dum2 = CDM[spin][ig][k][q][p];
	    dum = 0.50*(dum1 + dum2);
	    CDM[spin][ct_AN][i][p][q] = dum;
	    CDM[spin][ig][k][q][p]    = dum;
	    dum1 = EDM[spin][ct_AN][i][p][q];
	    dum2 = EDM[spin][ig][k][q][p];
	    dum = 0.50*(dum1 + dum2);
	    EDM[spin][ct_AN][i][p][q] = dum;
	    EDM[spin][ig][k][q][p]    = dum;
	  }
        }
      } 
    }
  }

  /****************************************************
                   calc Uele_IS
  ****************************************************/

  Uele_IS[0] = 0.0;
  Uele_IS[1] = 0.0;
  
  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    wanA = WhatSpecies[ct_AN];
    tnoA = Spe_Total_CNO[wanA];
    for (j=0; j<=FNAN[ct_AN]; j++){
      wanB = WhatSpecies[natn[ct_AN][j]];
      tnoB = Spe_Total_CNO[wanB];
      for (k=0; k<tnoA; k++){
        for (l=0; l<tnoB; l++){
          for (spin=0; spin<=SpinP_switch; spin++){
            Uele_IS[spin] += CDM[spin][ct_AN][j][k][l]*Ham[spin][ct_AN][j][k][l];
          }
	}
      }
    }
  }

  /****************************************************
    freeing of arrays:

    static double Uele_temp[YOUSO23][YOUSO1];
    static double UG[YOUSO3][YOUSO7][YOUSO7];
    static double UGUt[YOUSO8][YOUSO7][YOUSO7];
  ****************************************************/

  for (spin=0; spin<=SpinP_switch; spin++){
    free(Uele_temp[spin]);
  }
  free(Uele_temp);

  for (i=0; i<List_YOUSO[3]; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(UG[i][j]);
    }
    free(UG[i]);
  }
  free(UG);

  for (i=0; i<List_YOUSO[8]; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(UGUt[i][j]);
    }
    free(UGUt[i]);
  }
  free(UGUt);

}

void Search_ChemP(double ****OLP0)
{
  static int firsttime=1;
  static int ct_AN,spin,wan,tno,N,i,j,po,m,loop_num;
  static int rl1,rl2,h_AN,k,Gh_AN,hwan,n,htno;
  static double ChemP_MAX,ChemP_MIN,Tnum,sum;
  static double TZ,Dnum,FermiF,x,Lnum,tmp;
  static double ***UG;
  static double ***UGUt;
  static double ***Nele_E;

  /****************************************************
    allocation of arrays:

    static double UG[YOUSO3][YOUSO7][YOUSO7];
    static double UGUt[YOUSO8][YOUSO7][YOUSO7];
    static double Nele_E[YOUSO23][YOUSO1][YOUSO7*YOUSO3];
  ****************************************************/

  UG = (double***)malloc(sizeof(double**)*List_YOUSO[3]); 
  for (i=0; i<List_YOUSO[3]; i++){
    UG[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
    for (j=0; j<List_YOUSO[7]; j++){
      UG[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
    }
  }

  UGUt = (double***)malloc(sizeof(double**)*List_YOUSO[8]); 
  for (i=0; i<List_YOUSO[8]; i++){
    UGUt[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
    for (j=0; j<List_YOUSO[7]; j++){
      UGUt[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
    }
  }

  Nele_E = (double***)malloc(sizeof(double**)*(SpinP_switch+1)); 
  for (spin=0; spin<=SpinP_switch; spin++){
    Nele_E[spin] = (double**)malloc(sizeof(double*)*List_YOUSO[1]); 
    for (i=0; i<List_YOUSO[1]; i++){
      Nele_E[spin][i] = (double*)malloc(sizeof(double)*List_YOUSO[7]*List_YOUSO[3]); 
    }
  }

  /****************************************************
                      PrintMemory
  ****************************************************/

  if (firsttime) {
    PrintMemory("Search_ChemP: UG",
                 sizeof(double)*List_YOUSO[3]*List_YOUSO[7]*List_YOUSO[7],NULL);
    PrintMemory("Search_ChemP: UGUt",
                 sizeof(double)*List_YOUSO[8]*List_YOUSO[7]*List_YOUSO[7],NULL);

    PrintMemory("Search_ChemP: Nele_E",
                 sizeof(double)*List_YOUSO[23]*List_YOUSO[1]*
                                List_YOUSO[7]*List_YOUSO[3],NULL);

    firsttime=0;
  }


  TZ = 0.0;
  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    wan = WhatSpecies[ct_AN];
    TZ = TZ + Spe_Core_Charge[wan];
  }

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    wan = WhatSpecies[ct_AN];
    tno = Spe_Total_CNO[wan];
    N = tno*(RNUM[ct_AN] + 1);

    for (spin=0; spin<=SpinP_switch; spin++){
      for (i=0; i<N; i++){

	/* U*G */
         
	for (rl1=0; rl1<=RNUM[ct_AN]; rl1++){
	  for (rl2=0; rl2<=RNUM[ct_AN]; rl2++){
	    for (m=0; m<tno; m++){  
	      for (n=0; n<RCN_Rank[ct_AN][rl1]; n++){
		sum = 0.0;
		for (k=0; k<RCN_Rank[ct_AN][rl2]; k++){
		  sum = sum + LanU[spin][ct_AN][rl2][0][m][k]*
 		                REVec[spin][ct_AN][i][rl2][k]*
		                REVec[spin][ct_AN][i][rl1][n];
		}

		if (rl2==0) UG[rl1][m][n] = sum;
		else        UG[rl1][m][n] = UG[rl1][m][n] + sum;
	      }
	    }
	  }                    
	}            

	/* U*G*U^t */

	for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
	  Gh_AN = natn[ct_AN][h_AN];
	  hwan = WhatSpecies[Gh_AN];
	  htno = Spe_Total_CNO[hwan];
	  for (rl2=0; rl2<=RNUM[ct_AN]; rl2++){
	    for (m=0; m<tno; m++){  
	      for (n=0; n<htno; n++){
		sum = 0.0;
		for (k=0; k<RCN_Rank[ct_AN][rl2]; k++){
		  sum = sum + UG[rl2][m][k]*LanU[spin][ct_AN][rl2][h_AN][n][k];
		}
		if (rl2==0) UGUt[h_AN][m][n] = sum;
		else        UGUt[h_AN][m][n] = UGUt[h_AN][m][n] + sum;
	      }
	    }            
	  }
	}

	/* U*G*U^t*OLP0 */

	for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
	  Gh_AN = natn[ct_AN][h_AN];
	  hwan = WhatSpecies[Gh_AN];
	  htno = Spe_Total_CNO[hwan];

          sum = 0.0;
	  for (m=0; m<tno; m++){  
	    for (k=0; k<htno; k++){
	      sum = sum + UGUt[h_AN][m][k]*OLP0[ct_AN][h_AN][m][k];
	    }
	  }              

          if (h_AN==0) Nele_E[spin][ct_AN][i] = sum;
	  else         Nele_E[spin][ct_AN][i] = Nele_E[spin][ct_AN][i] + sum;
	}          
      } 
    }
  }

  po = 0;
  loop_num = 0;
  ChemP_MAX = 15.0;
  ChemP_MIN =-15.0;

  do {
    ChemP = 0.50*(ChemP_MAX + ChemP_MIN);
    Tnum = 0.0;

    loop_num++;

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      wan = WhatSpecies[ct_AN];
      tno = Spe_Total_CNO[wan];
      N = tno*(RNUM[ct_AN] + 1);

      Lnum = 0.0;
      for (spin=0; spin<=SpinP_switch; spin++){
        for (i=0; i<N; i++){
          x = (EigenValue[spin][ct_AN][i] - ChemP)*Beta;
          if (x<=-30.0) x = -30.0;
          if (30.0<=x)  x =  30.0;
          if (SpinP_switch==0) FermiF = 2.0/(1.0 + exp(x));
          else                 FermiF = 1.0/(1.0 + exp(x));

          Lnum = Lnum + FermiF*Nele_E[spin][ct_AN][i];
        } 
      }

      Tnum = Tnum + Lnum;

      /*
      printf("num=%2d  ct_AN=%2d  Lnum=%15.12f\n",loop_num,ct_AN,Lnum);
      */
    }

    Dnum = TZ - Tnum;

    if (2<=level_stdout){
      printf("  ChemP=%16.13f TZ=%16.13f Tnum=%16.13f\n",eV2Hartree*ChemP,TZ,Tnum);
    }

    if (0.0<=Dnum) ChemP_MIN = ChemP;
    else           ChemP_MAX = ChemP;
    if (fabs(Dnum)<1.0e-12) po = 1;

    if (100<loop_num) po = 1;

    if (po==1){

      for (spin=0; spin<=SpinP_switch; spin++) Uele_OS[spin] = 0.0;
      for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
        wan = WhatSpecies[ct_AN];
        tno = Spe_Total_CNO[wan];
        N = tno*(RNUM[ct_AN] + 1);

        for (spin=0; spin<=SpinP_switch; spin++){
          for (i=0; i<N; i++){
            x = (EigenValue[spin][ct_AN][i] - ChemP)*Beta;
            if (x<=-30.0) x = -30.0;
            if (30.0<=x)  x =  30.0;
            FermiF = 1.0/(1.0 + exp(x));
            tmp = FermiF*Nele_E[spin][ct_AN][i]*EigenValue[spin][ct_AN][i];
            Uele_OS[spin] = Uele_OS[spin] + tmp;       
	  }
	}
      }
    }

  }
  while (po==0); 

  /****************************************************
    freeing of arrays:

    static double UG[YOUSO3][YOUSO7][YOUSO7];
    static double UGUt[YOUSO8][YOUSO7][YOUSO7];
    static double Nele_E[YOUSO23][YOUSO1][YOUSO7*YOUSO3];
  ****************************************************/

  for (i=0; i<List_YOUSO[3]; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(UG[i][j]);
    }
    free(UG[i]);
  }
  free(UG);

  for (i=0; i<List_YOUSO[8]; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(UGUt[i][j]);
    }
    free(UGUt[i]);
  }
  free(UGUt);

  for (spin=0; spin<=SpinP_switch; spin++){
    for (i=0; i<List_YOUSO[1]; i++){
      free(Nele_E[spin][i]);
    }
    free(Nele_E[spin]);
  }
  free(Nele_E);

}


void Solve_EigenValue(int ct_AN, int spin)
{

  static int rl,i,j,i0,j0,ig,jg,k,m,wan,tno,tno1;
  static int ig2,jg2;

  static int k0,rl0,l0,l;

  static double tmp0,tmp1,tmp2;

  static double *A,*W,*WORK;
  static INTEGER N,LDA,LWORK,INFO;
  static char jobz = 'V';
  static char uplo = 'L';

  wan = WhatSpecies[ct_AN];
  tno = Spe_Total_CNO[wan];
  tno1 = tno - 1;
  N = tno*(RNUM[ct_AN] + 1);
  LDA = N;
  LWORK = 3*N;

  A    = (double*)malloc(sizeof(double)*N*N);
  W    = (double*)malloc(sizeof(double)*N);
  WORK = (double*)malloc(sizeof(double)*3*N);

  for (i=0; i<N*N; i++){
    A[i] = 0.0;
  }

  /* Setting of al */

  for (rl=0; rl<=RNUM[ct_AN]; rl++){
    i0 = rl*tno;
    j0 = rl*tno;
    for (i=0; i<=tno1; i++){
      for (j=0; j<=tno1; j++){
        ig = i0 + i;
        jg = j0 + j;
        k = N*jg + ig;
        A[k] = al[spin][ct_AN][rl][i][j]; 
        if (i==j && (RCN_Rank[ct_AN][rl]-1)<i)
          A[k] = 1.0e+4;
      }
    }
  }

  /* Setting of be */

  for (rl=1; rl<=RNUM[ct_AN]; rl++){
    i0 = rl*tno;
    j0 = (rl-1)*tno;
    for (i=0; i<=tno1; i++){
      for (j=0; j<=tno1; j++){
        ig = i0 + i;
        jg = j0 + j;
        k = N*jg + ig;
        A[k] = be[spin][ct_AN][rl][i][j]; 
      }
    }
  }
  
  /* Call dsbev_() */

  /*
  F77_NAME(dsyev,DSYEV)( &jobz, &uplo, &N, A, &LDA, W, WORK, &LWORK, &INFO);
  */

  if (INFO!=0){
    printf("Errors in dsbev_() info=%2d\n",INFO);
  }

  /*
  for (i=0; i<N; i++){
    printf("i=%2d  W=%18.15f\n",i,W[i]);
  }
  */

  for (i=0; i<N; i++){
    EigenValue[spin][ct_AN][i] = W[i];
  }

  for (j=0; j<N; j++){
    for (i=0; i<N; i++){
      k = N*j + i;
      m  = i%tno;
      rl = (i-m)/tno;
      REVec[spin][ct_AN][j][rl][m] = A[k];
    }
  }

  /* free arraies */

  free(A);  
  free(W);  
  free(WORK);  

}

void Solve_EigenValue2(int ct_AN, int spin)
{

  static int rl,i,j,i0,j0,ig,jg,k,m,wan,tno,tno1;
  static int ig2,jg2;

  static int k0,rl0,l0,l;

  static double tmp0,tmp1,tmp2;
  static double *AB,*W,*Z,*WORK;

  static INTEGER N,KD,LDAB,LDZ,INFO;
  static char jobz = 'V';
  static char uplo = 'L';

  wan = WhatSpecies[ct_AN];
  tno = Spe_Total_CNO[wan];
  tno1 = tno - 1;
  N = tno*(RNUM[ct_AN] + 1);

  KD = tno*2 - 1;
  LDAB = KD + 1;
  LDZ = N;

  AB = (double*)malloc(sizeof(double)*LDAB*(N+3));
  W = (double*)malloc(sizeof(double)*(N+3));
  Z = (double*)malloc(sizeof(double)*LDZ*(N+3));
  WORK = (double*)malloc(sizeof(double)*4*(N+3));

  for (i=0; i<LDAB*N; i++){
    AB[i] = 0.0;
  }

  /* Setting of al */

  for (rl=0; rl<=RNUM[ct_AN]; rl++){
    i0 = rl*tno;
    j0 = rl*tno;
    for (i=0; i<=tno1; i++){
      for (j=0; j<=tno1; j++){
        ig = i0 + i + 1;
        jg = j0 + j + 1;

        if (jg<=ig){
          ig2 = 1 + ig - jg - 1;
          jg2 = jg - 1;
          k = LDAB*jg2 + ig2;
          AB[k] = al[spin][ct_AN][rl][i][j]; 
          if (i==j && (RCN_Rank[ct_AN][rl]-1)<i) AB[k] = 1.0e+4;
	}

      }
    }
  }

  /* Setting of be */
  
  for (rl=1; rl<=RNUM[ct_AN]; rl++){
    i0 = rl*tno;
    j0 = (rl-1)*tno;
    for (i=0; i<=tno1; i++){
      for (j=0; j<=tno1; j++){
        ig = i0 + i + 1;
        jg = j0 + j + 1;
          
        ig2 = 1 + ig - jg - 1;
        jg2 = jg - 1;
        k = LDAB*jg2 + ig2;
        AB[k] = be[spin][ct_AN][rl][i][j]; 
      }
    }
  }
  
  /* Call dsbev_() */

  /*
  F77_NAME(dsbev,DSBEV)( &jobz, &uplo, &N, &KD, AB, &LDAB, W, Z, &LDZ, WORK, &INFO);
  */

  if (INFO!=0){
    printf("Errors in dsbev_() info=%2d\n",INFO);
  }

  /*
  for (i=0; i<N; i++){
    printf("ct_AN=%i i=%2d  W=%18.15f\n",ct_AN,i,W[i]);
  }
  */

  for (i=0; i<N; i++){
    EigenValue[spin][ct_AN][i] = W[i];
  }

  for (j=0; j<N; j++){
    for (i=0; i<N; i++){
      k = N*j + i;
      m  = i%tno;
      rl = (i-m)/tno;
      REVec[spin][ct_AN][j][rl][m] = Z[k];
    }
  }

  /* free arraies */

  free(AB);  
  free(W);  
  free(Z);  
  free(WORK);  

}


void SbyISH(double ****OLP0,  double ****HP,
            double ****nh, double ****DH)
{
  static int firsttime=1;
  static int ct_AN,fan,san,can,ig,ian,j,jg,jan;
  static int k,kg,kan,m,n,l,jl,size_AH;
  static int i,tno0,tno1,h_AN,Gh_AN,Cwan,Hwan;
  static double sum;
  static double ****AH;

  /****************************************************
    allocation of arrays:

    static double AH[YOUSO1][YOUSO2][YOUSO7][YOUSO7];
  ****************************************************/

  size_AH = 0;
  AH = (double****)malloc(sizeof(double***)*(atomnum+1)); 
  for (ct_AN=0; ct_AN<=atomnum; ct_AN++){

    if (ct_AN==0){
      tno0 = 1;
    }
    else{
      Cwan = WhatSpecies[ct_AN];
      tno0 = Spe_Total_NO[Cwan];  
    }    

    AH[ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN]+SNAN[ct_AN]+1));
    for (h_AN=0; h_AN<=(FNAN[ct_AN]+SNAN[ct_AN]); h_AN++){

      if (ct_AN==0){
	tno1 = 1;  
      }
      else{
	Gh_AN = natn[ct_AN][h_AN];
	Hwan = WhatSpecies[Gh_AN];
	tno1 = Spe_Total_NO[Hwan];
      } 

      AH[ct_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
      for (i=0; i<tno0; i++){
	AH[ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
      }
      size_AH += tno0*tno1; 
    }
  }

  /****************************************************
                      PrintMemory
  ****************************************************/

  if (firsttime) {
    PrintMemory("SbyISH: AH",sizeof(double)*size_AH,NULL);
    firsttime=0;
  }

  /****************************************************
                      S*HP -> AH
  ****************************************************/

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    fan = FNAN[ct_AN];
    san = SNAN[ct_AN];
    can = fan + san;
    ig = natn[ct_AN][0];
    ian = Spe_Total_CNO[WhatSpecies[ig]];
    for (j=0; j<=can; j++){
      jg = natn[ct_AN][j];
      jan = Spe_Total_CNO[WhatSpecies[jg]];
      for (k=0; k<=fan; k++){
        kg = natn[ct_AN][k];
        jl = RMI2[ct_AN][k][j];
        kan = Spe_Total_CNO[WhatSpecies[kg]];

        if (0<=jl){
	  for (m=0; m<=(ian-1); m++){
	    for (n=0; n<=(jan-1); n++){
	      sum = 0.0;
	      for (l=0; l<=(kan-1); l++){
		sum = sum + OLP0[ct_AN][k][m][l]*HP[kg][jl][l][n];
	      }
	      if (k==0){
		AH[ct_AN][j][m][n] = sum;
	      }
	      else {
		AH[ct_AN][j][m][n] = AH[ct_AN][j][m][n] + sum;
	      } 
	    }
	  }          
	}

      }
    }
  }

  /****************************************************
                        DH = AH
  ****************************************************/

  sum = 0.0;
  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    fan = FNAN[ct_AN];
    san = SNAN[ct_AN];
    can = fan + san;
    ig = natn[ct_AN][0];
    ian = Spe_Total_CNO[WhatSpecies[ig]];
    for (j=0; j<=fan; j++){
      jg = natn[ct_AN][j];
      jan = Spe_Total_CNO[WhatSpecies[jg]];
      for (m=0; m<=(ian-1); m++){
        for (n=0; n<=(jan-1); n++){
          DH[ct_AN][j][m][n] = AH[ct_AN][j][m][n] - nh[ct_AN][j][m][n];
          sum = sum + fabs(DH[ct_AN][j][m][n]);
        }
      }
    }
  }

  printf("sum=%15.12f\n",sum);

  /****************************************************
    freeing of arrays:

    static double AH[YOUSO1][YOUSO2][YOUSO7][YOUSO7];
  ****************************************************/

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
	free(AH[ct_AN][h_AN][i]);
      }
      free(AH[ct_AN][h_AN]);
    }
    free(AH[ct_AN]);
  }
  free(AH);

}


void ISH(double ****IOLP, double ****nh, double ****HP)
{
  static int i,j,k;
  static int Gc_AN,Lh_AN,Gh_AN,Li_AN,Gi_AN,Lh_AN2;
  static int wan,c_no,h_no,i_no;
  static double dum,sum;

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    wan = WhatSpecies[Gc_AN];
    c_no = Spe_Total_CNO[wan];
    for (Lh_AN=0; Lh_AN<=(FNAN[Gc_AN]+SNAN[Gc_AN]); Lh_AN++){
      Gh_AN = natn[Gc_AN][Lh_AN];
      wan = WhatSpecies[Gh_AN];
      h_no = Spe_Total_CNO[wan];
      for (i=0; i<=(c_no-1); i++){
       for (j=0; j<=(h_no-1); j++){
         HP[Gc_AN][Lh_AN][i][j] = 0.0; 
         HP[Gc_AN][Lh_AN][i][j] = 0.0; 
       }
      } 
      for (Li_AN=0; Li_AN<=(FNAN[Gc_AN]+SNAN[Gc_AN]); Li_AN++){
        Gi_AN = natn[Gc_AN][Li_AN];
        wan = WhatSpecies[Gi_AN];
        i_no = Spe_Total_CNO[wan];
        Lh_AN2 = RMI1[Gc_AN][Li_AN][Lh_AN];
        if (0<=Lh_AN2){
          for (i=0; i<=(c_no-1); i++){
            for (j=0; j<=(h_no-1); j++){
              sum = 0.0;
              for (k=0; k<=(i_no-1); k++){
                dum = IOLP[Gc_AN][Li_AN][i][k]*nh[Gi_AN][Lh_AN2][k][j];
                sum = sum + dum;
              }
              HP[Gc_AN][Lh_AN][i][j] = HP[Gc_AN][Lh_AN][i][j] + sum;       
            }
          } 
        }
      }    
    }
  }

}


void B_Lanczos(int ct_AN, int spin, 
               double ****hopH, double ****OLP0, double ****HP)
{
  static int firsttime=1;
  static int i,j,k,l,m,n,ie,i1,j1,po,ig,jg,Rni,Rnj,Rn;
  static int kl,l1,l2,l3,ian,jan;
  static int m1,m2,m3,ZeroNum,Rloop_END;
  static int fan,san,can,rl,rl0,wan,ct_on;
  static int tno0,tno1,Cwan,Hwan,h_AN,Gh_AN;
  static double sum,sum1,sum2,dum,dum1,sumx,sumy,sumz,xn,xa;
  static double tmp0,tmp1;
  static double **BeGa;
  static double **UmRn;
  static double ***RCpre;
  static double ***RCcnt;
  static double ***RRcnt;
  static double ***RRcnt2;
  static double ***RR0;
  static double ****RU;

  /****************************************************
    allocation of arrays:

    static double BeGa[YOUSO7][YOUSO7];
    static double UmRn[YOUSO7][YOUSO7];
    static double RCpre[YOUSO2+1][YOUSO7][YOUSO7];
    static double RCcnt[YOUSO2+1][YOUSO7][YOUSO7];
    static double RRcnt[YOUSO2+1][YOUSO7][YOUSO7];
    static double RRcnt2[YOUSO2+1][YOUSO7][YOUSO7];
    static double RR0[YOUSO2+1][YOUSO7][YOUSO7];
    static double RU[YOUSO3+1][YOUSO2][YOUSO7][YOUSO7];
  ****************************************************/

  /* BeGa */
  BeGa = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    BeGa[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
  }

  /* UmRn */
  UmRn = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    UmRn[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
  }

  /* RCpre */
  RCpre = (double***)malloc(sizeof(double**)*(List_YOUSO[2]+1)); 
  for (i=0; i<(List_YOUSO[2]+1); i++){
    RCpre[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
    for (j=0; j<List_YOUSO[7]; j++){
      RCpre[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
    }
  }

  /* RCcnt */
  RCcnt = (double***)malloc(sizeof(double**)*(List_YOUSO[2]+1)); 
  for (i=0; i<(List_YOUSO[2]+1); i++){
    RCcnt[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
    for (j=0; j<List_YOUSO[7]; j++){
      RCcnt[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
    }
  }

  /* RRcnt */
  RRcnt = (double***)malloc(sizeof(double**)*(List_YOUSO[2]+1)); 
  for (i=0; i<(List_YOUSO[2]+1); i++){
    RRcnt[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
    for (j=0; j<List_YOUSO[7]; j++){
      RRcnt[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
    }
  }

  /* RRcnt2 */
  RRcnt2 = (double***)malloc(sizeof(double**)*(List_YOUSO[2]+1)); 
  for (i=0; i<(List_YOUSO[2]+1); i++){
    RRcnt2[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
    for (j=0; j<List_YOUSO[7]; j++){
      RRcnt2[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
    }
  }

  /* RR0 */
  RR0 = (double***)malloc(sizeof(double**)*(List_YOUSO[2]+1)); 
  for (i=0; i<(List_YOUSO[2]+1); i++){
    RR0[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
    for (j=0; j<List_YOUSO[7]; j++){
      RR0[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
    }
  }

  /* RU */
  RU = (double****)malloc(sizeof(double***)*(List_YOUSO[3]+1)); 
  for (i=0; i<(List_YOUSO[3]+1); i++){
    RU[i] = (double***)malloc(sizeof(double**)*List_YOUSO[2]); 
    for (j=0; j<List_YOUSO[2]; j++){
      RU[i][j] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
      for (k=0; k<List_YOUSO[7]; k++){
        RU[i][j][k] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
      }
    }
  }

  /****************************************************
                       PrintMemory
  ****************************************************/

  if (firsttime) {
    PrintMemory("B_Lanczos: RCpre",
                 sizeof(double)*(List_YOUSO[2]+1)*List_YOUSO[7]*List_YOUSO[7],NULL);
    PrintMemory("B_Lanczos: RCcnt",
                 sizeof(double)*(List_YOUSO[2]+1)*List_YOUSO[7]*List_YOUSO[7],NULL);
    PrintMemory("B_Lanczos: RRcnt",
                 sizeof(double)*(List_YOUSO[2]+1)*List_YOUSO[7]*List_YOUSO[7],NULL);
    PrintMemory("B_Lanczos: RRcnt2",
                 sizeof(double)*(List_YOUSO[2]+1)*List_YOUSO[7]*List_YOUSO[7],NULL);
    PrintMemory("B_Lanczos: RR0",
                 sizeof(double)*(List_YOUSO[2]+1)*List_YOUSO[7]*List_YOUSO[7],NULL);
    PrintMemory("B_Lanczos: RU",
                 sizeof(double)*(List_YOUSO[3]+1)*List_YOUSO[2]*
                                 List_YOUSO[7]*List_YOUSO[7],NULL);
    firsttime=0;
  }


  fan = FNAN[ct_AN];
  san = SNAN[ct_AN];
  can = fan + san;
  wan = WhatSpecies[ct_AN];
  ct_on = Spe_Total_CNO[wan];

  for (i=0; i<=can; i++){
    for (i1=0; i1<List_YOUSO[7]; i1++){
      for (j1=0; j1<List_YOUSO[7]; j1++){
	RCpre[i][i1][j1] = 0.0;
	RCcnt[i][i1][j1] = 0.0;
	RRcnt[i][i1][j1] = 0.0;
      }
    }
  }

  /****************************************************
           Consider LanU as column vectors.
  ****************************************************/

  for (i=0; i<=can; i++){
    for (j=0; j<=rlmax; j++){
      for (i1=0; i1<ct_on; i1++){
	for (j1=0; j1<ct_on; j1++){
	  RU[j][i][i1][j1] = 0.0;
	}
      } 
    }
  }

  for (i=0; i<ct_on; i++){
    RCcnt[0][i][i] = 1.0;
    LanU[spin][ct_AN][0][0][i][i] = 1.0;
    RU[0][0][i][i] = 1.0;
  }

  for (i=0; i<ct_on; i++){
    for (j=0; j<ct_on; j++){
      be[spin][ct_AN][0][i][j] = 0.0;
      ga[spin][ct_AN][0][i][j] = 0.0;
    }
  }

  /****************************************************
     Recursion of one-sided block Lanczos algorithm
  ****************************************************/

  n = 0; 
  for (i=1; i<=can; i++){
    ig = natn[ct_AN][i];
    ian = Spe_Total_CNO[WhatSpecies[ig]];
    n = n + ian;        
  }

  if (n%ct_on==0) RNUM[ct_AN] = n/ct_on;
  else            RNUM[ct_AN] = n/ct_on + 1;
  if (rlmax<RNUM[ct_AN]) RNUM[ct_AN] = rlmax;

  Rloop_END = 0;
  RCN_Rank[ct_AN][0] = ct_on;

  rl = 0;
  do {

    /****************************************************
      H'|WRn} (=RCcnt[j][m][n])
       and
      H|WRn}  (=RCcnt2[j][m][n])
    
      RCcnt[j][m][n] and RCcnt2[j][m][n]
        are a column vector.
    
      j  ->  local index for atom
      m  ->  index for basis
      n  ->  index for basis
    ****************************************************/

    /* H'|WRn} */
 
    for (i=0; i<=can; i++){
      ig = natn[ct_AN][i];
      ian = Spe_Total_CNO[WhatSpecies[ig]];
      po = 0;
      for (j=0; j<=can; j++){
        kl = RMI2[ct_AN][i][j];
        jg = natn[ct_AN][j];
        jan = Spe_Total_CNO[WhatSpecies[jg]];

        if (0<=kl){
	  for (m=0; m<=(ian-1); m++){
	    for (n=0; n<RCN_Rank[ct_AN][rl]; n++){
	      sum = 0.0;
	      for (k=0; k<=(jan-1); k++){
		sum = sum + HP[ig][kl][m][k]*RCcnt[j][k][n];
	      }
	      if (po==0)
		RRcnt[i][m][n] = sum;
	      else
		RRcnt[i][m][n] = RRcnt[i][m][n] + sum;
	    }
	  } 
	  po++;
	}
      }
    } 

    /* H|WRn} */

    for (i=0; i<=can; i++){
      ig = natn[ct_AN][i];
      ian = Spe_Total_CNO[WhatSpecies[ig]];
      po = 0;
      for (j=0; j<=can; j++){
        kl = RMI1[ct_AN][i][j];
        jg = natn[ct_AN][j];
        jan = Spe_Total_CNO[WhatSpecies[jg]];
        if (0<=kl){
	  for (m=0; m<=(ian-1); m++){
	    for (n=0; n<RCN_Rank[ct_AN][rl]; n++){
	      sum = 0.0;
	      for (k=0; k<=(jan-1); k++){
		sum = sum + hopH[ig][kl][m][k]*RCcnt[j][k][n];
	      }
	      if (po==0)
		RRcnt2[i][m][n] = sum;
	      else
		RRcnt2[i][m][n] = RRcnt2[i][m][n] + sum;
	    }
	  } 
	  po++;
	}
      }
    }
    
    /****************************************************
       Alpha_n = {WRn|H|WRn}

                Alpha_n is calculated as
              {WRn|H|WRn} = <RCcnt|RRcnt2>.
    ****************************************************/

    for (m=0; m<ct_on; m++){
      for (n=0; n<ct_on; n++){
        al[spin][ct_AN][rl][m][n] = 0.0;
      }
    }

    for (i=0; i<=can; i++){
      ig = natn[ct_AN][i];
      ian = Spe_Total_CNO[WhatSpecies[ig]];
      for (m=0; m<RCN_Rank[ct_AN][rl]; m++){
        for (n=0; n<RCN_Rank[ct_AN][rl]; n++){
          sum = 0.0;
          for (k=0; k<=(ian-1); k++){
            sum = sum + RCcnt[i][k][m]*RRcnt2[i][k][n];
          }
          if (i==0)
            al[spin][ct_AN][rl][m][n] = sum;
          else
            al[spin][ct_AN][rl][m][n] = al[spin][ct_AN][rl][m][n] + sum;
        }
      }
    }

    for (m=0; m<RCN_Rank[ct_AN][rl]; m++){
      for (n=m+1; n<RCN_Rank[ct_AN][rl]; n++){
	tmp0 = al[spin][ct_AN][rl][m][n];
	tmp1 = al[spin][ct_AN][rl][n][m];
	al[spin][ct_AN][rl][n][m] = 0.5*(tmp0 + tmp1); 
	al[spin][ct_AN][rl][m][n] = 0.5*(tmp0 + tmp1); 
      }
    }

    if (rl!=RNUM[ct_AN]){

      /****************************************************
             |RRcnt} = H'|WRn} - |WRn-1}(Bn)^t - |WRn}An
      ****************************************************/

      for (i=0; i<=can; i++){
	ig = natn[ct_AN][i];
	ian = Spe_Total_CNO[WhatSpecies[ig]];    
	for (m=0; m<=(ian-1); m++){
	  for (n=0; n<RCN_Rank[ct_AN][rl]; n++){
	    sum = 0.0;
            if (rl!=0){
  	      for (k=0; k<RCN_Rank[ct_AN][rl-1]; k++){
	        sum = sum + RCpre[i][m][k]*be[spin][ct_AN][rl][n][k];
	      }
	    }
	    for (k=0; k<RCN_Rank[ct_AN][rl]; k++){
	      sum = sum + RCcnt[i][m][k]*al[spin][ct_AN][rl][k][n];
	    }
	    RRcnt[i][m][n] = RRcnt[i][m][n] - sum;
	  }
	}
      }

      /****************************************************
         Re-orthogonalization by the Gram-Schmidt method
      ****************************************************/

      /* |RRcnt) := |RRcnt) - (|RR0) := |RR0) + |RU_rl0)(RU_rl0|S|RRcnt)) */

      for (rl0=0; rl0<=rl; rl0++){

        /* (RU_rl0|S|RRcnt) */
    
        for (i=0; i<=can; i++){
          ig = natn[ct_AN][i];
          ian = Spe_Total_CNO[WhatSpecies[ig]];
          po = 0;
          for (j=0; j<=can; j++){
            kl = RMI1[ct_AN][i][j];
            jg = natn[ct_AN][j];
            jan = Spe_Total_CNO[WhatSpecies[jg]];
            if (0<=kl){
	      for (m=0; m<=(ian-1); m++){
		for (n=0; n<RCN_Rank[ct_AN][rl]; n++){
		  sum = 0.0;
		  for (k=0; k<=(jan-1); k++){
		    sum = sum + OLP0[ig][kl][m][k]*RRcnt[j][k][n];
		  }
		  if (po==0)
		    RRcnt2[i][m][n] = sum;
		  else
		    RRcnt2[i][m][n] = RRcnt2[i][m][n] + sum;
		}
	      } 
	      po++;
	    }
          }
        }

        for (i=0; i<=can; i++){
	  ig = natn[ct_AN][i];
	  ian = Spe_Total_CNO[WhatSpecies[ig]];
	  for (m=0; m<RCN_Rank[ct_AN][rl0]; m++){
	    for (n=0; n<RCN_Rank[ct_AN][rl]; n++){
	      sum = 0.0;
	      for (k=0; k<=(ian-1); k++){
		sum = sum + RU[rl0][i][k][m]*RRcnt2[i][k][n];
	      }
	      if (i==0)
		UmRn[m][n] = sum;
	      else
		UmRn[m][n] = UmRn[m][n] + sum;
	    }
	  }
	}

        /* |RR0) := |RR0) + |RU_rl0)(RU_rl0|S|RRcnt)) */
        
	for (i=0; i<=can; i++){
          ig = natn[ct_AN][i];
	  ian = Spe_Total_CNO[WhatSpecies[ig]];    
	  for (m=0; m<=(ian-1); m++){
	    for (n=0; n<RCN_Rank[ct_AN][rl]; n++){
	      sum = 0.0;
	      for (k=0; k<RCN_Rank[ct_AN][rl0]; k++){
		sum = sum + RU[rl0][i][m][k]*UmRn[k][n];
	      }
	      if (rl0==0) RR0[i][m][n] = sum;
              else        RR0[i][m][n] = RR0[i][m][n] + sum;
	    }
	  }
	}
      }

      /* |RRcnt) := |RRcnt) - |RR0) */

      for (i=0; i<=can; i++){
	ig = natn[ct_AN][i];
	ian = Spe_Total_CNO[WhatSpecies[ig]]; 
	for (m=0; m<=(ian-1); m++){
	  for (n=0; n<RCN_Rank[ct_AN][rl]; n++){
            RRcnt[i][m][n] = RRcnt[i][m][n] - RR0[i][m][n];           
	  } 
	}
      }

      /****************************************************
                (Bn+1)^t * Bn+1 = {RRcnt|S|RRcnt}
      ****************************************************/

      for (i=0; i<=can; i++){
	ig = natn[ct_AN][i];
	ian = Spe_Total_CNO[WhatSpecies[ig]];
	po = 0;
	for (j=0; j<=can; j++){
	  kl = RMI1[ct_AN][i][j];
	  jg = natn[ct_AN][j];
	  jan = Spe_Total_CNO[WhatSpecies[jg]];
          if (0<=kl){
	    for (m=0; m<=(ian-1); m++){
	      for (n=0; n<RCN_Rank[ct_AN][rl]; n++){
		sum = 0.0;
		for (k=0; k<=(jan-1); k++){
		  sum = sum + OLP0[ig][kl][m][k]*RRcnt[j][k][n];
		}
		if (po==0)
		  RRcnt2[i][m][n] = sum;
		else
		  RRcnt2[i][m][n] = RRcnt2[i][m][n] + sum;
	      }
	    } 
	    po++;
	  }

	}
      }

      for (i=0; i<=can; i++){
	ig = natn[ct_AN][i];
	ian = Spe_Total_CNO[WhatSpecies[ig]];
	for (m=0; m<RCN_Rank[ct_AN][rl]; m++){
	  for (n=0; n<RCN_Rank[ct_AN][rl]; n++){
	    sum = 0.0;
	    for (k=0; k<=(ian-1); k++){
	      sum = sum + RRcnt[i][k][m]*RRcnt2[i][k][n];
	    }
	    if (i==0)
	      BeGa[m][n] = sum;
	    else
	      BeGa[m][n] = BeGa[m][n] + sum;
	  }
	}
      }

      for (m=0; m<RCN_Rank[ct_AN][rl]; m++){
        for (n=m+1; n<RCN_Rank[ct_AN][rl]; n++){
          tmp0 = BeGa[m][n];
          tmp1 = BeGa[n][m];
          BeGa[m][n] = 0.5*(tmp0 + tmp1); 
          BeGa[n][m] = 0.5*(tmp0 + tmp1); 
	}
      }

      /****************************************************
	1) Perform the singular value decomposition of 
           {LRcnt|RRcnt} as U*W*V^t.
        2) Then, B_{n+1} = UW^{1/2} and C_{n+1}=W^{1/2}V^t.
        3) It is easy to calculate the inverses of B and C.
      ****************************************************/

      printf("al  ct_AN=%2d spin=%2d  rl=%2d\n",ct_AN,spin,rl);
      for (m=0; m<ct_on; m++){
        for (n=0; n<ct_on; n++){
          printf("%15.12f ",al[spin][ct_AN][rl][m][n]);
        }
        printf("\n");
      }

      printf("BeGa  ct_AN=%2d  rl=%2d\n",ct_AN,rl);
      for (m=0; m<RCN_Rank[ct_AN][rl]; m++){
        for (n=0; n<RCN_Rank[ct_AN][rl]; n++){
          printf("%15.12f ",BeGa[m][n]);
	}
        printf("\n");
      }

      /* Singular value decomposition */

      printf("A1\n");

      RCN_Rank[ct_AN][rl+1] = SVD_fact_inv(
                                rl,RCN_Rank[ct_AN][rl]-1,BeGa,
                                ga[spin][ct_AN][rl+1],iga[spin][ct_AN][rl+1],
                                be[spin][ct_AN][rl+1],ibe[spin][ct_AN][rl+1]);

      printf("A2\n");

      if (RCN_Rank[ct_AN][rl+1]==0){
        RNUM[ct_AN] = rl;
        Rloop_END = 1;
      }
      else if (RCN_Rank[ct_AN][rl+1]==1 && ct_on!=1 && Rloop_END!=-1){
        Rloop_END = -1;
      }
      else if (Rloop_END==-1){
        RNUM[ct_AN] = rl;
        Rloop_END = 1;
      }

      /*
      printf("ct_AN=%i  rl=%i  RNUM=%i  RCN_Rank=%i Rloop_END=%i  ct_on=%i\n",
              ct_AN,rl,RNUM[ct_AN],RCN_Rank[ct_AN][rl+1],Rloop_END,ct_on);
      */

      /****************************************************
                     |WRn+1} = |RRcnt}IBn+1
      ****************************************************/

      printf("A3\n");

      for (i=0; i<=can; i++){
	ig = natn[ct_AN][i];
	ian = Spe_Total_CNO[WhatSpecies[ig]];
	for (m=0; m<=(ian-1); m++){
          for (n=0; n<RCN_Rank[ct_AN][rl]; n++){
            RCpre[i][m][n] = RCcnt[i][m][n];
	  }
          for (n=0; n<RCN_Rank[ct_AN][rl+1]; n++){
            sum = 0.0;
            for (k=0; k<RCN_Rank[ct_AN][rl]; k++){
              sum = sum + RRcnt[i][m][k]*ibe[spin][ct_AN][rl+1][k][n];
	    }
            RCcnt[i][m][n] = sum;
            RU[rl+1][i][m][n] = sum;

            if (i<=FNAN[ct_AN])
              LanU[spin][ct_AN][rl+1][i][m][n] = sum;
	  }
	}
      }

      printf("A4\n");

      /****************************************************
                    End of if (rl!=RNUM[ct_AN])
      ****************************************************/

    }
    rl++;

    /****************************************************
                       End of rl loop
    ****************************************************/

  } 
  while (rl<=RNUM[ct_AN] && Rloop_END<=0);

  if (2<=level_stdout){
    printf("  ct_AN=%2d  RNUM=%2d\n",ct_AN,RNUM[ct_AN]);
  }


  /****************************************************
    freeing of arrays:

    static double BeGa[YOUSO7][YOUSO7];
    static double UmRn[YOUSO7][YOUSO7];
    static double RCpre[YOUSO2+1][YOUSO7][YOUSO7];
    static double RCcnt[YOUSO2+1][YOUSO7][YOUSO7];
    static double RRcnt[YOUSO2+1][YOUSO7][YOUSO7];
    static double RRcnt2[YOUSO2+1][YOUSO7][YOUSO7];
    static double RR0[YOUSO2+1][YOUSO7][YOUSO7];
    static double RU[YOUSO3+1][YOUSO2][YOUSO7][YOUSO7];
  ****************************************************/

  for (i=0; i<List_YOUSO[7]; i++){
    free(BeGa[i]);
  }
  free(BeGa);

  for (i=0; i<List_YOUSO[7]; i++){
    free(UmRn[i]);
  }
  free(UmRn);

  for (i=0; i<(List_YOUSO[2]+1); i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(RCpre[i][j]);
    }
    free(RCpre[i]);
  }
  free(RCpre);

  for (i=0; i<(List_YOUSO[2]+1); i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(RCcnt[i][j]);
    }
    free(RCcnt[i]);
  }
  free(RCcnt);

  for (i=0; i<(List_YOUSO[2]+1); i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(RRcnt[i][j]);
    }
    free(RRcnt[i]);
  }
  free(RRcnt);

  for (i=0; i<(List_YOUSO[2]+1); i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(RRcnt2[i][j]);
    }
    free(RRcnt2[i]);
  }
  free(RRcnt2);

  for (i=0; i<(List_YOUSO[2]+1); i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(RR0[i][j]);
    }
    free(RR0[i]);
  }
  free(RR0);

  for (i=0; i<(List_YOUSO[3]+1); i++){
    for (j=0; j<List_YOUSO[2]; j++){
      for (k=0; k<List_YOUSO[7]; k++){
        free(RU[i][j][k]);
      }
      free(RU[i][j]);
    }
    free(RU[i]);
  }
  free(RU);

}


int SVD_fact_inv(int rl, int n, double **a,
                 double **B, double **IB,
                 double **C, double **IC)
{

  static int i,j,j1,k,rankN,po,num;
  static double sum,rank_criterion,rank_criterion2;
  static double *a1,*vt,*u,*w,*work,*w2,*iw2;
  static double **vt1,**u1;
  static char jobu  = 'A';
  static char jobvt = 'A';
  static INTEGER ROW,COL,lda,ldu,ldvt,lwork,info;

  /****************************************************
    allocation of arrays:

    static double a1[YOUSO7*YOUSO7];
    static double vt[YOUSO7*YOUSO7];
    static double u[YOUSO7*YOUSO7];
    static double work[YOUSO7*YOUSO7];
    static double w[YOUSO7+2];
    static double w2[YOUSO7+2];
    static double iw2[YOUSO7+2];
    static double vt1[YOUSO7][YOUSO7];
    static double u1[YOUSO7][YOUSO7];
  ****************************************************/

  /* a1 */
  a1 = (double*)malloc(sizeof(double)*(List_YOUSO[7]+2)*(List_YOUSO[7]+2)); 

  /* vt */
  vt = (double*)malloc(sizeof(double)*(List_YOUSO[7]+2)*(List_YOUSO[7]+2)); 

  /* u */
  u = (double*)malloc(sizeof(double)*(List_YOUSO[7]+2)*(List_YOUSO[7]+2)); 

  /* work */
  work = (double*)malloc(sizeof(double)*(List_YOUSO[7]+2)*(List_YOUSO[7]+2)); 

  /* w */
  w = (double*)malloc(sizeof(double)*(List_YOUSO[7]+2)); 

  /* w2 */
  w2 = (double*)malloc(sizeof(double)*(List_YOUSO[7]+2)); 

  /* iw2 */
  iw2 = (double*)malloc(sizeof(double)*(List_YOUSO[7]+2)); 

  /* vt1 */
  vt1 = (double**)malloc(sizeof(double*)*(List_YOUSO[7]+2)); 
  for (i=0; i<(List_YOUSO[7]+2); i++){
    vt1[i] = (double*)malloc(sizeof(double)*(List_YOUSO[7]+2)); 
  }

  /* u1 */
  u1 = (double**)malloc(sizeof(double*)*(List_YOUSO[7]+2)); 
  for (i=0; i<(List_YOUSO[7]+2); i++){
    u1[i] = (double*)malloc(sizeof(double)*(List_YOUSO[7]+2)); 
  }

  /****************************************************
                      From 0 to n
  ****************************************************/

  num = -1;  
  for (j=0; j<=n; j++){
    for (i=0; i<=n; i++){
      num++;
      a1[num] = a[i][j];
    }
  }

  /****************************************************
             singular value decomposition
  ****************************************************/

  ROW  = n + 1;
  COL  = n + 1;
  lda  = n + 1;
  ldu  = n + 1;
  ldvt = n + 1;
  lwork = 5*List_YOUSO[7];
  /*
  F77_NAME(dgesvd,DGESVD)(&jobu, &jobvt, &ROW, &COL, a1, &lda, w, u,
          &ldu, vt, &ldvt, work, &lwork, &info);
  */

  num = -1;  
  for (j=0; j<=n; j++){
    for (i=0; i<=n; i++){
      num++;
      u1[i][j]  =  u[num];
      vt1[i][j] = vt[num];
    }
  }

  /****************************************************
             calculate rank of the BeGa matrix
  ****************************************************/

  if (rl<5){
    rank_criterion  = 1.0e-6;
    rank_criterion2 = 1.0e-7;
  }
  else if (rl<10){
    rank_criterion  = 1.0e-5;
    rank_criterion2 = 1.0e-6;
  }
  else if (rl<13) {
    rank_criterion  = 1.0e-4;
    rank_criterion2 = 1.0e-5;
  }
  else {
    rank_criterion  = 1.0e-3;
    rank_criterion2 = 1.0e-2;
  }

  rankN = -1;
  po = 0;
  do{
    rankN++;
    if (  (fabs(w[rankN]/w[0])<rank_criterion
        && fabs(w[rankN])<100.0*rank_criterion)
        ||
         fabs(w[rankN])<rank_criterion2 && rankN==0) po = 1; 
  } while(rankN<n && po==0);
  if (po==0) rankN++;

  /*
  printf("rankN=%2d\n",rankN);
  */

  for (i=0; i<=n; i++){
    if (i<rankN){
       w2[i] = sqrt(w[i]);
      iw2[i] = 1.0/w2[i];
    }
    else{
       w2[i] = 0.0;
      iw2[i] = 0.0;
    }
  }

  /****************************************************
             calculate beta and inverse beta
  ****************************************************/

  for (i=0; i<=n; i++){
    for (j=0; j<=n; j++){
      B[i][j] = u1[i][j]*w2[j];
    }
  }

  for (i=0; i<=n; i++){
    for (j=0; j<=n; j++){
      IB[i][j] = u1[j][i]*iw2[i];
    }
  }

  /****************************************************
            calculate gamma and inverse gamma
  ****************************************************/

  for (i=0; i<=n; i++){
    for (j=0; j<=n; j++){
      C[i][j] = vt1[i][j]*w2[i];
    }
  }

  for (i=0; i<=n; i++){
    for (j=0; j<=n; j++){
      IC[i][j] = vt1[j][i]*iw2[j];
    }
  }

  /****************************************************
    freeing of arrays:

    static double a1[YOUSO7*YOUSO7];
    static double vt[YOUSO7*YOUSO7];
    static double u[YOUSO7*YOUSO7];
    static double work[YOUSO7*YOUSO7];
    static double w[YOUSO7+2];
    static double w2[YOUSO7+2];
    static double iw2[YOUSO7+2];
    static double vt1[YOUSO7][YOUSO7];
    static double u1[YOUSO7][YOUSO7];
  ****************************************************/

  /* a1 */
  free(a1);

  /* vt */
  free(vt);

  /* u */
  free(u);

  /* work */
  free(work);

  /* w */
  free(w);

  /* w2 */
  free(w2);

  /* iw2 */
  free(iw2);

  /* vt1 */
  for (i=0; i<(List_YOUSO[7]+2); i++){
    free(vt1[i]);
  }
  free(vt1);

  /* u1 */
  for (i=0; i<(List_YOUSO[7]+2); i++){
    free(u1[i]);
  }
  free(u1);

  /* return */
  return rankN;
}


void Save_Recursion()
{
  static char fileRN[YOUSO10] = ".rcn";
  static FILE *fp_Rcn;

  strcpy(fileRN,".rcn");
  fnjoint(filepath,filename,fileRN);
  if ((fp_Rcn = fopen(fileRN,"w")) != NULL){
    Output_RcnCof(fp_Rcn);
    fclose(fp_Rcn);
  }
  else
    printf("Failure of saving the recusion logfile.\n");
}


void Output_RcnCof(FILE *fp_Rcn)
{
  static int spin,ct_AN,rl,i,j,wan1,TNO1;

  for (spin=0; spin<=SpinP_switch; spin++){

    if (spin==0){
      fprintf(fp_Rcn,"**********************************\n");
      fprintf(fp_Rcn,"************ Alpha UP ************\n");
      fprintf(fp_Rcn,"**********************************\n");
    } 
    else if (spin==1){
      fprintf(fp_Rcn,"**********************************\n");
      fprintf(fp_Rcn,"*********** Alpha DOWN ***********\n");
      fprintf(fp_Rcn,"**********************************\n");
    }

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      wan1 = WhatSpecies[ct_AN];
      TNO1 = Spe_Total_CNO[wan1];
      fprintf(fp_Rcn,"%i\n",ct_AN);
      for (rl=0; rl<=RNUM[ct_AN]; rl++){
        for (i=0; i<TNO1; i++){
          for (j=0; j<TNO1; j++){
            fprintf(fp_Rcn,"%19.14f ",al[spin][ct_AN][rl][i][j]);
          }
          fprintf(fp_Rcn,"\n");
        } 
      }
    }

    if (spin==0){
      fprintf(fp_Rcn,"**********************************\n");
      fprintf(fp_Rcn,"************ Beta UP *************\n");
      fprintf(fp_Rcn,"**********************************\n");
    }
    else if (spin==1){
      fprintf(fp_Rcn,"**********************************\n");
      fprintf(fp_Rcn,"*********** Beta DOWN ************\n");
      fprintf(fp_Rcn,"**********************************\n");
    }

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      wan1 = WhatSpecies[ct_AN];
      TNO1 = Spe_Total_CNO[wan1];
      fprintf(fp_Rcn,"%i\n",ct_AN);
      for (rl=0; rl<=RNUM[ct_AN]; rl++){
        for (i=0; i<TNO1; i++){
          for (j=0; j<TNO1; j++){
            fprintf(fp_Rcn,"%19.14f ",be[spin][ct_AN][rl][i][j]);
          }
          fprintf(fp_Rcn,"\n");
        } 
      }
    }

    if (spin==0){
      fprintf(fp_Rcn,"**********************************\n");
      fprintf(fp_Rcn,"******* Inverse of Beta UP *******\n");
      fprintf(fp_Rcn,"**********************************\n");
    }
    else if (spin==1){
      fprintf(fp_Rcn,"**********************************\n");
      fprintf(fp_Rcn,"****** Inverse of Beta DOWN ******\n");
      fprintf(fp_Rcn,"**********************************\n");
    }

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      wan1 = WhatSpecies[ct_AN];
      TNO1 = Spe_Total_CNO[wan1];
      fprintf(fp_Rcn,"%i\n",ct_AN);
      for (rl=0; rl<=RNUM[ct_AN]; rl++){
        for (i=0; i<TNO1; i++){
          for (j=0; j<TNO1; j++){
            fprintf(fp_Rcn,"%19.14f ",ibe[spin][ct_AN][rl][i][j]);
          }
          fprintf(fp_Rcn,"\n");
        } 
      }
    }

    if (spin==0){
      fprintf(fp_Rcn,"**********************************\n");
      fprintf(fp_Rcn,"************* Gamma UP ***********\n");
      fprintf(fp_Rcn,"**********************************\n");
    }
    else if (spin==1){
      fprintf(fp_Rcn,"**********************************\n");
      fprintf(fp_Rcn,"************ Gamma DOWN **********\n");
      fprintf(fp_Rcn,"**********************************\n");
    }


    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      wan1 = WhatSpecies[ct_AN];
      TNO1 = Spe_Total_CNO[wan1];
      fprintf(fp_Rcn,"%i\n",ct_AN);
      for (rl=0; rl<=RNUM[ct_AN]; rl++){
        for (i=0; i<TNO1; i++){
          for (j=0; j<TNO1; j++){
            fprintf(fp_Rcn,"%19.14f ",ga[spin][ct_AN][rl][i][j]);
          }
          fprintf(fp_Rcn,"\n");
        } 
      }
    }

    if (spin==0){
      fprintf(fp_Rcn,"**********************************\n");
      fprintf(fp_Rcn,"****** Inverse of Gamma UP *******\n");
      fprintf(fp_Rcn,"**********************************\n");
    }
    else if (spin==1){
      fprintf(fp_Rcn,"**********************************\n");
      fprintf(fp_Rcn,"***** Inverse of Gamma DOWN ******\n");
      fprintf(fp_Rcn,"**********************************\n");
    }

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      wan1 = WhatSpecies[ct_AN];
      TNO1 = Spe_Total_CNO[wan1];
      fprintf(fp_Rcn,"%i\n",ct_AN);
      for (rl=0; rl<=RNUM[ct_AN]; rl++){
        for (i=0; i<TNO1; i++){
          for (j=0; j<TNO1; j++){
            fprintf(fp_Rcn,"%19.14f ",iga[spin][ct_AN][rl][i][j]);
          }
          fprintf(fp_Rcn,"\n");
        } 
      }
    }

  }

}
