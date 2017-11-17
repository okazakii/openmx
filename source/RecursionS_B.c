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
#include <string.h>
#include <math.h>
#include <time.h>
#include "openmx_common.h"
#include "lapack_prototypes.h"

static void TSB_Lanczos(int ct_AN, int spin, double ****HP);
static void T_Conserve_MP(int T_switch);
static void BO_Calc_S(int T_switch, double *****CDM, double *****EDM, double ****IOLP);
static void GBTD00(int T_switch, int ct_AN, int spin,
		   dcomplex EpP, dcomplex **G00S);
static void RecurG(int ct_AN, int spin, dcomplex EpP, dcomplex **G00S, dcomplex ***G0S);

static int SVD_fact_inv(int rl, int n, double **a,
                        double **B, double **IB,
                        double **C, double **IC);

static void ISH(double ****IOLP, double ****nh, double ****HP);
static void SbyISH(double ****S, double ****HP,
		   double ****nh, double ****DH);
static void Save_Recursion();
static void Output_RcnCof(FILE *fp);
static void Output_LU(FILE *fp);

static double *****HPS;
static double *****DH;
static double *****al;
static double *****be;
static double *****ibe;
static double *****ga;
static double *****iga;
static double ******LU;
static int **RCN_Rank;
static double Uele_OS[4],Uele_IS[4];











double RecursionS_B(int SCF_iter,
                    double *****Ham,
                    double ****S,
                    double *****CDM,
                    double *****EDM,
                    double Eele0[2], double Eele1[2])
{
  static int firsttime=1;
  static int i,ct_AN,spin;
  static double time0;
  static int m,n,j,Rn,k,l,rl;
  static double sum;
  static double TStime,TEtime;
  static int size_HPS;
  static int tno0,tno1,Cwan,h_AN; 
  static int Gh_AN,Hwan;

  dtime(&TStime);
  Timetool_times("RecursionS_B","start");

  /****************************************************
    allocation of arrays:
  allocate_arrays_recursion();
  ****************************************************/

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


  al = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
  for (i=0; i<=SpinP_switch; i++){
    al[i] = (double****)malloc(sizeof(double***)*(atomnum+1)); 
    for (j=0; j<=atomnum; j++){
      al[i][j] = (double***)malloc(sizeof(double**)*List_YOUSO[3]); 
      for (k=0; k<List_YOUSO[3]; k++){
        al[i][j][k] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
        for (l=0; l<List_YOUSO[7]; l++){
          al[i][j][k][l] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
	}
      }
    }
  }

  be = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
  for (i=0; i<=SpinP_switch; i++){
    be[i] = (double****)malloc(sizeof(double***)*(atomnum+1)); 
    for (j=0; j<=atomnum; j++){
      be[i][j] = (double***)malloc(sizeof(double**)*List_YOUSO[3]); 
      for (k=0; k<List_YOUSO[3]; k++){
        be[i][j][k] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
        for (l=0; l<List_YOUSO[7]; l++){
          be[i][j][k][l] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
	}
      }
    }
  }

  ibe = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
  for (i=0; i<=SpinP_switch; i++){
    ibe[i] = (double****)malloc(sizeof(double***)*(atomnum+1)); 
    for (j=0; j<=atomnum; j++){
      ibe[i][j] = (double***)malloc(sizeof(double**)*List_YOUSO[3]); 
      for (k=0; k<List_YOUSO[3]; k++){
        ibe[i][j][k] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
        for (l=0; l<List_YOUSO[7]; l++){
          ibe[i][j][k][l] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
	}
      }
    }
  }

  ga = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
  for (i=0; i<=SpinP_switch; i++){
    ga[i] = (double****)malloc(sizeof(double***)*(atomnum+1)); 
    for (j=0; j<=atomnum; j++){
      ga[i][j] = (double***)malloc(sizeof(double**)*List_YOUSO[3]); 
      for (k=0; k<List_YOUSO[3]; k++){
        ga[i][j][k] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
        for (l=0; l<List_YOUSO[7]; l++){
          ga[i][j][k][l] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
	}
      }
    }
  }

  iga = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
  for (i=0; i<=SpinP_switch; i++){
    iga[i] = (double****)malloc(sizeof(double***)*(atomnum+1)); 
    for (j=0; j<=atomnum; j++){
      iga[i][j] = (double***)malloc(sizeof(double**)*List_YOUSO[3]); 
      for (k=0; k<List_YOUSO[3]; k++){
        iga[i][j][k] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
        for (l=0; l<List_YOUSO[7]; l++){
          iga[i][j][k][l] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
	}
      }
    }
  }

  LU = (double******)malloc(sizeof(double*****)*(SpinP_switch+1)); 
  for (i=0; i<=SpinP_switch; i++){
    LU[i] = (double*****)malloc(sizeof(double****)*(atomnum+1)); 
    for (j=0; j<=atomnum; j++){
      LU[i][j] = (double****)malloc(sizeof(double***)*List_YOUSO[3]); 
      for (k=0; k<List_YOUSO[3]; k++){
        LU[i][j][k] = (double***)malloc(sizeof(double**)*List_YOUSO[2]); 
        for (l=0; l<List_YOUSO[2]; l++){
          LU[i][j][k][l] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
          for (m=0; m<List_YOUSO[7]; m++){
            LU[i][j][k][l][m] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
	  }
	}
      }
    }
  }

  RCN_Rank = (int**)malloc(sizeof(int*)*(atomnum+1)); 
  for (i=0; i<=atomnum; i++){
    RCN_Rank[i] = (int*)malloc(sizeof(int)*List_YOUSO[3]); 
  }



  if (firsttime) {
  PrintMemory("RecursionS_B: al",sizeof(al),NULL);
  PrintMemory("RecursionS_B: be",sizeof(be),NULL);
  PrintMemory("RecursionS_B: ibe",sizeof(ibe),NULL);
  PrintMemory("RecursionS_B: ga",sizeof(ga),NULL);
  PrintMemory("RecursionS_B: iga",sizeof(iga),NULL);
  PrintMemory("RecursionS_B: LU",sizeof(LU),NULL);
  PrintMemory("RecursionS_B: IOLP",sizeof(IOLP),NULL);
  PrintMemory("RecursionS_B: HPS",sizeof(HPS),NULL);
  PrintMemory("RecursionS_B: DH",sizeof(DH),NULL);
  firsttime=0;
  }


  /****************************************************
     Calculate the inverse matrix of overlap matrix
  ****************************************************/
  
  printf("A1\n");

  if (IS_switch==1) IS_Lanczos(S,IOLP,rlmax_IS);
  if (IS_switch==2) IS_LU(S,IOLP);
  if (IS_switch==3) IS_Hotelling(S,IOLP,rlmax_IS);
  if (IS_switch==4) IS_Taylor(S,IOLP,rlmax_IS);

  /****************************************************
                   Calculate S^{-1}H
  ****************************************************/

  printf("A2\n");
  
  for (spin=0; spin<=SpinP_switch; spin++){
    ISH(IOLP,Ham[spin],HPS[spin]);
  }

  /****************************************************
          Two-sided block Lanczos transformation
  ****************************************************/





  printf("A3\n");

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    for (spin=0; spin<=SpinP_switch; spin++){





  printf("A4 ct_AN=%2d spin=%2d\n",ct_AN,spin);

      TSB_Lanczos(ct_AN,spin,HPS[spin]);
    }
  }

  /****************************************************
          Search the chemical potential to maintain
                the number of electrons
  ****************************************************/
  
  T_Conserve_MP(T_switch);
  
  /****************************************************
             Calculate the density matrix
  ****************************************************/

  BO_Calc_S(T_switch,CDM,EDM,IOLP);



  printf("Y1\n");

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
  

  printf("Y2\n");

  /****************************************************
   Output several informations of recursion algorithm
  ****************************************************/

  
  if (2<=level_fileout) Save_Recursion();


  printf("Y3\n");

  /****************************************************
    freeing of arrays:
  ****************************************************/


  /* HPS */  

  size_HPS = 0;
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

  printf("Y4\n");

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

  printf("Y5\n");


  for (i=0; i<=SpinP_switch; i++){
    for (j=0; j<=atomnum; j++){
      for (k=0; k<List_YOUSO[3]; k++){
        for (l=0; l<List_YOUSO[7]; l++){
          free(al[i][j][k][l]);
	}
        free(al[i][j][k]);
      }
      free(al[i][j]);
    }
    free(al[i]);
  }
  free(al);

  printf("Y6\n");


  for (i=0; i<=SpinP_switch; i++){
    for (j=0; j<=atomnum; j++){
      for (k=0; k<List_YOUSO[3]; k++){
        for (l=0; l<List_YOUSO[7]; l++){
          free(be[i][j][k][l]);
	}
        free(be[i][j][k]);
      }
      free(be[i][j]);
    }
    free(be[i]);
  }
  free(be);

  for (i=0; i<=SpinP_switch; i++){
    for (j=0; j<=atomnum; j++){
      for (k=0; k<List_YOUSO[3]; k++){
        for (l=0; l<List_YOUSO[7]; l++){
          free(ibe[i][j][k][l]);
	}
        free(ibe[i][j][k]);
      }
      free(ibe[i][j]);
    }
    free(ibe[i]);
  }
  free(ibe);

  for (i=0; i<=SpinP_switch; i++){
    for (j=0; j<=atomnum; j++){
      for (k=0; k<List_YOUSO[3]; k++){
        for (l=0; l<List_YOUSO[7]; l++){
          free(ga[i][j][k][l]);
	}
        free(ga[i][j][k]);
      }
      free(ga[i][j]);
    }
    free(ga[i]);
  }
  free(ga);

  for (i=0; i<=SpinP_switch; i++){
    for (j=0; j<=atomnum; j++){
      for (k=0; k<List_YOUSO[3]; k++){
        for (l=0; l<List_YOUSO[7]; l++){
          free(iga[i][j][k][l]);
	}
        free(iga[i][j][k]);
      }
      free(iga[i][j]);
    }
    free(iga[i]);
  }
  free(iga);

  for (i=0; i<=SpinP_switch; i++){
    for (j=0; j<=atomnum; j++){
      for (k=0; k<List_YOUSO[3]; k++){
        for (l=0; l<List_YOUSO[2]; l++){
          for (m=0; m<List_YOUSO[7]; m++){
            free(LU[i][j][k][l][m]);
	  }
          free(LU[i][j][k][l]);
	}
        free(LU[i][j][k]);
      }
      free(LU[i][j]);
    }
    free(LU[i]);
  }
  free(LU);

  for (i=0; i<=atomnum; i++){
    free(RCN_Rank[i]);
  }
  free(RCN_Rank);
  
  dtime(&TEtime);
  Timetool_times("RecursionS_B","end");
  time0 = TEtime - TStime;
  return time0;
}


void SbyISH(double ****S,
            double ****HP,
            double ****nh,
            double ****DH)
{
  static int firsttime=1;
  static int ct_AN,fan,san,can,ig,ian,j,jg,jan;
  static int k,kg,kan,m,n,l,jl,i;
  static double sum;
  static double ****AH;

  AH = (double****)malloc(sizeof(double***)*(atomnum+1));
  for (i=0; i<=atomnum; i++){
    AH[i] = (double***)malloc(sizeof(double**)*List_YOUSO[2]);
    for (j=0; j<List_YOUSO[2]; j++){
      AH[i][j] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
      for (k=0; k<List_YOUSO[7]; k++){
        AH[i][j][k] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
      }
    }
  }



  if (firsttime) {
  PrintMemory("SbyISH: AH",sizeof(AH),NULL);
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
        for (m=0; m<=(ian-1); m++){
          for (n=0; n<=(jan-1); n++){
            sum = 0.0;
            for (l=0; l<=(kan-1); l++){
              sum = sum + S[ct_AN][k][m][l]*HP[kg][jl][l][n];
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

  /****************************************************
                        DH = AH
  ****************************************************/

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
        }
      }
    }
  }

  for (i=0; i<=atomnum; i++){
    for (j=0; j<List_YOUSO[2]; j++){
      for (k=0; k<List_YOUSO[7]; k++){
        free(AH[i][j][k]);
      }
      free(AH[i][j]);
    }
    free(AH[i]);
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








void TSB_Lanczos(int ct_AN, int spin, double ****HP)
{
  static int firsttime=1;
  static int i,j,k,l,m,n,ie,i1,j1,po,ig,jg,Rni,Rnj,Rn;
  static int kl,l1,l2,l3,ian,jan;
  static int m1,m2,m3,ZeroNum,Rloop_END;
  static int fan,san,can,rl,rl0,wan,ct_on;
  static double sum,sum1,sum2,dum,dum1,sumx,sumy,sumz,xn,xa;
  static double **BeGa;
  static double **UmRn;
  static double **RnUm;
  static double ***RCpre;
  static double ***RCcnt;
  static double ***RRcnt;
  static double ***LCpre;
  static double ***LCcnt;
  static double ***LRcnt;
  static double ***RR0;
  static double ***LR0;
  static double ****RU;
  static int size_BeGa,size_RCpre,size_RU;

  /****************************************************
    allocation of arrays:
  ****************************************************/

  printf("B1\n");

  size_BeGa = List_YOUSO[7]*List_YOUSO[7];

  BeGa = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    BeGa[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
  }

  UmRn = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    UmRn[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
  }

  RnUm = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    RnUm[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
  }

  size_RCpre = List_YOUSO[2]*List_YOUSO[7]*List_YOUSO[7];

  RCpre = (double***)malloc(sizeof(double**)*List_YOUSO[2]);
  for (i=0; i<List_YOUSO[2]; i++){
    RCpre[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
    for (j=0; j<List_YOUSO[7]; j++){
      RCpre[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
      for (k=0; k<List_YOUSO[7]; k++){
        RCpre[i][j][k] = 0.0;
      }
    }
  }

  RCcnt = (double***)malloc(sizeof(double**)*List_YOUSO[2]);
  for (i=0; i<List_YOUSO[2]; i++){
    RCcnt[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
    for (j=0; j<List_YOUSO[7]; j++){
      RCcnt[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
      for (k=0; k<List_YOUSO[7]; k++){
        RCcnt[i][j][k] = 0.0;
      }
    }
  }

  RRcnt = (double***)malloc(sizeof(double**)*List_YOUSO[2]);
  for (i=0; i<List_YOUSO[2]; i++){
    RRcnt[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
    for (j=0; j<List_YOUSO[7]; j++){
      RRcnt[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
      for (k=0; k<List_YOUSO[7]; k++){
        RRcnt[i][j][k] = 0.0;
      }
    }
  }

  LCpre = (double***)malloc(sizeof(double**)*List_YOUSO[2]);
  for (i=0; i<List_YOUSO[2]; i++){
    LCpre[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
    for (j=0; j<List_YOUSO[7]; j++){
      LCpre[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
      for (k=0; k<List_YOUSO[7]; k++){
        LCpre[i][j][k] = 0.0;
      }
    }
  }

  LCcnt = (double***)malloc(sizeof(double**)*List_YOUSO[2]);
  for (i=0; i<List_YOUSO[2]; i++){
    LCcnt[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
    for (j=0; j<List_YOUSO[7]; j++){
      LCcnt[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
      for (k=0; k<List_YOUSO[7]; k++){
        LCcnt[i][j][k] = 0.0;
      }
    }
  }

  LRcnt = (double***)malloc(sizeof(double**)*List_YOUSO[2]);
  for (i=0; i<List_YOUSO[2]; i++){
    LRcnt[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
    for (j=0; j<List_YOUSO[7]; j++){
      LRcnt[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
      for (k=0; k<List_YOUSO[7]; k++){
        LRcnt[i][j][k] = 0.0;
      }
    }
  }

  RR0 = (double***)malloc(sizeof(double**)*List_YOUSO[2]);
  for (i=0; i<List_YOUSO[2]; i++){
    RR0[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
    for (j=0; j<List_YOUSO[7]; j++){
      RR0[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
      for (k=0; k<List_YOUSO[7]; k++){
        RR0[i][j][k] = 0.0;
      }
    }
  }

  LR0 = (double***)malloc(sizeof(double**)*List_YOUSO[2]);
  for (i=0; i<List_YOUSO[2]; i++){
    LR0[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
    for (j=0; j<List_YOUSO[7]; j++){
      LR0[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
      for (k=0; k<List_YOUSO[7]; k++){
        LR0[i][j][k] = 0.0;
      }
    }
  }

  size_RU = List_YOUSO[3]*List_YOUSO[2]*List_YOUSO[7]*List_YOUSO[7];

  RU = (double****)malloc(sizeof(double***)*List_YOUSO[3]);
  for (i=0; i<List_YOUSO[3]; i++){
    RU[i] = (double***)malloc(sizeof(double**)*List_YOUSO[2]);
    for (j=0; j<List_YOUSO[2]; j++){
      RU[i][j] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
      for (k=0; k<List_YOUSO[7]; k++){
        RU[i][j][k] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
        for (l=0; l<List_YOUSO[7]; l++){
          RU[i][j][k][l] = 0.0;
	}
      }
    }
  }  

  if (firsttime) {
  PrintMemory("TSB_Lanczos: BeGa",sizeof(double)*size_BeGa,NULL);
  PrintMemory("TSB_Lanczos: UmRn",sizeof(double)*size_BeGa,NULL);
  PrintMemory("TSB_Lanczos: RnUm",sizeof(double)*size_BeGa,NULL);

  PrintMemory("TSB_Lanczos: RCpre",sizeof(double)*size_RCpre,NULL);
  PrintMemory("TSB_Lanczos: RCcnt",sizeof(double)*size_RCpre,NULL);
  PrintMemory("TSB_Lanczos: RRcnt",sizeof(double)*size_RCpre,NULL);
  PrintMemory("TSB_Lanczos: LCpre",sizeof(double)*size_RCpre,NULL);
  PrintMemory("TSB_Lanczos: LCcnt",sizeof(double)*size_RCpre,NULL);
  PrintMemory("TSB_Lanczos: LRcnt",sizeof(double)*size_RCpre,NULL);
  PrintMemory("TSB_Lanczos: RR0",  sizeof(double)*size_RCpre,NULL);
  PrintMemory("TSB_Lanczos: LR0",  sizeof(double)*size_RCpre,NULL);
  PrintMemory("TSB_Lanczos: RU",   sizeof(double)*size_RU,   NULL);
  firsttime=0;
  }


  printf("B2\n");

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

	LCpre[i][i1][j1] = 0.0;
	LCcnt[i][i1][j1] = 0.0;
	LRcnt[i][i1][j1] = 0.0;
      }
    }
  }

  printf("B3\n");

  /****************************************************
           Consider RU as column vectors.
           Consider LU as row vectors
  ****************************************************/

  for (i=0; i<=can; i++){
    for (j=0; j<=rlmax; j++){
      for (i1=0; i1<=(ct_on-1); i1++){
	for (j1=0; j1<=(ct_on-1); j1++){
	  LU[spin][ct_AN][j][i][i1][j1] = 0.0;
	}
      } 
    }
  }

  for (i=0; i<=(ct_on-1); i++){
    RCcnt[0][i][i] = 1.0; 
    LCcnt[0][i][i] = 1.0;
    LU[spin][ct_AN][0][0][i][i] = 1.0;
  }

  for (i=0; i<List_YOUSO[7]; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      be[spin][ct_AN][0][i][j] = 0.0;
      ga[spin][ct_AN][0][i][j] = 0.0;
    }
  }

  /****************************************************
     Recursion of two-sided block Lanczos algorithm
  ****************************************************/

  Rloop_END = 0;
  RNUM[ct_AN] = rlmax;
  RCN_Rank[ct_AN][0] = ct_on;

  printf("B4\n");

  rl = 0;
  do {

    /****************************************************
      H'|WRn}

      RCcnt[j][m][n] is a column vector.

      j  ->  local index for atom
      m  ->  index for basis
      n  ->  index for basis
    ****************************************************/

  printf("B5\n");


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

    /****************************************************
      {WLn|H'

      LCcnt[j][m][n] is a row vector.

      j  ->  local index for atom
      m  ->  index for basis
      n  ->  index for basis
    ****************************************************/

    for (i=0; i<=can; i++){
      ig = natn[ct_AN][i];
      ian = Spe_Total_CNO[WhatSpecies[ig]];
      po = 0;
      for (j=0; j<=can; j++){
        kl = RMI2[ct_AN][j][i];
        jg = natn[ct_AN][j];
        jan = Spe_Total_CNO[WhatSpecies[jg]];

        if (0<=kl){
	  for (m=0; m<RCN_Rank[ct_AN][rl]; m++){
	    for (n=0; n<=(ian-1); n++){
	      sum = 0.0;
	      for (k=0; k<=(jan-1); k++){
		sum = sum + LCcnt[j][m][k]*HP[jg][kl][k][n];
	      }
	      if (po==0)
		LRcnt[i][m][n] = sum;
	      else
		LRcnt[i][m][n] += sum;
	    }
	  } 
	  po++;
	}
      }
    }


    /*
  printf("B7\n");
    for (i=0; i<=can; i++){
      printf("i=%3d LRcnt=%15.12f\n",i,LRcnt[i][0][0]); 
    }
    */


    /****************************************************
       Alpha_n = {WLn|H'|WRn}

                  Alpha_n is calculated as
                {WLn|H'|WRn} = <LCcnt|RRcnt>.
    ****************************************************/

    for (i=0; i<=can; i++){
      ig = natn[ct_AN][i];
      ian = Spe_Total_CNO[WhatSpecies[ig]];
      for (m=0; m<RCN_Rank[ct_AN][rl]; m++){
        for (n=0; n<RCN_Rank[ct_AN][rl]; n++){
          sum = 0.0;
          for (k=0; k<=(ian-1); k++){
            sum = sum + LCcnt[i][m][k]*RRcnt[i][k][n];
          }
          if (i==0)
            al[spin][ct_AN][rl][m][n] = sum;
          else
            al[spin][ct_AN][rl][m][n] = al[spin][ct_AN][rl][m][n] + sum;
        }
      }
    }


  printf("B4 al\n");

    for (m=0; m<RCN_Rank[ct_AN][rl]; m++){
      for (n=0; n<RCN_Rank[ct_AN][rl]; n++){
        printf("%9.5f ",al[spin][ct_AN][rl][m][n]);
      }
      printf("\n");
    }




    if (rl!=RNUM[ct_AN]){

      /****************************************************
             |RRcnt} = H'|WRn} - |WRn-1}Bn - |WRn}An
      ****************************************************/

      for (i=0; i<=can; i++){
	ig = natn[ct_AN][i];
	ian = Spe_Total_CNO[WhatSpecies[ig]];    
	for (m=0; m<ian; m++){
	  for (n=0; n<RCN_Rank[ct_AN][rl]; n++){
	    sum = 0.0;


	    if (1<=rl){ 
	      for (k=0; k<RCN_Rank[ct_AN][rl-1]; k++){
		sum += RCpre[i][m][k]*be[spin][ct_AN][rl][k][n];
	      }
	    }
	    for (k=0; k<RCN_Rank[ct_AN][rl]; k++){
	      sum += RCcnt[i][m][k]*al[spin][ct_AN][rl][k][n];
	    }
	    RRcnt[i][m][n] -= sum;
	  }
	}
      }

      /****************************************************
             {LRcnt| = {WLn|H' - Cn{WLn-1| - An{WLn|
      ****************************************************/

    /*
    printf("B8\n");
    for (i=0; i<=can; i++){
      printf("i=%3d LRcnt=%15.12f\n",i,LRcnt[i][0][0]); 
    }
    */

      for (i=0; i<=can; i++){
	ig = natn[ct_AN][i];
	ian = Spe_Total_CNO[WhatSpecies[ig]];    
	for (m=0; m<RCN_Rank[ct_AN][rl]; m++){
	  for (n=0; n<=(ian-1); n++){
	    sum = 0.0;
            if (1<=rl){ 
	      for (k=0; k<RCN_Rank[ct_AN][rl-1]; k++){
	       sum += ga[spin][ct_AN][rl][m][k]*LCpre[i][k][n];
	      }
	    }
	    for (k=0; k<RCN_Rank[ct_AN][rl]; k++){
	      sum += al[spin][ct_AN][rl][m][k]*LCcnt[i][k][n];
	    }
	    LRcnt[i][m][n] -= sum;
	  }
	}
      }

      /*
  printf("B9\n");
    for (i=0; i<=can; i++){
      printf("i=%3d LRcnt=%15.12f\n",i,LRcnt[i][0][0]); 
    }
      */

      /****************************************************
                     Rebiorthogonalization
           by the two-sided block Gram-Schmidt method
      ****************************************************/


      /* |RRcnt) := |RRcnt) - (|RR0) := |RR0) + |RU_rl0)(LU_rl0|RRcnt)) */

      for (rl0=0; rl0<=rl; rl0++){

        /* (LU_rl0|RRcnt) */

        for (i=0; i<=can; i++){
	  ig = natn[ct_AN][i];
	  ian = Spe_Total_CNO[WhatSpecies[ig]];
	  for (m=0; m<RCN_Rank[ct_AN][rl0]; m++){
	    for (n=0; n<RCN_Rank[ct_AN][rl]; n++){
	      sum = 0.0;
	      for (k=0; k<=(ian-1); k++){
		sum = sum + LU[spin][ct_AN][rl0][i][m][k]*RRcnt[i][k][n];
	      }
	      if (i==0)
		UmRn[m][n] = sum;
	      else
		UmRn[m][n] = UmRn[m][n] + sum;
	    }
	  }
	}

        /* |RR0) := |RR0) + |RU_rl0)(LU_rl0|RRn) */

        
	for (i=0; i<=can; i++){
	  ig = natn[ct_AN][i];
	  ian = Spe_Total_CNO[WhatSpecies[ig]];    
	  for (m=0; m<=(ian-1); m++){
	    for (n=0; n<RCN_Rank[ct_AN][rl]; n++){
	      sum = 0.0;
	      for (k=0; k<RCN_Rank[ct_AN][rl0]; k++){
		sum += RU[rl0][i][m][k]*UmRn[k][n];
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

      /* (LRcnt| := (LRcnt| - ((LR0| := (LR0| + (LRcnt|RU_rl0)(LU_rl0|) */







     



      for (rl0=0; rl0<=rl; rl0++){


	/*
    for (i=0; i<=can; i++){
      printf("i=%3d LRcnt=%15.12f RU=%15.12f\n",i,LRcnt[i][0][0],RU[rl0][i][0][0]); 
    }
	*/

        /* (LRcnt|RU_rl0) */

        for (i=0; i<=can; i++){
	  ig = natn[ct_AN][i];
	  ian = Spe_Total_CNO[WhatSpecies[ig]];
	  for (m=0; m<RCN_Rank[ct_AN][rl]; m++){
	    for (n=0; n<RCN_Rank[ct_AN][rl0]; n++){
	      sum = 0.0;
	      for (k=0; k<ian; k++){
                sum = sum + LRcnt[i][m][k]*RU[rl0][i][k][n];
	      }
	      if (i==0)
                RnUm[m][n]  = sum;
              else
	        RnUm[m][n] += sum;
	    } 
	  }
	}

        /* (LR0| := (LR0| + (LRcnt|RU_rl0)(LU_rl0| */


	/*
  printf("B14 RnUm=%15.12f\n",RnUm[0][0]);
    for (i=0; i<=can; i++){
      printf("i=%3d LU=%15.12f\n",i,LU[0][ct_AN][rl0][i][0][0]); 
    }
	*/


	for (i=0; i<=can; i++){
	  ig = natn[ct_AN][i];
	  ian = Spe_Total_CNO[WhatSpecies[ig]];
	  for (m=0; m<RCN_Rank[ct_AN][rl]; m++){
	    for (n=0; n<ian; n++){
	      sum = 0.0;
	      for (k=0; k<RCN_Rank[ct_AN][rl0]; k++){
                sum += RnUm[m][k]*LU[spin][ct_AN][rl0][i][k][n];
              }
	      if (rl0==0) LR0[i][m][n]  = sum;
              else        LR0[i][m][n] += sum;
	    }
	  }
	}
      }

      /* (LRcnt| := (LRcnt| - (LR0| */

      /*
  printf("B15\n");
    for (i=0; i<=can; i++){
      printf("i=%3d LRcnt=%15.12f\n",i,LRcnt[i][0][0]); 
    }
      */

      for (i=0; i<=can; i++){
	ig = natn[ct_AN][i];
	ian = Spe_Total_CNO[WhatSpecies[ig]]; 
	for (m=0; m<RCN_Rank[ct_AN][rl]; m++){
	  for (n=0; n<=(ian-1); n++){
            LRcnt[i][m][n] = LRcnt[i][m][n] - LR0[i][m][n];           
	  } 
	}
      }

      /*
  printf("B16\n");
    for (i=0; i<=can; i++){
      printf("i=%3d LRcnt=%15.12f\n",i,LRcnt[i][0][0]); 
    }
      */

      /****************************************************
                    Bn+1 * Cn+1 = {LRcnt|RRcnt}
      ****************************************************/


      for (i=0; i<=can; i++){
	ig = natn[ct_AN][i];
	ian = Spe_Total_CNO[WhatSpecies[ig]];    
	for (m=0; m<RCN_Rank[ct_AN][rl]; m++){
	  for (n=0; n<RCN_Rank[ct_AN][rl]; n++){
	    sum = 0.0;
	    for (k=0; k<=(ian-1); k++){
	      sum = sum + LRcnt[i][m][k]*RRcnt[i][k][n];
	    }
	    if (i==0)
	      BeGa[m][n] = sum;
	    else
	      BeGa[m][n] = BeGa[m][n] + sum;
	  } 
	}
      }

      /****************************************************
	1) Perform the singular value decomposition of 
           {LRcnt|RRcnt} as U*W*V^t.
        2) Then, B_{n+1} = UW^{1/2} and C_{n+1}=W^{1/2}V^t.
        3) It is easy to calculate the inverses of B and C.
      ****************************************************/



      printf("BeGa  ct_AN=%2d  rl=%2d\n",ct_AN,rl);
      for (m=0; m<RCN_Rank[ct_AN][rl]; m++){
        for (n=0; n<RCN_Rank[ct_AN][rl]; n++){
          printf("%15.12f ",BeGa[m][n]);
	}
        printf("\n");
      }


      /*
      MPI_Finalize();
      exit(0);
      */



      /* Singular value decomposition */


      RCN_Rank[ct_AN][rl+1] = SVD_fact_inv(
                                rl,RCN_Rank[ct_AN][rl]-1,BeGa,
                                be[spin][ct_AN][rl+1],ibe[spin][ct_AN][rl+1],
                                ga[spin][ct_AN][rl+1],iga[spin][ct_AN][rl+1]);


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
      printf("RCN_Rank[ct_AN][rl+1]=%i  Rloop_END=%i  ct_on=%i\n",
              RCN_Rank[ct_AN][rl+1],Rloop_END,ct_on);
      */

      /****************************************************
                     |WRn+1} = |RRcnt}ICn+1
      ****************************************************/


      for (k=0; k<RCN_Rank[ct_AN][rl]; k++){
        for (i=0; i<=can; i++){
  	  ig = natn[ct_AN][i];
	  ian = Spe_Total_CNO[WhatSpecies[ig]];    
  	  for (m=0; m<ian; m++){
            printf("1 RRcnt k=%2d %15.12f\n",k,RRcnt[i][m][k]);
	  }
	}
      }



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
              sum = sum + RRcnt[i][m][k]*iga[spin][ct_AN][rl+1][k][n];
	    }
            RCcnt[i][m][n] = sum;
            RU[rl+1][i][m][n] = sum;
	  }
	}
      }


      for (k=0; k<RCN_Rank[ct_AN][rl]; k++){
        for (i=0; i<=can; i++){
  	  ig = natn[ct_AN][i];
	  ian = Spe_Total_CNO[WhatSpecies[ig]];    
  	  for (m=0; m<ian; m++){
            printf("2 RRcnt k=%2d %15.12f\n",k,RRcnt[i][m][k]);
	  }
	}
      }


      /*
      for (i=0; i<=can; i++){
        printf("i=%3d RCcnt=%15.12f\n",i,RCcnt[i][0][0]); 
      }

  printf("B20\n");
      */

      /****************************************************
                      {WLn+1| = IBn+1{LRcnt|
      ****************************************************/

      /*
      for (i=0; i<=can; i++){
        printf("i=%3d LRcnt=%15.12f\n",i,LRcnt[i][0][0]); 
      }
      */


      /*
      for (k=0; k<RCN_Rank[ct_AN][rl]; k++){
        for (i=0; i<=can; i++){
  	  ig = natn[ct_AN][i];
	  ian = Spe_Total_CNO[WhatSpecies[ig]];    
  	  for (m=0; m<ian; m++){
            printf("1 LRcnt k=%2d %15.12f\n",k,LRcnt[i][k][m]);
	  }
	}
      }
      */


      for (i=0; i<=can; i++){
	ig = natn[ct_AN][i];
	ian = Spe_Total_CNO[WhatSpecies[ig]];

	for (m=0; m<RCN_Rank[ct_AN][rl]; m++){
	  for (n=0; n<ian; n++){
	    LCpre[i][m][n] = LCcnt[i][m][n];
	  }
	}

	for (m=0; m<RCN_Rank[ct_AN][rl+1]; m++){
	  for (n=0; n<ian; n++){
	    sum = 0.0;
	    for (k=0; k<RCN_Rank[ct_AN][rl]; k++){
	      sum = sum + ibe[spin][ct_AN][rl+1][m][k]*LRcnt[i][k][n];
	    }
	    LCcnt[i][m][n] = sum;
            LU[spin][ct_AN][rl+1][i][m][n] = sum;
	  }
	}
      }


      /*
      for (k=0; k<RCN_Rank[ct_AN][rl]; k++){
        for (i=0; i<=can; i++){
  	  ig = natn[ct_AN][i];
	  ian = Spe_Total_CNO[WhatSpecies[ig]];    
  	  for (m=0; m<ian; m++){
            printf("2 LRcnt k=%2d %15.12f\n",k,LRcnt[i][k][m]);
	  }
	}
      }
      */


      /*
      for (i=0; i<=can; i++){
        printf("i=%3d LCcnt=%15.12f\n",i,LCcnt[i][0][0]); 
      }
      */

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


  /****************************************************
    freeing of arrays:
  ****************************************************/

  for (i=0; i<List_YOUSO[7]; i++){
    free(BeGa[i]);
  }
  free(BeGa);

  for (i=0; i<List_YOUSO[7]; i++){
    free(UmRn[i]);
  }
  free(UmRn);

  for (i=0; i<List_YOUSO[7]; i++){
    free(RnUm[i]);
  }
  free(RnUm);

  for (i=0; i<List_YOUSO[2]; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(RCpre[i][j]);
    }
    free(RCpre[i]);
  }
  free(RCpre);

  for (i=0; i<List_YOUSO[2]; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(RCcnt[i][j]);
    }
    free(RCcnt[i]);
  }
  free(RCcnt);

  for (i=0; i<List_YOUSO[2]; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(RRcnt[i][j]);
    }
    free(RRcnt[i]);
  }
  free(RRcnt);

  for (i=0; i<List_YOUSO[2]; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(LCpre[i][j]);
    }
    free(LCpre[i]);
  }
  free(LCpre);

  for (i=0; i<List_YOUSO[2]; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(LCcnt[i][j]);
    }
    free(LCcnt[i]);
  }
  free(LCcnt);

  for (i=0; i<List_YOUSO[2]; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(LRcnt[i][j]);
    }
    free(LRcnt[i]);
  }
  free(LRcnt);

  for (i=0; i<List_YOUSO[2]; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(RR0[i][j]);
    }
    free(RR0[i]);
  }
  free(RR0);

  for (i=0; i<List_YOUSO[2]; i++){
    for (j=0; j<List_YOUSO[7]; j++){
      free(LR0[i][j]);
    }
    free(LR0[i]);
  }
  free(LR0);

  for (i=0; i<List_YOUSO[3]; i++){
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
  static double scale;
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
    if ( (fabs(w[rankN]/w[0])<rank_criterion
	  && fabs(w[rankN])<100.0*rank_criterion )
        ||
         (fabs(w[rankN])<rank_criterion2 && rankN==0) ) po = 1; 
  } while(rankN<n && po==0);
  if (po==0) rankN++;



  printf("rankN=%2d\n",rankN);
  for (i=0; i<=n; i++){
    printf("i=%2d w=%15.12f\n",i,w[i]);
  }


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

  scale = 2.0;

  /*
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
  */



  for (i=0; i<=n; i++){
    for (j=0; j<=n; j++){
      B[i][j] = u1[i][j]*w[j]/scale;
    }
  }

  for (i=0; i<=n; i++){
    for (j=0; j<=n; j++){
      IB[i][j] = u1[j][i]/w[i]*scale;
    }
  }



  /****************************************************
            calculate gamma and inverse gamma
  ****************************************************/

  /*
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
  */


  for (i=0; i<=n; i++){
    for (j=0; j<=n; j++){
      C[i][j] = vt1[i][j]*scale;
    }
  }

  for (i=0; i<=n; i++){
    for (j=0; j<=n; j++){
      IC[i][j] = vt1[j][i]/scale;
    }
  }




  printf("IB \n");

  for (i=0; i<=n; i++){
    for (j=0; j<=n; j++){
      printf("%10.5f ",IB[i][j]);
    }
    printf("\n");
  }

  printf("IC \n");

  for (i=0; i<=n; i++){
    for (j=0; j<=n; j++){
      printf("%10.5f ",IC[i][j]);
    }
    printf("\n");
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










void T_Conserve_MP(int T_switch)
{

  /****************************************************
    T_switch is a parameter to select a terminator
    in the Green's function.

    T_switch = 1 is no terminator.
    T_switch = 2 is a square root terminator.
  ****************************************************/

  static int firsttime=1;
  static int spin,i,j,k,n,p,po,ct_AN,rl_loop,doko;
  static int wan,tno,tno1,itnum,Num_loop;
  static double TN,TX,TZ,TQ,pTQ,Delta,pDelta;
  static double E_Temp0,rc,ac,DChemP,dtnum;
  static double dp,dM,dum,dum1;
  static dcomplex EpC,EpP,CN,CX;
  static dcomplex ***G00;
  static dcomplex *CNd,*CXd;
  static dcomplex cdum,Csum3;

  /****************************************************
    allocation of arrays:
  ****************************************************/

  G00 = (dcomplex***)malloc(sizeof(dcomplex**)*(SpinP_switch+1)); 
  for (k=0; k<=SpinP_switch; k++){
    G00[k] = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
    for (i=0; i<List_YOUSO[7]; i++){
      G00[k][i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
    }
  }  

  CNd = (dcomplex*)malloc(sizeof(dcomplex)*(atomnum+1)); 
  CXd = (dcomplex*)malloc(sizeof(dcomplex)*(atomnum+1)); 

  /****************************************************
    start calculation
  ****************************************************/

  TZ = 0.0;
  for (i=1; i<=atomnum; i++){
    wan = WhatSpecies[i];
    TZ = TZ + Spe_Core_Charge[wan];
  }

  printf("<Recursion>  Ideal num.of.electrons = %15.10f\n",TZ);

  po = 150;
  Num_loop = 0;
  TQ = 10.0;
  ac = 1.0;
  rc = 1.0;
  DChemP = 0.0040*eV2Hartree;
  pTQ = 10000;
  itnum = Av_num;
  dtnum = Av_num;

  while(0<po){

    Num_loop++;

    /****************************************************
                     Total_N and Total_X
    ****************************************************/

    CN = Complex(0.0,0.0);
    CX = Complex(0.0,0.0);

    for (p=0; p<=POLES-1; p++){
      EpC = RCadd(ChemP,Ep[p]);
      EpP.r = EpC.r;
      EpP.i = EpC.i;

      for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
	wan = WhatSpecies[ct_AN];
	tno = Spe_Total_CNO[wan];
	tno1 = tno - 1;

        /****************************************************
                          Multiple inverse
        ****************************************************/

        for (spin=0; spin<=SpinP_switch; spin++){
          GBTD00(T_switch,ct_AN,spin,EpP,G00[spin]);
        }

        /****************************************************
                         Mulliken population
        ****************************************************/

	Csum3.r = 0.0;
	Csum3.i = 0.0;
        for (spin=0; spin<=SpinP_switch; spin++){
          for (j=0; j<=tno1; j++){
 	    Csum3.r = Csum3.r + G00[spin][j][j].r;
	    Csum3.i = Csum3.i + G00[spin][j][j].i;
	  }
        }
	CNd[ct_AN].r = Csum3.r;
	CNd[ct_AN].i = Csum3.i;
        
        /****************************************************
                          Response functions
        ****************************************************/
        
	Csum3.r = 0.0;
	Csum3.i = 0.0;
        for (spin=0; spin<=SpinP_switch; spin++){
          for (j=0; j<=tno1; j++){
	    cdum.r = G00[spin][j][j].r*G00[spin][j][j].r
                   - G00[spin][j][j].i*G00[spin][j][j].i;
	    cdum.i = 2.0*G00[spin][j][j].i*G00[spin][j][j].r;
	    Csum3.r = Csum3.r + cdum.r;
	    Csum3.i = Csum3.i + cdum.i;
          }
        }
	CXd[ct_AN].r = Csum3.r;
	CXd[ct_AN].i = Csum3.i;

        /****************************************************
                         End of ct_AN loop 
        ****************************************************/

      }

      /****************************************************
               Summation of Mulliken populations and
               responce functions over atoms
      ****************************************************/

      Csum3.r = 0.0;
      Csum3.i = 0.0;
      for (i=1; i<=atomnum; i++){
	cdum.r  = zp[p].r*CNd[i].r - zp[p].i*CNd[i].i;
	cdum.i  = zp[p].i*CNd[i].r + zp[p].r*CNd[i].i;
	Csum3.r = Csum3.r + cdum.r;
	Csum3.i = Csum3.i + cdum.i;
      }
      CN.r = CN.r + Csum3.r;
      CN.i = CN.i + Csum3.i;
      Csum3.r = 0.0;
      Csum3.i = 0.0;
      for (i=1; i<=atomnum; i++){
	cdum.r  = zp[p].r*CXd[i].r - zp[p].i*CXd[i].i;
	cdum.i  = zp[p].i*CXd[i].r + zp[p].r*CXd[i].i;
	Csum3.r = Csum3.r + cdum.r;
	Csum3.i = Csum3.i + cdum.i;
      }
      CX.r = CX.r + Csum3.r;
      CX.i = CX.i + Csum3.i;

      /****************************************************
                       End of POLE loop
      ****************************************************/

    }

    if (SpinP_switch==0){
      TN = 4.0*CN.r/Beta;
      TX = 4.0*CX.r/Beta;
    }
    else if (SpinP_switch==1){
      TN = 2.0*CN.r/Beta;
      TX = 2.0*CX.r/Beta;
    }

    /****************************************************
                Correction of Chemical Potential
    ****************************************************/


    TQ = TZ - TN - Given_Total_Charge;

    if (fabs(TQ)<(CN_Error*(double)atomnum)) po=0;
    else                                     po--;

    if (0<po){
      if (0.1<fabs(TQ/TZ) || 1.0<fabs(TQ)){
	if (0.0<TQ){
	  ChemP = ChemP + DChemP;
          doko = 1;
	}
	else{
	  ChemP = ChemP - DChemP;
          doko = 2;
	}
	if (fabs(pTQ)<fabs(TQ)) DChemP = 0.5*DChemP;
      }
      else{
	if (fabs(pTQ)<fabs(TQ) && pTQ*TQ<0.0){
	  rc = rc + 3.0;
	  Delta = -0.5*pDelta;
	  pDelta = Delta;
          doko = 3;
	}
	else{
	  Delta = -ac*(TZ-TN-Given_Total_Charge)/TX/rc;
	  pDelta = Delta;
          doko = 4;          

	  if ((0.5<fabs(TQ/pTQ)) && (fabs(TQ/TZ)<0.05)){
	    if (0.0<(TQ/pTQ)){
	      ac = ac + 0.5;
              doko = 5;          
	    }
	    else{
	      ac = ac/3.0;
              doko = 6;
	    }
	  }
	}
	ChemP = ChemP + Delta;
      }
      pTQ = TQ;
    }

    printf("<Recursion> %3d doko=%2d  Num.of.eles= %15.10f   ChemP= %15.10f\n",
            Num_loop,doko,TN,eV2Hartree*ChemP);
  }

  /****************************************************
    freeing of arrays:
  ****************************************************/

  for (k=0; k<=SpinP_switch; k++){
    for (i=0; i<List_YOUSO[7]; i++){
      free(G00[k][i]);
    }
    free(G00[k]);
  }  
  free(G00);

  free(CNd);
  free(CXd);

}













void BO_Calc_S(int T_switch, double *****CDM, double *****EDM, double ****IOLP)
{
  /****************************************************
    T_switch is a parameter to select a terminator
    in the Green's function.
    T_switch = 1 is no terminator.
    T_switch = 2 is a square root terminator.
  ****************************************************/
  static int firsttime=1;
  static int ct_AN,spin,rl_loop;
  static int wan,tno,tno1,itnum,ig,ino1,jg,jno1,kg,kj,kno1;
  static int i,j,k,p,l,m,n,q,rl,i1,j1,h_AN;
  static double sum,sum2,sumx,sumy,sumz,dtnum;

  static dcomplex EpC,EpP,CE,CRho,CRes,CERho;
  static dcomplex ***G00;
  static dcomplex ****G0;

  static double ****Rho;
  static double ****ERho;
  static double **E00;
  static double *Fx,*Fy,*Fz;
  static double **Uele_temp;
  static double **Dmat,**Dmat2;
  static double dum,dum1,dum2;

  /****************************************************
    allocation of arrays:
  ****************************************************/

  G00 = (dcomplex***)malloc(sizeof(dcomplex**)*(SpinP_switch+1)); 
  for (k=0; k<=SpinP_switch; k++){
    G00[k] = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
    for (i=0; i<List_YOUSO[7]; i++){
      G00[k][i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
    }
  }  

  G0 = (dcomplex****)malloc(sizeof(dcomplex***)*(SpinP_switch+1)); 
  for (k=0; k<=SpinP_switch; k++){
    G0[k] = (dcomplex***)malloc(sizeof(dcomplex**)*List_YOUSO[3]);
    for (i=0; i<List_YOUSO[3]; i++){
      G0[k][i] = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
      for (j=0; j<List_YOUSO[7]; j++){
        G0[k][i][j] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
      }
    }
  }  

  Rho = (double****)malloc(sizeof(double***)*(SpinP_switch+1)); 
  for (k=0; k<=SpinP_switch; k++){
    Rho[k] = (double***)malloc(sizeof(double**)*List_YOUSO[3]);
    for (i=0; i<List_YOUSO[3]; i++){
      Rho[k][i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
      for (j=0; j<List_YOUSO[7]; j++){
        Rho[k][i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
      }
    }
  }  

  ERho = (double****)malloc(sizeof(double***)*(SpinP_switch+1)); 
  for (k=0; k<=SpinP_switch; k++){
    ERho[k] = (double***)malloc(sizeof(double**)*List_YOUSO[3]);
    for (i=0; i<List_YOUSO[3]; i++){
      ERho[k][i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
      for (j=0; j<List_YOUSO[7]; j++){
        ERho[k][i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
      }
    }
  }  

  E00 = (double**)malloc(sizeof(double*)*(SpinP_switch+1)); 
  for (k=0; k<=SpinP_switch; k++){
    E00[k] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
  }

  Fx = (double*)malloc(sizeof(double)*(atomnum+1));
  Fy = (double*)malloc(sizeof(double)*(atomnum+1));
  Fz = (double*)malloc(sizeof(double)*(atomnum+1));

  Uele_temp = (double**)malloc(sizeof(double*)*(SpinP_switch+1)); 
  for (k=0; k<=SpinP_switch; k++){
    Uele_temp[k] = (double*)malloc(sizeof(double)*(atomnum+1)); 
  }

  Dmat = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
  for (j=0; j<List_YOUSO[7]; j++){
    Dmat[j] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
  }

  Dmat2 = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
  for (j=0; j<List_YOUSO[7]; j++){
    Dmat2[j] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
  }




  if (firsttime) {
  PrintMemory("BO_Calc_S: Rho",sizeof(Rho),NULL);
  PrintMemory("BO_Calc_S: ERho",sizeof(ERho),NULL);
  firsttime=0;
  }


  itnum = Av_num;
  dtnum = Av_num;

  for (spin=0; spin<=SpinP_switch; spin++){
    Uele_IS[spin] = 0.0;
  }

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    wan = WhatSpecies[ct_AN];
    tno = Spe_Total_CNO[wan];
    tno1 = tno - 1;

    for (spin=0; spin<=SpinP_switch; spin++){
      for (i=0; i<=tno1; i++) E00[spin][i] = 0.0;
    }

    for (p=0; p<POLES; p++){
      EpC.r = ChemP + Ep[p].r;
      EpC.i = Ep[p].i;
      EpP.r = EpC.r;
      EpP.i = EpC.i;

      /****************************************************
                        Multiple inverse
      ****************************************************/

      for (spin=0; spin<=SpinP_switch; spin++){
        GBTD00(T_switch,ct_AN,spin,EpP,G00[spin]);
      }

      /****************************************************
                      Recurrence relation
      ****************************************************/

      for (spin=0; spin<=SpinP_switch; spin++){
        RecurG(ct_AN,spin,EpP,G00[spin],G0[spin]);
      } 

      for (spin=0; spin<=SpinP_switch; spin++){
        for (rl=0; rl<=RNUM[ct_AN]; rl++){
	  for (i=0; i<=tno1; i++){
	    for (j=0; j<RCN_Rank[ct_AN][rl]; j++){

	      CRho.r = zp[p].r*G0[spin][rl][i][j].r
                     - zp[p].i*G0[spin][rl][i][j].i;
	      CRho.i = zp[p].i*G0[spin][rl][i][j].r
                     + zp[p].r*G0[spin][rl][i][j].i;

	      if (p==0)
	        Rho[spin][rl][i][j] = CRho.r;
	      else
	        Rho[spin][rl][i][j] = Rho[spin][rl][i][j] + CRho.r;

              CERho.r = CRho.r*EpC.r - CRho.i*EpC.i;
	      if (p==0)
	        ERho[spin][rl][i][j] = CERho.r;
	      else
	        ERho[spin][rl][i][j] = ERho[spin][rl][i][j] + CERho.r;

              if (rl==0 && i==j) E00[spin][i] = E00[spin][i] + CERho.r;

            }
	  }
        }
      }

      /****************************************************
                         End of p loop
      ****************************************************/
    } 

    for (spin=0; spin<=SpinP_switch; spin++){
      for (rl=0; rl<=RNUM[ct_AN]; rl++){
        for (i=0; i<=tno1; i++){
	  for (j=0; j<RCN_Rank[ct_AN][rl]; j++){
            Rho[spin][rl][i][j]  = 2.0*Rho[spin][rl][i][j]/Beta;
            ERho[spin][rl][i][j] = 2.0*ERho[spin][rl][i][j]/Beta;
	  }
        }
      }
    }

    for (spin=0; spin<=SpinP_switch; spin++){
      for (i=0; i<=tno1; i++){
        E00[spin][i] = 2.0*E00[spin][i]/Beta;
      }
    }

    for (spin=0; spin<=SpinP_switch; spin++){
      sum = 0.0;
      for (i=0; i<=tno1; i++) sum += E00[spin][i];
      Uele_temp[spin][ct_AN] = sum;
    }

    /****************************************************
                       Intersite Uele
    ****************************************************/

    for (spin=0; spin<=SpinP_switch; spin++){
      for (i=0; i<=tno1; i++){
        for (k=0; k<=tno1; k++){
          Uele_IS[spin] = Uele_IS[spin]
                        + Rho[spin][0][i][k]*al[spin][ct_AN][0][k][i];
        }
      }
      for (i=0; i<=tno1; i++){
        for (k=0; k<=tno1; k++){
          Uele_IS[spin] = Uele_IS[spin]
                        + Rho[spin][1][i][k]*ga[spin][ct_AN][1][k][i];
        }
      }
    }
 




    for (rl=0; rl<=RNUM[ct_AN]; rl++){

      /*
      printf("Rho  ct_AN=%2d rl=%2d\n",ct_AN,rl);
      */

      sum = 0.0;   
      for (i=0; i<=tno1; i++){
        for (k=0; k<=tno1; k++){
          sum += fabs(Rho[0][rl][i][k]);
	  /*
          printf("%15.12f ",Rho[0][rl][i][k]);
	  */

	}
	/*
        printf("\n");
	*/
      }

      printf("ct_AN=%2d rl=%2d  Sum of Rho = %18.15f\n",ct_AN,rl,sum); 
    }





    /****************************************************
          Transformation of DM and EDM from L basis
                  into the original basis
    ****************************************************/

  printf("R1\n");

    for (spin=0; spin<=SpinP_switch; spin++){
      for (j=0; j<=FNAN[ct_AN]; j++){
        jg = natn[ct_AN][j];
        jno1 = Spe_Total_CNO[WhatSpecies[jg]] - 1;
        for (rl=0; rl<=RNUM[ct_AN]; rl++){
          for (k=0; k<=(FNAN[ct_AN]+SNAN[ct_AN]); k++){
            kg = natn[ct_AN][k];
            kj = RMI2[ct_AN][k][j];
            kno1 = Spe_Total_CNO[WhatSpecies[kg]] - 1;

	    if (0<=kj){

	      for (p=0; p<=tno1; p++){
		for (q=0; q<=kno1; q++){
		  sum = 0.0;
		  sum2 = 0.0;
		  for (m=0; m<RCN_Rank[ct_AN][rl]; m++){
		    sum  += Rho[spin][rl][p][m]*LU[spin][ct_AN][rl][k][m][q]; 
		    sum2 += ERho[spin][rl][p][m]*LU[spin][ct_AN][rl][k][m][q]; 
		  }
		  Dmat[p][q] = sum;
		  Dmat2[p][q] = sum2;
		}
	      }

	      for (p=0; p<=tno1; p++){
		for (q=0; q<=jno1; q++){                
		  sum = 0.0;
		  sum2 = 0.0;
		  for (m=0; m<=kno1; m++){
		    sum  += Dmat[p][m]*IOLP[kg][kj][m][q];
		    sum2 += Dmat2[p][m]*IOLP[kg][kj][m][q];
		  }
		  if (rl==0 && k==0){
		    CDM[spin][ct_AN][j][p][q] = sum;    
		    EDM[spin][ct_AN][j][p][q] = sum2;
		  }
		  else{
		    CDM[spin][ct_AN][j][p][q] += sum;
		    EDM[spin][ct_AN][j][p][q] += sum2;
		  }
		}
	      }

	    }
          }
        }
      }
    }

  printf("R2\n");

    /****************************************************
                      End of ct_AN loop
    ****************************************************/

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
	    CDM[spin][ig][k][q][p] = dum;
	    dum1 = EDM[spin][ct_AN][i][p][q];
	    dum2 = EDM[spin][ig][k][q][p];
	    dum = 0.50*(dum1 + dum2);
	    EDM[spin][ct_AN][i][p][q] = dum;
	    EDM[spin][ig][k][q][p] = dum;
	  }
        }
      } 
    }
  }


  printf("R2\n");


  /****************************************************
                Summation of Uele_temp
  ****************************************************/

  for (spin=0; spin<=SpinP_switch; spin++){
    sum = 0.0;
    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      sum = sum + Uele_temp[spin][ct_AN];
    }
    Uele_OS[spin] = sum;
  }

  if (SpinP_switch==0){
    printf("Uele_OS[0]=%15.12f   Uele_OS[1]=%15.12f\n",Uele_OS[0],Uele_OS[0]);
    printf("Uele_IS[0]=%15.12f   Uele_IS[1]=%15.12f\n",Uele_IS[0],Uele_IS[0]);
  }
  else if (SpinP_switch==1){
    printf("Uele_OS[0]=%15.12f   Uele_OS[1]=%15.12f\n",Uele_OS[0],Uele_OS[1]);
    printf("Uele_IS[0]=%15.12f   Uele_IS[1]=%15.12f\n",Uele_IS[0],Uele_IS[1]);
  }

  printf("R3\n");

  /****************************************************
    freeing of arrays:
  ****************************************************/

  for (k=0; k<=SpinP_switch; k++){
    for (i=0; i<List_YOUSO[7]; i++){
      free(G00[k][i]);
    }
    free(G00[k]);
  }  
  free(G00);

  for (k=0; k<=SpinP_switch; k++){
    for (i=0; i<List_YOUSO[3]; i++){
      for (j=0; j<List_YOUSO[7]; j++){
        free(G0[k][i][j]);
      }
      free(G0[k][i]);
    }
    free(G0[k]);
  }  
  free(G0);

  for (k=0; k<=SpinP_switch; k++){
    for (i=0; i<List_YOUSO[3]; i++){
      for (j=0; j<List_YOUSO[7]; j++){
        free(Rho[k][i][j]);
      }
      free(Rho[k][i]);
    }
    free(Rho[k]);
  }  
  free(Rho);

  for (k=0; k<=SpinP_switch; k++){
    for (i=0; i<List_YOUSO[3]; i++){
      for (j=0; j<List_YOUSO[7]; j++){
        free(ERho[k][i][j]);
      }
      free(ERho[k][i]);
    }
    free(ERho[k]);
  }  
  free(ERho);

  for (k=0; k<=SpinP_switch; k++){
    free(E00[k]);
  }
  free(E00);

  free(Fx);
  free(Fy);
  free(Fz);

  for (k=0; k<=SpinP_switch; k++){
    free(Uele_temp[k]);
  }
  free(Uele_temp);

  for (j=0; j<List_YOUSO[7]; j++){
    free(Dmat[j]);
  }
  free(Dmat);

  for (j=0; j<List_YOUSO[7]; j++){
    free(Dmat2[j]);
  }
  free(Dmat2);

  printf("R4\n");

}











void RecurG(int ct_AN, int spin, dcomplex EpP, dcomplex **G00S, dcomplex ***G0S)
{

  static int i,j,k,rl,wan;
  static int tno,tno1;
  static int rlm,rlm2;
  static double dum,dum1,dum2;
  static dcomplex Csum,Csum2,**Ctp,**Ctp1;
  static dcomplex **Ctp2,**Ctp3;
  static dcomplex Cdum,Cdum2,Cdum3;

  /****************************************************
    allocation of arrays:
  ****************************************************/

  Ctp = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  Ctp1 = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp1[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  Ctp2 = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp2[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  Ctp3 = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    Ctp3[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }


  wan = WhatSpecies[ct_AN];
  tno = Spe_Total_CNO[wan];
  tno1 = tno - 1;

  for (i=0; i<=tno1; i++){
    for (j=0; j<=tno1; j++){
      G0S[0][i][j] = G00S[i][j];
    }
  }
  for (rl=1; rl<=RNUM[ct_AN]; rl++){
    rlm = rl - 1;
    rlm2 = rl - 2;

    /****************************************************
                        G0(n-2)*B(n-1)
    ****************************************************/

    if (rl==1){
      for (i=0; i<=tno1; i++){
	for (j=0; j<=tno1; j++){
	  Ctp[i][j].r = 0.0;
	  Ctp[i][j].i = 0.0;
	}
      }
    }
    else{
      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[ct_AN][rlm]; j++){
	  Csum.r = 0.0;
	  Csum.i = 0.0;
	  for (k=0; k<RCN_Rank[ct_AN][rlm2]; k++){
	    Cdum.r = G0S[rlm2][i][k].r*be[spin][ct_AN][rlm][k][j];
	    Cdum.i = G0S[rlm2][i][k].i*be[spin][ct_AN][rlm][k][j];
	    Csum.r = Csum.r + Cdum.r;
	    Csum.i = Csum.i + Cdum.i;
	  }
	  Ctp[i][j].r = Csum.r;
	  Ctp[i][j].i = Csum.i;
	}
      }
    }

    /****************************************************
                       G0(n-1)*(Z-A(n-1))
    ****************************************************/

    for (i=0; i<RCN_Rank[ct_AN][rlm]; i++){
      for (j=0; j<RCN_Rank[ct_AN][rlm]; j++){
        dum = -al[spin][ct_AN][rlm][i][j];
        if (i==j){
          Ctp1[i][i].r = dum + EpP.r;
          Ctp1[i][i].i = EpP.i;
        } 
        else{
          Ctp1[i][j].r = dum;
          Ctp1[i][j].i = 0.0;
        } 
      }
    }

    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[ct_AN][rlm]; j++){
	Csum.r = 0.0;
	Csum.i = 0.0;
	for (k=0; k<RCN_Rank[ct_AN][rlm]; k++){
	  Cdum.r = G0S[rlm][i][k].r*Ctp1[k][j].r
                 - G0S[rlm][i][k].i*Ctp1[k][j].i;
          Cdum.i = G0S[rlm][i][k].i*Ctp1[k][j].r
                 + G0S[rlm][i][k].r*Ctp1[k][j].i;
          Csum.r = Csum.r + Cdum.r;
	  Csum.i = Csum.i + Cdum.i;
        }
        Ctp2[i][j].r = Csum.r;
        Ctp2[i][j].i = Csum.i;
      }
    }

    if (rl==1){
      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[ct_AN][rlm]; j++){
	  if (i==j){
	    Ctp3[i][j].r = Ctp2[i][j].r - Ctp[i][j].r - 1.0;
	    Ctp3[i][j].i = Ctp2[i][j].i - Ctp[i][j].i;
	  }
	  else{
	    Ctp3[i][j].r = Ctp2[i][j].r - Ctp[i][j].r;
	    Ctp3[i][j].i = Ctp2[i][j].i - Ctp[i][j].i;
	  }
	}
      }
    }
    else{
      for (i=0; i<=tno1; i++){
	for (j=0; j<RCN_Rank[ct_AN][rlm]; j++){
	  Ctp3[i][j].r = Ctp2[i][j].r - Ctp[i][j].r;
	  Ctp3[i][j].i = Ctp2[i][j].i - Ctp[i][j].i;
	}
      }
    }

    /****************************************************
                        [ ]*(Cn)^{-1}
    ****************************************************/

    for (i=0; i<=tno1; i++){
      for (j=0; j<RCN_Rank[ct_AN][rl]; j++){
	Csum.r = 0.0;
	Csum.i = 0.0;
	for (k=0; k<RCN_Rank[ct_AN][rlm]; k++){
	  Cdum.r = Ctp3[i][k].r*iga[spin][ct_AN][rl][k][j];
	  Cdum.i = Ctp3[i][k].i*iga[spin][ct_AN][rl][k][j];
	  Csum.r = Csum.r + Cdum.r;
	  Csum.i = Csum.i + Cdum.i;
	}
	G0S[rl][i][j].r = Csum.r;
	G0S[rl][i][j].i = Csum.i;
      }
    }
  }

  /****************************************************
    freeing of arrays:
  ****************************************************/

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp[i]);
  }
  free(Ctp);

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp1[i]);
  }
  free(Ctp1);

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp2[i]);
  }
  free(Ctp2);

  for (i=0; i<List_YOUSO[7]; i++){
    free(Ctp3[i]);
  }
  free(Ctp3);

}





void GBTD00(int T_switch, int ct_AN, int spin, dcomplex EpP, dcomplex **G00S)
{
  static int i,j,k,l,m,n,rl,wan,Avnum;
  static int tno1,tno,itnum,rl_loop;
  static double dtnum,Av_al,d1,d2;
  static double dum,sum,xd,yd;
  static double c1,s1,c2,s2;
  static double tr,ti,ai,bi;
  static double *be2;
  static dcomplex *cter,**cinv;
  static dcomplex **ctp1,**ctp2;
  static dcomplex Csum3,cdum;

  /****************************************************
    allocation of arrays:
  ****************************************************/

  be2 = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
  cter = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 

  cinv = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    cinv[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  ctp1 = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    ctp1[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  ctp2 = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    ctp2[i] = (dcomplex*)malloc(sizeof(dcomplex)*List_YOUSO[7]); 
  }

  wan = WhatSpecies[ct_AN];
  tno = Spe_Total_CNO[wan];
  tno1 = tno - 1;
  itnum = Av_num;
  dtnum = Av_num;

  /****************************************************
                       Non-Terminator 
  ****************************************************/

  if (T_switch==1){
    for (i=0; i<=tno1; i++){
      cter[i].r = 0.0;
      cter[i].i = 0.0;
    }
  }

  /****************************************************
                  Square Root Terminator
  ****************************************************/

  else{

    /****************************************************
                        Sum of Bn*Cn
    ****************************************************/

    Avnum = 0;    
    for (i=0; i<=tno1; i++) be2[i] = 0.0;
    for (rl=RNUM[ct_AN]; (RNUM[ct_AN]-itnum+1)<=rl; rl--){
      Avnum = Avnum + RCN_Rank[ct_AN][rl-1];
      for (i=0; i<RCN_Rank[ct_AN][rl-1]; i++){
	sum = 0.0;
	for (k=0; k<RCN_Rank[ct_AN][rl]; k++){
	  sum = sum + be[spin][ct_AN][rl][i][k]*ga[spin][ct_AN][rl][k][i];
	}
	be2[i] = be2[i] + 2.0*sum;
      }
    }

    /****************************************************
                         Averaging
    ****************************************************/

    sum = 0.0;
    for (i=0; i<=tno1; i++){
      sum = sum + be2[i];
    }
    sum = sum/(double)Avnum;
    for (i=0; i<=tno1; i++){
      be2[i] = sum;
    }

    Avnum = 0;    
    sum = 0.0;
    for (rl=RNUM[ct_AN]; (RNUM[ct_AN]-itnum+1)<=rl; rl--){
      Avnum = Avnum + RCN_Rank[ct_AN][rl];
      for (i=0; i<RCN_Rank[ct_AN][rl]; i++){  
	sum = sum + al[spin][ct_AN][rl][i][i];
      }
    }
    Av_al = sum/(double)Avnum;

    /****************************************************
         Calculation of the square root terminator
    ****************************************************/

    for (i=0; i<RCN_Rank[ct_AN][RNUM[ct_AN]]; i++){
      d1 = EpP.r - Av_al;
      d2 = EpP.i;
      xd = d1*d1 - d2*d2 - 2.0*be2[i];
      yd = 2.0*d1*d2;
      dum = sqrt(xd*xd+yd*yd);
      c1 = xd/dum;
      s1 = yd/dum;
      if (yd<0.0){
	c2 = -sqrt((1.0+c1)*0.5);
	s2 = sqrt(1.0-c2*c2);
      }
      else{
	c2 = sqrt((1.0+c1)*0.5);
	s2 = sqrt(1.0-c2*c2);
      }
          
      dum = sqrt(dum);
      d1 = d1 - dum*c2;
      d2 = d2 - dum*s2;
      cter[i].r = d1/be2[i];
      cter[i].i = d2/be2[i];
    }          
  }

  /****************************************************
            Calculation of multiple inverse
  ****************************************************/

  for (i=0; i<=tno1; i++){
    for (j=0; j<=tno1; j++){
      if (i==j)
	cinv[i][i] = cter[i];
      else{
	cinv[i][j].r = 0.0;
	cinv[i][j].i = 0.0;
      }
    }
  }

  if (T_switch==1)
    rl_loop = RNUM[ct_AN];
  else
    rl_loop = RNUM[ct_AN] - 1;

  for (rl=rl_loop; 0<=rl; rl--){
    if (rl==RNUM[ct_AN]){
      for (i=0; i<RCN_Rank[ct_AN][rl]; i++){
	for (j=0; j<RCN_Rank[ct_AN][rl]; j++){
	  if (i==j){
	    cinv[i][i].r = EpP.r - al[spin][ct_AN][rl][i][i];
	    cinv[i][i].i = EpP.i;
	  }
	  else{
	    cinv[i][j].r = -al[spin][ct_AN][rl][i][j];
	    cinv[i][j].i = 0.0;
	  }
	}
      }
    }
    else{

      /****************************************************
                           Bn+1*[]^-1
      ****************************************************/

      for (i=0; i<RCN_Rank[ct_AN][rl]; i++){
	for (j=0; j<RCN_Rank[ct_AN][rl+1]; j++){
	  Csum3.r = 0.0;
	  Csum3.i = 0.0;
	  for (k=0; k<RCN_Rank[ct_AN][rl+1]; k++){
	    cdum.r = be[spin][ct_AN][rl+1][i][k]*cinv[k][j].r;
	    cdum.i = be[spin][ct_AN][rl+1][i][k]*cinv[k][j].i;
	    Csum3.r = Csum3.r + cdum.r;
	    Csum3.i = Csum3.i + cdum.i;
	  } 
	  ctp1[i][j].r = Csum3.r;
	  ctp1[i][j].i = Csum3.i;
	}
      }

      /****************************************************
                        Bn+1 * []^-1 * Cn+1
      ****************************************************/

      for (i=0; i<RCN_Rank[ct_AN][rl]; i++){
	for (j=0; j<RCN_Rank[ct_AN][rl]; j++){
	  Csum3.r = 0.0;
	  Csum3.i = 0.0;
	  for (k=0; k<RCN_Rank[ct_AN][rl+1]; k++){
	    cdum.r = ctp1[i][k].r*ga[spin][ct_AN][rl+1][k][j];
            cdum.i = ctp1[i][k].i*ga[spin][ct_AN][rl+1][k][j];
            Csum3.r = Csum3.r + cdum.r;
	    Csum3.i = Csum3.i + cdum.i;
          }
	  ctp2[i][j].r = Csum3.r;
	  ctp2[i][j].i = Csum3.i;
	}
      }

      /****************************************************
                    Z - Alpha - Bn+1*[]^-1*Cn+1
      ****************************************************/

      for (i=0; i<RCN_Rank[ct_AN][rl]; i++){
	for (j=0; j<RCN_Rank[ct_AN][rl]; j++){
	  if (i==j){
	    cinv[i][i].r = EpP.r - al[spin][ct_AN][rl][i][j] - ctp2[i][j].r;
	    cinv[i][i].i = EpP.i - ctp2[i][j].i;
	  }
	  else{
	    cinv[i][j].r = -al[spin][ct_AN][rl][i][j] - ctp2[i][j].r;
	    cinv[i][j].i = -ctp2[i][j].i;
	  }
	}
      }
    }

    /****************************************************
                 [Z-Alpha-Beta*[]^-1*Beta]^-1
    ****************************************************/

    LU_inverse(RCN_Rank[ct_AN][rl]-1,cinv);

    /****************************************************
                         rl loop
    ****************************************************/

  }

  for (i=0; i<=tno1; i++){
    for (j=0; j<=tno1; j++){
      G00S[i][j] = cinv[i][j];
    }
  }

  /****************************************************
    allocation of arrays:
  ****************************************************/

  free(be2);
  free(cter);

  for (i=0; i<List_YOUSO[7]; i++){
    free(cinv[i]);
  }
  free(cinv);

  for (i=0; i<List_YOUSO[7]; i++){
    free(ctp1[i]);
  }
  free(ctp1);

  for (i=0; i<List_YOUSO[7]; i++){
    free(ctp2[i]);
  }
  free(ctp2);

}


void Save_Recursion()
{
  static char fileRN[YOUSO10] = ".rcn";
  static FILE *fp_Rcn;

  strcpy(fileRN,".rcn");
  fnjoint(filepath,filename,fileRN);
  if ((fp_Rcn = fopen(fileRN,"w")) != NULL){
    Output_RcnCof(fp_Rcn);
    Output_LU(fp_Rcn);
    fclose(fp_Rcn);
  }
  else
    printf("Failure of saving the recusion logfile.\n");
}

void Output_LU(FILE *fp_Rcn)
{
  static int ct_AN,rl,i,j,wan1,TNO1,wan2,TNO2;
  static int spin,Gh_AN,h_AN;

  for (spin=0; spin<=SpinP_switch; spin++){

    if (spin==0){
      fprintf(fp_Rcn,"******************************************\n");
      fprintf(fp_Rcn,"   Up spin  - Left side Lanczos vectors   \n");
      fprintf(fp_Rcn,"******************************************\n");
    }
    else if (spin==1){
      fprintf(fp_Rcn,"*****************************************\n");
      fprintf(fp_Rcn,"  Down spin - Left side Lanczos vectors  \n");
      fprintf(fp_Rcn,"*****************************************\n");
    }

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      wan1 = WhatSpecies[ct_AN];
      TNO1 = Spe_Total_CNO[wan1];
      for (rl=0; rl<=RNUM[ct_AN]; rl++){
        fprintf(fp_Rcn,"%i %i\n",ct_AN,rl);
        for (h_AN=0; h_AN<=(FNAN[ct_AN]+SNAN[ct_AN]); h_AN++){
          Gh_AN = natn[ct_AN][h_AN];
          wan2 = WhatSpecies[Gh_AN];
          TNO2 = Spe_Total_CNO[wan2];

          for (j=0; j<TNO2; j++){
            for (i=0; i<TNO1; i++){
              fprintf(fp_Rcn,"%19.14f ",LU[spin][ct_AN][rl][h_AN][i][j]);
            }
            fprintf(fp_Rcn,"\n");
          }

	  /*
          for (i=0; i<TNO1; i++){
            for (j=0; j<TNO2; j++){
              fprintf(fp_Rcn,"%19.14f ",LU[spin][ct_AN][rl][h_AN][i][j]);
            }
            fprintf(fp_Rcn,"\n");
          }
	  */

        }
      }   
    }

  }

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
