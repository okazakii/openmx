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
#include "abred_common.h"

static void TSB_Lanczos(int ct_AN,int spin,double HP[YOUSO1][YOUSO2][YOUSO7][YOUSO7]);
static void Solve_EigenValue(int ct_AN, int spin);
static void Search_ChemP();
static void Calc_DM(double CDM[YOUSO23][YOUSO1][YOUSO8][YOUSO7][YOUSO7],
		    double EDM[YOUSO23][YOUSO1][YOUSO8][YOUSO7][YOUSO7],
		    double IOLP[YOUSO1][YOUSO2][YOUSO7][YOUSO7]);
static int SVD_fact_inv(int rl, int n, double a[YOUSO7][YOUSO7],
                        double B[YOUSO7][YOUSO7],double IB[YOUSO7][YOUSO7],
                        double C[YOUSO7][YOUSO7],double IC[YOUSO7][YOUSO7]);

static void T_Conserve_MP(int T_switch);
static void BO_Calc_S(int T_switch,
		      double CDM[YOUSO23][YOUSO1][YOUSO8][YOUSO7][YOUSO7],
		      double EDM[YOUSO23][YOUSO1][YOUSO8][YOUSO7][YOUSO7],
		      double IOLP[YOUSO1][YOUSO2][YOUSO7][YOUSO7]);
static void GBTD00(int T_switch, int ct_AN, int spin,
		   dcomplex EpP, dcomplex G00S[YOUSO7][YOUSO7]);
static void RecurG(int ct_AN, int spin, dcomplex EpP,
		   dcomplex G00S[YOUSO7][YOUSO7],
		   dcomplex G0S[YOUSO3][YOUSO7][YOUSO7]);
static void ReLU_fact_inv(int n,double a[YOUSO7][YOUSO7],
			  double B[YOUSO7][YOUSO7],double IB[YOUSO7][YOUSO7],
			  double C[YOUSO7][YOUSO7],double IC[YOUSO7][YOUSO7]);
static void MoLU_fact_inv(int n,double a[YOUSO7][YOUSO7],
			  double B[YOUSO7][YOUSO7],double IB[YOUSO7][YOUSO7],
			  double C[YOUSO7][YOUSO7],double IC[YOUSO7][YOUSO7]);

static void ISH(double IOLP[YOUSO1][YOUSO2][YOUSO7][YOUSO7],
		double nh[YOUSO1][YOUSO8][YOUSO7][YOUSO7],
		double HP[YOUSO1][YOUSO2][YOUSO7][YOUSO7]);
static void SbyISH(double S[YOUSO1][YOUSO8][YOUSO7][YOUSO7],
		   double HP[YOUSO1][YOUSO2][YOUSO7][YOUSO7],
		   double nh[YOUSO1][YOUSO8][YOUSO7][YOUSO7],
		   double DH[YOUSO1][YOUSO2][YOUSO7][YOUSO7]);
static void Save_Recursion();
static void Output_RcnCof(FILE *fp);
static void Output_LU(FILE *fp);

static double al[YOUSO23][YOUSO1][YOUSO3][YOUSO7][YOUSO7];
static double be[YOUSO23][YOUSO1][YOUSO3][YOUSO7][YOUSO7];
static double ibe[YOUSO23][YOUSO1][YOUSO3][YOUSO7][YOUSO7];
static double ga[YOUSO23][YOUSO1][YOUSO3][YOUSO7][YOUSO7];
static double iga[YOUSO23][YOUSO1][YOUSO3][YOUSO7][YOUSO7];
static double LU[YOUSO23][YOUSO1][YOUSO3+1][YOUSO2][YOUSO7][YOUSO7];
static double LEVec[YOUSO23][YOUSO1][(YOUSO3+1)*YOUSO7][YOUSO3+1][YOUSO7];
static double REVec[YOUSO23][YOUSO1][(YOUSO3+1)*YOUSO7][YOUSO3+1][YOUSO7];
static double EigenValue[YOUSO23][YOUSO1][(YOUSO3+1)*YOUSO7];
static double Uele_OS[YOUSO23],Uele_IS[YOUSO23];
static int RCN_Rank[YOUSO1][YOUSO3];


static double TmpVec[YOUSO3+1][YOUSO7];


double RecursionS_D(int SCF_iter,
                 double Ham[YOUSO23][YOUSO1][YOUSO8][YOUSO7][YOUSO7],
                 double S[YOUSO1][YOUSO8][YOUSO7][YOUSO7],
                 double CDM[YOUSO23][YOUSO1][YOUSO8][YOUSO7][YOUSO7],
                 double EDM[YOUSO23][YOUSO1][YOUSO8][YOUSO7][YOUSO7],
                 double Eele0[2], double Eele1[2])
{
  static int firsttime=1;
  static int i,ct_AN,spin;
  static double time0;
  static int m,n,j,Rn,k,l,rl;
  static double sum;
  static double IOLP[YOUSO1][YOUSO2][YOUSO7][YOUSO7];
  static double HPS[YOUSO23][YOUSO1][YOUSO2][YOUSO7][YOUSO7];
  static double DH[YOUSO23][YOUSO1][YOUSO2][YOUSO7][YOUSO7];
  static double TStime,TEtime;

  dtime(&TStime);
  Timetool_times("RecursionS_D","start");
  if (firsttime) {
  PrintMemory("RecursionS_D: al",sizeof(al),NULL);
  PrintMemory("RecursionS_D: be",sizeof(be),NULL);
  PrintMemory("RecursionS_D: ibe",sizeof(ibe),NULL);
  PrintMemory("RecursionS_D: ga",sizeof(ga),NULL);
  PrintMemory("RecursionS_D: iga",sizeof(iga),NULL);
  PrintMemory("RecursionS_D: LU",sizeof(LU),NULL);
  PrintMemory("RecursionS_D: LEVec",sizeof(LEVec),NULL);
  PrintMemory("RecursionS_D: REVec",sizeof(REVec),NULL);
  PrintMemory("RecursionS_D: EigenValue",sizeof(EigenValue),NULL);
  PrintMemory("RecursionS_D: RCN_Rank",sizeof(RCN_Rank),NULL);
  PrintMemory("RecursionS_D: IOLP",sizeof(IOLP),NULL);
  PrintMemory("RecursionS_D: HPS",sizeof(HPS),NULL);
  PrintMemory("RecursionS_D: DH",sizeof(DH),NULL);
  firsttime=0;
  }



  /****************************************************
     Calculate the inverse matrix of overlap matrix
  ****************************************************/

  if (IS_switch==1) IS_Lanczos(S,IOLP,rlmax_IS);
  if (IS_switch==2) IS_LU(S,IOLP);
  if (IS_switch==3) IS_Hotelling(S,IOLP,rlmax_IS);
  if (IS_switch==4) IS_Taylor(S,IOLP,rlmax_IS);

  /****************************************************
                   Calculate S^{-1}H
  ****************************************************/
  
  for (spin=0; spin<=SpinP_switch; spin++){
    ISH(IOLP,Ham[spin],HPS[spin]);
  }

  /****************************************************
          Two-sided block Lanczos transformation
  ****************************************************/

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    for (spin=0; spin<=SpinP_switch; spin++){
      TSB_Lanczos(ct_AN,spin,HPS[spin]);
    }
  }

  /****************************************************
     Solve the eigenvalue problem of the matrix block
              tri-diagonalized by the TSB
  ****************************************************/

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){

    printf("Solve_EigenValue ct_AN=%2d\n",ct_AN); 

    for (spin=0; spin<=SpinP_switch; spin++){
      Solve_EigenValue(ct_AN,spin);
    }
  }
  
  /****************************************************
          Search the chemical potential to maintain
                the number of electrons
  ****************************************************/

  Search_ChemP();

  /****************************************************
             Calculate the density matrix
  ****************************************************/

  Calc_DM(CDM,EDM,IOLP);

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
  
  dtime(&TEtime);
  Timetool_times("RecursionS_D","end");
  time0 = TEtime - TStime;
  return time0;
}

void Calc_DM(double CDM[YOUSO23][YOUSO1][YOUSO8][YOUSO7][YOUSO7],
	     double EDM[YOUSO23][YOUSO1][YOUSO8][YOUSO7][YOUSO7],
             double IOLP[YOUSO1][YOUSO2][YOUSO7][YOUSO7])
{
  static int firsttime=1;
  static int ct_AN,spin,N;
  static int wan,tno,tno1,itnum,ig,ino1,jg,jno1,kg,kj,kno1;
  static int i,j,k,p,l,m,n,q,rl,i1,j1,h_AN;
  static double tmp,sum,sum2,sumx,sumy,sumz,dtnum;
  static double Rho[YOUSO23][YOUSO3][YOUSO7][YOUSO7];
  static double ERho[YOUSO23][YOUSO3][YOUSO7][YOUSO7]; 
  static double Uele_temp[YOUSO23][YOUSO1];
  static double Dmat[YOUSO7][YOUSO7],Dmat2[YOUSO7][YOUSO7];
  static double x,FermiF,dum,dum1,dum2;

  if (firsttime) {
  PrintMemory("Calc_DM: Rho",sizeof(Rho),NULL);
  PrintMemory("Calc_DM: ERho",sizeof(ERho),NULL);
  firsttime=0;
  }


  for (spin=0; spin<=SpinP_switch; spin++){
    Uele_IS[spin] = 0.0;
  }

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    wan = WhatSpecies[ct_AN];
    tno = Spe_Total_CNO[wan];
    tno1 = tno - 1;
    N = tno*(RNUM[ct_AN] + 1);
    for (spin=0; spin<=SpinP_switch; spin++){
      for (rl=0; rl<=RNUM[ct_AN]; rl++){
        for (i=0; i<=tno1; i++){
	  for (j=0; j<=tno1; j++){

            sum  = 0.0;
            sum2 = 0.0;
            for (k=0; k<N; k++){
              x = (EigenValue[spin][ct_AN][k] - ChemP)*Beta;
              if (x<=-30.0) x = -30.0;
              if (30.0<=x)  x =  30.0;
              FermiF = 1.0/(1.0 + exp(x));

              tmp = REVec[spin][ct_AN][k][0][i]*LEVec[spin][ct_AN][k][rl][j];
              sum  = sum  + FermiF*tmp;
              sum2 = sum2 + FermiF*EigenValue[spin][ct_AN][k]*tmp;
	    }

            Rho[spin][rl][i][j]  = sum;
            ERho[spin][rl][i][j] = sum2;

	  }
        }
      }
    }

    for (spin=0; spin<=SpinP_switch; spin++){
      sum = 0.0;
      for (i=0; i<=tno1; i++) sum = sum + ERho[spin][0][i][i];
      Uele_temp[spin][ct_AN] = sum;
    }

    for (rl=0; rl<=RNUM[ct_AN]; rl++){
      printf("Rho  ct_AN=%2d rl=%2d\n",ct_AN,rl);
      sum = 0.0;   
      for (i=0; i<=tno1; i++){
        for (k=0; k<=tno1; k++){
          sum = sum + fabs(Rho[0][rl][i][k]);
          printf("%15.12f ",Rho[0][rl][i][k]);
	}
        printf("\n");
      }
      printf("ct_AN=%2d rl=%2d  Sum of Rho = %18.15f\n",ct_AN,rl,sum); 
    }

    /****************************************************
          Transformation of DM and EDM from L basis
                  into the original basis
    ****************************************************/

    for (spin=0; spin<=SpinP_switch; spin++){
      for (j=0; j<=FNAN[ct_AN]; j++){
        jg = natn[ct_AN][j];
        jno1 = Spe_Total_CNO[WhatSpecies[jg]] - 1;
        for (rl=0; rl<=RNUM[ct_AN]; rl++){
          for (k=0; k<=(FNAN[ct_AN]+SNAN[ct_AN]); k++){
            kg = natn[ct_AN][k];
            kj = RMI2[ct_AN][k][j];
            kno1 = Spe_Total_CNO[WhatSpecies[kg]] - 1;

            for (p=0; p<=tno1; p++){
              for (q=0; q<=kno1; q++){
	        sum = 0.0;
                sum2 = 0.0;
	        for (m=0; m<RCN_Rank[ct_AN][rl]; m++){
                  sum  = sum  +  Rho[spin][rl][p][m]*LU[spin][ct_AN][rl][k][m][q]; 
                  sum2 = sum2 + ERho[spin][rl][p][m]*LU[spin][ct_AN][rl][k][m][q]; 
                }
                 Dmat[p][q] = sum;
                Dmat2[p][q] = sum2;
              }
            }

            if (0<=kg){
	      for (p=0; p<=tno1; p++){
		for (q=0; q<=jno1; q++){                
		  sum = 0.0;
		  sum2 = 0.0;
		  for (m=0; m<=kno1; m++){
		    sum  = sum  +  Dmat[p][m]*IOLP[kg][kj][m][q];
		    sum2 = sum2 + Dmat2[p][m]*IOLP[kg][kj][m][q];
		  }
		  if (rl==0 && k==0){
		    CDM[spin][ct_AN][j][p][q] = sum;    
		    EDM[spin][ct_AN][j][p][q] = sum2;
		  }
		  else{
		    CDM[spin][ct_AN][j][p][q] = CDM[spin][ct_AN][j][p][q] + sum;
		    EDM[spin][ct_AN][j][p][q] = EDM[spin][ct_AN][j][p][q] + sum2;
		  }
		}
	      }
	    }

          }
        }
      }
    }

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

  /*
  if (SpinP_switch==0){
    printf("Uele_OS[0]=%15.12f   Uele_OS[1]=%15.12f\n",Uele_OS[0],Uele_OS[0]);
    printf("Uele_IS[0]=%15.12f   Uele_IS[1]=%15.12f\n",Uele_IS[0],Uele_IS[0]);
  }
  else if (SpinP_switch==1){
    printf("Uele_OS[0]=%15.12f   Uele_OS[1]=%15.12f\n",Uele_OS[0],Uele_OS[1]);
    printf("Uele_IS[0]=%15.12f   Uele_IS[1]=%15.12f\n",Uele_IS[0],Uele_IS[1]);
  }
  */

}

void Search_ChemP()
{
  static int ct_AN,spin,wan,tno,N,i,po,m;
  static double ChemP_MAX,ChemP_MIN,Tnum;
  static double TZ,Dnum,FermiF,x,Lnum;

  ChemP_MAX = 10.0;  
  ChemP_MIN =-10.0;

  TZ = 0.0;
  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    wan = WhatSpecies[ct_AN];
    TZ = TZ + Spe_Core_Charge[wan];
  }

  po = 0;
  do {
    ChemP = 0.50*(ChemP_MAX + ChemP_MIN);
    Tnum = 0.0;

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
          for (m=0; m<tno; m++){
            Lnum = Lnum + FermiF*REVec[spin][ct_AN][i][0][m]
                                *LEVec[spin][ct_AN][i][0][m];  
	  }                            
        } 
      }
      Tnum = Tnum + Lnum;
      printf("ct_AN=%2d  Lnum=%15.12f\n",ct_AN,Lnum);

    }

    Dnum = TZ - Tnum;
    printf("ChemP=%16.13f  Tnum=%16.13f  TZ=%16.13f\n",eV2Hartree*ChemP,Tnum,TZ);

    if (0.0<=Dnum) ChemP_MIN = ChemP;
    else           ChemP_MAX = ChemP;
    if (fabs(Dnum)<1.0e-12) po = 1;
  }
  while (po==0); 

}


void Solve_EigenValue(int ct_AN, int spin)
{

  static int rl,i,j,i0,j0,ig,jg,k,m,N,wan,tno,tno1;

  static int k0,rl0,l0,l;

  static double tmp0,tmp1,tmp2;
  static double *A,*A2;
  static double *wr;
  static double *wr2;
  static double *wi;
  static double *vl;
  static double *vr;
  static double *work;
  static double *jun;

  static long int n,lda,ldvl,ldvr,lwork,info;
  static char jobvl = 'V';
  static char jobvr = 'V';

  wan = WhatSpecies[ct_AN];
  tno = Spe_Total_CNO[wan];
  tno1 = tno - 1;
  N = tno*(RNUM[ct_AN] + 1);

  A = (double*)malloc(sizeof(double)*N*N);

  A2 = (double*)malloc(sizeof(double)*N*N);

  wr = (double*)malloc(sizeof(double)*N);
  wr2 = (double*)malloc(sizeof(double)*(N+1));
  wi = (double*)malloc(sizeof(double)*N);
  vl = (double*)malloc(sizeof(double)*N*N);
  vr = (double*)malloc(sizeof(double)*N*N);
  work = (double*)malloc(sizeof(double)*N*4);
  jun = (double*)malloc(sizeof(double)*(N+1));

  n = N;
  lda = N;
  ldvl = N;
  ldvr = N;
  lwork = 4*N;

  for (i=0; i<N*N; i++){
    A[i] = 0.0;
    A2[i] = 0.0;
  }

  /* Setting of al (=A) */

  for (rl=0; rl<=RNUM[ct_AN]; rl++){
    i0 = rl*tno;
    j0 = rl*tno;
    for (i=0; i<=tno1; i++){
      for (j=0; j<=tno1; j++){
        ig = i0 + i;
        jg = j0 + j;
        k = N*jg + ig;
        A[k] = al[spin][ct_AN][rl][i][j]; 

        A2[k] = al[spin][ct_AN][rl][i][j]; 

      }
    }
  }

  /* Setting of be (=B) */

  for (rl=1; rl<=RNUM[ct_AN]; rl++){
    i0 = (rl-1)*tno;
    j0 = rl*tno;
    for (i=0; i<=tno1; i++){
      for (j=0; j<=tno1; j++){
        ig = i0 + i;
        jg = j0 + j;
        k = N*jg + ig;
        A[k] = be[spin][ct_AN][rl][i][j]; 

        A2[k] = be[spin][ct_AN][rl][i][j]; 


      }
    }
  }

  /* Setting of ga (=C) */

  for (rl=1; rl<=RNUM[ct_AN]; rl++){
    i0 = rl*tno;
    j0 = (rl-1)*tno;
    for (i=0; i<=tno1; i++){
      for (j=0; j<=tno1; j++){
        ig = i0 + i;
        jg = j0 + j;
        k = N*jg + ig;
        A[k] = ga[spin][ct_AN][rl][i][j]; 

        A2[k] = ga[spin][ct_AN][rl][i][j]; 


      }
    }
  }
  
  /* Call dgeev_() */

  F77_NAME(dgeev,DGEEV)( &jobvl, &jobvr, &n, A, &lda, wr, wi, vl,
          &ldvl, vr, &ldvr, work, &lwork, &info);

  if (info!=0){
    printf("Errors in dgeev_() info=%2d\n",info);
  }

  for (i=0; i<N; i++){
    jun[i+1] = (double)(i+1);
    if (fabs(wr[i])<1.0e-14) wr2[i+1] = 1.0e+10;
    else                     wr2[i+1] = wr[i]; 
  }
  
  QuickSort2(N,wr2,jun);  
  
  for (i=0; i<N; i++){
    EigenValue[spin][ct_AN][i] = wr2[i+1];
    wr[i] = wr2[i+1];
  }
  
  for (j=0; j<N; j++){
    j0 = (int)(jun[j+1] + 0.00001) - 1;
    
    for (i=0; i<N; i++){
      k = N*j0 + i;
      m  = i%tno;  
      rl = (i-m)/tno;         
      LEVec[spin][ct_AN][j][rl][m] = vl[k]; 
      REVec[spin][ct_AN][j][rl][m] = vr[k]; 
    }
  }

  /*
  for (j=0; j<N; j++){
    for (i=0; i<N; i++){
      k = N*j + i;
      printf("j0= j=%3d i=%3d  L=%18.15f  R=%18.15f\n",j,i,vl[k],vr[k]);
    }
  }
  */

  /* Correction for a set of conjugate complex eigenvectors */

  for (i=0; i<(N-1); i++){
    i0 = (int)(jun[i+1] + 0.000001) - 1; 
    j0 = (int)(jun[i+2] + 0.000001) - 1;

    if ((wi[i0]<0.0 || wi[j0]<0.0) && fabs(wr[i]-wr[i+1])<1.0e-14){
      for (rl=0; rl<=RNUM[ct_AN]; rl++){
        for (k=0; k<=tno1; k++){
          tmp0 = REVec[spin][ct_AN][i][rl][k];
          tmp1 = REVec[spin][ct_AN][i+1][rl][k];
          REVec[spin][ct_AN][i][rl][k]   = tmp1;
          REVec[spin][ct_AN][i+1][rl][k] = tmp0;
        }
      }
    }
  }

  /* Normalize the eigenvectors */

  for (i=0; i<N; i++){
    tmp0 = 0.0;
    for (rl=0; rl<=RNUM[ct_AN]; rl++){
      for (k=0; k<=tno1; k++){
        tmp0 = tmp0 + LEVec[spin][ct_AN][i][rl][k]*REVec[spin][ct_AN][i][rl][k];
      }

    }
    if (0.0<tmp0){
      tmp1 = 1.0/sqrt(tmp0); 
      tmp2 = tmp1;
    }
    else{
      tmp1 = 1.0/sqrt(fabs(tmp0)); 
      tmp2 = -tmp1;
    } 
 
    for (rl=0; rl<=RNUM[ct_AN]; rl++){
      for (k=0; k<=tno1; k++){
        LEVec[spin][ct_AN][i][rl][k] = tmp1*LEVec[spin][ct_AN][i][rl][k];
        REVec[spin][ct_AN][i][rl][k] = tmp2*REVec[spin][ct_AN][i][rl][k];
      }
    }
  }


  for (i=0; i<N; i++){
    i0 = (int)(jun[i+1] + 0.000001) - 1; 
    printf("i=%2d  wr=%18.15f  wi=%18.15f\n",i,wr[i],wi[i0]);  
  }

  /*
  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      tmp0 = 0.0;
      for (rl=0; rl<=RNUM[ct_AN]; rl++){
        for (k=0; k<=tno1; k++){
          tmp0 = tmp0 + LEVec[spin][ct_AN][i][rl][k]*REVec[spin][ct_AN][j][rl][k];
        }
      }
      printf("i=%3d j=%3d Norm=%18.15f\n",i,j,tmp0);
    }
  }
  */


  /*
  for (rl=0; rl<=RNUM[ct_AN]; rl++){
    for (i=0; i<=tno1; i++){
      printf("rl=%2d i=%2d L46=%18.15f L47=%18.15f R46=%18.15f R47=%18.15f\n",
               rl,i,
               LEVec[spin][ct_AN][46][rl][i],
               LEVec[spin][ct_AN][47][rl][i],
               REVec[spin][ct_AN][46][rl][i],
               REVec[spin][ct_AN][47][rl][i]
               );
    }
  }
  */

  /*
  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      k0 = -1;
      for (rl0=0; rl0<=RNUM[ct_AN]; rl0++){
        for (k=0; k<=tno1; k++){
          k0++;
          tmp0 = 0.0;
          l0 = -1;
          for (rl=0; rl<=RNUM[ct_AN]; rl++){
            for (l=0; l<=tno1; l++){
              l0++;
              tmp0 = tmp0 + LEVec[spin][ct_AN][i][rl][l]*A2[N*k0+l0];
	    }
	  }
          TmpVec[rl0][k] = tmp0;
	}
      }
      tmp0 = 0.0;
      for (rl=0; rl<=RNUM[ct_AN]; rl++){
        for (l=0; l<=tno1; l++){
          tmp0 = tmp0 + TmpVec[rl][l]*REVec[spin][ct_AN][j][rl][l];
	}
      }       
      printf("i=%3d j=%3d <L|H|R>=%18.15f\n",i,j,tmp0);
    }
  }
  */
  
  /* free arraies */

  free(A2);  


  free(A);  
  free(wr);  
  free(wr2);  
  free(wi);  
  free(vl);  
  free(vr);  
  free(work);  
  free(jun);  

}


void SbyISH(double S[YOUSO1][YOUSO8][YOUSO7][YOUSO7],
            double HP[YOUSO1][YOUSO2][YOUSO7][YOUSO7],
            double nh[YOUSO1][YOUSO8][YOUSO7][YOUSO7],
            double DH[YOUSO1][YOUSO2][YOUSO7][YOUSO7])
{
  static int firsttime=1;
  static int ct_AN,fan,san,can,ig,ian,j,jg,jan;
  static int k,kg,kan,m,n,l,jl;
  static double sum;
  static double AH[YOUSO1][YOUSO2][YOUSO7][YOUSO7];

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

        if (0<=jl){
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

}


void ISH(double IOLP[YOUSO1][YOUSO2][YOUSO7][YOUSO7],
         double nh[YOUSO1][YOUSO8][YOUSO7][YOUSO7],
         double HP[YOUSO1][YOUSO2][YOUSO7][YOUSO7])
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




void TSB_Lanczos(int ct_AN, int spin, 
                 double HP[YOUSO1][YOUSO2][YOUSO7][YOUSO7])
{
  static int firsttime=1;
  static int i,j,k,l,m,n,ie,i1,j1,po,ig,jg,Rni,Rnj,Rn;
  static int kl,l1,l2,l3,ian,jan;
  static int m1,m2,m3,ZeroNum,Rloop_END;
  static int fan,san,can,rl,rl0,wan,ct_on;
  static double sum,sum1,sum2,dum,dum1,sumx,sumy,sumz,xn,xa;
  static double BeGa[YOUSO7][YOUSO7];
  static double RCpre[YOUSO2+1][YOUSO7][YOUSO7];
  static double RCcnt[YOUSO2+1][YOUSO7][YOUSO7];
  static double RRcnt[YOUSO2+1][YOUSO7][YOUSO7];
  static double LCpre[YOUSO2+1][YOUSO7][YOUSO7];
  static double LCcnt[YOUSO2+1][YOUSO7][YOUSO7];
  static double LRcnt[YOUSO2+1][YOUSO7][YOUSO7];
  static double RU[YOUSO3+1][YOUSO2][YOUSO7][YOUSO7];
  static double RR0[YOUSO2+1][YOUSO7][YOUSO7];
  static double LR0[YOUSO2+1][YOUSO7][YOUSO7];
  static double UmRn[YOUSO7][YOUSO7],RnUm[YOUSO7][YOUSO7];

  if (firsttime) {
  PrintMemory("TSB_Lanczos: RCpre",sizeof(RCpre),NULL);
  PrintMemory("TSB_Lanczos: RCcnt",sizeof(RCcnt),NULL);
  PrintMemory("TSB_Lanczos: RRcnt",sizeof(RRcnt),NULL);
  PrintMemory("TSB_Lanczos: LCpre",sizeof(LCpre),NULL);
  PrintMemory("TSB_Lanczos: LCcnt",sizeof(LCcnt),NULL);
  PrintMemory("TSB_Lanczos: LRcnt",sizeof(LRcnt),NULL);
  PrintMemory("TSB_Lanczos: RU",sizeof(RU),NULL);
  PrintMemory("TSB_Lanczos: RR0",sizeof(RR0),NULL);
  PrintMemory("TSB_Lanczos: LR0",sizeof(LR0),NULL);
  firsttime=0;
  }


  fan = FNAN[ct_AN];
  san = SNAN[ct_AN];
  can = fan + san;
  wan = WhatSpecies[ct_AN];
  ct_on = Spe_Total_CNO[wan];

  for (i=0; i<=can; i++){
    for (i1=0; i1<=(YOUSO7-1); i1++){
      for (j1=0; j1<=(YOUSO7-1); j1++){
	RCpre[i][i1][j1] = 0.0;
	RCcnt[i][i1][j1] = 0.0;
	RRcnt[i][i1][j1] = 0.0;

	LCpre[i][i1][j1] = 0.0;
	LCcnt[i][i1][j1] = 0.0;
	LRcnt[i][i1][j1] = 0.0;
      }
    }
  }

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

  for (i=0; i<=(YOUSO7-1); i++){
    for (j=0; j<=(YOUSO7-1); j++){
      be[spin][ct_AN][0][i][j] = 0.0;
      ga[spin][ct_AN][0][i][j] = 0.0;
    }
  }

  /****************************************************
     Recursion of two-sided block Lanczos algorithm
  ****************************************************/

  n = 0; 
  for (i=0; i<=can; i++){
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
      H'|WRn}

      RCcnt[j][m][n] is a column vector.

      j  ->  local index for atom
      m  ->  index for basis
      n  ->  index for basis
    ****************************************************/

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
		LRcnt[i][m][n] = LRcnt[i][m][n] + sum;
	    }
	  } 
	  po++;
	}
      }
    }

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

    if (rl!=RNUM[ct_AN]){

      /****************************************************
             |RRcnt} = H'|WRn} - |WRn-1}Bn - |WRn}An
      ****************************************************/

      for (i=0; i<=can; i++){
	ig = natn[ct_AN][i];
	ian = Spe_Total_CNO[WhatSpecies[ig]];    
	for (m=0; m<=(ian-1); m++){
	  for (n=0; n<RCN_Rank[ct_AN][rl]; n++){
	    sum = 0.0;
	    for (k=0; k<RCN_Rank[ct_AN][rl-1]; k++){
	      sum = sum + RCpre[i][m][k]*be[spin][ct_AN][rl][k][n];
	    }
	    for (k=0; k<RCN_Rank[ct_AN][rl]; k++){
	      sum = sum + RCcnt[i][m][k]*al[spin][ct_AN][rl][k][n];
	    }
	    RRcnt[i][m][n] = RRcnt[i][m][n] - sum;
	  }
	}
      }

      /****************************************************
             {LRcnt| = {WLn|H' - Cn{WLn-1| - An{WLn|
      ****************************************************/

      for (i=0; i<=can; i++){
	ig = natn[ct_AN][i];
	ian = Spe_Total_CNO[WhatSpecies[ig]];    
	for (m=0; m<RCN_Rank[ct_AN][rl]; m++){
	  for (n=0; n<=(ian-1); n++){
	    sum = 0.0;
	    for (k=0; k<RCN_Rank[ct_AN][rl-1]; k++){
	      sum = sum + ga[spin][ct_AN][rl][m][k]*LCpre[i][k][n];
	    }
	    for (k=0; k<RCN_Rank[ct_AN][rl]; k++){
	      sum = sum + al[spin][ct_AN][rl][m][k]*LCcnt[i][k][n];
	    }
	    LRcnt[i][m][n] = LRcnt[i][m][n] - sum;
	  }
	}
      }

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

      /* (LRcnt| := (LRcnt| - ((LR0| := (LR0| + (LRcnt|RU_rl0)(LU_rl0|) */

      for (rl0=0; rl0<=rl; rl0++){

        /* (LRcnt|RU_rl0) */

        for (i=0; i<=can; i++){
	  ig = natn[ct_AN][i];
	  ian = Spe_Total_CNO[WhatSpecies[ig]];
	  for (m=0; m<RCN_Rank[ct_AN][rl]; m++){
	    for (n=0; n<RCN_Rank[ct_AN][rl0]; n++){
	      sum = 0.0;
	      for (k=0; k<=(ian-1); k++){
                sum = sum + LRcnt[i][m][k]*RU[rl0][i][k][n];
	      }
	      if (i==0)
                RnUm[m][n] = sum;
              else
	        RnUm[m][n] = RnUm[m][n] + sum;
	    } 
	  }
	}

        /* (LR0| := (LR0| + (LRcnt|RU_rl0)(LU_rl0| */

	for (i=0; i<=can; i++){
	  ig = natn[ct_AN][i];
	  ian = Spe_Total_CNO[WhatSpecies[ig]];
	  for (m=0; m<RCN_Rank[ct_AN][rl]; m++){
	    for (n=0; n<=(ian-1); n++){
	      sum = 0.0;
	      for (k=0; k<RCN_Rank[ct_AN][rl0]; k++){
                sum = sum + RnUm[m][k]*LU[spin][ct_AN][rl0][i][k][n];
              }
	      if (rl0==0) LR0[i][m][n] = sum;
              else        LR0[i][m][n] = LR0[i][m][n] + sum;
	    }
	  }
	}
      }

      /* (LRcnt| := (LRcnt| - (LR0| */

      for (i=0; i<=can; i++){
	ig = natn[ct_AN][i];
	ian = Spe_Total_CNO[WhatSpecies[ig]]; 
	for (m=0; m<RCN_Rank[ct_AN][rl]; m++){
	  for (n=0; n<=(ian-1); n++){
            LRcnt[i][m][n] = LRcnt[i][m][n] - LR0[i][m][n];           
	  } 
	}
      }

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

      /*
      printf("BeGa  ct_AN=%2d  rl=%2d\n",ct_AN,rl);
      for (m=0; m<RCN_Rank[ct_AN][rl]; m++){
        for (n=0; n<RCN_Rank[ct_AN][rl]; n++){
          printf("%15.12f ",BeGa[m][n]);
	}
        printf("\n");
      }
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

      /****************************************************
                      {WLn+1| = IBn+1{LRcnt|
      ****************************************************/

      for (i=0; i<=can; i++){
	ig = natn[ct_AN][i];
	ian = Spe_Total_CNO[WhatSpecies[ig]];

	for (m=0; m<RCN_Rank[ct_AN][rl]; m++){
	  for (n=0; n<=(ian-1); n++){
	    LCpre[i][m][n] = LCcnt[i][m][n];
	  }
	}

	for (m=0; m<RCN_Rank[ct_AN][rl+1]; m++){
	  for (n=0; n<=(ian-1); n++){
	    sum = 0.0;
	    for (k=0; k<RCN_Rank[ct_AN][rl]; k++){
	      sum = sum + ibe[spin][ct_AN][rl+1][m][k]*LRcnt[i][k][n];
	    }
	    LCcnt[i][m][n] = sum;
            LU[spin][ct_AN][rl+1][i][m][n] = sum;
	  }
	}
      }

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

  /*
  printf("ct_AN=%2d  RNUM=%2d\n",ct_AN,RNUM[ct_AN]);
  */

}


void MoLU_fact_inv(int n,double a[YOUSO7][YOUSO7],
                   double B[YOUSO7][YOUSO7],double IB[YOUSO7][YOUSO7],
                   double C[YOUSO7][YOUSO7],double IC[YOUSO7][YOUSO7])
{
  static int i,j,k,l,m;
  static double w,sum,dum;
  static double x[YOUSO7],y[YOUSO7],ia[YOUSO7][YOUSO7],bb[YOUSO7][YOUSO7];
  static double da[YOUSO7][YOUSO7];

  /****************************************************
                      From 0 to n
  ****************************************************/

  if (n==-1){
    for (i=0; i<=(YOUSO7-1); i++){
      for (j=0; j<=(YOUSO7-1); j++){
	a[i][j] = 0.0;
	B[i][j] = 0.0;
	C[i][j] = 0.0;
      }
    }
  }
  else{
    for (i=0; i<=n; i++){
      for (j=0; j<=n; j++){
	da[i][j] = a[i][j];
	B[i][j] = 0.0;
	C[i][j] = 0.0;
      }
    }

    /****************************************************
            Modified LU factorization, BC = LU 
    ****************************************************/

    for (k=0; k<=n; k++){

      sum = a[k][k];
      for (i=0; i<=(k-1); i++){
        sum = sum - B[k][i]*C[i][k];
      }
      C[k][k] = sqrt(fabs(sum));

      for (i=k; i<=n; i++){
        sum = a[i][k];
        for (l=0; l<=(k-1); l++){
          sum = sum - B[i][l]*C[l][k];
        }
        B[i][k] = sum/C[k][k];
      }

      for (j=(k+1); j<=n; j++){
        sum = a[k][j];
        for (l=0; l<=(k-1); l++){
          sum = sum - B[k][l]*C[l][j];
        }
        C[k][j] = sum/B[k][k];
      }
      
    }


    /****************************************************
      Calculate the inverse IB of B using the lower
      triangle matrix. 

      1). b is set as a unit vector in the Bx = b
      2). Solve the linear equation.
    ****************************************************/

    for (k=0; k<=n; k++){
      for (i=0; i<=n; i++){
	if (i==k)
	  y[i] = 1.0;
	else
	  y[i] = 0.0;
	for (j=0; j<=(i-1); j++){
	  dum = B[i][j]*y[j];
	  y[i] = y[i] - dum;
	}
        y[i] = y[i]/B[i][i];
	IB[i][k] = y[i];
      }
    }
     
    /****************************************************
      Calculate the inverse IC of C using the upper
      triangle matrix. 

      1). c is set as a unit vector in the Cy = c
      2). Solve the linear equation.

      Note: the diagonal elements of the upper triangle
      matrix are not 1.
    ****************************************************/

    for (k=n; 0<=k; k--){
      for (i=n; 0<=i; i--){
	if (i==k)
	  y[i] = 1.0;
	else
	  y[i] = 0.0;
	for (j=n; (i+1)<=j; j--){  
	  dum = C[i][j]*y[j];
	  y[i] = y[i] - dum;
	}
	y[i] = y[i]/C[i][i];
	IC[i][k] = y[i];
      }
    }
  }

}

int SVD_fact_inv(int rl, int n, double a[YOUSO7][YOUSO7],
                 double B[YOUSO7][YOUSO7],double IB[YOUSO7][YOUSO7],
                 double C[YOUSO7][YOUSO7],double IC[YOUSO7][YOUSO7])
{

  static int i,j,j1,k,rankN,po,num;
  static double sum,rank_criterion,rank_criterion2;
  static double a1[YOUSO7*YOUSO7];
  static double vt[YOUSO7*YOUSO7],u[YOUSO7*YOUSO7];
  static double vt1[YOUSO7][YOUSO7],u1[YOUSO7][YOUSO7];
  static double w[YOUSO7+2];
  static double work[YOUSO7*YOUSO7];
  static double w2[YOUSO7+2],iw2[YOUSO7+2];
  static char jobu  = 'A';
  static char jobvt = 'A';
  static long int ROW,COL,lda,ldu,ldvt,lwork,info;

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
  lwork = 5*YOUSO7;
  F77_NAME(dgesvd,DGESVD)(&jobu, &jobvt, &ROW, &COL, a1, &lda, w, u,
          &ldu, vt, &ldvt, work, &lwork, &info);

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

  rank_criterion  = 1.0e-6;
  rank_criterion2 = 1.0e-11;

  /*
  if (rl<5){
    rank_criterion  = 1.0e-6;
    rank_criterion2 = 1.0e-11;
  }
  else if (rl<10){
    rank_criterion  = 1.0e-5;
    rank_criterion2 = 1.0e-8;
  }
  else if (rl<13) {
    rank_criterion  = 1.0e-4;
    rank_criterion2 = 1.0e-7;
  }
  else {
    rank_criterion  = 1.0e-2;
    rank_criterion2 = 1.0e-3;
  }
  */

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

  for (i=0; i<YOUSO7; i++){
    for (j=0; j<YOUSO7; j++){
      B[i][j]  = 0.0;
      IB[i][j] = 0.0;
    }
  }

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

  for (i=0; i<YOUSO7; i++){
    for (j=0; j<YOUSO7; j++){
      C[i][j]  = 0.0;
      IC[i][j] = 0.0;
    }
  }

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

  /*
  for (i=0; i<=n; i++){
    printf("i=%2d  w2=%15.12f w=%15.12f\n",i,w2[i],w[i]); 
  }
  */

  /*
  printf("B\n");
  for (i=0; i<=n; i++){
    for (j=0; j<=n; j++){
      printf("%10.7f ",B[i][j]);
    }
    printf("\n");
  }

  printf("IB\n");
  for (i=0; i<=n; i++){
    for (j=0; j<=n; j++){
      printf("%10.7f ",IB[i][j]);
    }
    printf("\n");
  }

  printf("C\n");
  for (i=0; i<=n; i++){
    for (j=0; j<=n; j++){
      printf("%10.7f ",C[i][j]);
    }
    printf("\n");
  }

  printf("IC\n");
  for (i=0; i<=n; i++){
    for (j=0; j<=n; j++){
      printf("%10.7f ",IC[i][j]);
    }
    printf("\n");
  }
  */

  return rankN;

  /*
  printf("B*IB\n");
  for (i=0; i<=n; i++){
    for (j=0; j<=n; j++){
      sum = 0.0;
      for (k=0; k<=n; k++){
        sum = sum + B[i][k]*IB[k][j];
      }
      printf("%10.7f ",sum);
    }
    printf("\n");
  }

  printf("IB*B\n");
  for (i=0; i<=n; i++){
    for (j=0; j<=n; j++){
      sum = 0.0;
      for (k=0; k<=n; k++){
        sum = sum + IB[i][k]*B[k][j];
      }
      printf("%10.7f ",sum);
    }
    printf("\n");
  }


  printf("C*IC\n");
  for (i=0; i<=n; i++){
    for (j=0; j<=n; j++){
      sum = 0.0;
      for (k=0; k<=n; k++){
        sum = sum + C[i][k]*IC[k][j];
      }
      printf("%10.7f ",sum);
    }
    printf("\n");
  }

  printf("IC*C\n");
  for (i=0; i<=n; i++){
    for (j=0; j<=n; j++){
      sum = 0.0;
      for (k=0; k<=n; k++){
        sum = sum + IC[i][k]*C[k][j];
      }
      printf("%10.7f ",sum);
    }
    printf("\n");
  }

  printf("B*C\n");
  for (i=0; i<=n; i++){
    for (j=0; j<=n; j++){
      sum = 0.0;
      for (k=0; k<=n; k++){
        sum = sum + B[i][k]*C[k][j];
      }
      printf("%10.7f ",sum);
    }
    printf("\n");
  }
  */

}


void ReLU_fact_inv(int n,double a[YOUSO7][YOUSO7],
                   double B[YOUSO7][YOUSO7],double IB[YOUSO7][YOUSO7],
                   double C[YOUSO7][YOUSO7],double IC[YOUSO7][YOUSO7])
{
  static int i,j,k,l,m;
  static double w,sum,dum;
  static double x[YOUSO7],y[YOUSO7],ia[YOUSO7][YOUSO7],bb[YOUSO7][YOUSO7];
  static double da[YOUSO7][YOUSO7];

  /****************************************************
                      From 0 to n
  ****************************************************/

  if (n==-1){
    for (i=0; i<=(YOUSO7-1); i++){
      for (j=0; j<=(YOUSO7-1); j++){
	a[i][j] = 0.0;
      }
    }
  }
  else{
    for (i=0; i<=n; i++){
      for (j=0; j<=n; j++){
	da[i][j] = a[i][j];
      }
    }

    /****************************************************
                       LU factorization
    ****************************************************/

    for (k=0; k<=(n-1); k++){
      w = 1.0/a[k][k];
      for (i=k+1; i<=n; i++){
        a[i][k] = w*a[i][k];
        for (j=k+1; j<=n; j++){
          dum = a[i][k]*a[k][j];
          a[i][j] = a[i][j] - dum;
        }
      }
    }

    /****************************************************
      The lower part of array a[][] is L, where the
      diagonal elements are 1. The upper part of array
      a[][] is U, where the diagonal elements belong to U.
      B = L and C = U because of B*C = L*U.
    ****************************************************/

    for (i=0; i<=n; i++){
      for (j=0; j<=n; j++){
	if (i==j){
	  B[i][i] = 1.0;
	  C[i][i] = a[i][i];
	}
	else if (i<j){

          /****************************************************
                                 C = U
          ****************************************************/

	  B[i][j] = 0.0;
	  C[i][j] = a[i][j];
	}
	else{

          /****************************************************
                                 B = L
          ****************************************************/

	  B[i][j] = a[i][j];
	  C[i][j] = 0.0;
	}
      }
    }

    /****************************************************
      Calculate the inverse IB of B using the lower
      triangle matrix. 

      1). b is set as a unit vector in the Bx = b
      2). Solve the linear equation.
    ****************************************************/

    for (k=0; k<=n; k++){
      for (i=0; i<=n; i++){
	if (i==k)
	  y[i] = 1.0;
	else
	  y[i] = 0.0;
	for (j=0; j<=(i-1); j++){
	  dum = a[i][j]*y[j];
	  y[i] = y[i] - dum;
	}
	IB[i][k] = y[i];
      }
    }
     
    /****************************************************
      Calculate the inverse IC of C using the upper
      triangle matrix. 

      1). c is set as a unit vector in the Cy = c
      2). Solve the linear equation.

      Note: the diagonal elements of the upper triangle
      matrix are not 1.
    ****************************************************/

    for (k=n; 0<=k; k--){
      for (i=n; 0<=i; i--){
	if (i==k)
	  y[i] = 1.0;
	else
	  y[i] = 0.0;
	for (j=n; (i+1)<=j; j--){  
	  dum = a[i][j]*y[j];
	  y[i] = y[i] - dum;
	}
	y[i] = y[i]/a[i][i];
	IC[i][k] = y[i];
      }
    }
  }

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
  static dcomplex G00[YOUSO23][YOUSO7][YOUSO7];
  static dcomplex CNd[YOUSO1],CXd[YOUSO1];
  static dcomplex cdum,Csum3;

  if (firsttime) {
  PrintMemory("T_Conserve_MP: G00",sizeof(G00);
  firsttime=0;
  }

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
}



void BO_Calc_S(int T_switch,
               double CDM[YOUSO23][YOUSO1][YOUSO8][YOUSO7][YOUSO7],
	       double EDM[YOUSO23][YOUSO1][YOUSO8][YOUSO7][YOUSO7],
               double IOLP[YOUSO1][YOUSO2][YOUSO7][YOUSO7])
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
  static double Rho[YOUSO23][YOUSO3][YOUSO7][YOUSO7];
  static double E00[YOUSO23][YOUSO7];
  static double ERho[YOUSO23][YOUSO3][YOUSO7][YOUSO7]; 
  static double Fx[YOUSO1],Fy[YOUSO1],Fz[YOUSO1];
  static double Uele_temp[YOUSO23][YOUSO1];
  static double Dmat[YOUSO7][YOUSO7],Dmat2[YOUSO7][YOUSO7];
  static double dum,dum1,dum2;

  static dcomplex G00[YOUSO23][YOUSO7][YOUSO7];
  static dcomplex EpC,EpP,CE,CRho,CRes,CERho;
  static dcomplex G0[YOUSO23][YOUSO3][YOUSO7][YOUSO7];

  if (firsttime) {
  PrintMemory("BO_Calc_S: Rho",sizeof(Rho),NULL);
  PrintMemory("BO_Calc_S: E00",sizeof(E00),NULL);
  PrintMemory("BO_Calc_S: ERho",sizeof(ERho),NULL);
  PrintMemory("BO_Calc_S: G00",sizeof(G00),NULL);
  PrintMemory("BO_Calc_S: G0",sizeof(G0),NULL);
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
      for (i=0; i<=tno1; i++) sum = sum + E00[spin][i];
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
          sum = sum + fabs(Rho[0][rl][i][k]);
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

    for (spin=0; spin<=SpinP_switch; spin++){
      for (j=0; j<=FNAN[ct_AN]; j++){
        jg = natn[ct_AN][j];
        jno1 = Spe_Total_CNO[WhatSpecies[jg]] - 1;
        for (rl=0; rl<=RNUM[ct_AN]; rl++){
          for (k=0; k<=(FNAN[ct_AN]+SNAN[ct_AN]); k++){
            kg = natn[ct_AN][k];
            kj = RMI2[ct_AN][k][j];
            kno1 = Spe_Total_CNO[WhatSpecies[kg]] - 1;

            for (p=0; p<=tno1; p++){
              for (q=0; q<=kno1; q++){
	        sum = 0.0;
                sum2 = 0.0;
	        for (m=0; m<RCN_Rank[ct_AN][rl]; m++){
                  sum  = sum  +  Rho[spin][rl][p][m]*LU[spin][ct_AN][rl][k][m][q]; 
                  sum2 = sum2 + ERho[spin][rl][p][m]*LU[spin][ct_AN][rl][k][m][q]; 
                }
                 Dmat[p][q] = sum;
                Dmat2[p][q] = sum2;
              }
            }

	    if (0<=kj){
	      for (p=0; p<=tno1; p++){
		for (q=0; q<=jno1; q++){                
		  sum = 0.0;
		  sum2 = 0.0;
		  for (m=0; m<=kno1; m++){
		    sum  = sum  +  Dmat[p][m]*IOLP[kg][kj][m][q];
		    sum2 = sum2 + Dmat2[p][m]*IOLP[kg][kj][m][q];
		  }
		  if (rl==0 && k==0){
		    CDM[spin][ct_AN][j][p][q] = sum;    
		    EDM[spin][ct_AN][j][p][q] = sum2;
		  }
		  else{
		    CDM[spin][ct_AN][j][p][q] = CDM[spin][ct_AN][j][p][q] + sum;
		    EDM[spin][ct_AN][j][p][q] = EDM[spin][ct_AN][j][p][q] + sum2;
		  }
		}
	      }
	    }

          }
        }
      }
    }

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

}


void RecurG(int ct_AN, int spin, dcomplex EpP,
            dcomplex G00S[YOUSO7][YOUSO7],
            dcomplex G0S[YOUSO3][YOUSO7][YOUSO7])
{

  static int i,j,k,rl,wan;
  static int tno,tno1;
  static int rlm,rlm2;
  static double dum,dum1,dum2;
  static dcomplex Csum,Csum2,Ctp[YOUSO7][YOUSO7],Ctp1[YOUSO7][YOUSO7];
  static dcomplex Ctp2[YOUSO7][YOUSO7],Ctp3[YOUSO7][YOUSO7];
  static dcomplex Cdum,Cdum2,Cdum3;

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

}



void GBTD00(int T_switch, int ct_AN, int spin,
            dcomplex EpP, dcomplex G00S[YOUSO7][YOUSO7])
{
  static int i,j,k,l,m,n,rl,wan,Avnum;
  static int tno1,tno,itnum,rl_loop;
  static double dtnum,Av_al,d1,d2;
  static double dum,sum,xd,yd;
  static double c1,s1,c2,s2;
  static double tr,ti,ai,bi;
  static double be2[YOUSO7];
  static dcomplex cter[YOUSO7],cinv[YOUSO7][YOUSO7];
  static dcomplex ctp1[YOUSO7][YOUSO7],ctp2[YOUSO7][YOUSO7];
  static dcomplex Csum3,cdum;

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
