/**********************************************************************
  init.c:

     init.c is a subroutine to initialize several parameters at
     the starting point of calculations.

  Log of init.c:

     22/Nov/2001  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "openmx_common.h"
#include "mpi.h"


static void InitV();
static void HGF();
static void Matsubara();
static void MM();
static dcomplex G(dcomplex z);



void init()
{
  int i,p,wsp,wan;
  double dp,dM,dum,dum1;
  double kappa,R,R1,R2;

  /****************************************************
                 Correct_Position_flag
  ****************************************************/

  Correct_Position_flag = 0;

  /****************************************************
                      force flag
  ****************************************************/

  F_Kin_flag    = 1;
  F_NL_flag     = 1;
  F_VNA_flag    = 1;
  F_VEF_flag    = 1;
  F_dVHart_flag = 1;
  F_Vxc_flag    = 1;
  F_U_flag      = 1;
  F_dftD_flag   = 1; /* okuno */

  /****************************************************
                Setting of atomic weight
  ****************************************************/

  for (i=1; i<=atomnum; i++){
    wsp = WhatSpecies[i];

    if (0.0<Spe_AtomicMass[wsp]){
     Gxyz[i][20] = Spe_AtomicMass[wsp];
    }
    else{
      wan = Spe_WhatAtom[wsp];
      Gxyz[i][20] = Atom_Weight[wan];
    }
  }

  InitV();

  /* set Beta */

  Beta = 1.0/kB/E_Temp;

  /****************************************************
   (1+(1+x/n)^n)^{-1}
   Complex poles for the modified Matsubara summation
  ****************************************************/

  /*
  dM = POLES*1.000;
  dp = -1.0;
  dum1 = 0.5*PI/dM;
  Beta = 1.0/kB/E_Temp;

  for (p=0; p<POLES; p++){
    dp = dp + 1.0;
    dum = (2.0*dp+1.0)*dum1;
    zp[p] = Complex(cos(dum),sin(dum));
    Ep[p] = Complex(2.0*dM/Beta*(zp[p].r-1.0),2.0*dM/Beta*zp[p].i);
  }
  */


  /****************************************************
  CONTINUED FRACTION
  zero points and residues for the continued fraction
  expansion of exponential function used in integration
  of Green's functions with Fermi-Dirac distribution
  ****************************************************/

  if (Solver==1) {

    zero_cfrac( POLES, zp, Rp ); 

    /*
    for (p=0; p<(POLES-1); p++){
      printf("%5d %15.12f %15.12f %15.12f\n",p+1,zp[p].i,zp[p+1].i-zp[p].i,Rp[p].r); 
    }
    */

    for (p=0; p<POLES; p++){
      Ep[p].r = 0.0;
      Ep[p].i = zp[p].i/Beta;
    }
  }

  /**********************************************************
                 O(N^2) cluster diagonalization 
  **********************************************************/

  if (Solver==9) {

    ON2_Npoles_f = ON2_Npoles; 

    /************************************
      for calculation of density matrix
    ************************************/

    /* set up poles */

    zero_cfrac( ON2_Npoles, ON2_zp, ON2_Rp ); 

    for (p=0; p<ON2_Npoles; p++){
      ON2_method[p] = 1; 
    }

    /* for the zero-th moment */

    R = 1.0e+10;

    ON2_zp[ON2_Npoles].r = 0.0;
    ON2_zp[ON2_Npoles].i = R;

    ON2_Rp[ON2_Npoles].r = 0.0;
    ON2_Rp[ON2_Npoles].i = 0.5*R;

    ON2_method[ON2_Npoles] = 2; 

    ON2_Npoles++;

    /*
    for (p=0; p<ON2_Npoles; p++){
      printf("p=%5d zp %15.12f %15.12f Rp %15.12f %15.12f\n",
              p,ON2_zp[p].r,ON2_zp[p].i,ON2_Rp[p].r,ON2_Rp[p].i); 
    }
    */

    /*
    for (p=0; p<Npoles_ON2; p++){
      Ep[p].r = 0.0;
      Ep[p].i = zp[p].i/Beta;
    }
    */

    /*********************************************
      for calculation of energy density matrix
    *********************************************/

    /* calculation of kappa */

    kappa = 0.0; 
    for (p=0; p<ON2_Npoles_f; p++){
      kappa += ON2_Rp[p].r; 
    }
    kappa = 4.0*kappa/Beta;

    /* set up poles */

    for (p=0; p<ON2_Npoles_f; p++){
      ON2_zp_f[p] = ON2_zp[p];
      ON2_Rp_f[p] = ON2_Rp[p];
      ON2_method_f[p] = 1; 
    }

    /* for the zero-th order moment */

    if (0){

    R = 1.0e+8;

    ON2_zp_f[ON2_Npoles_f].r = 0.0;
    ON2_zp_f[ON2_Npoles_f].i = R;

    ON2_Rp_f[ON2_Npoles_f].r = 0.25*R*(kappa - R);
    ON2_Rp_f[ON2_Npoles_f].i = 0.25*R*(kappa + R);

    ON2_method_f[ON2_Npoles_f] = 2; 
    ON2_Npoles_f++;

    /* for the 1st order moment */

    ON2_zp_f[ON2_Npoles_f].r = -R;
    ON2_zp_f[ON2_Npoles_f].i = 0.0;

    ON2_Rp_f[ON2_Npoles_f].r = 0.25*R*(R - kappa);
    ON2_Rp_f[ON2_Npoles_f].i = 0.25*R*(R - kappa);

    ON2_method_f[ON2_Npoles_f] = 2; 
    ON2_Npoles_f++;

    }


    if (0){

    R1 = 1.0e+21;
    R2 = 1.0e+7;

    ON2_zp_f[ON2_Npoles_f].r = 0.0;
    ON2_zp_f[ON2_Npoles_f].i = R1;

    ON2_Rp_f[ON2_Npoles_f].r = 0.5*R1*R2;
    ON2_Rp_f[ON2_Npoles_f].i = 0.5*R1*kappa;

    ON2_method_f[ON2_Npoles_f] = 2; 
    ON2_Npoles_f++;

    /* for the 1st order moment */

    ON2_zp_f[ON2_Npoles_f].r = 0.0;
    ON2_zp_f[ON2_Npoles_f].i = R2;

    ON2_Rp_f[ON2_Npoles_f].r =-0.5*R2*R2;
    ON2_Rp_f[ON2_Npoles_f].i = 0.0;

    ON2_method_f[ON2_Npoles_f] = 2; 
    ON2_Npoles_f++;
    }


    R = 1.0e+7;

    ON2_zp_f[ON2_Npoles_f].r = 0.0;
    ON2_zp_f[ON2_Npoles_f].i = R*R;

    ON2_Rp_f[ON2_Npoles_f].r = 0.5*R*R*R*R/(R-1.0);
    ON2_Rp_f[ON2_Npoles_f].i = 0.5*R*R*R*kappa/(R-1.0);

    ON2_method_f[ON2_Npoles_f] = 2; 
    ON2_Npoles_f++;

    /* for the 1st order moment */

    ON2_zp_f[ON2_Npoles_f].r = 0.0;
    ON2_zp_f[ON2_Npoles_f].i = R;

    ON2_Rp_f[ON2_Npoles_f].r =-0.5*R*R*R/(R-1.0);
    ON2_Rp_f[ON2_Npoles_f].i =-0.5*R*kappa/(R-1.0);

    ON2_method_f[ON2_Npoles_f] = 2; 
    ON2_Npoles_f++;

  }


  /*
  {

    double sum;
    double complex alpha,weight;

    sum = 0.0;
 
    for (p=0; p<ON2_Npoles_f; p++){

      if (ON2_method_f[p]==1){
        alpha  = 0.0 + I*(ON2_zp_f[p].i/Beta);
        weight = -2.0*ON2_Rp_f[p].r/Beta*alpha;
      }
      else {
        alpha  = ON2_zp_f[p].r + I*ON2_zp_f[p].i;
        weight = ON2_Rp_f[p].r + I*ON2_Rp_f[p].i;
      }

      sum += creal(weight*1.0/(alpha+2.0));

      printf("%2d %18.15f\n",p,sum);

    }

    printf("\nsum = %18.15f\n",sum);

  }

  MPI_Finalize();
  exit(0);
  */


  /****************************************************
      Matsubara
  ****************************************************/

  /*
  Beta = 1.0/kB/E_Temp;

  for (p=0; p<POLES; p++){

    Rp[p].r = -1.0; 
    Rp[p].i = 0.0;

    Ep[p].r = 0.0;
    Ep[p].i = PI*(2.0*(double)p+1.0)/Beta;
  }
  */


  /*
  
  {
    double f[10],e[10];
    double ChemP,sum,x;

    ChemP = 0.0;

    e[1] = -10.0/eV2Hartree;
    e[2] = -5.0/eV2Hartree;
    e[3] = -2.0/eV2Hartree;
    e[4] =  5.0/eV2Hartree;

    sum = 0.0;

    for (i=1; i<=4; i++){
      x = Beta*(e[i]-ChemP); 

      printf("x=%15.12f f=%15.12f\n",x,1.0/(1.0+exp(x)));

      sum += 1.0/(1.0+exp(x)); 
    }

    printf("Analytic %18.15f\n",sum); 

  }

  printf("\n\n");    
  printf("POLES %2d\n",POLES);    
  */

  /*
  HGF();
  */

  /*
  MM();
  */

  /*
  Matsubara();
  */

  /*
  MPI_Finalize();
  exit(0);
  */

}






void MM()
{
  int p;
  double ChemP,N;
  dcomplex EpC,EpP,CN,CX;
  dcomplex G0; 

  ChemP = 0.0;

  CN.r = 0.0;
  CN.i = 0.0;

  for (p=0; p<POLES; p++){
    EpP.r = ChemP + Ep[p].r;
    EpP.i = Ep[p].i;

    G0 = G(EpP);
    CN.r += (G0.r*zp[p].r - G0.i*zp[p].i)*2.0/Beta;
    CN.i += (G0.r*zp[p].i + G0.i*zp[p].r)*2.0/Beta;
  }

  N = CN.r;

  printf("CN.r=%18.15f  CN.i=%18.15f\n",CN.r,CN.i);
  printf("N=%18.15f\n",N);

}



void Matsubara()
{
  int p;
  double ChemP,N,R,mu0;
  dcomplex EpC,EpP,CN,CX;
  dcomplex G0; 

  ChemP = 0.0;
  R = 1.0e+12;

  CN.r = 0.0;
  CN.i = 0.0;

  for (p=0; p<POLES; p++){
    EpP.r = ChemP + Ep[p].r;
    EpP.i = Ep[p].i;

    G0 = G(EpP);
    CN.r += (G0.r*Rp[p].r - G0.i*Rp[p].i)*2.0/Beta;
    CN.i += (G0.r*Rp[p].i + G0.i*Rp[p].r)*2.0/Beta;
  }

  

  EpP.r = 0.0;
  EpP.i = R;

  mu0 = -R*G(EpP).i;
  N = 0.5*mu0 - CN.r;

  printf("CN.r=%18.15f  CN.i=%18.15f\n",CN.r,CN.i);
  printf("mu0=%18.15f  N=%18.15f\n",mu0,N);
}





void HGF()
{
  int p;
  double ChemP,N,R,mu0;
  dcomplex EpC,EpP,CN,CX;
  dcomplex G0; 

  ChemP = 0.0;
  R = 1.0e+12;

  CN.r = 0.0;
  CN.i = 0.0;

  for (p=0; p<POLES; p++){
    EpP.r = ChemP + Ep[p].r;
    EpP.i = Ep[p].i;

    G0 = G(EpP);
    CN.r += (G0.r*Rp[p].r - G0.i*Rp[p].i)*2.0/Beta;
    CN.i += (G0.r*Rp[p].i + G0.i*Rp[p].r)*2.0/Beta;
  }

  

  EpP.r = 0.0;
  EpP.i = R;

  mu0 = -R*G(EpP).i;
  N = 0.5*mu0 - CN.r;

  printf("CN.r=%18.15f  CN.i=%18.15f\n",CN.r,CN.i);
  printf("mu0=%18.15f  N=%18.15f\n",mu0,N);

}



dcomplex G(dcomplex z)
{
  int i;
  dcomplex e[10];
  dcomplex d,sum,dum;
  
  e[1].r =-10.0/eV2Hartree;  e[1].i = 0.0;
  e[2].r = -5.0/eV2Hartree;  e[2].i = 0.0;
  e[3].r = -2.0/eV2Hartree;  e[3].i = 0.0;
  e[4].r =  5.0/eV2Hartree;  e[4].i = 0.0;

  sum.r = 0.0;
  sum.i = 0.0;
 
  for (i=1; i<=4; i++){
    d = Csub(z,e[i]);
    dum = RCdiv(1.0,d);
    sum.r += dum.r;
    sum.i += dum.i;
  }

  return sum; 
}








void InitV()
{
  /******************************************************* 
   1 a.u.=2.4189*10^-2 fs, 1fs=41.341105 a.u. 
  ********************************************************/

  double ax,ay,az,bx,by,bz,tmp,Wscale,sum,v;
  int j,Mc_AN,Gc_AN,ID,myid,numprocs;

  Wscale = unified_atomic_mass_unit/electron_mass;

  /* MPI */
  MPI_Comm_size(MPI_COMM_WORLD1,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD1,&myid);


  if ( (MD_switch==2 || MD_switch==9 || MD_switch==11 || MD_switch==14 || MD_switch==15)
       && MD_Init_Velocity==0 ){

    if (myid==Host_ID){
      for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){

        v = sqrt(3.0*kB*TempPara[1][2]/(Gxyz[Gc_AN][20]*Wscale*eV2Hartree));

	ax = rnd(1.0);
	ay = rnd(1.0);
	az = rnd(1.0);

        tmp = 1.0/sqrt(ax*ax+ay*ay+az*az);

        ax *= tmp; 
        ay *= tmp; 
        az *= tmp; 
        
	Gxyz[Gc_AN][24] = v*ax;
	Gxyz[Gc_AN][25] = v*ay;
	Gxyz[Gc_AN][26] = v*az;

	Gxyz[Gc_AN][27] = Gxyz[Gc_AN][24];
	Gxyz[Gc_AN][28] = Gxyz[Gc_AN][25];
	Gxyz[Gc_AN][29] = Gxyz[Gc_AN][26];
      }
    }

    /****************
    MPI:  Gxyz
    *****************/

    for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
      ID = G2ID[Gc_AN];
      MPI_Bcast(&Gxyz[Gc_AN][24], 3, MPI_DOUBLE, Host_ID, MPI_COMM_WORLD1);
    }
  }

  else if (MD_Init_Velocity==0) {

    for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
      Gxyz[Gc_AN][24] = 0.0;
      Gxyz[Gc_AN][25] = 0.0;
      Gxyz[Gc_AN][26] = 0.0;
      Gxyz[Gc_AN][27] = 0.0;
      Gxyz[Gc_AN][28] = 0.0;
      Gxyz[Gc_AN][29] = 0.0;
    }
  }

  /*********************************************************** 
   calculate the initial Ukc:

   1 a.u.=2.4189*10^-2 fs, 1fs=41.341105 a.u. 
  ***********************************************************/

  Ukc = 0.0;

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    sum = 0.0;
    for (j=1; j<=3; j++){
      if (atom_Fixed_XYZ[Gc_AN][j]==0){
	sum += Gxyz[Gc_AN][j+23]*Gxyz[Gc_AN][j+23];
      }
    }
    Ukc += 0.5*Gxyz[Gc_AN][20]*Wscale*sum;
  }

}
