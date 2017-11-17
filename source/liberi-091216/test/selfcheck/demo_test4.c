/*----------------------------------------------------------------------
  demo.c
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "eri.h"
#include "demo_sub.h"


static void demo(
  ERI_t *eri,
  ERI_Orbital_t *orb1,
  ERI_Orbital_t *orb2,
  ERI_Orbital_t *orb3,
  ERI_Orbital_t *orb4,
  double out_I4[10][2],
  double *sec
)
{
  int i;
  double scr;
  clock_t clk0, clk1;
  double I4[2];

  clk0 = clock(); 
  ERI_Integral(eri, I4, NULL, orb1, orb2, orb3, orb4, ERI_NOSCREEN);
  printf("    NOSCREEN            INT = %20.15e + I %20.15e\n", I4[0], I4[1]);  
  for (i=0; i<10; i++) {
    scr = 0.05*(double)(i+1);
    ERI_Integral(eri, I4, NULL, orb1, orb2, orb3, orb4, scr);
    printf("    SCR = %12.4f  INT = %20.15e + I %20.15e\n", scr, I4[0], I4[1]);  
    out_I4[i][0] = I4[0];
    out_I4[i][1] = I4[1];
  }
  clk1 = clock(); 
  *sec = (double)(clk1 - clk0)/(double)CLOCKS_PER_SEC;
}


void demo_test4(double I4[10][3], double out_sec[2])
{
  int i;
  double I4_lin[10][2], I4_log[10][2];

  double sec, sec_sum[2];
  clock_t clk0, clk1;
  size_t mem;

  ERI_t *eri_lin, *eri_log;

  /* parameters */
  int lmax, lmax_gl, ngrid, ngl;
  double rmax;
  ERI_Init_Misc info;

  /* basis functions */
  ERI_Orbital_t *orb1, *orb2, *orb3, *orb4;
  const char* pao_file[4] = {
    "./pao/H5.0.pao", 
    "./pao/C4.0.pao",
    "./pao/O5.0.pao", 
    "./pao/Fe5.5.pao" 
  };
  int         pao_l[4]   = { 0, 1, 2, 3 };
  int         pao_m[4]   = { 0,-1, 0, 2 };
  int         pao_n[4]   = { 0, 1, 2, 3 };


  printf("----------------------------------------------------------------\n");
  printf("  TEST #4: SCREENING\n");
  printf("----------------------------------------------------------------\n");
  printf("\n");

  {
    printf("INITIALIZING LIBERI WITH LINEAR-SBT:\n");
  
    lmax    = 15;
    lmax_gl = 10;
    ngrid   = 1024;
    ngl     = 100;
   
    info.sbttype = ERI_SBT_LINEAR;
    info.rmax    = 100.0;
 
    printf("  LMAX    = %5d\n"   , lmax   );
    printf("  LMAX_GL = %5d\n"   , lmax_gl);
    printf("  NGRID   = %5d\n"   , ngrid  );
    printf("  NGL     = %5d\n"   , ngl    );
    printf("  RMAX    = %12.4f\n", info.rmax);
    printf("  SBTTYPE = LINEAR-FSBT\n");
    printf("\n");

    mem = ERI_Required_Size(lmax, lmax_gl, ngrid, ngl, ERI_SH_COMPLEX, &info);

    eri_lin = ERI_Init(lmax, lmax_gl, ngrid, ngl, ERI_SH_COMPLEX, &info);
    if (NULL == eri_lin) {
      fprintf(stderr, "*** ERROR in %s (%d)\n", __FILE__, __LINE__);
      abort();
    }
    printf("  INITIALIZED SUCCESSFULLY.\n");
    printf("  MEMORY USAGE = %d BYTE\n", mem);
    printf("\n");
  }

 
  {
    printf("INITIALIZEING LIBERI WITH LOG-SBT:\n");
  
    lmax    = 15;
    lmax_gl = 10;
    ngrid   = 1024;
    ngl     = 100;

    info.sbttype = ERI_SBT_LOG;
    info.rmax    = 700.0;
    info.rho0    = -10.0;
    info.nq      = 3;
    info.qmin    = 1e-4;
    info.qmax    = 1e-2;
 
    printf("  LMAX    = %5d\n"   , lmax   );
    printf("  LMAX_GL = %5d\n"   , lmax_gl);
    printf("  NGRID   = %5d\n"   , ngrid  );
    printf("  NGL     = %5d\n"   , ngl    );
    printf("  RMAX    = %12.4f\n", info.rmax);
    printf("  SBTTYPE = LOG-FSBT\n");
    printf("  RHO0    = %12.4f\n", info.rho0);
    printf("  NQ      = %5d\n",    info.nq); 
    printf("\n");

    mem = ERI_Required_Size(lmax, lmax_gl, ngrid, ngl, ERI_SH_COMPLEX, &info);

    eri_log = ERI_Init(lmax, lmax_gl, ngrid, ngl, ERI_SH_COMPLEX, &info);
    if (NULL == eri_log) {
      fprintf(stderr, "*** ERROR in %s (%d)\n", __FILE__, __LINE__);
      abort();
    }
    printf("  INITIALIZED SUCCESSFULLY.\n");
    printf("  MEMORY USAGE = %d BYTE\n", mem);
    printf("\n");
  }


  if (0)
  {
    orb1 = Orbital_New_PAO(pao_file[0], pao_l[0], pao_m[0], pao_n[0]);
    orb2 = Orbital_New_PAO(pao_file[1], pao_l[1], pao_m[1], pao_n[1]);
    orb3 = Orbital_New_PAO(pao_file[2], pao_l[2], pao_m[2], pao_n[2]);
    orb4 = Orbital_New_PAO(pao_file[3], pao_l[3], pao_m[3], pao_n[3]);

    if (NULL == orb1 || NULL == orb2 || NULL == orb3 || NULL == orb4) {
      fprintf(stderr, "*** ERROR in %s (%d)\n", __FILE__, __LINE__);
      abort();
    }

    orb1->c[0] = 0.0000; orb1->c[1] = 0.0000; orb1->c[2] = 0.0000;
    orb2->c[0] = 0.0000; orb2->c[1] = 0.5000; orb2->c[2] = 0.0000;
    orb3->c[0] = 0.0000; orb3->c[1] = 0.0000; orb3->c[2] = 1.5000;
    orb4->c[0] = 1.1549; orb4->c[1] = 0.7500; orb4->c[2] = 2.0865;
  

    printf("BASIS FUNCTIONS:\n");
    printf("  1 : FILE=%10s   NODE=%d  L=%d  M=%d\n", pao_file[0], pao_l[0], pao_m[0], pao_n[0]);
    printf("  2 : FILE=%10s   NODE=%d  L=%d  M=%d\n", pao_file[1], pao_l[1], pao_m[1], pao_n[1]);
    printf("  3 : FILE=%10s   NODE=%d  L=%d  M=%d\n", pao_file[2], pao_l[2], pao_m[2], pao_n[2]);
    printf("  4 : FILE=%10s   NODE=%d  L=%d  M=%d\n", pao_file[3], pao_l[3], pao_m[3], pao_n[3]);
    printf("\n");

    printf("CENTERS:\n");
    printf("  1 : (%12.8f, %12.8f, %12.8f)\n", orb1->c[0], orb1->c[1], orb1->c[2]);
    printf("  2 : (%12.8f, %12.8f, %12.8f)\n", orb2->c[0], orb2->c[1], orb2->c[2]);
    printf("  3 : (%12.8f, %12.8f, %12.8f)\n", orb3->c[0], orb3->c[1], orb3->c[2]);
    printf("  4 : (%12.8f, %12.8f, %12.8f)\n", orb4->c[0], orb4->c[1], orb4->c[2]);
    printf("\n");
  }
  
  {
    ngrid = 8192;
    rmax  = 1000.0; 
    orb1 = Orbital_New_GTO(rmax, ngrid, 0, 0, 1.0, 1.0);    
    orb2 = Orbital_New_GTO(rmax, ngrid, 0, 0, 1.0, 1.0);    
    orb3 = Orbital_New_GTO(rmax, ngrid, 0, 0, 1.0, 1.0);    
    orb4 = Orbital_New_GTO(rmax, ngrid, 0, 0, 1.0, 1.0);    
    if (NULL == orb1 || NULL == orb2 || NULL == orb3 || NULL == orb4) {
      fprintf(stderr, "*** ERROR in %s (%d)\n", __FILE__, __LINE__);
      abort();
    }

    orb1->c[0] =  0.5; orb1->c[1] =  0.5; orb1->c[2] = 0.0;
    orb2->c[0] =  0.5; orb2->c[1] = -0.5; orb2->c[2] = 0.0;
    orb3->c[0] = -0.5; orb3->c[1] =  0.5; orb3->c[2] = 0.0;
    orb4->c[0] = -0.5; orb4->c[1] = -0.5; orb4->c[2] = 0.0;

    printf("BASIS FUNCTIONS (ssss):\n");
    printf("  (ssss) 1 : GTO  L=0  M=0   EXPONENT=1.0  NORMALIZATION=1.0\n");
    printf("         2 : GTO  L=0  M=0   EXPONENT=1.0  NORMALIZATION=1.0\n");
    printf("         3 : GTO  L=0  M=0   EXPONENT=1.0  NORMALIZATION=1.0\n");
    printf("         4 : GTO  L=0  M=0   EXPONENT=1.0  NORMALIZATION=1.0\n");
    printf("\n");

    printf("CENTERS:\n");
    printf("  1 : (%12.8f, %12.8f, %12.8f)\n", orb1->c[0], orb1->c[1], orb1->c[2]);
    printf("  2 : (%12.8f, %12.8f, %12.8f)\n", orb2->c[0], orb2->c[1], orb2->c[2]);
    printf("  3 : (%12.8f, %12.8f, %12.8f)\n", orb3->c[0], orb3->c[1], orb3->c[2]);
    printf("  4 : (%12.8f, %12.8f, %12.8f)\n", orb4->c[0], orb4->c[1], orb4->c[2]);
    printf("\n");
  }


  printf("  LINEAR-SBT:\n");
  demo(eri_lin, orb1, orb2, orb3, orb4, I4_lin, &sec);
  sec_sum[0] += sec;
  printf("  LOG-SBT:\n");
  demo(eri_log, orb1, orb2, orb3, orb4, I4_log, &sec);
  sec_sum[1] += sec;
  printf("\n");

  Orbital_Free(orb1);
  Orbital_Free(orb2);
  Orbital_Free(orb3);
  Orbital_Free(orb4);

  ERI_Free(eri_lin);
  ERI_Free(eri_log);

  for (i=0; i<10; i++) {
    I4[i][0] = I4_lin[i][0]; 
    I4[i][1] = I4_log[i][0];
  }
 
  I4[0][2] = 0.00710049115924762;
  I4[1][2] = 0.00659841154368471;
  I4[2][2] = 0.00610852657985739;
  I4[3][2] = 0.00563608579775074;
  I4[4][2] = 0.00518543604762238;
  I4[5][2] = 0.00475988669200648;
  I4[6][2] = 0.00436167497938318;
  I4[7][2] = 0.00399201835529064;
  I4[8][2] = 0.00365123100961191;
  I4[9][2] = 0.00333887870509403;
  
  out_sec[0] = sec_sum[0];
  out_sec[1] = sec_sum[1];
}


/* EOF */
