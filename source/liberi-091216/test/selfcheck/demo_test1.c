/*----------------------------------------------------------------------
  demo_test1.c
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "eri.h"
#include "demo_sub.h"


void demo_test1(double I4[3][3], double out_sec[2])
{
  int i;
  double I4_ssss_lin[2], I4_ssss_log[2];
  double I4_pssp_lin[2], I4_pssp_log[2];
  double I4_dddd_lin[2], I4_dddd_log[2];

  double sec, sec_sum[2];
  clock_t clk0, clk1;
  size_t mem;

  ERI_t *eri_lin, *eri_log;
  ERI_Orbital_t *orb1, *orb2, *orb3, *orb4;

  /* parameters */
  int lmax, lmax_gl, ngrid, ngl;
  double rmax;
  ERI_Init_Misc info;


  printf("----------------------------------------------------------------\n");
  printf("  TEST #1: GTO\n");
  printf("----------------------------------------------------------------\n");
  printf("\n");

  {
    printf("INITIALIZING LIBERI WITH LINEAR-SBT:\n");
  
    lmax    = 15;
    lmax_gl = 8;
    ngrid   = 1024;
    ngl     = 100;
   
    info.sbttype = ERI_SBT_LINEAR;
    info.rmax = 80.0; 
 
    printf("  LMAX    = %5d\n"   , lmax   );
    printf("  LMAX_GL = %5d\n"   , lmax_gl);
    printf("  NGRID   = %5d\n"   , ngrid  );
    printf("  NGL     = %5d\n"   , ngl    );
    printf("  RMAX    = %12.4f\n", info.rmax   );
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
    lmax_gl = 8;
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
    printf("  RMAX    = %12.4f\n", info.rmax   );
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

  
  {
    printf("BASIS FUNCTIONS (ssss):\n");
    printf("  (ssss) 1 : GTO  L=0  M=0   EXPONENT=1.0  NORMALIZATION=1.0\n");
    printf("         2 : GTO  L=0  M=0   EXPONENT=1.0  NORMALIZATION=1.0\n");
    printf("         3 : GTO  L=0  M=0   EXPONENT=1.0  NORMALIZATION=1.0\n");
    printf("         4 : GTO  L=0  M=0   EXPONENT=1.0  NORMALIZATION=1.0\n");
    printf("  (pssp) 1 : GTO  L=1  M=-1  EXPONENT=1.0  NORMALIZATION=1.0\n");
    printf("         2 : GTO  L=0  M=0   EXPONENT=1.0  NORMALIZATION=1.0\n");
    printf("         3 : GTO  L=0  M=0   EXPONENT=1.0  NORMALIZATION=1.0\n");
    printf("         4 : GTO  L=1  M=+1  EXPONENT=1.0  NORMALIZATION=1.0\n");
    printf("  (dddd) 1 : GTO  L=2  M=-2  EXPONENT=0.5  NORMALIZATION=1.0\n");
    printf("         2 : GTO  L=2  M=-1  EXPONENT=0.5  NORMALIZATION=1.0\n");
    printf("         3 : GTO  L=2  M=+1  EXPONENT=0.5  NORMALIZATION=1.0\n");
    printf("         4 : GTO  L=2  M=+2  EXPONENT=0.5  NORMALIZATION=1.0\n");
    printf("\n");

    printf("CENTERS:\n");
    printf("  1 : (%12.8f, %12.8f, %12.8f)\n",  0.5,  0.5, 0.0);
    printf("  2 : (%12.8f, %12.8f, %12.8f)\n",  0.5, -0.5, 0.0);
    printf("  3 : (%12.8f, %12.8f, %12.8f)\n", -0.5,  0.5, 0.0);
    printf("  4 : (%12.8f, %12.8f, %12.8f)\n", -0.5, -0.5, 0.0);
    printf("\n");
  }

  sec_sum[0] = 0.0;
  sec_sum[1] = 0.0;

  {
    printf("COMBINATION (ssss):\n");
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


    {
      clk0 = clock(); 
      ERI_Integral(eri_lin, I4_ssss_lin, NULL, orb1, orb2, orb3, orb4, ERI_NOSCREEN);
      clk1 = clock(); 
      sec = (double)(clk1 - clk0)/(double)CLOCKS_PER_SEC;

      printf("  LINEAR-SBT:\n");
      printf("    INT  = %20.15e + I %20.15e\n", I4_ssss_lin[0], I4_ssss_lin[1]);  
      printf("    TIME = %12.8f SEC\n", sec);
      
      sec_sum[0] += sec;
    }

    {
      clk0 = clock();
      ERI_Integral(eri_log, I4_ssss_log, NULL, orb1, orb2, orb3, orb4, ERI_NOSCREEN);
      clk1 = clock();
      sec = (double)(clk1 - clk0)/(double)CLOCKS_PER_SEC;

      printf("  LOG-SBT:\n");
      printf("    INT  = %20.15e + I %20.15e\n", I4_ssss_log[0], I4_ssss_log[1]);  
      printf("    TIME = %12.8f SEC\n", sec);
      printf("\n");
      
      sec_sum[1] += sec;
    }
  
    Orbital_Free(orb1);
    Orbital_Free(orb2);
    Orbital_Free(orb3);
    Orbital_Free(orb4);
  }
  

  {
    printf("COMBINATION (pssp):\n");
    ngrid = 8192;
    rmax  = 100.0;
    orb1 = Orbital_New_GTO(rmax, ngrid, 1, -1, 1.0, 1.0);
    orb2 = Orbital_New_GTO(rmax, ngrid, 0,  0, 1.0, 1.0);    
    orb3 = Orbital_New_GTO(rmax, ngrid, 0,  0, 1.0, 1.0);    
    orb4 = Orbital_New_GTO(rmax, ngrid, 1,  1, 1.0, 1.0);
    if (NULL == orb1 || NULL == orb2 || NULL == orb3 || NULL == orb4) {
      fprintf(stderr, "*** ERROR in %s (%d)\n", __FILE__, __LINE__);
      abort();
    }
 
    orb1->c[0] =  0.5; orb1->c[1] =  0.5; orb1->c[2] = 0.0;
    orb2->c[0] =  0.5; orb2->c[1] = -0.5; orb2->c[2] = 0.0;
    orb3->c[0] = -0.5; orb3->c[1] =  0.5; orb3->c[2] = 0.0;
    orb4->c[0] = -0.5; orb4->c[1] = -0.5; orb4->c[2] = 0.0;
 
    {
      clk0 = clock(); 
      ERI_Integral(eri_lin, I4_pssp_lin, NULL, orb1, orb2, orb3, orb4, ERI_NOSCREEN);
      clk1 = clock(); 
      sec = (double)(clk1 - clk0)/(double)CLOCKS_PER_SEC;

      printf("  LINEAR-SBT:\n");
      printf("    INT  = %20.15e + I %20.15e\n", I4_pssp_lin[0], I4_pssp_lin[1]);  
      printf("    TIME = %12.8f SEC\n", sec);
      
      sec_sum[0] += sec;
    }

    { 
      clk0 = clock();
      ERI_Integral(eri_log, I4_pssp_log, NULL, orb1, orb2, orb3, orb4, ERI_NOSCREEN);
      clk1 = clock();
      sec = (double)(clk1 - clk0)/(double)CLOCKS_PER_SEC;

      printf("  LOG-SBT:\n");
      printf("    INT  = %20.15e + I %20.15e\n", I4_pssp_log[0], I4_pssp_log[1]);  
      printf("    TIME = %12.8f SEC\n", sec);
      printf("\n");
      
      sec_sum[1] += sec;
    } 
    Orbital_Free(orb1);
    Orbital_Free(orb2);
    Orbital_Free(orb3);
    Orbital_Free(orb4);
  }


  {
    printf("COMBINATION (dddd):\n");
    ngrid = 8192;
    rmax  = 100.0;
    orb1 = Orbital_New_GTO(rmax, ngrid, 2, -2, 0.5, 1.0);
    orb2 = Orbital_New_GTO(rmax, ngrid, 2, -1, 0.5, 1.0);    
    orb3 = Orbital_New_GTO(rmax, ngrid, 2,  1, 0.5, 1.0);    
    orb4 = Orbital_New_GTO(rmax, ngrid, 2,  2, 0.5, 1.0);
    if (NULL == orb1 || NULL == orb2 || NULL == orb3 || NULL == orb4) {
      fprintf(stderr, "*** ERROR in %s (%d)\n", __FILE__, __LINE__);
      abort();
    }
 
    orb1->c[0] =  0.5; orb1->c[1] =  0.5; orb1->c[2] = 0.0;
    orb2->c[0] =  0.5; orb2->c[1] = -0.5; orb2->c[2] = 0.0;
    orb3->c[0] = -0.5; orb3->c[1] =  0.5; orb3->c[2] = 0.0;
    orb4->c[0] = -0.5; orb4->c[1] = -0.5; orb4->c[2] = 0.0;

    {
      clk0 = clock(); 
      ERI_Integral(eri_lin, I4_dddd_lin, NULL, orb1, orb2, orb3, orb4, ERI_NOSCREEN);
      clk1 = clock(); 
      sec = (double)(clk1 - clk0)/(double)CLOCKS_PER_SEC;

      printf("  LINEAR-SBT:\n");
      printf("    INT  = %20.15e + I %20.15e\n", I4_dddd_lin[0], I4_dddd_lin[1]);  
      printf("    TIME = %12.8f SEC\n", sec);
      
      sec_sum[0] += sec;
    }

    {
      clk0 = clock();
      ERI_Integral(eri_log, I4_dddd_log, NULL, orb1, orb2, orb3, orb4, ERI_NOSCREEN);
      clk1 = clock();
      sec = (double)(clk1 - clk0)/(double)CLOCKS_PER_SEC;

      printf("  LOG-SBT:\n");
      printf("    INT  = %20.15e + I %20.15e\n", I4_dddd_log[0], I4_dddd_log[1]);  
      printf("    TIME = %12.8f SEC\n", sec);
      printf("\n");
      
      sec_sum[1] += sec;
    }
 
    Orbital_Free(orb1);
    Orbital_Free(orb2);
    Orbital_Free(orb3);
    Orbital_Free(orb4);
  }
  printf("\n\n");

  ERI_Free(eri_lin);
  ERI_Free(eri_log);

  I4[0][0] = I4_ssss_lin[0];
  I4[0][1] = I4_ssss_log[0];
  I4[0][2] = 0.00760884651642564;

  I4[1][0] = I4_pssp_lin[0];
  I4[1][1] = I4_pssp_log[0];
  I4[1][2] = 0.00251250623063648;

  I4[2][0] = I4_dddd_lin[0];
  I4[2][1] = I4_dddd_log[0];
  I4[2][2] = -0.00196689855434091;

  out_sec[0] = sec_sum[0];
  out_sec[1] = sec_sum[1];
}


/* EOF */
