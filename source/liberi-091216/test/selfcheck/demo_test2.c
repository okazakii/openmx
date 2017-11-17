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
  double *I4,
  double *sec
)
{
  clock_t clk0, clk1;

  clk0 = clock(); 
  ERI_Integral(eri, I4, NULL, orb1, orb2, orb3, orb4, ERI_NOSCREEN);
  clk1 = clock(); 
  *sec = (double)(clk1 - clk0)/(double)CLOCKS_PER_SEC;

  printf("    INT  = %20.15e + I %20.15e\n", I4[0], I4[1]);  
  printf("    TIME = %12.8f SEC\n", sec);
}


void demo_test2(double I4[5][3], double out_sec[2])
{
  int i;
  double I4_1234_lin[2], I4_1234_log[2];
  double I4_1342_lin[2], I4_1342_log[2];
  double I4_1423_lin[2], I4_1423_log[2];
  double I4_2134_lin[2], I4_2134_log[2];
  double I4_3124_lin[2], I4_3124_log[2];

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
  printf("  TEST #2: PAO\n");
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
    ngrid = 8192;
    rmax  = 1000.0; 
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


  printf("COMBINATION (1234):\n");
  printf("  LINEAR-SBT:\n");
  demo(eri_lin, orb1, orb2, orb3, orb4, I4_1234_lin, &sec);
  sec_sum[0] += sec;
  printf("  LOG-SBT:\n");
  demo(eri_log, orb1, orb2, orb3, orb4, I4_1234_log, &sec);
  sec_sum[1] += sec;
  printf("\n");
    
  printf("COMBINATION (1342):\n");
  printf("  LINEAR-SBT:\n");
  demo(eri_lin, orb1, orb3, orb4, orb2, I4_1342_lin, &sec);
  sec_sum[0] += sec;
  printf("  LOG-SBT:\n");
  demo(eri_log, orb1, orb3, orb4, orb2, I4_1342_log, &sec);
  sec_sum[1] += sec;
  printf("\n");
 
  printf("COMBINATION (1423):\n");
  printf("  LINEAR-SBT:\n");
  demo(eri_lin, orb1, orb4, orb2, orb3, I4_1423_lin, &sec);
  sec_sum[0] += sec;
  printf("  LOG-SBT:\n");
  demo(eri_log, orb1, orb4, orb2, orb3, I4_1423_log, &sec);
  sec_sum[1] += sec;
  printf("\n");

  printf("COMBINATION (2134):\n");
  printf("  LINEAR-SBT:\n");
  demo(eri_lin, orb2, orb1, orb3, orb4, I4_2134_lin, &sec);
  sec_sum[0] += sec;
  printf("  LOG-SBT:\n");
  demo(eri_log, orb2, orb1, orb3, orb4, I4_2134_log, &sec);
  sec_sum[1] += sec;
  printf("\n");
    
  printf("COMBINATION (3124):\n");
  printf("  LINEAR-SBT:\n");
  demo(eri_lin, orb3, orb1, orb2, orb4, I4_3124_lin, &sec);
  sec_sum[0] += sec;
  printf("  LOG-SBT:\n");
  demo(eri_log, orb3, orb1, orb2, orb4, I4_3124_log, &sec);
  sec_sum[1] += sec;
  printf("\n\n");

  Orbital_Free(orb1);
  Orbital_Free(orb2);
  Orbital_Free(orb3);
  Orbital_Free(orb4);

  ERI_Free(eri_lin);
  ERI_Free(eri_log);
  
  I4[0][0] = I4_1234_lin[0];
  I4[0][1] = I4_1234_log[0];
  I4[0][2] = -0.00179763; /*  0.00051067 */
  
  I4[1][0] = I4_1342_lin[0];
  I4[1][1] = I4_1342_log[0];
  I4[1][2] = -0.00351254; /* -0.00039812 */

  I4[2][0] = I4_1423_lin[0];
  I4[2][1] = I4_1423_log[0];
  I4[2][2] = -0.00238960; /*  0.00045194 */
  
  I4[3][0] = I4_2134_lin[0];
  I4[3][1] = I4_2134_log[0];
  I4[3][2] = -0.00179763; /*  0.00051067 */
  
  I4[4][0] = I4_3124_lin[0];
  I4[4][1] = I4_3124_log[0];
  I4[4][2] = -0.00351254; /* -0.00039812 */

  out_sec[0] = sec_sum[0];
  out_sec[1] = sec_sum[1];
}


/* EOF */
