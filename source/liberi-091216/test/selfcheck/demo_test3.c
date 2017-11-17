/*----------------------------------------------------------------------
  demo_test3.c
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "eri.h"
#include "demo_sub.h"


static void demo_dirivative(
  ERI_t *eri, 
  const ERI_Orbital_t *in_o1,
  const ERI_Orbital_t *in_o2,
  const ERI_Orbital_t *in_o3,
  const ERI_Orbital_t *in_o4,
  double *out_err
)
{
  int i, j, k, ixyz;

  ERI_Orbital_t* qt[4];
  double I4[2], dI4[4][3][2], dif[12][2], err[12][2];
  double tmp, maxerr;

  const char *card[12] = {
    "dI4/dx1", "dI4/dy1", "dI4/dz1", "dI4/dx2", 
    "dI4/dy2", "dI4/dz2", "dI4/dx3", "dI4/dy3", 
    "dI4/dz3", "dI4/dx4", "dI4/dy4", "dI4/dz4" 
  };
  const char xyz[3] = "xyz";

  double delta = 1e-3;

  qt[0] = Orbital_Duplicate(in_o1);
  qt[1] = Orbital_Duplicate(in_o2);
  qt[2] = Orbital_Duplicate(in_o3);
  qt[3] = Orbital_Duplicate(in_o4);
  
  ERI_Integral(eri, I4, dI4, qt[0], qt[1], qt[2], qt[3], ERI_NOSCREEN);
  printf("  I4 = %12.8e  %12.8e\n", I4[0], I4[1]);
  printf("\n");
  
  /* differential by finite difference delta */
  for (i=0; i<4; i++) {
    for (ixyz=0; ixyz<3; ixyz++) {
      j = i*3+ixyz;
      tmp = qt[i]->c[ixyz]; 
      qt[i]->c[ixyz] += delta;
      ERI_Integral(eri, dif[j], NULL, qt[0], qt[1], qt[2], qt[3], ERI_NOSCREEN);
      printf("  %c%1d = ", (char)((int)('x')+ixyz), i);
      printf("%12.8e  %12.8e\n", dif[j][0], dif[j][1]);
      dif[j][0] = (dif[j][0] - I4[0])/delta;
      dif[j][1] = (dif[j][1] - I4[1])/delta;
      err[j][0]   = fabs(dif[j][0] - dI4[i][ixyz][0]);
      err[j][1]   = fabs(dif[j][1] - dI4[i][ixyz][1]);
      qt[i]->c[ixyz] = tmp;
    }
  } 
  printf("\n");
  
  printf("DERIVATIVES:\n");
  for (i=0; i<12; i++) {
    j = i%3;
    k = i/3;
    printf("  %s = %12.8e  %12.8e\n", card[i], dI4[k][j][0], dI4[k][j][1]);
  }
  printf("\n");
  
  printf("FINITE DIFFERENCE (DELTA=%12.8e):\n", delta);
  for (i=0; i<12; i++) {
    printf("  %s = %12.8e  %12.8e\n", card[i], dif[i][0], dif[i][1]);
  }
  printf("\n");

  printf("ERROR:\n");
  maxerr = 0.0;
  for (i=0; i<12; i++) {
    printf("  %s = %12.8e  %12.8e\n", card[i], err[i][0], err[i][1]);
    if (maxerr<err[i][0]) maxerr = err[i][0];
    if (maxerr<err[i][1]) maxerr = err[i][1];
  }
  printf("\n");

  *out_err = maxerr;
}




void demo_test3(double out_err[2], double out_sec[2])
{
  ERI_t *eri_lin, *eri_log;
  size_t mem;
  clock_t clk0, clk1;

  /* parameters */
  int lmax, lmax_gl, ngrid, ngl;
  double rmax;
  ERI_Init_Misc info;

  /* basis functions */
  ERI_Orbital_t *orb1, *orb2, *orb3, *orb4;
  const char* pao_file[4] = {"./pao/H5.0.pao", "./pao/C4.0.pao",
                             "./pao/O5.0.pao", "./pao/Fe5.5.pao" };
  int         pao_l[4]   = { 0, 1, 2, 3 };
  int         pao_m[4]   = { 0, -1, 0, 2 };
  int         pao_n[4]   = { 0, 1, 2, 3 };

  printf("----------------------------------------------------------------\n");
  printf("  TEST #3: DERIVATIVES\n");
  printf("----------------------------------------------------------------\n");
  printf("\n");

  {
    printf("INITIALIZING LIBERI WITH LINEAR-SBT:\n");
  
    lmax    = 15;
    lmax_gl = 8;
    ngrid   = 2048;
    ngl     = 100;

    info.sbttype = ERI_SBT_LINEAR;
    info.rmax    = 100.0;
 
    printf("  LMAX    = %5d\n"   , lmax   );
    printf("  LMAX_GL = %5d\n"   , lmax_gl);
    printf("  NGRID   = %5d\n"   , ngrid  );
    printf("  NGL     = %5d\n"   , ngl    );
    printf("  RMAX    = %12.4f\n", info.rmax );
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
    ngrid   = 2048;
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

  /* basis functions */
  {
    orb1 = Orbital_New_PAO(pao_file[0], pao_l[0], pao_m[0], pao_n[0]);
    orb2 = Orbital_New_PAO(pao_file[1], pao_l[1], pao_m[1], pao_n[1]);
    orb3 = Orbital_New_PAO(pao_file[2], pao_l[2], pao_m[2], pao_n[2]);
    orb4 = Orbital_New_PAO(pao_file[3], pao_l[3], pao_m[3], pao_n[3]);
    if (NULL==orb1 || NULL==orb2 || NULL==orb3 || NULL==orb4) {
      fprintf(stderr, "***ERROR in %s (%d)\n", __FILE__, __LINE__);
      abort();
    }

    orb1->c[0] = 0.0000; orb1->c[1] = 0.0000; orb1->c[2] = 0.0000;
    orb2->c[0] = 0.0000; orb2->c[1] = 0.5000; orb2->c[2] = 0.0000;
    orb3->c[0] = 0.0000; orb3->c[1] = 0.0000; orb3->c[2] = 1.5000;
    orb4->c[0] = 1.1549; orb4->c[1] = 0.7500; orb4->c[2] = 2.0865;
  
    printf("BASIS FUNCTIONS:\n");
    printf("  1 : FILE=%10s   NODE=%d  L=%d  M=%d\n", pao_file[0], pao_n[0], pao_l[0], pao_m[0]);
    printf("  2 : FILE=%10s   NODE=%d  L=%d  M=%d\n", pao_file[1], pao_n[1], pao_l[1], pao_m[1]);
    printf("  3 : FILE=%10s   NODE=%d  L=%d  M=%d\n", pao_file[2], pao_n[2], pao_l[2], pao_m[2]);
    printf("  4 : FILE=%10s   NODE=%d  L=%d  M=%d\n", pao_file[3], pao_n[3], pao_l[3], pao_m[3]);
    printf("\n");

    printf("CENTERS:\n");
    printf("  1 : (%12.8f, %12.8f, %12.8f)\n", orb1->c[0], orb1->c[1], orb1->c[2]);
    printf("  2 : (%12.8f, %12.8f, %12.8f)\n", orb2->c[0], orb2->c[1], orb2->c[2]);
    printf("  3 : (%12.8f, %12.8f, %12.8f)\n", orb3->c[0], orb3->c[1], orb3->c[2]);
    printf("  4 : (%12.8f, %12.8f, %12.8f)\n", orb4->c[0], orb4->c[1], orb4->c[2]);
    printf("\n");
  }

  clk0 = clock();
  demo_dirivative(eri_lin, orb1, orb2, orb3, orb4, &out_err[0]);
  clk1 = clock();
  out_sec[0] = (double)(clk1-clk0)/(double)CLOCKS_PER_SEC;

  clk0 = clock();
  demo_dirivative(eri_log, orb1, orb2, orb3, orb4, &out_err[1]);
  clk1 = clock();
  out_sec[1] = (double)(clk1-clk0)/(double)CLOCKS_PER_SEC;
}


/* EOF */
