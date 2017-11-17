/*----------------------------------------------------------------------
  demo.c

  Sample program of LIBERI

  Coded by TOYODA MAsayuki, June 2009.
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "eri.h"

#define NMESH 8192


static void integral(
  ERI_t         *solver,
  ERI_Orbital_t p[4],
  int           nitr,
  double        I4[3][2], /* (OUT) integral */
  double        *sec
)
{
  int lmax_gl, itr;
  double *glF, *glG;
  double cx12, cx34, R[3], a1[3], a2[3], a3[3], a4[3];
  double c1[3], c2[3], c3[3], c4[3];
  double Rd[3], scr;
  clock_t t1, t2;
 
  lmax_gl = ERI_lmax_gl(solver);

  /* allocate workspace memory */
  glF = (double*)malloc( ERI_Size_of_GLF(solver) );
  glG = (double*)malloc( ERI_Size_of_GLF(solver) );

  c1[0] = p[0].c[0]; c1[1] = p[0].c[1]; c1[2] = p[0].c[2];
  c2[0] = p[1].c[0]; c2[1] = p[1].c[1]; c2[2] = p[1].c[2];
  c3[0] = p[2].c[0]; c3[1] = p[2].c[1]; c3[2] = p[2].c[2];
  c4[0] = p[3].c[0]; c4[1] = p[3].c[1]; c4[2] = p[3].c[2];

  /* expansion center */
  cx12 = ERI_Center_r2(&p[0], &p[1]);
  cx34 = ERI_Center_r2(&p[2], &p[3]);

  scr = ERI_NOSCREEN;

  ERI_Overlap(solver, glF, NULL, &p[0], &p[1], cx12); 
  ERI_Overlap(solver, glG, NULL, &p[2], &p[3], cx34);

  /* R */
  ERI_Coordinate_Transform(R, a1, a2, a3, a4, c1, c2, c3, c4, cx12, cx34);

  t1 = clock();
  for (itr=0; itr<nitr; itr++) {
    Rd[0] = R[0];
    Rd[1] = R[1];
    Rd[2] = R[2];
    ERI_Integral_GL(solver, I4[0], glF, glG, 
      Rd[0], Rd[1], Rd[2], scr, lmax_gl);
  
    Rd[0] = R[0]+1.0;
    Rd[1] = R[1];
    Rd[2] = R[2];
    ERI_Integral_GL(solver, I4[1], glF, glG, 
      Rd[0], Rd[1], Rd[2], scr, lmax_gl);
  
    Rd[0] = R[0];
    Rd[1] = R[1]+0.4;
    Rd[2] = R[2];
    ERI_Integral_GL(solver, I4[2], glF, glG, 
      Rd[0], Rd[1], Rd[2], scr, lmax_gl);
  }
  t2 = clock();

  *sec = (double)(t2-t1)/(double)CLOCKS_PER_SEC;

  /* release workspace memory */
  free(glF);
  free(glG);
}




static void integral_pp(
  ERI_t         *solver,
  ERI_Orbital_t p[4],
  int           nitr,
  double        I4[3][2], 
  double        *sec
)
{
  int lmax_gl, ngl, itr, jmax_gl;
  double *glF, *glG;
  double cx12, cx34, R[3], a1[3], a2[3], a3[3], a4[3];
  double c1[3], c2[3], c3[3], c4[3], scr, Rd[3], Rd_t[3], Rd_p[3];
  double tmp_I4[6];
  clock_t t1, t2;
  
  int numR; 
  double *mul_gc, *prej, *preY;
  int *mul_j2, mul_n;

  int num_minimalR, *minimalR;

  int i,j;
 
  lmax_gl = ERI_lmax_gl(solver);
  jmax_gl = lmax_gl*lmax_gl;
  ngl     = ERI_ngl(solver);
  scr     = ERI_NOSCREEN;
  numR    = 3;

  /* allocate workspace memory */
  glF = (double*)malloc( ERI_Size_of_GLF(solver) );
  glG = (double*)malloc( ERI_Size_of_GLF(solver) );
  mul_j2 = (int*)malloc(sizeof(int)*jmax_gl*jmax_gl*jmax_gl);
  mul_gc = (double*)malloc(sizeof(double)*jmax_gl*jmax_gl*jmax_gl);
  mul_n  = (int*)malloc(sizeof(double)*jmax_gl*jmax_gl);
  prej = (double*)malloc(sizeof(double)*ngl*lmax_gl*numR); 
  preY = (double*)malloc(sizeof(double)*numR*jmax_gl);
  minimalR= (int*)malloc(sizeof(int)*numR);

  c1[0] = p[0].c[0]; c1[1] = p[0].c[1]; c1[2] = p[0].c[2];
  c2[0] = p[1].c[0]; c2[1] = p[1].c[1]; c2[2] = p[1].c[2];
  c3[0] = p[2].c[0]; c3[1] = p[2].c[1]; c3[2] = p[2].c[2];
  c4[0] = p[3].c[0]; c4[1] = p[3].c[1]; c4[2] = p[3].c[2];

  /* expansion center */
  cx12 = ERI_Center_r2(&p[0], &p[1]);
  cx34 = ERI_Center_r2(&p[2], &p[3]);

  /* R */
  ERI_Coordinate_Transform(R, a1, a2, a3, a4, c1, c2, c3, c4, cx12, cx34);

  ERI_Overlap(solver, glF, NULL, &p[0], &p[1], cx12); 
  ERI_Overlap(solver, glG, NULL, &p[2], &p[3], cx34);

  Rd  [0] = R[0];
  Rd_t[0] = R[1];
  Rd_p[0] = R[2];
  
  Rd  [1] = R[0]+1.0;
  Rd_t[1] = R[1];
  Rd_p[1] = R[2];
   
  Rd  [2] = R[0];
  Rd_t[2] = R[1]+0.4;
  Rd_p[2] = R[2];

  ERI_Integral_GL_PrejY(solver, Rd, Rd_t, Rd_p, numR,
    scr, prej, preY, mul_j2, mul_gc, mul_n, minimalR, &num_minimalR);

  t1 = clock();
  for (itr=0; itr<nitr; itr++) {
    ERI_Integral_GL_Post(solver, tmp_I4, glF, glG, numR,
      prej, preY, mul_j2, mul_gc, mul_n, minimalR, num_minimalR);
  }
  t2 = clock();

  *sec = (double)(t2-t1)/(double)CLOCKS_PER_SEC;

  I4[0][0] = tmp_I4[0];
  I4[0][1] = 0.0;
  I4[1][0] = tmp_I4[1];
  I4[1][1] = 0.0;
  I4[2][0] = tmp_I4[2];
  I4[2][1] = 0.0;

  /* release workspace memory */
  free(glF);
  free(glG);

  free(mul_j2);
  free(mul_gc);
  free(mul_n);
  free(prej);
  free(preY);
  free(minimalR);
}




int main(void)
{
  int i;
  double r, I4[3][2];
  const double rmax = 300.0;
  double s1, s2, s3;

  ERI_t *eri;
  ERI_Orbital_t p[4];

  double xr[NMESH], gto_s[NMESH];

  /* radial mesh and GTO */
  for (i=0; i<NMESH; i++) { 
    r = rmax*(double)i/(double)NMESH; 
    xr[i] = r;
    gto_s[i] = exp(-r*r); 
  } 

  /*  prepare orbital informations */
  p[0].fr = gto_s; p[0].xr = xr; p[0].ngrid = NMESH; p[0].l = 1; p[0].m = -1;
  p[1].fr = gto_s; p[1].xr = xr; p[1].ngrid = NMESH; p[1].l = 0; p[1].m = 0;
  p[2].fr = gto_s; p[2].xr = xr; p[2].ngrid = NMESH; p[2].l = 0; p[2].m = 0;
  p[3].fr = gto_s; p[3].xr = xr; p[3].ngrid = NMESH; p[3].l = 0; p[3].m = 0;

  /* positions */
  p[0].c[0] =  0.2; p[0].c[1] =  0.2; p[0].c[2] = 0.0;
  p[1].c[0] =  0.2; p[1].c[1] = -0.2; p[1].c[2] = 0.0;
  p[2].c[0] = -0.2; p[2].c[1] =  0.3; p[2].c[2] = 0.0;
  p[3].c[0] = -0.2; p[3].c[1] = -0.3; p[3].c[2] = 0.0;

  /* initialize LIBERI */
  //eri = ERI_Init(15, 8, 1024, 100, ERI_SH_COMPLEX, NULL);
  eri = ERI_Init(15, 8, 1024, 100, ERI_SH_REAL, NULL);
  if (NULL == eri) {
    fprintf(stderr, "*** Failed to initialize LIBERI.\n");
    return -1;
  }

  /* calculate ERI */
  integral(eri, p, 100, I4, &s1);
  printf("  %10.6f  %10.6f\n", I4[0][0], I4[0][1]);
  printf("  %10.6f  %10.6f\n", I4[1][0], I4[1][1]);
  printf("  %10.6f  %10.6f\n", I4[2][0], I4[2][1]);
  printf("  TIME = %10.6f\n", s1);
  printf("\n");

  integral_pp(eri, p, 100, I4, &s2);
  printf("  %10.6f  %10.6f\n", I4[0][0], I4[0][1]);
  printf("  %10.6f  %10.6f\n", I4[1][0], I4[1][1]);
  printf("  %10.6f  %10.6f\n", I4[2][0], I4[2][1]);
  printf("  TIME = %10.6f\n", s2);
  printf("\n");

  /* release LIBERI */
  ERI_Free(eri);

  return 0;
}


/* EOF */
