/*----------------------------------------------------------------------
  eri.c

  High-level interface of LIBERI 

  Coded by TOYODA Masayuki, 19 Jun. 2009.
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include "eri_def.h"
#include "eri.h"




static void coord_spherical(
  double  x,      /* (IN) Cartesian coordinate */ 
  double  y, 
  double  z,
  double *r,      /* (OUT) spherical coordinate */
  double *theta, 
  double *phi
)
{
  *r     = sqrt(x*x + y*y + z*z);
  *theta = atan2(sqrt(x*x+y*y),z);
  *phi   = atan2(y,x);
}




void ERI_Overlap(
  ERI_t               *solver,  /* (INOUT) ERI object */
  double              *glF,     /* (OUT) F term */
  double              *dglF[3], /* (OUT) derivatives of F */
  const ERI_Orbital_t *orb1,    /* (IN) orbital 1 */
  const ERI_Orbital_t *orb2,    /* (IN) orbital 2 */
  double               cx       /* (IN) expansion center */
)
{
  int ngrid, lmax, jmax, j;
  
  const double *fr1, *fr2, *xr1, *xr2;
  int ngrid1, ngrid2, l1, l2, m1, m2;
  double c[3], c1[3], c2[3];
  double ar1, at1, ap1, ar2, at2, ap2;
  
  double *fk1, *fk2, *gam1, *gam2, *alp1, *alp2, *p, *F;
  double *dgam1, *dgam2, *dalp1[3], *dalp2[3], *dp[3], *dF[3];

  ngrid = ERI_ngrid(solver);
  lmax  = ERI_lmax(solver);
  jmax  = lmax*lmax;
 
  /* WF 1 */ 
  fr1    = orb1->fr;
  xr1    = orb1->xr;
  ngrid1 = orb1->ngrid;
  l1     = orb1->l;
  m1     = orb1->m;
  c1[0]  = orb1->c[0]; 
  c1[1]  = orb1->c[1]; 
  c1[2]  = orb1->c[2];
  
  /* WF 2 */
  fr2    = orb2->fr;
  xr2    = orb2->xr;
  ngrid2 = orb2->ngrid;
  l2     = orb2->l;
  m2     = orb2->m;
  c2[0]  = orb2->c[0]; 
  c2[1]  = orb2->c[1]; 
  c2[2]  = orb2->c[2];

  STEPTRACE("ERI_Overlap: in\n");

  /* expansion centers */
  c[0] = cx*c1[0] + (1.0-cx)*c2[0]; 
  c[1] = cx*c1[1] + (1.0-cx)*c2[1]; 
  c[2] = cx*c1[2] + (1.0-cx)*c2[2]; 
 
  coord_spherical(c1[0]-c[0], c1[1]-c[1], c1[2]-c[2], &ar1, &at1, &ap1);
  coord_spherical(c2[0]-c[0], c2[1]-c[1], c2[2]-c[2], &ar2, &at2, &ap2);

  /* allocate workspace memory */
  fk1  = (double*)malloc( ERI_Size_of_Orbital(solver) );
  fk2  = (double*)malloc( ERI_Size_of_Orbital(solver) );
  gam1 = (double*)malloc( ERI_Size_of_Gamma(solver)   );
  gam2 = (double*)malloc( ERI_Size_of_Gamma(solver)   );
  alp1 = (double*)malloc( ERI_Size_of_Alpha(solver)   );
  alp2 = (double*)malloc( ERI_Size_of_Alpha(solver)   );
  p    = (double*)malloc( ERI_Size_of_Overlap(solver) );
  F    = (double*)malloc( ERI_Size_of_Overlap(solver) );

  if (dglF != NULL) {
    dgam1    = (double*)malloc( ERI_Size_of_Gamma(solver)    );
    dgam2    = (double*)malloc( ERI_Size_of_Gamma(solver)    );
    dalp1[0] = (double*)malloc( ERI_Size_of_Alpha(solver)    );
    dalp1[1] = (double*)malloc( ERI_Size_of_Alpha(solver)    );
    dalp1[2] = (double*)malloc( ERI_Size_of_Alpha(solver)    );
    dalp2[0] = (double*)malloc( ERI_Size_of_Alpha(solver)    );
    dalp2[1] = (double*)malloc( ERI_Size_of_Alpha(solver)    );
    dalp2[2] = (double*)malloc( ERI_Size_of_Alpha(solver)    );
    dp[0]    = (double*)malloc( ERI_Size_of_Overlap(solver) );
    dp[1]    = (double*)malloc( ERI_Size_of_Overlap(solver) );
    dp[2]    = (double*)malloc( ERI_Size_of_Overlap(solver) );
    dF[0]    = (double*)malloc( ERI_Size_of_Overlap(solver) );
    dF[1]    = (double*)malloc( ERI_Size_of_Overlap(solver) );
    dF[2]    = (double*)malloc( ERI_Size_of_Overlap(solver) );
  }
  
  ERI_Transform_Orbital(solver, fk1, fr1, xr1, ngrid1, l1);
  ERI_Transform_Orbital(solver, fk2, fr2, xr2, ngrid2, l2);
   
  if (NULL == dglF) {
    /* without derivatives */
    ERI_LL_Gamma(solver, gam1, NULL, fk1, ar1); 
    ERI_LL_Gamma(solver, gam2, NULL, fk2, ar2); 
    ERI_LL_Alpha(solver, alp1, gam1, at1, ap1, l1, m1);
    ERI_LL_Alpha(solver, alp2, gam2, at2, ap2, l2, m2);
    ERI_LL_Overlap(solver, p, alp1, alp2);
    ERI_Transform_Overlap(solver, F, p );
    ERI_GL_Interpolate(solver, glF, F);
  } else {
    ERI_LL_Gamma(solver, gam1, dgam1, fk1, ar1); 
    ERI_LL_Gamma(solver, gam2, dgam2, fk2, ar2); 
    ERI_LL_Alpha_d(solver, alp1, dalp1, gam1, dgam1, ar1, at1, ap1, l1, m1);
    ERI_LL_Alpha_d(solver, alp2, dalp2, gam2, dgam2, ar2, at2, ap2, l2, m2);
    ERI_LL_Overlap_d(solver, p, dp, alp1, dalp1, alp2, dalp2, cx);
    ERI_Transform_Overlap(solver, F, p );
    ERI_Transform_Overlap(solver, dF[0], dp[0] );
    ERI_Transform_Overlap(solver, dF[1], dp[1] );
    ERI_Transform_Overlap(solver, dF[2], dp[2] );
    ERI_GL_Interpolate(solver, glF, F);
    ERI_GL_Interpolate(solver, dglF[0], dF[0]);
    ERI_GL_Interpolate(solver, dglF[1], dF[1]);
    ERI_GL_Interpolate(solver, dglF[2], dF[2]);
  }

  /* release workspace memory */
  free(fk1);
  free(fk2);
  free(gam1);
  free(gam2);
  free(alp1);
  free(alp2);
  free(p);
  free(F);

  if (NULL != dglF) {
    free(dgam1); 
    free(dgam2); 
    free(dalp1[0]); 
    free(dalp1[1]); 
    free(dalp1[2]); 
    free(dalp2[0]); 
    free(dalp2[1]); 
    free(dalp2[2]); 
    free(dp[0]); 
    free(dp[1]); 
    free(dp[2]); 
    free(dF[0]); 
    free(dF[1]); 
    free(dF[2]); 
  }
  
  STEPTRACE("ERI_Overlap: out\n");
}



void ERI_Integral(
  ERI_t               *solver,
  double               I4[2], /* (OUT) integral */
  double              *dI4,   /* (OUT) derivatives */
  const ERI_Orbital_t *orb1, 
  const ERI_Orbital_t *orb2, 
  const ERI_Orbital_t *orb3, 
  const ERI_Orbital_t *orb4,
  double               scr
)
{
  int lmax_gl;
  double *glF, *glG, *dglF[3], *dglG[3];
  double cx12, cx34, R[3], a1[3], a2[3], a3[3], a4[3];
  double c1[3], c2[3], c3[3], c4[3];
  clock_t clk1, clk2, clk3;
 
  lmax_gl = ERI_lmax_gl(solver);

  STEPTRACE("ERI_Integral: in");

  /* allocate workspace memory */
  glF = (double*)malloc( ERI_Size_of_GLF(solver) );
  glG = (double*)malloc( ERI_Size_of_GLF(solver) );
  if (dI4) {
    dglF[0] = (double*)malloc( ERI_Size_of_GLF(solver) );
    dglF[1] = (double*)malloc( ERI_Size_of_GLF(solver) );
    dglF[2] = (double*)malloc( ERI_Size_of_GLF(solver) );
    dglG[0] = (double*)malloc( ERI_Size_of_GLF(solver) );
    dglG[1] = (double*)malloc( ERI_Size_of_GLF(solver) );
    dglG[2] = (double*)malloc( ERI_Size_of_GLF(solver) );
  }

  c1[0] = orb1->c[0]; c1[1] = orb1->c[1]; c1[2] = orb1->c[2];
  c2[0] = orb2->c[0]; c2[1] = orb2->c[1]; c2[2] = orb2->c[2];
  c3[0] = orb3->c[0]; c3[1] = orb3->c[1]; c3[2] = orb3->c[2];
  c4[0] = orb4->c[0]; c4[1] = orb4->c[1]; c4[2] = orb4->c[2];

  /* expansion center */
  cx12 = ERI_Center_r2(orb1, orb2);
  cx34 = ERI_Center_r2(orb3, orb4);

  /* R */
  ERI_Coordinate_Transform(R, a1, a2, a3, a4, c1, c2, c3, c4, cx12, cx34);

  if (NULL==dI4) {
    ERI_Overlap(solver, glF, NULL, orb1, orb2, cx12); 
    ERI_Overlap(solver, glG, NULL, orb3, orb4, cx34);
    ERI_Integral_GL(solver, I4, glF, glG, R[0], R[1], R[2], scr, lmax_gl);
  } else {
    ERI_Overlap(solver, glF, dglF, orb1, orb2, cx12); 
    ERI_Overlap(solver, glG, dglG, orb3, orb4, cx34);
    ERI_Integral_GL_d(solver, I4, dI4, glF, glG, dglF, dglG,
                      R[0], R[1], R[2], cx12, cx34, 1e-10, scr, lmax_gl);
  }
 
  /* release workspace memory */
  free(glF);
  free(glG);
  if (dI4) {
    free(dglF[0]); 
    free(dglF[1]); 
    free(dglF[2]);
    free(dglG[0]); 
    free(dglG[1]); 
    free(dglG[2]);
  }
  
  STEPTRACE("ERI_Integral: out");
}



  
/*----------------------------------------------------------------------
  orbital_T

  kinetic energy
----------------------------------------------------------------------*/
static double orbital_T(
  ERI_t        *solver,
  const double *fr,     
  const double *xr,     
  int           nrgrid, 
  int           l       
)
{
  int i;
  double k, af, ke, dk;
  double *fk;
  
  int nkgrid;
  const double *kmesh;

  nkgrid = ERI_ngrid(solver);
  kmesh  = ERI_Mesh_Array_k(solver);
 
  fk = (double*)malloc( ERI_Size_of_Orbital(solver) );

  ERI_Transform_Orbital(solver, fk, fr, xr, nrgrid, l);

  ke = 0.0;
  for (i=0; i<nkgrid; i++) {
    k  = ERI_Mesh_k(solver, i);
    dk = ERI_Mesh_dk(solver, i);
    af = fk[2*i+0]*fk[2*i+0]+fk[2*i+1]*fk[2*i+1];
    ke += 8.0*(k*k*k*k)*af*dk;
  }
  
  free(fk);

  return ke;
}




/*----------------------------------------------------------------------
  orbital_r2

  <r^2>
----------------------------------------------------------------------*/
static double orbital_r2(
  const double *fr,
  const double *xr,
  int           ngrid
)
{
  int ir;
  double norm, br2k, x0, x1, r2, r, f;

  norm = 0.0;
  br2k = 0.0;
  for (ir=0; ir<ngrid-1; ir++) {
    x0 = xr[ir];
    x1 = xr[ir+1];
    r  = 0.5*(x0+x1);
    r2 = r*r;
    f  = 0.5*(fr[ir]+fr[ir+1]);
    norm += r2*f*f*(x1-x0);
    br2k += r2*r2*f*f*(x1-x0);
  }

  return br2k/norm;
}




/*----------------------------------------------------------------------
  orbital_r

  <r>
----------------------------------------------------------------------*/
static double orbital_r(
  const double *fr,
  const double *xr,
  int           ngrid
)
{
  int ir;
  double norm, brk, x0, x1, r, r2, f;

  norm = 0.0;
  brk  = 0.0;
  for (ir=0; ir<ngrid-1; ir++) {
    x0 = xr[ir];
    x1 = xr[ir+1];
    r  = 0.5*(x0+x1);
    r2 = r*r;
    f  = 0.5*(fr[ir]+fr[ir+1]);
    norm += r2*f*f*(x1-x0);
    brk  += r2*r*f*f*(x1-x0);
  }

  return brk/norm;
}




void ERI_Coordinate_Transform(
  double out_R[3],  /* (OUT) displacement in spherical coord. */
  double out_a1[3], /* (OUT) translated center 1 in spherical coord. */
  double out_a2[3], /* (OUT) translated center 2 in spherical coord. */
  double out_a3[3], /* (OUT) translated center 3 in spherical coord. */
  double out_a4[3], /* (OUT) translated center 4 in spherical coord. */
  const double in_a1[3],  /* (IN) center 1 in Cartesian coord. */
  const double in_a2[3],  /* (IN) center 2 in Cartesian coord. */
  const double in_a3[3],  /* (IN) center 3 in Cartesian coord. */
  const double in_a4[3],  /* (IN) center 4 in Cartesian coord. */
  double x12,   /* (IN) ratio */
  double x34    /* (IN) ratio */
)
{
  double c12[3], c34[3];
  double a1[3], a2[3], a3[3], a4[3];

  c12[0] = x12*in_a1[0] + (1.0-x12)*in_a2[0];
  c12[1] = x12*in_a1[1] + (1.0-x12)*in_a2[1];
  c12[2] = x12*in_a1[2] + (1.0-x12)*in_a2[2];
  c34[0] = x34*in_a3[0] + (1.0-x34)*in_a4[0];
  c34[1] = x34*in_a3[1] + (1.0-x34)*in_a4[1];
  c34[2] = x34*in_a3[2] + (1.0-x34)*in_a4[2];

  a1[0] = in_a1[0]-c12[0]; a1[1] = in_a1[1]-c12[1]; a1[2] = in_a1[2]-c12[2];
  a2[0] = in_a2[0]-c12[0]; a2[1] = in_a2[1]-c12[1]; a2[2] = in_a2[2]-c12[2];
  a3[0] = in_a3[0]-c34[0]; a3[1] = in_a3[1]-c34[1]; a3[2] = in_a3[2]-c34[2];
  a4[0] = in_a4[0]-c34[0]; a4[1] = in_a4[1]-c34[1]; a4[2] = in_a4[2]-c34[2];

  coord_spherical(a1[0], a1[1], a1[2], &out_a1[0], &out_a1[1], &out_a1[2]);
  coord_spherical(a2[0], a2[1], a2[2], &out_a2[0], &out_a2[1], &out_a2[2]);
  coord_spherical(a3[0], a3[1], a3[2], &out_a3[0], &out_a3[1], &out_a3[2]);
  coord_spherical(a4[0], a4[1], a4[2], &out_a4[0], &out_a4[1], &out_a4[2]);

  coord_spherical(c34[0]-c12[0], c34[1]-c12[1], c34[2]-c12[2], 
                  &out_R[0], &out_R[1], &out_R[2]);
}

 


/*----------------------------------------------------------------------
  ERI_Center_r2
----------------------------------------------------------------------*/
double ERI_Center_r2(
  const ERI_Orbital_t *pA,
  const ERI_Orbital_t *pB
)
{
   double d1 = orbital_r2(pA->fr, pA->xr, pA->ngrid);
   double d2 = orbital_r2(pB->fr, pB->xr, pB->ngrid);

   return d2/(d1+d2);
}




/*----------------------------------------------------------------------
  ERI_Center_DK
----------------------------------------------------------------------*/
double ERI_Center_DK(
  ERI_t *solver,
  const ERI_Orbital_t *pA,
  const ERI_Orbital_t *pB
)
{
  double d1 = orbital_r2(pA->fr, pA->xr, pA->ngrid);
  double d2 = orbital_r2(pB->fr, pB->xr, pB->ngrid);
  double k1 = orbital_T(solver, pA->fr, pA->xr, pA->ngrid, pA->l);
  double k2 = orbital_T(solver, pB->fr, pB->xr, pB->ngrid, pB->l);
   
  return (d2*k1)/(d2*k1+d1*k2);
}


/* EOF */
