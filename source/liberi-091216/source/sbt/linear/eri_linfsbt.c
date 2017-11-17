/*----------------------------------------------------------------------
  eri_linfsbt.c

  Coded by M. Toyoda, June 2009, JAIST/RCIS.
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <assert.h>
#include "eri_def.h"
#include "eri_linfsbt.h"
#include "eri_def.h"


static double g_tmp_in[ERI_NGRIDMAX*2];
static double g_tmp_out[ERI_NGRIDMAX];
static double g_tmp_work[ERI_NGRIDMAX*32];


struct ERI_LinFSBT_struct {
  int     nmesh;       /* number of mesh points */
  double *xr;          /* r-mesh */
  double *xk;          /* k-mesh */
  int     max_order;
  double *powk;        /* powers of k */
  double *powr;        /* powers of r */
  int     pownmax;    
  double *fbuf;
  double *F;          /* Fourier cosie/sine transform  */
  double  gamma[12];    /* Integral */
  double *I;           /* definite integrals */
};


/*----------------------------------------------------------------------
  Evaluation of starting value 
----------------------------------------------------------------------*/
static double recA(double r, double k, unsigned int n)
{
  double ret;

  if (0 == n) { 
    ret = sin(k*r)/r; 
  } else if (0 == n%2) {
    /* even */
    ret = (-(double)n*recA(r, k, n-1) + pow(k,(double)n)*sin(k*r))/r;
  } else {
    /* odd */
    ret = ((double)n*recA(r, k, n-1) - pow(k,(double)n)*cos(k*r))/r;
  }

  return ret;
}


static void  S_start(
  double       *out_r,
  double       *out_i,
  const double *fr,
  const double *xr,
  int           nmesh,
  double        k,
  int           n
)
{
  int i;
  double r, dr;

  dr = xr[1]-xr[0];
  *out_r = 0.0;
  *out_i = 0.0;

  for (i=0; i<nmesh; i++) {
    r = xr[i];
    *out_r += recA(r, k, n)*fr[2*i+0]*r*r*dr;
    *out_i += recA(r, k, n)*fr[2*i+1]*r*r*dr;
  }
}



ERI_INLINE double S0_asymptotic(
  double       k,
  const double gm[6],
  int          n
)
{
  double ret;

  if (0 == n%2) {
    /* even */
    ret = pow(k,(double)(n+1))*( gm[0]/(double)(n+1) 
            - gm[2]*k*k/2.0/(double)(n+3) +gm[4]*k*k*k*k/24.0/(double)(n+5) );
  } else {
    /* odd */
    ret = pow(k,(double)(n+2))*( gm[1]/(double)(n+2) 
            - gm[3]*k*k/6.0/(double)(n+4) +gm[5]*k*k*k*k/120.0/(double)(n+6) );
  }

  return ret;
}



static void segment_cubic_even(
  double       *work,
  const double *F,
  int           nmesh,
  const double *xk
)
{
  int i;
 
  double ddx, dyCrddx, dyCiddx;
  double k0, k0k0, k0k0k0, k1, k1k1, k1k1k1;
  double yC0r, yC1r, dyC0r, dyC1r;
  double yC0i, yC1i, dyC0i, dyC1i;
  double c2_cr, c3_cr, c2_ci, c3_ci;
  double A_cr, B_cr, C_cr, A_ci, B_ci, C_ci;
  
  double *tmp_c0r = &work[0];
  double *tmp_c0i = &work[nmesh*4];
  double *tmp_c1r = &work[nmesh*8];
  double *tmp_c1i = &work[nmesh*12];

  const double *C = F;
  const double *dC = &F[4*nmesh];

  k1     = xk[0];
  k1k1   = k1*k1;
  k1k1k1 = k1k1*k1;

  yC1r   = C[0];
  yC1i   = C[1];
  dyC1r  = dC[0];
  dyC1i  = dC[1]; 

  for (i=1; i<nmesh; i++) {
    k0 = k1;
    k1 = xk[i];

    yC0r  = yC1r;
    yC0i  = yC1i;
    yC1r  = C[2*i+0];
    yC1i  = C[2*i+1];

    dyC0r = dyC1r;
    dyC0i = dyC1i;
    dyC1r = dC[2*i+0];
    dyC1i = dC[2*i+1];
  
    ddx = 1.0/(k1-k0);
    dyCrddx = (yC1r-yC0r)*ddx;
    dyCiddx = (yC1i-yC0i)*ddx;

    k0k0   = k1k1;
    k0k0k0 = k1k1k1;

    k1k1   = k1*k1;
    k1k1k1 = k1k1*k1;
 
    c2_cr = (3.0*dyCrddx - (dyC1r+2.0*dyC0r))*ddx;
    c2_ci = (3.0*dyCiddx - (dyC1i+2.0*dyC0i))*ddx;
    c3_cr = (-2.0*dyCrddx + (dyC1r+dyC0r))*ddx*ddx;
    c3_ci = (-2.0*dyCiddx + (dyC1i+dyC0i))*ddx*ddx;

    A_cr = yC0r - dyC0r*k0 + c2_cr*k0k0 - c3_cr*k0k0k0;
    A_ci = yC0i - dyC0i*k0 + c2_ci*k0k0 - c3_ci*k0k0k0;
    B_cr = dyC0r - 2.0*c2_cr*k0 + 3.0*c3_cr*k0k0;
    B_ci = dyC0i - 2.0*c2_ci*k0 + 3.0*c3_ci*k0k0;
    C_cr = c2_cr - 3.0*c3_cr*k0;
    C_ci = c2_ci - 3.0*c3_ci*k0;
    
    tmp_c0r[i*4+0] = A_cr;
    tmp_c0i[i*4+0] = A_ci;
    tmp_c0r[i*4+1] = B_cr*k0;
    tmp_c0i[i*4+1] = B_ci*k0;
    tmp_c0r[i*4+2] = C_cr*k0k0;
    tmp_c0i[i*4+2] = C_ci*k0k0;
    tmp_c0r[i*4+3] = c3_cr*k0k0k0;
    tmp_c0i[i*4+3] = c3_ci*k0k0k0;

    tmp_c1r[i*4+0] = A_cr;
    tmp_c1i[i*4+0] = A_ci;
    tmp_c1r[i*4+1] = B_cr*k1;
    tmp_c1i[i*4+1] = B_ci*k1;
    tmp_c1r[i*4+2] = C_cr*k1k1;
    tmp_c1i[i*4+2] = C_ci*k1k1;
    tmp_c1r[i*4+3] = c3_cr*k1k1k1;
    tmp_c1i[i*4+3] = c3_ci*k1k1k1;
  }
}


static void segment_cubic_odd(
  double       *work,
  const double *F,
  int           nmesh,
  const double *xk
)
{
  int i;
 
  double ddx, dySrddx, dySiddx;
  double k0, k0k0, k0k0k0, k1, k1k1, k1k1k1;
  double yS0r, yS1r, dyS0r, dyS1r;
  double yS0i, yS1i, dyS0i, dyS1i;
  double c2_sr,c3_sr, c2_si,c3_si;
  double A_sr, B_sr, C_sr, A_si, B_si, C_si;
  
  double *tmp_s0r = &work[nmesh*16];
  double *tmp_s0i = &work[nmesh*20];
  double *tmp_s1r = &work[nmesh*24];
  double *tmp_s1i = &work[nmesh*28];
  
  const double *S  = &F[2*nmesh];
  const double *dS = &F[6*nmesh];

  k1     = xk[0];
  k1k1   = k1*k1;
  k1k1k1 = k1k1*k1;

  yS1r   = S[0];
  yS1i   = S[1];
  dyS1r  = dS[0];
  dyS1i  = dS[1]; 

  for (i=1; i<nmesh; i++) {
    k0 = k1;
    k1 = xk[i];

    yS0r  = yS1r;
    yS0i  = yS1i;
    yS1r  = S[2*i+0];
    yS1i  = S[2*i+1];

    dyS0r = dyS1r;
    dyS0i = dyS1i;
    dyS1r = dS[2*i+0];
    dyS1i = dS[2*i+1];
     
    ddx = 1.0/(k1-k0);
    dySrddx = (yS1r-yS0r)*ddx;
    dySiddx = (yS1i-yS0i)*ddx;

    k0k0   = k1k1;
    k0k0k0 = k1k1k1;

    k1k1   = k1*k1;
    k1k1k1 = k1k1*k1;
 
    c2_sr = (3.0*dySrddx - (dyS1r+2.0*dyS0r))*ddx;
    c2_si = (3.0*dySiddx - (dyS1i+2.0*dyS0i))*ddx;
    c3_sr = (-2.0*dySrddx + (dyS1r+dyS0r))*ddx*ddx;
    c3_si = (-2.0*dySiddx + (dyS1i+dyS0i))*ddx*ddx;

    A_sr = yS0r - dyS0r*k0 + c2_sr*k0k0 - c3_sr*k0k0k0;
    A_si = yS0i - dyS0i*k0 + c2_si*k0k0 - c3_si*k0k0k0;
    B_sr = dyS0r - 2.0*c2_sr*k0 + 3.0*c3_sr*k0k0;
    B_si = dyS0i - 2.0*c2_si*k0 + 3.0*c3_si*k0k0;
    C_sr = c2_sr - 3.0*c3_sr*k0;
    C_si = c2_si - 3.0*c3_si*k0;

    tmp_s0r[i*4+0] = A_sr;
    tmp_s0i[i*4+0] = A_si;
    tmp_s0r[i*4+1] = B_sr*k0;
    tmp_s0i[i*4+1] = B_si*k0;
    tmp_s0r[i*4+2] = C_sr*k0k0;
    tmp_s0i[i*4+2] = C_si*k0k0;
    tmp_s0r[i*4+3] = c3_sr*k0k0k0;
    tmp_s0i[i*4+3] = c3_si*k0k0k0;

    tmp_s1r[i*4+0] = A_sr;
    tmp_s1i[i*4+0] = A_si;
    tmp_s1r[i*4+1] = B_sr*k1;
    tmp_s1i[i*4+1] = B_si*k1;
    tmp_s1r[i*4+2] = C_sr*k1k1;
    tmp_s1i[i*4+2] = C_si*k1k1;
    tmp_s1r[i*4+3] = c3_sr*k1k1k1;
    tmp_s1i[i*4+3] = c3_si*k1k1k1;
  }
}


static void segment_cubic(
  double       *work,
  const double *F,
  int           nmesh,
  const double *xk
)
{
  int i;
  double dx, dy, c0r, c0i, c1r, c1i, c2r, c2i, c3r, c3i;
 
  double ddx, dyCrddx, dyCiddx, dySrddx, dySiddx;
  double k0, k0k0, k0k0k0, k1, k1k1, k1k1k1;
  double yC0r, yC1r, yS0r, yS1r, dyC0r, dyC1r, dyS0r, dyS1r;
  double yC0i, yC1i, yS0i, yS1i, dyC0i, dyC1i, dyS0i, dyS1i;
  double c2_cr, c2_sr, c3_cr, c3_sr;
  double c2_ci, c2_si, c3_ci, c3_si;
  double A_cr, B_cr, C_cr, A_sr, B_sr, C_sr;
  double A_ci, B_ci, C_ci, A_si, B_si, C_si;
  
  double *tmp_c0r = &work[0];
  double *tmp_c0i = &work[nmesh*4];
  double *tmp_c1r = &work[nmesh*8];
  double *tmp_c1i = &work[nmesh*12];
  double *tmp_s0r = &work[nmesh*16];
  double *tmp_s0i = &work[nmesh*20];
  double *tmp_s1r = &work[nmesh*24];
  double *tmp_s1i = &work[nmesh*28];

  const double *C  = F;
  const double *S  = &F[2*nmesh];
  const double *dC = &F[4*nmesh];
  const double *dS = &F[6*nmesh];

  k1     = xk[0];
  k1k1   = k1*k1;
  k1k1k1 = k1k1*k1;

  yC1r   = C[0];
  yC1i   = C[1];
  yS1r   = S[0];
  yS1i   = S[1];
  dyC1r  = dC[0];
  dyC1i  = dC[1];
  dyS1r  = dS[0];
  dyS1i  = dS[1]; 

  for (i=1; i<nmesh; i++) {
    k0 = k1;
    k1 = xk[i];

    yC0r  = yC1r;
    yC0i  = yC1i;
    yC1r  = C[2*i+0];
    yC1i  = C[2*i+1]; 

    yS0r  = yS1r;
    yS0i  = yS1i;
    yS1r  = S[2*i+0];
    yS1i  = S[2*i+1];

    dyC0r = dyC1r;
    dyC0i = dyC1i;
    dyC1r = dC[2*i+0];
    dyC1i = dC[2*i+1]; 
  
    dyS0r = dyS1r;
    dyS0i = dyS1i;
    dyS1r = dS[2*i+0];
    dyS1i = dS[2*i+1]; 
   
    ddx = 1.0/(k1-k0);
    dyCrddx = (yC1r-yC0r)*ddx;
    dyCiddx = (yC1i-yC0i)*ddx;
    dySrddx = (yS1r-yS0r)*ddx;
    dySiddx = (yS1i-yS0i)*ddx;

    k0k0   = k1k1;
    k0k0k0 = k1k1k1;

    k1k1   = k1*k1;
    k1k1k1 = k1k1*k1;
 
    c2_cr = (3.0*dyCrddx - (dyC1r+2.0*dyC0r))*ddx;
    c2_ci = (3.0*dyCiddx - (dyC1i+2.0*dyC0i))*ddx;
    c3_cr = (-2.0*dyCrddx + (dyC1r+dyC0r))*ddx*ddx;
    c3_ci = (-2.0*dyCiddx + (dyC1i+dyC0i))*ddx*ddx;

    c2_sr = (3.0*dySrddx - (dyS1r+2.0*dyS0r))*ddx;
    c2_si = (3.0*dySiddx - (dyS1i+2.0*dyS0i))*ddx;
    c3_sr = (-2.0*dySrddx + (dyS1r+dyS0r))*ddx*ddx;
    c3_si = (-2.0*dySiddx + (dyS1i+dyS0i))*ddx*ddx;

    A_cr = yC0r - dyC0r*k0 + c2_cr*k0k0 - c3_cr*k0k0k0;
    A_ci = yC0i - dyC0i*k0 + c2_ci*k0k0 - c3_ci*k0k0k0;
    B_cr = dyC0r - 2.0*c2_cr*k0 + 3.0*c3_cr*k0k0;
    B_ci = dyC0i - 2.0*c2_ci*k0 + 3.0*c3_ci*k0k0;
    C_cr = c2_cr - 3.0*c3_cr*k0;
    C_ci = c2_ci - 3.0*c3_ci*k0;
    
    A_sr = yS0r - dyS0r*k0 + c2_sr*k0k0 - c3_sr*k0k0k0;
    A_si = yS0i - dyS0i*k0 + c2_si*k0k0 - c3_si*k0k0k0;
    B_sr = dyS0r - 2.0*c2_sr*k0 + 3.0*c3_sr*k0k0;
    B_si = dyS0i - 2.0*c2_si*k0 + 3.0*c3_si*k0k0;
    C_sr = c2_sr - 3.0*c3_sr*k0;
    C_si = c2_si - 3.0*c3_si*k0;

    tmp_c0r[i*4+0] = A_cr;
    tmp_c0i[i*4+0] = A_ci;
    tmp_c0r[i*4+1] = B_cr*k0;
    tmp_c0i[i*4+1] = B_ci*k0;
    tmp_c0r[i*4+2] = C_cr*k0k0;
    tmp_c0i[i*4+2] = C_ci*k0k0;
    tmp_c0r[i*4+3] = c3_cr*k0k0k0;
    tmp_c0i[i*4+3] = c3_ci*k0k0k0;

    tmp_c1r[i*4+0] = A_cr;
    tmp_c1i[i*4+0] = A_ci;
    tmp_c1r[i*4+1] = B_cr*k1;
    tmp_c1i[i*4+1] = B_ci*k1;
    tmp_c1r[i*4+2] = C_cr*k1k1;
    tmp_c1i[i*4+2] = C_ci*k1k1;
    tmp_c1r[i*4+3] = c3_cr*k1k1k1;
    tmp_c1i[i*4+3] = c3_ci*k1k1k1;

    tmp_s0r[i*4+0] = A_sr;
    tmp_s0i[i*4+0] = A_si;
    tmp_s0r[i*4+1] = B_sr*k0;
    tmp_s0i[i*4+1] = B_si*k0;
    tmp_s0r[i*4+2] = C_sr*k0k0;
    tmp_s0i[i*4+2] = C_si*k0k0;
    tmp_s0r[i*4+3] = c3_sr*k0k0k0;
    tmp_s0i[i*4+3] = c3_si*k0k0k0;

    tmp_s1r[i*4+0] = A_sr;
    tmp_s1i[i*4+0] = A_si;
    tmp_s1r[i*4+1] = B_sr*k1;
    tmp_s1i[i*4+1] = B_si*k1;
    tmp_s1r[i*4+2] = C_sr*k1k1;
    tmp_s1i[i*4+2] = C_si*k1k1;
    tmp_s1r[i*4+3] = c3_sr*k1k1k1;
    tmp_s1i[i*4+3] = c3_si*k1k1k1;
  }
}



static void recsum_cubic_twoway_new(
  double       *outI,
  const double *fr,
  const double *F,
  const double  gm[12],
  int           nmesh,
  int           nmax,
  const double *xr,
  const double *xk,
  int           ik0,
  const double  *powy,
  int            pownmax,
  double        *workspace /* nmesh*32 */
)
{
  int n, i;
  double sr, si, s0r, s0i;
 
  const double *tmp_cs0r, *tmp_cs1r; 
  const double *tmp_cs0i, *tmp_cs1i; 
  double powk0, powk1;

  double dn1, dn2, dn3, dn4, dn5, dn6, dn7;
  
  segment_cubic(workspace, F, nmesh, xk); 

  dn1 = 1.0; 
  dn2 = 1.0/2.0; 
  dn3 = 1.0/3.0; 
  dn4 = 1.0/4.0; 
  dn5 = 1.0/5.0; 
  dn6 = 1.0/6.0; 

  for (n=0; n<nmax; n++) {
    if (0==n%2) {
      tmp_cs0r = &workspace[0];
      tmp_cs0i = &workspace[nmesh*4];
      tmp_cs1r = &workspace[nmesh*8];
      tmp_cs1i = &workspace[nmesh*12];
    } else {
      tmp_cs0r = &workspace[nmesh*16];
      tmp_cs0i = &workspace[nmesh*20];
      tmp_cs1r = &workspace[nmesh*24];
      tmp_cs1i = &workspace[nmesh*28];
    }

    /* I_{n0} */
    if (ik0==0) {
      if (0 == n%2) {
        /* even */
        s0r = powy[n*nmesh]*gm[0]*dn1
              -powy[(n+2)*nmesh]*gm[4]/2.0*dn3
              +powy[(n+4)*nmesh]*gm[8]/24.0*dn5;
        s0i = powy[n*nmesh]*gm[1]*dn1
              -powy[(n+2)*nmesh]*gm[5]/2.0*dn3
              +powy[(n+4)*nmesh]*gm[9]/24.0*dn5;
      } else {
        /* odd */
        s0r = powy[(n+1)*nmesh]*gm[2]*dn2
             -powy[(n+3)*nmesh]*gm[6]/6.0*dn4
             +powy[(n+5)*nmesh]*gm[10]/120.0*dn6;
        s0i = powy[(n+1)*nmesh]*gm[3]*dn2
             -powy[(n+3)*nmesh]*gm[7]/6.0*dn4
             +powy[(n+5)*nmesh]*gm[11]/120.0*dn6;
      }
    } else {
      S_start(&s0r, &s0i, fr, xr, nmesh, xk[ik0], n);
    }
    outI[2*(n*nmesh + ik0)+0] = s0r/powy[n*nmesh+ik0];
    outI[2*(n*nmesh + ik0)+1] = s0i/powy[n*nmesh+ik0];

    /* Downward */
    sr = s0r;
    si = s0i;
    
    powk0 = powy[n*nmesh+ik0];
    for (i=ik0; i>0; i--) {
      powk1 = powk0;
      powk0 = powy[n*nmesh+i-1];

      sr -= powk1*(
              tmp_cs1r[i*4+0]*dn1
              + tmp_cs1r[i*4+1]*dn2
              + tmp_cs1r[i*4+2]*dn3
              + tmp_cs1r[i*4+3]*dn4 )
            -powk0*(
              tmp_cs0r[i*4+0]*dn1
              + tmp_cs0r[i*4+1]*dn2
              + tmp_cs0r[i*4+2]*dn3
              + tmp_cs0r[i*4+3]*dn4 );
 
      si -= powk1*(
              tmp_cs1i[i*4+0]*dn1
              + tmp_cs1i[i*4+1]*dn2
              + tmp_cs1i[i*4+2]*dn3
              + tmp_cs1i[i*4+3]*dn4 )
            -powk0*(
              tmp_cs0i[i*4+0]*dn1
              + tmp_cs0i[i*4+1]*dn2
              + tmp_cs0i[i*4+2]*dn3
              + tmp_cs0i[i*4+3]*dn4 );
  
      outI[2*(n*nmesh+i-1)+0] = sr/powk0;
      outI[2*(n*nmesh+i-1)+1] = si/powk0;
    }
 
    /* Upward */
    sr = s0r;
    si = s0i;

    powk1 = powy[n*nmesh+ik0];
    for (i=ik0+1; i<nmesh; i++) {
      powk0 = powk1;
      powk1 = powy[n*nmesh+i];

      sr += powk1*(
              tmp_cs1r[i*4+0]*dn1
              + tmp_cs1r[i*4+1]*dn2
              + tmp_cs1r[i*4+2]*dn3
              + tmp_cs1r[i*4+3]*dn4 )
            -powk0*(
              tmp_cs0r[i*4+0]*dn1
              + tmp_cs0r[i*4+1]*dn2
              + tmp_cs0r[i*4+2]*dn3
              + tmp_cs0r[i*4+3]*dn4 );

      si += powk1*(
              tmp_cs1i[i*4+0]*dn1
              + tmp_cs1i[i*4+1]*dn2
              + tmp_cs1i[i*4+2]*dn3
              + tmp_cs1i[i*4+3]*dn4 )
            -powk0*(
              tmp_cs0i[i*4+0]*dn1
              + tmp_cs0i[i*4+1]*dn2
              + tmp_cs0i[i*4+2]*dn3
              + tmp_cs0i[i*4+3]*dn4 );

      outI[2*(n*nmesh+i)+0] = sr/powk1;
      outI[2*(n*nmesh+i)+1] = si/powk1;
    }

    dn1 = dn2;
    dn2 = dn3; 
    dn3 = dn4;
    dn4 = dn5;
    dn5 = dn6;
    dn6 = 1.0/(double)(n+7.0);
  } /* end loop of n */
  
}



static void recsum_cubic_upward_new(
  double       *outI,
  const double *fr,
  const double *F,
  const double  gm[12],
  int           nmesh,
  int           nmax,
  const double *xr,
  const double *xk,
  const double  *powy,
  int            pownmax,
  double        *workspace, /* nmesh*16 */
  int            parity
)
{
  int n, i, j;
  double sr, si, s0r, s0i;
  double dn1, dn2, dn3, dn4, dn5, dn6, dn7;
 
  const double *tmp_cs0r, *tmp_cs1r; 
  const double *tmp_cs0i, *tmp_cs1i; 
  double powk0, powk1;


  dn1 = 1.0; 
  dn2 = 1.0/2.0; 
  dn3 = 1.0/3.0; 
  dn4 = 1.0/4.0; 
  dn5 = 1.0/5.0; 
  dn6 = 1.0/6.0; 


  if (1==parity) { /* --- odd only */
    dn1 = 1.0/2.0; 
    dn2 = 1.0/3.0; 
    dn3 = 1.0/4.0; 
    dn4 = 1.0/5.0; 
    dn5 = 1.0/6.0; 
    dn6 = 1.0/7.0; 

    segment_cubic_odd(workspace, F, nmesh, xk);
    tmp_cs0r = &workspace[nmesh*16];
    tmp_cs0i = &workspace[nmesh*20];
    tmp_cs1r = &workspace[nmesh*24];
    tmp_cs1i = &workspace[nmesh*28];
    for (n=1; n<nmax; n+=2) {
      s0r = powy[(n+1)*nmesh]*gm[2]*dn2
            -powy[(n+3)*nmesh]*gm[6]/6.0*dn4
            +powy[(n+5)*nmesh]*gm[10]/120.0*dn6;
      s0i = powy[(n+1)*nmesh]*gm[3]*dn2
            -powy[(n+3)*nmesh]*gm[7]/6.0*dn4
            +powy[(n+5)*nmesh]*gm[11]/120.0*dn6;
      sr = s0r;
      si = s0i;
      outI[2*(n*nmesh + 0)+0] = s0r/powy[n*nmesh+0];
      outI[2*(n*nmesh + 0)+1] = s0i/powy[n*nmesh+0];

      powk1 = powy[n*nmesh+0];
      for (i=1; i<nmesh; i++) {
        powk0 = powk1;
        powk1 = powy[n*nmesh+i];

        sr += powk1*( tmp_cs1r[i*4+0]*dn1 + tmp_cs1r[i*4+1]*dn2
                      + tmp_cs1r[i*4+2]*dn3 + tmp_cs1r[i*4+3]*dn4 )
             -powk0*( tmp_cs0r[i*4+0]*dn1 + tmp_cs0r[i*4+1]*dn2
                      + tmp_cs0r[i*4+2]*dn3 + tmp_cs0r[i*4+3]*dn4 );
        si += powk1*( tmp_cs1i[i*4+0]*dn1 + tmp_cs1i[i*4+1]*dn2
                      + tmp_cs1i[i*4+2]*dn3 + tmp_cs1i[i*4+3]*dn4 )
             -powk0*( tmp_cs0i[i*4+0]*dn1 + tmp_cs0i[i*4+1]*dn2
                      + tmp_cs0i[i*4+2]*dn3 + tmp_cs0i[i*4+3]*dn4 );

        outI[2*(n*nmesh+i)+0] = sr/powk1;
        outI[2*(n*nmesh+i)+1] = si/powk1;
      }

      dn1 = dn3;
      dn2 = dn4; 
      dn3 = dn5;
      dn4 = dn6;
      dn5 = 1.0/(double)(n+7.0);
      dn6 = 1.0/(double)(n+8.0);
    } /* end loop of n */
  } else if (2==parity) { /* --- even only */
    segment_cubic_even(workspace, F, nmesh, xk);
    tmp_cs0r = &workspace[nmesh*0];
    tmp_cs0i = &workspace[nmesh*4];
    tmp_cs1r = &workspace[nmesh*8];
    tmp_cs1i = &workspace[nmesh*12];
    for (n=0; n<nmax; n+=2) {
      s0r = powy[n*nmesh]*gm[0]*dn1
            -powy[(n+2)*nmesh]*gm[4]/2.0*dn3
            +powy[(n+4)*nmesh]*gm[8]/24.0*dn5;
      s0i = powy[n*nmesh]*gm[1]*dn1
            -powy[(n+2)*nmesh]*gm[5]/2.0*dn3
            +powy[(n+4)*nmesh]*gm[9]/24.0*dn5;
      sr = s0r;
      si = s0i;
      outI[2*(n*nmesh + 0)+0] = s0r/powy[n*nmesh+0];
      outI[2*(n*nmesh + 0)+1] = s0i/powy[n*nmesh+0];

      powk1 = powy[n*nmesh+0];
      for (i=1; i<nmesh; i++) {
        powk0 = powk1;
        powk1 = powy[n*nmesh+i];

        sr += powk1*( tmp_cs1r[i*4+0]*dn1 + tmp_cs1r[i*4+1]*dn2
                      + tmp_cs1r[i*4+2]*dn3 + tmp_cs1r[i*4+3]*dn4 )
             -powk0*( tmp_cs0r[i*4+0]*dn1 + tmp_cs0r[i*4+1]*dn2
                      + tmp_cs0r[i*4+2]*dn3 + tmp_cs0r[i*4+3]*dn4 );
        si += powk1*( tmp_cs1i[i*4+0]*dn1 + tmp_cs1i[i*4+1]*dn2
                      + tmp_cs1i[i*4+2]*dn3 + tmp_cs1i[i*4+3]*dn4 )
             -powk0*( tmp_cs0i[i*4+0]*dn1 + tmp_cs0i[i*4+1]*dn2
                      + tmp_cs0i[i*4+2]*dn3 + tmp_cs0i[i*4+3]*dn4 );

        outI[2*(n*nmesh+i)+0] = sr/powk1;
        outI[2*(n*nmesh+i)+1] = si/powk1;
      }

      dn1 = dn3;
      dn2 = dn4; 
      dn3 = dn5;
      dn4 = dn6;
      dn5 = 1.0/(double)(n+7.0);
      dn6 = 1.0/(double)(n+8.0);
    } /* end loop of n */
  } else if (3==parity) { /* --- both even and odd */
    segment_cubic(workspace, F, nmesh, xk);
    for (n=0; n<nmax; n++) {
      if (0==n%2) {
        tmp_cs0r = &workspace[nmesh*0];
        tmp_cs0i = &workspace[nmesh*4];
        tmp_cs1r = &workspace[nmesh*8];
        tmp_cs1i = &workspace[nmesh*12];

        s0r = powy[n*nmesh]*gm[0]*dn1
              -powy[(n+2)*nmesh]*gm[4]/2.0*dn3
              +powy[(n+4)*nmesh]*gm[8]/24.0*dn5;
        s0i = powy[n*nmesh]*gm[1]*dn1
              -powy[(n+2)*nmesh]*gm[5]/2.0*dn3
              +powy[(n+4)*nmesh]*gm[9]/24.0*dn5;
      } else {
        tmp_cs0r = &workspace[nmesh*16];
        tmp_cs0i = &workspace[nmesh*20];
        tmp_cs1r = &workspace[nmesh*24];
        tmp_cs1i = &workspace[nmesh*28];

        s0r = powy[(n+1)*nmesh]*gm[2]*dn2
              -powy[(n+3)*nmesh]*gm[6]/6.0*dn4
              +powy[(n+5)*nmesh]*gm[10]/120.0*dn6;
        s0i = powy[(n+1)*nmesh]*gm[3]*dn2
              -powy[(n+3)*nmesh]*gm[7]/6.0*dn4
              +powy[(n+5)*nmesh]*gm[11]/120.0*dn6;
      }
      sr = s0r;
      si = s0i;
      outI[2*(n*nmesh + 0)+0] = s0r/powy[n*nmesh+0];
      outI[2*(n*nmesh + 0)+1] = s0i/powy[n*nmesh+0];

      powk1 = powy[n*nmesh+0];
      for (i=1; i<nmesh; i++) {
        powk0 = powk1;
        powk1 = powy[n*nmesh+i];

        sr += powk1*( tmp_cs1r[i*4+0]*dn1 + tmp_cs1r[i*4+1]*dn2
                      + tmp_cs1r[i*4+2]*dn3 + tmp_cs1r[i*4+3]*dn4 )
             -powk0*( tmp_cs0r[i*4+0]*dn1 + tmp_cs0r[i*4+1]*dn2
                      + tmp_cs0r[i*4+2]*dn3 + tmp_cs0r[i*4+3]*dn4 );
        si += powk1*( tmp_cs1i[i*4+0]*dn1 + tmp_cs1i[i*4+1]*dn2
                      + tmp_cs1i[i*4+2]*dn3 + tmp_cs1i[i*4+3]*dn4 )
             -powk0*( tmp_cs0i[i*4+0]*dn1 + tmp_cs0i[i*4+1]*dn2
                      + tmp_cs0i[i*4+2]*dn3 + tmp_cs0i[i*4+3]*dn4 );

        outI[2*(n*nmesh+i)+0] = sr/powk1;
        outI[2*(n*nmesh+i)+1] = si/powk1;
      }

      dn1 = dn2;
      dn2 = dn3; 
      dn3 = dn4;
      dn4 = dn5;
      dn5 = dn6;
      dn6 = 1.0/(double)(n+7.0);
    } /* end loop of n */
  } else {
    abort();
  }
}




/*----------------------------------------------------------------------
  TOSBT_Calculat_Gamma

  Calculate Gamma factors which are defined as

  gamma[n] := int_0^infty fr(r) * r^n dr
----------------------------------------------------------------------*/
static void calculate_gamma(
  double        g[12],
  const double *fr,
  const double *xr,
  int           nmesh,
  const double *powr,
  int           parity
)
{
  int i;
  double x, dx, tmp, g2, g3, g4, g5, g6, g7;

  dx = xr[1]-xr[0];
  g[0] = 0.0; g[1] = 0.0; g[2]  = 0.0; g[3]  = 0.0;
  g[4] = 0.0; g[5] = 0.0; g[6]  = 0.0; g[7]  = 0.0;
  g[8] = 0.0; g[9] = 0.0; g[10] = 0.0; g[11] = 0.0;
 
  if (1==parity) { /* --- odd */
    for (i=0; i<nmesh; i++) {
      g[2]  += fr[2*i+0]*powr[2*nmesh+i]*dx;
      g[3]  += fr[2*i+1]*powr[2*nmesh+i]*dx;
      g[6]  += fr[2*i+0]*powr[4*nmesh+i]*dx;
      g[7]  += fr[2*i+1]*powr[4*nmesh+i]*dx;
      g[10] += fr[2*i+0]*powr[6*nmesh+i]*dx;
      g[11] += fr[2*i+1]*powr[6*nmesh+i]*dx;
    }   
  } else if (2==parity) { /* --- even */
    for (i=0; i<nmesh; i++) {
      g[0]  += fr[2*i+0]*powr[1*nmesh+i]*dx;
      g[1]  += fr[2*i+1]*powr[1*nmesh+i]*dx;
      g[4]  += fr[2*i+0]*powr[3*nmesh+i]*dx;
      g[5]  += fr[2*i+1]*powr[3*nmesh+i]*dx;
      g[8]  += fr[2*i+0]*powr[5*nmesh+i]*dx;
      g[9]  += fr[2*i+1]*powr[5*nmesh+i]*dx;
    }   
  } else if (3==parity) { /* --- both even and odd */
    for (i=0; i<nmesh; i++) {
      g[0]  += fr[2*i+0]*powr[1*nmesh+i]*dx;
      g[1]  += fr[2*i+1]*powr[1*nmesh+i]*dx;
      g[2]  += fr[2*i+0]*powr[2*nmesh+i]*dx;
      g[3]  += fr[2*i+1]*powr[2*nmesh+i]*dx;
      g[4]  += fr[2*i+0]*powr[3*nmesh+i]*dx;
      g[5]  += fr[2*i+1]*powr[3*nmesh+i]*dx;
      g[6]  += fr[2*i+0]*powr[4*nmesh+i]*dx;
      g[7]  += fr[2*i+1]*powr[4*nmesh+i]*dx;
      g[8]  += fr[2*i+0]*powr[5*nmesh+i]*dx;
      g[9]  += fr[2*i+1]*powr[5*nmesh+i]*dx;
      g[10] += fr[2*i+0]*powr[6*nmesh+i]*dx;
      g[11] += fr[2*i+1]*powr[6*nmesh+i]*dx;
    }   
  }
}




static void transform_f(
  double       *F,
  const double *fr,
  const double *meshx,
  unsigned int  nmesh,
  int           parity      /* 1: even only, 2: odd only, 3: both */
)
{
  int i, ieven, iodd;
  double x, dx;
  fftw_plan plan;

  double *tmp_in  = g_tmp_in;
  double *tmp_out = g_tmp_out;
 
  double *C  = F;
  double *S  = &F[2*nmesh];
  double *dC = &F[4*nmesh];
  double *dS = &F[6*nmesh];

  iodd  = 1 & parity; 
  ieven = 2 & parity; 

  dx = meshx[1]-meshx[0];
  for (i=0; i<nmesh; i++) { 
    x = meshx[i]; 
    tmp_in[i]       = 0.5*fr[2*i+0]*x*x*dx; 
    tmp_in[nmesh+i] = 0.5*fr[2*i+1]*x*x*dx; 
  }

  /* C and S */
  if (ieven) {
    plan = fftw_plan_r2r_1d(nmesh, &tmp_in[0], tmp_out, 
                            FFTW_REDFT11, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    for (i=0; i<nmesh; i++) { C[2*i+0] = tmp_out[i]; }
    plan = fftw_plan_r2r_1d(nmesh, &tmp_in[nmesh], tmp_out, 
                            FFTW_REDFT11, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    for (i=0; i<nmesh; i++) { C[2*i+1] = tmp_out[i]; }
  }
  if (iodd) {
    plan = fftw_plan_r2r_1d(nmesh, &tmp_in[0], tmp_out,
                            FFTW_RODFT11, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    for (i=0; i<nmesh; i++) { S[2*i+0] = tmp_out[i]; }
    plan = fftw_plan_r2r_1d(nmesh, &tmp_in[nmesh], tmp_out,
                            FFTW_RODFT11, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    for (i=0; i<nmesh; i++) { S[2*i+1] = tmp_out[i]; }
  }
 
  /* derivatives */
  for (i=0; i<nmesh; i++) { tmp_in[i] *= meshx[i]; }
  if (ieven) {
    plan = fftw_plan_r2r_1d(nmesh, &tmp_in[0], tmp_out,
                            FFTW_RODFT11, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    for (i=0; i<nmesh; i++) { dC[2*i+0] = -tmp_out[i]; }
    plan = fftw_plan_r2r_1d(nmesh, &tmp_in[nmesh], tmp_out,
                            FFTW_RODFT11, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    for (i=0; i<nmesh; i++) { dC[2*i+1] = -tmp_out[i]; }
  }
  if (iodd) {
    plan = fftw_plan_r2r_1d(nmesh, &tmp_in[0], tmp_out, 
                            FFTW_REDFT11, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    for (i=0; i<nmesh; i++) { dS[2*i+0] = tmp_out[i]; }
    plan = fftw_plan_r2r_1d(nmesh, &tmp_in[nmesh], tmp_out, 
                            FFTW_REDFT11, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    for (i=0; i<nmesh; i++) { dS[2*i+1] = tmp_out[i]; }
  }
}




static void transform_switch(
  double              *I,
  const double        *F,
  const double         gamma[12],
  const double        *fr,
  const double        *meshx,
  const double        *meshy,
  unsigned int         nmesh,
  unsigned int         nmax,
  unsigned int         ik0,
  int                  parity,     /* 1:even only, 2:odd only, 3:both */
  const double        *powy,
  int                  pownmax
)
{
  int n, ieven, iodd;
  double *tmp_work = g_tmp_work;

  if (ik0==0) {
    recsum_cubic_upward_new(I, fr, F, gamma, nmesh, nmax,  
                            meshx, meshy, powy, pownmax, tmp_work, parity);
  } else {
    recsum_cubic_twoway_new(I, fr, F, gamma, nmesh, nmax, 
                            meshx, meshy, ik0, powy, pownmax, tmp_work);
  }
}



static void transform_in(
  ERI_LinFSBT_t *ptr,
  const double  *in,
  int            direction, 
  int            parity)
{
  int i;
  const double *meshx, *meshy, *powy, *powx;

  int          nmax    = ptr->max_order + 1;
  int          nmesh   = ptr->nmesh;
  int          pownmax = ptr->pownmax;

  /* direction of transform */
  if (ERI_SBT_FORWARD == direction) {
    meshx = ptr->xr;
    meshy = ptr->xk;
    powx  = ptr->powr;
    powy  = ptr->powk;
  } else if (ERI_SBT_BACKWARD == direction) {
    meshx = ptr->xk;
    meshy = ptr->xr;
    powx  = ptr->powk;
    powy  = ptr->powr;
  } else {
    abort();
  }

  transform_f(ptr->F, in, meshx, nmesh, parity);

  /* gamma terms */
  calculate_gamma(ptr->gamma, in, meshx, nmesh, powx, parity);
 
  transform_switch(ptr->I, ptr->F, ptr->gamma, in, 
                   meshx, meshy, nmesh, nmax, 0, parity, 
                   powy, pownmax);
}




static void transform_out(ERI_LinFSBT_t *ptr, double *out, int l)
{
  int i, n, j;
  double fact, twop, twon;

  int           nmesh = ptr->nmesh;
  const double *I     = ptr->I;
  const double *fbuf  = ptr->fbuf;
  int           nmax  = ptr->max_order;

  
  for (i=0; i<2*nmesh; i++) { out[i] = 0.0; }

  if (0==l%2) { /* --- even */
    n = l/2;
    for (j=0; j<=n; j++) {
      for (i=0; i<nmesh; i++) {
        out[2*i+0] += fbuf[l*nmax+j]*I[2*(2*j*nmesh+i)+0]; 
        out[2*i+1] += fbuf[l*nmax+j]*I[2*(2*j*nmesh+i)+1]; 
      }
    }
  } else { /* --- odd */
    n = (l-1)/2;
    for (j=0; j<=n; j++) {
      for (i=0; i<nmesh; i++) {
        out[2*i+0] += fbuf[l*nmax+j]*I[2*((2*j+1)*nmesh+i)+0]; 
        out[2*i+1] += fbuf[l*nmax+j]*I[2*((2*j+1)*nmesh+i)+1]; 
      }
    }
  }
}



/*----------------------------------------------------------------------
                                                               interface 
----------------------------------------------------------------------*/
void ERI_LinFSBT_Free(void *p)
{
  ERI_LinFSBT_t *ptr = (ERI_LinFSBT_t*)p;
  if (ptr) {
    if (ptr->xr)   { free(ptr->xr);       }
    if (ptr->xk)   { free(ptr->xk);       }
    if (ptr->powr) { free(ptr->powr);     }
    if (ptr->powk) { free(ptr->powk);     }
    if (ptr->fbuf) { free(ptr->fbuf);     }
    if (ptr->I)    { free(ptr->I);        }
    if (ptr->F)    { free(ptr->F);        }
    free(ptr);
  }
}


double ERI_LinFSBT_Mesh_r(void *p, int i)
{
  ERI_LinFSBT_t *ptr = (ERI_LinFSBT_t*)p;
  assert( i < ptr->nmesh );
  return ptr->xr[i];
}


double ERI_LinFSBT_Mesh_k(void *p, int i)
{
  ERI_LinFSBT_t *ptr = (ERI_LinFSBT_t*)p;
  assert( i < ptr->nmesh );
  return ptr->xk[i];
}


int ERI_LinFSBT_Mesh_nmesh(void *p)
{
  ERI_LinFSBT_t *ptr = (ERI_LinFSBT_t*)p;
  return ptr->nmesh;
}


int ERI_LinFSBT_lmax(void *p)
{
  ERI_LinFSBT_t *ptr = (ERI_LinFSBT_t*)p;
  return ptr->max_order;
}


const double* ERI_LinFSBT_Mesh_Array_r(void *p)
{
  ERI_LinFSBT_t *ptr = (ERI_LinFSBT_t*)p;
  return ptr->xr;
}


const double* ERI_LinFSBT_Mesh_Array_k(void *p)
{
  ERI_LinFSBT_t *ptr = (ERI_LinFSBT_t*)p;
  return ptr->xk;
}


double ERI_LinFSBT_Mesh_dr(void *p, int i)
{
  ERI_LinFSBT_t *ptr = (ERI_LinFSBT_t*)p;
  return ptr->xr[1] - ptr->xr[0];
}


double ERI_LinFSBT_Mesh_dk(void *p, int i)
{
  ERI_LinFSBT_t *ptr = (ERI_LinFSBT_t*)p;
  return ptr->xk[1] - ptr->xk[0];
}


void ERI_LinFSBT_Transform_Input(void *p, const double *in, int dir)
{
  ERI_LinFSBT_t *ptr = (ERI_LinFSBT_t*)p;
  transform_in(ptr, in, dir, 3);
}


void ERI_LinFSBT_Transform_Output(void *p, double *out, int l)
{
  ERI_LinFSBT_t *ptr = (ERI_LinFSBT_t*)p;
  transform_out(ptr, out, l);
}


void ERI_LinFSBT_Transform(
  void         *p,
  double       *out, 
  const double *in, 
  int           l,
  int           dir
)
{
  ERI_LinFSBT_t *ptr = (ERI_LinFSBT_t*)p;
  int parity = (0==l%2) ? 2 : 1;

  transform_in(ptr, in, dir, parity);
  transform_out(ptr, out, l);
}




ERI_SBT_t* ERI_LinFSBT_Init(
  int    nmesh,      
  int    max_order,
  double range_r
)
{
  int i, n, p, hlfn;
  double dr, dk, fact, twon, twop;
  ERI_SBT_t *ptr;
  ERI_LinFSBT_t *pinst;
  
  int pownmax = max_order+10;

  assert( 0 < nmesh );
  assert( 0 < max_order );
  

  /* allocation */ 
  pinst = (ERI_LinFSBT_t*)malloc(sizeof(ERI_LinFSBT_t));
  if (NULL == pinst) return NULL;

  pinst->xr   = (double*)malloc(sizeof(double)*nmesh);
  pinst->xk   = (double*)malloc(sizeof(double)*nmesh);
  pinst->powk = (double*)malloc(sizeof(double)*nmesh*pownmax);
  pinst->powr = (double*)malloc(sizeof(double)*nmesh*pownmax);
  pinst->fbuf = (double*)malloc(sizeof(double)*max_order*max_order);
  pinst->I    = (double*)malloc(sizeof(double)*nmesh*(max_order+1)*2);
  pinst->F    = (double*)malloc(sizeof(double)*nmesh*8);
  if ( NULL == pinst->xr   || NULL == pinst->xk 
#if 0
    || NULL == pinst->in   || NULL == pinst->out 
    || NULL == pinst->work
#endif
    || NULL == pinst->powr 
    || NULL == pinst->powk || NULL == pinst->I || NULL == pinst->F ) {
    ERI_LinFSBT_Free(pinst);
    return NULL;
  }


  /* mesh */
  dr = range_r/(double)nmesh;  
  dk = PI/(double)nmesh/dr;
  for (i=0; i<nmesh; i++) {
    pinst->xr[i] = dr*(0.5+(double)i); 
    pinst->xk[i] = dk*(0.5+(double)i);  
    for (n=0; n<pownmax; n++) {
      pinst->powr[n*nmesh+i] = pow(pinst->xr[i],(double)n+1);
      pinst->powk[n*nmesh+i] = pow(pinst->xk[i],(double)n+1);
    }
  }

  /* factorial */
  for (n=0; n<max_order; n++) {
    if (n%2==0) {
      /* even */
      hlfn = n/2;
      twon = 2.0*(double)hlfn;
    
      /* (-1)^n*(4n-1)!!/(2n)! */
      fact = 1.0;
      for (p=0; p<2*hlfn; p++) { fact *= (double)(2*p+1)/(double)(p+1); }
      if (hlfn%2==1) fact = -fact;

      for (p=hlfn; p>=0; p--) {
        pinst->fbuf[n*max_order+p] = fact; 
        twop = 2.0*(double)p;
        fact *= -(twop*(twop-1.0))/((twon+twop-1.0)*(twon-twop+2.0));
      }
    } else { 
      /* odd */
      hlfn = (n-1)/2;
      twon = 2.0*(double)hlfn;
    
      /* (2n+1)!!/(2n)!! */
      fact = 1.0;
      for (p=1; p<=hlfn; p++) { fact *= (double)(2*p+1)/(double)(2*p); }

      for (p=0; p<=hlfn; p++) {
        pinst->fbuf[n*max_order+p] = fact;
        twop = 2.0*(double)p;
        fact *= -(twon+twop+3.0)*(twon-twop)/((twop+3.0)*(twop+2.0));
      }
    }
  }
 
  pinst->nmesh     = nmesh;
  pinst->max_order = max_order;
  pinst->pownmax   = pownmax;

  ptr = ERI_SBT_Init(
    ERI_LinFSBT_Free,
    ERI_LinFSBT_Mesh_nmesh,
    ERI_LinFSBT_lmax,
    ERI_LinFSBT_Mesh_r,
    ERI_LinFSBT_Mesh_k,
    ERI_LinFSBT_Mesh_Array_r,
    ERI_LinFSBT_Mesh_Array_k,
    ERI_LinFSBT_Mesh_dr,
    ERI_LinFSBT_Mesh_dk,
    ERI_LinFSBT_Transform_Input,
    ERI_LinFSBT_Transform_Output,
    ERI_LinFSBT_Transform,
    pinst);
  
  if (NULL==ptr) {
    ERI_LinFSBT_Free(pinst);
    return NULL;
  }
 
  return ptr;
}



size_t ERI_LinFSBT_Required_Size(
  int nmesh,      
  int max_order
)
{
  int pownmax = max_order+10;

  assert( 0 < nmesh );
  assert( 0 < max_order );

  return ERI_SBT_Required_Size()
         + sizeof(ERI_LinFSBT_t)
         + sizeof(double)*nmesh                 /* xr */
         + sizeof(double)*nmesh                 /* xk */
         + sizeof(double)*nmesh*pownmax         /* powk */
         + sizeof(double)*nmesh*pownmax         /* powr */
         + sizeof(double)*max_order*max_order   /* fbuf */
         + sizeof(double)*nmesh*(max_order+1)*2 /* I */
         + sizeof(double)*nmesh*8               /* F */
    ;
}


/* EOF */
