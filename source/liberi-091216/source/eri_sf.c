/*----------------------------------------------------------------------
  eri_sf.c

  Special functions used in LIBERI.
  
  This is based on Spherical_Bessel.c and openmx_common.c in OpenMX 
  3.4 package.

  22 Nov. 2001, opnemx_common.c by T. Ozaki
  08 Nov. 2005  Spherical_Bessel.c by T.Ozaki
  16 Feb. 2009, M. Toyoda
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <assert.h>
#include "eri_sf.h"
#include "eri_def.h"


#define LMAX_MAX 512
#define MAX_ITR 50
#define ERR_MAX 1e-12

#define CG_JMAX 30


static void sph_bessel_asym(double x, int l, double *sb, double *dsb);
static void sph_bessel_drec(double x, int l, double *sb, double *dsb);
static void sph_bessel_arec(double x, int l, double *sb, double *dsb);

static double CG(int j1, int m1, int j2, int m2, int j, int m);



/*----------------------------------------------------------------------
  rotation matrix for Complex -> Real conversion of SH

  Z(l, m) = a(m)*Y(l, m) + b(m)*Y(l, -m)

  where
    
    a(m) = 1/sqrt(2) (m>0), 1 (m=0), and i(-1)^m/sqrt(2) (m<0)
    b(m) = (-1)^m/sqrt(2) (m>0), 0 (m=0), and -i/sqrt(2) (m<0)
    
    a(m) = (-1)^m/sqrt(2) (m>0), 1 (m=0), and i/sqrt(2) (m<0)
    b(m) = 1/sqrt(2) (m>0), 0 (m=0), and -i(-1)^m/sqrt(2) (m<0)
----------------------------------------------------------------------*/
static void ERI_RCSH_Coeff(int m, double a[2], double b[2])
{
  double oost;

  oost = 1.0/sqrt(2.0);
  
  if (m>0) {
    if (m%2==0) {
      a[0] =  oost; a[1] = 0.0;
      b[0] =  oost; b[1] = 0.0;
    } else {
      a[0] = -oost; a[1] = 0.0;
      b[0] =  oost; b[1] = 0.0;
    }
  } else if (m<0) {
    if (abs(m)%2==0) {
      a[0] = 0.0; a[1] =  oost;
      b[0] = 0.0; b[1] = -oost;
    } else {
      a[0] = 0.0; a[1] =  oost;
      b[0] = 0.0; b[1] =  oost;
    }
  } else {
    /* m == 0 */
    a[0] = 1.0; a[1] = 0.0;
    b[0] = 0.0; b[1] = 0.0;
  }

#if 0
  if (m>0) {
    if (m%2==0) {
      a[0] =  oost; a[1] = 0.0;
      b[0] =  oost; b[1] = 0.0;
    } else {
      a[0] =  oost; a[1] = 0.0;
      b[0] = -oost; b[1] = 0.0;
    }
  } else if (m<0) {
    if (m%2==0) {
      a[0] = 0.0; a[1] =  oost;
      b[0] = 0.0; b[1] = -oost;
    } else {
      a[0] = 0.0; a[1] = -oost;
      b[0] = 0.0; b[1] = -oost;
    }
  } else {
    /* m == 0 */
    a[0] = 1.0; a[1] = 0.0;
    b[0] = 0.0; b[1] = 0.0;
  }
#endif
}


/*----------------------------------------------------------------------
  rotation matrix for Real -> Complex conversion of SH

  Y(l, m) = c(m)*Z(l, m) + d(m)*Z(l, -m)

  where
    
    c(m) = 1/sqrt(2) (m>0), 1 (m=0), and -i(-1)^m/sqrt(2) (m<0)
    f(m) = i/sqrt(2) (m>0), 0 (m=0), and (-1)^m/sqrt(2) (m<0)
----------------------------------------------------------------------*/
static void ERI_RCSH_Coeff_Inverse(int m, double c[2], double d[2])
{
  double na, nb, a[2], b[2];

  if (0==m) {
    c[0] = 1.0; c[1] = 0.0;
    d[0] = 0.0; d[1] = 0.0;
  } else { 
    ERI_RCSH_Coeff(m, a, b);
    na = a[0]*a[0] + a[1]*a[1];
    nb = b[0]*b[0] + b[1]*b[1];
    c[0] =  0.5*a[0]/na;
    c[1] = -0.5*a[1]/na;
    d[0] =  0.5*b[0]/nb;
    d[1] = -0.5*b[1]/nb;
  }
}




/*----------------------------------------------------------------------
  ERI_Spherical_Bessel

  Routine to calculate the spherical Bessel function of the first kind
  j_l(x) of the order from 0 to l and their derivatives.
  
  Descending recurrence algoritm is used for small x in order to avoid 
  accumulation of rounding errors, whereas ascending recurrence is used 
  for large x because rounding error can be negligible and the 
  descending recurrence becomes too slow for large x.

  Reference:
  A. jablonski, "Numerical Evaluation of Spherical Bessel Functions of 
  the First Kind", Journal of Computational Physics Vol. 111, pp. 256
  (1994).
----------------------------------------------------------------------*/
void ERI_Spherical_Bessel(
  double  x,  /* (IN) argument x >= 0*/
  int     l,  /* (IN) positive integer */
  double *sb, /* (OUT) function values of order from 0 to l */
  double *dsb /* (OUT) derivatives of order from 0 to l */
) 
{
  const double xmin = 1e-10;
  const double xthresh = 100.0;
  
  /* negative x value is invalid */
  assert( x >= 0.0 );
  assert( l >= 0 );

  if (x<xmin) {
    /* x is near zero */
    sph_bessel_asym(x, l, sb, dsb);
  } else if (x<xthresh) {
    /* descending recurrence */
    sph_bessel_drec(x, l, sb, dsb);
  } else {
    /* ascending recurrence */
    sph_bessel_arec(x, l, sb, dsb);
  }
}




static void sph_bessel_asym(
  double  x,
  int     l,
  double *sb,
  double *dsb
)
{
  int i;
  double f = 1.0;
  
  sb[0] = 1.0;
  for (i=1; i<=l+1; i++) { 
    f += 2.0; 
    sb[i] = sb[i-1]*x/f; 
  }

  dsb[0] = 0.0;
  for (i=1; i<=l; i++) {
    dsb[i] = ((double)i*sb[i-1]-(double)(i+1)*sb[i+1])
              /(double)(2*i+1);
  }
}




static void sph_bessel_drec(
  double  x,
  int     l,
  double *sb,
  double *dsb
)
{
  int i, j, lmax;
  double j0, j1, sf, tmp, si, co, ix, ix2, tsb[1024], huge;
  
  huge = sqrt(DBL_MAX);

  /* find an appropriate nmax */
  lmax = l + 3*(int)x + 20;
  assert( lmax+1 < LMAX_MAX );
  if (lmax<LMAX_MAX) lmax = LMAX_MAX; 

  /* initial values*/
  tsb[lmax]   = 0.0;
  tsb[lmax-1] = 1.0e-14;

  /* downward recurrence from nmax-2 to lmax+2 */
  for (i=lmax-1; (l+2)<i; i--) {
    tsb[i-1] = (2.0*(double)i+1.0)/x*tsb[i]-tsb[i+1];
    if (huge<tsb[i-1]){
      tsb[i]   /= tsb[i-1];
      tsb[i-1] = 1.0;
    }
  }

  /* downward recurrence from lmax+1 to 0 */
  tsb[l+3] /= tsb[l+2];
  tsb[l+2] = 1.0;
  for (i=l+2; 0<i; i--) {
    tsb[i-1] = (2.0*(double)i+1.0)/x*tsb[i] - tsb[i+1];
    if (huge<tsb[i-1]){
      tmp = tsb[i-1];
      for (j=i-1; j<=l+1; j++) { tsb[j] /= tmp; }
    }
  }

  /* normalization */
  si = sin(x);
  co = cos(x);
  ix = 1.0/x;
  ix2 = ix*ix;
  j0 = si*ix;
  j1 = si*ix*ix - co*ix;

  if (fabs(tsb[1])<fabs(tsb[0])) {
    sf = j0/tsb[0];
  } else {
    sf = j1/tsb[1];
  }

  /* tsb to sb */
  for (i=0; i<=l+1; i++) {
    sb[i] = tsb[i]*sf;
  }

  /* derivative of sb */
  dsb[0] = co*ix - si*ix*ix;
  for (i=1; i<=l; i++) {
    dsb[i] = ((double)i*sb[i-1] - (double)(i+1.0)*sb[i+1])
              /(2.0*(double)i + 1.0);
  }
}




static void sph_bessel_arec(
  double  x,
  int     l,
  double *sb,
  double *dsb
)
{
  int i;
  double si = sin(x), co = cos(x);

  sb[0]  = si/x;
  dsb[0] = co/x-si/x/x;
  if (l==0) { return; }

  sb[1] = -dsb[0];
  dsb[1] = sb[0] - 2.0*sb[1]/x;
  if (l==1) { return; }

  for (i=2; i<=l; i++) {
    sb[i]  = (double)(2*i-1)*sb[i-1]/x-sb[i-2];
    dsb[i] = sb[i-1]-(double)(i+1)*sb[i]/x;
  }
} /* end if */




/*----------------------------------------------------------------------
  ERI_Spherical_Harmonics
----------------------------------------------------------------------*/
void ERI_Spherical_Harmonics(
  int    l,       /* (IN)  azimuthal quantum number */
  int    m,       /* (IN)  magnetuc quantum number */
  double theta,   /* (IN)  argument */
  double phi,     /* (IN)  argument */
  double Y[2],    /* (OUT) function value (complex) */
  double dYdt[2], /* (OUT) derivative w.r.t. theta (complex) */
  double dYdp[2]  /* (OUT) derivative w.r.t. phi (complex) */
)
{
  unsigned int i, am;
  long double f0, f1;
  double norm, co, si, P, dP;

  am = abs(m);

  /* (l-|m|)! */
  f0 = 1.0;
  for (i=2; i<=(l-am); i++) { f0 *= (long double)i; }
  
  /* (l+|m|)! */
  f1 = 1.0;
  for (i=2; i<=(l+am); i++) { f1 *= (long double)i; }

  /* sqrt((2*l+1)/(4*PI)*(l-|m|)!/(l+|m|)!) */
  norm = sqrt((2.0*(double)l+1.0)/(4.0*PI)*f0/f1);
 
  /* i^(m-|m|) */ 
  if (m<0 && am%2==1) { norm = -norm; }

  /* P_l^|m| */
  ERI_Associated_Legendre_Polynomial(l, am, cos(theta), &P, &dP);

  /* exp(i*m*phi) */
  co = cos((double)m*phi);
  si = sin((double)m*phi);

  Y[0]   = norm*P*co;
  Y[1]   = norm*P*si;
  dYdt[0] = norm*dP*co;
  dYdt[1] = norm*dP*si;
  dYdp[0] = -(double)m*Y[1];
  dYdp[1] =  (double)m*Y[0];
}




/*----------------------------------------------------------------------
  ERI_Real_Spherical_Harmonics

  Real-valued spherical harmonics which is defined as:

  Z(l, m; theta, phi) = a(m)*Y(l, m; theta, phi)
                          + b(m)*Y(l, -m; theta, phi)
 
  where a(m) and b(m) are defined by real_coeff.
----------------------------------------------------------------------*/
void ERI_Real_Spherical_Harmonics(
  int    l,     /* (IN)  azimuthal quantum number */
  int    m,     /* (IN)  magnetuc quantum number */
  double theta, /* (IN)  argument */
  double phi,   /* (IN)  argument */
  double *Z,    /* (OUT) function value (optional) */
  double *dZdt, /* (OUT) derivative w.r.t. theta (optional)*/
  double *dZdp  /* (OUT) derivative w.r.t. phi (optional) */
)
{
  double Y1[2], dY1dt[2], dY1dp[2];
  double Y2[2], dY2dt[2], dY2dp[2];
  double a[2], b[2];

  ERI_Spherical_Harmonics(l,  m, theta, phi, Y1, dY1dt, dY1dp);
  ERI_Spherical_Harmonics(l, -m, theta, phi, Y2, dY2dt, dY2dp);

  ERI_RCSH_Coeff(m, a, b);
 
  if (Z) { 
    *Z    = a[0]*Y1[0]    - a[1]*Y1[1]    + b[0]*Y2[0]    - b[1]*Y2[1];
  }
  if (dZdt) { 
    *dZdt = a[0]*dY1dt[0] - a[1]*dY1dt[1] + b[0]*dY2dt[0] - b[1]*dY2dt[1]; 
  }
  if (dZdp) { 
    *dZdp = a[0]*dY1dp[0] - a[1]*dY1dp[1] + b[0]*dY2dp[0] - b[1]*dY2dp[1]; 
  }
}




/*----------------------------------------------------------------------
  ERI_Associated_Legendre_Polynomial 

  Associated Legendre Polynomial and its derivative.
----------------------------------------------------------------------*/
void ERI_Associated_Legendre_Polynomial(
  int     l, /* (IN)  positive integer */
  int     m, /* (IN)  integer 0 < m < l */
  double  x, /* (IN)  argument -1 <= x <= 1 */
  double *P, /* (OUT) function value */ 
  double *dP /* (OUT) derivative */ 
)
{
  const double cut0 = 1.0e-24;
  const double cut1 = 1.0e-12;
  double Pm, Pm1, f, p0, p1, tmp; 
  int i, ll;
   
  assert( fabs(x) <= 1.0 );
  assert( l >= 0 && m >= 0 && m <= l );
 
  if ((1.0-cut0)<fabs(x)){
    /* x = sgn(x)*(1.0-cut0); */
    if (x < 0.0) {
      x = -(1.0-cut0); 
    } else {
      x = (1.0-cut0); 
    }
  }

  /* calculate Pm */

  /* (-1)^m (2m-1)!! (1-x^2)^(m/2) */
  Pm = 1.0; 
  if (m>0){
    f = 1.0;
    tmp = -sqrt((1.0-x)*(1.0+x));
    for (i=1; i<=m; i++){
      Pm *= f*tmp;
      f += 2.0;
    }
  }
    
  if (l==m){
    /* m = l */
    p0 = Pm;
    p1 = 0.0;
  } else{
    /* m < l */
    Pm1 = x*(2.0*(double)m + 1.0)*Pm;
    if (m<l-1) {
      /* m < l-1 */
      for (ll=m+2; ll<=l; ll++) {
        tmp = (x*(2.0*(double)ll-1.0)*Pm1 - 
                ((double)ll+(double)m-1.0)*Pm)/(double)(ll-m);
        Pm  = Pm1;
        Pm1 = tmp;
      }
    } /* end if */
    p0 = Pm1;
    p1 = Pm;
  } /* end if */

  *P = p0;

  *dP = 0.0;
  tmp = sqrt(1.0-x*x);
  if (cut1<tmp) { 
    *dP = ((double)l*x*p0-(double)(l+m)*p1)/tmp;
  }
}




/*----------------------------------------------------------------------
  ERI_Associated_Laguerre_Polynomial

  The associated Laguerre polynomials in the Rodridues representation.
----------------------------------------------------------------------*/
void ERI_Associated_Laguerre_Polynomial(
  double  x,     /* (IN)  argument */
  int     n,     /* (IN)  zero or positive integer */
  double  alpha, /* (IN)  */
  double *L,     /* (OUT) function value */
  double *dL     /* (OUT) derivative */
)
{
  int i;
  double L0, L1, L2;

  assert( n >= 0 );
 
  if (n==0) {
    *L  = 1.0;
    *dL = 0.0;
    return;
  } else if (n==1) {
    *L  = (1.0 + alpha)-x;
    *dL = -1.0;
    return;
  } 

  /* n > 1 */
  L1 = 1.0;           /* L_0 */
  L2 = (1.0+alpha)-x; /* L_1 */

  /* recurrence */
  for (i=1; i<n; i++) {
    L0 = L1; /* L_{i-1} */
    L1 = L2; /* L_i */
    L2 = (((double)(2*i+1)+(alpha-x))*L1 
            - ((double)i+alpha)*L0)/(double)(i+1); /* L_{i+1} */
  }
  *L = L2;
  *dL = ((double)n*L2-((double)n+alpha)*L1)/x;
}




/*----------------------------------------------------------------------
  ERI_Gaunt

  Routine to calculate the Gaunt coefficients.
  
  See Eq. (3.7.73) in Modern Quantum Mechanics by J. J. Sakurai.
----------------------------------------------------------------------*/
double ERI_Gaunt(int l,  int m, int l1, int m1, int l2, int m2)
{
  double tmp0, tmp1, tmp2, tmp3;
  double result, cg1, cg2;

  cg1 = CG(l1, 0,  l2,  0, l, 0);
  cg2 = CG(l1, m1, l2, m2, l, m);

  return cg1*cg2*sqrt( (double)((2*l1+1)*(2*l2+1))/(double)(8*l+4)/PI );
}


/*----------------------------------------------------------------------
  ERI_Gaunt_R

  Gaunt coefficients for real spherical harmonics.

  GC_R(L,L1,L2) = int Z(L) Z(L1) Z(L2)
----------------------------------------------------------------------*/
double ERI_Gaunt_R(int l,  int m, int l1, int m1, int l2, int m2)
{
  double a[2], b[2], a1[2], b1[2], a2[2], b2[2];
  double aaa, aab, aba, abb, baa, bab, bba, bbb, s;
 
  ERI_RCSH_Coeff( m,  a,  b);
  ERI_RCSH_Coeff(m1, a1, b1);
  ERI_RCSH_Coeff(m2, a2, b2);

  s = (0==m%2) ? 1.0 : -1.0;

  aaa = (a[0]*a1[0]-a[1]*a1[1])*a2[0] - (a[0]*a1[1]+a[1]*a1[0])*a2[1];
  aab = (a[0]*a1[0]-a[1]*a1[1])*b2[0] - (a[0]*a1[1]+a[1]*a1[0])*b2[1];
  aba = (a[0]*b1[0]-a[1]*b1[1])*a2[0] - (a[0]*b1[1]+a[1]*b1[0])*a2[1];
  abb = (a[0]*b1[0]-a[1]*b1[1])*b2[0] - (a[0]*b1[1]+a[1]*b1[0])*b2[1];
  baa = (b[0]*a1[0]-b[1]*a1[1])*a2[0] - (b[0]*a1[1]+b[1]*a1[0])*a2[1];
  bab = (b[0]*a1[0]-b[1]*a1[1])*b2[0] - (b[0]*a1[1]+b[1]*a1[0])*b2[1];
  bba = (b[0]*b1[0]-b[1]*b1[1])*a2[0] - (b[0]*b1[1]+b[1]*b1[0])*a2[1];
  bbb = (b[0]*b1[0]-b[1]*b1[1])*b2[0] - (b[0]*b1[1]+b[1]*b1[0])*b2[1];

  return s*(   aaa*ERI_Gaunt(l, -m, l1,  m1, l2,  m2)
             + aab*ERI_Gaunt(l, -m, l1,  m1, l2, -m2)
             + aba*ERI_Gaunt(l, -m, l1, -m1, l2,  m2)
             + abb*ERI_Gaunt(l, -m, l1, -m1, l2, -m2)
             + baa*ERI_Gaunt(l,  m, l1,  m1, l2,  m2)
             + bab*ERI_Gaunt(l,  m, l1,  m1, l2, -m2)
             + bba*ERI_Gaunt(l,  m, l1, -m1, l2,  m2)
             + bbb*ERI_Gaunt(l,  m, l1, -m1, l2, -m2) );
}





/*----------------------------------------------------------------------
  Clebsch-Gordan coefficients

  <j1 j2; m1 m2 | j1 j2; j m>
  =
  delta(m,m1+m2)*sqrt(2j+1)
    *sqrt{ (j1+j2-j)! (j+j1-j2)! (j+j2-j1)! / (j+j1+j2+1)! } 
    *sqrt{ (j1+m1)! (j1-m1)! (j2+m2)! (j2-m2)! (j+m)! (j-m!) }
    *sum_k (-1)^k / ( k! (j1+j2-j-k)! (j1-m1-k)! (j2+m2-k)! 
                       (j-j2+m1+k)!  (j-j1-m2+k)! )

  See M. Avbramowitz and I.A. Stegun (eds.), "Handbook of Mathematical 
  Functions with Formulas, Graphs, and Tables," 
  Dover Publications, Inc. (Mineola, NY), 1972; an electronic copy of 
  this book is available at http://www.math.sfu.ca/~cbm/aands/.
----------------------------------------------------------------------*/
static double CG(int j1, int m1, int j2, int m2, int j, int m)
{
  int kmin, kmax, k;
  double sgn, cg;
  const int fmax = 3*CG_JMAX;
  double f[3*CG_JMAX];

  /* this routine safely calculates when j1, j2, j < CG_JMAX */
  if (j1>CG_JMAX || j2>CG_JMAX || j>CG_JMAX) {
    fprintf(stderr, "*** out of bound (%s, %d)\n", __FILE__, __LINE__);
    abort();
  }

  /* factorials */
  f[0] = f[1] = 1.0;
  for (k=2; k<fmax; k++) { f[k] = f[k-1]*(double)k; }

  /* delta(m,m1+m2) */
  if (m != m1+m2) return 0.0;

  /* conditions for j1, j2, and j */
  if (j1+j2 < j || j2+j < j1 || j+j1 < j2) return 0.0;

  /* conditions for m1, m2 and m */
  if (abs(m1)>j1 || abs(m2)>j2 || abs(m)>j) return 0.0;

  /* determin the range of k 
     max(0, -j+j2-m1, -j+j1+m2) <= k <= min(j1+j2-j, j1-m1, j2+m2) */
  kmin = 0;
  if (kmin < -j+j2-m1) kmin = -j+j2-m1;
  if (kmin < -j+j1+m2) kmin = -j+j1+m2;

  kmax = j1+j2-j;
  if (kmax > j1-m1) kmax = j1-m1;
  if (kmax > j2+m2) kmax = j2+m2;

  if (kmin>kmax) return 0.0;

  cg = 0.0;
  sgn = kmin%2==0 ? 1.0 : -1.0;
  for (k=kmin; k<=kmax; k++) {
    cg += sgn/(f[k]*f[j1+j2-j-k]*f[j1-m1-k]*f[j2+m2-k]
          *f[j-j2+m1+k]*f[j-j1-m2+k]);
    sgn = -sgn;
  }
  
  cg *= sqrt( (double)(2*j+1)*f[j1+j2-j]*f[j+j1-j2]*f[j+j2-j1]
    *f[j1+m1]*f[j1-m1]*f[j2+m2]*f[j2-m2]*f[j+m]*f[j-m]
    /f[j+j1+j2+1] );

  return cg;
}



void ERI_GLQ(double *x, double *w, int n)
{
  int i, itr;
  double err, xi, i1, L, dL;

  for (i=0; i<n; i++) {

    /* Initial guess for the roots of Laguerre polynomials 
       See A. Stroud and D. Secrest, "Gaussian Quadrature
       Formulas," Prentice-Hall, Englewood, NJ, 1966. */
    if (i==0) {
      xi = 3.0/(1.0+2.4*(double)n);
    } else if (i==1) {
      xi = x[i-1] + 15.0/(1.0+2.5*(double)n);
    } else {
      i1 = (double)(i-1);
      xi = x[i-1] + (x[i-1]-x[i-2])*(1.0+2.55*i1)/(1.9*i1);
    }

    /* Find the roots by Newton's method */ 
    err = DBL_MAX;
    for (itr=0; (itr<MAX_ITR) && (fabs(err)>ERR_MAX); itr++) {
      ERI_Associated_Laguerre_Polynomial(xi, n, 0.0, &L, &dL);
      err = L/dL;
      xi -= err;
    }

#ifndef NDEBUG
    if (fabs(err)>ERR_MAX) {
      fprintf(stderr,  "*** error in ERI_GLQ: no convergence\n");
      fprintf(stderr,  "    err = %12.4e\n", err);
    }
#endif

    x[i] = xi;
    w[i] = 1.0/dL/dL/xi;
  } /* end loop of i */
}


