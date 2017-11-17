/*----------------------------------------------------------------------
  eri_sf.h

  Special functions used in LIBERI.

  16 Feb. 2009, M. Toyoda
----------------------------------------------------------------------*/
#ifndef LIBERI_ERI_SF_H_INCLUDED
#define LIBERI_ERI_SF_H_INCLUDED

#if 0
void ERI_RCSH_Coeff(int m, double a[2], double b[2]);
void ERI_RCSH_Coeff_Inverse(int m, double c[2], double d[2]);
#endif

void ERI_Spherical_Bessel(
  double  x,  /* (IN) argument */
  int     l,  /* (IN) order */
  double *sb, /* (OUT) function values of order from 0 to l */
  double *dsb /* (OUT) derivatives of order from 0 to l */
);
 

void ERI_Spherical_Harmonics(
  int    l,       /* (IN)  azimuthal quantum number */
  int    m,       /* (IN)  magnetuc quantum number */
  double theta,   /* (IN)  argument */
  double phi,     /* (IN)  argument */
  double Y[2],    /* (OUT) function value (complex) */
  double dYdt[2], /* (OUT) derivative w.r.t. theta (complex) */
  double dYdp[2]  /* (OUT) derivative w.r.t. phi (complex) */
);


void ERI_Real_Spherical_Harmonics(
  int    l,     /* (IN)  azimuthal quantum number */
  int    m,     /* (IN)  magnetuc quantum number */
  double theta, /* (IN)  argument */
  double phi,   /* (IN)  argument */
  double *Z,    /* (OUT) function value */
  double *dZdt, /* (OUT) derivative w.r.t. theta */
  double *dZdp  /* (OUT) derivative w.r.t. phi */
);


void ERI_Associated_Legendre_Polynomial(
  int     l, /* (IN)  positive integer */
  int     m, /* (IN)  integer 0 < m < l */
  double  x, /* (IN)  argument -1 <= x <= 1 */
  double *P, /* (OUT) function value */ 
  double *dP /* (OUT) derivative */ 
);


void ERI_Associated_Laguerre_Polynomial(
  double  x,     /* (IN)  argument */
  int     n,     /* (IN)  zero or positive integer */
  double  alpha, /* (IN)  */
  double *L,     /* (OUT) function value */
  double *dL     /* (OUT) derivative */
);


double ERI_Gaunt(
  int l,  
  int m, 
  int l1, 
  int m1, 
  int l2, 
  int m2
);


double ERI_Gaunt_R(
  int l,
  int m,
  int l1,
  int m1,
  int l2,
  int m2
);


void ERI_GLQ(
  double *x, 
  double *w, 
  int     n
);


#endif /* LIBERI_ERI_SF_H_INCLUDED */
