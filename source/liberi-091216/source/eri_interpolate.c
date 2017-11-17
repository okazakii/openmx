/*----------------------------------------------------------------------
  eri_interpolate.c

  Coded by M. Toyoda, May 2009, JAIST/RCIS.

  * This uses CLAPACK, so that you must link with appropreate libraries
    by specifing, for example, the following link options:

    cc *.o -llapack -lblas -lf2c (or, -lg2c for gcc).
----------------------------------------------------------------------*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "eri_interpolate.h"


/* LAPACK driver to solve the linear eq. A*x = b, 
   where A is an N-by-N tri-diagonal matrix */
extern long dgtsv_(long*, long*, double*, double*, double*, double*, long*, long*);



/*----------------------------------------------------------------------
                     ***  LINEAR INTERPOLATION ***
----------------------------------------------------------------------*/
double ERI_Interpolate_Linear_Real(
  double        x, 
  const double *sx,
  const double *sy,
  int           ns
)
{
  int i;
  double x0, x1, y0, y1;

  assert( 0 < ns );

  /* find the smallest one of the mesh points which are larger than x */
  for (i=0; i<ns && sx[i]<x; i++)
    ;

  /* boundary check */
  if (i==0)  { return sy[0];    }
  if (i==ns) { return sy[ns-1]; }
  
  /* linear interpolatoion between i-th and (i-1)th mesh points */
  x0 = sx[i-1];
  x1 = sx[i];
  y0 = sy[i-1];
  y1 = sy[i];

  return y0 + (y1-y0)*(x-x0)/(x1-x0);
}




void ERI_Interpolate_Linear_Complex(
  double        out[2],
  double        x, 
  const double *sx,
  const double *sy,
  int           ns
)
{
  int i;
  double x0, x1, y0, y1;

  assert( 0 < ns );

  /* find the smallest one of the mesh points which are larger than x */
  for (i=0; i<ns && sx[i]<x; i++)
    ;

  /* boundary check */
  if (i==0) { 
    out[0] = sy[0]; 
    out[1] = sy[1]; 
    return;
  }
  if (i==ns) { 
    out[0] = sy[2*ns-2]; 
    out[1] = sy[2*ns-1];
    return;
  }

  /* linear interpoation between i-th and (i-1)th mesh points */
  x0 = sx[i-1];
  x1 = sx[i];
  y0 = sy[2*i-2];
  y1 = sy[2*i+0];
  out[0] = y0 + (y1-y0)*(x-x0)/(x1-x0);
  y0 = sy[2*i-1];
  y1 = sy[2*i+1];
  out[1] = y0 + (y1-y0)*(x-x0)/(x1-x0);
}


/*----------------------------------------------------------------------
  CUBIC SPLINE
----------------------------------------------------------------------*/
static double cspline_eval_core(
  double x,
  double x0, 
  double y0, 
  double u0,
  double x1, 
  double y1, 
  double u1
)
{
  double h, w, a, b, c, d, t;

  h = y1-y0;
  w = x1-x0;

  a = (u1-u0)/w/6.0; 
  b = 0.5*u0;
  c = h/w - (u1+2.0*u0)*w/6.0;
  d = y0;

  t = x - x0;

  return d +(c +( b + a*t)*t)*t;
}




struct ERI_CSpline_Struct {
  double *x;
  double *y;
  double *u;
  int     n;
};




ERI_CSpline_t* ERI_CSpline_Init(
  const double *sx,
  const double *sy,
  int           ns
)
{
  long i, n, nrhs, ldb, info;
  double *D, *DL, *DU, *b;
  ERI_CSpline_t *ptr;
 
  assert( 2 < ns );

  ptr    = (ERI_CSpline_t*)malloc(sizeof(ERI_CSpline_t));
  ptr->x = (double*)malloc(sizeof(double)*ns);
  ptr->y = (double*)malloc(sizeof(double)*ns);
  ptr->u = (double*)malloc(sizeof(double)*ns);
  ptr->n = ns;
 
  /* solve the linear equation */ 
  {
    n    = ns-2;
    nrhs = 1;
    ldb  = n;

    D  = (double*)malloc(sizeof(double)*n);
    DL = (double*)malloc(sizeof(double)*n);
    DU = (double*)malloc(sizeof(double)*n);
    b  = (double*)malloc(sizeof(double)*n);

    for (i=0; i<n; i++) {
      D[i]  = 2.0*(sx[i+2]-sx[i]);
      DL[i] = DU[i] = sx[i+2]-sx[i+1];
      b[i]  = 6.0*((sy[i+2]-sy[i+1])/(sx[i+2]-sx[i+1])
                    -(sy[i+1]-sy[i])/(sx[i+1]-sx[i]));
    }
    dgtsv_(&n, &nrhs, DL, D, DU, b, &ldb, &info);

    for (i=0; i<n; i++) { ptr->u[i+1] = b[i]; } 
    ptr->u[0] = 0.0;
    ptr->u[n+1] = 0.0;

    free(D); 
    free(DL); 
    free(DU); 
    free(b); 
  }
  
  /* copy samples */ 
  for (i=0; i<ns; i++) {
    ptr->x[i] = sx[i];
    ptr->y[i] = sy[i];
  }

  return ptr;
}






ERI_CSpline_t* ERI_CSpline_Init_Complex(
  const double *sx,
  const double *sy,
  int           ns
)
{
  long i, n, nrhs, ldb, info;
  double *D, *DL, *DU, *b;
  ERI_CSpline_t *ptr;
 
  assert( 2 < ns );

  ptr = (ERI_CSpline_t*)malloc(sizeof(ERI_CSpline_t));
  ptr->x = (double*)malloc(sizeof(double)*ns);
  ptr->y = (double*)malloc(sizeof(double)*ns*2);
  ptr->u = (double*)malloc(sizeof(double)*ns*2);
  ptr->n = ns;
 
  /* solve the linear equation */ 
  {
    n    = ns-2;
    nrhs = 1;
    ldb  = n;

    D  = (double*)malloc(sizeof(double)*n);
    DL = (double*)malloc(sizeof(double)*n);
    DU = (double*)malloc(sizeof(double)*n);
    b  = (double*)malloc(sizeof(double)*n);
    
    /* real part */
    for (i=0; i<n; i++) {
      D[i]  = 2.0*(sx[i+2]-sx[i]);
      DL[i] = DU[i] = sx[i+2]-sx[i+1];
      b[i]  = 6.0*((sy[2*i+4]-sy[2*i+2])/(sx[i+2]-sx[i+1])
                   -(sy[2*i+2]-sy[2*i])/(sx[i+1]-sx[i]));
    }
    dgtsv_(&n, &nrhs, DL, D, DU, b, &ldb, &info);
    for (i=0; i<n; i++) { ptr->u[2*i+2] = b[i]; } 

    /* imag part */
    for (i=0; i<n; i++) {
      D[i]  = 2.0*(sx[i+2]-sx[i]);
      DL[i] = DU[i] = sx[i+2]-sx[i+1];
      b[i]  = 6.0*((sy[2*i+5]-sy[2*i+3])/(sx[i+2]-sx[i+1])
                 -(sy[2*i+3]-sy[2*i+1])/(sx[i+1]-sx[i]));
    }
    dgtsv_(&n, &nrhs, DL, D, DU, b, &ldb, &info);
    for (i=0; i<n; i++) { ptr->u[2*i+3] = b[i]; } 

    ptr->u[0] = 0.0;
    ptr->u[1] = 0.0;
    ptr->u[2*n+2] = 0.0;
    ptr->u[2*n+3] = 0.0;

    free(D); 
    free(DL); 
    free(DU); 
    free(b); 
  }
  
  /* copy samples */ 
  for (i=0; i<ns; i++) {
    ptr->x[i] = sx[i];
    ptr->y[2*i]   = sy[2*i];
    ptr->y[2*i+1] = sy[2*i+1];
  }

  return ptr;
}





void ERI_CSpline_Free(ERI_CSpline_t *ptr)
{
  if (ptr) {
    if (ptr->x) { free(ptr->x); }
    if (ptr->y) { free(ptr->y); }
    if (ptr->u) { free(ptr->u); }
    free(ptr);
  }
}


static double extrapolate(
  double x,
  double x0,
  double x1,
  double y0,
  double y1,
  double u0,
  double u1
)
{
  double c;

  c = (y1-y0)/(x1-x0)*(x1-x0)*(2.0*u0+u1)/6.0;

  return y0 + x*(x-x0); 
}



double ERI_CSpline_Eval(
  double               x,
  const ERI_CSpline_t *ptr
)
{
  int i;

  const double *sx = ptr->x;
  const double *sy = ptr->y;
  const double *u  = ptr->u;
  int           ns = ptr->n;

  /* find the segment */
  for (i=0; sx[i]<x && i<ns; i++) 
    ;
 
  if (0==i) { 
    return extrapolate(x, sx[0], sx[1], sy[0], sy[1], u[0], u[1]);
  }
  if (ns==i) { return sy[ns-1]; }

  return cspline_eval_core(x, sx[i-1], sy[i-1], u[i-1], sx[i], sy[i], u[i]);
}




void ERI_CSpline_Eval_Complex(
  double               y[2],
  double               x,
  const ERI_CSpline_t *ptr
)
{
  int i;

  const double *sx = ptr->x;
  const double *sy = ptr->y;
  const double *u  = ptr->u;
  int           ns = ptr->n;

  /* find the segment */
  for (i=0; sx[i]<x && i<ns; i++) 
    ;
 
  if (0==i) { 
    y[0] = extrapolate(x, sx[0], sx[1], sy[0], sy[2], u[0], u[2]);
    y[1] = extrapolate(x, sx[0], sx[1], sy[1], sy[3], u[1], u[3]);
    return;
  }
  if (ns==i) { 
    y[0] = sy[2*ns-2];
    y[1] = sy[2*ns-1];
    return;
  }

  y[0] = cspline_eval_core(x, sx[i-1], sy[2*i-2], u[2*i-2], sx[i], sy[2*i], u[2*i]);
  y[1] = cspline_eval_core(x, sx[i-1], sy[2*i-1], u[2*i-1], sx[i], sy[2*i+1], u[2*i+1]);
}


/* EOF */
