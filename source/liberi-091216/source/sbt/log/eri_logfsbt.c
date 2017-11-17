/*----------------------------------------------------------------------
  eri_logfsbt.c
 
  Oct. 2008, M. Toyoda
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "eri_def.h"
#include "eri_logfsbt.h"
#include "eri_fsbt.h"
#include "eri_sf.h"

#define LOCAL_LMAX 100

    
static double g_tmp_work[ERI_NGRIDMAX*2];
static double g_tmp_fq[ERI_NQMAX*2];
static double g_tmp_dfq[ERI_NQMAX*2];


struct ERI_LogFSBT_struct {
  double     *xq;
  double     *jqr;
  double     *djqr;
  int         nq;
  ERI_FSBT_t *fsbt;
  double     *saved_in; 
};

typedef struct ERI_LogFSBT_struct ERI_LogFSBT_t;




/*----------------------------------------------------------------------
                                                      internal functions
----------------------------------------------------------------------*/
static double interpolation(
  double x,
  double x0, 
  double y0, 
  double dy0, 
  double x1,
  double y1, 
  double dy1
)
{
#if 1
  return y0 + (y1-y0)*(x-x0)/(x1-x0);
#else
  double p, q, b, c, d;
    
  p = x1-x0;
  q = y1-y0;
 
  b = dy0*p; 
  c = 3.0*q-(dy1+2.0*dy0)*p;
  d = (dy1+dy0)*p-2.0*q;

  x  = (x-x0)/p;

  return y0+(b+(c+d*x)*x)*x;
#endif
}


static void Integrate_Trapezoidal2(
  double        out[2],
  const double *in,
  const double *xmesh,
  int           ngrid
)
{
  int i;
  double dx;

  out[0] = 0.0;
  out[1] = 0.0;
  for (i=1; i<ngrid; i++) {
    dx = xmesh[i]-xmesh[i-1];
    out[0] += (in[2*i+0]+in[2*i-2])*dx;
    out[1] += (in[2*i+1]+in[2*i-1])*dx;
  }
  out[0] *= 0.5;
  out[1] *= 0.5;
}


static void correction(
  double       *out,
  const double *in,
  const double *mesh_r,
  const double *mesh_k, 
  const double *xq,
  const double *jqr,
  const double *djqr,
  int           nq,
  int           ngrid,
  int           lmax,
  int           l
)
{
  int iq, i, j, j0, j1; 
  double q, q0, q1, k, r, jq, djq;
  
  double *tmp_in  = g_tmp_work;
  double *tmp_fq  = g_tmp_fq;
  double *tmp_dfq = g_tmp_dfq;
 
  /* integrate */
  for (iq=0; iq<nq; iq++) {
    tmp_fq[2*iq+0] = 0.0;
    tmp_fq[2*iq+1] = 0.0;
    tmp_dfq[2*iq+0] = 0.0;
    tmp_dfq[2*iq+1] = 0.0;
 
    if (fabs(xq[iq])<1e-10) {
      if (l==0) {
        for (i=0; i<ngrid; i++) {
          r = mesh_r[i];
          tmp_in[2*i+0] = in[2*i+0]*r*r;
          tmp_in[2*i+1] = in[2*i+1]*r*r;
        } 
        Integrate_Trapezoidal2(&tmp_fq[2*iq], tmp_in, mesh_r, ngrid);
      } else if (l==1) {
        for (i=0; i<ngrid; i++) {
          r = mesh_r[i];
          tmp_in[2*i+0] = in[2*i+0]*r*r*r/3.0;
          tmp_in[2*i+1] = in[2*i+1]*r*r*r/3.0;
        } 
        Integrate_Trapezoidal2(&tmp_dfq[2*iq], tmp_in, mesh_r, ngrid);
      } 
    } else {
      for (i=0; i<ngrid; i++) {
        r = mesh_r[i];
        jq = jqr[(iq*lmax+l)*ngrid+i];
        tmp_in[2*i+0] = in[2*i+0]*r*r*jq;
        tmp_in[2*i+1] = in[2*i+1]*r*r*jq;
      } 
      Integrate_Trapezoidal2(&tmp_fq[2*iq], tmp_in, mesh_r, ngrid);
      for (i=0; i<ngrid; i++) {
        r = mesh_r[i];
        djq = djqr[(iq*lmax+l)*ngrid+i];
        tmp_in[2*i+0] = in[2*i+0]*r*r*r*djq;
        tmp_in[2*i+1] = in[2*i+1]*r*r*r*djq;
      } 
      Integrate_Trapezoidal2(&tmp_dfq[2*iq], tmp_in, mesh_r, ngrid);
    }
  }

  /* interporate */
  for (i=0; i<ngrid; i++) {
    k = mesh_k[i];
    j1 = -1;
    j0 = 0;
    for (j=0; j<nq; j++) {
      q = xq[j];
      if (k<q) { j1=j; j0=j-1; j=nq; }
    }
    if (j1==-1) {
     i=ngrid;
    } else {
      if (j0<0) {
        out[2*i+0] = tmp_fq[0];
        out[2*i+1] = tmp_fq[1];
      } else {
        /* interpolate by cubic spline */
        q0 = xq[j0]; 
        q1 = xq[j1];
        out[2*i+0] = interpolation(
          k, q0, tmp_fq[2*j0+0], tmp_dfq[2*j0+0], 
          q1, tmp_fq[2*j1+0], tmp_dfq[2*j1+0]);
        out[2*i+1] = interpolation(
          k, q0, tmp_fq[2*j0+1], tmp_dfq[2*j0+1], 
          q1, tmp_fq[2*j1+1], tmp_dfq[2*j1+1]);
      }
    }
  }         
}  




/*----------------------------------------------------------------------
                                                             interface 
----------------------------------------------------------------------*/
void ERI_LogFSBT_Free(void *p)
{
  ERI_LogFSBT_t *ptr = (ERI_LogFSBT_t*)p;
  if (ptr) {
    if (ptr->xq)   { free(ptr->xq);   }
    if (ptr->jqr)  { free(ptr->jqr);  }
    if (ptr->djqr) { free(ptr->djqr); }
    if (ptr->saved_in) { free(ptr->saved_in); }
    ERI_FSBT_Free(ptr->fsbt);
    free(ptr);
  }
}


double ERI_LogFSBT_Mesh_r(void *p, int i)
{
  ERI_LogFSBT_t *ptr = (ERI_LogFSBT_t*)p;
  return ERI_FSBT_Mesh_r(ptr->fsbt, i); 
}


double ERI_LogFSBT_Mesh_k(void *p, int i)
{
  ERI_LogFSBT_t *ptr = (ERI_LogFSBT_t*)p;
  return ERI_FSBT_Mesh_k(ptr->fsbt, i); 
}


int ERI_LogFSBT_Mesh_nmesh(void *p)
{
  ERI_LogFSBT_t *ptr = (ERI_LogFSBT_t*)p;
  return ERI_FSBT_ngrid(ptr->fsbt); 
}


int ERI_LogFSBT_lmax(void *p)
{
  ERI_LogFSBT_t *ptr = (ERI_LogFSBT_t*)p;
  return ERI_FSBT_lmax(ptr->fsbt);
}


const double* ERI_LogFSBT_Mesh_Array_r(void *p)
{
  ERI_LogFSBT_t *ptr = (ERI_LogFSBT_t*)p;
  return ERI_FSBT_Mesh_Array_r(ptr->fsbt); 
}


const double* ERI_LogFSBT_Mesh_Array_k(void *p)
{
  ERI_LogFSBT_t *ptr = (ERI_LogFSBT_t*)p;
  return ERI_FSBT_Mesh_Array_k(ptr->fsbt); 
}


double ERI_LogFSBT_Mesh_dr(void *p, int i)
{
  ERI_LogFSBT_t *ptr = (ERI_LogFSBT_t*)p;
  return ERI_FSBT_Mesh_dr(ptr->fsbt, i);
}


double ERI_LogFSBT_Mesh_dk(void *p, int i)
{
  ERI_LogFSBT_t *ptr = (ERI_LogFSBT_t*)p;
  return ERI_FSBT_Mesh_dk(ptr->fsbt, i);
}


void ERI_LogFSBT_Transform_Input(void *p, const double *in, int dir)
{
  int i;

  ERI_LogFSBT_t *ptr      = (ERI_LogFSBT_t*)p;
  int            ngrid    = ERI_FSBT_ngrid(ptr->fsbt);
  double        *saved_in = ptr->saved_in;
   
  for (i=0; i<2*ngrid; i++) { saved_in[i] = in[i]; }

  ERI_FSBT_Transform_Input(ptr->fsbt, in, 0);
}


void ERI_LogFSBT_Transform_Output(void *p, double *out, int l)
{
  ERI_LogFSBT_t *ptr   = (ERI_LogFSBT_t*)p;
  int            ngrid = ERI_FSBT_ngrid(ptr->fsbt);
  const double  *rmesh = ERI_FSBT_Mesh_Array_r(ptr->fsbt);
  const double  *kmesh = ERI_FSBT_Mesh_Array_k(ptr->fsbt);
  int            lmax  = ERI_FSBT_lmax(ptr->fsbt);
  const double  *in    = ptr->saved_in;
 
  ERI_FSBT_Transform_Output(ptr->fsbt, out, l);
  correction(out, in, rmesh, kmesh, ptr->xq, ptr->jqr, ptr->djqr, 
    ptr->nq, ngrid, lmax, l);
}


void ERI_LogFSBT_Transform(
  void         *p,
  double       *out, 
  const double *in,
  int           l, 
  int           dir
)
{
  ERI_LogFSBT_t *ptr   = (ERI_LogFSBT_t*)p;
  int            ngrid = ERI_FSBT_ngrid(ptr->fsbt);
  const double  *rmesh = ERI_FSBT_Mesh_Array_r(ptr->fsbt);
  const double  *kmesh = ERI_FSBT_Mesh_Array_k(ptr->fsbt);
  int            lmax  = ERI_FSBT_lmax(ptr->fsbt);

  ERI_FSBT_Transform_Input(ptr->fsbt, in, 0);
  ERI_FSBT_Transform_Output(ptr->fsbt, out, l);
  correction(out, in, rmesh, kmesh, ptr->xq, ptr->jqr, ptr->djqr, 
    ptr->nq, ngrid, lmax, l);
}


ERI_SBT_t* ERI_LogFSBT_Init(
  int    lmax,  /* (IN) maximum number of angular momentum */
  int    ngrid, /* (IN) number of radial logarithm mesh */
  double rho0,  /* (IN) lower bound of rho (radial) mesh */
  double dt,    /* (IN) interval of t-mesh */
  double qmin,
  double qmax,
  int    nq
)
{
  int i, j, l;
  double dq, r, tmp_sb[LOCAL_LMAX], tmp_dsb[LOCAL_LMAX];
  const double *rmesh;
  ERI_LogFSBT_t *pinst;
  ERI_SBT_t *ptr;

  STEPTRACE( "ERI_LogFSBT_Init: in" );
    
  if (lmax >LOCAL_LMAX) {
    fprintf(stderr, "error in liberi: too large lmax\n");
    return NULL;
  }

  {
    pinst = (ERI_LogFSBT_t*)malloc(sizeof(ERI_LogFSBT_t));
    if (NULL==pinst) { return NULL; }

    pinst->fsbt = ERI_FSBT_Init(lmax, ngrid, rho0, dt);
    rmesh = ERI_FSBT_Mesh_Array_r(pinst->fsbt);
  
    pinst->xq   = (double*)malloc(sizeof(double)*nq);
    pinst->jqr  = (double*)malloc(sizeof(double)*nq*ngrid*lmax);
    pinst->djqr = (double*)malloc(sizeof(double)*nq*ngrid*lmax);
    pinst->saved_in = (double*)malloc(sizeof(double)*ngrid*2);
    if (NULL==pinst->fsbt || NULL==pinst->xq || NULL==pinst->jqr 
      || NULL==pinst->djqr || NULL==pinst->saved_in) { 
      ERI_LogFSBT_Free(pinst);
      return NULL; 
    }

    /* points to be corrected */
    dq = (log(qmax)-log(qmin))/(double)(nq-1);
    for (i=1; i<nq; i++) {
      pinst->xq[i] = qmin*exp(dq*(double)i);
    }
    pinst->xq[0] = 0.0;

    /*--- pre-calculate sph. Bessel functions ---*/
    for (i=0; i<nq; i++) {
      for (j=0; j<ngrid; j++) {
        r = rmesh[j];
        ERI_Spherical_Bessel(pinst->xq[i]*r, lmax, tmp_sb, tmp_dsb);
        for (l=0; l<lmax; l++) {
          pinst->jqr[(i*lmax+l)*ngrid+j] = tmp_sb[l];
          pinst->djqr[(i*lmax+l)*ngrid+j] = tmp_dsb[l];
        }
      }
    }
    pinst->nq = nq;
  }

  ptr = ERI_SBT_Init(ERI_LogFSBT_Free, ERI_LogFSBT_Mesh_nmesh,
    ERI_LogFSBT_lmax, ERI_LogFSBT_Mesh_r, ERI_LogFSBT_Mesh_k,
    ERI_LogFSBT_Mesh_Array_r, ERI_LogFSBT_Mesh_Array_k,
    ERI_LogFSBT_Mesh_dr, ERI_LogFSBT_Mesh_dk,
    ERI_LogFSBT_Transform_Input, ERI_LogFSBT_Transform_Output,
    ERI_LogFSBT_Transform,
    pinst
  );

  if (NULL==ptr) {
    ERI_LogFSBT_Free(pinst);
    return NULL;
  }
  
  STEPTRACE( "ERI_LogFSBT_Init: out" );

  return ptr;
}


size_t ERI_LogFSBT_Required_Size(
  int lmax,  /* (IN) maximum number of angular momentum */
  int ngrid, /* (IN) number of radial logarithm mesh */
  int nq
)
{
  return ERI_SBT_Required_Size()
       + sizeof(ERI_LogFSBT_t)
       + ERI_FSBT_Required_Size(lmax, ngrid)
       + sizeof(double)*nq            /* xq */ 
       + sizeof(double)*nq*ngrid*lmax /* jqr */
       + sizeof(double)*nq*ngrid*lmax /* djqr */ 
       + sizeof(double)*ngrid*2       /* saved_in */
  ;
}


/* EOF */
