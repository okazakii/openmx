/*----------------------------------------------------------------------
  eri_fsbt.h

  Routine to perform fast spherical Bessel transform (FSBT).
  This also manages the logarithmic radial mesh information.
 
  Reference:
 
  A. E. Siegman, "Quasi fast Hankel transform", Opt. Lett., vol. 1,
  no. 1, pp.13 (1977). 
  
  J. D. Talman, "Numerical Fourier and Bessel Transforms in Logarithmic
  Variables", J. Comp. Phys., vol. 29, pp.35 (1978).
 
  Jul. 2008, M. Toyoda
----------------------------------------------------------------------*/
#ifndef LIBERI_ERI_FSBT_H_INCLUDED
#define LIBERI_ERI_FSBT_H_INCLUDED

typedef struct ERI_FSBT_Struct ERI_FSBT_t;


ERI_FSBT_t* ERI_FSBT_Init(
  int    lmax,
  int    ngrid,
  double rho0, 
  double drho
);


void ERI_FSBT_Free(ERI_FSBT_t *ptr);


size_t ERI_FSBT_Required_Size(int lmax, int ngrid);


int           ERI_FSBT_ngrid(const ERI_FSBT_t *ptr);
int           ERI_FSBT_lmax(const ERI_FSBT_t *ptr);
double        ERI_FSBT_Mesh_r(const ERI_FSBT_t *ptr, int i);
double        ERI_FSBT_Mesh_k(const ERI_FSBT_t *ptr, int i);
const double* ERI_FSBT_Mesh_Array_r(const ERI_FSBT_t *ptr);
const double* ERI_FSBT_Mesh_Array_k(const ERI_FSBT_t *ptr);
double        ERI_FSBT_Mesh_dr(const ERI_FSBT_t *ptr, int i);
double        ERI_FSBT_Mesh_dk(const ERI_FSBT_t *pre, int i);


void ERI_FSBT_Transform(
  ERI_FSBT_t   *ptr,
  double       *out,
  const double *in,
  int           l,
  int           m
); 


void ERI_FSBT_Transform_Input(
  ERI_FSBT_t   *ptr,
  const double *in,
  int           m
); 


void ERI_FSBT_Transform_Output(
  ERI_FSBT_t *ptr,
  double     *out,
  int        l
); 

#if 0
void ERI_FSBT_Transform_Real(
  ERI_FSBT_t   *ptr,
  double       *out,
  const double *in,
  int           l,
  int           m
); 


void ERI_FSBT_Transform_Input_Real(
  ERI_FSBT_t   *ptr,
  const double *in,
  int           m
); 


void ERI_FSBT_Transform_Output_Real(
  ERI_FSBT_t *ptr,
  double     *out,
  int        l
); 
#endif

#endif

/* EOF */
