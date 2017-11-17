/*----------------------------------------------------------------------
  eri_sbt.h
----------------------------------------------------------------------*/
#ifndef LIBERI_ERI_SBT_H_INCLUDED
#define LIBERI_ERI_SBT_H_INCLUDED

#define ERI_SBT_FORWARD  1
#define ERI_SBT_BACKWARD 2


/* prototypes */
typedef void          (*fnc_free)(void *);
typedef int           (*fnc_ngrid)(void* );
typedef int           (*fnc_lmax)(void* );
typedef double        (*fnc_mesh)(void*, int );
typedef const double* (*fnc_array)(void* );
typedef double        (*fnc_mesh_dr)(void*, int );
typedef void          (*fnc_transform_in)(void*, const double*, int);
typedef void          (*fnc_transform_out)(void*, double*, int);
typedef void          (*fnc_transform)(void*, double*, const double*, int, int);


typedef struct ERI_SBT_struct ERI_SBT_t;


ERI_SBT_t* ERI_SBT_Init(
  fnc_free          free,
  fnc_ngrid         ngrid,
  fnc_lmax          lmax,
  fnc_mesh          mesh_r,
  fnc_mesh          mesh_k,
  fnc_array         mesh_array_r,
  fnc_array         mesh_array_k,
  fnc_mesh_dr       mesh_dr,
  fnc_mesh_dr       mesh_dk,
  fnc_transform_in  transform_in,
  fnc_transform_out transform_out,
  fnc_transform     transform,
  void              *instance
);


size_t ERI_SBT_Required_Size(void);

void ERI_SBT_Free(ERI_SBT_t *ptr);

int           ERI_SBT_ngrid(ERI_SBT_t *ptr);
int           ERI_SBT_lmax(ERI_SBT_t *ptr);
double        ERI_SBT_Mesh_r(ERI_SBT_t *ptr, int i);
double        ERI_SBT_Mesh_k(ERI_SBT_t *ptr, int i);
const double* ERI_SBT_Mesh_Array_r(ERI_SBT_t *ptr);
const double* ERI_SBT_Mesh_Array_k(ERI_SBT_t *ptr);
double        ERI_SBT_Mesh_dr(ERI_SBT_t *ptr, int i);
double        ERI_SBT_Mesh_dk(ERI_SBT_t *ptr, int i);

void ERI_SBT_Transform(
  ERI_SBT_t    *ptr, 
  double       *out, 
  const double *in, 
  int l, 
  int direction
);

void ERI_SBT_Transform_Input(
  ERI_SBT_t *ptr, 
  const double *in, 
  int direction
);

void ERI_SBT_Transform_Output(
  ERI_SBT_t *ptr, 
  double *out, 
  int l
);

#if 0
void ERI_SBT_Transform_Real(
  ERI_SBT_t *ptr,
  double *out,
  const double *in,
  int l,
  int direction
);

void ERI_SBT_Transform_Input_Real(
  ERI_SBT_t *ptr, 
  const double *in,
  int direction
);

void ERI_SBT_Transform_Output_Real(
  ERI_SBT_t *ptr, 
  double *out, 
  int l
);
#endif
#endif
