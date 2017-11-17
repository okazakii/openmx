/*----------------------------------------------------------------------
  eri_sbt.c

  Coded by M. Toyoda, June 2009, JAIST/RCIS.
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include "eri_sbt.h"


struct ERI_SBT_struct {
  /* functions table */
  fnc_free          free;
  fnc_ngrid         ngrid;
  fnc_lmax          lmax;
  fnc_mesh          mesh_r;
  fnc_mesh          mesh_k;
  fnc_array         array_r;
  fnc_array         array_k;
  fnc_mesh_dr       mesh_dr;
  fnc_mesh_dr       mesh_dk;
  fnc_transform_in  transform_in;
  fnc_transform_out transform_out;
  fnc_transform     transform;
  /* instance */
  void              *instance;
};


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
)
{
  ERI_SBT_t *ptr = (ERI_SBT_t*)malloc(sizeof(ERI_SBT_t));
  if (ptr) {
    ptr->free          = free;
    ptr->ngrid         = ngrid;
    ptr->lmax          = lmax;
    ptr->mesh_r        = mesh_r;
    ptr->mesh_k        = mesh_k;
    ptr->array_r       = mesh_array_r;
    ptr->array_k       = mesh_array_k;
    ptr->mesh_dr       = mesh_dr;
    ptr->mesh_dk       = mesh_dk;
    ptr->transform_in  = transform_in;
    ptr->transform_out = transform_out;
    ptr->transform     = transform;
    ptr->instance      = instance;
  }
 
  return ptr;
}


size_t ERI_SBT_Required_Size(void)
{
  return sizeof(ERI_SBT_t);
}


void ERI_SBT_Free(ERI_SBT_t *ptr)
{
  ptr->free(ptr->instance);
  free(ptr);
}


int ERI_SBT_ngrid(ERI_SBT_t *ptr)
{
  return ptr->ngrid(ptr->instance);
}


int ERI_SBT_lmax(ERI_SBT_t *ptr)
{
  return ptr->lmax(ptr->instance);
}


double ERI_SBT_Mesh_r(ERI_SBT_t *ptr, int i)
{
  return ptr->mesh_r(ptr->instance, i);
}


double ERI_SBT_Mesh_k(ERI_SBT_t *ptr, int i)
{
  return ptr->mesh_k(ptr->instance, i);
}


const double* ERI_SBT_Mesh_Array_r(ERI_SBT_t *ptr)
{
  return ptr->array_r(ptr->instance);
}


const double* ERI_SBT_Mesh_Array_k(ERI_SBT_t *ptr)
{
  return ptr->array_k(ptr->instance);
}


double ERI_SBT_Mesh_dr(ERI_SBT_t *ptr, int i)
{
  return ptr->mesh_dr(ptr->instance, i);
}


double ERI_SBT_Mesh_dk(ERI_SBT_t *ptr, int i)
{
  return ptr->mesh_dk(ptr->instance, i);
}


void ERI_SBT_Transform(
  ERI_SBT_t    *ptr,
  double       *out,
  const double *in,
  int           l,
  int           direction
)
{
  if (ptr->transform) {
    ptr->transform(ptr->instance, out, in, l, direction);
  } else {
    ERI_SBT_Transform_Input(ptr, in, direction);
    ERI_SBT_Transform_Output(ptr, out, l);
  }
}


void ERI_SBT_Transform_Input(
  ERI_SBT_t    *ptr,
  const double *in,
  int           direction
)
{
  ptr->transform_in(ptr->instance, in, direction);
}


void ERI_SBT_Transform_Output(
  ERI_SBT_t *ptr,
  double    *out,
  int        l
)
{
  ptr->transform_out(ptr->instance, out, l);
}



#if 0
void ERI_SBT_Transform_Real(
  ERI_SBT_t    *ptr,
  double       *out,
  const double *in,
  int           l,
  int           direction
)
{
  if (ptr->transform) {
    ptr->transform_real(ptr->instance, out, in, l, direction);
  } else {
    ERI_SBT_Transform_Input_Real(ptr, in, direction);
    ERI_SBT_Transform_Output_Real(ptr, out, l);
  }
}

void ERI_SBT_Transform_Input_Real(
  ERI_SBT_t    *ptr,
  const double *in,
  int           direction
)
{
  ptr->transform_in_real(ptr->instance, in, direction);
}


void ERI_SBT_Transform_Output_Real(
  ERI_SBT_t *ptr,
  double    *out,
  int        l
)
{
  ptr->transform_out_real(ptr->instance, out, l);
}
#endif

/* EOF */
