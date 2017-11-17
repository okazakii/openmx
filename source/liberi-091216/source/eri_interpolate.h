/*----------------------------------------------------------------------
  eri_interpolate.h

  Coded by M. Toyoda, May 2009, JAIST/RCIS.

  * This uses CLAPACK, so that you must link with appropreate libraries
    by specifing, for example, the following link options:

    cc *.o -llapack -lblas -lf2c (or, -lg2c for gcc).
----------------------------------------------------------------------*/
#ifndef LIBERI_ERI_INTERPOLATE_H_INCLUDED
#define LIBERI_ERI_INTERPOLATE_H_INCLUDED


double ERI_Interpolate_Linear_Real(
  double        x, 
  const double *in,
  const double *xmesh, 
  int           nmesh
);


void ERI_Interpolate_Linear_Complex(
  double        out[2],
  double        x, 
  const double *in,
  const double *xmesh, 
  int           nmesh
);



typedef struct ERI_CSpline_Struct ERI_CSpline_t;


ERI_CSpline_t* ERI_CSpline_Init(
  const double *sx,
  const double *sy,
  int           ns
);


ERI_CSpline_t* ERI_CSpline_Init_Complex(
  const double *sx,
  const double *sy,
  int           ns
);


ERI_CSpline_t* ERI_CSpline_Init_Plus2(
  const double *sx,
  const double *sy,
  int           ns
);


ERI_CSpline_t* ERI_CSpline_Init_Complex_Plus2(
  const double *sx,
  const double *sy,
  int           ns
);


void ERI_CSpline_Free(ERI_CSpline_t *ptr);


double ERI_CSpline_Eval(
  double               x,
  const ERI_CSpline_t *ptr
);


void ERI_CSpline_Eval_Complex(
  double               y[2],
  double               x,
  const ERI_CSpline_t *ptr
);


#endif
