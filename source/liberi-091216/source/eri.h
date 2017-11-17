/*----------------------------------------------------------------------
  erl.h

  16 Feb. 2009, M. Toyoda
----------------------------------------------------------------------*/
#ifndef LIBERI_ERI_H_INCLUDED
#define LIBERI_ERI_H_INCLUDED

#include <stdlib.h>
#include <float.h>


#define ERI_SBT_LOG    1 /* Siegman-Talman method */
#define ERI_SBT_LINEAR 2 /* Toyoda-Ozaki method */
#define ERI_NOSCREEN -1.0

#define ERI_SH_REAL    1  /* Real-valued SH is used */
#define ERI_SH_COMPLEX 2  /* Complex-valued SH is used */


#if defined(__cplusplus)
extern "C" {
#endif


/*----------------------------------------------------------------------
  Initializer and Finalizer
----------------------------------------------------------------------*/
typedef struct ERI_Struct ERI_t;

const char* ERI_Version(void);


typedef struct {
  int    sbttype;
  double rmax;
  double rho0;
  double qmin;
  double qmax;
  int    nq;
} ERI_Init_Misc; 


/*----------------------------------------------------------------------
  ERI_Init

  Initializes LIBERI.

  IN:
    lmax    : maximum angular momentum 
    lmax_gl : maximum angular momenrum for the final summation
    ngrid   : number of radial mesh points
    ngl     : number of abscissas for GL quadrature
    sbttyp  : type of SBT routine (ERI_SBT_LOG or ERI_SBT_LINEAR)
    rmax    : range of r-mesh
    info    : additional information 
              (this is optional. if this is NULL, the default settings
               are used).
    rcsh    : real- or complex-valued SH is used.
  RETURN:
    pointer to an ERI_Struct object.
----------------------------------------------------------------------*/
ERI_t* ERI_Init(
  int lmax,
  int lmax_gl,
  int ngrid,
  int ngl,
  int rcsh,
  const ERI_Init_Misc* info
);


/*----------------------------------------------------------------------
  ERI_Required_Size
 
  Estimates the byte size which is used internally by LIBERI.

  See ERI_Init for the explantation of the arguments.
----------------------------------------------------------------------*/
size_t ERI_Required_Size(
  int lmax,
  int lmax_gl,
  int ngrid,
  int ngl,
  int rcsh,
  const ERI_Init_Misc* info
);


/*----------------------------------------------------------------------
  ERI_Free

  Finalizes LIBERI.

  IN:
    ptr : pointer returned by ERI_Init.
----------------------------------------------------------------------*/
void ERI_Free(ERI_t *ptr);


/*----------------------------------------------------------------------
  ERI_lmax, ERI_lmax_gl, ERI_ngird, ERI_ngl, ERI_Misc

  Returns the initialization parameters.
----------------------------------------------------------------------*/
int                  ERI_lmax    (const ERI_t *ptr);
int                  ERI_lmax_gl (const ERI_t *ptr);
int                  ERI_ngl     (const ERI_t *ptr);
const ERI_Init_Misc* ERI_Misc    (const ERI_t *ptr);
int                  ERI_ngrid   (const ERI_t *ptr);
int                  ERI_rcsh    (const ERI_t *ptr);

/*----------------------------------------------------------------------
  ERI_Mesh_r, ERI_Mesh_k, ERI_Mesh_dr, ERI_Mesh_dk

  Returns the radial mesh points or the intervals of the mesh.
----------------------------------------------------------------------*/
double        ERI_Mesh_r       (const ERI_t *ptr, int i);
double        ERI_Mesh_k       (const ERI_t *ptr, int i);
double        ERI_Mesh_dr      (const ERI_t *ptr, int i);
double        ERI_Mesh_dk      (const ERI_t *ptr, int i);


/*----------------------------------------------------------------------
  ERI_Mesh_Array_r, ERI_Mesh_Array_k

  Returns the all radial mesh points.
----------------------------------------------------------------------*/
const double* ERI_Mesh_Array_r (const ERI_t *ptr);
const double* ERI_Mesh_Array_k (const ERI_t *ptr);


/*----------------------------------------------------------------------
  ERI_Mesh_Array_glx, ERI_Mesh_Array_glw

  Returns the abscissas and the weight factors for the GL quadrature.
----------------------------------------------------------------------*/
const double* ERI_Mesh_Array_glx(const ERI_t *ptr);
const double* ERI_Mesh_Array_glw(const ERI_t *ptr);




/*----------------------------------------------------------------------
  High-level `black box' functions
----------------------------------------------------------------------*/

typedef struct {
  double       *fr;   /* radial function */
  double       *xr;   /* radial grid */
  unsigned int  ngrid; /* number of radial grid */
  int           l;     /* angular momentum */
  int           m;     /* angular momentum */
  double        c[3];  /* position */
} ERI_Orbital_t;



/*----------------------------------------------------------------------
  ERI_Overlap

  Calculated the transformed overlap function for given orbital 
  functions.
 
  IN:
    ptr  : pointer returned by ERI_Init
    orb1 : an orbital function
    orb2 : another orbital fuction
    cx   : parameter for expansion center 

  OUT:
    F    : claculated overlap function
    dF   : derivatives of F
           (this is optional. if you do not need the derivatives, 
            specify NULL here.)
----------------------------------------------------------------------*/
void ERI_Overlap(
  ERI_t  *ptr,
  double *F,
  double *dF[3],
  const ERI_Orbital_t *orb1, 
  const ERI_Orbital_t *orb2,
  double cx
);



/*----------------------------------------------------------------------
  ERI_Integral

  Calculates ERI for given four orbital functions.

  IN:
    ptr  : pointer returned by ERI_Init.
    orb1 : orbital function #1
    orb2 : orbital function #2
    orb3 : orbital function #3
    orb4 : orbital function #4
   
  OUT:
    I4   : calculated ERI (complex number)
    dI4  : derivatives of ERI
           (this is optional. if you do not need the derivatives,
            specify NULL here.)
----------------------------------------------------------------------*/
void ERI_Integral(
  ERI_t  *ptr,
  double I4[2],
  double *dI4,
  const ERI_Orbital_t *orb1,
  const ERI_Orbital_t *orb2,
  const ERI_Orbital_t *orb3,
  const ERI_Orbital_t *orb4,
  double               scr
);


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
);



/*----------------------------------------------------------------------
  ERI_Center_r2

  IN:
    orb1 : orbital #1 
    orb2 : orbital #2 

  RETURN:
    ratio for the expansion center.
----------------------------------------------------------------------*/
double ERI_Center_r2(
  const ERI_Orbital_t *orb1,
  const ERI_Orbital_t *orb2
);


/*----------------------------------------------------------------------
  ERI_Center_DK


  IN:
    ptr  : pointer returned by ERI_Init
    orb1 : orbital #1
    orb1 : orbital #2

  RETURN:
    ratio for the expansion center.
----------------------------------------------------------------------*/
double ERI_Center_DK(
  ERI_t        *ptr,
  const ERI_Orbital_t *orb1,
  const ERI_Orbital_t *orb2
);




/*----------------------------------------------------------------------
  Low-level functions
----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
  ERI_Size_of_*

  Returns required byte size for the matrices.
----------------------------------------------------------------------*/
size_t ERI_Size_of_Orbital (const ERI_t *ptr);
size_t ERI_Size_of_Gamma   (const ERI_t *ptr);
size_t ERI_Size_of_Alpha   (const ERI_t *ptr);
size_t ERI_Size_of_Overlap (const ERI_t *ptr);
size_t ERI_Size_of_GLF     (const ERI_t *ptr);


/*----------------------------------------------------------------------
  ERI_Transform_Orbital

  Given input function, this performs the SBT.

  Note that the radial mesh for input function (fr) is NOT the r-mesh 
  defined by LIBERI. It is the user-defined mesh (xr). The user-defined
  radial mesh (xr) can be of any type, as long as it is in ascenfing 
  order.

  The k-mesh for the output function (fk) is defined by LIBERI.

  IN:
    ptr : pointer returned by ERI_Init
    fr  : values of input function at each 
    xr  : points of the radial mesh for fr.
    n   : number of the radial mesh for fr.
    l   : angular momenrum

  OUT:
    fk  : transformed orbital
          (The size of the array should be no smaller then the size 
           returned by ERI_Size_of_Orbital.)
----------------------------------------------------------------------*/
void ERI_Transform_Orbital(
  ERI_t        *ptr,
  double       *fk,  /* (OUT) transformed orbital */
  const double *fr,  /* (IN) radial function */
  const double *xr,  /* (IN) radial mesh */
  int           n,   /* (IN) number of mesh poitns */
  int           l    /* (IN) angular momentum */
);


#if 0
/*----------------------------------------------------------------------
  ERI_LL_Gamma

  Calculates the gamma-functions.

  IN:
    ptr : pointer returned by ERI_Init
    fk  : transformed function (returned by ERI_Transform_Orbital)
    rd  : displacement from the center
  
  OUT: 
    g   : gamma-functions.
          (The size of the array should be no smaller than the size 
           returned by ERI_Size_of_Gamma.)
----------------------------------------------------------------------*/
void ERI_LL_Gamma(
  ERI_t        *ptr,
  double       *g,   /* (OUT) Gamma term */
  const double *fk,  /* (IN)  Transformed radial function */
  double        rd   /* (IN)  Displacemenent from the center */
);
#endif


/*----------------------------------------------------------------------
  ERI_LL_Alpha

  Calculates the alpha-functions.

  IN:
    ptr        : pointer returned by ERI_Init
    g          : gamma-functions (returned by ERI_LL_Gamma)
    theta, phi : angle displacement from the center
    m          : magnetic angular momentum of the orbital

  OIT:
    alp        : alpha-functions
                 (The size of the array shuold be no smaller than the 
                  size returned by ERI_Size_of_Alpha.)
----------------------------------------------------------------------*/
void ERI_LL_Alpha(
  ERI_t        *ptr,
  double       *alp,   /* (OUT) Alpha terms */
  const double *g,     /* (IN)  Gamma terms */
  double        theta,
  double        phi, 
  int           l,     /* (IN)  Angular momentum of the orbital */
  int           m
);


/*----------------------------------------------------------------------
  ERI_LL_Overlap
  
  Calculatets the overlap function.

  IN:
    ptr : pointer returned by ERI_Init
    a1  : alpha-function of orbital 1
    a2  : alpha-function of orbital 2
  
  OUT:
    p   : overlap function
          (The size of the array should be no smaller than the size 
           returned by ERI_Size_of_Overlap.) 
----------------------------------------------------------------------*/
void ERI_LL_Overlap(
  ERI_t        *ptr,
  double       *p,
  const double *a1,
  const double *a2 
);


/*----------------------------------------------------------------------
  ERI_Transform_Overlap
 
  Performes SBT of the overlap functions.
 
  If you do not need the derivatives, specify NULL for dp and df.
 
  IN:
    ptr : pointer returned by ERI_Init
    p   : overlap function
    dp  : derivatives of overlap function (optional)

  OUT:
    f   : transformed overlap functions 
    df  : derivatives wrt nucleus position in x, y and z direction,
          respectively (optional).
          (The size of each array should be no smaller than the size 
           returned by ERI_Size_of_Overlap.)
----------------------------------------------------------------------*/
void ERI_Transform_Overlap(
  ERI_t        *ptr,
  double       *f,
  const double *p
);
#if 0
void ERI_Transform_Overlap(
  ERI_t        *ptr,
  double       *f,
  double       *df[3],
  const double *p,
  const double *dp[3]
);
#endif

/*----------------------------------------------------------------------
  ERI_GL_Interpolate

  Interpolates the transformed overlap function from SBT-mesh to 
  GL-mesh.

  IN:
    ptr : poitner returned by ERI_Init
    f   : overlap function on SBT-mesh

  OUT:
    glf : overlap function on GL-mesh
----------------------------------------------------------------------*/
void ERI_GL_Interpolate(
  ERI_t        *ptr,
  double       *glf,
  const double *f
);


/*----------------------------------------------------------------------
  ERI_LL_Gamma

  Calculates the gamma functions and the derivatives.

  IN:
    ptr : poitner returned by ERI_Init
    fk  : transformed function returned by ERI_Transform_Orbital.
    rd  : displacement from the center

  OUT:
    g   : gamma-functions
    dg  : derivatives of the gamma-functions
          (The size of each array should be no smaller than the size 
           returned by ERI_Size_of_Gamma)
----------------------------------------------------------------------*/
void ERI_LL_Gamma(
  ERI_t        *ptr,
  double       *g,   /* (OUT) Gamma terms */
  double       *dg,  /* (OUT) Derivatives of the Gamma terms */
  const double *fk,  /* (IN) Transformed radidal function */
  double        rd   /* (IN) Displacement from the center */
);


/*----------------------------------------------------------------------
  ERI_LL_Alpha_d

  Calculates the alpha-funtions and the derivatives.

  IN:
    ptr            : pointer returned by ERI_Init
    g              : gamma-functions
    dg             : derivatives of gamma-functions
    rd, theta, phi : displacement from the center in spherical coord.
    l, m           : angular momentum or the orbital. 

  OUT:
    alp  : alpha-functions
    dalp : derivatives of alpha-functions wrt nuclear position
           in x, y and z directions, respectively.
           (The size of each array should be no smaller than the size 
            returned by ERI_Size_of_Alpha.)
----------------------------------------------------------------------*/
void ERI_LL_Alpha_d(
  ERI_t        *ptr,
  double       *alp,
  double       *dalp[3],
  const double *g,
  const double *dg,
  double        rd,
  double        theta,
  double        phi,
  int           l,
  int           m
);


/*----------------------------------------------------------------------
  ERI_LL_Overlap_d
 
  Calculates the overlap functions and the derivatives.

  IN:
    ptr : pointer returned by ERI_Init
    a1  : alpha-function of orbtial 1
    da1 : derivatives of alpha-function of orbital 1 
    a2  : alpha-function of orbtial 2
    da2 : derivatives of alpha-function of orbital 2 
    cx  : ratio for expansion center 

  OUT:
    p   : overlap function
    dp  : derivatives of overlap function
          (The size of each array should be no smaller than the size 
           returned by ERI_Size_of_Overlap.)
----------------------------------------------------------------------*/
void ERI_LL_Overlap_d(
  ERI_t        *ptr,
  double       *p,
  double       *dp[3],
  const double *a1,
  const double *da1[3],
  const double *a2,
  const double *da2[3],
  double        x
);




/*----------------------------------------------------------------------
  ERI_Integral_GL

  Calculates ERI for given overlap functions.

  IN:
    ptr : pointer returned by ERI_Init
    p12 : transformed overlap function for orbiatal 1 and 2
    p34 : transformed overlap function for orbiatal 3 and 4
    rd, theta, phi : displacement of expansion centers in sph. coord. 
    omega          : screening factor 
                     (speficy ERI_NOSCREEN if you do not need it)
    lmax           : cut-off of the angular momentum summation

  OUT:
    I4  : ERI (complex number)
----------------------------------------------------------------------*/
void ERI_Integral_GL(
  ERI_t        *ptr,
  double        I4[2],
  const double *p12,
  const double *p34,
  double        rd,
  double        theta,
  double        phi,
  double        omega,
  int           lmax1
);


void ERI_Integral_GL_d(
  ERI_t        *ptr,
  double        I4[2],  /* (OUT) */
  double        dI4[4][3][2],
  const double *F1,     /* (IN) Overlap matrix */
  const double *F2,     /* (IN) */
  const double *dF1[3],    /* (IN) Overlap matrix */
  const double *dF2[3],     /* (IN) */
  double        R,      /* (IN) Displacement of two expansion centers */
  double        theta,      
  double        phi,
  double        cx12,
  double        cx34,
  double        delta,
  double        omega,  /* (IN) screening parameter */
  int           lmax1
);



#if 0
void ERI_Integral_GL_PrejY(
  ERI_t        *solver,
  const double *F1,    /* (IN) Overlap matrix */
  const double *F2,    /* (IN) */
  double       *R,     /* (IN) Displacement of two expansion centers */
  double       *theta,      
  double       *phi,
  int           numR,
  double        omega,  /* (IN) screening parameter */
  int          *mul_j2, /* [jmax1*jmax1*jmax1] */
  double       *mul_gc, /* [jmax1*jmax1*jmax1] */
  int          *mul_n,  /* [jmax1*jmax1] */
  double       *prej,   /* [lmax*ngl*numR] */
  double       *preY,   /* [numR*jmax1] */
  int          *num_minimalR,
  int          *normalR,      /* [numR*numR] */
  int          *num_normalR    /* [numR] */
);




void ERI_Integral_GL_Post(
  ERI_t        *solver,
  double       *I4,    /* (OUT) [numR] */
  const double *F1,    /* (IN) Overlap matrix */
  const double *F2,    /* (IN) */
  int           numR,
  int          *mul_j2, /* [jmax1*jmax1*jmax1] */
  double       *mul_gc, /* [jmax1*jmax1*jmax1] */
  int          *mul_n,  /* [jmax1*jmax1] */
  double       *prej,   /* [lmax*ngl*numR] */
  double       *preY,    /* [numR*jmax1] */
  int           num_minimalR,
  int          *normalR,      /* [numR*numR] */
  int          *num_normalR    /* [numR] */
);
#endif


void ERI_Integral_GL_PrejY(
  ERI_t        *solver,
  double       *R,     /* (IN) Displacement of two expansion centers */
  double       *theta,      
  double       *phi,
  int           numR,
  double        omega,  /* (IN) screening parameter */
  double       *prej,   /* [lmax*ngl*numR] */
  double       *preY,   /* [numR*jmax1] */
  int          *mul_j2, /* [jmax1*jmax1*jmax1] */
  double       *mul_gc, /* [jmax1*jmax1*jmax1] */
  int          *mul_n,  /* [jmax1*jmax1] */
  int          *minimalR,     /* [numR] */
  int          *num_minimalR
);


void ERI_Integral_GL_Post(
  ERI_t        *solver,
  double       *I4,    /* (OUT) [numR] */
  const double *F1,    /* (IN) Overlap matrix */
  const double *F2,    /* (IN) */
  /*const double *Q1,*/    /* (IN) Overlap matrix */
  /*const double *Q2,*/    /* (IN) */
  int           numR,
  double       *prej,   /* [lmax*ngl*numR] */
  double       *preY,    /* [numR*jmax1] */
  int          *mul_j2, /* [jmax1*jmax1*jmax1] */
  double       *mul_gc, /* [jmax1*jmax1*jmax1] */
  int          *mul_n,  /* [jmax1*jmax1] */
  int          *minimalR,     /* [numR] */
  int           num_minimalR
);

void ERI_Integral_GL_Post2(
  ERI_t        *solver,
  double       *I4,    /* (OUT) [numR] */
  const double *F1,    /* (IN) Overlap matrix */
  const double *F2,    /* (IN) */
  /*const double *Q1,*/    /* (IN) Overlap matrix */
  /*const double *Q2,*/    /* (IN) */
  int           numR,
  double       *prej,   /* [lmax*ngl*numR] */
  double       *preY,    /* [numR*jmax1] */
  int          *mul_j2, /* [jmax1*jmax1*jmax1] */
  double       *mul_gc, /* [jmax1*jmax1*jmax1] */
  int          *mul_n,  /* [jmax1*jmax1] */
  int          *minimalR,     /* [numR] */
  int           num_minimalR
);



void ERI_Integral_GL_X(
  ERI_t        *solver,
  double       *X4,    /* (OUT) [numR] */
  const double *F1,    /* (IN) Overlap matrix */
  const double *F2,    /* (IN) */
  int           numR,
  double       *prej,   /* [lmax*ngl*numR] */
  int          *mul_j2, /* [jmax1*jmax1*jmax1] */
  double       *mul_gc, /* [jmax1*jmax1*jmax1] */
  int          *mul_n,  /* [jmax1*jmax1] */
  int          *minimalR,    /* [numR] */
  int           num_minimalR
);

void ERI_Integral_GL_X_Post(
  ERI_t        *solver,
  double       *I4,       /* (OUT) [numR] */
  const double *X,        /* (IN)  [numR] */
  int           numR,
  const double *preY,     /* (IN)  [numR*jmax1] */
  const int    *minimalR  /* (IN)  [numR] */
);


#if defined(__cplusplus)
}
#endif

#endif /* LIBERI_ERI_H_INCLUDED */
/* EOF */
