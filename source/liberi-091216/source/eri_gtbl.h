/*----------------------------------------------------------------------
  eri_gtbl.h
----------------------------------------------------------------------*/
#ifndef LIBERI_ERI_GTBL_H_INCLUDED
#define LIBERI_ERI_GTBL_H_INCLUDED

#include <stdlib.h>
#include "eri.h"

typedef struct ERI_Gtbl_Struct ERI_Gtbl_t;

#define ERI_GTBL_PACK_J1J2(j1,j2) ((long)(j1%32768)*65536+(long)(j2%32768))
#define ERI_GTBL_UNPACK_J1(j1j2) ((j1j2/65536)%32768)
#define ERI_GTBL_UNPACK_J2(j1j2) (j1j2%32768)


/*----------------------------------------------------------------------
  ERI_Gtbl_Init

  Initialize Gaunt coefficient table. 
  
  IN:
    lmax : maximum value of angular momentum 
    type : ERI_SH_COMPLEX or ERI_SH_REAL
  
  OUT:
    return : pointer to a new object, or NULL on error.
----------------------------------------------------------------------*/
ERI_Gtbl_t* ERI_Gtbl_Init(int lmax, int type);


/*----------------------------------------------------------------------
  ERI_Gtbl_Free

  Destroy Gaunt coefficient table 

  IN:
    p : pointer to object
----------------------------------------------------------------------*/
void ERI_Gtbl_Free(ERI_Gtbl_t *p);


/*----------------------------------------------------------------------
  ERI_Gtbl_Required_Size

  IN:
    lmax : maximum value of angular momentum 
    type : ERI_SH_COMPLEX or ERI_SH_REAL
  
  OUT:
    return : required byte size of an object
----------------------------------------------------------------------*/
size_t ERI_Gtbl_Required_Size(int lmax, int type);


int ERI_Gtbl_index(ERI_Gtbl_t* p, int j);

int ERI_Gtbl_jmax(ERI_Gtbl_t* p);
const int* ERI_Gtbl_n(ERI_Gtbl_t* p);
const long* ERI_Gtbl_j1j2(ERI_Gtbl_t *p);
const double* ERI_Gtbl_gc(ERI_Gtbl_t *p);


#endif /* LIBERI_ERI_GTBL_H_INCLUDED */
