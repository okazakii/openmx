/*----------------------------------------------------------------------
  eri_linear.h

  Coded by M. Toyoda, May 2009, JAIST/RCIS
----------------------------------------------------------------------*/
#ifndef LIBERI_ERI_LINFSBT_H_INCLUDED
#define LIBERI_ERI_LINFSBT_H_INCLUDED

#include "sbt/eri_sbt.h"
#include <stdlib.h>


typedef struct ERI_LinFSBT_struct ERI_LinFSBT_t;


ERI_SBT_t* ERI_LinFSBT_Init(
  int    nmesh,     /* (IN) number of mesh points */
  int    max_order, /* (IN) maximum order */
  double range_r    /* (IN) range of k-mesh */
);


size_t ERI_LinFSBT_Required_Size(
  int    nmesh,
  int    max_order
);

#endif /* LIBERI_ERI_LINFSBT_H_INCLUDED */
