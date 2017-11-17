/*----------------------------------------------------------------------
  eri_logfsbt.h
 
  Oct. 2008, M. Toyoda
----------------------------------------------------------------------*/
#ifndef LIBERI_ERI_LOGFSBT_H_INCLUDED
#define LIBERI_ERI_LOGFSBT_H_INCLUDED

#include "sbt/eri_sbt.h"


ERI_SBT_t* ERI_LogFSBT_Init(
  int    lmax,  /* (IN) maximum number of angular momentum */
  int    ngrid, /* (IN) number of radial logarithm mesh */
  double rho0,  /* (IN) lower bound of rho (radial) mesh */
  double dt,    /* (IN) interval of t-mesh */
  double qmin,
  double qmax,
  int    nq
);


size_t ERI_LogFSBT_Required_Size(
  int lmax,  /* (IN) maximum number of angular momentum */
  int ngrid, /* (IN) number of radial logarithm mesh */
  int nq
);


#endif

/* EOF */
