/*----------------------------------------------------------------------
  demo.c

  Sample program of LIBERI

  Coded by TOYODA MAsayuki, June 2009.
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "eri.h"

#define NMESH 8192


int main(void)
{
  int i;
  double r, I4[2];
  const double rmax = 300.0;

  ERI_t *eri;
  ERI_Orbital_t p[4];

  double xr[NMESH], gto_s[NMESH];

  /* radial mesh and GTO */
  for (i=0; i<NMESH; i++) { 
    r = rmax*(double)i/(double)NMESH; 
    xr[i] = r;
    gto_s[i] = exp(-r*r); 
  } 

  /*  prepare orbital informations */
  p[0].fr = gto_s; p[0].xr = xr; p[0].ngrid = NMESH; p[0].l = 0; p[0].m = 0;
  p[1].fr = gto_s; p[1].xr = xr; p[1].ngrid = NMESH; p[1].l = 0; p[1].m = 0;
  p[2].fr = gto_s; p[2].xr = xr; p[2].ngrid = NMESH; p[2].l = 0; p[2].m = 0;
  p[3].fr = gto_s; p[3].xr = xr; p[3].ngrid = NMESH; p[3].l = 0; p[3].m = 0;

  /* positions */
  p[0].c[0] =  0.5; p[0].c[1] =  0.5; p[0].c[2] = 0.0;
  p[1].c[0] =  0.5; p[1].c[1] = -0.5; p[1].c[2] = 0.0;
  p[2].c[0] = -0.5; p[2].c[1] =  0.5; p[2].c[2] = 0.0;
  p[3].c[0] = -0.5; p[3].c[1] = -0.5; p[3].c[2] = 0.0;

  /* initialize LIBERI */
  eri = ERI_Init(15, 8, 1024, 100, ERI_SH_COMPLEX, NULL);
  if (NULL == eri) {
    fprintf(stderr, "*** Failed to initialize LIBERI.\n");
    return -1;
  }

  /* calculate ERI */
  ERI_Integral(eri, I4, NULL, &p[0], &p[1], &p[2], &p[3], ERI_NOSCREEN);
  printf("  %10.6f  %10.6f\n", I4[0], I4[1]);
 
  /* release LIBERI */
  ERI_Free(eri);

  return 0;
}


/* EOF */
