/**********************************************************************
  AtomicDenF.c:

     AtomicDenF.c is a subroutine to calculate the atomic charge
     density of one atom specified "Gensi" at R.

  Log of AtomicDenF.c:

     22/Nov/2001  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <math.h>
#include "openmx_common.h"

double AtomicDenF(int spe, double R)
{
  int mp_min,mp_max,m;
  double h1,h2,h3,f1,f2,f3,f4;
  double g1,g2,x1,x2,y1,y2,f,df;
  double rm,y12,y22,a,b;
  double result;

  if (Spe_WhatAtom[spe]==0) return 0.0;

  mp_min = 0;
  mp_max = Spe_Num_Mesh_PAO[spe] - 1;

  if (Spe_PAO_RV[spe][Spe_Num_Mesh_PAO[spe]-1]<R){
    result = 0.0;
  }

  else if (R<Spe_PAO_RV[spe][0]){
    m = 4;
    rm = Spe_PAO_RV[spe][m];

    h1 = Spe_PAO_RV[spe][m-1] - Spe_PAO_RV[spe][m-2];
    h2 = Spe_PAO_RV[spe][m]   - Spe_PAO_RV[spe][m-1];
    h3 = Spe_PAO_RV[spe][m+1] - Spe_PAO_RV[spe][m];

    f1 = Spe_Atomic_Den[spe][m-2];
    f2 = Spe_Atomic_Den[spe][m-1];
    f3 = Spe_Atomic_Den[spe][m];
    f4 = Spe_Atomic_Den[spe][m+1];

    g1 = ((f3-f2)*h1/h2 + (f2-f1)*h2/h1)/(h1+h2);
    g2 = ((f4-f3)*h2/h3 + (f3-f2)*h3/h2)/(h2+h3);

    x1 = rm - Spe_PAO_RV[spe][m-1];
    x2 = rm - Spe_PAO_RV[spe][m];
    y1 = x1/h2;
    y2 = x2/h2;
    y12 = y1*y1;
    y22 = y2*y2;

    f =  y22*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
       + y12*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1);

    df = 2.0*y2/h2*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
       + y22*(2.0*f2 + h2*g1)/h2
       + 2.0*y1/h2*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1)
       - y12*(2.0*f3 - h2*g2)/h2;

    a = 0.5*df/rm;
    b = f - a*rm*rm;      
    result = a*R*R + b;
  }

  else{
    do{
      m = (mp_min + mp_max)/2;
      if (Spe_PAO_RV[spe][m]<R)
        mp_min = m;
      else 
        mp_max = m;
    }
    while((mp_max-mp_min)!=1);
    m = mp_max;

    if (m<2)
      m = 2;
    else if (Spe_Num_Mesh_PAO[spe]<=m)
      m = Spe_Num_Mesh_PAO[spe] - 2;

    /****************************************************
                   Spline like interpolation
    ****************************************************/
   
    if (m==1){
      h2 = Spe_PAO_RV[spe][m]   - Spe_PAO_RV[spe][m-1];
      h3 = Spe_PAO_RV[spe][m+1] - Spe_PAO_RV[spe][m];

      f2 = Spe_Atomic_Den[spe][m-1];
      f3 = Spe_Atomic_Den[spe][m];
      f4 = Spe_Atomic_Den[spe][m+1];

      h1 = -(h2+h3);
      f1 = f4;
    }
    else if (m==(Spe_Num_Mesh_PAO[spe]-1)){
      h1 = Spe_PAO_RV[spe][m-1] - Spe_PAO_RV[spe][m-2];
      h2 = Spe_PAO_RV[spe][m]   - Spe_PAO_RV[spe][m-1];

      f1 = Spe_Atomic_Den[spe][m-2];
      f2 = Spe_Atomic_Den[spe][m-1];
      f3 = Spe_Atomic_Den[spe][m];

      h3 = -(h1+h2);
      f4 = f1;
    }
    else{
      h1 = Spe_PAO_RV[spe][m-1] - Spe_PAO_RV[spe][m-2];
      h2 = Spe_PAO_RV[spe][m]   - Spe_PAO_RV[spe][m-1];
      h3 = Spe_PAO_RV[spe][m+1] - Spe_PAO_RV[spe][m];

      f1 = Spe_Atomic_Den[spe][m-2];
      f2 = Spe_Atomic_Den[spe][m-1];
      f3 = Spe_Atomic_Den[spe][m];
      f4 = Spe_Atomic_Den[spe][m+1];
    } 

    /****************************************************
                Calculate the value at R
    ****************************************************/

    g1 = ((f3-f2)*h1/h2 + (f2-f1)*h2/h1)/(h1+h2);
    g2 = ((f4-f3)*h2/h3 + (f3-f2)*h3/h2)/(h2+h3);

    x1 = R - Spe_PAO_RV[spe][m-1];
    x2 = R - Spe_PAO_RV[spe][m];
    y1 = x1/h2;
    y2 = x2/h2;

    result = y2*y2*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
           + y1*y1*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1);
  }

  return fabs(result);
}



