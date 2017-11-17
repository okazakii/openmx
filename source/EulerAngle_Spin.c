/**********************************************************************
  EulerAngle_Spin.c:

     EulerAngle_Spin.c is a subroutine to find the Euler angle of spin
     orientation from the density matrix.

  Log of EulerAngle_Spin.c:

     15/Feb/2006  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <math.h>
#include "openmx_common.h"

#ifdef c_complex
#include <complex.h>
#endif

#define MIND  1.0e-8 


void EulerAngle_Spin( int quickcalc_flag,
                      double Re11, double Re22,
                      double Re12, double Im12,
                      double Re21, double Im21,
                      double Nup[2], double Ndown[2],
                      double t[2], double p[2] )
{
  double phi,theta,d1,d2,d3,d4;
  double cop,sip,sit,cot,tmp,tmp1,tmp2,prod;
  double mx,my,mz,tn,absm;
  double S_coordinate[3];
  dcomplex bunbo,bunsi;
  dcomplex cd1,cd2,cd3,cd4,cNup,cNdown;

#ifdef c_complex
  double complex ctmp1,ctheta,csit,ccot;
#endif


#ifdef c_complex

  if (fabs(Re12)<1.0e-14){
    phi = PI*90.0/180.0;
  }
  else{

    bunbo.r = Re12 + Re21;       
    bunbo.i = Im12 + Im21;       
    bunsi.r = Re12 - Re21;       
    bunsi.i = Im12 - Im21;       

    tmp = -(bunsi.i*bunbo.r - bunsi.r*bunbo.i)/(bunbo.r*bunbo.r+bunbo.i*bunbo.i);
    phi = atan(tmp);
  }

  cop = cos(phi);
  sip = sin(phi);

  if (fabs(Re11 - Re22)<1.0e-14){
    ctheta = PI*90.0/180.0 + 0.0*I;
  }
  else {
    tmp1 = (Re12*cop - Im12*sip + Re21*cop + Im21*sip)/(Re11 - Re22); 
    tmp2 = (Re12*sip + Im12*cop - Re21*sip + Im21*cop)/(Re11 - Re22); 

    ctmp1 = tmp1 + tmp2*I;
    ctheta = catan(ctmp1);
  }

  csit = csin(ctheta); 
  ccot = ccos(ctheta); 

  cd1.r = 0.5*(Re11 + Re22);
  cd1.i = 0.0; 
  
  cd2.r = 0.5*creal(ccot)*(Re11 - Re22);
  cd2.i = 0.5*cimag(ccot)*(Re11 - Re22);

  cd3.r = 0.5*( (Re12*creal(csit)-Im12*cimag(csit))*cop
		-(Re12*cimag(csit)+Im12*creal(csit))*sip );

  cd3.i = 0.5*( (Re12*creal(csit)-Im12*cimag(csit))*sip
		+(Re12*cimag(csit)+Im12*creal(csit))*cop );

  cd4.r = 0.5*( (Re21*creal(csit)-Im21*cimag(csit))*cop
		+(Re21*cimag(csit)+Im21*creal(csit))*sip );

  cd4.i = 0.5*(-(Re21*creal(csit)-Im21*cimag(csit))*sip
	       +(Re21*cimag(csit)+Im21*creal(csit))*cop );

  cNup.r   = cd1.r + cd2.r + cd3.r + cd4.r;
  cNup.i   = cd1.i + cd2.i + cd3.i + cd4.i;
  cNdown.r = cd1.r - cd2.r - cd3.r - cd4.r;
  cNdown.i = cd1.i - cd2.i - cd3.i - cd4.i;

  Nup[0] = cNup.r;
  Nup[1] = cNup.i;

  Ndown[0] = cNdown.r;
  Ndown[1] = cNdown.i;

  t[0] = creal(ctheta);
  t[1] = cimag(ctheta);

  p[0] = phi;
  p[1] = 0.0;

#else


  /*
  if (fabs(Re12)<MIND){

    mx =  2.0*Re12;
    my = -2.0*Im12;
    phi = acos( mx/sqrt(mx*mx + my*my + MIND*MIND) );
  }
  else{
    phi = atan(-Im12/Re12);
  }

  cop = cos(phi);
  sip = sin(phi);

  if (fabs(Re11 - Re22)<MIND){

    mx =  2.0*Re12;
    my = -2.0*Im12;
    mz = Re11 - Re22;

    theta = acos( mz/sqrt(mx*mx + my*my + mz*mz + MIND*MIND) );
  }
  else {
    theta = atan( 2.0*(Re12*cop - Im12*sip) / (Re11 - Re22) );
  }

  sit = sin(theta); 
  cot = cos(theta); 

  d1 = 0.5*(Re11 + Re22);
  d2 = 0.5*cot*(Re11 - Re22);
  d3 = (Re12*cop - Im12*sip)*sit;

  Nup[0]   = d1 + d2 + d3;
  Nup[1]   = 0.0;

  Ndown[0] = d1 - d2 - d3;
  Ndown[1] = 0.0;

  t[0] = theta;
  t[1] = 0.0;

  p[0] = phi;
  p[1] = 0.0;
  */



  /*
  printf("theta=%10.5f phi=%10.5f Re11=%10.5f Re22=%10.5f Re12=%10.5f Im12=%10.5f d1=%10.5f d2=%10.5f d3=%10.5f\n",
          180.0*theta/PI,180.0*phi/PI,Re11,Re22,Re12,Im12,d1,d2,d3);
  */



  mx =  2.0*Re12;
  my = -2.0*Im12;
  mz = Re11 - Re22;
  
  xyz2spherical(mx,my,mz,0.0,0.0,0.0,S_coordinate);
  
  absm = S_coordinate[0];
  theta = S_coordinate[1];
  phi = S_coordinate[2]; 
  tn = Re11 + Re22;
  
  Nup[0]   = 0.5*(tn + absm);
  Nup[1]   = 0.0;
    
  Ndown[0] = 0.5*(tn - absm);
  Ndown[1] = 0.0;
  
  t[0] = theta;
  t[1] = 0.0;
  
  p[0] = phi;
  p[1] = 0.0;

#endif

}
