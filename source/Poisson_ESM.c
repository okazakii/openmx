/**********************************************************************
  Poisson_ESM.c:

     Poisson_ESM.c is a subroutine to solve Poisson's equations
     for effective screening medium (ESM) method calculations.
     This subroutine is written based on "Poisson.c".

  References for ESM method:

     [1] M.Otani and O.Sugino, PRB 73, 115407 (2006).
     [2] O.Sugino et al., Surf.Sci., 601, 5237 (2007).
     [3] M.Otani et al., J.Phys.Soc.Jp., 77, 024802 (2008).
     [4] T. Ohwaki et al., submitted to J. Chem. Phys. (2011).

  Log of Poisson_ESM.c:

     01/10/2011  Released by T.Ohwaki and M.Otani

***********************************************************************/

#define  measure_time   0

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "openmx_common.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

#ifdef fftw2 
#include <fftw.h>
#else
#include <fftw3.h> 
#endif

static void One_dim_FFT(FILE *fp, 
			int sgn1, int sgn2, 
			int kn2,  int kn3, 
			double ***ReF2, double ***ImF2,
                        double prefac);


double Poisson_ESM(int fft_charge_flag,
                   double ***ReV1, double ***ImV1,
                   double ***ReV2, double ***ImV2)
{ 
  int i,j,k,k1,k2,k3,kk2,kz1,kk1;
  double time0;
  double GridArea,tmpv[4];
  double tmp0,sk1,sk2,sk3,skz1,tot_den;
  double x,y,z,Gx,Gy,Gz,Eden[2],DenA,G2;
  double ReTmp,ImTmp,c_coe,p_coe;
  double *tmp_array0;
  double *tmp_array1;
  double TStime,TEtime;

  double Gz2,Gp,Gp2,z0,z1,zz;
  double tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
  double tmp1r,tmp1i,tmp2r,tmp2i,tmp3r,tmp3i,tmp4r,tmp4i,tmp5r,tmp5i;
  int numprocs,myid,tag=999,ID,IDS,IDR;

  FILE *fp;
  char fname[YOUSO10];

  MPI_Status stat;
  MPI_Request request;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (myid==Host_ID) { 
    printf("<Poisson_ESM> Poisson's equation using FFT & ESM...\n");
    if (ESM_switch == 1) printf("<Poisson_ESM> Boundary condition = vacuum|vacuum|vacuum\n");
    if (ESM_switch == 2) printf("<Poisson_ESM> Boundary condition = metal|vacuum|metal \n");
    if (ESM_switch == 3) printf("<Poisson_ESM> Boundary condition = vacuum|vacuum|metal \n");
    if (ESM_switch == 4) printf("<Poisson_ESM> Boundary condition = uniform electric field \n");
  }

  MPI_Barrier(mpi_comm_level1);
  dtime(&TStime);

  Cross_Product(gtv[2],gtv[3],tmpv);
  GridArea = sqrt(tmpv[1]*tmpv[1] + tmpv[2]*tmpv[2] + tmpv[3]*tmpv[3]) * 0.529177249;

  /****************************************************
                 FFT of charge density 
  ****************************************************/

  if (fft_charge_flag==1) FFT_Density(0,ReV1,ImV1,ReV2,ImV2);

  if (myid==Host_ID){ 
    printf("<Poisson_ESM> Total number of electrons = %12.9f \n",ReV2[0][0][0]*GridVol);fflush(stdout);
  }

  if (myid==Host_ID){

    sprintf(fname,"%s%s.ESM.dcharge",filepath,filename);

    if ( (fp=fopen(fname,"w"))!=NULL ) {

      fprintf(fp,"## Check of difference charge density (e/Ang)      ##\n");fflush(stdout);
      fprintf(fp,"## Grid : x-corrdinate (Ang) : delta-rho(G_||=0,z) ##\n");fflush(stdout);
      One_dim_FFT(fp,1,1,0,0,ReV2,ImV2,GridArea);

      fprintf(fp,"\n");
      fclose(fp);
    }
    else {
      printf("Failure in saving *.ESM.dcharge.\n");
    }
  }

  /************************************************************

      One_dim_FFT(fp,n1,n2,m1,m2,ReV,ImV,prefac)

       n1= 1: for check of output (ReV & ImV are not changed)
         = 2: for calculation (ReV & ImV are changed)
       n2= 1: FFT
         =-1: inverse FFT
       m1,m2: Gx,Gy
       prefac: prefactor of output data

  *************************************************************/

  /***************************************************************
      ESM calculation for vacuum|slab|vacuum boundary condition
  ***************************************************************/

  tmp0 =1.0/(double)Ngrid1/(double)Ngrid2/(double)Ngrid3;

  /*******************************************************************************/
  /*  taking phase factor coming from grid-origin shift ( rho{Gz}*exp{i*Gz*dz} ) */
  /*******************************************************************************/

  tmp2 = tv[1][1]/2.0; 

  for (k2=0; k2<My_NGrid2_Poisson; k2++){
    for (k3=0; k3<Ngrid3; k3++){
      for (k1=0; k1<Ngrid1; k1++){

	if (k1<Ngrid1/2) sk1 = (double)k1;
	else             sk1 = (double)(k1 - Ngrid1);

	Gz  = sk1*rtv[1][1];

	tmp1r = ReV2[k2][k1][k3]*cos(Gz*tmp2) + ImV2[k2][k1][k3]*sin(Gz*tmp2);
	tmp1i = ImV2[k2][k1][k3]*cos(Gz*tmp2) - ReV2[k2][k1][k3]*sin(Gz*tmp2);

	ReV2[k2][k1][k3] = tmp1r;
	ImV2[k2][k1][k3] = tmp1i;

      }
    }
  }

  /*
    if (myid==Host_ID){
    printf(" $$$ check for delta-rho(G_||=0,z) (2) $$$ \n");
    One_dim_FFT(fp,1,1,0,0,ReV2,ImV2,GridArea); 
    }
  */

  /***************************************************************************/
  /*  "vacuum|vacuum|vacuum" boundary condition (for isolated slab systems)  */
  /***************************************************************************/

  if (ESM_switch==1){

    z0 = tv[1][1]/2.0;

    /* * * *  G_|| != 0  * * * */

    for (k2=0; k2<My_NGrid2_Poisson; k2++){

      kk2 = k2 + Start_Grid2[myid];

      if (kk2<Ngrid2/2) sk2 = (double)kk2;
      else              sk2 = (double)(kk2 - Ngrid2);

      for (k3=0; k3<Ngrid3; k3++){

	if (k3<Ngrid3/2) sk3 = (double)k3;
	else             sk3 = (double)(k3 - Ngrid3);

	Gx  = sk2*rtv[2][2] + sk3*rtv[3][2];  /* original Gy,Gz -> Gx,Gy */
	Gy  = sk2*rtv[2][3] + sk3*rtv[3][3];
	Gp2 = Gx*Gx + Gy*Gy;
	Gp  = sqrt(Gp2);

	if (kk2!=0 || k3!=0){

	  tmp1r = 0.0; tmp1i = 0.0; tmp2r = 0.0; tmp2i = 0.0;

	  for (k1=0; k1<Ngrid1; k1++){  /* Gz-loop */

	    if (k1<Ngrid1/2) sk1 = (double)k1;
	    else             sk1 = (double)(k1 - Ngrid1);

	    Gz  = sk1*rtv[1][1];  /* original Gx -> Gz */
	    Gz2 = Gz*Gz;
	    G2  = Gp2 + Gz2;

	    tmp1r +=  (ReV2[k2][k1][k3]*(Gp*cos(Gz*z0) - Gz*sin(Gz*z0))
		       - ImV2[k2][k1][k3]*(Gp*sin(Gz*z0) + Gz*cos(Gz*z0)))/G2*tmp0;
	    tmp1i +=  (ImV2[k2][k1][k3]*(Gp*cos(Gz*z0) - Gz*sin(Gz*z0))
		       + ReV2[k2][k1][k3]*(Gp*sin(Gz*z0) + Gz*cos(Gz*z0)))/G2*tmp0;

	    tmp2r +=  (ReV2[k2][k1][k3]*(Gp*cos(Gz*z0) - Gz*sin(Gz*z0))
		       + ImV2[k2][k1][k3]*(Gp*sin(Gz*z0) + Gz*cos(Gz*z0)))/G2*tmp0;
	    tmp2i +=  (ImV2[k2][k1][k3]*(Gp*cos(Gz*z0) - Gz*sin(Gz*z0))
		       - ReV2[k2][k1][k3]*(Gp*sin(Gz*z0) + Gz*cos(Gz*z0)))/G2*tmp0;

	    ReV2[k2][k1][k3] = 4.0*PI*ReV2[k2][k1][k3]/G2*tmp0;
	    ImV2[k2][k1][k3] = 4.0*PI*ImV2[k2][k1][k3]/G2*tmp0;

	  } /* end of Gz-loop */

	  One_dim_FFT(fp,2,1,k2,k3,ReV2,ImV2,0.0);  /*FFT for Gz -> z */

	  for (k1=0; k1<Ngrid1; k1++){  /* z-loop */

	    if (k1<Ngrid1/2) sk1 = (double)k1;
	    else             sk1 = (double)(k1 - Ngrid1);

	    zz = sk1*tv[1][1]/(double)Ngrid1;

	    ReV2[k2][k1][k3] += -2.0*PI*exp( Gp*(zz-z0))*tmp1r/Gp
	      -2.0*PI*exp(-Gp*(zz+z0))*tmp2r/Gp;
	    ImV2[k2][k1][k3] += -2.0*PI*exp( Gp*(zz-z0))*tmp1i/Gp
	      -2.0*PI*exp(-Gp*(zz+z0))*tmp2i/Gp;

	  } /* end of z-loop */

	} /* end of if */
      } /* end of Gy-loop */
    } /* end of Gx-loop */

    /* * * *  End of G_|| != 0 case  * * * */


    /* * * *  G_|| = 0  * * * */

    if(myid==Host_ID){

      tmp5r = ReV2[0][0][0];
      tmp5i = ImV2[0][0][0];

      ReV2[0][0][0] = -2.0*PI*(z0*z0)*ReV2[0][0][0]*tmp0;
      ImV2[0][0][0] = -2.0*PI*(z0*z0)*ImV2[0][0][0]*tmp0;

      tmp1r = 0.0; tmp1i = 0.0; tmp2r = 0.0; tmp2i = 0.0; tmp3r = 0.0; tmp3i = 0.0;

      for (k1=1; k1<Ngrid1; k1++){  /* Gz-loop */

	if (k1<Ngrid1/2) sk1 = (double)k1;
	else             sk1 = (double)(k1 - Ngrid1);

	Gz  = sk1*rtv[1][1];  /* original Gx -> Gz */
	Gz2 = Gz*Gz;

	tmp1r += (-ReV2[0][k1][0]*sin(Gz*z0) - ImV2[0][k1][0]*cos(Gz*z0))/Gz*tmp0;
	tmp1i += (-ImV2[0][k1][0]*sin(Gz*z0) + ReV2[0][k1][0]*cos(Gz*z0))/Gz*tmp0;

	tmp2r += ( ReV2[0][k1][0]*sin(Gz*z0) - ImV2[0][k1][0]*cos(Gz*z0))/Gz*tmp0;
	tmp2i += ( ImV2[0][k1][0]*sin(Gz*z0) + ReV2[0][k1][0]*cos(Gz*z0))/Gz*tmp0;

	tmp3r += ReV2[0][k1][0]*cos(Gz*z0)/Gz2*tmp0;
	tmp3i += ImV2[0][k1][0]*cos(Gz*z0)/Gz2*tmp0;

	ReV2[0][k1][0] = 4.0*PI*ReV2[0][k1][0]/Gz2*tmp0;
	ImV2[0][k1][0] = 4.0*PI*ImV2[0][k1][0]/Gz2*tmp0;

      } /* end of Gz-loop */

      One_dim_FFT(fp,2,1,0,0,ReV2,ImV2,0.0);  /* FFT for Gz -> z */

      for (k1=0; k1<Ngrid1; k1++){  /* z-loop */

	if (k1<Ngrid1/2) sk1 = (double)k1;
	else             sk1 = (double)(k1 - Ngrid1);

	zz = sk1*tv[1][1]/(double)Ngrid1;

	ReV2[0][k1][0] += -2.0*PI*(zz*zz)*tmp5r*tmp0
	  -2.0*PI*(zz-z0)*tmp1r    
	  -2.0*PI*(zz+z0)*tmp2r    
	  -4.0*PI*tmp3r;

	ImV2[0][k1][0] += -2.0*PI*(zz*zz)*tmp5i*tmp0
	  -2.0*PI*(zz-z0)*tmp1i    
	  -2.0*PI*(zz+z0)*tmp2i    
	  -4.0*PI*tmp3i;

      } /* end of z-loop */

    } /* end of if myid==Host_ID */

    /* * * *  End of G_|| = 0 case  * * * */

  } /* ESM_switch==1 */  


  /************************************************************/
  /*  "metal|vacuum|metal" boundary condition (ESM_switch==2) */
  /*  "inside-capacitor"   boundary condition (ESM_switch==4) */
  /************************************************************/

  else if (ESM_switch==2 || ESM_switch==4){

    z0 = tv[1][1]/2.0;
    z1 = z0*1.5;  /* z1: position of ideal metal surface */

    /* * * *  G_|| != 0  * * * */

    for (k2=0; k2<My_NGrid2_Poisson; k2++){

      kk2 = k2 + Start_Grid2[myid];

      if (kk2<Ngrid2/2) sk2 = (double)kk2;
      else              sk2 = (double)(kk2 - Ngrid2);

      for (k3=0; k3<Ngrid3; k3++){

	if (k3<Ngrid3/2) sk3 = (double)k3;
	else             sk3 = (double)(k3 - Ngrid3);

	Gx  = sk2*rtv[2][2] + sk3*rtv[3][2];  /* original Gy,Gz -> Gx,Gy */
	Gy  = sk2*rtv[2][3] + sk3*rtv[3][3];
	Gp2 = Gx*Gx + Gy*Gy;
	Gp  = sqrt(Gp2);

	if (kk2!=0 || k3!=0){

	  tmp1r = 1.0; tmp1i = 0.0; tmp2r = 0.0; tmp2i = 0.0; 
	  tmp3r = 0.0; tmp3i = 0.0; tmp4r = 0.0; tmp4i = 0.0;

	  for (k1=0; k1<Ngrid1; k1++){  /* Gz-loop */

	    if (k1<Ngrid1/2) sk1 = (double)k1;
	    else             sk1 = (double)(k1 - Ngrid1);

	    Gz  = sk1*rtv[1][1];  /* original Gx -> Gz */
	    Gz2 = Gz*Gz;
	    G2  = Gp2 + Gz2;

	    tmp1r +=  (ReV2[k2][k1][k3]*(Gp*cos(Gz*z0) - Gz*sin(Gz*z0))
		       - ImV2[k2][k1][k3]*(Gp*sin(Gz*z0) + Gz*cos(Gz*z0)))/G2*exp(Gp*(z1-z0))*tmp0;
	    tmp1i +=  (ImV2[k2][k1][k3]*(Gp*cos(Gz*z0) - Gz*sin(Gz*z0))
		       + ReV2[k2][k1][k3]*(Gp*sin(Gz*z0) + Gz*cos(Gz*z0)))/G2*exp(Gp*(z1-z0))*tmp0;

	    tmp2r +=  (ReV2[k2][k1][k3]*(Gp*cos(Gz*z0) + Gz*sin(Gz*z0))
		       - ImV2[k2][k1][k3]*(Gp*sin(Gz*z0) - Gz*cos(Gz*z0)))/G2*exp(Gp*(z0-z1))*tmp0;
	    tmp2i +=  (ImV2[k2][k1][k3]*(Gp*cos(Gz*z0) + Gz*sin(Gz*z0))
		       + ReV2[k2][k1][k3]*(Gp*sin(Gz*z0) - Gz*cos(Gz*z0)))/G2*exp(Gp*(z0-z1))*tmp0;

	    tmp3r +=  (ReV2[k2][k1][k3]*(Gp*cos(Gz*z0) - Gz*sin(Gz*z0))
		       + ImV2[k2][k1][k3]*(Gp*sin(Gz*z0) + Gz*cos(Gz*z0)))/G2*exp(Gp*(z1-z0))*tmp0;
	    tmp3i +=  (ImV2[k2][k1][k3]*(Gp*cos(Gz*z0) - Gz*sin(Gz*z0))
		       - ReV2[k2][k1][k3]*(Gp*sin(Gz*z0) + Gz*cos(Gz*z0)))/G2*exp(Gp*(z1-z0))*tmp0;

	    tmp4r +=  (ReV2[k2][k1][k3]*(Gp*cos(Gz*z0) + Gz*sin(Gz*z0))
		       + ImV2[k2][k1][k3]*(Gp*sin(Gz*z0) - Gz*cos(Gz*z0)))/G2*exp(Gp*(z0-z1))*tmp0;
	    tmp4i +=  (ImV2[k2][k1][k3]*(Gp*cos(Gz*z0) + Gz*sin(Gz*z0))
		       - ReV2[k2][k1][k3]*(Gp*sin(Gz*z0) - Gz*cos(Gz*z0)))/G2*exp(Gp*(z0-z1))*tmp0;

	    ReV2[k2][k1][k3] = 4.0*PI*ReV2[k2][k1][k3]/G2*tmp0;
	    ImV2[k2][k1][k3] = 4.0*PI*ImV2[k2][k1][k3]/G2*tmp0;

	  } /* end of Gz-loop */

	  One_dim_FFT(fp,2,1,k2,k3,ReV2,ImV2,0.0);  /*FFT for Gz -> z */

	  for (k1=0; k1<Ngrid1; k1++){  /* z-loop */

	    if (k1<Ngrid1/2) sk1 = (double)k1;
	    else             sk1 = (double)(k1 - Ngrid1);

	    zz = sk1*tv[1][1]/(double)Ngrid1;

	    ReV2[k2][k1][k3] += -4.0*PI/Gp/(1.0-exp(-4.0*Gp*z1))*
	      ((exp(Gp*(zz-z1)) - exp(-Gp*(zz+3.0*z1)))*(tmp1r + tmp2r)
	       +(exp(Gp*(zz-3.0*z1)) - exp(-Gp*(zz+z1)))*(tmp3r + tmp4r));
	    ImV2[k2][k1][k3] += -4.0*PI/Gp/(1.0-exp(-4.0*Gp*z1))*
	      ((exp(Gp*(zz-z1)) - exp(-Gp*(zz+3.0*z1)))*(tmp1i + tmp2i)
	       +(exp(Gp*(zz-3.0*z1)) - exp(-Gp*(zz+z1)))*(tmp3i + tmp4i));

	  } /* end of z-loop */

	} /* end of if */
      } /* end of Gy-loop */
    } /* end of Gx-loop */

    /* * * *  End of G_|| != 0 case  * * * */


    /* * * *  G_|| = 0  * * * */

    if(myid==Host_ID){

      tmp5r = ReV2[0][0][0];
      tmp5i = ImV2[0][0][0];

      ReV2[0][0][0] = -2.0*PI*(z0*z0 - 2.0*z0*z1)*ReV2[0][0][0]*tmp0;
      ImV2[0][0][0] = -2.0*PI*(z0*z0 - 2.0*z0*z1)*ImV2[0][0][0]*tmp0;

      tmp1r = 0.0; tmp1i = 0.0; tmp2r = 0.0; tmp2i = 0.0; 
      tmp3r = 0.0; tmp3i = 0.0; tmp4r = 0.0; tmp4i = 0.0;

      for (k1=1; k1<Ngrid1; k1++){  /* Gz-loop */

	if (k1<Ngrid1/2) sk1 = (double)k1;
	else             sk1 = (double)(k1 - Ngrid1);

	Gz  = sk1*rtv[1][1];  /* original Gx -> Gz */
	Gz2 = Gz*Gz;

	tmp1r += ( ReV2[0][k1][0]*cos(Gz*z0) - ImV2[0][k1][0]*sin(Gz*z0))/Gz2*tmp0;
	tmp1i += ( ImV2[0][k1][0]*cos(Gz*z0) + ReV2[0][k1][0]*sin(Gz*z0))/Gz2*tmp0;

	tmp2r += ( ReV2[0][k1][0]*cos(Gz*z0) + ImV2[0][k1][0]*sin(Gz*z0))/Gz2*tmp0;
	tmp2i += ( ImV2[0][k1][0]*cos(Gz*z0) - ReV2[0][k1][0]*sin(Gz*z0))/Gz2*tmp0;

	tmp3r += ( ReV2[0][k1][0]*sin(Gz*z0))/Gz*tmp0;
	tmp3i += ( ImV2[0][k1][0]*sin(Gz*z0))/Gz*tmp0;

	tmp4r += (-ImV2[0][k1][0]*cos(Gz*z0))/Gz*tmp0;
	tmp4i += ( ReV2[0][k1][0]*cos(Gz*z0))/Gz*tmp0;

	ReV2[0][k1][0] = 4.0*PI*ReV2[0][k1][0]/Gz2*tmp0;
	ImV2[0][k1][0] = 4.0*PI*ImV2[0][k1][0]/Gz2*tmp0;

      } /* end of Gz-loop */

      One_dim_FFT(fp,2,1,0,0,ReV2,ImV2,0.0);  /* FFT for Gz -> z */

      for (k1=0; k1<Ngrid1; k1++){  /* z-loop */

	if (k1<Ngrid1/2) sk1 = (double)k1;
	else             sk1 = (double)(k1 - Ngrid1);

	zz = sk1*tv[1][1]/(double)Ngrid1;

	ReV2[0][k1][0] += -2.0*PI*(zz*zz)*tmp5r*tmp0
	  -2.0*PI*(zz+z1)/z1*tmp1r
	  +2.0*PI*(zz-z1)/z1*tmp2r
	  +4.0*PI*(z1-z0)*tmp3r
	  -4.0*PI*(z1-z0)*tmp4r*zz/z1
	  -0.5*V_ESM*(zz-z1)/z1; /* <-- for the case ESM_switch==4 */

	ImV2[0][k1][0] += -2.0*PI*(zz*zz)*tmp5i*tmp0
	  -2.0*PI*(zz+z1)/z1*tmp1i
	  +2.0*PI*(zz-z1)/z1*tmp2i
	  +4.0*PI*(z1-z0)*tmp3i
	  -4.0*PI*(z1-z0)*tmp4i*zz/z1;

      } /* end of z-loop */

    } /* end of if myid==Host_ID */

    /* * * *  End of G_|| = 0 case  * * * */

  } /* ESM_switch==2 or 4 */ 


  /****************************************************************************/
  /*  "vacuum|vacuum|metal" boundary condition (for electrochemical systems)  */
  /****************************************************************************/

  else if (ESM_switch==3){

    z0 = tv[1][1]/2.0;   
    z1 = z0;  /* z1: position of ideal metal surface */

    /* * * *  G_|| != 0  * * * */

    for (k2=0; k2<My_NGrid2_Poisson; k2++){

      kk2 = k2 + Start_Grid2[myid];

      if (kk2<Ngrid2/2) sk2 = (double)kk2;
      else              sk2 = (double)(kk2 - Ngrid2);

      for (k3=0; k3<Ngrid3; k3++){

	if (k3<Ngrid3/2) sk3 = (double)k3;
	else             sk3 = (double)(k3 - Ngrid3);

	Gx  = sk2*rtv[2][2] + sk3*rtv[3][2];  /* original Gy,Gz -> Gx,Gy */
	Gy  = sk2*rtv[2][3] + sk3*rtv[3][3];
	Gp2 = Gx*Gx + Gy*Gy;
	Gp  = sqrt(Gp2);

	if (kk2!=0 || k3!=0){

	  tmp1r = 0.0; tmp1i = 0.0; tmp2r = 0.0; tmp2i = 0.0; tmp3r = 0.0; tmp3i = 0.0;

	  for (k1=0; k1<Ngrid1; k1++){  /* Gz-loop */

	    if (k1<Ngrid1/2) sk1 = (double)k1;
	    else             sk1 = (double)(k1 - Ngrid1);

	    Gz  = sk1*rtv[1][1];  /* original Gx -> Gz */
	    Gz2 = Gz*Gz;
	    G2  = Gp2 + Gz2;

	    tmp1r +=  (ReV2[k2][k1][k3]*(Gp*cos(Gz*z0) - Gz*sin(Gz*z0))
		       - ImV2[k2][k1][k3]*(Gp*sin(Gz*z0) + Gz*cos(Gz*z0)))/G2*exp(Gp*(z1-z0))*tmp0;
	    tmp1i +=  (ImV2[k2][k1][k3]*(Gp*cos(Gz*z0) - Gz*sin(Gz*z0))
		       + ReV2[k2][k1][k3]*(Gp*sin(Gz*z0) + Gz*cos(Gz*z0)))/G2*exp(Gp*(z1-z0))*tmp0;

	    tmp2r +=  (ReV2[k2][k1][k3]*(Gp*cos(Gz*z0) + Gz*sin(Gz*z0))
		       - ImV2[k2][k1][k3]*(Gp*sin(Gz*z0) - Gz*cos(Gz*z0)))/G2*exp(Gp*(z0-z1))*tmp0;
	    tmp2i +=  (ImV2[k2][k1][k3]*(Gp*cos(Gz*z0) + Gz*sin(Gz*z0))
		       + ReV2[k2][k1][k3]*(Gp*sin(Gz*z0) - Gz*cos(Gz*z0)))/G2*exp(Gp*(z0-z1))*tmp0;

	    tmp3r +=  (ReV2[k2][k1][k3]*(Gp*cos(Gz*z0) - Gz*sin(Gz*z0))
		       + ImV2[k2][k1][k3]*(Gp*sin(Gz*z0) + Gz*cos(Gz*z0)))/G2*tmp0;
	    tmp3i +=  (ImV2[k2][k1][k3]*(Gp*cos(Gz*z0) - Gz*sin(Gz*z0))
		       - ReV2[k2][k1][k3]*(Gp*sin(Gz*z0) + Gz*cos(Gz*z0)))/G2*tmp0;

	    ReV2[k2][k1][k3] = 4.0*PI*ReV2[k2][k1][k3]/G2*tmp0;
	    ImV2[k2][k1][k3] = 4.0*PI*ImV2[k2][k1][k3]/G2*tmp0;

	  } /* end of Gz-loop */

	  One_dim_FFT(fp,2,1,k2,k3,ReV2,ImV2,0.0);  /*FFT for Gz -> z */

	  for (k1=0; k1<Ngrid1; k1++){  /* z-loop */

	    if (k1<Ngrid1/2) sk1 = (double)k1;
	    else             sk1 = (double)(k1 - Ngrid1);

	    zz = sk1*tv[1][1]/(double)Ngrid1 /*- Grid_Origin[1]*/;

	    ReV2[k2][k1][k3] += -2.0*PI*exp(Gp*(zz-z1))*(tmp1r + tmp2r)/Gp
	      +2.0*PI*(exp(Gp*(zz-z0-2.0*z1))-exp(-Gp*(zz+z0)))*tmp3r/Gp;
	    ImV2[k2][k1][k3] += -2.0*PI*exp(Gp*(zz-z1))*(tmp1i + tmp2i)/Gp
	      +2.0*PI*(exp(Gp*(zz-z0-2.0*z1))-exp(-Gp*(zz+z0)))*tmp3i/Gp;

	  } /* end of z-loop */

	} /* end of if */
      } /* end of Gy-loop */
    } /* end of Gx-loop */

    /* * * *  End of G_|| != 0 case  * * * */


    /* * * *  G_|| = 0  * * * */

    if(myid==Host_ID){

      tmp5r = ReV2[0][0][0];
      tmp5i = ImV2[0][0][0];

      ReV2[0][0][0] = -2.0*PI*(z0*z0 - 4.0*z0*z1)*ReV2[0][0][0]*tmp0;
      ImV2[0][0][0] = -2.0*PI*(z0*z0 - 4.0*z0*z1)*ImV2[0][0][0]*tmp0;

      tmp1r = 0.0; tmp1i = 0.0; tmp2r = 0.0; tmp2i = 0.0; tmp3r = 0.0; tmp3i = 0.0;

      for (k1=1; k1<Ngrid1; k1++){  /* Gz-loop */

	if (k1<Ngrid1/2) sk1 = (double)k1;
	else             sk1 = (double)(k1 - Ngrid1);

	Gz  = sk1*rtv[1][1];  /* original Gx -> Gz */
	Gz2 = Gz*Gz;

	tmp1r += ( ReV2[0][k1][0]*cos(Gz*z0) - ImV2[0][k1][0]*sin(Gz*z0))/Gz2*tmp0;
	tmp1i += ( ImV2[0][k1][0]*cos(Gz*z0) + ReV2[0][k1][0]*sin(Gz*z0))/Gz2*tmp0;

	tmp2r += ( ReV2[0][k1][0]*sin(Gz*z0) - ImV2[0][k1][0]*cos(Gz*z0))/Gz*tmp0;
	tmp2i += ( ImV2[0][k1][0]*sin(Gz*z0) + ReV2[0][k1][0]*cos(Gz*z0))/Gz*tmp0;

	tmp3r += (-ReV2[0][k1][0]*sin(Gz*z0) - ImV2[0][k1][0]*cos(Gz*z0))/Gz*tmp0;
	tmp3i += (-ImV2[0][k1][0]*sin(Gz*z0) + ReV2[0][k1][0]*cos(Gz*z0))/Gz*tmp0;

	ReV2[0][k1][0] = 4.0*PI*ReV2[0][k1][0]/Gz2*tmp0;
	ImV2[0][k1][0] = 4.0*PI*ImV2[0][k1][0]/Gz2*tmp0;

      } /* end of Gz-loop */

      One_dim_FFT(fp,2,1,0,0,ReV2,ImV2,0.0);  /* FFT for Gz -> z */

      for (k1=0; k1<Ngrid1; k1++){  /* z-loop */

	if (k1<Ngrid1/2) sk1 = (double)k1;
	else             sk1 = (double)(k1 - Ngrid1);

	zz = sk1*tv[1][1]/(double)Ngrid1 /*- Grid_Origin[1]*/;

	ReV2[0][k1][0] += -2.0*PI*(zz*zz + 2.0*z0*zz)*tmp5r*tmp0
	  -4.0*PI*tmp1r
	  -4.0*PI*(zz-z1)*tmp2r
	  -4.0*PI*(z1-z0)*tmp3r;

	ImV2[0][k1][0] += -2.0*PI*(zz*zz + 2.0*z0*zz)*tmp5i*tmp0
	  -4.0*PI*tmp1i 
	  -4.0*PI*(zz-z1)*tmp2i
	  -4.0*PI*(z1-z0)*tmp3i;

      } /* end of z-loop */

    } /* end of if myid==Host_ID */

    /* * * *  End of G_|| = 0 case  * * * */

  } /* ESM_switch==3 */ 

/*
  if(myid==Host_ID){
    printf(" $$$ check for delta-V_H(G_||=0,z) (3) $$$ \n");
    for (kz1=0; kz1<Ngrid1; kz1++){
      printf("## ReV2,ImV2[0][%4d][0]= %16.9f %16.9f \n",kz1,ReV2[0][kz1][0],ImV2[0][kz1][0]);
    }
  }
*/

  /*************************************************************************
     V(G_||,z') -> V(G_||,z) change the order of z-coordinate mesh points  
  *************************************************************************/

    tmp2 = tv[1][1]/2.0; 

    for (k2=0; k2<My_NGrid2_Poisson; k2++){
      for (k3=0; k3<Ngrid3; k3++){

	One_dim_FFT(fp,2,-1,k2,k3,ReV2,ImV2,0.0);

	for (k1=0; k1<Ngrid1; k1++){

	  if (k1<Ngrid1/2) sk1 = (double)k1;
	  else             sk1 = (double)(k1 - Ngrid1);

	  Gz  = sk1*rtv[1][1];

	  tmp1r = ReV2[k2][k1][k3]*cos(Gz*tmp2) - ImV2[k2][k1][k3]*sin(Gz*tmp2);
	  tmp1i = ImV2[k2][k1][k3]*cos(Gz*tmp2) + ReV2[k2][k1][k3]*sin(Gz*tmp2);

	  ReV2[k2][k1][k3] = tmp1r/(double)Ngrid1;
	  ImV2[k2][k1][k3] = tmp1i/(double)Ngrid1;
  
	}
  
      }
    }

    if (myid==Host_ID){

      sprintf(fname,"%s%s.ESM.dhart",filepath,filename);

      if ( (fp=fopen(fname,"w"))!=NULL ) {

        fprintf(fp,"## check of difference Hartree potential (Hartree) ##\n");fflush(stdout);
        fprintf(fp,"## Grid : x-corrdinate (Ang) : delta-V_H(G_||=0,z) ##\n");fflush(stdout);

        One_dim_FFT(fp,1,1,0,0,ReV2,ImV2,1.0); 

        fprintf(fp,"\n");
        fclose(fp);
      }
      else {
        printf("Failure in saving *.ESM.dhart.\n");
      }
    }

  /****************************************************
        find the Hartree potential in real space
  ****************************************************/
  
  Get_Value_inReal(0,ReV2,ImV2,ReV1,ImV1,dVHart_Grid,dVHart_Grid); 

/*
  if(myid==Host_ID){
    printf(" $$$ check for delta-V_H(r_||=0,z) (5) $$$ \n");
    for (kz1=0; kz1<My_NGrid1_Poisson; kz1++){
      printf("## ReV1,ImV1[%4d][0][0]= %16.9f %16.9f \n",kz1,ReV1[kz1][0][0],ImV1[kz1][0][0]);
    }
  }
*/

  /* for time */
  MPI_Barrier(mpi_comm_level1);
  dtime(&TEtime);
  time0 = TEtime - TStime;
  return time0;

}



void One_dim_FFT(FILE *fp, 
                 int sgn1, int sgn2, 
                 int kn2,  int kn3, 
                 double ***ReF2, double ***ImF2,
                 double prefac)
{
  int k1,n1;
  double xcord,tv11 = tv[1][1];

  /**************************************************************

      One_dim_FFT(fp,sgn1,sgn2,kn2,kn3,ReF3,ImF3,prefac)

       sgn1= 1: for check of output (ReV & ImV are not changed)
           = 2: for calculation (ReV & ImV are changed)
       sgn2= 1: FFT
           =-1: inverse FFT
       kn2,kn3: Gx,Gy
       prefac: prefactor of output data

  ***************************************************************/

  fftw_complex *in, *out;
  fftw_plan p;

  /****************************************************
    allocation of arrays:

    fftw_complex  in[List_YOUSO[17]];
    fftw_complex out[List_YOUSO[17]];
  ****************************************************/

#ifdef fftw2
  in  = (fftw_complex*)malloc(sizeof(fftw_complex)*List_YOUSO[17]);
  out = (fftw_complex*)malloc(sizeof(fftw_complex)*List_YOUSO[17]);
#else
  in  = fftw_malloc(sizeof(fftw_complex)*List_YOUSO[17]);
  out = fftw_malloc(sizeof(fftw_complex)*List_YOUSO[17]);
#endif


#ifdef fftw2
  p = fftw_create_plan(Ngrid1, sgn2, FFTW_ESTIMATE);
#else
  p = fftw_plan_dft_1d(Ngrid1, in, out, sgn2, FFTW_ESTIMATE);
#endif


  for (k1=0; k1<Ngrid1; k1++){

#ifdef fftw2
    c_re(in[k1]) = ReF2[kn2][k1][kn3];
    c_im(in[k1]) = ImF2[kn2][k1][kn3];
#else
    in[k1][0] = ReF2[kn2][k1][kn3];
    in[k1][1] = ImF2[kn2][k1][kn3];
#endif

  }

#ifdef fftw2
  fftw_one(p, in, out);
#else
  fftw_execute(p);
#endif

  for (n1=0; n1<Ngrid1; n1++){

    if (sgn1 == 1){

      xcord = tv11 / Ngrid1 * n1 * 0.529177249;

#ifdef fftw2
      fprintf(fp,"   %4d  %12.9f  %12.9f \n",n1,xcord,c_re(out[n1])*prefac);
#else
      fprintf(fp,"   %4d  %12.9f  %12.9f \n",n1,xcord,out[n1][0]*prefac);
#endif

    }
    else if (sgn1 == 2){
#ifdef fftw2
      ReF2[kn2][n1][kn3] = c_re(out[n1]);
      ImF2[kn2][n1][kn3] = c_im(out[n1]);
#else
      ReF2[kn2][n1][kn3] = out[n1][0];
      ImF2[kn2][n1][kn3] = out[n1][1];
#endif
    }

  }

  fftw_destroy_plan(p);
  fftw_cleanup();

  /****************************************************
    freeing of arrays:

    fftw_complex  in[List_YOUSO[17]];
    fftw_complex out[List_YOUSO[17]];
  ****************************************************/

#ifdef fftw2
  free(in);
  free(out);
#else
  fftw_free(in);
  fftw_free(out);
#endif

}



