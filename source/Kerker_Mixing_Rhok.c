/**********************************************************************
  Kerker_Mixing_Rhok.c:

     Kerker_Mixing_Rhok.c is a subroutine to achieve self-consistent
     field using the Kerker mixing in k-space. 

  Log of Kerker_Mixing_Rhok.c

     30/Dec/2004  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "openmx_common.h"
#include "tran_prototypes.h"
#include "tran_variables.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

#ifdef noomp
#include "mimic_omp.h"
#else
#include <omp.h>
#endif

#define  maxima_step  1000000.0


void Kerker_Mixing_Rhok_Normal(int Change_switch,
			       double Mix_wgt,
			       double *****ReRhok,
			       double *****ImRhok,
			       double ****ReBestRhok,
			       double ****ImBestRhok,
			       double *****Residual_ReRhok,
			       double *****Residual_ImRhok,
			       double ***ReV1,
			       double ***ImV1,
			       double ***ReV2,
			       double ***ImV2,
			       double ***ReRhoAtomk,
			       double ***ImRhoAtomk);

void Kerker_Mixing_Rhok_NEGF(int Change_switch,
			     double Mix_wgt,
			     double *****ReRhok,
			     double *****ImRhok,
			     double ****ReBestRhok,
			     double ****ImBestRhok,
			     double *****Residual_ReRhok,
			     double *****Residual_ImRhok,
			     double ***ReV1,
			     double ***ImV1,
			     double ***ReV2,
			     double ***ImV2,
			     double ***ReRhoAtomk,
			     double ***ImRhoAtomk);



void Kerker_Mixing_Rhok(int Change_switch,
                        double Mix_wgt,
                        double *****ReRhok,
                        double *****ImRhok,
                        double ****ReBestRhok,
                        double ****ImBestRhok,
                        double *****Residual_ReRhok,
                        double *****Residual_ImRhok,
                        double ***ReV1,
                        double ***ImV1,
                        double ***ReV2,
                        double ***ImV2,
                        double ***ReRhoAtomk,
                        double ***ImRhoAtomk)
{

  if (Solver!=4 || TRAN_Poisson_flag==4){

    Kerker_Mixing_Rhok_Normal(Change_switch,
			      Mix_wgt,
			      ReRhok,
			      ImRhok,
			      ReBestRhok,
			      ImBestRhok,
			      Residual_ReRhok,
			      Residual_ImRhok,
			      ReV1,
			      ImV1,
			      ReV2,
			      ImV2,
			      ReRhoAtomk,
			      ImRhoAtomk);
  }

  else if (Solver==4){

    Kerker_Mixing_Rhok_NEGF(Change_switch,
			    Mix_wgt,
			    ReRhok,
			    ImRhok,
			    ReBestRhok,
			    ImBestRhok,
			    Residual_ReRhok,
			    Residual_ImRhok,
			    ReV1,
			    ImV1,
			    ReV2,
			    ImV2,
			    ReRhoAtomk,
			    ImRhoAtomk);
  }

}










void Kerker_Mixing_Rhok_Normal(int Change_switch,
			       double Mix_wgt,
			       double *****ReRhok,
			       double *****ImRhok,
			       double ****ReBestRhok,
			       double ****ImBestRhok,
			       double *****Residual_ReRhok,
			       double *****Residual_ImRhok,
			       double ***ReV1,
			       double ***ImV1,
			       double ***ReV2,
			       double ***ImV2,
			       double ***ReRhoAtomk,
			       double ***ImRhoAtomk)
{
  static int firsttime=1;
  int ian,jan,Mc_AN,Gc_AN,spin,spinmax;
  int h_AN,Gh_AN,m,n,i,j,k,k1,k2,k3,kk2;
  int knum,knum_full,MN,pSCF_iter;
  double Mix_wgt2,Norm,My_Norm,tmp0,tmp1;
  double Min_Weight,Max_Weight,wgt0,wgt1;
  double Gx,Gy,Gz,G2,size_Kerker_weight;
  double sk1,sk2,sk3,G12,G22,G32;
  double G0,G02,weight;
  int numprocs,myid,ID;
  double ***Kerker_weight;
  /* for OpenMP */
  int OMPID,Nthrds,Nprocs,Nloop,Nthrds0;
  double *My_Norm_threads;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /* find an optimum G0 */

  G12 = rtv[1][1]*rtv[1][1] + rtv[1][2]*rtv[1][2] + rtv[1][3]*rtv[1][3]; 
  G22 = rtv[2][1]*rtv[2][1] + rtv[2][2]*rtv[2][2] + rtv[2][3]*rtv[2][3]; 
  G32 = rtv[3][1]*rtv[3][1] + rtv[3][2]*rtv[3][2] + rtv[3][3]*rtv[3][3]; 

  if (G12<G22) G0 = G12;
  else         G0 = G22;
  if (G32<G0)  G0 = G32;

  if (Change_switch==0) G0 = 0.0;
  else                  G0 = sqrt(G0);

  G02 = Kerker_factor*Kerker_factor*G0*G0;

  if      (SpinP_switch==0)  spinmax = 1;
  else if (SpinP_switch==1)  spinmax = 2;
  else if (SpinP_switch==3)  spinmax = 3;

  /***********************************
            set Kerker_weight 
  ************************************/

  size_Kerker_weight = My_NGrid2_Poisson*Ngrid1*Ngrid3;
  knum = My_NGrid2_Poisson*Ngrid1*Ngrid3;
  knum_full = spinmax*My_NGrid2_Poisson*Ngrid1*Ngrid3;

  Kerker_weight = (double***)malloc(sizeof(double**)*My_NGrid2_Poisson); 
  for (i=0; i<My_NGrid2_Poisson; i++){
    Kerker_weight[i] = (double**)malloc(sizeof(double*)*Ngrid1); 
    for (j=0; j<Ngrid1; j++){
      Kerker_weight[i][j] = (double*)malloc(sizeof(double)*Ngrid3); 
    }
  }

  if (firsttime)
  PrintMemory("Kerker_Mixing_Rhok: Kerker_weight",sizeof(double)*size_Kerker_weight,NULL);
  firsttime=0;

#pragma omp parallel shared(myid,Start_Grid2,knum,G02,Kerker_weight,rtv,Ngrid3,Ngrid2,Ngrid1) private(OMPID,Nthrds,Nprocs,k2,kk2,sk2,k,k1,k3,sk1,sk3,Gx,Gy,Gz,G2,weight)
  {

    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    /* one-dimensionalized loop */

    for (k=OMPID*knum/Nthrds; k<(OMPID+1)*knum/Nthrds; k++){

      /* get k2, k1, and k3 */
 
      k2 = k/(Ngrid1*Ngrid3);
      k1 = (k - k2*(Ngrid1*Ngrid3))/Ngrid3;
      k3 = k - k2*(Ngrid1*Ngrid3) - k1*Ngrid3;

      kk2 = k2 + Start_Grid2[myid];

      if (kk2<Ngrid2/2) sk2 = (double)kk2;
      else              sk2 = (double)(kk2 - Ngrid2);

      if (k1<Ngrid1/2)  sk1 = (double)k1;
      else              sk1 = (double)(k1 - Ngrid1);

      if (k3<Ngrid3/2)  sk3 = (double)k3;
      else              sk3 = (double)(k3 - Ngrid3);

      Gx = sk1*rtv[1][1] + sk2*rtv[2][1] + sk3*rtv[3][1];
      Gy = sk1*rtv[1][2] + sk2*rtv[2][2] + sk3*rtv[3][2]; 
      Gz = sk1*rtv[1][3] + sk2*rtv[2][3] + sk3*rtv[3][3];
      G2 = Gx*Gx + Gy*Gy + Gz*Gz;

      if (k1==0 && kk2==0 && k3==0)  weight = 1.0;
      else                           weight = (G2 + G02)/G2;

      Kerker_weight[k2][k1][k3] = weight;

    } /* k */
  } /* #pragma omp parallel */

  /* start... */

  Min_Weight = Min_Mixing_weight;
  if (SCF_RENZOKU==-1){
    Max_Weight = Max_Mixing_weight;
    Max_Mixing_weight2 = Max_Mixing_weight;
  }
  else if (SCF_RENZOKU==1000){  /* past 3 */
    Max_Mixing_weight2 = 2.0*Max_Mixing_weight2;
    if (0.7<Max_Mixing_weight2) Max_Mixing_weight2 = 0.7;
    Max_Weight = Max_Mixing_weight2;
    SCF_RENZOKU = 0;
  }
  else{
    Max_Weight = Max_Mixing_weight2;
  }

  /****************************************************
       norm of residual charge density in k-space
  ****************************************************/

  /* get Nthrds0 */  
#pragma omp parallel shared(Nthrds0)
  {
    Nthrds0 = omp_get_num_threads();
  }
  /* allocation of array */
  My_Norm_threads = (double*)malloc(sizeof(double)*Nthrds0);
 for (Nloop=0; Nloop<Nthrds0; Nloop++) My_Norm_threads[Nloop] = 0.0;

#pragma omp parallel shared(Residual_ReRhok,Residual_ImRhok,ImRhok,ReRhok,Kerker_weight,My_Norm_threads,knum_full,My_NGrid2_Poisson,Ngrid1,Ngrid3) private(OMPID,Nthrds,Nprocs,k,spin,k2,k1,k3,Nloop,weight,tmp0,tmp1)
  {
    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    /* one-dimensionalized loop */

    for (k=OMPID*knum_full/Nthrds; k<(OMPID+1)*knum_full/Nthrds; k++){

      /* get spin, k2, k1, and k3 */

      spin = k/(My_NGrid2_Poisson*Ngrid1*Ngrid3);
      k2 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3))/(Ngrid1*Ngrid3);
      k1 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3))/Ngrid3;  
      k3 = k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3) - k1*Ngrid3;

      weight = Kerker_weight[k2][k1][k3];
      tmp0 = ReRhok[0][spin][k2][k1][k3] - ReRhok[1][spin][k2][k1][k3];
      tmp1 = ImRhok[0][spin][k2][k1][k3] - ImRhok[1][spin][k2][k1][k3];
      Residual_ReRhok[1][spin][k2][k1][k3] = tmp0;
      Residual_ImRhok[1][spin][k2][k1][k3] = tmp1;
      My_Norm_threads[OMPID] += (tmp0*tmp0 + tmp1*tmp1)*weight;

    } /* k */
  } /* #pragma omp parallel */

  /* sum of My_Norm_threads */
  My_Norm = 0.0;
  for (Nloop=0; Nloop<Nthrds0; Nloop++){
    My_Norm += My_Norm_threads[Nloop];
  }

  /* freeing of array */
  free(My_Norm_threads);

  /****************************************************
    MPI: 

    My_Norm
  ****************************************************/

  MPI_Allreduce(&My_Norm, &Norm, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

  /* normalization by the number of grids */
  Norm = Norm/(double)Ngrid1/(double)Ngrid2/(double)Ngrid3;

  /****************************************************
    find an optimum mixing weight
  ****************************************************/

  for (i=4; 1<=i; i--){
    NormRD[i] = NormRD[i-1];
    History_Uele[i] = History_Uele[i-1];
  }
  NormRD[0] = Norm;
  History_Uele[0] = Uele;

  if (Change_switch==1){

    if ((int)sgn(History_Uele[0]-History_Uele[1])
	  ==(int)sgn(History_Uele[1]-History_Uele[2])
       && NormRD[0]<NormRD[1]){

      /* tmp0 = 1.6*Mixing_weight; */

      tmp0 = NormRD[1]/(largest(NormRD[1]-NormRD[0],10e-10))*Mixing_weight;

      if (tmp0<Max_Weight){
        if (Min_Weight<tmp0){
          Mixing_weight = tmp0;
	}
        else{ 
          Mixing_weight = Min_Weight;
	}
      }
      else{ 
        Mixing_weight = Max_Weight;
        SCF_RENZOKU++;  
      }
    }
   
    else if ((int)sgn(History_Uele[0]-History_Uele[1])
             ==(int)sgn(History_Uele[1]-History_Uele[2])
             && NormRD[1]<NormRD[0]){

      tmp0 = NormRD[1]/(largest(NormRD[1]+NormRD[0],10e-10))*Mixing_weight;

      /* tmp0 = Mixing_weight/1.6; */

      if (tmp0<Max_Weight){
        if (Min_Weight<tmp0)
          Mixing_weight = tmp0;
        else 
          Mixing_weight = Min_Weight;
      }
      else 
        Mixing_weight = Max_Weight;

      SCF_RENZOKU = -1;  
    }

    else if ((int)sgn(History_Uele[0]-History_Uele[1])
	     !=(int)sgn(History_Uele[1]-History_Uele[2])
             && NormRD[0]<NormRD[1]){

      /* tmp0 = Mixing_weight*1.2; */

      tmp0 = NormRD[1]/(largest(NormRD[1]-NormRD[0],10e-10))*Mixing_weight;

      if (tmp0<Max_Weight){
        if (Min_Weight<tmp0)
          Mixing_weight = tmp0;
        else 
          Mixing_weight = Min_Weight;
      }
      else{ 
        Mixing_weight = Max_Weight;
        SCF_RENZOKU++;
      }
    }

    else if ((int)sgn(History_Uele[0]-History_Uele[1])
	     !=(int)sgn(History_Uele[1]-History_Uele[2])
             && NormRD[1]<NormRD[0]){

      /* tmp0 = Mixing_weight/2.0; */

      tmp0 = NormRD[1]/(largest(NormRD[1]+NormRD[0],10e-10))*Mixing_weight;

      if (tmp0<Max_Weight){
        if (Min_Weight<tmp0)
          Mixing_weight = tmp0;
        else 
          Mixing_weight = Min_Weight;
      }
      else 
        Mixing_weight = Max_Weight;

      SCF_RENZOKU = -1;
    }

    Mix_wgt = Mixing_weight;
  }

  /****************************************************
                        Mixing
  ****************************************************/

#pragma omp parallel shared(ImRhok,ReRhok,Mix_wgt,Kerker_weight,Ngrid3,Ngrid1,My_NGrid2_Poisson,knum_full) private(OMPID,Nthrds,Nprocs,k,spin,k2,k1,k3,weight,wgt0,wgt1)
  {
    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    /* one-dimensionalized loop */

    for (k=OMPID*knum_full/Nthrds; k<(OMPID+1)*knum_full/Nthrds; k++){

      /* get spin, k2, k1, and k3 */

      spin = k/(My_NGrid2_Poisson*Ngrid1*Ngrid3);
      k2 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3))/(Ngrid1*Ngrid3);
      k1 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3))/Ngrid3;  
      k3 = k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3) - k1*Ngrid3;

      weight = 1.0/Kerker_weight[k2][k1][k3];
      wgt0  = Mix_wgt*weight;
      wgt1 =  1.0 - wgt0;

      ReRhok[0][spin][k2][k1][k3] = wgt0*ReRhok[0][spin][k2][k1][k3]
 	                          + wgt1*ReRhok[1][spin][k2][k1][k3];

      ImRhok[0][spin][k2][k1][k3] = wgt0*ImRhok[0][spin][k2][k1][k3]
	                          + wgt1*ImRhok[1][spin][k2][k1][k3];

      /* correction to largely changing components */

      tmp0 = ReRhok[0][spin][k2][k1][k3] - ReRhok[1][spin][k2][k1][k3];  
      tmp1 = ImRhok[0][spin][k2][k1][k3] - ImRhok[1][spin][k2][k1][k3];  

      if ( maxima_step<(fabs(tmp0)+fabs(tmp1)) ){
	ReRhok[0][spin][k2][k1][k3] = sgn(tmp0)*maxima_step + ReRhok[1][spin][k2][k1][k3]; 
	ImRhok[0][spin][k2][k1][k3] = sgn(tmp1)*maxima_step + ImRhok[1][spin][k2][k1][k3]; 
      }

    } /* k */
  } /* #pragma omp parallel */

  /****************************************************
                         shift of rho
  ****************************************************/

  for (pSCF_iter=(List_YOUSO[38]-1); 0<pSCF_iter; pSCF_iter--){
    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
	for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){
	    ReRhok[pSCF_iter][spin][k2][k1][k3] = ReRhok[pSCF_iter-1][spin][k2][k1][k3]; 
	    ImRhok[pSCF_iter][spin][k2][k1][k3] = ImRhok[pSCF_iter-1][spin][k2][k1][k3]; 
	  }
	}
      }
    }
  }

  /****************************************************
                    shift of residual rho
  ****************************************************/

  for (pSCF_iter=(List_YOUSO[38]-1); 0<pSCF_iter; pSCF_iter--){
    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
	for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){
	    Residual_ReRhok[pSCF_iter][spin][k2][k1][k3] = Residual_ReRhok[pSCF_iter-1][spin][k2][k1][k3];
	    Residual_ImRhok[pSCF_iter][spin][k2][k1][k3] = Residual_ImRhok[pSCF_iter-1][spin][k2][k1][k3]; 
	  }
	}
      }
    }
  }

  /****************************************************
        find the charge density in real space 
  ****************************************************/

  tmp0 = 1.0/(double)Ngrid1/(double)Ngrid2/(double)Ngrid3;

  for (spin=0; spin<spinmax; spin++){

    for (k2=0; k2<My_NGrid2_Poisson; k2++){
      for (k1=0; k1<Ngrid1; k1++){
        for (k3=0; k3<Ngrid3; k3++){
          ReV2[k2][k1][k3] = ReRhok[0][spin][k2][k1][k3];
          ImV2[k2][k1][k3] = ImRhok[0][spin][k2][k1][k3];
        }
      }
    }

    if (spin==0 || spin==1){
      Get_Value_inReal(0,ReV2,ImV2,ReV1,ImV1,Density_Grid[spin],Density_Grid[spin]);

#pragma omp parallel shared(tmp0,spin,Density_Grid) private(OMPID,Nthrds,Nprocs,MN)
      {
	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	for (MN=OMPID*My_NumGrid1/Nthrds; MN<(OMPID+1)*My_NumGrid1/Nthrds; MN++){
	  Density_Grid[spin][MN] = Density_Grid[spin][MN]*tmp0;
	}

      } /* #pragma omp parallel */

    }

    else if (spin==2){
      Get_Value_inReal(1,ReV2,ImV2,ReV1,ImV1,Density_Grid[2],Density_Grid[3]);


#pragma omp parallel shared(tmp0,spin,Density_Grid) private(OMPID,Nthrds,Nprocs,MN)
      {
	/* get info. on OpenMP */ 

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	for (MN=OMPID*My_NumGrid1/Nthrds; MN<(OMPID+1)*My_NumGrid1/Nthrds; MN++){
	  Density_Grid[2][MN] = Density_Grid[2][MN]*tmp0;
	  Density_Grid[3][MN] = Density_Grid[3][MN]*tmp0;
	}

      } /* #pragma omp parallel */

    }
  }

  if (SpinP_switch==0){
    for (MN=0; MN<My_NumGrid1; MN++){
      Density_Grid[1][MN] = Density_Grid[0][MN];
    }
  }

  /****************************************************
     set ReV2 and ImV2 which are used in Poisson.c
  ****************************************************/

  if (SpinP_switch==0){

#pragma omp parallel shared(ImRhoAtomk,ReRhoAtomk,ImRhok,ReRhok,ImV2,ReV2,Ngrid3,Ngrid1,knum) private(OMPID,Nthrds,Nprocs,k,k2,k1,k3)
    {

      /* get info. on OpenMP */ 

      OMPID = omp_get_thread_num();
      Nthrds = omp_get_num_threads();
      Nprocs = omp_get_num_procs();

      /* one-dimensionalized loop */

      for (k=OMPID*knum/Nthrds; k<(OMPID+1)*knum/Nthrds; k++){

	/* get k2, k1, and k3 */
 
	k2 = k/(Ngrid1*Ngrid3);
	k1 = (k - k2*(Ngrid1*Ngrid3))/Ngrid3;
	k3 = k - k2*(Ngrid1*Ngrid3) - k1*Ngrid3;

	ReV2[k2][k1][k3] = 2.0*ReRhok[0][0][k2][k1][k3] - ReRhoAtomk[k2][k1][k3];
	ImV2[k2][k1][k3] = 2.0*ImRhok[0][0][k2][k1][k3] - ImRhoAtomk[k2][k1][k3];
      }

    } /* #pragma omp parallel */

  }
  
  else {

#pragma omp parallel shared(ImRhoAtomk,ReRhoAtomk,ImRhok,ReRhok,ImV2,ReV2,Ngrid3,Ngrid1,knum) private(OMPID,Nthrds,Nprocs,k,k2,k1,k3)
    {

      /* get info. on OpenMP */ 

      OMPID = omp_get_thread_num();
      Nthrds = omp_get_num_threads();
      Nprocs = omp_get_num_procs();

      /* one-dimensionalized loop */

      for (k=OMPID*knum/Nthrds; k<(OMPID+1)*knum/Nthrds; k++){

	/* get k2, k1, and k3 */
 
	k2 = k/(Ngrid1*Ngrid3);
	k1 = (k - k2*(Ngrid1*Ngrid3))/Ngrid3;
	k3 = k - k2*(Ngrid1*Ngrid3) - k1*Ngrid3;

	ReV2[k2][k1][k3] = ReRhok[0][0][k2][k1][k3] + ReRhok[0][1][k2][k1][k3] - ReRhoAtomk[k2][k1][k3];
	ImV2[k2][k1][k3] = ImRhok[0][0][k2][k1][k3] + ImRhok[0][1][k2][k1][k3] - ImRhoAtomk[k2][k1][k3];

      } /* k */
    } /* #pragma omp parallel */
  }

  /****************************************************
          store the best input charge density 
  ****************************************************/

  if (Change_switch==0 || NormRD[0]<BestNormRD ){

    BestNormRD = NormRD[0];

    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
	for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){
            ReBestRhok[spin][k2][k1][k3] = ReRhok[2][spin][k2][k1][k3];
            ImBestRhok[spin][k2][k1][k3] = ImRhok[2][spin][k2][k1][k3];
	  }
	}
      }
    }
  }

  /* freeing of arrays */

  for (i=0; i<My_NGrid2_Poisson; i++){
    for (j=0; j<Ngrid1; j++){
      free(Kerker_weight[i][j]);
    }
    free(Kerker_weight[i]);
  }
  free(Kerker_weight);

}



void Kerker_Mixing_Rhok_NEGF(int Change_switch,
			     double Mix_wgt,
			     double *****ReRhok,
			     double *****ImRhok,
			     double ****ReBestRhok,
			     double ****ImBestRhok,
			     double *****Residual_ReRhok,
			     double *****Residual_ImRhok,
			     double ***ReV1,
			     double ***ImV1,
			     double ***ReV2,
			     double ***ImV2,
			     double ***ReRhoAtomk,
			     double ***ImRhoAtomk)
{
  static int firsttime=1;
  int ian,jan,Mc_AN,Gc_AN,spin,spinmax;
  int h_AN,Gh_AN,m,n,n1,i,j,k,k1,k2,k3,kk2;
  int knum,knum_full,MN,pSCF_iter;
  double Mix_wgt2,Norm,My_Norm,tmp0,tmp1;
  double Min_Weight,Max_Weight,wgt0,wgt1;
  double Gx,Gy,Gz,G2,size_Kerker_weight;
  double sk1,sk2,sk3,G12,G22,G32;
  double G0,G02,G02p,weight;
  double Q,sumL,sumR;
  int numprocs,myid,ID;
  double ***Kerker_weight;
  /* for OpenMP */
  int OMPID,Nthrds,Nprocs,Nloop,Nthrds0;
  double *My_Norm_threads;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /* find an optimum G0 */

  G22 = rtv[2][1]*rtv[2][1] + rtv[2][2]*rtv[2][2] + rtv[2][3]*rtv[2][3]; 
  G32 = rtv[3][1]*rtv[3][1] + rtv[3][2]*rtv[3][2] + rtv[3][3]*rtv[3][3]; 

  if (G22<G32) G0 = G22;
  else         G0 = G32;

  if (Change_switch==0) G0 = 0.0;
  else                  G0 = sqrt(G0);

  G02 = Kerker_factor*Kerker_factor*G0*G0;
  G02p = (0.3*Kerker_factor*G0)*(0.3*Kerker_factor*G0);

  if      (SpinP_switch==0)  spinmax = 1;
  else if (SpinP_switch==1)  spinmax = 2;
  else if (SpinP_switch==3)  spinmax = 3;

  /***********************************
            set Kerker_weight 
  ************************************/

  size_Kerker_weight = My_NGrid2_Poisson*Ngrid1*Ngrid3;
  knum = My_NGrid2_Poisson*Ngrid1*Ngrid3;
  knum_full = spinmax*My_NGrid2_Poisson*Ngrid1*Ngrid3;

  Kerker_weight = (double***)malloc(sizeof(double**)*My_NGrid2_Poisson); 
  for (i=0; i<My_NGrid2_Poisson; i++){
    Kerker_weight[i] = (double**)malloc(sizeof(double*)*Ngrid1); 
    for (j=0; j<Ngrid1; j++){
      Kerker_weight[i][j] = (double*)malloc(sizeof(double)*Ngrid3); 
    }
  }
  
  if (firsttime)
  PrintMemory("Kerker_Mixing_Rhok: Kerker_weight",sizeof(double)*size_Kerker_weight,NULL);
  firsttime=0;

#pragma omp parallel shared(myid,Start_Grid2,knum,G02,Kerker_weight,rtv,Ngrid3,Ngrid2,Ngrid1) private(OMPID,Nthrds,Nprocs,k2,kk2,sk2,k,k1,k3,sk1,sk3,Gx,Gy,Gz,G2,weight)
  {

    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    /* one-dimensionalized loop */

    for (k=OMPID*knum/Nthrds; k<(OMPID+1)*knum/Nthrds; k++){

      /* get k2, k1, and k3 */
 
      k2 = k/(Ngrid1*Ngrid3);
      k1 = (k - k2*(Ngrid1*Ngrid3))/Ngrid3;
      k3 = k - k2*(Ngrid1*Ngrid3) - k1*Ngrid3;

      kk2 = k2 + Start_Grid2[myid];

      if (kk2<Ngrid2/2) sk2 = (double)kk2;
      else              sk2 = (double)(kk2 - Ngrid2);

      if (k3<Ngrid3/2)  sk3 = (double)k3;
      else              sk3 = (double)(k3 - Ngrid3);

      Gx = sk2*rtv[2][1] + sk3*rtv[3][1];
      Gy = sk2*rtv[2][2] + sk3*rtv[3][2]; 
      Gz = sk2*rtv[2][3] + sk3*rtv[3][3];
      G2 = Gx*Gx + Gy*Gy + Gz*Gz;

      if (Change_switch==0){
        weight = (G2 + G02 + 0.1)/(G2 + 0.01);
      }
      else{

        if (kk2==0 && k3==0){

	  sumL = 0.0;
          sumR = 0.0; 

          if (SpinP_switch==0){

            for (n1=0; n1<k1; n1++){
              sumL += 2.0*ReRhok[0][0][k2][n1][k3];
            }

            for (n1=k1+1; n1<Ngrid1; n1++){
              sumR += 2.0*ReRhok[0][0][k2][n1][k3];
            }
	  }
        
          else if (SpinP_switch==1){

            for (n1=0; n1<k1; n1++){
              sumL += ReRhok[0][0][k2][n1][k3] + ReRhok[0][1][k2][n1][k3];
            }

            for (n1=k1+1; n1<Ngrid1; n1++){
              sumR += ReRhok[0][0][k2][n1][k3] + ReRhok[0][1][k2][n1][k3];
            }
	  }

          Q = 4.0*fabs(sumL - sumR)*GridVol + 1.0;

	}    
        else {
          Q = 1.0;
        }

        weight = (G2 + G02)/(G2 + G02p)*Q;

      }

      Kerker_weight[k2][k1][k3] = weight;

    } /* k */
  } /* #pragma omp parallel */

  /* start... */

  Min_Weight = Min_Mixing_weight;
  if (SCF_RENZOKU==-1){
    Max_Weight = Max_Mixing_weight;
    Max_Mixing_weight2 = Max_Mixing_weight;
  }
  else if (SCF_RENZOKU==1000){  /* past 3 */
    Max_Mixing_weight2 = 2.0*Max_Mixing_weight2;
    if (0.7<Max_Mixing_weight2) Max_Mixing_weight2 = 0.7;
    Max_Weight = Max_Mixing_weight2;
    SCF_RENZOKU = 0;
  }
  else{
    Max_Weight = Max_Mixing_weight2;
  }

  /****************************************************
       norm of residual charge density in k-space
  ****************************************************/

  /* get Nthrds0 */  
#pragma omp parallel shared(Nthrds0)
  {
    Nthrds0 = omp_get_num_threads();
  }
  /* allocation of array */
  My_Norm_threads = (double*)malloc(sizeof(double)*Nthrds0);
 for (Nloop=0; Nloop<Nthrds0; Nloop++) My_Norm_threads[Nloop] = 0.0;

#pragma omp parallel shared(Residual_ReRhok,Residual_ImRhok,ImRhok,ReRhok,Kerker_weight,My_Norm_threads,knum_full,My_NGrid2_Poisson,Ngrid1,Ngrid3) private(OMPID,Nthrds,Nprocs,k,spin,k2,k1,k3,Nloop,weight,tmp0,tmp1)
  {
    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    /* one-dimensionalized loop */

    for (k=OMPID*knum_full/Nthrds; k<(OMPID+1)*knum_full/Nthrds; k++){

      /* get spin, k2, k1, and k3 */

      spin = k/(My_NGrid2_Poisson*Ngrid1*Ngrid3);
      k2 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3))/(Ngrid1*Ngrid3);
      k1 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3))/Ngrid3;  
      k3 = k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3) - k1*Ngrid3;

      weight = Kerker_weight[k2][k1][k3];
      tmp0 = ReRhok[0][spin][k2][k1][k3] - ReRhok[1][spin][k2][k1][k3];
      tmp1 = ImRhok[0][spin][k2][k1][k3] - ImRhok[1][spin][k2][k1][k3];
      Residual_ReRhok[1][spin][k2][k1][k3] = tmp0;
      Residual_ImRhok[1][spin][k2][k1][k3] = tmp1;
      My_Norm_threads[OMPID] += (tmp0*tmp0 + tmp1*tmp1)*weight;

    } /* k */
  } /* #pragma omp parallel */

  /* sum of My_Norm_threads */
  My_Norm = 0.0;
  for (Nloop=0; Nloop<Nthrds0; Nloop++){
    My_Norm += My_Norm_threads[Nloop];
  }

  /* freeing of array */
  free(My_Norm_threads);

  /****************************************************
    MPI: 

    My_Norm
  ****************************************************/

  MPI_Allreduce(&My_Norm, &Norm, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

  /****************************************************
    find an optimum mixing weight
  ****************************************************/

  for (i=4; 1<=i; i--){
    NormRD[i] = NormRD[i-1];
    History_Uele[i] = History_Uele[i-1];
  }
  NormRD[0] = Norm;
  History_Uele[0] = Uele;

  if (Change_switch==1){

    if ((int)sgn(History_Uele[0]-History_Uele[1])
	  ==(int)sgn(History_Uele[1]-History_Uele[2])
       && NormRD[0]<NormRD[1]){

      /* tmp0 = 1.6*Mixing_weight; */

      tmp0 = NormRD[1]/(largest(NormRD[1]-NormRD[0],10e-10))*Mixing_weight;

      if (tmp0<Max_Weight){
        if (Min_Weight<tmp0){
          Mixing_weight = tmp0;
	}
        else{ 
          Mixing_weight = Min_Weight;
	}
      }
      else{ 
        Mixing_weight = Max_Weight;
        SCF_RENZOKU++;  
      }
    }
   
    else if ((int)sgn(History_Uele[0]-History_Uele[1])
             ==(int)sgn(History_Uele[1]-History_Uele[2])
             && NormRD[1]<NormRD[0]){

      tmp0 = NormRD[1]/(largest(NormRD[1]+NormRD[0],10e-10))*Mixing_weight;

      /* tmp0 = Mixing_weight/1.6; */

      if (tmp0<Max_Weight){
        if (Min_Weight<tmp0)
          Mixing_weight = tmp0;
        else 
          Mixing_weight = Min_Weight;
      }
      else 
        Mixing_weight = Max_Weight;

      SCF_RENZOKU = -1;  
    }

    else if ((int)sgn(History_Uele[0]-History_Uele[1])
	     !=(int)sgn(History_Uele[1]-History_Uele[2])
             && NormRD[0]<NormRD[1]){

      /* tmp0 = Mixing_weight*1.2; */

      tmp0 = NormRD[1]/(largest(NormRD[1]-NormRD[0],10e-10))*Mixing_weight;

      if (tmp0<Max_Weight){
        if (Min_Weight<tmp0)
          Mixing_weight = tmp0;
        else 
          Mixing_weight = Min_Weight;
      }
      else{ 
        Mixing_weight = Max_Weight;
        SCF_RENZOKU++;
      }
    }

    else if ((int)sgn(History_Uele[0]-History_Uele[1])
	     !=(int)sgn(History_Uele[1]-History_Uele[2])
             && NormRD[1]<NormRD[0]){

      /* tmp0 = Mixing_weight/2.0; */

      tmp0 = NormRD[1]/(largest(NormRD[1]+NormRD[0],10e-10))*Mixing_weight;

      if (tmp0<Max_Weight){
        if (Min_Weight<tmp0)
          Mixing_weight = tmp0;
        else 
          Mixing_weight = Min_Weight;
      }
      else 
        Mixing_weight = Max_Weight;

      SCF_RENZOKU = -1;
    }

    Mix_wgt = Mixing_weight;
  }

  /****************************************************
                        Mixing
  ****************************************************/

#pragma omp parallel shared(ImRhok,ReRhok,Mix_wgt,Kerker_weight,Ngrid3,Ngrid1,My_NGrid2_Poisson,knum_full) private(OMPID,Nthrds,Nprocs,k,spin,k2,k1,k3,weight,wgt0,wgt1)
  {
    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    /* one-dimensionalized loop */

    for (k=OMPID*knum_full/Nthrds; k<(OMPID+1)*knum_full/Nthrds; k++){

      /* get spin, k2, k1, and k3 */

      spin = k/(My_NGrid2_Poisson*Ngrid1*Ngrid3);
      k2 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3))/(Ngrid1*Ngrid3);
      k1 = (k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3))/Ngrid3;  
      k3 = k - spin*(My_NGrid2_Poisson*Ngrid1*Ngrid3) - k2*(Ngrid1*Ngrid3) - k1*Ngrid3;

      weight = 1.0/Kerker_weight[k2][k1][k3];
      wgt0  = Mix_wgt*weight;
      wgt1 =  1.0 - wgt0;

      /***********************************
       [0]: output at n -> input at n+1
       [1]: input  at n
      ***********************************/

      ReRhok[0][spin][k2][k1][k3] = wgt0*ReRhok[0][spin][k2][k1][k3]
 	                          + wgt1*ReRhok[1][spin][k2][k1][k3];

      ImRhok[0][spin][k2][k1][k3] = wgt0*ImRhok[0][spin][k2][k1][k3]
	                          + wgt1*ImRhok[1][spin][k2][k1][k3];

      /* correction to large changing components */

      tmp0 = ReRhok[0][spin][k2][k1][k3] - ReRhok[1][spin][k2][k1][k3];  
      tmp1 = ImRhok[0][spin][k2][k1][k3] - ImRhok[1][spin][k2][k1][k3];  

      if ( maxima_step<(fabs(tmp0)+fabs(tmp1)) ){
	ReRhok[0][spin][k2][k1][k3] = sgn(tmp0)*maxima_step + ReRhok[1][spin][k2][k1][k3]; 
	ImRhok[0][spin][k2][k1][k3] = sgn(tmp1)*maxima_step + ImRhok[1][spin][k2][k1][k3]; 
      }

    } /* k */
  } /* #pragma omp parallel */

  /****************************************************
                       shift of rho
  ****************************************************/

  for (pSCF_iter=(List_YOUSO[38]-1); 0<pSCF_iter; pSCF_iter--){
    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
	for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){
	    ReRhok[pSCF_iter][spin][k2][k1][k3] = ReRhok[pSCF_iter-1][spin][k2][k1][k3]; 
	    ImRhok[pSCF_iter][spin][k2][k1][k3] = ImRhok[pSCF_iter-1][spin][k2][k1][k3]; 
	  }
	}
      }
    }
  }

  /****************************************************
                    shift of residual rho
  ****************************************************/

  for (pSCF_iter=(List_YOUSO[38]-1); 0<pSCF_iter; pSCF_iter--){
    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
	for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){
	    Residual_ReRhok[pSCF_iter][spin][k2][k1][k3] = Residual_ReRhok[pSCF_iter-1][spin][k2][k1][k3];
	    Residual_ImRhok[pSCF_iter][spin][k2][k1][k3] = Residual_ImRhok[pSCF_iter-1][spin][k2][k1][k3]; 
	  }
	}
      }
    }
  }

  /****************************************************
        find the charge density in real space 
  ****************************************************/

  for (spin=0; spin<spinmax; spin++){

    for (k2=0; k2<My_NGrid2_Poisson; k2++){
      for (k1=0; k1<Ngrid1; k1++){
        for (k3=0; k3<Ngrid3; k3++){
          ReV2[k2][k1][k3] = ReRhok[0][spin][k2][k1][k3];
          ImV2[k2][k1][k3] = ImRhok[0][spin][k2][k1][k3];
        }
      }
    }

    Get_Value_inReal2D(ReV2,ImV2,ReV1,ImV1,Density_Grid[spin]);
  }

  if (SpinP_switch==0){
    for (MN=0; MN<My_NumGrid1; MN++){
      Density_Grid[1][MN] = Density_Grid[0][MN];
    }
  }

  /****************************************************
          store the best input charge density 
  ****************************************************/

  if (Change_switch==0 || NormRD[0]<BestNormRD ){

    BestNormRD = NormRD[0];

    for (spin=0; spin<spinmax; spin++){
      for (k2=0; k2<My_NGrid2_Poisson; k2++){
	for (k1=0; k1<Ngrid1; k1++){
	  for (k3=0; k3<Ngrid3; k3++){
            ReBestRhok[spin][k2][k1][k3] = ReRhok[2][spin][k2][k1][k3];
            ImBestRhok[spin][k2][k1][k3] = ImRhok[2][spin][k2][k1][k3];
	  }
	}
      }
    }
  }
  
  /* freeing of arrays */
  
  for (i=0; i<My_NGrid2_Poisson; i++){
    for (j=0; j<Ngrid1; j++){
      free(Kerker_weight[i][j]);
    }
    free(Kerker_weight[i]);
  }
  free(Kerker_weight);

}
