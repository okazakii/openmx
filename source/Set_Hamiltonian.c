/**********************************************************************
  Set_Hamiltonian.c:

     Set_Hamiltonian.c is a subroutine to make Hamiltonian matrix
     within LDA or GGA.

  Log of Set_Hamiltonian.c:

     24/April/2002  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "openmx_common.h"

#define  measure_time   0

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


void simple_mixing_H(double *****Hcnt);


double Set_Hamiltonian(char *mode,
                       int SCF_iter,
                       int SucceedReadingDMfile,
                       int Cnt_kind,
                       double *****H0,
                       double *****HNL,
                       double *****CDM,
		       double *****H)
{
  /***************************************************************
      Cnt_kind
        0:  Uncontracted Hamiltonian    
        1:  Contracted Hamiltonian    
  ***************************************************************/

  double time0;
  int Mc_AN,Gc_AN,Mh_AN,h_AN,h_ANs,h_ANe,Gh_AN;
  int i,j,Cwan,Hwan,NO0,NO1,spinmax;
  int Rnh,Rnk,spin,N,NumC[4];
  int n1,n2,n3,L0,Mul0,M0,L1,Mul1,M1;
  int Nc,Ncs,GNc,GRc,Nog,Nh,MN,XC_P_switch;
  int Nc_0,Nc_1,Nc_2,Nc_3,Nh_0,Nh_1,Nh_2,Nh_3;
  int mm;

  double x,y,z,dx,dy,dz;
  double bc,dv,r,theta,phi,sum,tmp0,tmp1;
  double tmp0_0,tmp0_1,tmp0_2,tmp0_3;

  double xo,yo,zo,S_coordinate[3];
  double ***ChiV0;
  double **tmp_ChiV0_0,**tmp_ChiV0_1,**tmp_ChiV0_2,**tmp_ChiV0_3;
  double *tmp_Orbs_Grid_0,*tmp_Orbs_Grid_1,*tmp_Orbs_Grid_2,*tmp_Orbs_Grid_3;
  double *Chi0;

  double ***tmp_H;
  double TStime,TEtime;
  double TStime0,TEtime0;
  double TStime1,TEtime1;
  double TStime2,TEtime2;
  double TStime3,TEtime3;
  int numprocs,myid,tag=999,ID;
  double Stime_atom, Etime_atom;
  /* snlxxx */
  double time1,time2,mflops;

  int GNh_0,GNh_1,GNh_2,GNh_3;
  int GRh_0,GRh_1,GRh_2,GRh_3;
  double Cxyz_0[4],Cxyz_1[4],Cxyz_2[4],Cxyz_3[4];
  double x_0,y_0,z_0;
  double x_1,y_1,z_1;
  double x_2,y_2,z_2;
  double x_3,y_3,z_3;
  double dx_0,dy_0,dz_0;
  double dx_1,dy_1,dz_1;
  double dx_2,dy_2,dz_2;
  double dx_3,dy_3,dz_3;

  /* for OpenMP */
  int OMPID,Nthrds,Nprocs;
  int Nloop,OneD_Nloop,*OneD2Mc_AN,*OneD2h_AN;

  if (atomnum<=MYID_MPI_COMM_WORLD) return 0.0;
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if(myid==0 && measure_time){ 
    printf(" *******Set_Hamiltonian IN\n");
  }

  MPI_Barrier(mpi_comm_level1);
  dtime(&TStime);

  if (myid==Host_ID && strcasecmp(mode,"stdout")==0 && 0<level_stdout ){
    printf("<Set_Hamiltonian>  Hamiltonian matrix for VNA+dVH+Vxc...\n");fflush(stdout);
  }

  /*****************************************************
                  adding H0+HNL to H 
  *****************************************************/

  /* spin non-collinear */

  if (SpinP_switch==3){
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];    
      Cwan = WhatSpecies[Gc_AN];
      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
        Gh_AN = natn[Gc_AN][h_AN];
        Hwan = WhatSpecies[Gh_AN];
        for (i=0; i<Spe_Total_NO[Cwan]; i++){
          for (j=0; j<Spe_Total_NO[Hwan]; j++){

            if (ProExpn_VNA==0){
              H[0][Mc_AN][h_AN][i][j] = F_Kin_flag*H0[0][Mc_AN][h_AN][i][j]
		+ F_NL_flag*HNL[0][Mc_AN][h_AN][i][j];
              H[1][Mc_AN][h_AN][i][j] = F_Kin_flag*H0[0][Mc_AN][h_AN][i][j]
		+ F_NL_flag*HNL[1][Mc_AN][h_AN][i][j];
              H[2][Mc_AN][h_AN][i][j] = F_NL_flag*HNL[2][Mc_AN][h_AN][i][j];
              H[3][Mc_AN][h_AN][i][j] = 0.0;
	    }
            else{
              H[0][Mc_AN][h_AN][i][j] = F_Kin_flag*H0[0][Mc_AN][h_AN][i][j]
		+ F_VNA_flag*HVNA[Mc_AN][h_AN][i][j]
		+ F_NL_flag*HNL[0][Mc_AN][h_AN][i][j];
              H[1][Mc_AN][h_AN][i][j] = F_Kin_flag*H0[0][Mc_AN][h_AN][i][j]
		+ F_VNA_flag*HVNA[Mc_AN][h_AN][i][j]
		+ F_NL_flag*HNL[1][Mc_AN][h_AN][i][j];
              H[2][Mc_AN][h_AN][i][j] = F_NL_flag*HNL[2][Mc_AN][h_AN][i][j];
              H[3][Mc_AN][h_AN][i][j] = 0.0;
            }

            /* Effective Hubbard Hamiltonain --- added by MJ */

	    if( (Hub_U_switch==1 || Constraint_NCS_switch==1) && F_U_flag==1 && 2<=SCF_iter ){
	      H[0][Mc_AN][h_AN][i][j] += H_Hub[0][Mc_AN][h_AN][i][j];
	      H[1][Mc_AN][h_AN][i][j] += H_Hub[1][Mc_AN][h_AN][i][j];
	      H[2][Mc_AN][h_AN][i][j] += H_Hub[2][Mc_AN][h_AN][i][j];
	    }

          }
        }
      }
    }
  }

  /* spin collinear */

  else{
    /* snlxxx */
    dtime(&time1);
    /**/
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];    
      Cwan = WhatSpecies[Gc_AN];
      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
        Gh_AN = natn[Gc_AN][h_AN];
        Hwan = WhatSpecies[Gh_AN];
        for (i=0; i<Spe_Total_NO[Cwan]; i++){
          for (j=0; j<Spe_Total_NO[Hwan]; j++){
            for (spin=0; spin<=SpinP_switch; spin++){

              if (ProExpn_VNA==0){
                H[spin][Mc_AN][h_AN][i][j] = F_Kin_flag*H0[0][Mc_AN][h_AN][i][j]
		                           + F_NL_flag*HNL[spin][Mc_AN][h_AN][i][j];
	      }
              else{

                H[spin][Mc_AN][h_AN][i][j] = F_Kin_flag*H0[0][Mc_AN][h_AN][i][j]
		                           + F_VNA_flag*HVNA[Mc_AN][h_AN][i][j]
		                           + F_NL_flag*HNL[spin][Mc_AN][h_AN][i][j];
              }

	      /* Effective Hubbard Hamiltonain --- added by MJ */
	      if( (Hub_U_switch==1 || Constraint_NCS_switch==1) && F_U_flag==1 && 2<=SCF_iter ){
		H[spin][Mc_AN][h_AN][i][j] += H_Hub[spin][Mc_AN][h_AN][i][j];
	      }
            }
          }
        }
      }
    }
    /* snlxxx */
    dtime(&time2);
    if(myid==0 && measure_time) 
      printf("Time for Part1=%18.5f\n",time2-time1);fflush(stdout);
    /**/
  }

  if (Cnt_kind==1) {
    Contract_Hamiltonian( H, CntH, OLP, CntOLP );
    if (SO_switch==1) Contract_iHNL(iHNL,iCntHNL);
  }

  /*****************************************************
   calculation of matrix elements for dVH + Vxc (+ VNA)
  *****************************************************/

  XC_P_switch = 1;

  /* snlxxx */
  dtime(&time1);
  /**/

  Set_Vpot(SCF_iter,XC_P_switch,CDM);

  /* snlxxx */
  dtime(&time2);
  if(myid==0 && measure_time) 
    printf("Set_Vpot  time=%18.5f\n",time2-time1);fflush(stdout);
  /**/

  dtime(&time1);

  /* one-dimensionalize the Mc_AN and h_AN loops */

  OneD_Nloop = 0;
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];    
    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
      OneD_Nloop++;
    }
  }  

  OneD2Mc_AN = (int*)malloc(sizeof(int)*(OneD_Nloop+1));
  OneD2h_AN = (int*)malloc(sizeof(int)*(OneD_Nloop+1));

  OneD_Nloop = 0;
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];    
    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
      OneD2Mc_AN[OneD_Nloop] = Mc_AN; 
      OneD2h_AN[OneD_Nloop] = h_AN; 
      OneD_Nloop++;
    }
  }

  /* OpenMP */

#pragma omp parallel shared(gtv,Gxyz,atv,GridListAtom,CellListAtom,OneD2h_AN,OneD2Mc_AN,OneD_Nloop,GridN_Atom,H,CntH,List_YOUSO,GridVol,Orbs_Grid,Vpot_Grid,MGridListAtom,GListTAtoms2,GListTAtoms1,Matomnum,NumOLG,SpinP_switch,Spe_Total_CNO,Spe_Total_NO,Cnt_kind,WhatSpecies,ncn,natn,F_G2M,FNAN)  private(OMPID,Nthrds,Nprocs,h_AN,Gh_AN,Mh_AN,Rnh,Hwan,NO1,spin,i,j,tmp_H,Nog,Nc_0,Nc_1,Nc_2,Nc_3,Nh_0,Nh_1,Nh_2,Nh_3,tmp_ChiV0_0,tmp_ChiV0_1,tmp_ChiV0_2,tmp_ChiV0_3,tmp_Orbs_Grid_0,tmp_Orbs_Grid_1,tmp_Orbs_Grid_2,tmp_Orbs_Grid_3,mm,Nc,Nh,tmp0,tmp0_0,tmp0_1,tmp0_2,tmp0_3,Mc_AN,Gc_AN,Cwan,NO0,Ncs,ChiV0,Chi0,Nloop)
  {
    /* allocation of arrays */ 

    Chi0 = (double*)malloc(sizeof(double)*List_YOUSO[7]);

    ChiV0 = (double***)malloc(sizeof(double**)*4);
    for (i=0; i<4; i++){
      ChiV0[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
      for (j=0; j<List_YOUSO[7]; j++){
	ChiV0[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[11]);
      }
    }

    tmp_H = (double***)malloc(sizeof(double**)*4);
    for (i=0; i<4; i++){
      tmp_H[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
      for (j=0; j<List_YOUSO[7]; j++){
	tmp_H[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
      }
    }

    tmp_ChiV0_0 = (double**)malloc(sizeof(double*)*4);
    for (i=0; i<4; i++){
      tmp_ChiV0_0[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    }

    tmp_ChiV0_1 = (double**)malloc(sizeof(double*)*4);
    for (i=0; i<4; i++){
      tmp_ChiV0_1[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    }
    
    tmp_ChiV0_2 = (double**)malloc(sizeof(double*)*4);
    for (i=0; i<4; i++){
      tmp_ChiV0_2[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    }
    
    tmp_ChiV0_3 = (double**)malloc(sizeof(double*)*4);
    for (i=0; i<4; i++){
      tmp_ChiV0_3[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    }
    
    tmp_Orbs_Grid_0 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    tmp_Orbs_Grid_1 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    tmp_Orbs_Grid_2 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    tmp_Orbs_Grid_3 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    
    /* get info. on OpenMP */ 
     
    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    /* one-dimensionalized loop */

    for (Nloop=OMPID*OneD_Nloop/Nthrds; Nloop<(OMPID+1)*OneD_Nloop/Nthrds; Nloop++){

      /* get Mc_AN and h_AN */

      Mc_AN = OneD2Mc_AN[Nloop];
      h_AN  = OneD2h_AN[Nloop];

      /* set data on Mc_AN */

      Gc_AN = M2G[Mc_AN];    
      Cwan = WhatSpecies[Gc_AN];

      if (Cnt_kind==0) NO0 = Spe_Total_NO[Cwan];
      else             NO0 = Spe_Total_CNO[Cwan];
    
      for (spin=0; spin<=SpinP_switch; spin++){
	for (i=0; i<NO0; i++){

	  for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
	    ChiV0[spin][i][Nc] = Orbs_Grid[Mc_AN][i][Nc]*Vpot_Grid[spin][MGridListAtom[Mc_AN][Nc]];
	  }

	  Ncs = GridN_Atom[Gc_AN] - GridN_Atom[Gc_AN]%4; 
	  for (Nc=Ncs; Nc<GridN_Atom[Gc_AN]; Nc++){
	    ChiV0[spin][i][Nc] = Orbs_Grid[Mc_AN][i][Nc]*Vpot_Grid[spin][MGridListAtom[Mc_AN][Nc]];
	  }
	}
      }

      /* set data on h_AN */

      Gh_AN = natn[Gc_AN][h_AN];
      Mh_AN = F_G2M[Gh_AN];

      Rnh = ncn[Gc_AN][h_AN];
      Hwan = WhatSpecies[Gh_AN];

      if (Cnt_kind==0)
	NO1 = Spe_Total_NO[Hwan];
      else   
	NO1 = Spe_Total_CNO[Hwan];

      /* initialize tmp_H */

      for (spin=0; spin<=SpinP_switch; spin++){
	for (i=0; i<NO0; i++){
	  for (j=0; j<NO1; j++){
	    tmp_H[spin][i][j] = 0.0;                          
	  }
	}
      }

      /* summation of non-zero elements */

      for (Nog=0; Nog<NumOLG[Mc_AN][h_AN]-3; Nog+=4){

        Nc_0 = GListTAtoms1[Mc_AN][h_AN][Nog+0];
        Nc_1 = GListTAtoms1[Mc_AN][h_AN][Nog+1];
        Nc_2 = GListTAtoms1[Mc_AN][h_AN][Nog+2];
        Nc_3 = GListTAtoms1[Mc_AN][h_AN][Nog+3];

        Nh_0 = GListTAtoms2[Mc_AN][h_AN][Nog+0];
        Nh_1 = GListTAtoms2[Mc_AN][h_AN][Nog+1];
        Nh_2 = GListTAtoms2[Mc_AN][h_AN][Nog+2];
        Nh_3 = GListTAtoms2[Mc_AN][h_AN][Nog+3];

	/* store ChiV0 in tmp_ChiV0 */
 
	for (spin=0; spin<=SpinP_switch; spin++){
	  for (i=0; i<NO0; i++){
	    tmp_ChiV0_0[spin][i] = ChiV0[spin][i][Nc_0];
	    tmp_ChiV0_1[spin][i] = ChiV0[spin][i][Nc_1];
	    tmp_ChiV0_2[spin][i] = ChiV0[spin][i][Nc_2];
	    tmp_ChiV0_3[spin][i] = ChiV0[spin][i][Nc_3];
	  }
	}

	/* store Orbs_Grid in tmp_Orbs_Grid */

	for (j=0; j<NO1; j++){
	  tmp_Orbs_Grid_0[j] = Orbs_Grid[Mh_AN][j][Nh_0];
	  tmp_Orbs_Grid_1[j] = Orbs_Grid[Mh_AN][j][Nh_1];
	  tmp_Orbs_Grid_2[j] = Orbs_Grid[Mh_AN][j][Nh_2];
	  tmp_Orbs_Grid_3[j] = Orbs_Grid[Mh_AN][j][Nh_3];
	}

	/* integration */
 
	for (spin=0; spin<=SpinP_switch; spin++){
	  for (i=0; i<NO0; i++){

	    tmp0_0 = tmp_ChiV0_0[spin][i];
	    tmp0_1 = tmp_ChiV0_1[spin][i];
	    tmp0_2 = tmp_ChiV0_2[spin][i];
	    tmp0_3 = tmp_ChiV0_3[spin][i];

	    for (j=0; j<NO1; j++){
	      tmp_H[spin][i][j] += (  tmp0_0*tmp_Orbs_Grid_0[j]
			             +tmp0_1*tmp_Orbs_Grid_1[j]
				     +tmp0_2*tmp_Orbs_Grid_2[j]
				     +tmp0_3*tmp_Orbs_Grid_3[j]);
	    }

	  }
	}
      }

      mm = NumOLG[Mc_AN][h_AN]-(NumOLG[Mc_AN][h_AN]/4)*4;

      if(mm != 0){
 
	for (Nog=NumOLG[Mc_AN][h_AN]-mm; Nog<NumOLG[Mc_AN][h_AN]; Nog++){
 
	  Nc = GListTAtoms1[Mc_AN][h_AN][Nog];
	  Nh = GListTAtoms2[Mc_AN][h_AN][Nog];

	  /* store ChiV0 in tmp_ChiV0_0 */
 
	  for (spin=0; spin<=SpinP_switch; spin++){
	    for (i=0; i<NO0; i++){
	      tmp_ChiV0_0[spin][i] = ChiV0[spin][i][Nc];
	    }
	  }
 
	  /* store Orbs_Grid in tmp_Orbs_Grid */

	  for (j=0; j<NO1; j++){
	    tmp_Orbs_Grid_0[j] = Orbs_Grid[Mh_AN][j][Nh];
	  }

	  /* integration */

	  for (spin=0; spin<=SpinP_switch; spin++){

	    for (i=0; i<NO0; i++){
	      tmp0 = tmp_ChiV0_0[spin][i];
	      for (j=0; j<NO1; j++){
		tmp_H[spin][i][j] += tmp0*tmp_Orbs_Grid_0[j];
	      }
	    }
	  }
	}
      }

      /* add tmp_H to H or CntH */

      if (Cnt_kind==0){

	for (spin=0; spin<=SpinP_switch; spin++){
	  for (i=0; i<NO0; i++){
	    for (j=0; j<NO1; j++){
	      H[spin][Mc_AN][h_AN][i][j] += tmp_H[spin][i][j]*GridVol;
	    }
	  }
	}
      }
 
      else{
	for (spin=0; spin<=SpinP_switch; spin++){
	  for (i=0; i<NO0; i++){
	    for (j=0; j<NO1; j++){
	      CntH[spin][Mc_AN][h_AN][i][j] += tmp_H[spin][i][j]*GridVol;
	    }
	  }
	}
      }

    } /* Nloop */

    /* freeing of arrays */ 

    free(tmp_Orbs_Grid_3);
    free(tmp_Orbs_Grid_2);
    free(tmp_Orbs_Grid_1);
    free(tmp_Orbs_Grid_0);

    for (i=0; i<4; i++){
      free(tmp_ChiV0_3[i]);
    }
    free(tmp_ChiV0_3);

    for (i=0; i<4; i++){
      free(tmp_ChiV0_2[i]);
    }
    free(tmp_ChiV0_2);

    for (i=0; i<4; i++){
      free(tmp_ChiV0_1[i]);
    }
    free(tmp_ChiV0_1);

    for (i=0; i<4; i++){
      free(tmp_ChiV0_0[i]);
    }
    free(tmp_ChiV0_0);

    for (i=0; i<4; i++){
      for (j=0; j<List_YOUSO[7]; j++){
	free(tmp_H[i][j]);
      }
      free(tmp_H[i]);
    }
    free(tmp_H);

    for (i=0; i<4; i++){
      for (j=0; j<List_YOUSO[7]; j++){
	free(ChiV0[i][j]);
      }
      free(ChiV0[i]);
    }
    free(ChiV0);

    free(Chi0);

#pragma omp flush(H,CntH)

  } /* #pragma omp parallel */


  /*
 {
   double sum0,sum1;

   sum0 = 0.0;
   sum1 = 0.0;

   spin = 0;
   for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
     Gc_AN = M2G[Mc_AN];    
     Cwan = WhatSpecies[Gc_AN];

     for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

       Gh_AN = natn[Gc_AN][h_AN];
       Hwan = WhatSpecies[Gh_AN];

       printf("Gc_AN=%2d Gh_AN=%2d\n",Gc_AN,Gh_AN);

       for (i=0; i<Spe_Total_NO[Cwan]; i++){
	 for (j=0; j<Spe_Total_NO[Hwan]; j++){
	   printf("%10.5f ",H[spin][Mc_AN][h_AN][i][j]);

	   sum0 += fabs(H[spin][Mc_AN][h_AN][i][j]);
	   sum1 += H[spin][Mc_AN][h_AN][i][j];

	 }
	 printf("\n");
       }
     }
   }

   printf("myid=%2d sum0=%18.15f\n",myid,sum0);
   printf("myid=%2d sum1=%18.15f\n",myid,sum1);
 }
  */



  dtime(&time2);
  if(myid==0 && measure_time){
    printf("Time for Part2=%18.5f\n",(time2-time1));fflush(stdout);
  }
  
  /*****************************************************
   When the restart file is used, a simple mixing of the
   Hamiltonian matrix is made at SCF_iter==2 
  *****************************************************/

  /*
    if (SucceedReadingDMfile==1 && SCF_iter==2){
    if (Cnt_kind==0)  simple_mixing_H(H);
    else              simple_mixing_H(CntH);
    }
  */

  /* freeing of arrays */ 

  free(OneD2h_AN);
  free(OneD2Mc_AN);

  /* for time */

  MPI_Barrier(mpi_comm_level1);
  dtime(&TEtime);
  time0 = TEtime - TStime;
  /* snlxxx */
  if(myid==0 && measure_time) 
    printf(" *******Set_Hamiltonian RETURN time=%18.5f\n",time0);
  /**/

  return time0;
}





void simple_mixing_H(double *****Hcnt)
{
  int Mc_AN,Gc_AN,h_AN,i,j,can;
  int Gh_AN,pSCF,spin,Rn;
  int wan1,wan2,TNO1,TNO2;
  int h_AN0,Gh_AN0,Rn0,wan20,TNO20;
  int i_vec[20],*p_vec,po;
  int my_check,exit_flag;
  int pFNAN;
  int numprocs,myid;
  char fileHKS[YOUSO10];
  FILE *fp;
  double *tmpvec;
  double Uele0;
  char buf[fp_bsize];          /* setvbuf */

  /* allocation of array */

  tmpvec = (double*)malloc(sizeof(double)*List_YOUSO[7]);

  /****************************************************
   List_YOUSO[23] spin poralized     1
                  non spin poralized 0
   List_YOUSO[1]  atomnum
   List_YOUSO[8]  max # of atoms in a rcut-off cluster
   List_YOUSO[7]  max # of orbitals in an atom
  ****************************************************/

  my_check = 1;

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];
    wan1 = WhatSpecies[Gc_AN];
    TNO1 = Spe_Total_CNO[wan1];

    sprintf(fileHKS,"%s%s_rst/%s.rst%i",filepath,filename,filename,Gc_AN);

    if ((fp = fopen(fileHKS,"r")) != NULL){

#ifdef xt3
      setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      /****************************************************
       List_YOUSO[23] 0:  non spin poralized
                      1:  spin poralized
                      3:  spin non-collinear
       List_YOUSO[1]  atomnum
       List_YOUSO[8]  max # of atoms in a rcut-off cluster
       List_YOUSO[7]  max # of orbitals in an atom
      ****************************************************/

      fread(i_vec,sizeof(int),10,fp);

      pFNAN = i_vec[8];   

      if ( i_vec[0]!=SpinP_switch   ||
           i_vec[1]!=List_YOUSO[23] || 
           i_vec[2]!=List_YOUSO[1]  ||
           i_vec[4]!=List_YOUSO[7]  ||
           i_vec[5]!=atomnum        ||
           i_vec[6]!=wan1           ||
           i_vec[7]!=TNO1           
         )
      {
        printf("Failed (1) in reading the restart file %s\n",fileHKS); fflush(stdout);     
        my_check = 0;
      }

      if (my_check!=0){

        /****************************************************
                   read Gh_AN0, Rn0, wan20, TNO20
        ****************************************************/

        p_vec = (int*)malloc(sizeof(int)*(pFNAN+1)*4);
        fread(p_vec, sizeof(int), (pFNAN+1)*4, fp);
        fread(&Uele0,sizeof(double),1,fp);

        /****************************************************
          store Hamiltonian to appropriate position while 
          comparing Gh_AN, Rn, wan2, TNO2
        ****************************************************/

        for (spin=0; spin<=SpinP_switch; spin++){ 
	  for (h_AN0=0; h_AN0<=pFNAN; h_AN0++){
	    Gh_AN0 = p_vec[              h_AN0];
	    Rn0    = p_vec[(pFNAN+1)*1 + h_AN0];
	    wan20  = p_vec[(pFNAN+1)*2 + h_AN0];
	    TNO20  = p_vec[(pFNAN+1)*3 + h_AN0]; 
 
	    exit_flag = 0;
	    h_AN = 0;

	    do {

	      Gh_AN = natn[Gc_AN][h_AN];
	      Rn    = ncn[Gc_AN][h_AN];
	      wan2  = WhatSpecies[Gh_AN];
	      TNO2  = Spe_Total_CNO[wan2];

	      if ( Gh_AN==Gh_AN0 &&
		   Rn==Rn0       && 
		   wan2==wan20   &&
		   TNO2==TNO20 )
		{

		  for (i=0; i<TNO1; i++){
		    fread(&tmpvec[0],sizeof(double),TNO2,fp);

  		    for (j=0; j<TNO2; j++){
                      Hcnt[spin][Mc_AN][h_AN][i][j] = 0.9999*tmpvec[j] + 0.0001*Hcnt[spin][Mc_AN][h_AN][i][j];
		    }
		  }

		  exit_flag = 1;
		}

	      h_AN++;

	    } while (h_AN<=FNAN[Gc_AN] && exit_flag==0);

            /* In case appropriate one is not found, just read the Hamiltonian */ 

            if (exit_flag==0){
	      for (i=0; i<TNO1; i++){
	        fread(&tmpvec[0],sizeof(double),TNO20,fp);
	      }
            }
	  }
	}

        /* freeing of array */ 
        free(p_vec);

      }

      /* close the file */
      fclose(fp);
    }
    else{
      printf("Failed (2) in reading the restart file %s\n",fileHKS); fflush(stdout);     
      my_check = 0;
    }
  }

  /* freeing of array */
  free(tmpvec);
} 

