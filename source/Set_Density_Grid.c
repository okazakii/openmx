/**********************************************************************
  Set_Density_Grid.c:

     Set_Density_Grid.c is a subroutine to calculate a charge density 
     on grid by one-particle wave functions.

  Log of Set_Density_Grid.c:

     22/Nov/2001  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
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




double Set_Density_Grid(int Cnt_kind, int Calc_CntOrbital_ON, double *****CDM)
{
  int al,L0,Mul0,M0,p,size1,size2;
  int Gc_AN,Mc_AN,Mh_AN,MN;
  int n1,n2,n3,k1,k2,k3,N3[4];
  int Cwan,NO0,NO1,Rn,N,Hwan,i,j,k,n;
  int Max_Size,My_Max,top_num;
  double time0;
  int h_AN,Gh_AN,Rnh,spin,Nc,GNc,GRc,GN,Nh,Nog;
  int Nc_0,Nc_1,Nc_2,Nc_3,Nh_0,Nh_1,Nh_2,Nh_3;
  int mm;

  double tmp0,tmp1,sk1,sk2,sk3,tot_den,sum;
  double tmp0_0,tmp0_1,tmp0_2,tmp0_3;
  double sum_0,sum_1,sum_2,sum_3;

  double d1,d2,d3,cop,sip,sit,cot;
  double Re11,Re22,Re12,Im12,phi,theta;
  double x,y,z,Cxyz[4];
  double TStime,TEtime;
  double *tmp_array;
  double *tmp_array2;
  double ***Tmp_Den_Grid;
  double **Tmp_Den2_Grid;
  double *orbs0,*orbs1;
  double *orbs0_0,*orbs0_1,*orbs0_2,*orbs0_3;
  double *orbs1_0,*orbs1_1,*orbs1_2,*orbs1_3;

  double ***tmp_CDM;
  int *Snd_Size,*Rcv_Size;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  double Stime_atom, Etime_atom;

  MPI_Status stat;
  MPI_Request request;

  /* for OpenMP */
  int OMPID,Nthrds,Nprocs;

  /* snlxxx */
  double time1,time2,mflops;
  /**/
  /* MPI */
  if (atomnum<=MYID_MPI_COMM_WORLD) return 0.0;
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  /* snlxxx */
  if(myid==0 && measure_time)
    printf(" >>>>>>>Set_Density_Grid IN\n");
  /**/
  
  dtime(&TStime);

  /****************************************************
             initialize electron densities
  ****************************************************/

  for (spin=0; spin<=SpinP_switch; spin++){
    for (MN=0; MN<My_NumGrid1; MN++){
      Density_Grid[spin][MN] = 0.0;
    }
  }

  /****************************************************
                when orbital optimization
  ****************************************************/

  if (Calc_CntOrbital_ON==1 && Cnt_kind==0 && Cnt_switch==1){
      
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
       
      dtime(&Stime_atom);
      
      /* COrbs_Grid */
 
      Gc_AN = M2G[Mc_AN];
      Cwan = WhatSpecies[Gc_AN];
      NO0 = Spe_Total_CNO[Cwan]; 
      for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
        al = -1;
	for (L0=0; L0<=Spe_MaxL_Basis[Cwan]; L0++){
	  for (Mul0=0; Mul0<Spe_Num_CBasis[Cwan][L0]; Mul0++){
	    for (M0=0; M0<=2*L0; M0++){
	      al++;
	      tmp0 = 0.0;
	      for (p=0; p<Spe_Specified_Num[Cwan][al]; p++){
	        j = Spe_Trans_Orbital[Cwan][al][p];  
	        tmp0 += CntCoes[Mc_AN][al][p]*Orbs_Grid[Mc_AN][j][Nc];
	      }
	      COrbs_Grid[Mc_AN][al][Nc] = tmp0;
	    }
	  }
        }
      }

      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
    }

    /**********************************************
     MPI:

     COrbs_Grid    
     dCOrbs_Grid    
    ***********************************************/

    /* allocation of arrays  */
    Snd_Size = (int*)malloc(sizeof(int)*numprocs); 
    Rcv_Size = (int*)malloc(sizeof(int)*numprocs); 

    /* find data size for sending and receiving */

    My_Max = -10000;
    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      if (ID!=0){
        /*  sending size */
        if (F_Snd_Num[IDS]!=0){
          /* find data size */ 
          size1 = 0; 
          for (n=0; n<F_Snd_Num[IDS]; n++){
            Gc_AN = Snd_GAN[IDS][n];
            Cwan = WhatSpecies[Gc_AN];
            size1 += GridN_Atom[Gc_AN]*Spe_Total_CNO[Cwan];
          }

          Snd_Size[IDS] = size1;
          MPI_Isend(&size1, 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
        }
        else{
          Snd_Size[IDS] = 0;
        }

        /*  receiving size */
        if (F_Rcv_Num[IDR]!=0){
          MPI_Recv(&size2, 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
          Rcv_Size[IDR] = size2;
        }
        else{
          Rcv_Size[IDR] = 0;
        }
        if (F_Snd_Num[IDS]!=0) MPI_Wait(&request,&stat);
      } 
      else{
        Snd_Size[IDS] = 0;
        Rcv_Size[IDR] = 0;
      }

      if (My_Max<Snd_Size[IDS]) My_Max = Snd_Size[IDS];
      if (My_Max<Rcv_Size[IDR]) My_Max = Rcv_Size[IDR];

    }  

    MPI_Allreduce(&My_Max, &Max_Size, 1, MPI_INT, MPI_MAX, mpi_comm_level1);
    /* allocation of arrays */ 
    tmp_array  = (double*)malloc(sizeof(double)*Max_Size);
    tmp_array2 = (double*)malloc(sizeof(double)*Max_Size);

    /* send and receive COrbs_Grid */

    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      if (ID!=0){

        /* sending of data */ 

        if (F_Snd_Num[IDS]!=0){

          /* find data size */
          size1 = Snd_Size[IDS];

          /* multidimentional array to vector array */
          k = 0; 
          for (n=0; n<F_Snd_Num[IDS]; n++){
            Mc_AN = Snd_MAN[IDS][n];
            Gc_AN = Snd_GAN[IDS][n];
            Cwan = WhatSpecies[Gc_AN];
            NO0 = Spe_Total_CNO[Cwan]; 
            for (i=0; i<NO0; i++){
              for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
                tmp_array[k] = COrbs_Grid[Mc_AN][i][Nc];
                k++;
              }          
            }
          } 

          /* MPI_Isend */
          MPI_Isend(&tmp_array[0], size1, MPI_DOUBLE, IDS,
                    tag, mpi_comm_level1, &request);
        }

        /* receiving of block data */

        if (F_Rcv_Num[IDR]!=0){

          /* find data size */
          size2 = Rcv_Size[IDR]; 

          /* MPI_Recv */
          MPI_Recv(&tmp_array2[0], size2, MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

          k = 0;
          Mc_AN = F_TopMAN[IDR] - 1;
          for (n=0; n<F_Rcv_Num[IDR]; n++){
            Mc_AN++;
            Gc_AN = Rcv_GAN[IDR][n];
            Cwan = WhatSpecies[Gc_AN];
            NO0 = Spe_Total_CNO[Cwan]; 

            for (i=0; i<NO0; i++){
              for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
                COrbs_Grid[Mc_AN][i][Nc] = tmp_array2[k];
                k++;
              }          
            }
          }
        }
        if (F_Snd_Num[IDS]!=0) MPI_Wait(&request,&stat);
      } 
    }  

    /* freeing of arrays  */
    free(tmp_array);
    free(tmp_array2);
    free(Snd_Size);
    free(Rcv_Size);
  }

  /**********************************************
              calculate Tmp_Den_Grid
  ***********************************************/
    
  dtime(&time1);

  /* allocation of Tmp_Den_Grid */

  Tmp_Den_Grid = (double***)malloc(sizeof(double**)*(SpinP_switch+1)); 
  for (i=0; i<(SpinP_switch+1); i++){
    Tmp_Den_Grid[i] = (double**)malloc(sizeof(double*)*(Matomnum+MatomnumF+1)); 
    Tmp_Den_Grid[i][0] = (double*)malloc(sizeof(double)*1); 
    for (Mc_AN=1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
      Gc_AN = F_M2G[Mc_AN];
      Tmp_Den_Grid[i][Mc_AN] = (double*)malloc(sizeof(double)*GridN_Atom[Gc_AN]);
    }
  }

#pragma omp parallel shared(List_YOUSO,time_per_atom,Tmp_Den_Grid,Orbs_Grid,COrbs_Grid,Cnt_switch,Cnt_kind,GListTAtoms2,GListTAtoms1,NumOLG,CDM,SpinP_switch,WhatSpecies,ncn,F_G2M,natn,Spe_Total_CNO,M2G) private(OMPID,Nthrds,Nprocs,Mc_AN,h_AN,Stime_atom,Etime_atom,Gc_AN,Cwan,NO0,Gh_AN,Mh_AN,Rnh,Hwan,NO1,spin,i,j,tmp_CDM,Nog,Nc_0,Nc_1,Nc_2,Nc_3,Nh_0,Nh_1,Nh_2,Nh_3,orbs0_0,orbs0_1,orbs0_2,orbs0_3,orbs1_0,orbs1_1,orbs1_2,orbs1_3,sum_0,sum_1,sum_2,sum_3,tmp0_0,tmp0_1,tmp0_2,tmp0_3,mm,Nc,Nh,orbs0,orbs1,sum,tmp0)
  {

    orbs0 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    orbs1 = (double*)malloc(sizeof(double)*List_YOUSO[7]);

    orbs0_0 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    orbs0_1 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    orbs0_2 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    orbs0_3 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    orbs1_0 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    orbs1_1 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    orbs1_2 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    orbs1_3 = (double*)malloc(sizeof(double)*List_YOUSO[7]);

    tmp_CDM = (double***)malloc(sizeof(double**)*(SpinP_switch+1)); 
    for (i=0; i<(SpinP_switch+1); i++){
      tmp_CDM[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
      for (j=0; j<List_YOUSO[7]; j++){
	tmp_CDM[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
      }
    }

    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    for (Mc_AN=(OMPID*Matomnum/Nthrds+1); Mc_AN<((OMPID+1)*Matomnum/Nthrds+1); Mc_AN++){

      dtime(&Stime_atom);

      /* set data on Mc_AN */

      Gc_AN = M2G[Mc_AN];
      Cwan = WhatSpecies[Gc_AN];
      NO0 = Spe_Total_CNO[Cwan]; 

      for (spin=0; spin<=SpinP_switch; spin++){
	for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
	  Tmp_Den_Grid[spin][Mc_AN][Nc] = 0.0;
	}
      }

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	/* set data on h_AN */
    
	Gh_AN = natn[Gc_AN][h_AN];
	Mh_AN = F_G2M[Gh_AN];
	Rnh = ncn[Gc_AN][h_AN];
	Hwan = WhatSpecies[Gh_AN];
	NO1 = Spe_Total_CNO[Hwan];

	/* store CDM into tmp_CDM */

	for (spin=0; spin<=SpinP_switch; spin++){
	  for (i=0; i<NO0; i++){
	    for (j=0; j<NO1; j++){
	      tmp_CDM[spin][i][j] = CDM[spin][Mc_AN][h_AN][i][j];
	    }
	  }
	}

	/* summation of non-zero elements */

	for (Nog=0; Nog<NumOLG[Mc_AN][h_AN]-3; Nog+=4){

	  Nc_0 = GListTAtoms1[Mc_AN][h_AN][Nog];
	  Nc_1 = GListTAtoms1[Mc_AN][h_AN][Nog+1];
	  Nc_2 = GListTAtoms1[Mc_AN][h_AN][Nog+2];
	  Nc_3 = GListTAtoms1[Mc_AN][h_AN][Nog+3];

	  Nh_0 = GListTAtoms2[Mc_AN][h_AN][Nog];
	  Nh_1 = GListTAtoms2[Mc_AN][h_AN][Nog+1];
	  Nh_2 = GListTAtoms2[Mc_AN][h_AN][Nog+2];
	  Nh_3 = GListTAtoms2[Mc_AN][h_AN][Nog+3];

	  /* Now under the orbital optimization */
	  if (Cnt_kind==0 && Cnt_switch==1){
	    for (i=0; i<NO0; i++){
	      orbs0_0[i] = COrbs_Grid[Mc_AN][i][Nc_0];
	      orbs0_1[i] = COrbs_Grid[Mc_AN][i][Nc_1];
	      orbs0_2[i] = COrbs_Grid[Mc_AN][i][Nc_2];
	      orbs0_3[i] = COrbs_Grid[Mc_AN][i][Nc_3];
	    }
	    for (j=0; j<NO1; j++){
	      orbs1_0[j] = COrbs_Grid[Mh_AN][j][Nh_0];
	      orbs1_1[j] = COrbs_Grid[Mh_AN][j][Nh_1];
	      orbs1_2[j] = COrbs_Grid[Mh_AN][j][Nh_2];
	      orbs1_3[j] = COrbs_Grid[Mh_AN][j][Nh_3];
	    }
	  }
	  /* else if ! "now under the orbital optimization" */
	  else{
	    for (i=0; i<NO0; i++){
	      orbs0_0[i] = Orbs_Grid[Mc_AN][i][Nc_0];
	      orbs0_1[i] = Orbs_Grid[Mc_AN][i][Nc_1];
	      orbs0_2[i] = Orbs_Grid[Mc_AN][i][Nc_2];
	      orbs0_3[i] = Orbs_Grid[Mc_AN][i][Nc_3];
	    }
	    for (j=0; j<NO1; j++){
	      orbs1_0[j] = Orbs_Grid[Mh_AN][j][Nh_0];
	      orbs1_1[j] = Orbs_Grid[Mh_AN][j][Nh_1];
	      orbs1_2[j] = Orbs_Grid[Mh_AN][j][Nh_2];
	      orbs1_3[j] = Orbs_Grid[Mh_AN][j][Nh_3];
	    }
	  }

	  for (spin=0; spin<=SpinP_switch; spin++){

	    /* Tmp_Den_Grid */

	    sum_0 = 0.0;
	    sum_1 = 0.0;
	    sum_2 = 0.0;
	    sum_3 = 0.0;

	    for (i=0; i<NO0; i++){

	      tmp0_0 = 0.0;
	      tmp0_1 = 0.0;
	      tmp0_2 = 0.0;
	      tmp0_3 = 0.0;

	      for (j=0; j<NO1; j++){
		tmp0_0 += orbs1_0[j]*tmp_CDM[spin][i][j];
		tmp0_1 += orbs1_1[j]*tmp_CDM[spin][i][j];
		tmp0_2 += orbs1_2[j]*tmp_CDM[spin][i][j];
		tmp0_3 += orbs1_3[j]*tmp_CDM[spin][i][j];
	      }

	      sum_0 += orbs0_0[i]*tmp0_0;
	      sum_1 += orbs0_1[i]*tmp0_1;
	      sum_2 += orbs0_2[i]*tmp0_2;
	      sum_3 += orbs0_3[i]*tmp0_3;
	    }

	    Tmp_Den_Grid[spin][Mc_AN][Nc_0] += sum_0;
	    Tmp_Den_Grid[spin][Mc_AN][Nc_1] += sum_1;
	    Tmp_Den_Grid[spin][Mc_AN][Nc_2] += sum_2;
	    Tmp_Den_Grid[spin][Mc_AN][Nc_3] += sum_3;

	  }
	}

	mm = NumOLG[Mc_AN][h_AN]-(NumOLG[Mc_AN][h_AN]/4)*4;

	if (mm!=0){

	  for (Nog=NumOLG[Mc_AN][h_AN]-mm; Nog<NumOLG[Mc_AN][h_AN]; Nog++){
 
	    Nc = GListTAtoms1[Mc_AN][h_AN][Nog];
	    Nh = GListTAtoms2[Mc_AN][h_AN][Nog]; 
 
	    /* Now under the orbital optimization */
	    if (Cnt_kind==0 && Cnt_switch==1){
	      for (i=0; i<NO0; i++){
		orbs0[i] = COrbs_Grid[Mc_AN][i][Nc];
	      }
	      for (j=0; j<NO1; j++){
		orbs1[j] = COrbs_Grid[Mh_AN][j][Nh];
	      }
	    }
	    /* else if ! "now under the orbital optimization" */
	    else{
	      for (i=0; i<NO0; i++){
		orbs0[i] = Orbs_Grid[Mc_AN][i][Nc];
	      }
	      for (j=0; j<NO1; j++){
		orbs1[j] = Orbs_Grid[Mh_AN][j][Nh];
	      }
	    }

	    for (spin=0; spin<=SpinP_switch; spin++){
 
	      /* Tmp_Den_Grid */
 
	      sum = 0.0;
	      for (i=0; i<NO0; i++){
		tmp0 = 0.0;
		for (j=0; j<NO1; j++){
		  tmp0 += orbs1[j]*tmp_CDM[spin][i][j];
		}
		sum += orbs0[i]*tmp0;
	      }
 
	      Tmp_Den_Grid[spin][Mc_AN][Nc] += sum;
	    }
	  }
	}

      } /* h_AN */

      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;

    } /* Mc_AN */

    /* freeing of arrays */ 

    free(orbs0);
    free(orbs1);

    free(orbs0_0);
    free(orbs0_1);
    free(orbs0_2);
    free(orbs0_3);
    free(orbs1_0);
    free(orbs1_1);
    free(orbs1_2);
    free(orbs1_3);

    for (i=0; i<(SpinP_switch+1); i++){
      for (j=0; j<List_YOUSO[7]; j++){
	free(tmp_CDM[i][j]);
      }
      free(tmp_CDM[i]);
    }
    free(tmp_CDM);

#pragma omp flush(Tmp_Den_Grid)

  } /* #pragma omp parallel */

  dtime(&time2);
  if(myid==0 && measure_time){
    printf("Time for Part1=%18.5f\n",(time2-time1));fflush(stdout);
  }

  /******************************************************
   MPI:
        Tmp_Den_Grid 
  ******************************************************/
 
  /* allocation of arrays  */
  Snd_Size = (int*)malloc(sizeof(int)*numprocs); 
  Rcv_Size = (int*)malloc(sizeof(int)*numprocs); 
   
  /* find data size for sending and recieving */

  My_Max = -10000;
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){
      /*  sending size */
      if (F_Snd_Num[IDS]!=0){
        /* find data size */ 
        size1 = 0; 
        for (n=0; n<F_Snd_Num[IDS]; n++){
          Gc_AN = Snd_GAN[IDS][n];
          size1 += GridN_Atom[Gc_AN]*(SpinP_switch+1);
        }

        Snd_Size[IDS] = size1;
        MPI_Isend(&size1, 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
      }
      else{
        Snd_Size[IDS] = 0;
      }

      /*  receiving size */
      if (F_Rcv_Num[IDR]!=0){
        MPI_Recv(&size2, 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
        Rcv_Size[IDR] = size2;
      }
      else{
        Rcv_Size[IDR] = 0;
      }
      if (F_Snd_Num[IDS]!=0) MPI_Wait(&request,&stat);
    } 
    else{
      Snd_Size[IDS] = 0;
      Rcv_Size[IDR] = 0;
    }

    if (My_Max<Snd_Size[IDS]) My_Max = Snd_Size[IDS];
    if (My_Max<Rcv_Size[IDR]) My_Max = Rcv_Size[IDR];
  }  

  MPI_Allreduce(&My_Max, &Max_Size, 1, MPI_INT, MPI_MAX, mpi_comm_level1);
  tmp_array  = (double*)malloc(sizeof(double)*Max_Size);
  tmp_array2 = (double*)malloc(sizeof(double)*Max_Size);

  /* send and receive Tmp_Den_Grid */

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){

      /* sending of data  */

      if (F_Snd_Num[IDS]!=0){

        /* find data size */
        size1 = Snd_Size[IDS];

        /* multidimentional array to vector array */
        k = 0; 
	for (spin=0; spin<=SpinP_switch; spin++){
          for (n=0; n<F_Snd_Num[IDS]; n++){
            Mc_AN = Snd_MAN[IDS][n];
            Gc_AN = Snd_GAN[IDS][n];
            for (i=0; i<GridN_Atom[Gc_AN]; i++){
              tmp_array[k] = Tmp_Den_Grid[spin][Mc_AN][i];
              k++;
            }
  	  } 
	}

        /* MPI_Isend */
        MPI_Isend(&tmp_array[0], size1, MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);
      }

      /*****************************
         receiving of block data
      *****************************/

      if (F_Rcv_Num[IDR]!=0){
  
        /* find data size */
        size2 = Rcv_Size[IDR]; 

        /* MPI_Recv */
        MPI_Recv(&tmp_array2[0], size2, MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

        k = 0;
        for (spin=0; spin<=SpinP_switch; spin++){
          Mc_AN = F_TopMAN[IDR] - 1;
          for (n=0; n<F_Rcv_Num[IDR]; n++){
            Mc_AN++;
            Gc_AN = Rcv_GAN[IDR][n];

            for (i=0; i<GridN_Atom[Gc_AN]; i++){
              Tmp_Den_Grid[spin][Mc_AN][i] = tmp_array2[k];
              k++;
            }
          }
	}
      }
      if (F_Snd_Num[IDS]!=0) MPI_Wait(&request,&stat);
    } 
  }  

  /* freeing of arrays  */
  free(tmp_array);
  free(tmp_array2);
  free(Snd_Size);
  free(Rcv_Size);

  /**********************************************
                 calc Density_Grid
  ***********************************************/

  /* snlxxx */
  dtime(&time1);
  /**/

  for (Mc_AN=1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){

    dtime(&Stime_atom);

    Gc_AN = F_M2G[Mc_AN];
    Cwan = WhatSpecies[Gc_AN];

    for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
      MN = MGridListAtom[Mc_AN][Nc];
      GRc = CellListAtom[Mc_AN][Nc];

      if ( (0<=MN && Solver!=4) || (0<=MN && Solver==4 && atv_ijk[GRc][1]==0) ){ 

	for (spin=0; spin<=SpinP_switch; spin++){
	  Density_Grid[spin][MN] += Tmp_Den_Grid[spin][Mc_AN][Nc];
	}
      }
    }

    dtime(&Etime_atom);
    time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
  }

  /* snlxxx */
  dtime(&time2);
  if(myid==0 && measure_time)
    printf("Time for Part2=%18.5f\n",time2-time1);fflush(stdout);
  /**/

  /******************************************************
         add Tmp_Den2_Grid in terms of FNAN2 
  ******************************************************/

  /* allocation of array */
  Tmp_Den2_Grid = (double**)malloc(sizeof(double*)*(SpinP_switch+1));
  for (i=0; i<(SpinP_switch+1); i++){
    Tmp_Den2_Grid[i] = (double*)malloc(sizeof(double)*FNAN2_Grid);
  }

  for (spin=0; spin<=SpinP_switch; spin++){

    /* MPI */

    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      if (ID!=0){

        /*****************************
              sending of data 
        *****************************/

        if (Num_Snd_FNAN2_Grid[IDS]!=0){
        
          tmp_array = (double*)malloc(sizeof(double)*Num_Snd_FNAN2_Grid[IDS]);

          /* vector array */
          for (i=0; i<Num_Snd_FNAN2_Grid[IDS]; i++){
            Gc_AN = Snd_FNAN2_At[IDS][i];
            Mc_AN = F_G2M[Gc_AN];
            Nc    = Snd_FNAN2_Nc[IDS][i];
            tmp_array[i] = Tmp_Den_Grid[spin][Mc_AN][Nc];
          }

          /* MPI_Isend */
          MPI_Isend(&tmp_array[0], Num_Snd_FNAN2_Grid[IDS], MPI_DOUBLE,
                    IDS, tag, mpi_comm_level1, &request);

        }

        /*****************************
              receiving of data
        *****************************/

        if (Num_Rcv_FNAN2_Grid[IDR]!=0){
          top_num = TopMAN2_Grid[IDR];
          /* MPI_Recv */
          MPI_Recv(&Tmp_Den2_Grid[spin][top_num], Num_Rcv_FNAN2_Grid[IDR],
                   MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);
	}

        if (Num_Snd_FNAN2_Grid[IDS]!=0){
          MPI_Wait(&request,&stat);
          free(tmp_array);
        }
      }
    }
  } /* spin */

  for (spin=0; spin<=SpinP_switch; spin++){

    /* Density_Grid += Tmp_Den2_Grid */ 

    if (Solver!=4){
      for (i=0; i<FNAN2_Grid; i++){
	MN = Rcv_FNAN2_MN[i];
	Density_Grid[spin][MN] += Tmp_Den2_Grid[spin][i];
      }
    }
    else if (Solver==4){
      for (i=0; i<FNAN2_Grid; i++){
	MN = Rcv_FNAN2_MN[i];
        GRc = Rcv_FNAN2_GRc[i];

        if (atv_ijk[GRc][1]==0){
  	  Density_Grid[spin][MN] += Tmp_Den2_Grid[spin][i];
	}
      }
    }
  }

  /* if (SpinP_switch==0) */

  if (SpinP_switch==0){
    for (MN=0; MN<My_NumGrid1; MN++){
      Density_Grid[1][MN] = Density_Grid[0][MN]; 
    }
  }

  /****************************************************
   Conjugate complex of Density_Grid[3][MN] due to
   difference in the definition between density matrix
   and charge density
  ****************************************************/

  if (SpinP_switch==3){
    for (MN=0; MN<My_NumGrid1; MN++){
      Density_Grid[3][MN] = -Density_Grid[3][MN]; 
    }
  }


  /*
  {

    int nn1,nn0,MN1,MN2;

  n2 = Ngrid2/2;
  n3 = Ngrid3/2;

  for (n1=Start_Grid1[myid]; n1<=End_Grid1[myid]; n1++){
    nn1 = My_Cell0[n1];
    nn0 = n1 - Start_Grid1[myid]; 
    MN1 = nn1*Ngrid2*Ngrid3;
    MN2 = n2*Ngrid3;
    MN = MN1 + MN2 + n3;
    printf("%2d %18.15f %18.15f\n",n1,Density_Grid[0][MN],Density_Grid[1][MN]);
  }

  }
  */




  /* freeing of arrays */

  for (i=0; i<(SpinP_switch+1); i++){
    free(Tmp_Den2_Grid[i]);
  }
  free(Tmp_Den2_Grid);

  for (i=0; i<(SpinP_switch+1); i++){
    for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
      free(Tmp_Den_Grid[i][Mc_AN]);
    }
    free(Tmp_Den_Grid[i]);
  }
  free(Tmp_Den_Grid);

  /* elapsed time */
  dtime(&TEtime);
  time0 = TEtime - TStime;
  /* snlxxx */
  if(myid==0 && measure_time)
    printf(">>>>>>>Set_Density_Grid RETURN time=%18.5f\n",time0);
  /**/

  return time0;
}




void diagonalize_nc_density()
{
  int MN,Mc_AN,Gc_AN,Nog,GNc,GRc;
  double Re11,Re22,Re12,Im12;
  double phi[2],theta[2],sit,cot,sip,cop;
  double d1,d2,d3,x,y,z,Cxyz[4];
  double Nup[2],Ndown[2];
  /* for OpenMP */
  int OMPID,Nthrds,Nprocs;

  if (atomnum<=MYID_MPI_COMM_WORLD) return;

#pragma omp parallel shared(Density_Grid,My_NumGrid1) private(OMPID,Nthrds,Nprocs,MN,Re11,Re22,Re12,Im12,Nup,Ndown,theta,phi)
  {

    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    for (MN=OMPID*My_NumGrid1/Nthrds; MN<(OMPID+1)*My_NumGrid1/Nthrds; MN++){

      Re11 = Density_Grid[0][MN];
      Re22 = Density_Grid[1][MN];
      Re12 = Density_Grid[2][MN];
      Im12 = Density_Grid[3][MN];

      EulerAngle_Spin( 1, Re11, Re22, Re12, Im12, Re12, -Im12, Nup, Ndown, theta, phi );

      Density_Grid[0][MN] = Nup[0];
      Density_Grid[1][MN] = Ndown[0];
      Density_Grid[2][MN] = theta[0];
      Density_Grid[3][MN] = phi[0];
    }

#pragma omp flush(Density_Grid)

  } /* #pragma omp parallel */

}
