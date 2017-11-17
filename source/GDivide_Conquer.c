/**********************************************************************
  GDivide_Conquer.c:

     GDivide_Conquer.c is a subroutine to perform a generalized
     divide and conquer method for the eigenvalue problem

  Log of GDivide_Conquer.c:

     20/Dec/2004  Released by T.Ozaki

***********************************************************************/

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

static double GDC_Col(
              int SCF_iter,
              double *****Hks, double ****OLP0,
              double *****CDM,
              double *****EDM,
              double Eele0[2], double Eele1[2]);

static double GDC_NonCol(
              int SCF_iter,
              double *****Hks,
              double *****ImNL,
              double ****OLP0,
              double *****CDM,
              double *****EDM,
              double Eele0[2], double Eele1[2]);

double GDivide_Conquer(int SCF_iter,
                       double *****Hks,
                       double *****ImNL,
                       double ****OLP0,
                       double *****CDM,
                       double *****EDM,
                       double Eele0[2], double Eele1[2])
{
  double time0;

  /****************************************************
         collinear without spin-orbit coupling
  ****************************************************/

  if ( (SpinP_switch==0 || SpinP_switch==1) && SO_switch==0 ){
    time0 = GDC_Col(SCF_iter, Hks, OLP0, CDM, EDM, Eele0, Eele1);
  }

  /****************************************************
         collinear with spin-orbit coupling
  ****************************************************/

  else if ( (SpinP_switch==0 || SpinP_switch==1) && SO_switch==1 ){
    printf("Spin-orbit coupling is not supported for collinear DFT calculations.\n");
    MPI_Finalize();
    exit(1);
  }

  /****************************************************
   non-collinear with and without spin-orbit coupling
  ****************************************************/

  else if (SpinP_switch==3){
    time0 = GDC_NonCol(SCF_iter, Hks, ImNL, OLP0, CDM, EDM, Eele0, Eele1);
  }

  return time0;
}



static double GDC_Col(
              int SCF_iter,
              double *****Hks, double ****OLP0,
              double *****CDM,
              double *****EDM,
              double Eele0[2], double Eele1[2])
{
  static int firsttime=1;
  int Mc_AN,Mc_AN_GDC,Gc_AN,i,Gi,wan,wanA,wanB,Anum;
  int size1,size2,num,NUM,NUM1,n2,Cwan,Hwan;
  int ih,ig,ian,j,kl,jg,jan,Bnum,m,n,spin;
  int l,i1,j1,P_min,m_size;
  int Rni,Rnj,Gi_AN,Gj_AN;
  int h_AN_GDC,Mh_AN_GDC,Gh_AN_GDC;
  int po,po0,loopN,tno1,tno2,h_AN,Gh_AN;
  double My_TZ,TZ,sum,FermiF,time0,tmp1,tmp2;
  double My_Num_State,Num_State,x,Dnum;
  double TStime,TEtime;
  double My_Eele0[2],My_Eele1[2];
  double max_x=50.0;
  double ChemP_MAX,ChemP_MIN,spin_degeneracy;
  double **S_DC,***H_DC,*ko,*M1;
  double **C;
  double ***EVal;
  double ******Residues;
  double ***PDOS_DC;
  int **MP,*Msize;
  double *tmp_array;
  double *tmp_array2;
  int *Snd_H_Size,*Rcv_H_Size;
  int *Snd_S_Size,*Rcv_S_Size;
  int numprocs,myid,ID,IDS,IDR,tag=999;
  double Stime_atom, Etime_atom;
  double OLP_eigen_cut=Threshold_OLP_Eigen;

  MPI_Status stat;
  MPI_Request request;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  dtime(&TStime);

  /****************************************************
    allocation of arrays:

    int Msize[Matomnum+1];
    double EVal[SpinP_switch+1][Matomnum+1][n2];
  ****************************************************/

  Snd_H_Size = (int*)malloc(sizeof(int)*numprocs);
  Rcv_H_Size = (int*)malloc(sizeof(int)*numprocs);
  Snd_S_Size = (int*)malloc(sizeof(int)*numprocs);
  Rcv_S_Size = (int*)malloc(sizeof(int)*numprocs);

  Msize = (int*)malloc(sizeof(int)*(Matomnum+1));

  /* find Msize */

  Msize[0] = 1;
  for (Mc_AN_GDC=1; Mc_AN_GDC<=Matomnum_GDC; Mc_AN_GDC++){

    Mc_AN = Mnatn_GDC[Mc_AN_GDC][0];
    Gc_AN = M2G[Mc_AN];

    Anum = 1;
    for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
      Gi = natn[Gc_AN][i];
      wanA = WhatSpecies[Gi];
      Anum = Anum + Spe_Total_CNO[wanA];
    }
    NUM = Anum - 1;

    for (h_AN_GDC=0; h_AN_GDC<=MFNAN_GDC[Mc_AN_GDC]; h_AN_GDC++){
      Mh_AN_GDC = Mnatn_GDC[Mc_AN_GDC][h_AN_GDC];
      Msize[Mh_AN_GDC] = NUM;
    }
  }

  /* allocate EVal */

  m_size = 0;

  EVal = (double***)malloc(sizeof(double**)*(SpinP_switch+1));
  for (spin=0; spin<=SpinP_switch; spin++){
    EVal[spin] = (double**)malloc(sizeof(double*)*(Matomnum+1));

    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];

      if (Mc_AN==0){
        FNAN[0] = 1;
        SNAN[0] = 0;
        n2 = 1;
      }
      else{
        n2 = Msize[Mc_AN] + 3;
      }

      m_size += n2;
      EVal[spin][Mc_AN] = (double*)malloc(sizeof(double)*n2);
    }
  }

  if (firsttime)
  PrintMemory("GDivide_Conquer: EVal",sizeof(double)*m_size,NULL);

  if (2<=level_stdout){
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
        printf("<DC> myid=%i Mc_AN=%2d Gc_AN=%2d Msize=%3d\n",
        myid,Mc_AN,M2G[Mc_AN],Msize[Mc_AN]);
    }
  }

  /****************************************************
    allocation of arrays:

    double Residues[SpinP_switch+1]
                   [Matomnum+1]
                   [FNAN[Gc_AN]+1]
                   [Spe_Total_CNO[Gc_AN]] 
                   [Spe_Total_CNO[Gh_AN]] 
                   [NUM2]
     To reduce the memory size, the size of NUM2 is
     needed to be found in the loop.  
  ****************************************************/

  m_size = 0;
  Residues = (double******)malloc(sizeof(double*****)*(SpinP_switch+1));
  for (spin=0; spin<=SpinP_switch; spin++){
    Residues[spin] = (double*****)malloc(sizeof(double****)*(Matomnum+1));
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];

      if (Mc_AN==0){
        Gc_AN = 0;
        FNAN[0] = 1;
        tno1 = 1;
        n2 = 1;
      }
      else{
        wanA = WhatSpecies[Gc_AN];
        tno1 = Spe_Total_CNO[wanA];
        n2 = Msize[Mc_AN] + 2;
      }

      Residues[spin][Mc_AN] = (double****)malloc(sizeof(double***)*(FNAN[Gc_AN]+1));

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

        if (Mc_AN==0){
          tno2 = 1;
        }
        else {
          Gh_AN = natn[Gc_AN][h_AN];
          wanB = WhatSpecies[Gh_AN];
          tno2 = Spe_Total_CNO[wanB];
        }

        Residues[spin][Mc_AN][h_AN] = (double***)malloc(sizeof(double**)*tno1);
        for (i=0; i<tno1; i++){
          Residues[spin][Mc_AN][h_AN][i] = (double**)malloc(sizeof(double*)*tno2);
          for (j=0; j<tno2; j++){
            Residues[spin][Mc_AN][h_AN][i][j] = (double*)malloc(sizeof(double)*n2);
	  }
        }

        m_size += tno1*tno2*n2;
      }
    }
  }

  if (firsttime)
  PrintMemory("GDivide_Conquer: Residues",sizeof(double)*m_size,NULL);

  /****************************************************
    allocation of arrays:

    double PDOS_DOS[SpinP_switch+1]
                   [Matomnum+1]
                   [NUM]
  ****************************************************/

  m_size = 0;
  PDOS_DC = (double***)malloc(sizeof(double**)*(SpinP_switch+1));
  for (spin=0; spin<=SpinP_switch; spin++){
    PDOS_DC[spin] = (double**)malloc(sizeof(double*)*(Matomnum+1));
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];

      if (Mc_AN==0) n2 = 1;
      else          n2 = Msize[Mc_AN] + 2;

      m_size += n2;
      PDOS_DC[spin][Mc_AN] = (double*)malloc(sizeof(double)*n2);
    }
  }

  if (firsttime)
  PrintMemory("GDivide_Conquer: PDOS_DC",sizeof(double)*m_size,NULL);

  /****************************************************
   MPI

   Hks
  ****************************************************/

  /***********************************
             set data size
  ************************************/

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){
      tag = 999;

      /* find data size to send block data */
      if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0){

        size1 = 0;
        for (spin=0; spin<=SpinP_switch; spin++){
          for (n=0; n<(F_Snd_Num[IDS]+S_Snd_Num[IDS]); n++){
            Mc_AN = Snd_MAN[IDS][n];
            Gc_AN = Snd_GAN[IDS][n];
            Cwan = WhatSpecies[Gc_AN]; 
            tno1 = Spe_Total_CNO[Cwan];
            for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
              Gh_AN = natn[Gc_AN][h_AN];        
              Hwan = WhatSpecies[Gh_AN];
              tno2 = Spe_Total_CNO[Hwan];
              for (i=0; i<tno1; i++){
                for (j=0; j<tno2; j++){
                  size1++; 
                } 
              } 
	    }
          }
	}
 
        Snd_H_Size[IDS] = size1;
        MPI_Isend(&size1, 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
      }
      else{
        Snd_H_Size[IDS] = 0;
      }

      /* receiving of size of data */

      if ((F_Rcv_Num[IDR]+S_Rcv_Num[IDR])!=0){
        MPI_Recv(&size2, 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
        Rcv_H_Size[IDR] = size2;
      }
      else{
        Rcv_H_Size[IDR] = 0;
      }

      if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0) MPI_Wait(&request,&stat);

    }
    else{
      Snd_H_Size[IDS] = 0;
      Rcv_H_Size[IDR] = 0;
    }
  }

  /***********************************
             data transfer
  ************************************/

  tag = 999;
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){

      /*****************************
              sending of data 
      *****************************/

      if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0){

        size1 = Snd_H_Size[IDS];

        /* allocation of array */

        tmp_array = (double*)malloc(sizeof(double)*size1);

        /* multidimentional array to vector array */

        num = 0;
        for (spin=0; spin<=SpinP_switch; spin++){
          for (n=0; n<(F_Snd_Num[IDS]+S_Snd_Num[IDS]); n++){
            Mc_AN = Snd_MAN[IDS][n];
            Gc_AN = Snd_GAN[IDS][n];
            Cwan = WhatSpecies[Gc_AN]; 
            tno1 = Spe_Total_CNO[Cwan];
            for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
              Gh_AN = natn[Gc_AN][h_AN];
              Hwan = WhatSpecies[Gh_AN];
              tno2 = Spe_Total_CNO[Hwan];
              for (i=0; i<tno1; i++){
                for (j=0; j<tno2; j++){
                  tmp_array[num] = Hks[spin][Mc_AN][h_AN][i][j];
                  num++;
                } 
              } 
	    }
          }
	}

        MPI_Isend(&tmp_array[0], size1, MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);

      }

      /*****************************
         receiving of block data
      *****************************/

      if ((F_Rcv_Num[IDR]+S_Rcv_Num[IDR])!=0){

        size2 = Rcv_H_Size[IDR];
        
        /* allocation of array */
        tmp_array2 = (double*)malloc(sizeof(double)*size2);
        
        MPI_Recv(&tmp_array2[0], size2, MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

        num = 0;
        for (spin=0; spin<=SpinP_switch; spin++){
          Mc_AN = S_TopMAN[IDR] - 1;  /* S_TopMAN should be used. */
          for (n=0; n<(F_Rcv_Num[IDR]+S_Rcv_Num[IDR]); n++){
            Mc_AN++;
            Gc_AN = Rcv_GAN[IDR][n];
            Cwan = WhatSpecies[Gc_AN]; 
            tno1 = Spe_Total_CNO[Cwan];

            for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
              Gh_AN = natn[Gc_AN][h_AN];        
              Hwan = WhatSpecies[Gh_AN];
              tno2 = Spe_Total_CNO[Hwan];
              for (i=0; i<tno1; i++){
                for (j=0; j<tno2; j++){
                  Hks[spin][Mc_AN][h_AN][i][j] = tmp_array2[num];
                  num++;
		}
	      }
	    }
	  }        
	}

        /* freeing of array */
        free(tmp_array2);
      }

      if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0){
        MPI_Wait(&request,&stat);
        free(tmp_array); /* freeing of array */
      }
    }
  }

  /****************************************************
   MPI

   OLP0
  ****************************************************/

  /***********************************
             set data size
  ************************************/

  if (SCF_iter<=2){

    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      if (ID!=0){
        tag = 999;

        /* find data size to send block data */
        if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0){
          size1 = 0;
          for (n=0; n<(F_Snd_Num[IDS]+S_Snd_Num[IDS]); n++){
            Mc_AN = Snd_MAN[IDS][n];
            Gc_AN = Snd_GAN[IDS][n];
            Cwan = WhatSpecies[Gc_AN]; 
            tno1 = Spe_Total_CNO[Cwan];
            for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
              Gh_AN = natn[Gc_AN][h_AN];        
              Hwan = WhatSpecies[Gh_AN];
              tno2 = Spe_Total_CNO[Hwan];
              for (i=0; i<tno1; i++){
                for (j=0; j<tno2; j++){
                  size1++; 
                } 
              } 
	    }
          }

          Snd_S_Size[IDS] = size1;
          MPI_Isend(&size1, 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
        }
        else{
          Snd_S_Size[IDS] = 0;
        }

        /* receiving of size of data */
 
        if ((F_Rcv_Num[IDR]+S_Rcv_Num[IDR])!=0){
          MPI_Recv(&size2, 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
          Rcv_S_Size[IDR] = size2;
        }
        else{
          Rcv_S_Size[IDR] = 0;
        }

        if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0) MPI_Wait(&request,&stat);
      }
      else{
        Snd_S_Size[IDS] = 0;
        Rcv_S_Size[IDR] = 0;
      }
    }

    /***********************************
               data transfer
    ************************************/

    tag = 999;
    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      if (ID!=0){

        /*****************************
                sending of data 
        *****************************/

        if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0){

          size1 = Snd_S_Size[IDS];

          /* allocation of array */

          tmp_array = (double*)malloc(sizeof(double)*size1);

          /* multidimentional array to vector array */

          num = 0;

          for (n=0; n<(F_Snd_Num[IDS]+S_Snd_Num[IDS]); n++){
            Mc_AN = Snd_MAN[IDS][n];
            Gc_AN = Snd_GAN[IDS][n];
            Cwan = WhatSpecies[Gc_AN]; 
            tno1 = Spe_Total_CNO[Cwan];
            for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
              Gh_AN = natn[Gc_AN][h_AN];        
              Hwan = WhatSpecies[Gh_AN];
              tno2 = Spe_Total_CNO[Hwan];
              for (i=0; i<tno1; i++){
                for (j=0; j<tno2; j++){
                  tmp_array[num] = OLP0[Mc_AN][h_AN][i][j];
                  num++;
                } 
              } 
	    }
          }

          MPI_Isend(&tmp_array[0], size1, MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);
        }

        /*****************************
           receiving of block data
        *****************************/

        if ((F_Rcv_Num[IDR]+S_Rcv_Num[IDR])!=0){
          
          size2 = Rcv_S_Size[IDR];
        
         /* allocation of array */
          tmp_array2 = (double*)malloc(sizeof(double)*size2);
         
          MPI_Recv(&tmp_array2[0], size2, MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

          num = 0;
          Mc_AN = S_TopMAN[IDR] - 1; /* S_TopMAN should be used. */
          for (n=0; n<(F_Rcv_Num[IDR]+S_Rcv_Num[IDR]); n++){
            Mc_AN++;
            Gc_AN = Rcv_GAN[IDR][n];
            Cwan = WhatSpecies[Gc_AN]; 
            tno1 = Spe_Total_CNO[Cwan];

            for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
              Gh_AN = natn[Gc_AN][h_AN];        
              Hwan = WhatSpecies[Gh_AN];
              tno2 = Spe_Total_CNO[Hwan];
              for (i=0; i<tno1; i++){
                for (j=0; j<tno2; j++){
                  OLP0[Mc_AN][h_AN][i][j] = tmp_array2[num];
                  num++;
   	        }
	      }
	    }
	  }        

          /* freeing of array */
          free(tmp_array2);

        }

        if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0){
          MPI_Wait(&request,&stat);
          free(tmp_array); /* freeing of array */
	}
      }
    }
  }

  /****************************************************
            find the total number of electrons 
  ****************************************************/

  My_TZ = 0.0;
  for (i=1; i<=Matomnum; i++){
    Gc_AN = M2G[i];
    wan = WhatSpecies[Gc_AN];
    My_TZ = My_TZ + Spe_Core_Charge[wan];
  }

  /* MPI, My_TZ */

  MPI_Barrier(mpi_comm_level1);
  MPI_Allreduce(&My_TZ, &TZ, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

  /****************************************************
      Setting of Hamiltonian and overlap matrices
  ****************************************************/

  for (Mc_AN_GDC=1; Mc_AN_GDC<=Matomnum_GDC; Mc_AN_GDC++){

    dtime(&Stime_atom);

    Mc_AN = Mnatn_GDC[Mc_AN_GDC][0];
    Gc_AN = M2G[Mc_AN];
    wan = WhatSpecies[Gc_AN];

    /***********************************************
      find the size of matrix for the atom Mc_AN
                and set the MP vector

     Note:
         MP indicates the starting position of
              atom i in arraies H and S

     allocation: 
       int MP[(MFNAN_GDC[Mc_AN_GDC]+1)]
                    [List_YOUSO[2]];
    ***********************************************/

    MP = (int**)malloc(sizeof(int*)*(MFNAN_GDC[Mc_AN_GDC]+1));
    for (h_AN_GDC=0; h_AN_GDC<=MFNAN_GDC[Mc_AN_GDC]; h_AN_GDC++){
      MP[h_AN_GDC] = (int*)malloc(sizeof(int)*List_YOUSO[2]);
    }
    
    /* set MP[0][i] */
    Anum = 1;
    for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
      MP[0][i] = Anum;
      Gi = natn[Gc_AN][i];
      wanA = WhatSpecies[Gi];
      Anum = Anum + Spe_Total_CNO[wanA];
    }
    NUM = Anum - 1;
    n2 = NUM + 3;
    
    /* set MP[h_AN_GDC][j] */
    for (h_AN_GDC=1; h_AN_GDC<=MFNAN_GDC[Mc_AN_GDC]; h_AN_GDC++){

      Mh_AN_GDC = Mnatn_GDC[Mc_AN_GDC][h_AN_GDC];
      Gh_AN_GDC = M2G[Mh_AN_GDC];

      for (j=0; j<=FNAN[Gh_AN_GDC]; j++){

        Gj_AN = natn[Gh_AN_GDC][j];
        Rnj   =  ncn[Gh_AN_GDC][j];

        po0= 0;
        i = 0;
        do{

          Gi_AN = natn[Gc_AN][i];
          Rni   =  ncn[Gc_AN][i];

          if (Gj_AN==Gi_AN && Rnj==Rni){
            MP[h_AN_GDC][j] = MP[0][i];
            po0 = 1;
          }

          i++;
	} while ( i<=(FNAN[Gc_AN]+SNAN[Gc_AN]) && po0==0 );

        if (po0==0){
          printf("Error (1) in GDivide_Conquer.c  Mc_AN_GDC=%2d h_AN_GDC=%2d j=%2d\n",Mc_AN_GDC,h_AN_GDC,j); 
          MPI_Finalize();
          exit(0);
        }        
      }
    }

    /***********************************************
     allocation of arrays:
     
     double S_DC[n2][n2];     
     double H_DC[n2][n2];     
     double ko[n2];
    ***********************************************/

    S_DC = (double**)malloc(sizeof(double*)*n2);
    for (i=0; i<n2; i++){
      S_DC[i] = (double*)malloc(sizeof(double)*n2);
    }

    H_DC = (double***)malloc(sizeof(double**)*(SpinP_switch+1));
    for (spin=0; spin<=SpinP_switch; spin++){
      H_DC[spin] = (double**)malloc(sizeof(double*)*n2);
      for (i=0; i<n2; i++){
        H_DC[spin][i] = (double*)malloc(sizeof(double)*n2);
      }
    }

    ko = (double*)malloc(sizeof(double)*n2);
    M1 = (double*)malloc(sizeof(double)*n2);

    C = (double**)malloc(sizeof(double*)*n2);
    for (i=0; i<n2; i++){
      C[i] = (double*)malloc(sizeof(double)*n2);
    }

    /***********************************************
     construct cluster full matrices of Hamiltonian
           and overlap for the atom Mc_AN             
    ***********************************************/

    for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
      ig = natn[Gc_AN][i];
      ian = Spe_Total_CNO[WhatSpecies[ig]];
      Anum = MP[0][i];
      ih = S_G2M[ig];

      for (j=0; j<=(FNAN[Gc_AN]+SNAN[Gc_AN]); j++){

	kl = RMI1[Mc_AN][i][j];
	jg = natn[Gc_AN][j];
	jan = Spe_Total_CNO[WhatSpecies[jg]];
	Bnum = MP[0][j];

	if (0<=kl){

          if (SCF_iter<=2){
 	    for (m=0; m<ian; m++){
	      for (n=0; n<jan; n++){
	        S_DC[Anum+m][Bnum+n] = OLP0[ih][kl][m][n];
	      }
	    }
	  }

	  for (spin=0; spin<=SpinP_switch; spin++){
	    for (m=0; m<ian; m++){
	      for (n=0; n<jan; n++){
		H_DC[spin][Anum+m][Bnum+n] = Hks[spin][ih][kl][m][n];
	      }
	    }
	  }
	}

	else{

          if (SCF_iter<=2){
	    for (m=0; m<ian; m++){
	      for (n=0; n<jan; n++){
	        S_DC[Anum+m][Bnum+n] = 0.0;
	      }
	    }
	  }

	  for (spin=0; spin<=SpinP_switch; spin++){
	    for (m=0; m<ian; m++){
	      for (n=0; n<jan; n++){
		H_DC[spin][Anum+m][Bnum+n] = 0.0;
	      }
	    }
	  }
	}
      }
    }

    /****************************************************
     Solve the generalized eigenvalue problem
     HC = SCE

     1) diagonalize S
     2) search negative eigenvalues of S  
    ****************************************************/

    if (SCF_iter<=2){

      Eigen_lapack(S_DC,ko,NUM,NUM);

      /***********************************************
             Searching of negative eigenvalues
      ************************************************/

      P_min = 1;
      for (l=1; l<=NUM; l++){
        if (ko[l]<OLP_eigen_cut){
          P_min = l + 1;
          if (3<=level_stdout){
            printf("<DC>  Negative EV of OLP %2d %15.12f\n",l,ko[l]);
	  }
        }
      }

      S12[Mc_AN_GDC][0][0] = P_min;

      for (l=1; l<P_min; l++)     M1[l] = 0.0;
      for (l=P_min; l<=NUM; l++)  M1[l] = 1.0/sqrt(ko[l]);

      for (i1=1; i1<=NUM; i1++){
	for (j1=1; j1<=NUM; j1++){
	  S_DC[i1][j1]           = S_DC[i1][j1]*M1[j1];
          S12[Mc_AN_GDC][i1][j1] = S_DC[i1][j1];
	}
      }

    }

    else{

      P_min = (int)S12[Mc_AN_GDC][0][0];

      for (i1=1; i1<=NUM; i1++){
	for (j1=1; j1<=NUM; j1++){
	  S_DC[i1][j1] = S12[Mc_AN_GDC][i1][j1];
	}
      }
    }

    /***********************************************
      transform Hamiltonian matrix
    ************************************************/

    for (spin=0; spin<=SpinP_switch; spin++){

      /* transpose S */
      for (i1=1; i1<=NUM; i1++){
	for (j1=i1+1; j1<=NUM; j1++){
	  tmp1 = S_DC[i1][j1];
	  tmp2 = S_DC[j1][i1];
	  S_DC[i1][j1] = tmp2;
	  S_DC[j1][i1] = tmp1;
	}
      }

      /* H * U * M1 */

      for (j1=1; j1<=NUM; j1++){
        for (i1=1; i1<=NUM; i1++){
          sum = 0.0;
          for (l=1; l<=NUM; l++){
	    sum += H_DC[spin][i1][l]*S_DC[j1][l];
          }
          C[j1][i1] = sum;
        }
      }

      /* M1 * U^+ H * U * M1 */

      for (j1=1; j1<=NUM; j1++){
        for (i1=1; i1<=NUM; i1++){
          sum = 0.0;
          for (l=1; l<=NUM; l++){
	    sum += S_DC[i1][l]*C[j1][l];
          }
          H_DC[spin][j1][i1] = sum;
        }
      }

      /* H_DC to C (transposition) */

      for (i1=P_min; i1<=NUM; i1++){
        for (j1=P_min; j1<=NUM; j1++){
          C[j1-(P_min-1)][i1-(P_min-1)] = H_DC[spin][i1][j1];
        }
      }

      /***********************************************
       diagonalize the trasformed Hamiltonian matrix
      ************************************************/

      NUM1 = NUM - (P_min - 1);
      Eigen_lapack(C,ko,NUM1,NUM1);

      /* C to H (transposition) */

      for (i1=1; i1<=NUM; i1++){
	for (j1=1; j1<=NUM1; j1++){
	  H_DC[spin][j1][i1] = C[i1][j1];
	}
      }

      /***********************************************
        transformation to the original eigenvectors
                        NOTE 244P
      ***********************************************/

      /* transpose */

      for (i1=1; i1<=NUM; i1++){
	for (j1=i1+1; j1<=NUM; j1++){
	  tmp1 = S_DC[i1][j1];
	  tmp2 = S_DC[j1][i1];
	  S_DC[i1][j1] = tmp2;
	  S_DC[j1][i1] = tmp1;
	}
      }

      for (j1=1; j1<=NUM1; j1++){
	for (l=NUM; P_min<=l; l--){
	  H_DC[spin][j1][l] = H_DC[spin][j1][l-(P_min-1)];
	}
      }

      for (j1=1; j1<=NUM1; j1++){
        for (i1=1; i1<=NUM; i1++){
          sum = 0.0;
          for (l=P_min; l<=NUM; l++){
            sum += S_DC[i1][l]*H_DC[spin][j1][l];
          }
          C[i1][j1] = sum;
        }
      }

      /***********************************************
          store eigenvalues and residues of poles
      ***********************************************/

      for (h_AN_GDC=0; h_AN_GDC<=MFNAN_GDC[Mc_AN_GDC]; h_AN_GDC++){

        Mh_AN_GDC = Mnatn_GDC[Mc_AN_GDC][h_AN_GDC];
        Gh_AN_GDC = M2G[Mh_AN_GDC];
        wanA = WhatSpecies[Gh_AN_GDC];
        tno1 = Spe_Total_CNO[wanA];
        Anum = MP[h_AN_GDC][0];

        for (i1=1; i1<=NUM; i1++){
          EVal[spin][Mh_AN_GDC][i1-1] = 1000.0;
        }

        for (i1=1; i1<=NUM1; i1++){
          EVal[spin][Mh_AN_GDC][i1-1] = ko[i1];
        }

	for (i=0; i<tno1; i++){
	  for (h_AN=0; h_AN<=FNAN[Gh_AN_GDC]; h_AN++){

	    Gh_AN = natn[Gh_AN_GDC][h_AN];
	    wanB = WhatSpecies[Gh_AN];
	    tno2 = Spe_Total_CNO[wanB];
	    Bnum = MP[h_AN_GDC][h_AN];

	    for (j=0; j<tno2; j++){
	      for (i1=1; i1<=NUM1; i1++){
		Residues[spin][Mh_AN_GDC][h_AN][i][j][i1-1] = C[Anum+i][i1]*C[Bnum+j][i1];
	      }
	    }
	  }
	}      
      }
    } /* end of spin */

    /****************************************************
                        free arrays
    ****************************************************/

    for (i=0; i<n2; i++){
      free(S_DC[i]);
    }
    free(S_DC);

    for (spin=0; spin<=SpinP_switch; spin++){
      for (i=0; i<n2; i++){
        free(H_DC[spin][i]);
      }
      free(H_DC[spin]);
    }
    free(H_DC);

    free(ko);
    free(M1);

    for (i=0; i<n2; i++){
      free(C[i]);
    }
    free(C);

    for (h_AN_GDC=0; h_AN_GDC<=MFNAN_GDC[Mc_AN_GDC]; h_AN_GDC++){
      free(MP[h_AN_GDC]);
    }
    free(MP);

    dtime(&Etime_atom);
    time_per_atom[Gc_AN] += Etime_atom - Stime_atom;

  } /* end of Mc_AN_GDC */

  /****************************************************
                calculate projected DOS
  ****************************************************/

  for (spin=0; spin<=SpinP_switch; spin++){
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

      dtime(&Stime_atom);

      Gc_AN = M2G[Mc_AN];
      wanA = WhatSpecies[Gc_AN];
      tno1 = Spe_Total_CNO[wanA];

      for (i1=0; i1<Msize[Mc_AN]; i1++){
        PDOS_DC[spin][Mc_AN][i1] = 0.0;
      }

      for (i=0; i<tno1; i++){
        for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
          Gh_AN = natn[Gc_AN][h_AN];
          wanB = WhatSpecies[Gh_AN];
          tno2 = Spe_Total_CNO[wanB];
          for (j=0; j<tno2; j++){

            tmp1 = OLP0[Mc_AN][h_AN][i][j];
            for (i1=0; i1<Msize[Mc_AN]; i1++){
              PDOS_DC[spin][Mc_AN][i1] += Residues[spin][Mc_AN][h_AN][i][j][i1]*tmp1;
  	    }
          }            
        }        
      }

      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
    }
  }

  /****************************************************
    find chemical potential
  ****************************************************/

  po = 0;
  loopN = 0;

  ChemP_MAX = 15.0;  
  ChemP_MIN =-15.0;
  if      (SpinP_switch==0) spin_degeneracy = 2.0;
  else if (SpinP_switch==1) spin_degeneracy = 1.0;

  do {
    ChemP = 0.50*(ChemP_MAX + ChemP_MIN);

    My_Num_State = 0.0;
    for (spin=0; spin<=SpinP_switch; spin++){
      for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

        dtime(&Stime_atom);

        Gc_AN = M2G[Mc_AN];

        for (i=0; i<Msize[Mc_AN]; i++){
          x = (EVal[spin][Mc_AN][i] - ChemP)*Beta;
          if (x<=-max_x) x = -max_x;
          if (max_x<=x)  x = max_x;
          FermiF = 1.0/(1.0 + exp(x));
          My_Num_State += spin_degeneracy*FermiF*PDOS_DC[spin][Mc_AN][i];
	}

        dtime(&Etime_atom);
        time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
      }
    }

    /* MPI, My_Num_State */

    MPI_Barrier(mpi_comm_level1);
    MPI_Allreduce(&My_Num_State, &Num_State, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

    Dnum = (TZ - Num_State) - system_charge;
    if (0.0<=Dnum) ChemP_MIN = ChemP;
    else           ChemP_MAX = ChemP;
    if (fabs(Dnum)<1.0e-11) po = 1;


    if (myid==Host_ID && 2<=level_stdout){
      printf("ChemP=%15.12f TZ=%15.12f Num_state=%15.12f\n",ChemP,TZ,Num_State); 
    }

    loopN++;
  }
  while (po==0 && loopN<1000); 

  /****************************************************
        eigenenergy by summing up eigenvalues
  ****************************************************/

  My_Eele0[0] = 0.0;
  My_Eele0[1] = 0.0;
  for (spin=0; spin<=SpinP_switch; spin++){
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

      dtime(&Stime_atom);

      Gc_AN = M2G[Mc_AN];
      for (i=0; i<Msize[Mc_AN]; i++){
        x = (EVal[spin][Mc_AN][i] - ChemP)*Beta;
        if (x<=-max_x) x = -max_x;
        if (max_x<=x)  x = max_x;
        FermiF = 1.0/(1.0 + exp(x));
        My_Eele0[spin] += FermiF*EVal[spin][Mc_AN][i]*PDOS_DC[spin][Mc_AN][i];
      }

      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
    }
  }

  /* MPI, My_Eele0 */

  for (spin=0; spin<=SpinP_switch; spin++){
    MPI_Barrier(mpi_comm_level1);
    MPI_Allreduce(&My_Eele0[spin], &Eele0[spin], 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
  }

  if (SpinP_switch==0){
    Eele0[1] = Eele0[0];
  }

  /****************************************************
      calculate density and energy density matrices
  ****************************************************/

  for (spin=0; spin<=SpinP_switch; spin++){
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];
      wanA = WhatSpecies[Gc_AN];
      tno1 = Spe_Total_CNO[wanA];
      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	Gh_AN = natn[Gc_AN][h_AN];
	wanB = WhatSpecies[Gh_AN];
	tno2 = Spe_Total_CNO[wanB];
	for (i=0; i<tno1; i++){
	  for (j=0; j<tno2; j++){
	    CDM[spin][Mc_AN][h_AN][i][j] = 0.0;
	    EDM[spin][Mc_AN][h_AN][i][j] = 0.0;
	  }
	}
      }
    }
  }

  for (spin=0; spin<=SpinP_switch; spin++){
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

      dtime(&Stime_atom);

      Gc_AN = M2G[Mc_AN];
      wanA = WhatSpecies[Gc_AN];
      tno1 = Spe_Total_CNO[wanA];

      for (i1=0; i1<Msize[Mc_AN]; i1++){
	x = (EVal[spin][Mc_AN][i1] - ChemP)*Beta;
	if (x<=-max_x) x = -max_x;
	if (max_x<=x)  x = max_x;
	FermiF = 1.0/(1.0 + exp(x));

        for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
          Gh_AN = natn[Gc_AN][h_AN];
          wanB = WhatSpecies[Gh_AN];
          tno2 = Spe_Total_CNO[wanB];
          for (i=0; i<tno1; i++){
            for (j=0; j<tno2; j++){
              tmp1 = FermiF*Residues[spin][Mc_AN][h_AN][i][j][i1]; 
              CDM[spin][Mc_AN][h_AN][i][j] += tmp1;
              EDM[spin][Mc_AN][h_AN][i][j] += tmp1*EVal[spin][Mc_AN][i1];
  	    }
	  }
        }
      }

      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
    }
  }

  /****************************************************
    freeing of arrays:

  ****************************************************/

  free(Snd_H_Size);
  free(Rcv_H_Size);

  free(Snd_S_Size);
  free(Rcv_S_Size);

  free(Msize);

  for (spin=0; spin<=SpinP_switch; spin++){
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      free(EVal[spin][Mc_AN]);
    }
    free(EVal[spin]);
  }
  free(EVal);

  for (spin=0; spin<=SpinP_switch; spin++){
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];

      if (Mc_AN==0){
        Gc_AN = 0;
        FNAN[0] = 1;
        tno1 = 1;
      }
      else{
        wanA = WhatSpecies[Gc_AN];
        tno1 = Spe_Total_CNO[wanA];
      }

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
        if (Mc_AN==0){
          tno2 = 1;
        }
        else {
          Gh_AN = natn[Gc_AN][h_AN];
          wanB = WhatSpecies[Gh_AN];
          tno2 = Spe_Total_CNO[wanB];
        }

        for (i=0; i<tno1; i++){
          for (j=0; j<tno2; j++){
            free(Residues[spin][Mc_AN][h_AN][i][j]);
	  }
          free(Residues[spin][Mc_AN][h_AN][i]);
        }
        free(Residues[spin][Mc_AN][h_AN]);
      }
      free(Residues[spin][Mc_AN]);
    }
    free(Residues[spin]);
  }
  free(Residues);

  for (spin=0; spin<=SpinP_switch; spin++){
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      free(PDOS_DC[spin][Mc_AN]);
    }
    free(PDOS_DC[spin]);
  }
  free(PDOS_DC);

  /* for time */
  dtime(&TEtime);
  time0 = TEtime - TStime;

  /* for PrintMemory */
  firsttime=0;

  return time0;
}






static double GDC_NonCol(
              int SCF_iter,
              double *****Hks,
              double *****ImNL,
              double ****OLP0,
              double *****CDM,
              double *****EDM,
              double Eele0[2],
              double Eele1[2])
{
  static int firsttime=1;
  int Mc_AN,Mc_AN_GDC,Gc_AN,i,Gi,wan,wanA,wanB,Anum;
  int size1,size2,num,NUM,NUM1,n2,Cwan,Hwan,l1,l2,k1;
  int ih,ig,ian,j,kl,jg,jan,Bnum,m,n,spin,so;
  int l,i1,j1,P_min,m_size;
  int Rni,Rnj,Gi_AN,Gj_AN;
  int h_AN_GDC,Mh_AN_GDC,Gh_AN_GDC;
  int po,po0,loopN,tno1,tno2,h_AN,Gh_AN;
  int k,ii1,jj1;
  double My_TZ,TZ,sum,FermiF,time0;
  double My_Num_State,Num_State,x,Dnum;
  double sum_r,sum_i,tmp1,tmp2;
  double sum1_r,sum1_i,sum2_r,sum2_i;
  double TStime,TEtime;
  double My_Eele0[2],My_Eele1[2];
  double max_x=50.0;
  double ChemP_MAX,ChemP_MIN,spin_degeneracy;
  double **S_DC,*ko,*M1;
  dcomplex **C,**H_DC;
  double **EVal;
  dcomplex ******Residues;
  double **PDOS_DC;
  int **MP,*Msize;
  double *tmp_array;
  double *tmp_array2;
  int *Snd_H_Size,*Rcv_H_Size;
  int *Snd_iHNL_Size,*Rcv_iHNL_Size;
  int *Snd_S_Size,*Rcv_S_Size;
  int numprocs,myid,ID,IDS,IDR,tag=999;
  double Stime_atom, Etime_atom;
  double OLP_eigen_cut=Threshold_OLP_Eigen;

  MPI_Status stat;
  MPI_Request request;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  dtime(&TStime);

  /****************************************************
    allocation of arrays:

    int MP[List_YOUSO[2]];
    int Msize[Matomnum+1];
    double EVal[Matomnum+1][n2];
  ****************************************************/

  Snd_H_Size = (int*)malloc(sizeof(int)*numprocs);
  Rcv_H_Size = (int*)malloc(sizeof(int)*numprocs);
  Snd_iHNL_Size = (int*)malloc(sizeof(int)*numprocs);
  Rcv_iHNL_Size = (int*)malloc(sizeof(int)*numprocs);
  Snd_S_Size = (int*)malloc(sizeof(int)*numprocs);
  Rcv_S_Size = (int*)malloc(sizeof(int)*numprocs);
  Msize = (int*)malloc(sizeof(int)*(Matomnum+1));

  /* find Msize */

  Msize[0] = 1;
  for (Mc_AN_GDC=1; Mc_AN_GDC<=Matomnum_GDC; Mc_AN_GDC++){

    Mc_AN = Mnatn_GDC[Mc_AN_GDC][0];
    Gc_AN = M2G[Mc_AN];

    Anum = 1;
    for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
      Gi = natn[Gc_AN][i];
      wanA = WhatSpecies[Gi];
      Anum = Anum + Spe_Total_CNO[wanA];
    }
    NUM = Anum - 1;

    for (h_AN_GDC=0; h_AN_GDC<=MFNAN_GDC[Mc_AN_GDC]; h_AN_GDC++){
      Mh_AN_GDC = Mnatn_GDC[Mc_AN_GDC][h_AN_GDC];
      Msize[Mh_AN_GDC] = NUM;
    }
  }

  /* allocate EVal */

  m_size = 0;

  EVal = (double**)malloc(sizeof(double*)*(Matomnum+1));
  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];

    if (Mc_AN==0){
      Gc_AN = 0;
      FNAN[0] = 1;
      SNAN[0] = 0;
      n2 = 1;
    }
    else{
      n2 = 2*Msize[Mc_AN] + 3;
    }

    m_size += n2;
    EVal[Mc_AN] = (double*)malloc(sizeof(double)*n2);
  }

  if (firsttime)
  PrintMemory("GDivide_Conquer: EVal",sizeof(double)*m_size,NULL);

  if (2<=level_stdout){
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
        printf("<DC> myid=%i Mc_AN=%2d Gc_AN=%2d Msize=%3d\n",
        myid,Mc_AN,M2G[Mc_AN],Msize[Mc_AN]);
    }
  }

  /****************************************************
    allocation of arrays:

    dcomplex Residues[3]
                     [Matomnum+1]
                     [FNAN[Gc_AN]+1]
                     [Spe_Total_CNO[Gc_AN]] 
                     [Spe_Total_CNO[Gh_AN]] 
                     [n2]
     To reduce the memory size, the size of NUM2
     should be found in the loop.  
  ****************************************************/

  m_size = 0;

  Residues = (dcomplex******)malloc(sizeof(dcomplex*****)*3);
  for (spin=0; spin<3; spin++){
    Residues[spin] = (dcomplex*****)malloc(sizeof(dcomplex****)*(Matomnum+1));
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];

      if (Mc_AN==0){
	Gc_AN = 0;
	FNAN[0] = 1;
	tno1 = 1;
	n2 = 1;
      }
      else{
	wanA = WhatSpecies[Gc_AN];
	tno1 = Spe_Total_CNO[wanA];
	n2 = 2*Msize[Mc_AN] + 2;
      }

      Residues[spin][Mc_AN] = (dcomplex****)malloc(sizeof(dcomplex***)*(FNAN[Gc_AN]+1));

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	if (Mc_AN==0){
	  tno2 = 1;
	}
	else {
	  Gh_AN = natn[Gc_AN][h_AN];
	  wanB = WhatSpecies[Gh_AN];
	  tno2 = Spe_Total_CNO[wanB];
	}

	Residues[spin][Mc_AN][h_AN] = (dcomplex***)malloc(sizeof(dcomplex**)*tno1);
	for (i=0; i<tno1; i++){
	  Residues[spin][Mc_AN][h_AN][i] = (dcomplex**)malloc(sizeof(dcomplex*)*tno2);
	  for (j=0; j<tno2; j++){
	    Residues[spin][Mc_AN][h_AN][i][j] = (dcomplex*)malloc(sizeof(dcomplex)*n2);
	  }
	}

	m_size += tno1*tno2*n2;
      }
    }
  }

  if (firsttime)
  PrintMemory("GDivide_Conquer: Residues",sizeof(dcomplex)*m_size,NULL);

  /****************************************************
    allocation of arrays:

    double PDOS_DC[Matomnum+1]
                  [n2]
  ****************************************************/

  m_size = 0;

  PDOS_DC = (double**)malloc(sizeof(double*)*(Matomnum+1));
  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];

    if (Mc_AN==0) n2 = 1;
    else          n2 = 2*Msize[Mc_AN] + 2;

    m_size += n2;
    PDOS_DC[Mc_AN] = (double*)malloc(sizeof(double)*n2);
  }

  if (firsttime)
  PrintMemory("GDivide_Conquer: PDOS_DC",sizeof(double)*m_size,NULL);

  /****************************************************
   MPI

   Hks
  ****************************************************/

  /***********************************
            set data size
  ************************************/

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){
      tag = 999;

      /* find data size to send block data */
      if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0){

        size1 = 0;
        for (spin=0; spin<=SpinP_switch; spin++){
          for (n=0; n<(F_Snd_Num[IDS]+S_Snd_Num[IDS]); n++){
            Mc_AN = Snd_MAN[IDS][n];
            Gc_AN = Snd_GAN[IDS][n];
            Cwan = WhatSpecies[Gc_AN]; 
            tno1 = Spe_Total_CNO[Cwan];
            for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
              Gh_AN = natn[Gc_AN][h_AN];        
              Hwan = WhatSpecies[Gh_AN];
              tno2 = Spe_Total_CNO[Hwan];
              for (i=0; i<tno1; i++){
                for (j=0; j<tno2; j++){
                  size1++; 
                } 
              } 
	    }
          }
	}
        
        Snd_H_Size[IDS] = size1;
        MPI_Isend(&size1, 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
      }
      else{
        Snd_H_Size[IDS] = 0;
      }

      /* receiving of size of data */

      if ((F_Rcv_Num[IDR]+S_Rcv_Num[IDR])!=0){
        MPI_Recv(&size2, 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
        Rcv_H_Size[IDR] = size2;
      }
      else{
        Rcv_H_Size[IDR] = 0;
      }

      if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0) MPI_Wait(&request,&stat);

    }
    else{
      Snd_H_Size[IDS] = 0;
      Rcv_H_Size[IDR] = 0;
    }
  }

  /***********************************
             data transfer
  ************************************/

  tag = 999;
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){

      /*****************************
              sending of data 
      *****************************/

      if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0){

        size1 = Snd_H_Size[IDS];

        /* allocation of array */

        tmp_array = (double*)malloc(sizeof(double)*size1);

        /* multidimentional array to vector array */

        num = 0;
        for (spin=0; spin<=SpinP_switch; spin++){
          for (n=0; n<(F_Snd_Num[IDS]+S_Snd_Num[IDS]); n++){
            Mc_AN = Snd_MAN[IDS][n];
            Gc_AN = Snd_GAN[IDS][n];
            Cwan = WhatSpecies[Gc_AN]; 
            tno1 = Spe_Total_CNO[Cwan];
            for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
              Gh_AN = natn[Gc_AN][h_AN];
              Hwan = WhatSpecies[Gh_AN];
              tno2 = Spe_Total_CNO[Hwan];
              for (i=0; i<tno1; i++){
                for (j=0; j<tno2; j++){
                  tmp_array[num] = Hks[spin][Mc_AN][h_AN][i][j];
                  num++;
                } 
              } 
	    }
          }
	}

        MPI_Isend(&tmp_array[0], size1, MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);

      }

      /*****************************
         receiving of block data
      *****************************/

      if ((F_Rcv_Num[IDR]+S_Rcv_Num[IDR])!=0){

        size2 = Rcv_H_Size[IDR];
        
        /* allocation of array */
        tmp_array2 = (double*)malloc(sizeof(double)*size2);
        
        MPI_Recv(&tmp_array2[0], size2, MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

        num = 0;
        for (spin=0; spin<=SpinP_switch; spin++){
          Mc_AN = S_TopMAN[IDR] - 1;  /* S_TopMAN should be used. */
          for (n=0; n<(F_Rcv_Num[IDR]+S_Rcv_Num[IDR]); n++){
            Mc_AN++;
            Gc_AN = Rcv_GAN[IDR][n];
            Cwan = WhatSpecies[Gc_AN]; 
            tno1 = Spe_Total_CNO[Cwan];

            for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
              Gh_AN = natn[Gc_AN][h_AN];        
              Hwan = WhatSpecies[Gh_AN];
              tno2 = Spe_Total_CNO[Hwan];
              for (i=0; i<tno1; i++){
                for (j=0; j<tno2; j++){
                  Hks[spin][Mc_AN][h_AN][i][j] = tmp_array2[num];
                  num++;
		}
	      }
	    }
	  }        
	}

        /* freeing of array */
        free(tmp_array2);
      }

      if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0){
        MPI_Wait(&request,&stat);
        free(tmp_array); /* freeing of array */
      }
    }
  }

  /****************************************************
   MPI

   ImNL
  ****************************************************/

  /***********************************
             set data size
  ************************************/

  /* spin-orbit coupling or LDA+U */  
  if (SO_switch==1 || Hub_U_switch==1 || Constraint_NCS_switch==1 || Zeeman_NCS_switch==1 || Zeeman_NCO_switch==1){

    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      if (ID!=0){
	tag = 999;

	/* find data size to send block data */
	if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0){

	  size1 = 0;
	  for (so=0; so<List_YOUSO[5]; so++){
	    for (n=0; n<(F_Snd_Num[IDS]+S_Snd_Num[IDS]); n++){
	      Mc_AN = Snd_MAN[IDS][n];
	      Gc_AN = Snd_GAN[IDS][n];
	      Cwan = WhatSpecies[Gc_AN]; 
	      tno1 = Spe_Total_CNO[Cwan];
	      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
		Gh_AN = natn[Gc_AN][h_AN];        
		Hwan = WhatSpecies[Gh_AN];
		tno2 = Spe_Total_CNO[Hwan];
		for (i=0; i<tno1; i++){
		  for (j=0; j<tno2; j++){
		    size1++; 
		  } 
		} 
	      }
	    }
	  }
 
	  Snd_iHNL_Size[IDS] = size1;
	  MPI_Isend(&size1, 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
	}
	else{
	  Snd_iHNL_Size[IDS] = 0;
	}

	/* receiving of size of data */

	if ((F_Rcv_Num[IDR]+S_Rcv_Num[IDR])!=0){
	  MPI_Recv(&size2, 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
	  Rcv_iHNL_Size[IDR] = size2;
	}
	else{
	  Rcv_iHNL_Size[IDR] = 0;
	}

	if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0) MPI_Wait(&request,&stat);

      }
      else{
	Snd_iHNL_Size[IDS] = 0;
	Rcv_iHNL_Size[IDR] = 0;
      }
    }

    /***********************************
               data transfer
    ************************************/

    tag = 999;
    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      if (ID!=0){

	/*****************************
	       sending of data 
	*****************************/

	if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0){

	  size1 = Snd_iHNL_Size[IDS];

	  /* allocation of array */

	  tmp_array = (double*)malloc(sizeof(double)*size1);

	  /* multidimentional array to vector array */

	  num = 0;
	  for (so=0; so<List_YOUSO[5]; so++){
	    for (n=0; n<(F_Snd_Num[IDS]+S_Snd_Num[IDS]); n++){
	      Mc_AN = Snd_MAN[IDS][n];
	      Gc_AN = Snd_GAN[IDS][n];
	      Cwan = WhatSpecies[Gc_AN]; 
	      tno1 = Spe_Total_CNO[Cwan];
	      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
		Gh_AN = natn[Gc_AN][h_AN];
		Hwan = WhatSpecies[Gh_AN];
		tno2 = Spe_Total_CNO[Hwan];
		for (i=0; i<tno1; i++){
		  for (j=0; j<tno2; j++){
		    tmp_array[num] = ImNL[so][Mc_AN][h_AN][i][j];
		    num++;
		  } 
		} 
	      }
	    }
	  }

	  MPI_Isend(&tmp_array[0], size1, MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);

	}

	/*****************************
           receiving of block data
	*****************************/

	if ((F_Rcv_Num[IDR]+S_Rcv_Num[IDR])!=0){

	  size2 = Rcv_iHNL_Size[IDR];
        
	  /* allocation of array */
	  tmp_array2 = (double*)malloc(sizeof(double)*size2);
        
	  MPI_Recv(&tmp_array2[0], size2, MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

	  num = 0;
	  for (so=0; so<List_YOUSO[5]; so++){
	    Mc_AN = S_TopMAN[IDR] - 1;  /* S_TopMAN should be used. */
	    for (n=0; n<(F_Rcv_Num[IDR]+S_Rcv_Num[IDR]); n++){
	      Mc_AN++;
	      Gc_AN = Rcv_GAN[IDR][n];
	      Cwan = WhatSpecies[Gc_AN]; 
	      tno1 = Spe_Total_CNO[Cwan];

	      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
		Gh_AN = natn[Gc_AN][h_AN];        
		Hwan = WhatSpecies[Gh_AN];
		tno2 = Spe_Total_CNO[Hwan];
		for (i=0; i<tno1; i++){
		  for (j=0; j<tno2; j++){
		    ImNL[so][Mc_AN][h_AN][i][j] = tmp_array2[num];
		    num++;
		  }
		}
	      }
	    }        
	  }

	  /* freeing of array */
	  free(tmp_array2);
	}

	if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0){
	  MPI_Wait(&request,&stat);
	  free(tmp_array); /* freeing of array */
	}
      }
    }
  }

  /****************************************************
   MPI

   OLP0
  ****************************************************/

  /***********************************
             set data size
  ************************************/

  if (SCF_iter<=2){

    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      if (ID!=0){
        tag = 999;

        /* find data size to send block data */
        if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0){
          size1 = 0;
          for (n=0; n<(F_Snd_Num[IDS]+S_Snd_Num[IDS]); n++){
            Mc_AN = Snd_MAN[IDS][n];
            Gc_AN = Snd_GAN[IDS][n];
            Cwan = WhatSpecies[Gc_AN]; 
            tno1 = Spe_Total_CNO[Cwan];
            for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
              Gh_AN = natn[Gc_AN][h_AN];        
              Hwan = WhatSpecies[Gh_AN];
              tno2 = Spe_Total_CNO[Hwan];
              for (i=0; i<tno1; i++){
                for (j=0; j<tno2; j++){
                  size1++; 
                } 
              } 
	    }
          }

          Snd_S_Size[IDS] = size1;
          MPI_Isend(&size1, 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
        }
        else{
          Snd_S_Size[IDS] = 0;
        }

        /* receiving of size of data */
 
        if ((F_Rcv_Num[IDR]+S_Rcv_Num[IDR])!=0){
          MPI_Recv(&size2, 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
          Rcv_S_Size[IDR] = size2;
        }
        else{
          Rcv_S_Size[IDR] = 0;
        }

        if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0) MPI_Wait(&request,&stat);
      }
      else{
        Snd_S_Size[IDS] = 0;
        Rcv_S_Size[IDR] = 0;
      }
    }

    /***********************************
               data transfer
    ************************************/

    tag = 999;
    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      if (ID!=0){

        /*****************************
                sending of data 
        *****************************/

        if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0){

          size1 = Snd_S_Size[IDS];

          /* allocation of array */

          tmp_array = (double*)malloc(sizeof(double)*size1);

          /* multidimentional array to vector array */

          num = 0;

          for (n=0; n<(F_Snd_Num[IDS]+S_Snd_Num[IDS]); n++){
            Mc_AN = Snd_MAN[IDS][n];
            Gc_AN = Snd_GAN[IDS][n];
            Cwan = WhatSpecies[Gc_AN]; 
            tno1 = Spe_Total_CNO[Cwan];
            for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
              Gh_AN = natn[Gc_AN][h_AN];        
              Hwan = WhatSpecies[Gh_AN];
              tno2 = Spe_Total_CNO[Hwan];
              for (i=0; i<tno1; i++){
                for (j=0; j<tno2; j++){
                  tmp_array[num] = OLP0[Mc_AN][h_AN][i][j];
                  num++;
                } 
              } 
	    }
          }

          MPI_Isend(&tmp_array[0], size1, MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);
        }

        /*****************************
           receiving of block data
        *****************************/

        if ((F_Rcv_Num[IDR]+S_Rcv_Num[IDR])!=0){
          
          size2 = Rcv_S_Size[IDR];
        
         /* allocation of array */
          tmp_array2 = (double*)malloc(sizeof(double)*size2);
         
          MPI_Recv(&tmp_array2[0], size2, MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

          num = 0;
          Mc_AN = S_TopMAN[IDR] - 1; /* S_TopMAN should be used. */
          for (n=0; n<(F_Rcv_Num[IDR]+S_Rcv_Num[IDR]); n++){
            Mc_AN++;
            Gc_AN = Rcv_GAN[IDR][n];
            Cwan = WhatSpecies[Gc_AN]; 
            tno1 = Spe_Total_CNO[Cwan];

            for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
              Gh_AN = natn[Gc_AN][h_AN];        
              Hwan = WhatSpecies[Gh_AN];
              tno2 = Spe_Total_CNO[Hwan];
              for (i=0; i<tno1; i++){
                for (j=0; j<tno2; j++){
                  OLP0[Mc_AN][h_AN][i][j] = tmp_array2[num];
                  num++;
   	        }
	      }
	    }
	  }        

          /* freeing of array */
          free(tmp_array2);
        }

        if ((F_Snd_Num[IDS]+S_Snd_Num[IDS])!=0){
          MPI_Wait(&request,&stat);
          free(tmp_array); /* freeing of array */
	}
      }
    }
  }

  /****************************************************
            find the total number of electrons 
  ****************************************************/

  My_TZ = 0.0;
  for (i=1; i<=Matomnum; i++){
    Gc_AN = M2G[i];
    wan = WhatSpecies[Gc_AN];
    My_TZ = My_TZ + Spe_Core_Charge[wan];
  }

  /* MPI, My_TZ */

  MPI_Barrier(mpi_comm_level1);
  MPI_Allreduce(&My_TZ, &TZ, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

  /****************************************************
      Setting of Hamiltonian and overlap matrices
  ****************************************************/

  for (Mc_AN_GDC=1; Mc_AN_GDC<=Matomnum_GDC; Mc_AN_GDC++){

    dtime(&Stime_atom);

    Mc_AN = Mnatn_GDC[Mc_AN_GDC][0];
    Gc_AN = M2G[Mc_AN];
    wan = WhatSpecies[Gc_AN];

    /***********************************************
      find the size of matrix for the atom Mc_AN
                and set the MP vector

     Note:
         MP indicates the starting position of
              atom i in arraies H and S

     allocation: 
       int MP[(MFNAN_GDC[Mc_AN_GDC]+1)]
                    [List_YOUSO[2]];
    ***********************************************/

    MP = (int**)malloc(sizeof(int*)*(MFNAN_GDC[Mc_AN_GDC]+1));
    for (h_AN_GDC=0; h_AN_GDC<=MFNAN_GDC[Mc_AN_GDC]; h_AN_GDC++){
      MP[h_AN_GDC] = (int*)malloc(sizeof(int)*List_YOUSO[2]);
    }

    /* set MP[0][i] */

    Anum = 1;
    for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
      MP[0][i] = Anum;
      Gi = natn[Gc_AN][i];
      wanA = WhatSpecies[Gi];
      Anum = Anum + Spe_Total_CNO[wanA];
    }
    NUM = Anum - 1;
    n2 = 2*NUM + 3;

    /* set MP[h_AN_GDC][j] */
    for (h_AN_GDC=1; h_AN_GDC<=MFNAN_GDC[Mc_AN_GDC]; h_AN_GDC++){

      Mh_AN_GDC = Mnatn_GDC[Mc_AN_GDC][h_AN_GDC];
      Gh_AN_GDC = M2G[Mh_AN_GDC];

      for (j=0; j<=FNAN[Gh_AN_GDC]; j++){

        Gj_AN = natn[Gh_AN_GDC][j];
        Rnj   =  ncn[Gh_AN_GDC][j];

        po0= 0;
        i = 0;
        do{

          Gi_AN = natn[Gc_AN][i];
          Rni   =  ncn[Gc_AN][i];

          if (Gj_AN==Gi_AN && Rnj==Rni){
            MP[h_AN_GDC][j] = MP[0][i];
            po0 = 1;
          }

          i++;
	} while ( i<=(FNAN[Gc_AN]+SNAN[Gc_AN]) && po0==0 );
      }
    }

    /***********************************************
     allocation of arrays:
     
     double   S_DC[NUM+2][NUM+2];     
     dcomplex H_DC[n2][n2];     
     double   ko[n2];
     double   M1[n2];
     dcomplex C[n2][n2];     
    ***********************************************/

    S_DC = (double**)malloc(sizeof(double*)*(NUM+2));
    for (i=0; i<(NUM+2); i++){
      S_DC[i] = (double*)malloc(sizeof(double)*(NUM+2));
    }

    H_DC = (dcomplex**)malloc(sizeof(dcomplex*)*n2);
    for (i=0; i<n2; i++){
      H_DC[i] = (dcomplex*)malloc(sizeof(dcomplex)*n2);
    }

    ko = (double*)malloc(sizeof(double)*n2);
    M1 = (double*)malloc(sizeof(double)*n2);

    C = (dcomplex**)malloc(sizeof(dcomplex*)*n2);
    for (i=0; i<n2; i++){
      C[i] = (dcomplex*)malloc(sizeof(dcomplex)*n2);
    }

    /***********************************************
     construct cluster full matrices of Hamiltonian
             and overlap for the atom Mc_AN             
    ***********************************************/

    for (i=0; i<=(FNAN[Gc_AN]+SNAN[Gc_AN]); i++){
      ig = natn[Gc_AN][i];
      ian = Spe_Total_CNO[WhatSpecies[ig]];
      Anum = MP[0][i];
      ih = S_G2M[ig];

      for (j=0; j<=(FNAN[Gc_AN]+SNAN[Gc_AN]); j++){

	kl = RMI1[Mc_AN][i][j];
	jg = natn[Gc_AN][j];
	jan = Spe_Total_CNO[WhatSpecies[jg]];
	Bnum = MP[0][j];

	if (0<=kl){

          if (SCF_iter<=2){
  	    for (m=0; m<ian; m++){
	      for (n=0; n<jan; n++){
	        S_DC[Anum+m][Bnum+n] = OLP0[ih][kl][m][n];
	      }
	    }
	  }

          /* non-spin-orbit coupling and non-LDA+U */  
          if (SO_switch==0 && Hub_U_switch==0 && Constraint_NCS_switch==0 && Zeeman_NCS_switch==0 && Zeeman_NCO_switch==0){
  	    for (m=0; m<ian; m++){
	      for (n=0; n<jan; n++){
	        H_DC[Anum+m    ][Bnum+n    ].r =  Hks[0][ih][kl][m][n];
	        H_DC[Anum+m    ][Bnum+n    ].i =  0.0;
	        H_DC[Anum+m+NUM][Bnum+n+NUM].r =  Hks[1][ih][kl][m][n];
	        H_DC[Anum+m+NUM][Bnum+n+NUM].i =  0.0;
	        H_DC[Anum+m    ][Bnum+n+NUM].r =  Hks[2][ih][kl][m][n];
	        H_DC[Anum+m    ][Bnum+n+NUM].i =  Hks[3][ih][kl][m][n];
                H_DC[Bnum+n+NUM][Anum+m    ].r =  H_DC[Anum+m    ][Bnum+n+NUM].r; 
	        H_DC[Bnum+n+NUM][Anum+m    ].i = -H_DC[Anum+m    ][Bnum+n+NUM].i;
	      }
	    }
	  }

          /* spin-orbit coupling or LDA+U */  
          else {
  	    for (m=0; m<ian; m++){
	      for (n=0; n<jan; n++){
	        H_DC[Anum+m    ][Bnum+n    ].r =  Hks[0][ih][kl][m][n];
	        H_DC[Anum+m    ][Bnum+n    ].i =  ImNL[0][ih][kl][m][n];
	        H_DC[Anum+m+NUM][Bnum+n+NUM].r =  Hks[1][ih][kl][m][n];
	        H_DC[Anum+m+NUM][Bnum+n+NUM].i =  ImNL[1][ih][kl][m][n];
	        H_DC[Anum+m    ][Bnum+n+NUM].r =  Hks[2][ih][kl][m][n];
	        H_DC[Anum+m    ][Bnum+n+NUM].i =  Hks[3][ih][kl][m][n] + ImNL[2][ih][kl][m][n];
                H_DC[Bnum+n+NUM][Anum+m    ].r =  H_DC[Anum+m    ][Bnum+n+NUM].r; 
	        H_DC[Bnum+n+NUM][Anum+m    ].i = -H_DC[Anum+m    ][Bnum+n+NUM].i;
	      }
	    }
	  }

	}

	else{

          if (SCF_iter<=2){
	    for (m=0; m<ian; m++){
	      for (n=0; n<jan; n++){
	        S_DC[Anum+m][Bnum+n] = 0.0;
	      }
	    }
	  }

	  for (m=0; m<ian; m++){
	    for (n=0; n<jan; n++){
	      H_DC[Anum+m    ][Bnum+n    ].r = 0.0;
	      H_DC[Anum+m    ][Bnum+n    ].i = 0.0;
	      H_DC[Anum+m+NUM][Bnum+n+NUM].r = 0.0;
	      H_DC[Anum+m+NUM][Bnum+n+NUM].i = 0.0;
	      H_DC[Anum+m    ][Bnum+n+NUM].r = 0.0;
	      H_DC[Anum+m    ][Bnum+n+NUM].i = 0.0;
              H_DC[Bnum+n+NUM][Anum+m    ].r = 0.0;
	      H_DC[Bnum+n+NUM][Anum+m    ].i = 0.0;
	    }
	  }

	}
      }
    }

    /****************************************************
     Solve the generalized eigenvalue problem
     HC = SCE

     1) diagonalize S
     2) search negative eigenvalues of S  
    ****************************************************/

    if (SCF_iter<=2){

      Eigen_lapack(S_DC, ko, NUM, NUM);

      /***********************************************
             Searching of negative eigenvalues
      ************************************************/

      P_min = 1;
      for (l=1; l<=NUM; l++){
        if (ko[l]<OLP_eigen_cut){
          P_min = l + 1;
          if (3<=level_stdout){
            printf("<DC>  Negative EV of OLP %2d %15.12f\n",l,ko[l]);
	  }
        }
      }

      S12[Mc_AN_GDC][0][0] = P_min;

      for (l=1; l<P_min; l++)     M1[l] = 0.0;
      for (l=P_min; l<=NUM; l++)  M1[l] = 1.0/sqrt(ko[l]);

      for (i1=1; i1<=NUM; i1++){
	for (j1=1; j1<=NUM; j1++){
	  S_DC[i1][j1]       = S_DC[i1][j1]*M1[j1];
          S12[Mc_AN_GDC][i1][j1] = S_DC[i1][j1];
	}
      }

    }

    else{

      P_min = (int)S12[Mc_AN_GDC][0][0];

      for (i1=1; i1<=NUM; i1++){
	for (j1=1; j1<=NUM; j1++){
	  S_DC[i1][j1] = S12[Mc_AN_GDC][i1][j1];
	}
      }
    }

    /***********************************************
      transform Hamiltonian matrix
    ************************************************/

    /* transpose S */

    for (i1=1; i1<=NUM; i1++){
      for (j1=i1+1; j1<=NUM; j1++){
	tmp1 = S_DC[i1][j1];
	tmp2 = S_DC[j1][i1];
	S_DC[i1][j1] = tmp2;
	S_DC[j1][i1] = tmp1;
      }
    }

    /* H * U * M1 */

    for (j1=P_min; j1<=NUM; j1++){
      for (k=0; k<=1; k++){
        jj1 = 2*j1 - P_min + k;
        k1 = k*NUM;
        
        for (i1=1; i1<=2*NUM; i1++){

	  sum_r = 0.0;
	  sum_i = 0.0;

	  for (l=1; l<=NUM; l++){
            l1 = k1 + l;
	    sum_r += H_DC[i1][l1].r*S_DC[j1][l];
	    sum_i += H_DC[i1][l1].i*S_DC[j1][l];
	  }

	  C[jj1][i1].r = sum_r;
	  C[jj1][i1].i = sum_i;
	}
      }
    }

    /* M1 * U^+ H * U * M1 */

    for (j1=1; j1<=2*NUM; j1++){
      for (i1=P_min; i1<=NUM; i1++){

        sum1_r = 0.0;
	sum1_i = 0.0;
	sum2_r = 0.0;
	sum2_i = 0.0;

	for (l=1; l<=NUM; l++){
	  sum1_r += S_DC[i1][l]*C[j1][l].r;
	  sum1_i += S_DC[i1][l]*C[j1][l].i;
	}

	for (l=NUM+1; l<=2*NUM; l++){
          l1 = l - NUM; 
	  sum2_r += S_DC[i1][l1]*C[j1][l].r;
	  sum2_i += S_DC[i1][l1]*C[j1][l].i;
	}

	ii1 = 2*i1 - P_min;
	H_DC[j1][ii1].r = sum1_r;
	H_DC[j1][ii1].i = sum1_i;

	ii1 = 2*i1 - P_min + 1;
	H_DC[j1][ii1].r = sum2_r;
	H_DC[j1][ii1].i = sum2_i;
      }
    }

    /* H to C (transposition) */

    for (i1=P_min; i1<=2*NUM; i1++){
      for (j1=P_min; j1<=2*NUM; j1++){
	C[j1-(P_min-1)][i1-(P_min-1)].r = H_DC[i1][j1].r;
	C[j1-(P_min-1)][i1-(P_min-1)].i = H_DC[i1][j1].i;
      }
    }

    /***********************************************
     diagonalize the trasformed Hamiltonian matrix
    ************************************************/

    NUM1 = 2*NUM - (P_min - 1);
    EigenBand_lapack(C, ko, NUM1, 1);

    /* C to H (transposition) */

    for (i1=1; i1<=NUM1; i1++){
      for (j1=1; j1<=NUM1; j1++){
        H_DC[j1][i1].r = C[i1][j1].r;
        H_DC[j1][i1].i = C[i1][j1].i;
      }
    }

    /***********************************************
      transformation to the original eigen vectors.
      NOTE 244P    C = U * lambda^{-1/2} * D
    ***********************************************/

    /* transpose S */

    for (i1=1; i1<=NUM; i1++){
      for (j1=i1+1; j1<=NUM; j1++){
        tmp1 = S_DC[i1][j1];
        tmp2 = S_DC[j1][i1];
        S_DC[i1][j1] = tmp2;
        S_DC[j1][i1] = tmp1;
      }
    }

    for (j1=1; j1<=NUM1; j1++){
      for (i1=1; i1<=NUM; i1++){

	sum1_r = 0.0;
	sum1_i = 0.0;
	sum2_r = 0.0;
	sum2_i = 0.0;

	for (l=P_min; l<=NUM; l++){
	  l1 = 2*(l - P_min) + 1;
	  l2 = l1 + 1;
	  sum1_r += S_DC[i1][l]*H_DC[j1][l1].r;
	  sum2_r += S_DC[i1][l]*H_DC[j1][l2].r;
	  sum1_i += S_DC[i1][l]*H_DC[j1][l1].i;
	  sum2_i += S_DC[i1][l]*H_DC[j1][l2].i;
	}

	C[i1    ][j1].r = sum1_r;
	C[i1    ][j1].i = sum1_i;
	C[i1+NUM][j1].r = sum2_r;
	C[i1+NUM][j1].i = sum2_i;
      }
    }

    /***********************************************
        store eigenvalues and residues of poles
    ***********************************************/

    for (h_AN_GDC=0; h_AN_GDC<=MFNAN_GDC[Mc_AN_GDC]; h_AN_GDC++){

      Mh_AN_GDC = Mnatn_GDC[Mc_AN_GDC][h_AN_GDC];
      Gh_AN_GDC = M2G[Mh_AN_GDC];
      wanA = WhatSpecies[Gh_AN_GDC];
      tno1 = Spe_Total_CNO[wanA];
      Anum = MP[h_AN_GDC][0];

      for (i1=1; i1<=2*NUM; i1++){
        EVal[Mh_AN_GDC][i1-1] = 1000.0;
      }
      for (i1=1; i1<=NUM1; i1++){
        EVal[Mh_AN_GDC][i1-1] = ko[i1];
      }

      for (i=0; i<tno1; i++){
        for (h_AN=0; h_AN<=FNAN[Gh_AN_GDC]; h_AN++){
          Gh_AN = natn[Gh_AN_GDC][h_AN];
  	  wanB = WhatSpecies[Gh_AN];
	  tno2 = Spe_Total_CNO[wanB];
          Bnum = MP[h_AN_GDC][h_AN];

	  for (j=0; j<tno2; j++){
	    for (i1=1; i1<=NUM1; i1++){

              /* Re11 */
	      Residues[0][Mh_AN_GDC][h_AN][i][j][i1-1].r = C[Anum+i    ][i1].r*C[Bnum+j    ][i1].r
                                                          +C[Anum+i    ][i1].i*C[Bnum+j    ][i1].i;

              /* Re22 */
	      Residues[1][Mh_AN_GDC][h_AN][i][j][i1-1].r = C[Anum+i+NUM][i1].r*C[Bnum+j+NUM][i1].r
                                                          +C[Anum+i+NUM][i1].i*C[Bnum+j+NUM][i1].i;

              /* Re12 */
	      Residues[2][Mh_AN_GDC][h_AN][i][j][i1-1].r = C[Anum+i    ][i1].r*C[Bnum+j+NUM][i1].r
                                                          +C[Anum+i    ][i1].i*C[Bnum+j+NUM][i1].i;

              /* Im12 */
	      Residues[2][Mh_AN_GDC][h_AN][i][j][i1-1].i = C[Anum+i    ][i1].r*C[Bnum+j+NUM][i1].i 
                                                          -C[Anum+i    ][i1].i*C[Bnum+j+NUM][i1].r; 

              /* spin-orbit coupling or LDA+U */
              if (SO_switch==1 || Hub_U_switch==1 || Constraint_NCS_switch==1 
                 || Zeeman_NCS_switch==1 || Zeeman_NCO_switch==1){

                /* Im11 */
	        Residues[0][Mh_AN_GDC][h_AN][i][j][i1-1].i = C[Anum+i    ][i1].r*C[Bnum+j    ][i1].i
                                                            -C[Anum+i    ][i1].i*C[Bnum+j    ][i1].r;
                /* Im22 */
	        Residues[1][Mh_AN_GDC][h_AN][i][j][i1-1].i = C[Anum+i+NUM][i1].r*C[Bnum+j+NUM][i1].i
                                                            -C[Anum+i+NUM][i1].i*C[Bnum+j+NUM][i1].r;
              } 

	    }
	  }
        }
      }
    }      

    /****************************************************
                        free arrays
    ****************************************************/

    for (i=0; i<(NUM+2); i++){
      free(S_DC[i]);
    }
    free(S_DC);

    for (i=0; i<n2; i++){
      free(H_DC[i]);
    }
    free(H_DC);

    free(ko);
    free(M1);

    for (i=0; i<n2; i++){
      free(C[i]);
    }
    free(C);

    for (h_AN_GDC=0; h_AN_GDC<=MFNAN_GDC[Mc_AN_GDC]; h_AN_GDC++){
      free(MP[h_AN_GDC]);
    }
    free(MP);

    dtime(&Etime_atom);
    time_per_atom[Gc_AN] += Etime_atom - Stime_atom;

    /*
    for (i1=1; i1<=NUM1; i1++){
      printf("Mc_AN=%2d i1=%2d Eval=%15.12f\n",Mc_AN,i1,EVal[Mc_AN][i1-1]);
    }
    */

  } /* end of Mc_AN */

  /****************************************************
              calculate projected DOS
  ****************************************************/

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    dtime(&Stime_atom);

    Gc_AN = M2G[Mc_AN];
    wanA = WhatSpecies[Gc_AN];
    tno1 = Spe_Total_CNO[wanA];

    for (i1=0; i1<2*Msize[Mc_AN]; i1++){
      PDOS_DC[Mc_AN][i1] = 0.0;
    }

    for (i=0; i<tno1; i++){
      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	Gh_AN = natn[Gc_AN][h_AN];
	wanB = WhatSpecies[Gh_AN];
	tno2 = Spe_Total_CNO[wanB];
	for (j=0; j<tno2; j++){

          tmp1 = OLP0[Mc_AN][h_AN][i][j];
	  for (i1=0; i1<2*Msize[Mc_AN]; i1++){
	    PDOS_DC[Mc_AN][i1] += ( Residues[0][Mc_AN][h_AN][i][j][i1].r
			          + Residues[1][Mc_AN][h_AN][i][j][i1].r )*tmp1;
	  }

	}            
      }        
    }

    dtime(&Etime_atom);
    time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
  }

  /****************************************************
    find chemical potential
  ****************************************************/

  po = 0;
  loopN = 0;

  ChemP_MAX = 15.0;  
  ChemP_MIN =-15.0;

  do {
    ChemP = 0.50*(ChemP_MAX + ChemP_MIN);

    My_Num_State = 0.0;

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

      dtime(&Stime_atom);

      Gc_AN = M2G[Mc_AN];

      for (i=0; i<2*Msize[Mc_AN]; i++){
	x = (EVal[Mc_AN][i] - ChemP)*Beta;
	if (x<=-max_x) x = -max_x;
	if (max_x<=x)  x = max_x;
	FermiF = 1.0/(1.0 + exp(x));
	My_Num_State += FermiF*PDOS_DC[Mc_AN][i];
      }

      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
    }

    /* MPI, My_Num_State */

    MPI_Barrier(mpi_comm_level1);
    MPI_Allreduce(&My_Num_State, &Num_State, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

    Dnum = (TZ - Num_State) - system_charge;
    if (0.0<=Dnum) ChemP_MIN = ChemP;
    else           ChemP_MAX = ChemP;
    if (fabs(Dnum)<1.0e-11) po = 1;

    if (myid==Host_ID && 2<=level_stdout){
      printf("ChemP=%15.12f TZ=%15.12f Num_state=%15.12f\n",ChemP,TZ,Num_State); 
    }

    loopN++;
  }
  while (po==0 && loopN<1000); 

  /****************************************************
        eigenenergy by summing up eigenvalues
  ****************************************************/

  My_Eele0[0] = 0.0;
  My_Eele0[1] = 0.0;

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    dtime(&Stime_atom);

    Gc_AN = M2G[Mc_AN];
    for (i=0; i<2*Msize[Mc_AN]; i++){
      x = (EVal[Mc_AN][i] - ChemP)*Beta;
      if (x<=-max_x) x = -max_x;
      if (max_x<=x)  x = max_x;
      FermiF = 1.0/(1.0 + exp(x));
      My_Eele0[0] += FermiF*EVal[Mc_AN][i]*PDOS_DC[Mc_AN][i];
    }

    dtime(&Etime_atom);
    time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
  }

  /* MPI, My_Eele0 */

  MPI_Barrier(mpi_comm_level1);
  MPI_Allreduce(&My_Eele0[0], &Eele0[0], 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
  MPI_Allreduce(&My_Eele0[1], &Eele0[1], 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

  /****************************************************
      calculate density and energy density matrices

        CDM[0]  Re alpha alpha density matrix
        CDM[1]  Re beta  beta  density matrix
        CDM[2]  Re alpha beta  density matrix
        CDM[3]  Im alpha beta  density matrix
        iDM[0][0]  Im alpha alpha density matrix
        iDM[0][1]  Im beta  beta  density matrix

        EDM[0]  Re alpha alpha energy density matrix
        EDM[1]  Re beta  beta  energy density matrix
        EDM[2]  Re alpha beta  energy density matrix
        EDM[3]  Im alpha beta  energy density matrix
  ****************************************************/

  for (spin=0; spin<=SpinP_switch; spin++){
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];
      wanA = WhatSpecies[Gc_AN];
      tno1 = Spe_Total_CNO[wanA];
      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	Gh_AN = natn[Gc_AN][h_AN];
	wanB = WhatSpecies[Gh_AN];
	tno2 = Spe_Total_CNO[wanB];

	for (i=0; i<tno1; i++){
	  for (j=0; j<tno2; j++){
	    CDM[spin][Mc_AN][h_AN][i][j] = 0.0;
	    EDM[spin][Mc_AN][h_AN][i][j] = 0.0;
	  }
	}

        /* spin-orbit coupling or LDA+U */
        if ( (SO_switch==1 || Hub_U_switch==1 || Constraint_NCS_switch==1 
             || Zeeman_NCS_switch==1 || Zeeman_NCO_switch==1) && spin==0 ){
	  for (i=0; i<tno1; i++){
	    for (j=0; j<tno2; j++){
	      iDM[0][0][Mc_AN][h_AN][i][j] = 0.0;
	      iDM[0][1][Mc_AN][h_AN][i][j] = 0.0;
	    }
	  }
	}

      }
    }
  }

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    dtime(&Stime_atom);

    Gc_AN = M2G[Mc_AN];
    wanA = WhatSpecies[Gc_AN];
    tno1 = Spe_Total_CNO[wanA];

    for (i1=0; i1<2*Msize[Mc_AN]; i1++){
      x = (EVal[Mc_AN][i1] - ChemP)*Beta;
      if (x<=-max_x) x = -max_x;
      if (max_x<=x)  x = max_x;
      FermiF = 1.0/(1.0 + exp(x));

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	Gh_AN = natn[Gc_AN][h_AN];
	wanB = WhatSpecies[Gh_AN];
	tno2 = Spe_Total_CNO[wanB];
	for (i=0; i<tno1; i++){
	  for (j=0; j<tno2; j++){

            tmp1 = FermiF*Residues[0][Mc_AN][h_AN][i][j][i1].r;

	    CDM[0][Mc_AN][h_AN][i][j] += tmp1;
	    EDM[0][Mc_AN][h_AN][i][j] += tmp1*EVal[Mc_AN][i1];

            tmp1 = FermiF*Residues[1][Mc_AN][h_AN][i][j][i1].r;

	    CDM[1][Mc_AN][h_AN][i][j] += tmp1;
	    EDM[1][Mc_AN][h_AN][i][j] += tmp1*EVal[Mc_AN][i1];

            tmp1 = FermiF*Residues[2][Mc_AN][h_AN][i][j][i1].r;

	    CDM[2][Mc_AN][h_AN][i][j] += tmp1;
	    EDM[2][Mc_AN][h_AN][i][j] += tmp1*EVal[Mc_AN][i1];

            tmp1 = FermiF*Residues[2][Mc_AN][h_AN][i][j][i1].i; 

	    CDM[3][Mc_AN][h_AN][i][j] += tmp1;
	    EDM[3][Mc_AN][h_AN][i][j] += tmp1*EVal[Mc_AN][i1];

            /* spin-orbit coupling or LDA+U */
            if (SO_switch==1 || Hub_U_switch==1 || Constraint_NCS_switch==1 
             || Zeeman_NCS_switch==1 || Zeeman_NCO_switch==1){
	      iDM[0][0][Mc_AN][h_AN][i][j] += FermiF*Residues[0][Mc_AN][h_AN][i][j][i1].i;
	      iDM[0][1][Mc_AN][h_AN][i][j] += FermiF*Residues[1][Mc_AN][h_AN][i][j][i1].i;
            }

	  }
	}
      }
    }

    dtime(&Etime_atom);
    time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
  }

  /****************************************************
    freeing of arrays:

  ****************************************************/

  free(Snd_H_Size);
  free(Rcv_H_Size);

  free(Snd_iHNL_Size);
  free(Rcv_iHNL_Size);

  free(Snd_S_Size);
  free(Rcv_S_Size);

  free(Msize);

  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
    free(EVal[Mc_AN]);
  }
  free(EVal);

  for (spin=0; spin<3; spin++){
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];

      if (Mc_AN==0){
	Gc_AN = 0;
	FNAN[0] = 1;
	tno1 = 1;
      }
      else{
	wanA = WhatSpecies[Gc_AN];
	tno1 = Spe_Total_CNO[wanA];
      }

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	if (Mc_AN==0){
	  tno2 = 1;
	}
	else {
	  Gh_AN = natn[Gc_AN][h_AN];
	  wanB = WhatSpecies[Gh_AN];
	  tno2 = Spe_Total_CNO[wanB];
	}

	for (i=0; i<tno1; i++){
	  for (j=0; j<tno2; j++){
	    free(Residues[spin][Mc_AN][h_AN][i][j]);
	  }
	  free(Residues[spin][Mc_AN][h_AN][i]);
	}
	free(Residues[spin][Mc_AN][h_AN]);
      }
      free(Residues[spin][Mc_AN]);
    }
    free(Residues[spin]);
  }
  free(Residues);

  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
    free(PDOS_DC[Mc_AN]);
  }
  free(PDOS_DC);

  /* for time */
  dtime(&TEtime);
  time0 = TEtime - TStime;

  /* for PrintMemory */
  firsttime=0;

  return time0;
}

