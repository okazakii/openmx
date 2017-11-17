/**********************************************************************
  RestartFileDFT.c:

     RestartFileDFT.c is a subroutine to make a restart file
     as filename.rst.

  Log of RestartFileDFT.c:

     22/Nov/2001  Released by T.Ozaki 

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
/*  stat section */
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
/*  end stat section */
#include "openmx_common.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif




static int  Input_HKS( int MD_iter, double *Uele, double *****CH );
static void Output_HKS(int MD_iter, double *Uele, double *****CH );
static void Output_Charge_Density(int MD_iter);
static void Input_Charge_Density(int MD_iter, double *extpln_coes);
static void Inverse(int n, double **a, double **ia);
static void Extp_Charge(int MD_iter, double *extpln_coes);


int RestartFileDFT(char *mode, int MD_iter, double *Uele, double *****H, double *****CntH)
{
  int ret;
  int numprocs,myid;
  double *extpln_coes;

  /* MPI */
  if (atomnum<=MYID_MPI_COMM_WORLD) return 1;
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (ML_flag) return 0;

  ret = 1;

  /* write */
  if ( strcasecmp(mode,"write")==0 ) {

    MPI_Barrier(mpi_comm_level1);

    if      (Cnt_switch==0)  Output_HKS(MD_iter, Uele, H);
    else if (Cnt_switch==1)  Output_HKS(MD_iter, Uele, CntH);

    MPI_Barrier(mpi_comm_level1);

    Output_Charge_Density(MD_iter);
  }

  /* read */
  else if (strcasecmp(mode,"read")==0) {
    
    /* find coefficients for extrapolation of charge density or H */

    /* allocation */

    extpln_coes = (double*)malloc(sizeof(double)*(Extrapolated_Charge_History+2)); 

    Extp_Charge(MD_iter,extpln_coes);

    /* read H */

    if      (Cnt_switch==0)  ret = Input_HKS(MD_iter, Uele, H);
    else if (Cnt_switch==1)  ret = Input_HKS(MD_iter, Uele, CntH);

    /* read charge density */

    Input_Charge_Density(MD_iter,extpln_coes);

    /* free */

    free(extpln_coes);
  }

  return ret;
}

  

 
void Extp_Charge(int MD_iter, double *extpln_coes)
{
  int i,j,k;
  int NumHis;
  double sum;
  double **A,**IA,*B;

  for (i=0; i<Extrapolated_Charge_History; i++){
    extpln_coes[i] = 0.0;
  }
  extpln_coes[0] = 1.0;

  if (Cnt_switch==1) return;

  if (Extrapolated_Charge_History<MD_iter && 1<Extrapolated_Charge_History){

    /* allocation */

    A = (double**)malloc(sizeof(double*)*(Extrapolated_Charge_History+2));
    for (i=0; i<(Extrapolated_Charge_History+2); i++){
      A[i] = (double*)malloc(sizeof(double)*(Extrapolated_Charge_History+2));
    }

    IA = (double**)malloc(sizeof(double*)*(Extrapolated_Charge_History+2));
    for (i=0; i<(Extrapolated_Charge_History+2); i++){
      IA[i] = (double*)malloc(sizeof(double)*(Extrapolated_Charge_History+2));
    }
 
    B = (double*)malloc(sizeof(double)*(Extrapolated_Charge_History+2));

    /* find the current number of history */

    NumHis = MD_iter - 1;
    
    if (Extrapolated_Charge_History<NumHis){
      NumHis = Extrapolated_Charge_History;
    }

    /* make the matrix A */

    for (i=0; i<NumHis; i++){
      for (j=i; j<NumHis; j++){

        sum = rnd(1.0e-16);
        for (k=0; k<3*atomnum; k++){
          sum += His_Gxyz[i][k]*His_Gxyz[j][k];       
	}

        A[i][j] = sum;
        A[j][i] = sum;
      }
    }

    /* make the vector B */

    for (i=0; i<NumHis; i++){

      sum = 0.0;
      k = 0;

      for (j=1; j<=atomnum; j++){
	sum += His_Gxyz[i][k]*Gxyz[j][1]; k++;       
	sum += His_Gxyz[i][k]*Gxyz[j][2]; k++;       
	sum += His_Gxyz[i][k]*Gxyz[j][3]; k++;       
      }

      B[i] = sum; 
    }
  
    /* calculate the inverse of A */

    Inverse(NumHis-1,A,IA);

    /* calculate the coefficients */

    for (i=0; i<NumHis; i++){

      sum = 0.0;
      for (j=0; j<NumHis; j++){
        sum += IA[i][j]*B[j]; 
      }

      extpln_coes[i] = sum;
    }

    /* free */

    for (i=0; i<(Extrapolated_Charge_History+2); i++){
      free(A[i]);
    }
    free(A);
 
    for (i=0; i<(Extrapolated_Charge_History+2); i++){
      free(IA[i]);
    }
    free(IA);

    free(B);

  } /* if (2<MD_iter)  */   

  /* shift His_Gxyz */

  for (i=(Extrapolated_Charge_History-2); 0<=i; i--){
    for (k=0; k<3*atomnum; k++){
      His_Gxyz[i+1][k] = His_Gxyz[i][k];
    }
  }

  k = 0;
  for (i=1; i<=atomnum; i++){
    His_Gxyz[0][k] = Gxyz[i][1]; k++;       
    His_Gxyz[0][k] = Gxyz[i][2]; k++;       
    His_Gxyz[0][k] = Gxyz[i][3]; k++;       
  }
}



void Inverse(int n, double **a, double **ia)
{
  int i,j,k;
  double sum;
  double **a0,*ko,*iko;

  /***************************************************
    allocation of arrays: 
  ***************************************************/

  a0 = (double**)malloc(sizeof(double*)*(Extrapolated_Charge_History+3));
  for (i=0; i<(Extrapolated_Charge_History+3); i++){
    a0[i] = (double*)malloc(sizeof(double)*(Extrapolated_Charge_History+3));
  }

  ko = (double*)malloc(sizeof(double)*(Extrapolated_Charge_History+3));
  iko = (double*)malloc(sizeof(double)*(Extrapolated_Charge_History+3));

  /***************************************************
    calculate the inverse
  ***************************************************/

  for (i=0; i<=n; i++){
    for (j=0; j<=n; j++){
      a0[i+1][j+1] = a[i][j];
    }  
  }  

  Eigen_lapack(a0,ko,n+1,n+1);

  for (i=1; i<=(n+1); i++){
    if (fabs(ko[i])<1.0e-12) 
      iko[i] = 0.0;
    else    
      iko[i] = 1.0/ko[i];
  }

  for (i=1; i<=(n+1); i++){
    for (j=1; j<=(n+1); j++){

      sum = 0.0;

      for (k=1; k<=(n+1); k++){
        sum += a0[i][k]*iko[k]*a0[j][k]; 
      }
      ia[i-1][j-1] = sum;
    }
  }

  /***************************************************
    freeing of arrays: 
  ***************************************************/

  for (i=0; i<(Extrapolated_Charge_History+3); i++){
    free(a0[i]);
  }
  free(a0);

  free(ko);
  free(iko);
}




int Input_HKS( int MD_iter, double *Uele, double *****CH )
{
  int Mc_AN,Gc_AN,h_AN,i,j,can;
  int Gh_AN,pSCF,spin,Rn;
  int wan1,wan2,TNO1,TNO2;
  int h_AN0,Gh_AN0,Rn0,wan20,TNO20;
  int i_vec[20],*p_vec,po;
  int my_check,exit_flag;
  int pFNAN;
  int l1,l2,l3,l10,l20,l30;
  int numprocs,myid;
  char fileHKS[YOUSO10];
  FILE *fp;
  double *tmpvec;  
  char buf[fp_bsize];          /* setvbuf */

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

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

    if ((fp = fopen(fileHKS,"rb")) != NULL){

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
              read Gh_AN0, l10, l20, l30, wan20, TNO20
        ****************************************************/

        p_vec = (int*)malloc(sizeof(int)*(pFNAN+1)*6);
        fread(p_vec, sizeof(int), (pFNAN+1)*6, fp);
        fread(Uele,sizeof(double),1,fp);

        /****************************************************
          store Hamiltonian to appropriate position while 
          comparing Gh_AN, l1, l2, l3, wan2, TNO2
        ****************************************************/

        for (spin=0; spin<=SpinP_switch; spin++){ 
	  for (h_AN0=0; h_AN0<=pFNAN; h_AN0++){
	    Gh_AN0 = p_vec[              h_AN0];
	    l10    = p_vec[(pFNAN+1)*1 + h_AN0];
	    l20    = p_vec[(pFNAN+1)*2 + h_AN0];
	    l30    = p_vec[(pFNAN+1)*3 + h_AN0];
	    wan20  = p_vec[(pFNAN+1)*4 + h_AN0];
	    TNO20  = p_vec[(pFNAN+1)*5 + h_AN0]; 
 
	    exit_flag = 0;
	    h_AN = 0;

	    do {

	      Gh_AN = natn[Gc_AN][h_AN];
	      Rn    = ncn[Gc_AN][h_AN];
              l1 = atv_ijk[Rn][1];
              l2 = atv_ijk[Rn][2];
              l3 = atv_ijk[Rn][3];
	      wan2  = WhatSpecies[Gh_AN];
	      TNO2  = Spe_Total_CNO[wan2];

	      if ( Gh_AN==Gh_AN0 &&
		   l1==l10 &&
		   l2==l20 &&
		   l3==l30 &&
		   wan2==wan20   &&
		   TNO2==TNO20 )
		{

		  for (i=0; i<TNO1; i++){
		    fread(&CH[spin][Mc_AN][h_AN][i][0],sizeof(double),TNO2,fp);
		  }

		  /* add contribution of SO coupling */
		  if (i_vec[9]==0 && SO_switch==1 && spin==2){

		    for (i=0; i<TNO1; i++){
		      for (j=0; j<TNO2; j++){
			CH[2][Mc_AN][h_AN][i][j] += HNL[2][Mc_AN][h_AN][i][j];
		      }
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

	  } /* h_AN */
	} /* spin */

        /****************************************************
          store iHNL to appropriate position while 
          comparing Gh_AN, Rn, wan2, TNO2
        ****************************************************/

        if (SpinP_switch==3){

	  for (spin=0; spin<SpinP_switch; spin++){
	    for (h_AN0=0; h_AN0<=pFNAN; h_AN0++){

	      Gh_AN0 = p_vec[              h_AN0];
	      l10    = p_vec[(pFNAN+1)*1 + h_AN0];
	      l20    = p_vec[(pFNAN+1)*2 + h_AN0];
	      l30    = p_vec[(pFNAN+1)*3 + h_AN0];
	      wan20  = p_vec[(pFNAN+1)*4 + h_AN0];
	      TNO20  = p_vec[(pFNAN+1)*5 + h_AN0]; 
 
	      exit_flag = 0;
	      h_AN = 0;

	      do {

		Gh_AN = natn[Gc_AN][h_AN];
		Rn    = ncn[Gc_AN][h_AN];
                l1 = atv_ijk[Rn][1];
                l2 = atv_ijk[Rn][2];
                l3 = atv_ijk[Rn][3];
		wan2  = WhatSpecies[Gh_AN];
		TNO2  = Spe_Total_CNO[wan2];

		if ( Gh_AN==Gh_AN0 &&
		     l1==l10 &&
		     l2==l20 &&
		     l3==l30 &&
		     wan2==wan20   &&
		     TNO2==TNO20 )
		  {

		    for (i=0; i<TNO1; i++){
		      fread(&iHNL[spin][Mc_AN][h_AN][i][0],sizeof(double),TNO2,fp);
		    }

                    /* add contribution of SO coupling */
                    if (i_vec[9]==0 && SO_switch==1){

  	  	      for (i=0; i<TNO1; i++){
    	    	        for (j=0; j<TNO2; j++){
                          iHNL[spin][Mc_AN][h_AN][i][j] += iHNL0[spin][Mc_AN][h_AN][i][j];
			}
		      }
                    }

		    exit_flag = 1;
		  }

		h_AN++;

	      } while (h_AN<=FNAN[Gc_AN] && exit_flag==0);

	      /* In case appropriate one is not found, just read the iHNL */ 

	      if (exit_flag==0){
		for (i=0; i<TNO1; i++){
		  fread(&tmpvec[0],sizeof(double),TNO20,fp);
		}
	      }

	    } /* h_AN */
	  } /* spin */

	} /* if (SpinP_switch==3) */

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

  return my_check;
}





void Output_HKS(int MD_iter, double *Uele, double *****CH )
{
  int Mc_AN,Gc_AN,h_AN,i,j,can,Gh_AN;
  int wan1,wan2,TNO1,TNO2,spin,Rn,po;
  int i_vec[20],*p_vec;
  int numprocs,myid;
  char operate[300];
  char fileHKS[YOUSO10];
  FILE *fp;
  char buf[fp_bsize];          /* setvbuf */

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /* delete if there are */

  if (myid==Host_ID){

    if (MD_iter==1){
      sprintf(operate,"%s%s_rst",filepath,filename);
      mkdir(operate,0775); 
    }

    for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
      sprintf(operate,"%s%s_rst/%s.rst%i",filepath,filename,filename,Gc_AN);
      remove(operate);
    }
  }

  MPI_Barrier(mpi_comm_level1);

  /* write */

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];
    wan1 = WhatSpecies[Gc_AN];
    TNO1 = Spe_Total_CNO[wan1];

    sprintf(fileHKS,"%s%s_rst/%s.rst%i",filepath,filename,filename,Gc_AN);

    po = 0;

    do {

      if ((fp = fopen(fileHKS,"wb")) != NULL){

#ifdef xt3
         setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

	/****************************************************
         List_YOUSO[23] 0:  non spin poralized
                        1:  spin poralized
                        3:  spin non-collinear
         List_YOUSO[1]  atomnum
         List_YOUSO[8]  max # of atoms in a rcut-off cluster
         List_YOUSO[7]  max # of orbitals including an atom
	****************************************************/

	i_vec[0] = SpinP_switch;
	i_vec[1] = List_YOUSO[23];
	i_vec[2] = List_YOUSO[1];
	i_vec[3] = List_YOUSO[8];
	i_vec[4] = List_YOUSO[7];
	i_vec[5] = atomnum;
	i_vec[6] = wan1;
	i_vec[7] = TNO1;
	i_vec[8] = FNAN[Gc_AN];
	i_vec[9] = SO_switch;

	fwrite(i_vec,sizeof(int),10,fp);

	/****************************************************
                  # of orbitals in each FNAN atom
	****************************************************/

	p_vec = (int*)malloc(sizeof(int)*(FNAN[Gc_AN]+1)*6);
	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	  Gh_AN = natn[Gc_AN][h_AN];
	  Rn = ncn[Gc_AN][h_AN];
	  wan2 = WhatSpecies[Gh_AN];
	  TNO2 = Spe_Total_CNO[wan2];
	  p_vec[                    h_AN] = Gh_AN;
	  p_vec[(FNAN[Gc_AN]+1)*1 + h_AN] = atv_ijk[Rn][1];
	  p_vec[(FNAN[Gc_AN]+1)*2 + h_AN] = atv_ijk[Rn][2];
	  p_vec[(FNAN[Gc_AN]+1)*3 + h_AN] = atv_ijk[Rn][3];
	  p_vec[(FNAN[Gc_AN]+1)*4 + h_AN] = wan2;
	  p_vec[(FNAN[Gc_AN]+1)*5 + h_AN] = TNO2;
	}
	fwrite(p_vec,sizeof(int), (FNAN[Gc_AN]+1)*6, fp);
	free(p_vec);

	/****************************************************
                           Uele
	****************************************************/

	fwrite(Uele,sizeof(double),1,fp);

	/****************************************************
                        Kohn-Sham Hamiltonian
	****************************************************/

	for (spin=0; spin<=SpinP_switch; spin++){
	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	    Gh_AN = natn[Gc_AN][h_AN];
	    wan2 = WhatSpecies[Gh_AN];
	    TNO2 = Spe_Total_CNO[wan2];
	    for (i=0; i<TNO1; i++){
	      fwrite(CH[spin][Mc_AN][h_AN][i],sizeof(double),TNO2,fp);
	    }
	  }
	}

        if (SpinP_switch==3){
	  for (spin=0; spin<SpinP_switch; spin++){
	    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	      Gh_AN = natn[Gc_AN][h_AN];
	      wan2 = WhatSpecies[Gh_AN];
	      TNO2 = Spe_Total_CNO[wan2];
	      for (i=0; i<TNO1; i++){
	        fwrite(iHNL[spin][Mc_AN][h_AN][i],sizeof(double),TNO2,fp);
	      }
	    }
  	  }
        }

	fclose(fp);

	po = 1;

      }
      else{
	printf("Failed in saving the restart file %s\n",fileHKS);      
        /* openmx gives up to save the restart file. */
        po = 1;
      }

    } while (po==0);
  }

}





void Output_Charge_Density(int MD_iter)
{
  int i,j,k,i1,i2,i3,c;
  int GN,mul,n1,n2,n3,nn1,nn0;
  int cmd,MN,MN0,MN1,MN2,MN3;
  double ***V;
  double *tmp_array0;
  double *tmp_array1;
  int numprocs,myid,ID,IDS,IDR,tag=999;
  char operate[300];
  char fname1[300];
  char fname2[300];
  FILE *fp;
  char fileCD0[YOUSO10];
  char fileCD1[YOUSO10];
  char fileCD2[YOUSO10];
  char buf[fp_bsize];          /* setvbuf */

  MPI_Status stat;
  MPI_Request request;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if      (SpinP_switch==0) mul = 1;
  else if (SpinP_switch==1) mul = 2;
  else if (SpinP_switch==3) mul = 4;
  
  /****************************************************
   allocation of arrays:

   double V[mul][My_NGrid1_Poisson][Ngrid2*Ngrid3];
  ****************************************************/

  V = (double***)malloc(sizeof(double**)*mul); 
  for (k=0; k<mul; k++){
    V[k] = (double**)malloc(sizeof(double*)*My_NGrid1_Poisson); 
    for (i=0; i<My_NGrid1_Poisson; i++){
      V[k][i] = (double*)malloc(sizeof(double)*Ngrid2*Ngrid3); 
    }
  }

  /****************************************************
                    set V0 and V1
  ****************************************************/

  /* initialize */

  for (k=0; k<mul; k++){
    for (n1=0; n1<My_NGrid1_Poisson; n1++){
      for (n2=0; n2<Ngrid2*Ngrid3; n2++){
        V[k][n1][n2] = 0.0;
      }
    }
  }

  /* use their densities using MPI */

  for (k=0; k<mul; k++){

    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      /* Isend */
      if (Num_Snd_Grid1[IDS]!=0){

        tmp_array0 = (double*)malloc(sizeof(double)*Num_Snd_Grid1[IDS]*Ngrid2*Ngrid3); 
  
        for (i=0; i<Num_Snd_Grid1[IDS]; i++){ 

	  n1 = Snd_Grid1[IDS][i];
          nn1 = My_Cell0[n1];
          MN1 = nn1*Ngrid2*Ngrid3;
          MN0 = i*Ngrid2*Ngrid3;
          for (n2=0; n2<Ngrid2; n2++){
            MN2 = n2*Ngrid3;

            if (k<=1){
	      for (n3=0; n3<Ngrid3; n3++){
		MN = MN1 + MN2 + n3;
		tmp_array0[MN0+MN2+n3] = Density_Grid[k][MN] - ADensity_Grid[MN];
	      }
	    }

            else {
	      for (n3=0; n3<Ngrid3; n3++){
		MN = MN1 + MN2 + n3;
		tmp_array0[MN0+MN2+n3] = Density_Grid[k][MN];
	      }
            }

	  }
        }

        MPI_Isend(&tmp_array0[0], Num_Snd_Grid1[IDS]*Ngrid2*Ngrid3,
                  MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);
      }

      /* Recv */
      if (Num_Rcv_Grid1[IDR]!=0){

        tmp_array1 = (double*)malloc(sizeof(double)*Num_Rcv_Grid1[IDR]*Ngrid2*Ngrid3); 

        MPI_Recv(&tmp_array1[0], Num_Rcv_Grid1[IDR]*Ngrid2*Ngrid3,
                MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

        for (i=0; i<Num_Rcv_Grid1[IDR]; i++){ 
	  n1 = Rcv_Grid1[IDR][i];
          nn1 = My_Cell0[n1];
          nn0 = n1 - Start_Grid1[myid];
          MN0 = i*Ngrid2*Ngrid3;
          for (n2=0; n2<Ngrid2; n2++){
            MN2 = n2*Ngrid3;
	    for (n3=0; n3<Ngrid3; n3++){
	      MN = MN0 + MN2 + n3;
	      V[k][nn0][MN2+n3] = tmp_array1[MN];
	    }
	  }
        }

        free(tmp_array1);
      }

      if (Num_Snd_Grid1[IDS]!=0){
        MPI_Wait(&request,&stat);
        free(tmp_array0);
      }

    }

    /* use own densities */
    for (n1=Start_Grid1[myid]; n1<=End_Grid1[myid]; n1++){
      nn1 = My_Cell0[n1];
      nn0 = n1 - Start_Grid1[myid]; 
      if (nn1!=-1){
        MN1 = nn1*Ngrid2*Ngrid3;
        for (n2=0; n2<Ngrid2; n2++){
          MN2 = n2*Ngrid3;

          if (k<=1){
	    for (n3=0; n3<Ngrid3; n3++){
	      MN = MN1 + MN2 + n3;
	      V[k][nn0][MN2+n3] = Density_Grid[k][MN] - ADensity_Grid[MN];
	    }    
	  }
          else {
	    for (n3=0; n3<Ngrid3; n3++){
	      MN = MN1 + MN2 + n3;
	      V[k][nn0][MN2+n3] = Density_Grid[k][MN];
	    }    
	  }
        }    
      }
    }

  } /* mul */

  /****************************************************
                    output data 
  ****************************************************/


  for (k=0; k<mul; k++){
    for (n1=0; n1<My_NGrid1_Poisson; n1++){

      nn0 = n1 + Start_Grid1[myid]; 

      if (Cnt_switch!=1){ 

	for (i=(Extrapolated_Charge_History-2); 0<=i; i--){
 
	  sprintf(fileCD1,"%s%s_rst/%s.crst%i_%i_%i",filepath,filename,filename,k,nn0,i);
	  sprintf(fileCD2,"%s%s_rst/%s.crst%i_%i_%i",filepath,filename,filename,k,nn0,i+1);

	  if ((fp = fopen(fileCD1,"rb")) != NULL){
	    fclose(fp);
	    rename(fileCD1,fileCD2); 
	  }
	} 
      }

      sprintf(fileCD0,"%s%s_rst/%s.crst%i_%i_0",filepath,filename,filename,k,nn0);

      if ((fp = fopen(fileCD0,"wb")) != NULL){
        fwrite(V[k][n1],sizeof(double),Ngrid2*Ngrid3,fp);
	fclose(fp);
      }
      else{
        printf("Could not open a file %s\n",fileCD0);
      }
    }
  }

  /****************************************************
   freeing of arrays:

   double V[mul][My_NGrid1_Poisson][Ngrid2*Ngrid3];
  ****************************************************/

  for (k=0; k<mul; k++){
    for (i=0; i<My_NGrid1_Poisson; i++){
      free(V[k][i]);
    }
    free(V[k]);
  }
  free(V);

}








void Input_Charge_Density(int MD_iter, double *extpln_coes)
{
  int i,j,k,i1,i2,i3,c,spin,ri;
  int GN,mul,n1,n2,n3,nn1,nn0;
  int cmd,MN,MN0,MN1,MN2,MN3;
  int numprocs,myid,ID,IDS,IDR,tag=999;
  char operate[300];
  char fname1[300];
  char fname2[300];
  FILE *fp;
  char fileCD[YOUSO10];
  char buf[fp_bsize];          /* setvbuf */
  double *tmp_array;

  MPI_Status stat;
  MPI_Request request;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /* allocation */
  tmp_array = (double*)malloc(sizeof(double)*Ngrid2*Ngrid3);

  for (spin=0; spin<=SpinP_switch; spin++){

    for (i=0; i<Num_Cells0; i++){
      ri = My_Cell1[i];
      MN = i*Ngrid2*Ngrid3;

      for (k=0; k<Extrapolated_Charge_History; k++){

        sprintf(fileCD,"%s%s_rst/%s.crst%i_%i_%i",filepath,filename,filename,spin,ri,k);

	if ((fp = fopen(fileCD,"rb")) != NULL){

	  fread(tmp_array,sizeof(double),Ngrid2*Ngrid3,fp);
	  fclose(fp);

          if (k==0){

            if (spin<=1){
	      for (MN1=0; MN1<Ngrid2*Ngrid3; MN1++){
		Density_Grid[spin][MN+MN1] = ADensity_Grid[MN+MN1] + extpln_coes[k]*tmp_array[MN1];
	      }
	    }
            else{
	      for (MN1=0; MN1<Ngrid2*Ngrid3; MN1++){
		Density_Grid[spin][MN+MN1] = extpln_coes[k]*tmp_array[MN1];
	      }
	    }

          }

          else{
            for (MN1=0; MN1<Ngrid2*Ngrid3; MN1++){
              Density_Grid[spin][MN+MN1] += extpln_coes[k]*tmp_array[MN1];
	    }
          }

	}
	else{
	  /*
	  printf("Could not open a file %s\n",fileCD);
	  */
	}
      }

    }
  }

  if (SpinP_switch==0){
    for (MN=0; MN<My_NumGrid1; MN++){
      Density_Grid[1][MN] = Density_Grid[0][MN];
    }
  }

  /* free */
  free(tmp_array);

}
