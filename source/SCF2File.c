/**********************************************************************
  SCF2File.c:

     SCF2File.c is a subroutine to output connectivity, Hamiltonian,
     overlap, and etc. to a binary file, filename.scfout.

  Log of SCF2DFT.c:

     1/July/2003  Released by T.Ozaki 

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "openmx_common.h"
#include "mpi.h"

#define MAX_LINE_SIZE 256

static void Output( FILE *fp, char *inputfile );
static void Calc_OLPpo();

static double ****OLPpox;
static double ****OLPpoy;
static double ****OLPpoz;



void SCF2File(char *mode, char *inputfile)
{
  int Mc_AN,Gc_AN,tno0,tno1;
  int Cwan,h_AN,Gh_AN,Hwan;
  int i;
  char fname[YOUSO10];
  char operate[300];
  FILE *fp;
  int numprocs,myid,ID;
  int succeed_open;
  char buf[fp_bsize];          /* setvbuf */

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /* allocation of array */

  OLPpox = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
  FNAN[0] = 0;
  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

    if (Mc_AN==0){
      Gc_AN = 0;
      tno0 = 1;
    }
    else{
      Gc_AN = S_M2G[Mc_AN];
      Cwan = WhatSpecies[Gc_AN];
      tno0 = Spe_Total_CNO[Cwan];  
    }    

    OLPpox[Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

      if (Mc_AN==0){
	tno1 = 1;  
      }
      else{
	Gh_AN = natn[Gc_AN][h_AN];
	Hwan = WhatSpecies[Gh_AN];
	tno1 = Spe_Total_CNO[Hwan];
      } 

      OLPpox[Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
      for (i=0; i<tno0; i++){
	OLPpox[Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
      }
    }
  }

  OLPpoy = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
  FNAN[0] = 0;
  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

    if (Mc_AN==0){
      Gc_AN = 0;
      tno0 = 1;
    }
    else{
      Gc_AN = S_M2G[Mc_AN];
      Cwan = WhatSpecies[Gc_AN];
      tno0 = Spe_Total_CNO[Cwan];  
    }    

    OLPpoy[Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

      if (Mc_AN==0){
	tno1 = 1;  
      }
      else{
	Gh_AN = natn[Gc_AN][h_AN];
	Hwan = WhatSpecies[Gh_AN];
	tno1 = Spe_Total_CNO[Hwan];
      } 

      OLPpoy[Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
      for (i=0; i<tno0; i++){
	OLPpoy[Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
      }
    }
  }

  OLPpoz = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
  FNAN[0] = 0;
  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

    if (Mc_AN==0){
      Gc_AN = 0;
      tno0 = 1;
    }
    else{
      Gc_AN = S_M2G[Mc_AN];
      Cwan = WhatSpecies[Gc_AN];
      tno0 = Spe_Total_CNO[Cwan];  
    }    

    OLPpoz[Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

      if (Mc_AN==0){
	tno1 = 1;  
      }
      else{
	Gh_AN = natn[Gc_AN][h_AN];
	Hwan = WhatSpecies[Gh_AN];
	tno1 = Spe_Total_CNO[Hwan];
      } 

      OLPpoz[Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
      for (i=0; i<tno0; i++){
	OLPpoz[Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
      }
    }
  }

  /* write a file */

  if ( strcasecmp(mode,"write")==0 ) {
    sprintf(fname,"%s%s.scfout",filepath,filename);

    if (myid==Host_ID){

      if ((fp = fopen(fname,"w")) != NULL){

#ifdef xt3
        setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

        succeed_open = 1;
      }
      else {
        succeed_open = 0;
        printf("Failure of making the scfout file (%s).\n",fname);
      }
    }

    MPI_Bcast(&succeed_open, 1, MPI_INT, Host_ID, mpi_comm_level1);
 
    if (succeed_open==1){
      if (myid==Host_ID) printf("Save the scfout file (%s)\n",fname);

      Calc_OLPpo();
      Output( fp, inputfile );
      if (myid==Host_ID) fclose(fp);
    } 
  }

  /* free array */

  FNAN[0] = 0;
  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

    if (Mc_AN==0){
      Gc_AN = 0;
      tno0 = 1;
    }
    else{
      Gc_AN = S_M2G[Mc_AN];
      Cwan = WhatSpecies[Gc_AN];
      tno0 = Spe_Total_CNO[Cwan];  
    }    

    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

      if (Mc_AN==0){
	tno1 = 1;  
      }
      else{
	Gh_AN = natn[Gc_AN][h_AN];
	Hwan = WhatSpecies[Gh_AN];
	tno1 = Spe_Total_CNO[Hwan];
      } 

      for (i=0; i<tno0; i++){
	free(OLPpox[Mc_AN][h_AN][i]);
      }
      free(OLPpox[Mc_AN][h_AN]);
    }
    free(OLPpox[Mc_AN]);
  }
  free(OLPpox);

  FNAN[0] = 0;
  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

    if (Mc_AN==0){
      Gc_AN = 0;
      tno0 = 1;
    }
    else{
      Gc_AN = S_M2G[Mc_AN];
      Cwan = WhatSpecies[Gc_AN];
      tno0 = Spe_Total_CNO[Cwan];  
    }    

    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

      if (Mc_AN==0){
	tno1 = 1;  
      }
      else{
	Gh_AN = natn[Gc_AN][h_AN];
	Hwan = WhatSpecies[Gh_AN];
	tno1 = Spe_Total_CNO[Hwan];
      } 

      for (i=0; i<tno0; i++){
	free(OLPpoy[Mc_AN][h_AN][i]);
      }
      free(OLPpoy[Mc_AN][h_AN]);
    }
    free(OLPpoy[Mc_AN]);
  }
  free(OLPpoy);

  FNAN[0] = 0;
  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

    if (Mc_AN==0){
      Gc_AN = 0;
      tno0 = 1;
    }
    else{
      Gc_AN = S_M2G[Mc_AN];
      Cwan = WhatSpecies[Gc_AN];
      tno0 = Spe_Total_CNO[Cwan];  
    }    

    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

      if (Mc_AN==0){
	tno1 = 1;  
      }
      else{
	Gh_AN = natn[Gc_AN][h_AN];
	Hwan = WhatSpecies[Gh_AN];
	tno1 = Spe_Total_CNO[Hwan];
      } 

      for (i=0; i<tno0; i++){
	free(OLPpoz[Mc_AN][h_AN][i]);
      }
      free(OLPpoz[Mc_AN][h_AN]);
    }
    free(OLPpoz[Mc_AN]);
  }
  free(OLPpoz);
}




void Output( FILE *fp, char *inputfile )
{
  int Gc_AN,Mc_AN,ct_AN,h_AN,i,j,can,Gh_AN;
  int num,wan1,wan2,TNO1,TNO2,spin,Rn,count;
  int k,Mh_AN,Rnh,Hwan,NO0,NO1,Qwan;
  int q_AN,Mq_AN,Gq_AN,Rnq;
  int i_vec[20],*p_vec;
  int numprocs,myid,ID,tag=999;
  double *Tmp_Vec,d_vec[20];
  FILE *fp_inp;
  char strg[MAX_LINE_SIZE];
  char buf[fp_bsize];          /* setvbuf */

  MPI_Status stat;
  MPI_Request request;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /***************************************************************
                     information of connectivity
  ****************************************************************/

  if (myid==Host_ID){

    /****************************************************
       atomnum
       spinP_switch 
    ****************************************************/

    i = 0;
    i_vec[0] = atomnum;
    i_vec[1] = SpinP_switch;
    i_vec[2] = Catomnum;
    i_vec[3] = Latomnum;
    i_vec[4] = Ratomnum;
    i_vec[5] = TCpyCell;

    fwrite(i_vec,sizeof(int),6,fp);

    /****************************************************
                       atv[Rn][4]
    ****************************************************/

    for (Rn=0; Rn<=TCpyCell; Rn++){
      fwrite(atv[Rn],sizeof(double),4,fp);
    }

    /****************************************************
                       atv_ijk[Rn][4]
    ****************************************************/

    for (Rn=0; Rn<=TCpyCell; Rn++){
      fwrite(atv_ijk[Rn],sizeof(int),4,fp);
    }

    /****************************************************
                 # of orbitals in each atom
    ****************************************************/

    p_vec = (int*)malloc(sizeof(int)*atomnum);
    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      wan1 = WhatSpecies[ct_AN];
      TNO1 = Spe_Total_CNO[wan1];
      p_vec[ct_AN-1] = TNO1;
    }
    fwrite(p_vec,sizeof(int),atomnum,fp);
    free(p_vec);

    /****************************************************
     FNAN[]:
     number of first nearest neighbouring atoms
    ****************************************************/

    p_vec = (int*)malloc(sizeof(int)*atomnum);
    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      p_vec[ct_AN-1] = FNAN[ct_AN];
    }
    fwrite(p_vec,sizeof(int),atomnum,fp);
    free(p_vec);

    /****************************************************
      natn[][]:
      grobal index of neighboring atoms of an atom ct_AN
     ****************************************************/
  
    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      fwrite(natn[ct_AN],sizeof(int),FNAN[ct_AN]+1,fp);
    }  
  
    /****************************************************
      ncn[][]:
      grobal index for cell of neighboring atoms
      of an atom ct_AN
    ****************************************************/
  
    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      fwrite(ncn[ct_AN],sizeof(int),FNAN[ct_AN]+1,fp);
    }  

    /****************************************************
      tv[4][4]:
      unit cell vectors in Bohr
    ****************************************************/

    fwrite(tv[1],sizeof(double),4,fp);
    fwrite(tv[2],sizeof(double),4,fp);
    fwrite(tv[3],sizeof(double),4,fp);

    /****************************************************
      rtv[4][4]:
      reciprocal unit cell vectors in Bohr^{-1}
    ****************************************************/

    fwrite(rtv[1],sizeof(double),4,fp);
    fwrite(rtv[2],sizeof(double),4,fp);
    fwrite(rtv[3],sizeof(double),4,fp);

    /****************************************************
      Gxyz[][1-3]:
      atomic coordinates in Bohr
    ****************************************************/

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      fwrite(Gxyz[ct_AN],sizeof(double),4,fp);
    }  

  } /* if (myid==Host_ID) */

  /***************************************************************
                         Hamiltonian matrix 
  ****************************************************************/

  Tmp_Vec = (double*)malloc(sizeof(double)*List_YOUSO[8]*List_YOUSO[7]*List_YOUSO[7]);

  for (spin=0; spin<=SpinP_switch; spin++){

    for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
      ID = G2ID[Gc_AN];

      if (myid==ID){

        num = 0;

        Mc_AN = F_G2M[Gc_AN];
        wan1 = WhatSpecies[Gc_AN];
        TNO1 = Spe_Total_CNO[wan1];
        for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
          Gh_AN = natn[Gc_AN][h_AN];
          wan2 = WhatSpecies[Gh_AN];
          TNO2 = Spe_Total_CNO[wan2];

          if (Cnt_switch==0){
            for (i=0; i<TNO1; i++){
              for (j=0; j<TNO2; j++){
                Tmp_Vec[num] = H[spin][Mc_AN][h_AN][i][j];
                num++;
	      }
            }
          }
          else{
            for (i=0; i<TNO1; i++){
              for (j=0; j<TNO2; j++){
                Tmp_Vec[num] = CntH[spin][Mc_AN][h_AN][i][j];
                num++;
	      }
            }
          }
        }

        if (myid!=Host_ID){
          MPI_Isend(&num, 1, MPI_INT, Host_ID, tag, mpi_comm_level1, &request);
          MPI_Wait(&request,&stat);
          MPI_Isend(&Tmp_Vec[0], num, MPI_DOUBLE, Host_ID, tag, mpi_comm_level1, &request);
          MPI_Wait(&request,&stat);
	}
        else{
          fwrite(Tmp_Vec, sizeof(double), num, fp);
        }
      }

      else if (ID!=myid && myid==Host_ID){
        MPI_Recv(&num, 1, MPI_INT, ID, tag, mpi_comm_level1, &stat);
        MPI_Recv(&Tmp_Vec[0], num, MPI_DOUBLE, ID, tag, mpi_comm_level1, &stat);
        fwrite(Tmp_Vec, sizeof(double), num, fp);
      }

    }  
  }

  /***************************************************************
                                  iHNL
  ****************************************************************/

  if ( SpinP_switch==3 ){

    for (spin=0; spin<List_YOUSO[5]; spin++){

      for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
	ID = G2ID[Gc_AN];

	if (myid==ID){

	  num = 0;

	  Mc_AN = F_G2M[Gc_AN];
	  wan1 = WhatSpecies[Gc_AN];
	  TNO1 = Spe_Total_CNO[wan1];
	  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	    Gh_AN = natn[Gc_AN][h_AN];
	    wan2 = WhatSpecies[Gh_AN];
	    TNO2 = Spe_Total_CNO[wan2];

	    for (i=0; i<TNO1; i++){
	      for (j=0; j<TNO2; j++){
		Tmp_Vec[num] = iHNL[spin][Mc_AN][h_AN][i][j];
		num++;
	      }
	    }
	  }

	  if (myid!=Host_ID){
	    MPI_Isend(&num, 1, MPI_INT, Host_ID, tag, mpi_comm_level1, &request);
	    MPI_Wait(&request,&stat);
	    MPI_Isend(&Tmp_Vec[0], num, MPI_DOUBLE, Host_ID, tag, mpi_comm_level1, &request);
	    MPI_Wait(&request,&stat);
	  }
	  else{
	    fwrite(Tmp_Vec, sizeof(double), num, fp);
	  }
	}

	else if (ID!=myid && myid==Host_ID){
	  MPI_Recv(&num, 1, MPI_INT, ID, tag, mpi_comm_level1, &stat);
	  MPI_Recv(&Tmp_Vec[0], num, MPI_DOUBLE, ID, tag, mpi_comm_level1, &stat);
	  fwrite(Tmp_Vec, sizeof(double), num, fp);
	}

      }  
    }
  }

  /***************************************************************
                          overlap matrix 
  ****************************************************************/

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    ID = G2ID[Gc_AN];

    if (myid==ID){

      num = 0;

      Mc_AN = F_G2M[Gc_AN];
      wan1 = WhatSpecies[Gc_AN];
      TNO1 = Spe_Total_CNO[wan1];
      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
        Gh_AN = natn[Gc_AN][h_AN];
        wan2 = WhatSpecies[Gh_AN];
        TNO2 = Spe_Total_CNO[wan2];

        if (Cnt_switch==0){
          for (i=0; i<TNO1; i++){
            for (j=0; j<TNO2; j++){
              Tmp_Vec[num] = OLP[0][Mc_AN][h_AN][i][j];
              num++;
            }
          }
        }
        else{
          for (i=0; i<TNO1; i++){
            for (j=0; j<TNO2; j++){
              Tmp_Vec[num] = CntOLP[0][Mc_AN][h_AN][i][j];
              num++;
	    }
          }
        }
      }

      if (myid!=Host_ID){
        MPI_Isend(&num, 1, MPI_INT, Host_ID, tag, mpi_comm_level1, &request);
        MPI_Wait(&request,&stat);
        MPI_Isend(&Tmp_Vec[0], num, MPI_DOUBLE, Host_ID, tag, mpi_comm_level1, &request);
        MPI_Wait(&request,&stat);
      }
      else{
        fwrite(Tmp_Vec, sizeof(double), num, fp);
      }
    }

    else if (ID!=myid && myid==Host_ID){
      MPI_Recv(&num, 1, MPI_INT, ID, tag, mpi_comm_level1, &stat);
      MPI_Recv(&Tmp_Vec[0], num, MPI_DOUBLE, ID, tag, mpi_comm_level1, &stat);
      fwrite(Tmp_Vec, sizeof(double), num, fp);
    }

  }  

  /***************************************************************
                overlap matrix with position operator x
  ****************************************************************/

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    ID = G2ID[Gc_AN];

    if (myid==ID){

      num = 0;
      Mc_AN = F_G2M[Gc_AN];
      wan1 = WhatSpecies[Gc_AN];
      TNO1 = Spe_Total_CNO[wan1];

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
        Gh_AN = natn[Gc_AN][h_AN];
        wan2 = WhatSpecies[Gh_AN];
        TNO2 = Spe_Total_CNO[wan2];

	for (i=0; i<TNO1; i++){
	  for (j=0; j<TNO2; j++){
	    Tmp_Vec[num] = OLPpox[Mc_AN][h_AN][i][j];
	    num++;
	  }
	}
      }

      if (myid!=Host_ID){
        MPI_Isend(&num, 1, MPI_INT, Host_ID, tag, mpi_comm_level1, &request);
        MPI_Wait(&request,&stat);
        MPI_Isend(&Tmp_Vec[0], num, MPI_DOUBLE, Host_ID, tag, mpi_comm_level1, &request);
        MPI_Wait(&request,&stat);
      }
      else{
        fwrite(Tmp_Vec, sizeof(double), num, fp);
      }
    }

    else if (ID!=myid && myid==Host_ID){
      MPI_Recv(&num, 1, MPI_INT, ID, tag, mpi_comm_level1, &stat);
      MPI_Recv(&Tmp_Vec[0], num, MPI_DOUBLE, ID, tag, mpi_comm_level1, &stat);
      fwrite(Tmp_Vec, sizeof(double), num, fp);
    }
  }  

  /***************************************************************
                overlap matrix with position operator y
  ****************************************************************/

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    ID = G2ID[Gc_AN];

    if (myid==ID){

      num = 0;
      Mc_AN = F_G2M[Gc_AN];
      wan1 = WhatSpecies[Gc_AN];
      TNO1 = Spe_Total_CNO[wan1];

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
        Gh_AN = natn[Gc_AN][h_AN];
        wan2 = WhatSpecies[Gh_AN];
        TNO2 = Spe_Total_CNO[wan2];

	for (i=0; i<TNO1; i++){
	  for (j=0; j<TNO2; j++){
	    Tmp_Vec[num] = OLPpoy[Mc_AN][h_AN][i][j];
	    num++;
	  }
	}
      }

      if (myid!=Host_ID){
        MPI_Isend(&num, 1, MPI_INT, Host_ID, tag, mpi_comm_level1, &request);
        MPI_Wait(&request,&stat);
        MPI_Isend(&Tmp_Vec[0], num, MPI_DOUBLE, Host_ID, tag, mpi_comm_level1, &request);
        MPI_Wait(&request,&stat);
      }
      else{
        fwrite(Tmp_Vec, sizeof(double), num, fp);
      }
    }

    else if (ID!=myid && myid==Host_ID){
      MPI_Recv(&num, 1, MPI_INT, ID, tag, mpi_comm_level1, &stat);
      MPI_Recv(&Tmp_Vec[0], num, MPI_DOUBLE, ID, tag, mpi_comm_level1, &stat);
      fwrite(Tmp_Vec, sizeof(double), num, fp);
    }
  }  

  /***************************************************************
                overlap matrix with position operator z
  ****************************************************************/

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    ID = G2ID[Gc_AN];

    if (myid==ID){

      num = 0;
      Mc_AN = F_G2M[Gc_AN];
      wan1 = WhatSpecies[Gc_AN];
      TNO1 = Spe_Total_CNO[wan1];

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
        Gh_AN = natn[Gc_AN][h_AN];
        wan2 = WhatSpecies[Gh_AN];
        TNO2 = Spe_Total_CNO[wan2];

	for (i=0; i<TNO1; i++){
	  for (j=0; j<TNO2; j++){
	    Tmp_Vec[num] = OLPpoz[Mc_AN][h_AN][i][j];
	    num++;
	  }
	}
      }

      if (myid!=Host_ID){
        MPI_Isend(&num, 1, MPI_INT, Host_ID, tag, mpi_comm_level1, &request);
        MPI_Wait(&request,&stat);
        MPI_Isend(&Tmp_Vec[0], num, MPI_DOUBLE, Host_ID, tag, mpi_comm_level1, &request);
        MPI_Wait(&request,&stat);
      }
      else{
        fwrite(Tmp_Vec, sizeof(double), num, fp);
      }
    }

    else if (ID!=myid && myid==Host_ID){
      MPI_Recv(&num, 1, MPI_INT, ID, tag, mpi_comm_level1, &stat);
      MPI_Recv(&Tmp_Vec[0], num, MPI_DOUBLE, ID, tag, mpi_comm_level1, &stat);
      fwrite(Tmp_Vec, sizeof(double), num, fp);
    }
  }  

  /***************************************************************
                         density matrix
  ****************************************************************/

  for (spin=0; spin<=SpinP_switch; spin++){

    for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
      ID = G2ID[Gc_AN];

      if (myid==ID){

        num = 0;

        Mc_AN = F_G2M[Gc_AN];
        wan1 = WhatSpecies[Gc_AN];
        TNO1 = Spe_Total_CNO[wan1];
        for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
          Gh_AN = natn[Gc_AN][h_AN];
          wan2 = WhatSpecies[Gh_AN];
          TNO2 = Spe_Total_CNO[wan2];

          for (i=0; i<TNO1; i++){
            for (j=0; j<TNO2; j++){
              Tmp_Vec[num] = DM[0][spin][Mc_AN][h_AN][i][j];
              num++;
	    }
          }
        }

        if (myid!=Host_ID){
          MPI_Isend(&num, 1, MPI_INT, Host_ID, tag, mpi_comm_level1, &request);
          MPI_Wait(&request,&stat);
          MPI_Isend(&Tmp_Vec[0], num, MPI_DOUBLE, Host_ID, tag, mpi_comm_level1, &request);
          MPI_Wait(&request,&stat);
	}
        else{
          fwrite(Tmp_Vec, sizeof(double), num, fp);
        }
      }

      else if (ID!=myid && myid==Host_ID){
        MPI_Recv(&num, 1, MPI_INT, ID, tag, mpi_comm_level1, &stat);
        MPI_Recv(&Tmp_Vec[0], num, MPI_DOUBLE, ID, tag, mpi_comm_level1, &stat);
        fwrite(Tmp_Vec, sizeof(double), num, fp);
      }

    }  
  }

  /* freeing of Tmp_Vec */
  free(Tmp_Vec);

  if (myid==Host_ID){

    /****************************************************
        Solver
    ****************************************************/

    if (PeriodicGamma_flag==1)  i_vec[0] = 3;
    else                        i_vec[0] = Solver;
    fwrite(i_vec,sizeof(int),1,fp);

    /****************************************************
        ChemP
        Electronic Temp
        dipole moment (x, y, z) from core charge
        dipole moment (x, y, z) from back ground charge
    ****************************************************/

    d_vec[0] = ChemP;
    d_vec[1] = E_Temp*eV2Hartree;
    d_vec[2] = dipole_moment[1][1]; 
    d_vec[3] = dipole_moment[1][2];
    d_vec[4] = dipole_moment[1][3];
    d_vec[5] = dipole_moment[3][1];
    d_vec[6] = dipole_moment[3][2];
    d_vec[7] = dipole_moment[3][3];
    d_vec[8] = Total_Num_Electrons;
    d_vec[9] = Total_SpinS;
    fwrite(d_vec,sizeof(double),10,fp);
  }

  /****************************************************
      input file 
  ****************************************************/

  if (myid==Host_ID){

    /* find the number of lines in the input file */

    if ((fp_inp = fopen(inputfile,"r")) != NULL){

#ifdef xt3
      setvbuf(fp_inp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      count = 0;
      /* fgets gets a carriage return */
      while( fgets(strg, MAX_LINE_SIZE, fp_inp ) != NULL ){
	count++;
      }

      fclose(fp_inp);
    }
    else{
      printf("error #1 in saving *.scfout\n"); 
    }

    /* write the input file */

    if ((fp_inp = fopen(inputfile,"r")) != NULL){

#ifdef xt3
      setvbuf(fp_inp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      i_vec[0] = count;
      fwrite(i_vec, sizeof(int), 1, fp);

      /* fgets gets a carriage return */
      while( fgets(strg, MAX_LINE_SIZE, fp_inp ) != NULL ){
	fwrite(strg, sizeof(char), MAX_LINE_SIZE, fp);
      }

      fclose(fp_inp);
    }
    else{
      printf("error #1 in saving *.scfout\n"); 
    }
  }

}





void Calc_OLPpo()
{
  double time0;
  int Mc_AN,Gc_AN,Mh_AN,h_AN,Gh_AN;
  int i,j,Cwan,Hwan,NO0,NO1,spinmax;
  int Rnh,Rnk,spin,N,NumC[4];
  int n1,n2,n3,L0,Mul0,M0,L1,Mul1,M1;
  int Nc,GNc,GRc,Nog,Nh,MN,XC_P_switch;
  double x,y,z,dx,dy,dz,tmpx,tmpy,tmpz;
  double bc,dv,r,theta,phi,sum,tmp0,tmp1;
  double xo,yo,zo,S_coordinate[3];
  double Cxyz[4];
  double **ChiVx;
  double **ChiVy;
  double **ChiVz;
  double *tmp_ChiVx;
  double *tmp_ChiVy;
  double *tmp_ChiVz;
  double *tmp_Orbs_Grid;
  double **tmp_OLPpox;
  double **tmp_OLPpoy;
  double **tmp_OLPpoz;
  double TStime,TEtime;
  double TStime0,TEtime0;
  double TStime1,TEtime1;
  double TStime2,TEtime2;
  double TStime3,TEtime3;
  int numprocs,myid,tag=999,ID;
  double Stime_atom, Etime_atom;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  /****************************************************
    allocation of arrays:
  ****************************************************/

  ChiVx = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    ChiVx[i] = (double*)malloc(sizeof(double)*List_YOUSO[11]);
  }

  ChiVy = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    ChiVy[i] = (double*)malloc(sizeof(double)*List_YOUSO[11]);
  }

  ChiVz = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    ChiVz[i] = (double*)malloc(sizeof(double)*List_YOUSO[11]);
  }

  tmp_ChiVx = (double*)malloc(sizeof(double)*List_YOUSO[7]);
  tmp_ChiVy = (double*)malloc(sizeof(double)*List_YOUSO[7]);
  tmp_ChiVz = (double*)malloc(sizeof(double)*List_YOUSO[7]);

  tmp_Orbs_Grid = (double*)malloc(sizeof(double)*List_YOUSO[7]);

  tmp_OLPpox = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    tmp_OLPpox[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
  }

  tmp_OLPpoy = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    tmp_OLPpoy[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
  }

  tmp_OLPpoz = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (i=0; i<List_YOUSO[7]; i++){
    tmp_OLPpoz[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
  }

  /*****************************************************
      calculation of matrix elements for OLPpox,y,z
  *****************************************************/

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    dtime(&Stime_atom);

    Gc_AN = M2G[Mc_AN];    
    Cwan = WhatSpecies[Gc_AN];
    NO0 = Spe_Total_CNO[Cwan];
    
    for (i=0; i<NO0; i++){
      for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){

        GNc = GridListAtom[Mc_AN][Nc];
        GRc = CellListAtom[Mc_AN][Nc];

        Get_Grid_XYZ(GNc,Cxyz);
        x = Cxyz[1] + atv[GRc][1] - Gxyz[Gc_AN][1]; 
        y = Cxyz[2] + atv[GRc][2] - Gxyz[Gc_AN][2]; 
        z = Cxyz[3] + atv[GRc][3] - Gxyz[Gc_AN][3];

	ChiVx[i][Nc] = x*Orbs_Grid[Mc_AN][Nc][i];/* AITUNE */
	ChiVy[i][Nc] = y*Orbs_Grid[Mc_AN][Nc][i];/* AITUNE */
	ChiVz[i][Nc] = z*Orbs_Grid[Mc_AN][Nc][i];/* AITUNE */
      }
    }
    
    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

      Gh_AN = natn[Gc_AN][h_AN];
      Mh_AN = F_G2M[Gh_AN];

      Rnh = ncn[Gc_AN][h_AN];
      Hwan = WhatSpecies[Gh_AN];
      NO1 = Spe_Total_CNO[Hwan];

      /* initialize */

      for (i=0; i<NO0; i++){
	for (j=0; j<NO1; j++){
	  tmp_OLPpox[i][j] = 0.0;
	  tmp_OLPpoy[i][j] = 0.0;
	  tmp_OLPpoz[i][j] = 0.0;
	}
      }

      /* summation of non-zero elements */

      for (Nog=0; Nog<NumOLG[Mc_AN][h_AN]; Nog++){

        Nc = GListTAtoms1[Mc_AN][h_AN][Nog];
        Nh = GListTAtoms2[Mc_AN][h_AN][Nog];

        /* store ChiVx,y,z in tmp_ChiVx,y,z */

        for (i=0; i<NO0; i++){
          tmp_ChiVx[i] = ChiVx[i][Nc];
          tmp_ChiVy[i] = ChiVy[i][Nc];
          tmp_ChiVz[i] = ChiVz[i][Nc];
	}

        /* store Orbs_Grid in tmp_Orbs_Grid */

        if (G2ID[Gh_AN]==myid){
	  for (j=0; j<NO1; j++){
	    tmp_Orbs_Grid[j] = Orbs_Grid[Mh_AN][Nh][j];/* AITUNE */
	  }
	}
        else{
	  for (j=0; j<NO1; j++){
	    tmp_Orbs_Grid[j] = Orbs_Grid_FNAN[Mc_AN][h_AN][Nog][j];/* AITUNE */ 
	  }
        }

        /* integration */

        for (i=0; i<NO0; i++){
          tmpx = tmp_ChiVx[i]; 
          tmpy = tmp_ChiVy[i]; 
          tmpz = tmp_ChiVz[i]; 
          for (j=0; j<NO1; j++){
            tmp_OLPpox[i][j] += tmpx*tmp_Orbs_Grid[j];
            tmp_OLPpoy[i][j] += tmpy*tmp_Orbs_Grid[j];
            tmp_OLPpoz[i][j] += tmpz*tmp_Orbs_Grid[j];
	  }
	}
      }

      /* OLPpox,y,z */

      for (i=0; i<NO0; i++){
        for (j=0; j<NO1; j++){
          OLPpox[Mc_AN][h_AN][i][j] = tmp_OLPpox[i][j]*GridVol;
          OLPpoy[Mc_AN][h_AN][i][j] = tmp_OLPpoy[i][j]*GridVol;
          OLPpoz[Mc_AN][h_AN][i][j] = tmp_OLPpoz[i][j]*GridVol;
        }
      }

    }

    dtime(&Etime_atom);
    time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
  }
  
  /****************************************************
    freeing of arrays:
  ****************************************************/

  for (i=0; i<List_YOUSO[7]; i++){
    free(ChiVx[i]);
  }
  free(ChiVx);

  for (i=0; i<List_YOUSO[7]; i++){
    free(ChiVy[i]);
  }
  free(ChiVy);

  for (i=0; i<List_YOUSO[7]; i++){
    free(ChiVz[i]);
  }
  free(ChiVz);

  free(tmp_ChiVx);
  free(tmp_ChiVy);
  free(tmp_ChiVz);

  free(tmp_Orbs_Grid);

  for (i=0; i<List_YOUSO[7]; i++){
    free(tmp_OLPpox[i]);
  }
  free(tmp_OLPpox);

  for (i=0; i<List_YOUSO[7]; i++){
    free(tmp_OLPpoy[i]);
  }
  free(tmp_OLPpoy);

  for (i=0; i<List_YOUSO[7]; i++){
    free(tmp_OLPpoz[i]);
  }
  free(tmp_OLPpoz);
}
