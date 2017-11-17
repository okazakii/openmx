/**********************************************************************
  Output_Energy_Decomposition.c:

     Output_Energy_Decomposition.c is a subrutine to output 
     decomposed energies.
 
  Log of Output_Energy_Decomposition.c:

     03/Aug./2014  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "openmx_common.h"
#include "mpi.h"


void Output_Energy_Decomposition()
{ 
  int Mc_AN,Gc_AN,tno0,Cwan,num,l,m,n,mul,Nene;
  int wan1,wan2,i,j,k,spin,tag=999;
  double Stime_atom,Etime_atom;
  double sum0,sum1,sum,time0;
  double ***DecE,*tmp_array0;
  char *Name_Angular[Supported_MaxL+1][2*(Supported_MaxL+1)+1];
  char file_DecE[YOUSO10];
  FILE *fp_DecE;
  int numprocs,myid,ID;
  char buf[fp_bsize];          /* setvbuf */
  MPI_Status stat;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  Nene = 12;

  /* allocation of arrays */

  DecE = (double***)malloc(sizeof(double**)*(atomnum+1));
  for (Gc_AN=0; Gc_AN<=atomnum; Gc_AN++){
    DecE[Gc_AN] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
    for (i=0; i<List_YOUSO[7]; i++){
      DecE[Gc_AN][i] = (double*)malloc(sizeof(double)*Nene);
      for (n=0; n<Nene; n++) DecE[Gc_AN][i][n] = 0.0;
    }
  }

  tmp_array0 = (double*)malloc(sizeof(double)*List_YOUSO[7]*Nene);

  /****************************************
   MPI: DecE
  ****************************************/

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){

    ID = G2ID[Gc_AN];
    wan1 = WhatSpecies[Gc_AN];

    /* sending from ID to Host_ID */ 

    if (myid==ID){
  
      Mc_AN = F_G2M[Gc_AN];  

      num = 0; 

      if (SpinP_switch==0){

	for (i=0; i<Spe_Total_CNO[wan1]; i++){
	  tmp_array0[num] = 2.0*DecEkin[0][Mc_AN][i]; num++;
	  tmp_array0[num] = 2.0*DecEna[ 0][Mc_AN][i]; num++;
	  tmp_array0[num] = 2.0*DecEnl[ 0][Mc_AN][i]; num++;
	  tmp_array0[num] = 2.0*DecEdee[0][Mc_AN][i]; num++;
	  tmp_array0[num] = 2.0*DecExc[ 0][Mc_AN][i]; num++;
	  tmp_array0[num] = 2.0*DecEscc[0][Mc_AN][i]; num++;
	  tmp_array0[num] = 2.0*DecEef[ 0][Mc_AN][i]; num++;
	  tmp_array0[num] = 2.0*DecEhub[0][Mc_AN][i]; num++;
	  tmp_array0[num] = 2.0*DecEcs[ 0][Mc_AN][i]; num++;
	  tmp_array0[num] = 2.0*DecEzs[ 0][Mc_AN][i]; num++;
	  tmp_array0[num] = 2.0*DecEzo[ 0][Mc_AN][i]; num++;
	  tmp_array0[num] = 2.0*DecEvdw[0][Mc_AN][i]; num++;
	}
      }

      else{

	for (i=0; i<Spe_Total_CNO[wan1]; i++){
	  tmp_array0[num] = DecEkin[0][Mc_AN][i] + DecEkin[1][Mc_AN][i]; num++;
	  tmp_array0[num] = DecEna[ 0][Mc_AN][i] + DecEna[ 1][Mc_AN][i]; num++;
	  tmp_array0[num] = DecEnl[ 0][Mc_AN][i] + DecEnl[ 1][Mc_AN][i]; num++;
	  tmp_array0[num] = DecEdee[0][Mc_AN][i] + DecEdee[1][Mc_AN][i]; num++;
	  tmp_array0[num] = DecExc[ 0][Mc_AN][i] + DecExc[ 1][Mc_AN][i]; num++;
	  tmp_array0[num] = DecEscc[0][Mc_AN][i] + DecEscc[1][Mc_AN][i]; num++;
	  tmp_array0[num] = DecEef[ 0][Mc_AN][i] + DecEef[ 1][Mc_AN][i]; num++;
	  tmp_array0[num] = DecEhub[0][Mc_AN][i] + DecEhub[1][Mc_AN][i]; num++;
	  tmp_array0[num] = DecEcs[ 0][Mc_AN][i] + DecEcs[ 1][Mc_AN][i]; num++;
	  tmp_array0[num] = DecEzs[ 0][Mc_AN][i] + DecEzs[ 1][Mc_AN][i]; num++;
	  tmp_array0[num] = DecEzo[ 0][Mc_AN][i] + DecEzo[ 1][Mc_AN][i]; num++;
	  tmp_array0[num] = DecEvdw[0][Mc_AN][i] + DecEvdw[1][Mc_AN][i]; num++;
	}
      }

      if (myid!=Host_ID){
        MPI_Send(&tmp_array0[0], num, MPI_DOUBLE, Host_ID, tag, mpi_comm_level1);
      }
    }

    /* receiving at Host_ID from ID */ 

    if (myid==Host_ID){

      num = Nene*Spe_Total_CNO[wan1]; 

      if (ID!=Host_ID){
	MPI_Recv(&tmp_array0[0], num, MPI_DOUBLE, ID, tag, mpi_comm_level1, &stat);
      }

      num = 0; 
      for (i=0; i<Spe_Total_CNO[wan1]; i++){
        for (n=0; n<Nene; n++){
          DecE[Gc_AN][i][n] = tmp_array0[num];
	  num++;
        }
      }

    }
  }

  /****************************************
             making a file, *.DecE
  ****************************************/

  if ( myid==Host_ID ){

    sprintf(file_DecE,"%s%s.DecE",filepath,filename);

    if ((fp_DecE = fopen(file_DecE,"w")) != NULL){

#ifdef xt3
      setvbuf(fp_DecE,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      fprintf(fp_DecE,"\n");
      fprintf(fp_DecE,"***********************************************************\n");
      fprintf(fp_DecE,"***********************************************************\n");
      fprintf(fp_DecE,"            Decomposed energies in Hartree unit            \n");
      fprintf(fp_DecE,"***********************************************************\n");
      fprintf(fp_DecE,"***********************************************************\n\n");

      /* total energy */

      sum = 0.0;
      for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
	wan1 = WhatSpecies[Gc_AN];
	for (i=0; i<Spe_Total_CNO[wan1]; i++){
          for (n=0; n<Nene; n++){
	    sum += DecE[Gc_AN][i][n];
	  }
	}
      }

      fprintf(fp_DecE,"  Total energy (Hartree) = %18.15f\n",sum);

      /* Decomposed energies with respect to atom */

      fprintf(fp_DecE,"\n  Decomposed energies (Hartree) with respect to atom\n\n");
      fprintf(fp_DecE,"                Utot              Ukin       Una        Unl         UH1       Uxc        Ucore+UH0   Uef        Uhub       Ucs        Uzs        Uzo        Uvdw\n");
      
      for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){

	wan1 = WhatSpecies[Gc_AN];

        sum = 0.0;
        for (i=0; i<Spe_Total_CNO[wan1]; i++){
          for (n=0; n<Nene; n++){
	    sum += DecE[Gc_AN][i][n];
	  }
	}

	fprintf(fp_DecE,"  %4d %4s  %18.12f",Gc_AN,SpeName[wan1],sum);

        for (n=0; n<Nene; n++){
          sum = 0.0;
  	  for (i=0; i<Spe_Total_CNO[wan1]; i++){
	    sum += DecE[Gc_AN][i][n];
	  }

          fprintf(fp_DecE," %10.5f",sum);
	}
	fprintf(fp_DecE,"\n");
      }

      /* Decomposed energies with respect to atomic orbital */

      Name_Angular[0][0] = "s          ";
      Name_Angular[1][0] = "px         ";
      Name_Angular[1][1] = "py         ";
      Name_Angular[1][2] = "pz         ";
      Name_Angular[2][0] = "d3z^2-r^2  ";
      Name_Angular[2][1] = "dx^2-y^2   ";
      Name_Angular[2][2] = "dxy        ";
      Name_Angular[2][3] = "dxz        ";
      Name_Angular[2][4] = "dyz        ";
      Name_Angular[3][0] = "f5z^2-3r^2 ";
      Name_Angular[3][1] = "f5xz^2-xr^2";
      Name_Angular[3][2] = "f5yz^2-yr^2";
      Name_Angular[3][3] = "fzx^2-zy^2 ";
      Name_Angular[3][4] = "fxyz       ";
      Name_Angular[3][5] = "fx^3-3*xy^2";
      Name_Angular[3][6] = "f3yx^2-y^3 ";
      Name_Angular[4][0] = "g1         ";
      Name_Angular[4][1] = "g2         ";
      Name_Angular[4][2] = "g3         ";
      Name_Angular[4][3] = "g4         ";
      Name_Angular[4][4] = "g5         ";
      Name_Angular[4][5] = "g6         ";
      Name_Angular[4][6] = "g7         ";
      Name_Angular[4][7] = "g8         ";
      Name_Angular[4][8] = "g9         ";

      fprintf(fp_DecE,"\n\n  Decomposed energies (Hartree) with respect to atomic orbital\n");

      for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){

        wan1 = WhatSpecies[Gc_AN];

        fprintf(fp_DecE,"\n %4d %4s          Utot       Ukin       Una        Unl         UH1       Uxc        Ucore+UH0  Uef        Uhub       Ucs        Uzs        Uzo        Uvdw\n",
                Gc_AN,SpeName[wan1]);
	fprintf(fp_DecE,"            multiple\n");

	num = 0;
	for (l=0; l<=Supported_MaxL; l++){

	  for (mul=0; mul<Spe_Num_CBasis[wan1][l]; mul++){
	    for (m=0; m<(2*l+1); m++){

	      fprintf(fp_DecE,"  %s%2d ",Name_Angular[l][m],mul);
 
              sum = 0.0;
              for (n=0; n<Nene; n++) sum += DecE[Gc_AN][num][n];
              fprintf(fp_DecE," %10.5f",sum);

              for (n=0; n<Nene; n++){
                fprintf(fp_DecE," %10.5f",DecE[Gc_AN][num][n]);
	      }
              fprintf(fp_DecE,"\n");

	      num++;
	    }
	  } 
	}
      }

      fclose(fp_DecE);
    }
    else{
      printf("Failure of saving the DecE file.\n");
    }
  }

  /* freeing of arrays */

  for (Gc_AN=0; Gc_AN<=atomnum; Gc_AN++){
    for (i=0; i<List_YOUSO[7]; i++){
      free(DecE[Gc_AN][i]);
    }
    free(DecE[Gc_AN]);
  }
  free(DecE);

  free(tmp_array0);
}

