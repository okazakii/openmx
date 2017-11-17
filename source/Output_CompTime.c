/**********************************************************************
  Output_CompTime.c:

     Output_CompTime.c is a subrutine to write computational time
     to filename.TRN.

  Log of Output_CompTime.c:

     22/Nov/2001  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "openmx_common.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

#define Num_CompTime 16

void Output_CompTime()
{
  int ID,j;
  int ID_Mintime1,ID_Maxtime1;
  int MinID[Num_CompTime];
  int MaxID[Num_CompTime];
  double Mintime1,Maxtime1;
  double *time0;
  double *time1;
  double MinCompTime[Num_CompTime];
  double MaxCompTime[Num_CompTime];
  char file_CompTime[YOUSO10] = ".CompTime";
  FILE *fp;
  int numprocs,myid;

  /* MPI */
  if (atomnum<=MYID_MPI_COMM_WORLD) return;
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /* allocate arrays */

  time0 = (double*)malloc(sizeof(double)*numprocs);
  time1 = (double*)malloc(sizeof(double)*numprocs);

  /* MPI_Bcast CompTime */

  for (ID=0; ID<numprocs; ID++){
    MPI_Bcast(&CompTime[ID][0], Num_CompTime, MPI_DOUBLE, ID, mpi_comm_level1);
  }

  MPI_Barrier(mpi_comm_level1);

  /* find the min and max CompTime */
  for (j=0; j<Num_CompTime; j++){
    MaxCompTime[j] = -10.0;
    MinCompTime[j] = 1.0e+10;

    for (ID=0; ID<numprocs; ID++){ 
      if (CompTime[ID][j]<MinCompTime[j]){
        MinCompTime[j] = CompTime[ID][j];
        MinID[j] = ID;
      }         
      if (MaxCompTime[j]<CompTime[ID][j]){
        MaxCompTime[j] = CompTime[ID][j];
        MaxID[j] = ID;
      }         
    }
  }

  if (myid==Host_ID){

    fnjoint(filepath,filename,file_CompTime);

    for (ID=0; ID<numprocs; ID++){
      time0[ID] = CompTime[ID][5]
                + CompTime[ID][6]
                + CompTime[ID][7]
                + CompTime[ID][8]
                + CompTime[ID][9]
                + CompTime[ID][10]
                + CompTime[ID][11]
                + CompTime[ID][12]
                + CompTime[ID][13]
                + CompTime[ID][14]
                + CompTime[ID][15];
      time1[ID] = CompTime[ID][3] - time0[ID];
    }

    Maxtime1 = -100.0; 
    Mintime1 = 1.0e+10;

    for (ID=0; ID<numprocs; ID++){
      if (time1[ID]<Mintime1){
        Mintime1 = time1[ID];
        ID_Mintime1 = ID;
      }         
      if (Maxtime1<time1[ID]){
        Maxtime1 = time1[ID];
        ID_Maxtime1 = ID;
      }         
    }    

    if ((fp = fopen(file_CompTime,"w")) != NULL){

      /* write */ 

      fprintf(fp,"\n");
      fprintf(fp,"***********************************************************\n");
      fprintf(fp,"***********************************************************\n");
      fprintf(fp,"               Computational Time (second)                 \n");
      fprintf(fp,"***********************************************************\n");
      fprintf(fp,"***********************************************************\n\n");

      fprintf(fp,"   Elapsed.Time.  %12.3f\n\n",MaxCompTime[0]);

      fprintf(fp,"                            Min_ID    Min_Time    Max_ID    Max_Time\n");

      fprintf(fp,"   Total Computational Time = %2d %12.3f      %2d %12.3f\n",
              MinID[0],MinCompTime[0],MaxID[0],MaxCompTime[0]);
      fprintf(fp,"   readfile                 = %2d %12.3f      %2d %12.3f\n",
              MinID[1],MinCompTime[1],MaxID[1],MaxCompTime[1]);
      fprintf(fp,"   truncation               = %2d %12.3f      %2d %12.3f\n",
              MinID[2],MinCompTime[2],MaxID[2],MaxCompTime[2]);
      fprintf(fp,"   MD_pac                   = %2d %12.3f      %2d %12.3f\n",
              MinID[4],MinCompTime[4],MaxID[4],MaxCompTime[4]);
      fprintf(fp,"   DFT                      = %2d %12.3f      %2d %12.3f\n",
              MinID[3],MinCompTime[3],MaxID[3],MaxCompTime[3]);
      fprintf(fp,"\n");
      fprintf(fp,"*** In DFT ***\n\n");
      fprintf(fp,"   Set_OLP_Kin              = %2d %12.3f      %2d %12.3f\n",
              MinID[5],MinCompTime[5],MaxID[5],MaxCompTime[5]);
      fprintf(fp,"   Set_Nonlocal             = %2d %12.3f      %2d %12.3f\n",
              MinID[6],MinCompTime[6],MaxID[6],MaxCompTime[6]);
      fprintf(fp,"   Set_Hamiltonian          = %2d %12.3f      %2d %12.3f\n",
              MinID[7],MinCompTime[7],MaxID[7],MaxCompTime[7]);
      fprintf(fp,"   Poisson                  = %2d %12.3f      %2d %12.3f\n",
              MinID[8],MinCompTime[8],MaxID[8],MaxCompTime[8]);
      fprintf(fp,"   Diagonalization          = %2d %12.3f      %2d %12.3f\n",
              MinID[9],MinCompTime[9],MaxID[9],MaxCompTime[9]);
      fprintf(fp,"   Mixing_DM                = %2d %12.3f      %2d %12.3f\n",
              MinID[10],MinCompTime[10],MaxID[10],MaxCompTime[10]);
      fprintf(fp,"   Force                    = %2d %12.3f      %2d %12.3f\n",
              MinID[11],MinCompTime[11],MaxID[11],MaxCompTime[11]);
      fprintf(fp,"   Total_Energy             = %2d %12.3f      %2d %12.3f\n",
              MinID[12],MinCompTime[12],MaxID[12],MaxCompTime[12]);
      fprintf(fp,"   Set_Aden_Grid            = %2d %12.3f      %2d %12.3f\n",
              MinID[13],MinCompTime[13],MaxID[13],MaxCompTime[13]);
      fprintf(fp,"   Set_Orbitals_Grid        = %2d %12.3f      %2d %12.3f\n",
              MinID[14],MinCompTime[14],MaxID[14],MaxCompTime[14]);
      fprintf(fp,"   Set_Density_Grid         = %2d %12.3f      %2d %12.3f\n",
              MinID[15],MinCompTime[15],MaxID[15],MaxCompTime[15]);
      fprintf(fp,"   Others                   = %2d %12.3f      %2d %12.3f\n",
              ID_Mintime1,Mintime1,ID_Maxtime1,Maxtime1);
      fclose(fp);

      /* stdout */ 

      if (0<level_stdout){

	printf("\n");
	printf("***********************************************************\n");
	printf("***********************************************************\n");
	printf("               Computational Time (second)                 \n");
	printf("***********************************************************\n");
	printf("***********************************************************\n\n");

	printf("                            Min_ID    Min_Time    Max_ID    Max_Time\n");

	printf("   Total Computational Time = %2d %12.3f      %2d %12.3f\n",
	       MinID[0],MinCompTime[0],MaxID[0],MaxCompTime[0]);
	printf("   readfile                 = %2d %12.3f      %2d %12.3f\n",
	       MinID[1],MinCompTime[1],MaxID[1],MaxCompTime[1]);
	printf("   truncation               = %2d %12.3f      %2d %12.3f\n",
	       MinID[2],MinCompTime[2],MaxID[2],MaxCompTime[2]);
	printf("   MD_pac                   = %2d %12.3f      %2d %12.3f\n",
	       MinID[4],MinCompTime[4],MaxID[4],MaxCompTime[4]);
	printf("   DFT                      = %2d %12.3f      %2d %12.3f\n",
	       MinID[3],MinCompTime[3],MaxID[3],MaxCompTime[3]);
	printf("\n");
	printf("*** In DFT ***\n\n");
	printf("   Set_OLP_Kin              = %2d %12.3f      %2d %12.3f\n",
	       MinID[5],MinCompTime[5],MaxID[5],MaxCompTime[5]);
	printf("   Set_Nonlocal             = %2d %12.3f      %2d %12.3f\n",
	       MinID[6],MinCompTime[6],MaxID[6],MaxCompTime[6]);
	printf("   Set_Hamiltonian          = %2d %12.3f      %2d %12.3f\n",
	       MinID[7],MinCompTime[7],MaxID[7],MaxCompTime[7]);
	printf("   Poisson                  = %2d %12.3f      %2d %12.3f\n",
	       MinID[8],MinCompTime[8],MaxID[8],MaxCompTime[8]);
	printf("   Diagonalization          = %2d %12.3f      %2d %12.3f\n",
	       MinID[9],MinCompTime[9],MaxID[9],MaxCompTime[9]);
	printf("   Mixing_DM                = %2d %12.3f      %2d %12.3f\n",
	       MinID[10],MinCompTime[10],MaxID[10],MaxCompTime[10]);
	printf("   Force                    = %2d %12.3f      %2d %12.3f\n",
	       MinID[11],MinCompTime[11],MaxID[11],MaxCompTime[11]);
	printf("   Total_Energy             = %2d %12.3f      %2d %12.3f\n",
	       MinID[12],MinCompTime[12],MaxID[12],MaxCompTime[12]);
	printf("   Set_Aden_Grid            = %2d %12.3f      %2d %12.3f\n",
	       MinID[13],MinCompTime[13],MaxID[13],MaxCompTime[13]);
	printf("   Set_Orbitals_Grid        = %2d %12.3f      %2d %12.3f\n",
	       MinID[14],MinCompTime[14],MaxID[14],MaxCompTime[14]);
	printf("   Set_Density_Grid         = %2d %12.3f      %2d %12.3f\n",
	       MinID[15],MinCompTime[15],MaxID[15],MaxCompTime[15]);
	printf("   Others                   = %2d %12.3f      %2d %12.3f\n",
	       ID_Mintime1,Mintime1,ID_Maxtime1,Maxtime1);

      }

    } 
    else{
      printf("could not save the CompTime file.\n");
    }
  }

  /* free arrays */

  free(time0);
  free(time1);

}


