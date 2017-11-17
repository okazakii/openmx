/**********************************************************************
  Make_Comm_Worlds.c:

   Make_Comm_Worlds.c is a subroutine to make new COMM worlds.

  Log of Make_Comm_Worlds.c:

     16/Dec/2006  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"



void Make_Comm_Worlds(
   MPI_Comm MPI_Curret_Comm_WD,   
   int myid0,
   int numprocs0,
   int Num_Comm_World, 
   int *myworld1, 
   MPI_Comm *MPI_CommWD,     /* size: Num_Comm_World */
   int *NPROCS1_ID,          /* size: numprocs0 */
   int *Comm_World1,         /* size: numprocs0 */
   int *NPROCS1_WD,          /* size: Num_Comm_World */
   int *Comm_World_StartID   /* size: Num_Comm_World */
   )
{
  int i,j;
  int ID0;
  int numprocs1,myid1;
  int num;
  double avnum;
  int *new_ranks; 
  MPI_Group new_group,old_group; 

  /******************************************
     Set up informations to construct the 
     (numprocs0/Num_Comm_World)-th worlds 
  ******************************************/

  if (Num_Comm_World<=numprocs0){

    avnum = (double)numprocs0/(double)Num_Comm_World;
    for (i=0; i<Num_Comm_World; i++){

      if ( (int)((double)(i)*avnum+1.0e-12)<=myid0
         && myid0<(int)((double)(i+1)*avnum+1.0e-12)){

        numprocs1 = (int)((double)(i+1)*avnum+1.0e-12) - (int)((double)(i)*avnum+1.0e-12);
        *myworld1 = i;
      }
    }

    for (i=0; i<Num_Comm_World; i++){
      num = 0;
      for (ID0=0; ID0<numprocs0; ID0++){

        if ( (int)((double)(i)*avnum+1.0e-12)<=ID0
           && ID0<(int)((double)(i+1)*avnum+1.0e-12)){

          NPROCS1_ID[ID0] = (int)((double)(i+1)*avnum+1.0e-12) - (int)((double)(i)*avnum+1.0e-12);
          Comm_World1[ID0] = i;
          if (num==0) Comm_World_StartID[i] = ID0; 
          num++;

        }
      }
 
      NPROCS1_WD[i] = num;
    }

    /**************************
     make a set of MPI_CommWD
    **************************/

    for (i=0; i<Num_Comm_World; i++){

      new_ranks = (int*)malloc(sizeof(int)*NPROCS1_WD[i]);

      for (j=0; j<NPROCS1_WD[i]; j++) {
	new_ranks[j] = Comm_World_StartID[i] + j;
      }

      MPI_Comm_group(MPI_Curret_Comm_WD, &old_group);

      /* define a new group */
      MPI_Group_incl(old_group,NPROCS1_WD[i],new_ranks,&new_group);
      MPI_Comm_create(MPI_Curret_Comm_WD,new_group,&MPI_CommWD[i]);

      MPI_Group_free(&new_group);
      free(new_ranks); /* never forget cleaning! */

    }
  }

  else {

    numprocs1 = numprocs0;
    *myworld1 = 0;
    for (ID0=0; ID0<numprocs0; ID0++){
      NPROCS1_ID[ID0] = numprocs0;
      Comm_World1[ID0] = 0;
    }
    for (i=0; i<Num_Comm_World; i++){
      Comm_World_StartID[i] = 0;
      NPROCS1_WD[i] = numprocs0;
    }

    for (i=0; i<Num_Comm_World; i++){
      MPI_CommWD[i] = MPI_Curret_Comm_WD;
    }
  }
}


