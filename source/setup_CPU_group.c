#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "openmx_common.h"
#include "Inputtools.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif


void setup_CPU_group(char *file)
{
  char *s_vec[20];
  int  i_vec[20];
  int po=0;
  int numprocs,myid; 

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (input_open(file)==0){
    MPI_Finalize(); 
    exit(0);
  }
  else{
    input_close();
  }

  if (myid==Host_ID) {

    input_open(file);

    /* Solver */

    s_vec[0]="Recursion";     s_vec[1]="Cluster"; s_vec[2]="Band";
    s_vec[3]="NEGF";          s_vec[4]="DC";      s_vec[5]="GDC";
    s_vec[6]="Cluster-DIIS";  s_vec[7]="Krylov";  s_vec[8]="Cluster2";  

    i_vec[0]=1;  i_vec[1]=2;  i_vec[2]=3;
    i_vec[3]=4;  i_vec[4]=5;  i_vec[5]=6;
    i_vec[6]=7;  i_vec[7]=8;  i_vec[8]=9;

    input_string2int("scf.EigenvalueSolver", &Solver, 9, s_vec,i_vec);

    /* atomnum */

    input_int("Atoms.Number",&atomnum,0);
    input_double("scf.energycutoff",&Grid_Ecut,(double)150.0);
    input_logical("scf.Mixed.Basis",&Mixed_Basis_flag,0); /* default=off */

    if (Mixed_Basis_flag==1){
      Setup_Mixed_Basis(file,myid);
      atomnum += Ngrid1_FE*Ngrid2_FE*Ngrid3_FE;
    }

    if (atomnum<=0){
      printf("Atoms.Number may be wrong.\n");
      po++;
    }

    if (Solver==4) {
      /* if NEGF */
      /* left */
      input_int("LeftLeadAtoms.Number",&Latomnum,0);
      if (Latomnum<=0){
	printf("LeftLeadAtoms.Number may be wrong.\n");
	po++;
      }

      /* right */
      input_int("RightLeadAtoms.Number",&Ratomnum,0);
      if (Ratomnum<=0){
	printf("RightLeadAtoms.Number may be wrong.\n");
	po++;
      }

      atomnum += Latomnum+Ratomnum;
    }

    input_close();

  } /* myid==Host_ID */

  MPI_Bcast(&atomnum, 1, MPI_INT, Host_ID, mpi_comm_level1);
  MPI_Bcast(&po,      1, MPI_INT, Host_ID, mpi_comm_level1);

  if (po>0) {
     MPI_Finalize();
     exit(0);
  }

  /* separate CPUs */
  if ( atomnum < numprocs )  {
    int *new_ranks; 
    int i;
    MPI_Group new_group,old_group; 

    if (myid==Host_ID) {
       printf("******************************************\n");fflush(stdout);
       printf("Cut off CPU(s), New group contains %d CPUs\n",atomnum); fflush(stdout);
       printf("******************************************\n");fflush(stdout);
    }

    new_ranks = (int*)malloc(sizeof(int)*atomnum);
    for (i=0; i<atomnum; i++) {
      new_ranks[i]=i; /* a new group is made of original rank=0:atomnum-1 */
    }

    MPI_Comm_group(MPI_COMM_WORLD1, &old_group);

    /* define a new group */
    MPI_Group_incl(old_group,atomnum,new_ranks,&new_group);
    MPI_Comm_create(MPI_COMM_WORLD1,new_group,&mpi_comm_level1);

    MPI_Group_free(&new_group);
    free(new_ranks); /* never forget cleaning! */
  }
  /* default mpi_comm_level1 is already set in main() */
}

