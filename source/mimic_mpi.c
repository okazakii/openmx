#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef nompi
#include "mimic_mpi.h"

int MPI_Init(int *argc, char **argv[])
{}

int MPI_Comm_size(MPI_Comm comm, int *numprocs)
{
  *numprocs = 1; 
}

int MPI_Comm_rank(MPI_Comm comm, int *myid)
{
  *myid = 0; 
}


int MPI_Isend( void *buf, int count, MPI_Datatype datatype, int dest, int tag,
	       MPI_Comm comm, MPI_Request *request )
{}

int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request)
{}

int MPI_Waitall(int aaa, MPI_Request *request, MPI_Status *status)
{}


int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
{}


int MPI_Recv( void *buf, int count, MPI_Datatype datatype, int source, 
	      int tag, MPI_Comm comm, MPI_Status *status )
{}

int MPI_Wait(MPI_Request *request, MPI_Status *stat)
{}

int MPI_Finalize()
{}

int MPI_Abort(MPI_Comm comm, int count)
{}

int MPI_Bcast ( void *buffer, int count, MPI_Datatype datatype, int root, 
		MPI_Comm comm )
{}

int MPI_Reduce ( void *sendbuf, void *recvbuf, int count, 
		 MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm )
{

  static int i;
  static int *is,*ir;
  static long int *lis,*lir;
  static double *ds,*dr;

  if (datatype==6){ /* MPI_INT */
    is = sendbuf;
    ir = recvbuf;
    for (i=0; i<count; i++) *ir++ = *is++;
  }

  else if (datatype==8){ /* MPI_LONG */
    lis = sendbuf;
    lir = recvbuf;
    for (i=0; i<count; i++) *lir++ = *lis++;
  }

  else if (datatype==11){ /* MPI_DOUBLE */
    ds = sendbuf;
    dr = recvbuf;
    for (i=0; i<count; i++) *dr++ = *ds++;
  }

}

int MPI_Allreduce ( void *sendbuf, void *recvbuf, int count,
		    MPI_Datatype datatype, MPI_Op op, MPI_Comm comm )
{

  static int i;
  static int *is,*ir;
  static long int *lis,*lir;
  static double *ds,*dr;

  if (datatype==6){
    is = sendbuf;
    ir = recvbuf;
    for (i=0; i<count; i++) *ir++ = *is++;
  }

  else if (datatype==8){
    lis = sendbuf;
    lir = recvbuf;
    for (i=0; i<count; i++) *lir++ = *lis++;
  }

  else if (datatype==11){
    ds = sendbuf;
    dr = recvbuf;
    for (i=0; i<count; i++) *dr++ = *ds++;
  }

}

int MPI_Barrier ( MPI_Comm comm )
{}


int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *comm1)
{}

int MPI_Group_incl(MPI_Group group, int i1, int *j1, MPI_Group *group2)
{}

int MPI_Comm_group(MPI_Comm comm, MPI_Group *group)
{}

int MPI_Comm_free(MPI_Comm *comm)
{}

int MPI_Group_free(MPI_Comm *comm)
{}

#endif

