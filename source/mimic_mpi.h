/* Communicators */

#ifndef ___MPI_Comm_definition___
typedef int MPI_Comm;
#define ___MPI_Comm_definition___ 
#endif

#define MPI_COMM_WORLD 91
#define MPI_COMM_SELF  92
/* Groups */
typedef int MPI_Group;
#define MPI_GROUP_EMPTY 90

/* Data types
 * A more aggressive yet homogeneous implementation might want to 
 * make the values here the number of bytes in the basic type, with
 * a simple test against a max limit (e.g., 16 for long double), and
 * non-contiguous structures with indices greater than that.
 * 
 * Note: Configure knows these values for providing the Fortran optional
 * types (like MPI_REAL8).  Any changes here must be matched by changes
 * in configure.in
 */
typedef int MPI_Datatype;
#define MPI_CHAR           ((MPI_Datatype)1)
#define MPI_UNSIGNED_CHAR  ((MPI_Datatype)2)
#define MPI_BYTE           ((MPI_Datatype)3)
#define MPI_SHORT          ((MPI_Datatype)4)
#define MPI_UNSIGNED_SHORT ((MPI_Datatype)5)
#define MPI_INT            ((MPI_Datatype)6)
#define MPI_UNSIGNED       ((MPI_Datatype)7)
#define MPI_LONG           ((MPI_Datatype)8)
#define MPI_UNSIGNED_LONG  ((MPI_Datatype)9)
#define MPI_FLOAT          ((MPI_Datatype)10)
#define MPI_DOUBLE         ((MPI_Datatype)11)
#define MPI_LONG_DOUBLE    ((MPI_Datatype)12)
#define MPI_LONG_LONG_INT  ((MPI_Datatype)13)
/* MPI_LONG_LONG is in the complete ref 2nd edition, though not in the 
   standard.  Rather, MPI_LONG_LONG_INT is on page 40 in the HPCA version */
#define MPI_LONG_LONG      ((MPI_Datatype)13)
#define MPI_PACKED         ((MPI_Datatype)14)
#define MPI_LB             ((MPI_Datatype)15)
#define MPI_UB             ((MPI_Datatype)16)
/* 
   The layouts for the types MPI_DOUBLE_INT etc are simply
   struct { 
       double var;
       int    loc;
   }
   This is documented in the man pages on the various datatypes.   
 */
#define MPI_FLOAT_INT         ((MPI_Datatype)17)
#define MPI_DOUBLE_INT        ((MPI_Datatype)18)
#define MPI_LONG_INT          ((MPI_Datatype)19)
#define MPI_SHORT_INT         ((MPI_Datatype)20)
#define MPI_2INT              ((MPI_Datatype)21)
#define MPI_LONG_DOUBLE_INT   ((MPI_Datatype)22)

/* Collective operations */
typedef int MPI_Op;

#define MPI_MAX    (MPI_Op)(100)
#define MPI_MIN    (MPI_Op)(101)
#define MPI_SUM    (MPI_Op)(102)
#define MPI_PROD   (MPI_Op)(103)
#define MPI_LAND   (MPI_Op)(104)
#define MPI_BAND   (MPI_Op)(105)
#define MPI_LOR    (MPI_Op)(106)
#define MPI_BOR    (MPI_Op)(107)
#define MPI_LXOR   (MPI_Op)(108)
#define MPI_BXOR   (MPI_Op)(109)
#define MPI_MINLOC (MPI_Op)(110)
#define MPI_MAXLOC (MPI_Op)(111)

#ifndef ___MPI_Status_definition___
typedef struct MPIStatus{int i;}  MPI_Status;  
#define ___MPI_Status_definition___ 
#endif

#ifndef ___MPI_Request_definition___
typedef struct MPIRequest{int i;} MPI_Request;  
#define ___MPI_Request_definition___ 
#endif

int MPI_Init(int *, char ***);
int MPI_Comm_size(MPI_Comm comm, int *);
int MPI_Comm_rank(MPI_Comm comm, int *);
int MPI_Finalize(void);
int MPI_Abort(MPI_Comm comm, int count);
int MPI_Send(void*, int, MPI_Datatype, int, int, MPI_Comm);
int MPI_Recv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status *);
int MPI_Isend(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int MPI_Irecv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int MPI_Wait(MPI_Request *, MPI_Status *);
int MPI_Bcast(void*, int, MPI_Datatype, int, MPI_Comm );
int MPI_Reduce(void* , void*, int, MPI_Datatype, MPI_Op, int, MPI_Comm);
int MPI_Allreduce(void* , void*, int, MPI_Datatype, MPI_Op, MPI_Comm);
int MPI_Barrier(MPI_Comm );
int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *);
int MPI_Group_incl(MPI_Group group, int, int *, MPI_Group *);
int MPI_Comm_group(MPI_Comm comm, MPI_Group *);
int MPI_Waitall(int, MPI_Request *, MPI_Status *);
int MPI_Comm_free(MPI_Comm *comm);
int MPI_Group_free(MPI_Comm *comm);


MPI_Comm  mpi_comm_level1;
