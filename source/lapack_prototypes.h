#include "f77func.h"
#ifdef blaswrap
#define zgemm_ f2c_zgemm
#define zcopy_ f2c_zcopy
#endif

#ifndef ___INTEGER_definition___
typedef int INTEGER; /* for fortran integer */
#define ___INTEGER_definition___ 
#endif

#ifndef ___logical_definition___
typedef long int logical;
#define ___logical_definition___ 
#endif

#ifndef ___dcomplex_definition___
typedef struct { double r,i; } dcomplex;
#define ___dcomplex_definition___ 
#endif

int zgesvd_(char *jobu, char *jobvt, INTEGER *m, INTEGER *n, 
	    dcomplex *a, INTEGER *lda,
	    double *s, dcomplex *u, 
	    INTEGER *ldu, dcomplex *vt,
	    INTEGER *ldvt, dcomplex *work, 
	    INTEGER *lwork, double *rwork, INTEGER *info);

int dgesv_(INTEGER *n, INTEGER *nrhs, double *a, INTEGER *lda, INTEGER *ipiv, double *b,
           INTEGER *ldb, INTEGER *info);

void dsyev_(char *JOBZ, char *UPLO,INTEGER *N,double *A,INTEGER *LDA,double *W,
             double *WORK,INTEGER *LWORK, INTEGER *INFO);
void dsbev_(char *JOBZ, char *UPLO, INTEGER *N, INTEGER *KD, double *AB, INTEGER *LDAB,
         double *W, double *Z,
        INTEGER * LDZ,double *WORK, INTEGER *INFO );
void dgesvd_(char *JOBU, char *JOBVT, INTEGER *M, INTEGER *N, double *A, INTEGER *LDA, double *S, double *U,
        INTEGER *LDU, double *VT, INTEGER *LDVT, double *WORK, INTEGER *LWORK, INTEGER *INFO );

void dsygv_(INTEGER *itype, char *jobz, char *uplo, INTEGER *n, 
            double *a, INTEGER *lda, double *b, INTEGER *ldb, 
	    double *w, double *work, INTEGER *lwork, INTEGER *info);

void dsyevx_(char *JOBZ, char *RANGE, char *UPLO, INTEGER *N, double *A, INTEGER *LDA, 
          double *VL, double *VU,
          INTEGER *IL, INTEGER *IU, double *ABSTOL, INTEGER *M, double *W, double *Z, 
           INTEGER *LDZ, double *WORK,
          INTEGER *LWORK, INTEGER *IWORK, INTEGER *IFAIL, INTEGER *INFO);
void dsyevd_(char *JOBZ, char *UPLO, INTEGER *N, double *A, INTEGER *LDA, double *W, double *WORK, 
         INTEGER *LWORK, INTEGER *IWORK,
                       INTEGER  *LIWORK, INTEGER *INFO);
void dsyevr_(char *JOBZ, char *RANGE, char *UPLO, INTEGER *N, double *A, INTEGER *LDA, double *VL, double *VU, 
          INTEGER *IL, INTEGER *IU,
                       double  *ABSTOL, INTEGER *M, double *W, double *Z,INTEGER *LDZ, INTEGER *ISUPPZ, double *WORK, 
        INTEGER *LWORK, INTEGER *IWORK, INTEGER *LIWORK, INTEGER *INFO);

int dstevx_(char *JOBZ, char *RANGE, INTEGER *N, double *D, double *E, double *VL, double *VU, INTEGER *IL, INTEGER *IU, double *ABSTOL,
                        INTEGER *M, double *W, double *Z, INTEGER *LDZ, double *WORK, INTEGER *IWORK, INTEGER *IFAIL, INTEGER *INFO);

void dgetrf_(INTEGER *M, INTEGER *N,double *A, INTEGER *LDA, INTEGER *IPIV, INTEGER *INFO);
void dgetri_( INTEGER *N, double *A, INTEGER *LDA, INTEGER *IPIV, double *WORK, INTEGER *LWORK, INTEGER *INFO);

int dsysv_(char *UPLO, INTEGER *N, INTEGER *NRHS, double *A, INTEGER *LDA, INTEGER *IPIV, double *B, INTEGER *LDB, double *WORK,
                       INTEGER *LWORK, INTEGER *INFO);

int dstedc_(char *COMPZ , INTEGER *N , double *D , double *E , double *Z , 
       INTEGER *LDZ , double *WORK , INTEGER *LWORK , INTEGER *IWORK , INTEGER *LIWORK , 
       INTEGER *INFO);

int dstegr_(char *JOBZ , char *RANGE , INTEGER *N , double *D , double *E , 
       double *VL , double *VU , INTEGER *IL , INTEGER *IU , double *ABSTOL , 
       INTEGER *M , double *W , double *Z , INTEGER *LDZ , INTEGER *ISUPPZ , 
       double *WORK , INTEGER *LWORK , INTEGER *IWORK , INTEGER *LIWORK , INTEGER *INFO
       );

int dsteqr_(char *compz, INTEGER *n, double *d__, 
	    double *e, double *z__, INTEGER *ldz, double *work, INTEGER *info);


void zheev_(char *JOBZ, char *UPLO, INTEGER *N, dcomplex *A, INTEGER *LDA, double *W, dcomplex *WORK, INTEGER *LWORK, 
       double *RWORK, INTEGER *INFO);
void zheevx_(char *JOBZ, char *RANGE, char *UPLO, INTEGER *N, dcomplex *A, INTEGER *LDA, double *VL, double *VU, INTEGER *IL, INTEGER *IU,
                        double *ABSTOL, INTEGER *M, double *W, dcomplex *Z, INTEGER *LDZ, dcomplex *WORK, INTEGER *LWORK, double *RWORK,
                        INTEGER *IWORK, INTEGER *IFAIL, INTEGER *INFO);


void zheev_(char *JOBZ, char *UPLO, int *N, dcomplex *A, int *LDA, double *W, dcomplex *WORK, int *LWORK, double *RWORK,
                       int *INFO);
void zheevx_(char *JOBZ, char *RANGE, char *UPLO, int *N, dcomplex *A, int *LDA, double *VL, double *VU, int *IL, int *IU,
                        double *ABSTOL, int *M, double *W, dcomplex *Z, int *LDZ, dcomplex *WORK, int *LWORK, double *RWORK,
                        int *IWORK, int *IFAIL, int *INFO);

void zgemm_(char* TRANSA, char* TRANSB, int * M, int * N,int *K, dcomplex *alpha, 
         dcomplex *A, int *LDA, dcomplex *B, int*LDB, dcomplex *beta, dcomplex *C, int *LDC);
void zgetrf_(int *m, int *n, dcomplex *a,int *lda,int *ipvt, int *info );
void zgetri_(int *n,dcomplex *a,int *lda, int *ipvt, dcomplex *work, int *lwork, int *info);
int zcopy_(int *n, dcomplex *zx, int *incx, dcomplex *zy, int *incy);

void zgtsv_(INTEGER *n, INTEGER *nrhs, dcomplex *dl, 
	    dcomplex *d__, dcomplex *du, dcomplex *b, INTEGER *ldb,
	    INTEGER *info);

void dpotrf_(char *uplo, INTEGER *n, double *a, INTEGER *lda, INTEGER *info);
void dpotri_(char *uplo, INTEGER *n, double *a, INTEGER *lda, INTEGER *info);
void dggevx_(char *balanc, char *jobvl, char *jobvr, char *sense,
             INTEGER *n, double *a, INTEGER *lda, double *b, 
	     INTEGER *ldb, double *alphar, double *alphai, double *beta,
             double *vl, INTEGER *ldvl, double *vr, INTEGER *ldvr, 
	     INTEGER *ilo, INTEGER *ihi, double *lscale, double *rscale, 
	     double *abnrm, double *bbnrm, double *rconde, double *rcondv,
             double *work, INTEGER *lwork, INTEGER *iwork, logical *bwork,
             INTEGER *info);

void dsytrd_(char *uplo, INTEGER *n, double *a, INTEGER *lda, double *d__, double *e, 
             double *tau, double *work, INTEGER *lwork, INTEGER *info);
