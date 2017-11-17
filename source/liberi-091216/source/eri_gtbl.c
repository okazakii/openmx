/*----------------------------------------------------------------------
  eri_gtbl.c

  Table for Gaunt coefficients (GC).

  06 Nov. 2009, M. Toyoda.
----------------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include "eri_sf.h"
#include "eri_gtbl.h"

#define THRESHOLD 1e-10


struct ERI_Gtbl_Struct {
  int     jmax;
  int    *n;
  long   *j1j2;
  double *gc;
};




static int nonzero(int m, int m1, int m2, int type)
{
  switch (type) {
  case ERI_SH_COMPLEX: return m == (m1+m2);
  case ERI_SH_REAL:    return ( m == ( m1+m2)) || ( m == ( m1-m2))
                           || ( m == (-m1+m2)) || ( m == (-m1-m2)) ;
  }

  fprintf(stderr, "***ERROR in %s (%d)\n", __FILE__, __LINE__);
  fprintf(stderr, "   undefined type (%d)\n", type);
  abort();

  return 0;
}


static double gaunt(int l, int m, int l1, int m1, int l2, int m2, int type)
{
  switch (type) {
  case ERI_SH_COMPLEX: return ERI_Gaunt(l, m, l1, m1, l2, m2); break;
  case ERI_SH_REAL:    return ERI_Gaunt_R(l, m, l1, m1, l2, m2); break;
  }
  
  fprintf(stderr, "***ERROR in %s (%d)\n", __FILE__, __LINE__);
  fprintf(stderr, "   undefined type (%d)\n", type);
  abort();

  return 0.0;
}




/*----------------------------------------------------------------------
  gtbl_nmax

  Count number of non-zero values of GC where l, l1 and l2 are less than 
  lmax. 

  Note that GC is zero unless all of the following conditions are hold:
    (1) l<l1+l2 and l1<l+l2 and l2<l+l1
    (2) m = m1 + m2
    (3) l + l1 + l2 = even
----------------------------------------------------------------------*/
static int gtbl_nmax(int lmax, int type)
{
  int l, m, j, l1, m1, j1, l2, m2, j2, lsum, n, jmax;
  double g;

  jmax = lmax*lmax;

  n = 0;

  for (j=0; j<jmax; j++) {
    l = (int)sqrt((double)j);
    m = j-l*(l+1);
    for (l1=0; l1<lmax; l1++) {
      for (l2=0; l2<lmax; l2++) {
        lsum = l+l1+l2; 
        if (l<=l1+l2 && l1<=l+l2 && l2<=l+l1 && lsum%2==0) {
          for (m1=-l1; m1<=l1; m1++) {
            for (m2=-l2; m2<=l2; m2++) {
              if (0==nonzero(m, m1, m2, type)) { continue; }
              g = gaunt(l, m, l1, m1, l2, m2, type);
              if (fabs(g)>=THRESHOLD) { n++; }
            }
          } /* end loop of m1 */
        } /* end if */
      } /* end loop of l2 */
    } /* end loop of l1 */
  }
 
  return n;
}




ERI_Gtbl_t* ERI_Gtbl_Init(int lmax, int type)
{
  int l, m, j, l1, m1, j1, l2, m2, j2, i, n, lsum, jmax, nmax;
  double g;
  ERI_Gtbl_t* gt;

  jmax = lmax*lmax;
  nmax = gtbl_nmax(lmax, type);

  /* allocation */
  gt = (ERI_Gtbl_t*)malloc(sizeof(ERI_Gtbl_t));
  gt->jmax = jmax;
  gt->n    = (int*)malloc(sizeof(int)*jmax);
  gt->j1j2 = (long*)malloc(sizeof(long)*nmax);
  gt->gc   = (double*)malloc(sizeof(double)*nmax);

  i = 0;
  for (j=0; j<jmax; j++) {
    l = (int)sqrt((double)j);
    m = j-l*(l+1);
    n=0;
    for (l1=0; l1<lmax; l1++) {
      for (l2=0; l2<lmax; l2++) {
        lsum = l+l1+l2; 
        if (l<=l1+l2 && l1<=l+l2 && l2<=l+l1 && lsum%2==0) {
          for (m1=-l1; m1<=l1; m1++) {
            for (m2=-l2; m2<=l2; m2++) {
              if (0==nonzero(m, m1, m2, type)) { continue; }
              g = gaunt(l, m, l1, m1, l2, m2, type);
              if (fabs(g)>=THRESHOLD) {
                j1 = l1*(l1+1)+m1;
                j2 = l2*(l2+1)+m2;
                gt->j1j2[i] = ERI_GTBL_PACK_J1J2(j1, j2);
                gt->gc[i] = g;
                n++;
                i++;
              }
            } /* end loop of m2 */
          } /* end loop of m1 */
        } /* end if */
      } /* end loop of l2 */
    } /* end loop of l1 */
    gt->n[j] = n;
  }

  return gt;
}



void ERI_Gtbl_Free(ERI_Gtbl_t *p)
{
  if (p) {
    if (p->n) { free(p->n); }
    if (p->n) { free(p->j1j2); }
    if (p->n) { free(p->gc); }
    free(p); 
  }
}


size_t ERI_Gtbl_Required_Size(int lmax, int type)
{
  int jmax, nmax;

  jmax = lmax*lmax;
  nmax = gtbl_nmax(lmax, type);

  return sizeof(ERI_Gtbl_t)
    + sizeof(int)*jmax   /* n */
    + sizeof(long)*nmax  /* j1j2 */
    + sizeof(double)*nmax /* gc */
  ;
}



int ERI_Gtbl_index(ERI_Gtbl_t* p, int j) 
{
  int i, n;

  n = 0;
  for (i=0; i<j; i++) { n += p->n[i]; }

  return n;
}

int ERI_Gtbl_jmax(ERI_Gtbl_t* p) { return p->jmax; }
const int* ERI_Gtbl_n(ERI_Gtbl_t* p) { return p->n; }
const long* ERI_Gtbl_j1j2(ERI_Gtbl_t *p) { return p->j1j2; }
const double* ERI_Gtbl_gc(ERI_Gtbl_t *p) { return p->gc; }

/* EOF */
