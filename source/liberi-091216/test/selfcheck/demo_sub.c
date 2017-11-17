/*----------------------------------------------------------------------
  demo_sub.c
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "demo_sub.h"

 
#define BUFFER_SIZE 4096



ERI_Orbital_t* Orbital_New_GTO(double rmax, int ngrid, int l, int m, double a, double norm)
{
  int i, k, kmax;
  double dr, r, pi;
  ERI_Orbital_t *p;

  pi = acos(-1.0);

#if 0
  /* normalization factor */
  if (0 == l%2) {
    norm = sqrt(a/pi);
    kmax = l/2;
    for (k=0; k<kmax; k++) { norm *= 2.0*a/(double)(2*k+1); }
  } else {
    norm = 2.0*a;
    kmax = (l+1)/2;
    for (k=0; k<kmax; k++) { norm *= a/(double)(k+1); }
  }
  norm = sqrt(1.5); //0.5/sqrt(pi);
#endif

  /* memory allocation */ 
  p = (ERI_Orbital_t*)malloc(sizeof(ERI_Orbital_t));
  if (p) {
    p->xr = (double*)malloc(sizeof(double)*ngrid);
    p->fr = (double*)malloc(sizeof(double)*ngrid);
  }
  if (NULL==p || NULL==p->xr || NULL==p->fr) {
    Orbital_Free(p);
    return NULL;
  }

  dr = rmax/(double)ngrid; 
  for (i=0; i<ngrid; i++) {
    r = dr*(double)i;
    p->xr[i] = r;
    p->fr[i] = norm*pow(r, (double)l)*exp(-a*r*r);
  }
   
  p->ngrid = ngrid;
  p->l     = l;
  p->m     = m;

  return p;
}




ERI_Orbital_t* Orbital_New_STO(double rmax, int ngrid, int l, int m, double a)
{
  int i, k;
  double dr, r, norm;
  ERI_Orbital_t *p;
 
  /* normalization factor */
  norm = a;
  for (k=0; k<l; k++) { norm *= a/(double)(k+1); }

  /* memory allocation */ 
  p = (ERI_Orbital_t*)malloc(sizeof(ERI_Orbital_t));
  if (p) {
    p->xr = (double*)malloc(sizeof(double)*ngrid);
    p->fr = (double*)malloc(sizeof(double)*ngrid);
  }
  if (NULL==p || NULL==p->xr || NULL==p->fr) {
    Orbital_Free(p);
    return NULL;
  }

  dr = rmax/(double)ngrid;
  for (i=0; i<ngrid; i++) {
    r = dr*(double)i;
    p->xr[i] = r;
    p->fr[i] = norm*pow(r, (double)l)*exp(-a*r);
  }
  
  p->ngrid = ngrid;
  p->l     = l;
  p->m     = m;
  
  return p;
}




/*----------------------------------------------------------------------
  Orbital_New_PAO
----------------------------------------------------------------------*/
ERI_Orbital_t* Orbital_New_PAO(
  const char *path,
  int         l, 
  int         m,
  int         node
)
{
  int i, j, n, stat, ngrid, lmax, mul;
  double rcut;
  char buffer[BUFFER_SIZE], headline[256], pao_name[256];

  ERI_Orbital_t *p;
  FILE *fp;


  /*-------------------------------------------------- open .pao file */
  fp = fopen(path, "r");
  if (fp == NULL) { return NULL; }


  /*------------------------------------------------ read information */
  stat = 1;
  while (1 == stat) {
    n = fscanf(fp, "%s", buffer);

    /* name */
    if (strcmp(buffer, "System.Name")==0) {
      n = fscanf(fp, "%s", buffer);
      if (EOF != n) strncpy(pao_name, buffer, 256);
    }
    /* number of mesh points */
    if (strcmp(buffer, "grid.num.output")==0) {
      n = fscanf(fp, "%s", buffer);
      if (EOF != n) ngrid = atoi(buffer);
    }
    /* max value of l */
    if (strcmp(buffer, "maxL.pao")==0) {
      n = fscanf(fp, "%s", buffer);
      if (EOF != n) lmax = atoi(buffer);
    }
    /* multiplicity */
    if (strcmp(buffer, "num.pao")==0) {
      n = fscanf(fp, "%s", buffer);
      if (EOF != n) mul = atoi(buffer);
    }
    /* cutoff radius */
    if (strcmp(buffer, "radial.cutoff.pao")==0) {
      n = fscanf(fp, "%s", buffer);
      if (EOF != n) rcut = atof(buffer);
    }
    /* head of the first data set */
    if (strcmp(buffer, "<pseudo.atomic.orbitals.L=0")==0) {
      stat = 0;
    } 
    if (n == EOF) stat = -1;
  }
   
  if (0 != stat) {
    fprintf(stderr, "*** FAILED TO READ PAO FILE.\n");
    fclose(fp);
    return NULL;
  }


  /*------------------------------------------------- allocate memoty */
  p = (ERI_Orbital_t*)malloc(sizeof(ERI_Orbital_t));
  if (p) {
    p->fr    = (double*)malloc(sizeof(double)*ngrid);
    p->xr    = (double*)malloc(sizeof(double)*ngrid);
  }
  if (NULL == p || NULL == p->fr || NULL == p->xr) {
    fclose(fp);
    Orbital_Free(p);
    return NULL;
  }

  p->l     = l;
  p->m     = m;
  p->ngrid = ngrid;

 
  /*--------------------------------------------------- read PAO data */
  while (1) 
  {
    if (l>lmax || node>=mul) { stat = -1; break; }

    /* find the data set for given l value */
    stat = 0;
    sprintf(headline, "<pseudo.atomic.orbitals.L=%d", l);
    while (0 == stat && strcmp(headline, buffer)!=0) {
      n = fscanf(fp, "%s", buffer);
      if (n==EOF) stat = -1;
    }
    if (0 != stat) { break ; }
 
    for (i=0; i<ngrid; i++) {
      fscanf(fp, "%s", buffer); /* ln(r) */
      n = fscanf(fp, "%s", buffer); /* r */
      if (EOF != n) p->xr[i] = atof(buffer);
      for (j=0; j<mul; j++) {
        n = fscanf(fp, "%s", buffer);
        if (j==node && EOF != n) p->fr[i]=atof(buffer);
      }
      if (n==EOF) { break; }
    }

    stat = 0;
    break;
  }
 
  fclose(fp);
  
  if (0 != stat) {
    fprintf(stderr, "*** FAILED TO READ PAO FILE\n");
    Orbital_Free(p);
    return NULL;
  }
  
  return p;
}


ERI_Orbital_t* Orbital_Duplicate(const ERI_Orbital_t* p)
{
  int i;

  ERI_Orbital_t *new;
  int ngrid = p->ngrid;

  /* memory allocation */ 
  new = (ERI_Orbital_t*)malloc(sizeof(ERI_Orbital_t));
  if (new) {
    new->xr = (double*)malloc(sizeof(double)*ngrid);
    new->fr = (double*)malloc(sizeof(double)*ngrid);
  }
  if (NULL==new || NULL==new->xr || NULL==new->fr) {
    Orbital_Free(new);
    return NULL;
  }

  for (i=0; i<ngrid; i++) {
    new->xr[i] = p->xr[i];
    new->fr[i] = p->fr[i];
  }
   
  new->ngrid = p->ngrid;
  new->l     = p->l;
  new->m     = p->m;
  new->c[0]  = p->c[0];
  new->c[1]  = p->c[1];
  new->c[2]  = p->c[2];
  
  return new;
}


void Orbital_Free(ERI_Orbital_t *p)
{
  if (p) {
    if (p->xr) { free(p->xr); }
    if (p->fr) { free(p->fr); }
    free(p);
  }
}


/*----------------------------------------------------------------------
----------------------------------------------------------------------*/
Quartet_t* Quartet_GTO_ssss(int ngrid, double rmax)
{
  Quartet_t *p;

  p = (Quartet_t*)malloc(sizeof(Quartet_t));
  if (p) {
    p->o1 = Orbital_New_GTO(rmax, ngrid, 0, 0, 1.0, 1.0);    
    p->o2 = Orbital_New_GTO(rmax, ngrid, 0, 0, 1.0, 1.0);    
    p->o3 = Orbital_New_GTO(rmax, ngrid, 0, 0, 1.0, 1.0);    
    p->o4 = Orbital_New_GTO(rmax, ngrid, 0, 0, 1.0, 1.0);    
  }
  if (NULL==p || NULL == p->o1 || NULL == p->o2 
              || NULL == p->o3 || NULL == p->o4) {
    Quartet_Free(p);
    return NULL;
  }
 
  p->o1->c[0] =  0.5; p->o1->c[1] =  0.5; p->o1->c[2] = 0.0;
  p->o2->c[0] =  0.5; p->o2->c[1] = -0.5; p->o2->c[2] = 0.0;
  p->o3->c[0] = -0.5; p->o3->c[1] =  0.5; p->o3->c[2] = 0.0;
  p->o4->c[0] = -0.5; p->o4->c[1] = -0.5; p->o4->c[2] = 0.0;

  return p;
}


Quartet_t* Quartet_GTO_pssp(int ngrid, double rmax)
{
  Quartet_t *p;

  p = (Quartet_t*)malloc(sizeof(Quartet_t));
  if (p) {
    p->o1 = Orbital_New_GTO(rmax, ngrid, 1, -1, 1.0, 1.0);    
    p->o2 = Orbital_New_GTO(rmax, ngrid, 0,  0, 1.0, 1.0);    
    p->o3 = Orbital_New_GTO(rmax, ngrid, 0,  0, 1.0, 1.0);    
    p->o4 = Orbital_New_GTO(rmax, ngrid, 1, +1, 1.0, 1.0);    
  }
  if (NULL==p || NULL == p->o1 || NULL == p->o2 
              || NULL == p->o3 || NULL == p->o4) {
    Quartet_Free(p);
    return NULL;
  }
 
  p->o1->c[0] =  0.5; p->o1->c[1] =  0.5; p->o1->c[2] = 0.0;
  p->o2->c[0] =  0.5; p->o2->c[1] = -0.5; p->o2->c[2] = 0.0;
  p->o3->c[0] = -0.5; p->o3->c[1] =  0.5; p->o3->c[2] = 0.0;
  p->o4->c[0] = -0.5; p->o4->c[1] = -0.5; p->o4->c[2] = 0.0;

  return p;
}


Quartet_t* Quartet_GTO_dddd(int ngrid, double rmax)
{
  Quartet_t *p;

  p = (Quartet_t*)malloc(sizeof(Quartet_t));
  if (p) {
    p->o1 = Orbital_New_GTO(rmax, ngrid, 2, -2, 0.5, 1.0);    
    p->o2 = Orbital_New_GTO(rmax, ngrid, 2, -1, 0.5, 1.0);    
    p->o3 = Orbital_New_GTO(rmax, ngrid, 2, +1, 0.5, 1.0);    
    p->o4 = Orbital_New_GTO(rmax, ngrid, 2, +2, 0.5, 1.0);    
  }
  if (NULL==p || NULL == p->o1 || NULL == p->o2 
              || NULL == p->o3 || NULL == p->o4) {
    Quartet_Free(p);
    return NULL;
  }
 
  p->o1->c[0] =  0.5; p->o1->c[1] =  0.5; p->o1->c[2] = 0.0;
  p->o2->c[0] =  0.5; p->o2->c[1] = -0.5; p->o2->c[2] = 0.0;
  p->o3->c[0] = -0.5; p->o3->c[1] =  0.5; p->o3->c[2] = 0.0;
  p->o4->c[0] = -0.5; p->o4->c[1] = -0.5; p->o4->c[2] = 0.0;

  return p;
}


#if 0
Quartet_t* Quartet_PAO(void)
{
  Quartet_t *p;

  p = (Quartet_t*)malloc(sizeof(Quartet_t));
  if (p) {
    p->o1 = Orbital_New_PAO("./pao/H5.0.pao",  0, 0, 0);
    p->o2 = Orbital_New_PAO("./pao/C4.0.pao",  1, 0, 1);
    p->o3 = Orbital_New_PAO("./pao/O5.0.pao",  2, 0, 2);
    p->o4 = Orbital_New_PAO("./pao/Fe5.5.pao", 3, 0, 3);
  }
  if (NULL==p || NULL==p->o1 || NULL==p->o2 || NULL==p->o3 || NULL==p->o4) {
    Quartet_Free(p);
    return NULL;
  }

  p->o1->c[0] = 0.0000; p->o1->c[1] = 0.0000; p->o1->c[2] = 0.0000;
  p->o2->c[0] = 0.0000; p->o2->c[1] = 0.5000; p->o2->c[2] = 0.0000;
  p->p3->c[0] = 0.0000; p->o3->c[1] = 0.0000; p->o3->c[2] = 1.5000;
  p->o4->c[0] = 1.1549; p->o4->c[1] = 0.7500; p->o4->c[2] = 2.0865;
  
#if 0

    printf("BASIS FUNCTIONS:\n");
    printf("  1 : FILE=%10s   NODE=%d  L=%d  M=%d\n", pao_file[0], pao_l[0], pao_m[0], pao_n[0]);
    printf("  2 : FILE=%10s   NODE=%d  L=%d  M=%d\n", pao_file[1], pao_l[1], pao_m[1], pao_n[1]);
    printf("  3 : FILE=%10s   NODE=%d  L=%d  M=%d\n", pao_file[2], pao_l[2], pao_m[2], pao_n[2]);
    printf("  4 : FILE=%10s   NODE=%d  L=%d  M=%d\n", pao_file[3], pao_l[3], pao_m[3], pao_n[3]);
    printf("\n");

    printf("CENTERS:\n");
    printf("  1 : (%12.8f, %12.8f, %12.8f)\n", orb1->c[0], orb1->c[1], orb1->c[2]);
    printf("  2 : (%12.8f, %12.8f, %12.8f)\n", orb2->c[0], orb2->c[1], orb2->c[2]);
    printf("  3 : (%12.8f, %12.8f, %12.8f)\n", orb3->c[0], orb3->c[1], orb3->c[2]);
    printf("  4 : (%12.8f, %12.8f, %12.8f)\n", orb4->c[0], orb4->c[1], orb4->c[2]);
    printf("\n");
#endif
  return p;
}
#endif

Quartet_t* Quartet_Duplicate(const Quartet_t *p)
{
  Quartet_t *new;

  new = (Quartet_t*)malloc(sizeof(Quartet_t));
  if (new) {
    new->o1 = Orbital_Duplicate(p->o1);
    new->o2 = Orbital_Duplicate(p->o2);
    new->o3 = Orbital_Duplicate(p->o3);
    new->o4 = Orbital_Duplicate(p->o4);
  }
  if (NULL==new || NULL==new->o1 || NULL==new->o2
                || NULL==new->o3 || NULL==new->o4 ) {
    Quartet_Free(new);
    return NULL;
  }

  return new;
}


void Quartet_Free(Quartet_t* p)
{
  if (p) {
    if (p->o1) { Orbital_Free(p->o1); }
    if (p->o2) { Orbital_Free(p->o2); }
    if (p->o3) { Orbital_Free(p->o3); }
    if (p->o4) { Orbital_Free(p->o4); }
    free(p);
  }
}


/* EOF */
