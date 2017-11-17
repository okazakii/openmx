/*----------------------------------------------------------------------
  demo_sub.h
----------------------------------------------------------------------*/
#include "eri.h"

ERI_Orbital_t* Orbital_New_GTO(double rmax, int ngrid, int l, int m, double a, double norm);
ERI_Orbital_t* Orbital_New_STO(double rmax, int ngrid, int l, int m, double a);
ERI_Orbital_t* Orbital_New_PAO(const char *path, int l, int m, int node);

ERI_Orbital_t* Orbital_Duplicate(const ERI_Orbital_t *p);
void Orbital_Free(ERI_Orbital_t *p);


typedef struct {
  ERI_Orbital_t *o1;
  ERI_Orbital_t *o2;
  ERI_Orbital_t *o3;
  ERI_Orbital_t *o4;
} Quartet_t;

Quartet_t* Quartet_GTO_ssss(int ngrid, double rmax);
Quartet_t* Quartet_GTO_pssp(int ngrid, double rmax);
Quartet_t* Quartet_GTO_dddd(int ngrid, double rmax);

Quartet_t* Quartet_Duplicate(const Quartet_t *p);
void Quartet_Free(Quartet_t* );


