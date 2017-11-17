/**********************************************************************
  TRAN_Add_ADensity_Lead.c:

  TRAN_Add_ADensity_Lead.c is a subroutine to correct atomic charge
  density near the boundary region of the extended system
  The super position of atomic charge density from that of electrodes 
  is added to the regions [0:TRAN_grid_bound[0]] and 
  [TRAN_grid_bound[1]:Ngrid1-1].

  Log of TRAN_Add_ADensity_Lead.c:

     24/July/2008  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif

#include "tran_variables.h"
#include "tran_prototypes.h"


void TRAN_Add_ADensity_Lead(
            MPI_Comm comm1,
            int SpinP_switch,
            int Ngrid1,
            int Ngrid2,
            int Ngrid3,
            int Num_Cells0, 
            int *My_Cell0, 
            int *My_Cell1,
            double *ADensity_Grid)

#define grid_ref(i,j,k)    ( (i)*Ngrid2*Ngrid3+(j)*Ngrid3+(k) )
#define grid_e_ref(i,j,k)  ( ((i)-l1[0]) *Ngrid2*Ngrid3+(j)*Ngrid3+(k) )

{
  int side,l1[2];
  int i,j,k;
  int spin;
  int ie;
  int myid;

  MPI_Comm_rank(comm1,&myid);

  if (myid==Host_ID){
    printf("<TRAN_Add_ADensity_Lead>\n");
  }

  /* left lead */

  side = 0;
  l1[0] = 0;
  l1[1] = TRAN_grid_bound[0]; 

  for (i=0; i<Num_Cells0; i++) {

    ie = My_Cell1[i]; 

    if ( l1[0]<=ie && ie<=l1[1] ) {

      for (j=0; j<Ngrid2; j++) {
	for (k=0; k<Ngrid3; k++) {
	  ADensity_Grid[ grid_ref(i,j,k) ] += ElectrodeADensity_Grid[side][ grid_e_ref(ie,j,k) ];
	}
      }
    }
  }

  /* right lead */

  side = 1;
  l1[0] = TRAN_grid_bound[1];
  l1[1] = Ngrid1-1;
  
  for (i=0; i<Num_Cells0; i++) {

    ie = My_Cell1[i];

    if ( l1[0]<=ie && ie<=l1[1] ) {
      for (j=0; j<Ngrid2; j++) {
	for (k=0; k<Ngrid3; k++) {
	  ADensity_Grid[ grid_ref(i,j,k) ] += ElectrodeADensity_Grid[side][ grid_e_ref(ie,j,k) ];
	}
      }
    }
  }
    
}
