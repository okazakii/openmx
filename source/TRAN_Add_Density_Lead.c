/**********************************************************************
  TRAN_Add_Density_Lead.c:

  TRAN_Add_Density_Lead.c is a subroutine to correct charge density 
  near the boundary region of the extended system.
  The charge density from that of electrodes is added to the regions
  [0:TRAN_grid_bound[0]] and [TRAN_grid_bound[1]:Ngrid1-1].

  Log of TRAN_Add_Density_Lead.c:

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


void TRAN_Add_Density_Lead(
            MPI_Comm comm1,
            int SpinP_switch,
            int Ngrid1,
            int Ngrid2,
            int Ngrid3,
            int Num_Cells0, 
            int *My_Cell0, 
            int *My_Cell1,
            double **Density_Grid)

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
    printf("<TRAN_Add_Density_Lead>\n");
  }

  /* left lead */

  side = 0;
  l1[0] = 0;
  l1[1] = TRAN_grid_bound[0]; 

  for (spin=0; spin<=SpinP_switch; spin++){

    for (i=0; i<Num_Cells0; i++) {

      ie = My_Cell1[i]; 

      if ( l1[0]<=ie && ie<=l1[1] ) {

	for (j=0; j<Ngrid2; j++) {

	  for (k=0; k<Ngrid3; k++) {
	    Density_Grid[spin][ grid_ref(i,j,k) ] += ElectrodeDensity_Grid[side][spin][ grid_e_ref(ie,j,k) ];
	  }

          if (SpinP_switch==0){
	    for (k=0; k<Ngrid3; k++) {
	      Density_Grid[1][ grid_ref(i,j,k) ] += ElectrodeDensity_Grid[side][0][ grid_e_ref(ie,j,k) ];
  	    }
          }
	}
      }
    }

  }  /* spin */

  /* right lead */

  side = 1;
  l1[0] = TRAN_grid_bound[1];
  l1[1] = Ngrid1-1;
  
  for (spin=0; spin<=SpinP_switch; spin++){

    for (i=0; i<Num_Cells0; i++) {

      ie = My_Cell1[i];

      if ( l1[0]<=ie && ie<=l1[1] ) {
	for (j=0; j<Ngrid2; j++) {
	  for (k=0; k<Ngrid3; k++) {
	    Density_Grid[spin][ grid_ref(i,j,k) ] += ElectrodeDensity_Grid[side][spin][ grid_e_ref(ie,j,k) ];
	  }

          if (SpinP_switch==0){
	    for (k=0; k<Ngrid3; k++) {
	      Density_Grid[1][ grid_ref(i,j,k) ] += ElectrodeDensity_Grid[side][0][ grid_e_ref(ie,j,k) ];
  	    }
	  }

	}
      }
    }
    
  }  /* spin */

}
