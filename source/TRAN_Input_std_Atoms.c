/**********************************************************************
  TRAN_Input_std_Atoms.c:

  TRAN_Input_std_Atoms.c is a subroutine to read the input data.

  Log of TRAN_Input_std_Atoms.c:

     24/July/2008  Released by H.Kino and T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "Inputtools.h"
#include "openmx_common.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif

#include "tran_prototypes.h"
#include "tran_variables.h"

#define MAXBUF 1024

int OrbPol2int(char OrbPol[YOUSO10]);


void TRAN_Input_std_Atoms(  MPI_Comm comm1, int Solver )
{
  int po=0;
  int idim=1;
  int myid;

  FILE *fp;
  char *s_vec[20];
  int i_vec[20];
  double r_vec[20];

  int i,j,spe,spe_e; 
  char Species[YOUSO10];
  double GO_XR,GO_YR,GO_ZR;
  double GO_XL,GO_YL,GO_ZL;
  char OrbPol[YOUSO10];
  char buf[MAXBUF];

  if (Solver!=4) return; 

  MPI_Comm_rank(comm1,&myid);

  /* center */
  input_int("Atoms.Number",&Catomnum,0);
  if (Catomnum<=0){

    if (myid==Host_ID){
      printf("Atoms.Number may be wrong.\n");
    }
    MPI_Finalize();
    exit(1);
  }

  /* left */
  input_int("LeftLeadAtoms.Number",&Latomnum,0);
  if (Latomnum<=0){

    if (myid==Host_ID){
      printf("LeftLeadAtoms.Number may be wrong.\n");
    }
    MPI_Finalize();
    exit(1);
  }

  /* right */
  input_int("RightLeadAtoms.Number",&Ratomnum,0);
  if (Ratomnum<=0){

    if (myid==Host_ID){
      printf("RightLeadAtoms.Number may be wrong.\n");
    }
    MPI_Finalize();
    exit(1);
  }
    
  atomnum = Catomnum + Latomnum + Ratomnum;
  List_YOUSO[1] = atomnum + 1;

  /* memory allocation */

  Allocate_Arrays(1);

  /* memory allocation for TRAN_* */
  TRAN_Allocate_Atoms(atomnum);

  s_vec[0]="Ang";  s_vec[1]="AU";
  i_vec[0]= 0;     i_vec[1]=1;
  input_string2int("Atoms.SpeciesAndCoordinates.Unit",
		   &coordinates_unit,2,s_vec,i_vec);

  /* left lead */

  if (fp=input_find("<LeftLeadAtoms.SpeciesAndCoordinates") ) {

    for (i=1; i<=Latomnum; i++){

      fgets(buf,MAXBUF,fp);

      sscanf(buf,"%i %s %lf %lf %lf %lf %lf %s",&j,Species,
	     &Gxyz[i][1],&Gxyz[i][2],&Gxyz[i][3],
	     &InitN_USpin[i],&InitN_DSpin[i], OrbPol);

      WhatSpecies[i] = Species2int(Species);
      TRAN_region[i]=2;
      TRAN_Original_Id[i]=j;

      if (Hub_U_switch==1) OrbPol_flag[i] = OrbPol2int(OrbPol);

      /* check the consistency of basis set */

      spe   = WhatSpecies[i];
      spe_e = WhatSpecies_e[0][i];

      if (i!=j){

        if (myid==Host_ID){
	  printf("Error of sequential number %i in <LeftLeadAtoms.SpeciesAndCoordinates\n",j);
        }
        MPI_Finalize();
        exit(1);
      }

      if (2<=level_stdout && myid==Host_ID){

	printf("<Input_std> L_AN=%2d T_AN=%2d WhatSpecies=%2d USpin=%7.4f DSpin=%7.4f\n",
	       i,i,
	       WhatSpecies[i],
	       InitN_USpin[i],
	       InitN_DSpin[i]);
      }
    }

    ungetc('\n',fp);

    if (!input_last("LeftLeadAtoms.SpeciesAndCoordinates>")) {
      /* format error */

      if (myid==Host_ID){
        printf("Format error for LeftLeadAtoms.SpeciesAndCoordinates\n");
      }
      MPI_Finalize();
      exit(1);
    }
  }

  /* center */
  if (fp=input_find("<Atoms.SpeciesAndCoordinates") ) {
    for (i=1; i<=Catomnum; i++){

      fgets(buf,MAXBUF,fp);

      sscanf(buf,"%i %s %lf %lf %lf %lf %lf %s",
             &j,Species,
	     &Gxyz[Latomnum+i][1],
             &Gxyz[Latomnum+i][2],
             &Gxyz[Latomnum+i][3],
	     &InitN_USpin[Latomnum+i],
             &InitN_DSpin[Latomnum+i], 
             OrbPol);

      WhatSpecies[Latomnum+i] = Species2int(Species);
      TRAN_region[Latomnum+i]= 1;
      TRAN_Original_Id[Latomnum+i]= j;

      if (Hub_U_switch==1) OrbPol_flag[Latomnum+i] = OrbPol2int(OrbPol);

      if (i!=j){

        if (myid==Host_ID){
	  printf("Error of sequential number %i in <Atoms.SpeciesAndCoordinates\n",j);
        }
        MPI_Finalize();
        exit(1);
      }

      if (2<=level_stdout && myid==Host_ID){
	printf("<Input_std>  ct_AN=%2d WhatSpecies=%2d USpin=%7.4f DSpin=%7.4f\n",
	       Latomnum+i,WhatSpecies[Latomnum+i],InitN_USpin[Latomnum+i],InitN_DSpin[Latomnum+i]);
      }
    }

    ungetc('\n',fp);

    if (!input_last("Atoms.SpeciesAndCoordinates>")) {
      /* format error */
      if (myid==Host_ID){
        printf("Format error for Atoms.SpeciesAndCoordinates\n");
      }
      MPI_Finalize();
      exit(1);
 
    }
  }

  /* right */

  if (fp=input_find("<RightLeadAtoms.SpeciesAndCoordinates") ) {

    for (i=1; i<=Ratomnum; i++){

      fgets(buf,MAXBUF,fp);

      sscanf(buf,"%i %s %lf %lf %lf %lf %lf %s",
             &j,Species,
	     &Gxyz[Catomnum+Latomnum+i][1],
	     &Gxyz[Catomnum+Latomnum+i][2],
	     &Gxyz[Catomnum+Latomnum+i][3],
	     &InitN_USpin[Catomnum+Latomnum+i],
	     &InitN_DSpin[Catomnum+Latomnum+i], OrbPol);

      WhatSpecies[Catomnum+Latomnum+i] = Species2int(Species);

      TRAN_region[Catomnum+Latomnum+i]= 3;
      TRAN_Original_Id[Catomnum+Latomnum+i]= j;

      if (Hub_U_switch==1) OrbPol_flag[Catomnum+Latomnum+i] = OrbPol2int(OrbPol);

      if (i!=j){
        if (myid==Host_ID){
 	  printf("Error of sequential number %i in <RightLeadAtoms.SpeciesAndCoordinates\n",j);
        }
        MPI_Finalize();
        exit(1);
      }

      if (2<=level_stdout && myid==Host_ID){
	printf("<Input_std> R_AN=%2d T_AN=%2d WhatSpecies=%2d USpin=%7.4f DSpin=%7.4f\n",
	       i,Catomnum+Latomnum+i,
	       WhatSpecies[Catomnum+Latomnum+i],
	       InitN_USpin[Catomnum+Latomnum+i],
	       InitN_DSpin[Catomnum+Latomnum+i]);
      }
    }

    ungetc('\n',fp);

    if (!input_last("RightLeadAtoms.SpeciesAndCoordinates>")) {
      /* format error */

      if (myid==Host_ID){
        printf("Format error for RightLeadAtoms.SpeciesAndCoordinates\n");
      }
      MPI_Finalize();
      exit(1);

    }
  }

  if (coordinates_unit==0){
    for (i=1; i<=atomnum; i++){
      Gxyz[i][1] = Gxyz[i][1]/BohrR;
      Gxyz[i][2] = Gxyz[i][2]/BohrR;
      Gxyz[i][3] = Gxyz[i][3]/BohrR;
    }
  }

  /* compare the coordinates with those used for the band calculation of the left lead */

  for (i=1; i<=Latomnum; i++){

    if (    1.0e-12<fabs(Gxyz_e[0][i][1]-Gxyz[i][1]+Gxyz[1][1]-Gxyz_e[0][1][1])
	 || 1.0e-12<fabs(Gxyz_e[0][i][2]-Gxyz[i][2]+Gxyz[1][2]-Gxyz_e[0][1][2])
	 || 1.0e-12<fabs(Gxyz_e[0][i][3]-Gxyz[i][3]+Gxyz[1][3]-Gxyz_e[0][1][3]) ){  

      if (myid==Host_ID){
	printf("\n\n");
	printf("The LEFT lead cannot be superposed on the original cell even after the translation.\n");
	printf("Check your atomic coordinates of the LEFT lead.\n\n");
      }

      MPI_Finalize();
      exit(1);
    }
  }

  /* compare the coordinates with those used for the band calculation of the right lead */

  for (i=1; i<=Ratomnum; i++){

    if (    1.0e-12<fabs((Gxyz_e[1][i][1]+(Gxyz[Catomnum+Latomnum+1][1]-Gxyz_e[1][1][1]))-Gxyz[Catomnum+Latomnum+i][1])
         || 1.0e-12<fabs((Gxyz_e[1][i][2]+(Gxyz[Catomnum+Latomnum+1][2]-Gxyz_e[1][1][2]))-Gxyz[Catomnum+Latomnum+i][2])
         || 1.0e-12<fabs((Gxyz_e[1][i][3]+(Gxyz[Catomnum+Latomnum+1][3]-Gxyz_e[1][1][3]))-Gxyz[Catomnum+Latomnum+i][3]) ){  

      if (myid==Host_ID){
	printf("\n\n");
	printf("The RIGHT lead cannot be superposed on the original cell even after the translation.\n");
	printf("Check your atomic coordinates of the RIGHT lead.\n\n");
      }

      MPI_Finalize();
      exit(1);
    }

  }

  /****************************************************
                          Unit cell
  ****************************************************/

  for (i=1; i<=3; i++){

    Left_tv[i][1]  = tv_e[0][i][1];
    Left_tv[i][2]  = tv_e[0][i][2];
    Left_tv[i][3]  = tv_e[0][i][3];

    Right_tv[i][1] = tv_e[1][i][1];
    Right_tv[i][2] = tv_e[1][i][2];
    Right_tv[i][3] = tv_e[1][i][3];
  }

  for (i=2; i<=3; i++){
    tv[i][1]  = tv_e[0][i][1];
    tv[i][2]  = tv_e[0][i][2];
    tv[i][3]  = tv_e[0][i][3];
  }

  /*****************************************************
  set the grid origin and the unit vectors for the NEGF
  calculations so that the boundaries between leads and 
  the extended central region match with those in the 
  band calculations for the leads. 
  *****************************************************/

  GO_XL = Gxyz[1][1] - Gxyz_e[0][1][1] + Grid_Origin_e[0][1];
  GO_YL = Gxyz[1][2] - Gxyz_e[0][1][2] + Grid_Origin_e[0][2];
  GO_ZL = Gxyz[1][3] - Gxyz_e[0][1][3] + Grid_Origin_e[0][3];

  GO_XR = Gxyz[Catomnum+Latomnum+1][1] - Gxyz_e[1][1][1] + Grid_Origin_e[1][1];
  GO_YR = Gxyz[Catomnum+Latomnum+1][2] - Gxyz_e[1][1][2] + Grid_Origin_e[1][2];
  GO_ZR = Gxyz[Catomnum+Latomnum+1][3] - Gxyz_e[1][1][3] + Grid_Origin_e[1][3];

  /* large cell = left lead + central region + right lead */    

  tv[idim][1] = GO_XR + Right_tv[idim][1] - GO_XL;
  tv[idim][2] = GO_YR + Right_tv[idim][2] - GO_YL;
  tv[idim][3] = GO_ZR + Right_tv[idim][3] - GO_ZL;

  scf_fixed_origin[0] = GO_XL;
  scf_fixed_origin[1] = GO_YL;
  scf_fixed_origin[2] = GO_ZL;

}
