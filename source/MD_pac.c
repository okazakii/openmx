/**********************************************************************
  MD_pac.c:

     MD_pac.c is a subroutine to perform molecular dynamics
     simulations and geometry optimization.

  Log of MD_pac.c:

     22/Nov/2001  Released by T. Ozaki
     15/Dec/2003  DIISs are added by H. Kino
     14/May/2004  NVT_VS is added by M. Ohfuti
     25/May/2004  Modified by T. Ozaki
     14/Jul/2007  RF is added by H.M.Weng
     08/Jan/2010  NVT_VS2 is added by T. Ohwaki 

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "openmx_common.h"
#include "lapack_prototypes.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

#define Criterion_Max_Step       0.06
 
static void NoMD(int iter);
static void VerletXYZ(int iter);
static void NVT_VS(int iter);  /* added by mari */
static void NVT_VS2(int iter); /* added by Ohwaki */
static void NVT_NH(int iter); 
static void Steepest_Descent(int iter, int SD_scaling_flag);
static void GDIIS(int iter, int iter0);
static void GDIIS_BFGS(int iter, int iter0);
static void GDIIS_EF(int iter, int iter0);
static void Geometry_Opt_DIIS(int iter);
static void Geometry_Opt_DIIS_BFGS(int iter);
static void Geometry_Opt_DIIS_EF(int iter);
static void Correct_Position_In_First_Cell();
static void Geometry_Opt_RF(int iter);
static void RF(int iter, int iter0);
static void EvsLC(int iter);


double MD_pac(int iter, char *fname_input)
{
  double time0;
  double TStime,TEtime;
  int numprocs,myid;

  dtime(&TStime);
 
  /* MPI */
  if (atomnum<=MYID_MPI_COMM_WORLD) goto LAST_Proc;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (myid==Host_ID && 0<level_stdout){
    printf("\n*******************************************************\n"); 
    printf("             MD or geometry opt. at MD =%2d              \n",iter);
    printf("*******************************************************\n\n"); 
  }

  /* making of an input file with the final structure */
  if (Runtest_flag==0){
    Make_InputFile_with_FinalCoord(fname_input,iter);
  }

  switch (MD_switch) {
    case  0: NoMD(iter);                    break;
    case  1: VerletXYZ(iter);               break;
    case  2: NVT_VS(iter);                  break;  /* added by mari */
    case  3: Steepest_Descent(iter,1);      break;
    case  4: Geometry_Opt_DIIS_EF(iter);    break; 
    case  5: Geometry_Opt_DIIS_BFGS(iter);  break;  
    case  6: Geometry_Opt_RF(iter);         break;  /* added by hmweng */
    case  7: Geometry_Opt_DIIS(iter);       break;                    
    case  8:                                break;  /* not used */
    case  9: NVT_NH(iter);                  break;
    case 10:                                break;  /* not used */
    case 11: NVT_VS2(iter);                 break;  /* added by Ohwaki */
    case 12: EvsLC(iter);                   break; 
  }

  /***************************************************************
    correct atoms which are out of the first unit cell during 
    molecular dynamics simulations. The correction is not applied 
    for geometry optimization.
  ***************************************************************/

  if (   MD_switch==1 ||
	 MD_switch==2 ||
	 MD_switch==9 ||
	 MD_switch==11 ){

    Correct_Position_In_First_Cell();
  }

 LAST_Proc:

  MPI_Bcast(&MD_Opt_OK, 1, MPI_INT, Host_ID, MPI_COMM_WORLD1);

  dtime(&TEtime);
  time0 = TEtime - TStime;
  return time0;
}


void NoMD(int iter)
{
  char fileCoord[YOUSO10];
  FILE *fp_crd,*fp_SD;
  int i,j,k;
  char buf[fp_bsize];          /* setvbuf */
  int numprocs,myid,ID;
  char fileE[YOUSO10];

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  MD_Opt_OK = 1;

  if (myid==Host_ID){ 

    if (MD_Opt_OK==1 || iter==MD_IterNumber){

      strcpy(fileCoord,".crd");
      fnjoint(filepath,filename,fileCoord);
      if ((fp_crd = fopen(fileCoord,"w")) != NULL){

#ifdef xt3
        setvbuf(fp_crd,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

        fprintf(fp_crd,"\n");
        fprintf(fp_crd,"***********************************************************\n");
        fprintf(fp_crd,"***********************************************************\n");
        fprintf(fp_crd,"       xyz-coordinates (Ang) and forces (Hartree/Bohr)  \n");
        fprintf(fp_crd,"***********************************************************\n");
        fprintf(fp_crd,"***********************************************************\n\n");

        fprintf(fp_crd,"<coordinates.forces\n");
        fprintf(fp_crd,"  %i\n",atomnum);
        for (k=1; k<=atomnum; k++){
          i = WhatSpecies[k];
          j = Spe_WhatAtom[i];
          fprintf(fp_crd," %4d  %4s   %9.5f %9.5f %9.5f  %15.12f %15.12f %15.12f\n",
                  k,Atom_Symbol[j],
	          Gxyz[k][1]*BohrR,Gxyz[k][2]*BohrR,Gxyz[k][3]*BohrR,
	    	  -Gxyz[k][17],-Gxyz[k][18],-Gxyz[k][19]);
        }
        fprintf(fp_crd,"coordinates.forces>\n");
        fclose(fp_crd);
      }
      else
        printf("error(1) in MD_pac.c\n");
    }

  } /* if (myid==Host_ID) */

  /****************************************************
   write informatins to *.ene
  ****************************************************/

  if (myid==Host_ID){  
    sprintf(fileE,"%s%s.ene",filepath,filename);
    iterout_md(iter,MD_TimeStep*(iter-1),fileE);
  }

}

    


void Correct_Position_In_First_Cell()
{
  int i,Mc_AN,Gc_AN,ct_AN,k;
  int itmp,My_Correct_Position_flag;
  int numprocs,myid,ID,tag=999;
  double Cxyz[4],Frac[4];

  MPI_Status stat;
  MPI_Request request;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /* My_Correct_Position_flag  */

  My_Correct_Position_flag = 0; 

  /* loop for Mc_AN */

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];
    Cxyz[1] = Gxyz[Gc_AN][1] - Grid_Origin[1];
    Cxyz[2] = Gxyz[Gc_AN][2] - Grid_Origin[2];
    Cxyz[3] = Gxyz[Gc_AN][3] - Grid_Origin[3];
    Cell_Gxyz[Gc_AN][1] = Dot_Product(Cxyz,rtv[1])*0.5/PI;
    Cell_Gxyz[Gc_AN][2] = Dot_Product(Cxyz,rtv[2])*0.5/PI;
    Cell_Gxyz[Gc_AN][3] = Dot_Product(Cxyz,rtv[3])*0.5/PI;

    /* The fractional coordinates are kept within 0 to 1. */

    for (i=1; i<=3; i++){

      itmp = (int)Cell_Gxyz[Gc_AN][i]; 

      if (1.0<Cell_Gxyz[Gc_AN][i]){

        My_Correct_Position_flag = 1; 

	Cell_Gxyz[Gc_AN][i] = Cell_Gxyz[Gc_AN][i] - (double)itmp;

        if (0<level_stdout){ 
  	  if (i==1) printf("The fractional coordinate of a-axis for atom %d was translated within the range (0 to 1).\n",Gc_AN);
	  if (i==2) printf("The fractional coordinate of b-axis for atom %d was translated within the range (0 to 1).\n",Gc_AN);
	  if (i==3) printf("The fractional coordinate of c-axis for atom %d was translated within the range (0 to 1).\n",Gc_AN);
	}

        /* update His_Gxyz */ 

        for (k=0; k<Extrapolated_Charge_History; k++){

	  Cxyz[1] = His_Gxyz[k][(Gc_AN-1)*3+0] - Grid_Origin[1];
	  Cxyz[2] = His_Gxyz[k][(Gc_AN-1)*3+1] - Grid_Origin[2];
	  Cxyz[3] = His_Gxyz[k][(Gc_AN-1)*3+2] - Grid_Origin[3];

	  Frac[1] = Dot_Product(Cxyz,rtv[1])*0.5/PI;
	  Frac[2] = Dot_Product(Cxyz,rtv[2])*0.5/PI;
	  Frac[3] = Dot_Product(Cxyz,rtv[3])*0.5/PI;

          Frac[i] = Frac[i] - (double)itmp;
          
          His_Gxyz[k][(Gc_AN-1)*3+0] = 
                      Frac[1]*tv[1][1]
                    + Frac[2]*tv[2][1]
                    + Frac[3]*tv[3][1] + Grid_Origin[1];

          His_Gxyz[k][(Gc_AN-1)*3+1] = 
                      Frac[1]*tv[1][2]
                    + Frac[2]*tv[2][2]
                    + Frac[3]*tv[3][2] + Grid_Origin[2];

          His_Gxyz[k][(Gc_AN-1)*3+2] = 
                      Frac[1]*tv[1][3]
                    + Frac[2]*tv[2][3]
                    + Frac[3]*tv[3][3] + Grid_Origin[3];
	}

        /* update GxyzHistoryIn */ 

        for (k=0; k<(M_GDIIS_HISTORY+1); k++){

	  Cxyz[1] = GxyzHistoryIn[k][Gc_AN][1] - Grid_Origin[1];
	  Cxyz[2] = GxyzHistoryIn[k][Gc_AN][2] - Grid_Origin[2];
	  Cxyz[3] = GxyzHistoryIn[k][Gc_AN][3] - Grid_Origin[3];

	  Frac[1] = Dot_Product(Cxyz,rtv[1])*0.5/PI;
	  Frac[2] = Dot_Product(Cxyz,rtv[2])*0.5/PI;
	  Frac[3] = Dot_Product(Cxyz,rtv[3])*0.5/PI;

          Frac[i] = Frac[i] - (double)itmp;

          GxyzHistoryIn[k][Gc_AN][1] = 
                      Frac[1]*tv[1][1]
                    + Frac[2]*tv[2][1]
                    + Frac[3]*tv[3][1] + Grid_Origin[1];

          GxyzHistoryIn[k][Gc_AN][2] = 
                      Frac[1]*tv[1][2]
                    + Frac[2]*tv[2][2]
                    + Frac[3]*tv[3][2] + Grid_Origin[2];

          GxyzHistoryIn[k][Gc_AN][3] = 
                      Frac[1]*tv[1][3]
                    + Frac[2]*tv[2][3]
                    + Frac[3]*tv[3][3] + Grid_Origin[3];
	}

      }
      else if (Cell_Gxyz[Gc_AN][i]<0.0){

        My_Correct_Position_flag = 1; 

	Cell_Gxyz[Gc_AN][i] = Cell_Gxyz[Gc_AN][i] + (double)(abs(itmp)+1);

        if (0<level_stdout){ 
          if (i==1) printf("The fractional coordinate of a-axis for atom %d was translated within the range (0 to 1).\n",Gc_AN);
	  if (i==2) printf("The fractional coordinate of b-axis for atom %d was translated within the range (0 to 1).\n",Gc_AN);
	  if (i==3) printf("The fractional coordinate of c-axis for atom %d was translated within the range (0 to 1).\n",Gc_AN);
	}

        /* update His_Gxyz */ 

        for (k=0; k<Extrapolated_Charge_History; k++){

	  Cxyz[1] = His_Gxyz[k][(Gc_AN-1)*3+0] - Grid_Origin[1];
	  Cxyz[2] = His_Gxyz[k][(Gc_AN-1)*3+1] - Grid_Origin[2];
	  Cxyz[3] = His_Gxyz[k][(Gc_AN-1)*3+2] - Grid_Origin[3];

	  Frac[1] = Dot_Product(Cxyz,rtv[1])*0.5/PI;
	  Frac[2] = Dot_Product(Cxyz,rtv[2])*0.5/PI;
	  Frac[3] = Dot_Product(Cxyz,rtv[3])*0.5/PI;

          Frac[i] = Frac[i] + (double)(abs(itmp)+1);
          
          His_Gxyz[k][(Gc_AN-1)*3+0] = 
                      Frac[1]*tv[1][1]
                    + Frac[2]*tv[2][1]
                    + Frac[3]*tv[3][1] + Grid_Origin[1];

          His_Gxyz[k][(Gc_AN-1)*3+1] = 
                      Frac[1]*tv[1][2]
                    + Frac[2]*tv[2][2]
                    + Frac[3]*tv[3][2] + Grid_Origin[2];

          His_Gxyz[k][(Gc_AN-1)*3+2] = 
                      Frac[1]*tv[1][3]
                    + Frac[2]*tv[2][3]
                    + Frac[3]*tv[3][3] + Grid_Origin[3];
	}

        /* update GxyzHistoryIn */ 

        for (k=0; k<(M_GDIIS_HISTORY+1); k++){

	  Cxyz[1] = GxyzHistoryIn[k][Gc_AN][1] - Grid_Origin[1];
	  Cxyz[2] = GxyzHistoryIn[k][Gc_AN][2] - Grid_Origin[2];
	  Cxyz[3] = GxyzHistoryIn[k][Gc_AN][3] - Grid_Origin[3];

	  Frac[1] = Dot_Product(Cxyz,rtv[1])*0.5/PI;
	  Frac[2] = Dot_Product(Cxyz,rtv[2])*0.5/PI;
	  Frac[3] = Dot_Product(Cxyz,rtv[3])*0.5/PI;

          Frac[i] = Frac[i] + (double)(abs(itmp)+1);

          GxyzHistoryIn[k][Gc_AN][1] = 
                      Frac[1]*tv[1][1]
                    + Frac[2]*tv[2][1]
                    + Frac[3]*tv[3][1] + Grid_Origin[1];

          GxyzHistoryIn[k][Gc_AN][2] = 
                      Frac[1]*tv[1][2]
                    + Frac[2]*tv[2][2]
                    + Frac[3]*tv[3][2] + Grid_Origin[2];

          GxyzHistoryIn[k][Gc_AN][3] = 
                      Frac[1]*tv[1][3]
                    + Frac[2]*tv[2][3]
                    + Frac[3]*tv[3][3] + Grid_Origin[3];
	}

      }
    }

    Gxyz[Gc_AN][1] =  Cell_Gxyz[Gc_AN][1]*tv[1][1]
                    + Cell_Gxyz[Gc_AN][2]*tv[2][1]
                    + Cell_Gxyz[Gc_AN][3]*tv[3][1] + Grid_Origin[1];

    Gxyz[Gc_AN][2] =  Cell_Gxyz[Gc_AN][1]*tv[1][2]
                    + Cell_Gxyz[Gc_AN][2]*tv[2][2]
                    + Cell_Gxyz[Gc_AN][3]*tv[3][2] + Grid_Origin[2];

    Gxyz[Gc_AN][3] =  Cell_Gxyz[Gc_AN][1]*tv[1][3]
                    + Cell_Gxyz[Gc_AN][2]*tv[2][3]
                    + Cell_Gxyz[Gc_AN][3]*tv[3][3] + Grid_Origin[3];
  }

  /****************
    MPI:  Gxyz
  *****************/

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    ID = G2ID[ct_AN];

    MPI_Bcast(&Gxyz[ct_AN][0], 4, MPI_DOUBLE, ID, mpi_comm_level1);

    for (k=0; k<Extrapolated_Charge_History; k++){
      MPI_Bcast(&His_Gxyz[k][(ct_AN-1)*3], 3, MPI_DOUBLE, ID, mpi_comm_level1);
    }

    for (k=0; k<(M_GDIIS_HISTORY+1); k++){
      MPI_Bcast(&GxyzHistoryIn[k][ct_AN][0], 4, MPI_DOUBLE, ID, mpi_comm_level1);
    }

  }

  /*
  for (k=0; k<(M_GDIIS_HISTORY+1); k++){
    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      printf("myid=%2d k=%2d ct_AN=%2d %15.12f %15.12f %15.12f\n",
             myid,k,ct_AN,
             GxyzHistoryIn[k][ct_AN][1],
             GxyzHistoryIn[k][ct_AN][2],
             GxyzHistoryIn[k][ct_AN][3]);fflush(stdout);
    }
  }
  */

  /* MPI:  My_Correct_Position_flag = 0; */

  MPI_Allreduce(&My_Correct_Position_flag, &Correct_Position_flag, 1, MPI_INT, MPI_MAX, mpi_comm_level1);
}


void VerletXYZ(int iter)
{
  /***********************************************************
   NVE molecular dynamics with velocity-Verlet integrator
  ***********************************************************/
  /*********************************************************** 
   1 a.u.=2.4189*10^-2 fs, 1fs=41.341105 a.u. 
   Atom weight trasformation: proton = 1836.1526 a.u 
  ***********************************************************/
  /****************************************************
    Gxyz[][1] = x-coordinate at current step
    Gxyz[][2] = y-coordinate at current step
    Gxyz[][3] = z-coordinate at current step

    Gxyz[][14] = dEtot/dx at previous step
    Gxyz[][15] = dEtot/dy at previous step
    Gxyz[][16] = dEtot/dz at previous step

    Gxyz[][17] = dEtot/dx at current step
    Gxyz[][18] = dEtot/dy at current step
    Gxyz[][19] = dEtot/dz at current step

    Gxyz[][20] = atomic mass

    Gxyz[][21] = x-coordinate at previous step
    Gxyz[][22] = y-coordinate at previous step
    Gxyz[][23] = z-coordinate at previous step

    Gxyz[][24] = x-component of velocity at current step
    Gxyz[][25] = y-component of velocity at current step
    Gxyz[][26] = z-component of velocity at current step

    Gxyz[][27] = x-component of velocity at t+dt/2
    Gxyz[][28] = y-component of velocity at t+dt/2
    Gxyz[][29] = z-component of velocity at t+dt/2

    Gxyz[][30] = hx
    Gxyz[][31] = hy
    Gxyz[][32] = hz

  ****************************************************/

  char fileE[YOUSO10];
  double dt,dt2,back,sum,My_Ukc;
  double Wscale,scaled_force;
  int Mc_AN,Gc_AN,j,k,l;
  int numprocs,myid,ID;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  MD_Opt_OK = 0;
  dt = 41.3411*MD_TimeStep;
  dt2 = dt*dt;
  Wscale = 1836.1526;

  /****************************************************
                 velocity-Verlet algorithm
  ****************************************************/

  if (iter==1){

    /****************************************************
                       Kinetic Energy 
    ****************************************************/

    My_Ukc = 0.0;

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];
      sum = 0.0;
      for (j=1; j<=3; j++){
        if (atom_Fixed_XYZ[Gc_AN][j]==0){
          sum += Gxyz[Gc_AN][j+23]*Gxyz[Gc_AN][j+23];
	}
      }
      My_Ukc = My_Ukc + 0.5*Gxyz[Gc_AN][20]*Wscale*sum;
    }

    /****************************************************
     MPI, Ukc
    ****************************************************/

    MPI_Allreduce(&My_Ukc, &Ukc, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

    /* calculation of temperature (K) */
    Temp = Ukc/(1.5*kB*(double)atomnum)*eV2Hartree;

    /****************************************************
     write informatins to *.ene
    ****************************************************/

    if (myid==Host_ID){  

      sprintf(fileE,"%s%s.ene",filepath,filename);
      iterout_md(iter,MD_TimeStep*(iter-1),fileE);
    }

    /****************************************************
      first step in velocity Verlet 
    ****************************************************/

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];

      for (j=1; j<=3; j++){

        if (atom_Fixed_XYZ[Gc_AN][j]==0){

          scaled_force = -Gxyz[Gc_AN][16+j]/(Gxyz[Gc_AN][20]*Wscale);  

          /* v( t+0.5*dt ) */
          Gxyz[Gc_AN][26+j] = Gxyz[Gc_AN][23+j] + scaled_force*0.5*dt;
          Gxyz[Gc_AN][23+j] = Gxyz[Gc_AN][26+j];

          /* r( t+dt ) */
          Gxyz[Gc_AN][20+j] = Gxyz[Gc_AN][j];
 	  Gxyz[Gc_AN][j] =  Gxyz[Gc_AN][j] + Gxyz[Gc_AN][26+j]*dt;

	}
      }

    }
  }
  else{

    /****************************************************
      second step in velocity Verlet 
    ****************************************************/

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];
      for (j=1; j<=3; j++){

        if (atom_Fixed_XYZ[Gc_AN][j]==0){
          scaled_force = -Gxyz[Gc_AN][16+j]/(Gxyz[Gc_AN][20]*Wscale);  
          Gxyz[Gc_AN][23+j] = Gxyz[Gc_AN][26+j] + scaled_force*0.5*dt;
	}
      }
    }

    /****************************************************
                       Kinetic Energy 
    ****************************************************/

    My_Ukc = 0.0;

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];
      sum = 0.0;
      for (j=1; j<=3; j++){
        if (atom_Fixed_XYZ[Gc_AN][j]==0){
          sum += Gxyz[Gc_AN][j+23]*Gxyz[Gc_AN][j+23];
	}
      }
      My_Ukc = My_Ukc + 0.5*Gxyz[Gc_AN][20]*Wscale*sum;
    }

    /****************************************************
     MPI, Ukc
    ****************************************************/

    MPI_Allreduce(&My_Ukc, &Ukc, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

    /* calculation of temperature (K) */
    Temp = Ukc/(1.5*kB*(double)atomnum)*eV2Hartree;

    /****************************************************
     write informatins to *.ene
    ****************************************************/

    if (myid==Host_ID){  

      sprintf(fileE,"%s%s.ene",filepath,filename);
      iterout_md(iter,MD_TimeStep*(iter-1),fileE);
    } 

    /****************************************************
      first step in velocity Verlet 
    ****************************************************/

    if (iter!=MD_IterNumber){

      for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
	Gc_AN = M2G[Mc_AN];
	for (j=1; j<=3; j++){

	  if (atom_Fixed_XYZ[Gc_AN][j]==0){

	    scaled_force = -Gxyz[Gc_AN][16+j]/(Gxyz[Gc_AN][20]*Wscale);  
	    /* v( t+0.5*dt ) */
	    Gxyz[Gc_AN][26+j] = Gxyz[Gc_AN][23+j] + scaled_force*0.5*dt;
	    Gxyz[Gc_AN][23+j] = Gxyz[Gc_AN][26+j];

	    /* r( t+dt ) */
	    Gxyz[Gc_AN][20+j] = Gxyz[Gc_AN][j];
	    Gxyz[Gc_AN][j] =  Gxyz[Gc_AN][j] + Gxyz[Gc_AN][26+j]*dt;
	  }

	}

      }
    }
  }

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    ID = G2ID[Gc_AN];
    MPI_Bcast(&Gxyz[Gc_AN][1],  3, MPI_DOUBLE, ID, mpi_comm_level1);
    MPI_Bcast(&Gxyz[Gc_AN][17],13, MPI_DOUBLE, ID, mpi_comm_level1);
  }
}




void Steepest_Descent(int iter, int SD_scaling_flag)
{
  /* 1au=2.4189*10^-2 fs, 1fs=41.341105 au */
  int i,j,k,l,Mc_AN,Gc_AN;
  double dt,SD_max,SD_min,SD_init,Atom_W,tmp0,scale;
  double My_Max_Force;
  char fileCoord[YOUSO10];
  char fileSD[YOUSO10];
  FILE *fp_crd,*fp_SD;
  int numprocs,myid,ID;
  double tmp1,MaxStep;
  char buf[fp_bsize];          /* setvbuf */
  char fileE[YOUSO10];

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  MD_Opt_OK = 0;

  /****************************************************
   find the maximum value of force 
  ****************************************************/

  My_Max_Force = 0.0;
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];

    tmp0 = 0.0;
    for (j=1; j<=3; j++){
      if (atom_Fixed_XYZ[Gc_AN][j]==0){
        tmp0 += Gxyz[Gc_AN][16+j]*Gxyz[Gc_AN][16+j];
      }
    }
    tmp0 = sqrt(tmp0); 
    if (My_Max_Force<tmp0) My_Max_Force = tmp0;
  }

  /****************************************************
   MPI, Max_Force
  ****************************************************/

  MPI_Allreduce(&My_Max_Force, &Max_Force, 1, MPI_DOUBLE, MPI_MAX, mpi_comm_level1);
  if (Max_Force<MD_Opt_criterion) MD_Opt_OK = 1;

  /****************************************************
   set up SD_scaling
  ****************************************************/

  dt = 41.3411*2.0;
  SD_init = dt*dt/1836.1526;
  SD_max = SD_init*10.0;   /* default 10   */
  SD_min = SD_init*0.04;   /* default 0.02 */
  Atom_W = 12.0;

  if (iter==1 || SD_scaling_flag==0){

    SD_scaling_user = Max_Force/BohrR/2.0;
    SD_scaling = SD_scaling_user/(Max_Force+1.0e-10);

    if (SD_max<SD_scaling) SD_scaling = SD_max;
    if (SD_scaling<SD_min) SD_scaling = SD_min;
  }

  else{

    if (Past_Utot[1]<Utot && iter%3==1){ 
      SD_scaling = SD_scaling/5.0;
    }
    else if (Past_Utot[1]<Past_Utot[2] && Utot<Past_Utot[1] && iter%4==1){
      SD_scaling = SD_scaling*2.5;
    }

    if (SD_max<SD_scaling) SD_scaling = SD_max;
    if (SD_scaling<SD_min) SD_scaling = SD_min;

    Past_Utot[5] = Past_Utot[4];
    Past_Utot[4] = Past_Utot[3];
    Past_Utot[3] = Past_Utot[2];
    Past_Utot[2] = Past_Utot[1];
    Past_Utot[1] = Utot;
  }
  
  if (myid==Host_ID && 0<level_stdout) printf("<Steepest_Descent>  SD_scaling=%15.12f\n",SD_scaling);

  /****************************************************
   write informatins to *.ene
  ****************************************************/

  if (myid==Host_ID){  

    sprintf(fileE,"%s%s.ene",filepath,filename);
    iterout_md(iter,MD_TimeStep*(iter-1),fileE);
  }

  /****************************************************
    move atoms
  ****************************************************/

  if (MD_Opt_OK!=1 && iter!=MD_IterNumber){

    /* avoid too large movement */
    if ( Criterion_Max_Step<(Max_Force*SD_scaling) )
      scale = Criterion_Max_Step/(Max_Force*SD_scaling);
    else 
      scale = 1.0; 

    /* update coordinates */
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

      Gc_AN = M2G[Mc_AN];

      for (j=1; j<=3; j++){
	if (atom_Fixed_XYZ[Gc_AN][j]==0){
	  Gxyz[Gc_AN][j] = Gxyz[Gc_AN][j] - scale*SD_scaling*Gxyz[Gc_AN][16+j];
	}
      }
    }
  }

  /****************************************************
   MPI, Gxyz
  ****************************************************/

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    ID = G2ID[Gc_AN];
    MPI_Bcast(&Gxyz[Gc_AN][1],  3, MPI_DOUBLE, ID, mpi_comm_level1);
    MPI_Bcast(&Gxyz[Gc_AN][17], 3, MPI_DOUBLE, ID, mpi_comm_level1);
  }

  if (myid==Host_ID){ 

    if (0<level_stdout){ 

      printf("<Steepest_Descent>  |Maximum force| (Hartree/Bohr) =%15.12f\n",
             Max_Force);
      printf("<Steepest_Descent>  Criterion       (Hartree/Bohr) =%15.12f\n",
	     MD_Opt_criterion);

      printf("\n");
      for (i=1; i<=atomnum; i++){
	printf("  atom=%3d, XYZ(ang) Fxyz(a.u.)=%9.4f %9.4f %9.4f  %9.4f %9.4f %9.4f\n",
	       i,BohrR*Gxyz[i][1],BohrR*Gxyz[i][2],BohrR*Gxyz[i][3],
	       Gxyz[i][17],Gxyz[i][18],Gxyz[i][19] ); 
      }   
    }

    strcpy(fileSD,".SD");
    fnjoint(filepath,filename,fileSD);
    if ((fp_SD = fopen(fileSD,"a")) != NULL){

#ifdef xt3
      setvbuf(fp_SD,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      if ( Criterion_Max_Step<(Max_Force*SD_scaling) )
	MaxStep = Criterion_Max_Step;
      else 
	MaxStep = SD_scaling*Max_Force;

      if (iter==1){

        fprintf(fp_SD,"\n");
        fprintf(fp_SD,"***********************************************************\n");
        fprintf(fp_SD,"***********************************************************\n");
        fprintf(fp_SD,"              History of geometry optimization             \n");
        fprintf(fp_SD,"***********************************************************\n");
        fprintf(fp_SD,"***********************************************************\n\n");

        fprintf(fp_SD,"  MD_iter   SD_scaling     |Maximum force|   Maximum step        Utot\n");
        fprintf(fp_SD,"                           (Hartree/Bohr)        (Ang)         (Hartree)\n\n");
      }
      fprintf(fp_SD,"  %3d  %15.8f  %15.8f  %15.8f  %15.8f\n",
              iter,SD_scaling,Max_Force,MaxStep*BohrR, Utot);
      fclose(fp_SD);
    }
    else
      printf("Error(7) in MD_pac.c\n");

    if (MD_Opt_OK==1 || iter==MD_IterNumber){

      strcpy(fileCoord,".crd");
      fnjoint(filepath,filename,fileCoord);
      if ((fp_crd = fopen(fileCoord,"w")) != NULL){

#ifdef xt3
        setvbuf(fp_crd,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

        fprintf(fp_crd,"\n");
        fprintf(fp_crd,"***********************************************************\n");
        fprintf(fp_crd,"***********************************************************\n");
        fprintf(fp_crd,"       xyz-coordinates (Ang) and forces (Hartree/Bohr)  \n");
        fprintf(fp_crd,"***********************************************************\n");
        fprintf(fp_crd,"***********************************************************\n\n");

        fprintf(fp_crd,"<coordinates.forces\n");
        fprintf(fp_crd,"  %i\n",atomnum);
        for (k=1; k<=atomnum; k++){
          i = WhatSpecies[k];
          j = Spe_WhatAtom[i];
          fprintf(fp_crd," %4d  %4s   %9.5f %9.5f %9.5f  %15.12f %15.12f %15.12f\n",
                  k,Atom_Symbol[j],
	          Gxyz[k][1]*BohrR,Gxyz[k][2]*BohrR,Gxyz[k][3]*BohrR,
	    	  -Gxyz[k][17],-Gxyz[k][18],-Gxyz[k][19]);
        }
        fprintf(fp_crd,"coordinates.forces>\n");
        fclose(fp_crd);
      }
      else
        printf("error(1) in MD_pac.c\n");
    }

  } /* if (myid==Host_ID) */

}




/*BEGIN hmweng*/
void Geometry_Opt_RF(int iter)
{
  int i,iatom,k,diis_iter;
  double sMD_TimeStep;
  static int local_iter=1;
  static int SD_iter=0,GDIIS_iter=0;
  static int flag=0;
  static int Every_iter;
  int everyiter,buf_iter;
  int numprocs,myid,ID;  

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  Every_iter = OptEveryDIIS;

  if (iter<M_GDIIS_HISTORY)
    diis_iter = iter;
  else   
    diis_iter = M_GDIIS_HISTORY;

  /* increament of iter */

  if (iter<OptStartDIIS){
    flag = 0;     
  }

  else if (iter==OptStartDIIS){
    flag = 1;     
    GDIIS_iter++;
  }

  else if (flag==0){
    SD_iter++; 
  }  
  else if (flag==1){
    GDIIS_iter++;
  }

  /* SD */

  if (flag==0){

    if (SD_iter==1)
      Steepest_Descent(iter,0);
    else 
      Steepest_Descent(iter,1);

    /* shift one */

    for (i=(diis_iter-2); 0<=i; i--) {
      for (iatom=1; iatom<=atomnum; iatom++) {
	for (k=1; k<=3; k++) {
	  GxyzHistoryIn[i+1][iatom][k] = GxyzHistoryIn[i][iatom][k];
	  GxyzHistoryR[i+1][iatom][k]  = GxyzHistoryR[i][iatom][k];
	}
      }
    }

    /* add GxyzHistoryIn and GxyzHisotryR */

    for (iatom=1; iatom<=atomnum; iatom++) {
      for (k=1; k<=3; k++) {
        if (atom_Fixed_XYZ[iatom][k]==0){
	  GxyzHistoryIn[0][iatom][k] = Gxyz[iatom][k];
	  GxyzHistoryR[0][iatom][k]  = Gxyz[iatom][16+k];
	}
        else{
	  GxyzHistoryIn[0][iatom][k] = Gxyz[iatom][k];
	  GxyzHistoryR[0][iatom][k]  = 0.0;
        }
      }
    }

    /* initialize local_iter */

    local_iter = 1;

  }

  /* GDIIS */

  else {

    RF(local_iter,iter);
    local_iter++;
  }

  /* check the number of iterations */

  if (Every_iter<=SD_iter){
    flag = 1;
    SD_iter = 0;
    GDIIS_iter = 0;
  }

  else if (Every_iter<=GDIIS_iter){
    flag = 0;
    SD_iter = 0;
    GDIIS_iter = 0;
  }

}
/*END hmweng */


void Geometry_Opt_DIIS(int iter)
{
  int i,iatom,k,diis_iter;
  double sMD_TimeStep;
  static int local_iter=1;
  static int SD_iter=0,GDIIS_iter=0;
  static int flag=0;
  static int Every_iter;
  int everyiter,buf_iter;
  int numprocs,myid,ID;  

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  Every_iter = OptEveryDIIS;

  if (iter<M_GDIIS_HISTORY)
    diis_iter = iter;
  else   
    diis_iter = M_GDIIS_HISTORY;

  /* increament of iter */

  if (iter<OptStartDIIS){
    flag = 0;     
  }

  else if (iter==OptStartDIIS){
    flag = 1;     
    GDIIS_iter++;
  }

  else if (flag==0){
    SD_iter++; 
  }  
  else if (flag==1){
    GDIIS_iter++;
  }

  /* SD */

  if (flag==0){

    if (SD_iter==1)
      Steepest_Descent(iter,0);
    else 
      Steepest_Descent(iter,1);

    /* shift one */

    for (i=(diis_iter-2); 0<=i; i--) {
      for (iatom=1; iatom<=atomnum; iatom++) {
	for (k=1; k<=3; k++) {
	  GxyzHistoryIn[i+1][iatom][k] = GxyzHistoryIn[i][iatom][k];
	  GxyzHistoryR[i+1][iatom][k]  = GxyzHistoryR[i][iatom][k];
	}
      }
    }

    /* add GxyzHistoryIn and GxyzHisotryR */

    for (iatom=1; iatom<=atomnum; iatom++) {
      for (k=1; k<=3; k++) {
        if (atom_Fixed_XYZ[iatom][k]==0){
	  GxyzHistoryIn[0][iatom][k] = Gxyz[iatom][k];
	  GxyzHistoryR[0][iatom][k]  = Gxyz[iatom][16+k];
	}
        else{
	  GxyzHistoryIn[0][iatom][k] = Gxyz[iatom][k];
	  GxyzHistoryR[0][iatom][k]  = 0.0;
        }
      }
    }

    /* initialize local_iter */

    local_iter = 1;

  }

  /* GDIIS */

  else {

    GDIIS(local_iter,iter);
    local_iter++;
  }

  /* check the number of iterations */

  if (Every_iter<=SD_iter){
    flag = 1;
    SD_iter = 0;
    GDIIS_iter = 0;
  }

  else if (Every_iter<=GDIIS_iter){
    flag = 0;
    SD_iter = 0;
    GDIIS_iter = 0;
  }

}




void Geometry_Opt_DIIS_BFGS(int iter)
{
  int i,iatom,k,diis_iter;
  double sMD_TimeStep;
  static int local_iter=1;
  static int SD_iter=0,GDIIS_iter=0;
  static int flag=0;
  static int Every_iter;
  int everyiter,buf_iter;
  int numprocs,myid,ID;  

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  Every_iter = OptEveryDIIS;

  if (iter<M_GDIIS_HISTORY)
    diis_iter = iter;
  else   
    diis_iter = M_GDIIS_HISTORY;

  /* increament of iter */

  if (iter<OptStartDIIS){
    flag = 0;     
  }

  else if (iter==OptStartDIIS){
    flag = 1;     
    GDIIS_iter++;
  }

  else if (flag==0){
    SD_iter++; 
  }  
  else if (flag==1){
    GDIIS_iter++;
  }

  /* SD */

  if (flag==0){

    if (SD_iter==1)
      Steepest_Descent(iter,0);
    else 
      Steepest_Descent(iter,1);

    /* shift one */

    for (i=(diis_iter-2); 0<=i; i--) {
      for (iatom=1; iatom<=atomnum; iatom++) {
	for (k=1; k<=3; k++) {
	  GxyzHistoryIn[i+1][iatom][k] = GxyzHistoryIn[i][iatom][k];
	  GxyzHistoryR[i+1][iatom][k]  = GxyzHistoryR[i][iatom][k];
	}
      }
    }

    /* add GxyzHistoryIn and GxyzHisotryR */

    for (iatom=1; iatom<=atomnum; iatom++) {
      for (k=1; k<=3; k++) {
        if (atom_Fixed_XYZ[iatom][k]==0){
	  GxyzHistoryIn[0][iatom][k] = Gxyz[iatom][k];
	  GxyzHistoryR[0][iatom][k]  = Gxyz[iatom][16+k];
	}
        else{
	  GxyzHistoryIn[0][iatom][k] = Gxyz[iatom][k];
	  GxyzHistoryR[0][iatom][k]  = 0.0;
        }
      }
    }

    /* initialize local_iter */

    local_iter = 1;

  }

  /* GDIIS */

  else {
    GDIIS_BFGS(local_iter,iter);
    local_iter++;
  }

  /* check the number of iterations */

  if (Every_iter<=SD_iter){
    flag = 1;
    SD_iter = 0;
    GDIIS_iter = 0;
  }

  else if (Every_iter<=GDIIS_iter){
    flag = 0;
    SD_iter = 0;
    GDIIS_iter = 0;
  }

}


void Geometry_Opt_DIIS_EF(int iter)
{
  int i,iatom,k,diis_iter;
  double sMD_TimeStep;
  static int local_iter=1;
  static int SD_iter=0,GDIIS_iter=0;
  static int flag=0;
  static int Every_iter;
  int everyiter,buf_iter;
  int numprocs,myid,ID;  

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  Every_iter = OptEveryDIIS;

  if (iter<M_GDIIS_HISTORY)
    diis_iter = iter;
  else   
    diis_iter = M_GDIIS_HISTORY;

  /* increament of iter */

  if (iter<OptStartDIIS){
    flag = 0;     
  }

  else if (iter==OptStartDIIS){
    flag = 1;     
    GDIIS_iter++;
  }

  else if (flag==0){
    SD_iter++; 
  }  
  else if (flag==1){
    GDIIS_iter++;
  }

  /* SD */

  if (flag==0){

    if (SD_iter==1)
      Steepest_Descent(iter,0);
    else 
      Steepest_Descent(iter,1);

    /* shift one */

    for (i=(diis_iter-2); 0<=i; i--) {
      for (iatom=1; iatom<=atomnum; iatom++) {
	for (k=1; k<=3; k++) {
	  GxyzHistoryIn[i+1][iatom][k] = GxyzHistoryIn[i][iatom][k];
	  GxyzHistoryR[i+1][iatom][k]  = GxyzHistoryR[i][iatom][k];
	}
      }
    }

    /* add GxyzHistoryIn and GxyzHisotryR */

    for (iatom=1; iatom<=atomnum; iatom++) {
      for (k=1; k<=3; k++) {
        if (atom_Fixed_XYZ[iatom][k]==0){
	  GxyzHistoryIn[0][iatom][k] = Gxyz[iatom][k];
	  GxyzHistoryR[0][iatom][k]  = Gxyz[iatom][16+k];
	}
        else{
	  GxyzHistoryIn[0][iatom][k] = Gxyz[iatom][k];
	  GxyzHistoryR[0][iatom][k]  = 0.0;
        }
      }
    }

    /* initialize local_iter */

    local_iter = 1;

  }

  /* GDIIS */

  else {

    GDIIS_EF(local_iter,iter);
    local_iter++;
  }

  /* check the number of iterations */

  if (Every_iter<=SD_iter){
    flag = 1;
    SD_iter = 0;
    GDIIS_iter = 0;
  }

  else if (Every_iter<=GDIIS_iter){
    flag = 0;
    SD_iter = 0;
    GDIIS_iter = 0;
  }

}



void GDIIS(int iter, int iter0)
{
  /* 1au=2.4189*10^-2 fs, 1fs=41.341105 au */

  char *func_name="DIIS";
  char *JOBB="L";
  double *A,*B,sumB,max_A, RR,dRi[4],dRj[4];
  double *work;
  double mixing,force_Max;
  static double sMD_TimeStep,dx_max=0.05; 
  double diff_dx,diff,Max_Step;
  int *ipiv;
  INTEGER i,j,k,iatom,N,LDA,LWORK,info;
  int diis_iter;
  char fileCoord[YOUSO10];
  char fileSD[YOUSO10];
  FILE *fp_crd,*fp_SD;
  char buf[fp_bsize];          /* setvbuf */
  char fileE[YOUSO10];

  /* variables for MPI */
  int Gc_AN;
  int numprocs,myid,ID;  

  /* MPI myid */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  /* share Gxyz */
  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    ID = G2ID[Gc_AN];
    MPI_Bcast(&Gxyz[Gc_AN][1],  3, MPI_DOUBLE, ID, mpi_comm_level1);
    MPI_Bcast(&Gxyz[Gc_AN][17], 3, MPI_DOUBLE, ID, mpi_comm_level1);
  }

  if (iter<M_GDIIS_HISTORY)
    diis_iter = iter;
  else   
    diis_iter = M_GDIIS_HISTORY;

  /* shift one */

  for (i=(diis_iter-2); 0<=i; i--) {
    for (iatom=1; iatom<=atomnum; iatom++) {
      for (k=1; k<=3; k++) {
	GxyzHistoryIn[i+1][iatom][k] = GxyzHistoryIn[i][iatom][k];
	GxyzHistoryR[i+1][iatom][k]  = GxyzHistoryR[i][iatom][k];
      }
    }
  }

  /* add GxyzHistoryIn and GxyzHisotryR */

  for (iatom=1; iatom<=atomnum; iatom++)   {

    for (k=1;k<=3;k++) {
      if (atom_Fixed_XYZ[iatom][k]==0){
	GxyzHistoryIn[0][iatom][k] = Gxyz[iatom][k];
	GxyzHistoryR[0][iatom][k]  = Gxyz[iatom][16+k];
      }
      else{
	GxyzHistoryIn[0][iatom][k] = Gxyz[iatom][k];
	GxyzHistoryR[0][iatom][k]  = 0.0;
      }
    }
  }

  if (myid!=Host_ID)   goto Last_Bcast; 

  /*********************** for myid==Host_ID **************************/

  /* allocation of arrays */

  A = (double*)malloc(sizeof(double)*(diis_iter+1)*(diis_iter+1));
  for (i=0; i<(diis_iter+1)*(diis_iter+1); i++) A[i] = 0.0;
  B = (double*)malloc(sizeof(double)*(diis_iter+1));

  /* Max of force */

  force_Max=0.0;
  for (iatom=1;iatom<=atomnum;iatom++)   {
    for (k=1;k<=3;k++) {
      if (atom_Fixed_XYZ[iatom][k]==0){
	if (force_Max< fabs(Gxyz[iatom][16+k]) ) force_Max = fabs(Gxyz[iatom][16+k]);
      }
    }
  }

  sMD_TimeStep = 0.05/(0.01*41.341105);

  if (2<=level_stdout){
    printf("<%s>  |Maximum force| (Hartree/Bohr) = %15.12f tuned_dt= %f\n",func_name,force_Max, sMD_TimeStep);
    printf("<%s>  Criterion      (Hartree/Bohr) = %15.12f\n", func_name, MD_Opt_criterion);
  }

  for (i=0; i<diis_iter; i++) {
    for (j=0; j<diis_iter; j++) {

      RR = 0.0;

      for (iatom=1; iatom<=atomnum; iatom++)   {
	for (k=1; k<=3; k++) {
	  dRi[k] = GxyzHistoryR[i][iatom][k];
	  dRj[k] = GxyzHistoryR[j][iatom][k];
	}

	RR += dRi[1]*dRj[1] + dRi[2]*dRj[2] + dRi[3]*dRj[3];
      }

      A[ i*(diis_iter+1)+j ]= RR;
    }
  }

  /* find max of A */

  max_A = 0.0;

  for (i=0;i<diis_iter;i++) {
    for (j=0;j<diis_iter;j++) {
      RR = fabs(A[i*(diis_iter+1)+j]) ;
      if (max_A< RR ) max_A = RR;
    }
  }

  max_A = 1.0/max_A;

  for (i=0; i<diis_iter; i++) {
    for (j=0; j<diis_iter; j++) {
      A[ i*(diis_iter+1)+j ] *= max_A;
    }
  }

  for (i=0; i<diis_iter; i++) {
    A[ i*(diis_iter+1)+diis_iter ] = 1.0;
    A[ diis_iter*(diis_iter+1)+i ] = 1.0;
  }

  A[diis_iter*(diis_iter+1)+diis_iter] = 0.0;

  for (i=0; i<diis_iter; i++) B[i] = 0.0;
  B[diis_iter] = 1.0;

  if (2<=level_stdout){
    printf("<%s>  DIIS matrix\n",func_name);
    for (i=0; i<(diis_iter+1); i++) {
      printf("<%s> ",func_name);
      for (j=0; j<(diis_iter+1); j++) {
        printf("%6.3f ",A[i*(diis_iter+1)+j]);
      }
      printf("\n");
    }
  }

  /* lapack routine */

  N=diis_iter+1;
  LDA=diis_iter+1;
  LWORK=diis_iter+1;
  work=(double*)malloc(sizeof(double)*LWORK);
  ipiv = (int*)malloc(sizeof(int)*(diis_iter+1));

  i = 1; 

  if (2<=level_stdout){
    printf("M_GDIIS_HISTORY=%2d diis_iter=%2d\n",M_GDIIS_HISTORY,diis_iter);
  }

  F77_NAME(dsysv,DSYSV)( JOBB, &N, &i, A, &LDA,  ipiv, B, &LDA, work, &LWORK, &info);

  if (info!=0) {
    printf("<%s> dsysv_, info=%d\n",func_name,info);
    printf("<%s> \n",func_name);
    printf("<%s> ERROR, aborting\n",func_name);
    printf("<%s> \n",func_name);

    MD_Opt_OK =1; 
    /* no change */

    goto Last_Bcast ;
  }

  if (2<=level_stdout){
    printf("<%s> diis alpha=",func_name);
    sumB = 0;
    for (i=0; i<diis_iter; i++) {
      printf("%f ",B[i]);
      sumB += B[i];
    }
    printf("%lf\n",B[diis_iter]);
  }

  if (force_Max<MD_Opt_criterion )  MD_Opt_OK = 1;

  /****************************************************
   write informatins to *.ene
  ****************************************************/

  if (myid==Host_ID){  

    sprintf(fileE,"%s%s.ene",filepath,filename);
    iterout_md(iter,MD_TimeStep*(iter-1),fileE);
  }

  if (MD_Opt_OK!=1 && iter!=MD_IterNumber){

    /* initialize */

    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1; j<=3; j++) {
	Gxyz[iatom][j] = 0.0;
      }
    }

    /* add tilde{R} */

    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1; j<=3; j++) {
	if (atom_Fixed_XYZ[iatom][j]==0){
	  for (i=0; i<diis_iter; i++) {
	    Gxyz[iatom][j] += GxyzHistoryR[i][iatom][j]*B[i];
	  }
	}
      }
    }

    sumB = 0.0;
    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1; j<=3; j++) {
	sumB += Gxyz[iatom][j]*Gxyz[iatom][j] ;
      }
    }

    sumB = sqrt(sumB)/(double)atomnum;

    if (2<=level_stdout){
      printf("<%s> |tilde{R}|=%E\n",func_name, sumB);
    }

    if      (1.0e-2<sumB)  mixing =-0.2;
    else if (1.0e-3<sumB)  mixing =-0.3;
    else if (1.0e-4<sumB)  mixing =-0.4;
    else if (1.0e-5<sumB)  mixing =-0.5;
    else if (1.0e-6<sumB)  mixing =-0.6;
    else                   mixing =-0.7;

    if (2<=level_stdout){
      printf("<%s> mixing=%15.12f\n",func_name,mixing);
    }

    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1;j<=3;j++) {
	Gxyz[iatom][j] *= mixing;
      }
    }

    /* tilde{x} */

    for (iatom=1;iatom<=atomnum;iatom++) {
      for (j=1;j<=3;j++) {
	for (i=0; i<diis_iter; i++) {
	  Gxyz[iatom][j] += GxyzHistoryIn[i][iatom][j]*B[i]; 
	}
      }
    }

    /************************************************************
        In case of a too large updating, do a modest updating
    ************************************************************/  

    Max_Step = 0.0;
    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1; j<=3; j++) {

        diff = fabs(Gxyz[iatom][j] - GxyzHistoryIn[0][iatom][j]);
        if (Max_Step<diff) Max_Step = diff;
      }
    }

    if (Criterion_Max_Step<Max_Step){

      for (iatom=1; iatom<=atomnum; iatom++) {
	for (j=1; j<=3; j++) {

	  diff = Gxyz[iatom][j] - GxyzHistoryIn[0][iatom][j];
	  Gxyz[iatom][j] = GxyzHistoryIn[0][iatom][j] + diff/Max_Step*Criterion_Max_Step;
	}
      }
    }

    /* find Max_Step */
    Max_Step = 0.0;
    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1; j<=3; j++) {

        diff = fabs(Gxyz[iatom][j] - GxyzHistoryIn[0][iatom][j]);
        if (Max_Step<diff) Max_Step = diff;
      }
    }

    if (2<=level_stdout){

      printf("<%s> diff_x= %f , dE= %f\n",func_name, Max_Step, fabs(Utot-Past_Utot[1]) );

      /* print atomic positions */
      printf("<%s> atomnum= %d\n",func_name,atomnum);
      for (i=1; i<=atomnum; i++){
	j = Spe_WhatAtom[WhatSpecies[i]];
	printf("  %3d %s XYZ(ang) Fxyz(a.u.)= %9.4f %9.4f %9.4f  %9.4f %9.4f %9.4f\n",
	       i,Atom_Symbol[j],
	       BohrR*Gxyz[i][1],BohrR*Gxyz[i][2],BohrR*Gxyz[i][3],
	       Gxyz[i][17],Gxyz[i][18],Gxyz[i][19] ); 
      }   
    }

  } /* if (MD_Opt_OK!=1 && iter!=MD_IterNumber) */

  Past_Utot[1]=Utot;

  Max_Force = force_Max;
  SD_scaling_user = SD_scaling*Max_Force*0.2;

  /* free arrays */

  free(A);
  free(B);
  free(ipiv);
  free(work);

  /*********************** end of "myid==Host_ID" **************************/

 Last_Bcast: 

  MPI_Bcast(&MD_Opt_OK,1,MPI_INT, Host_ID, mpi_comm_level1);

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    /* ID = G2ID[Gc_AN]; */
    MPI_Bcast(&Gxyz[Gc_AN][1],  3, MPI_DOUBLE, Host_ID, mpi_comm_level1);
    /*    MPI_Bcast(&Gxyz[Gc_AN][17], 3, MPI_DOUBLE, Host_ID, mpi_comm_level1); */
  }


  if (myid==Host_ID){ 

    if (0<level_stdout){
      printf("<%s>  |Maximum force| (Hartree/Bohr) =%15.12f\n",
	     func_name,Max_Force);fflush(stdout);
      printf("<%s>  Criterion       (Hartree/Bohr) =%15.12f\n",
	     func_name,MD_Opt_criterion);fflush(stdout);

      printf("\n");
      for (i=1; i<=atomnum; i++){
	printf("     atom=%4d, XYZ(ang) Fxyz(a.u.)=%9.4f %9.4f %9.4f  %9.4f %9.4f %9.4f\n",
	       i,BohrR*Gxyz[i][1],BohrR*Gxyz[i][2],BohrR*Gxyz[i][3],
	       Gxyz[i][17],Gxyz[i][18],Gxyz[i][19] ); fflush(stdout);
      }   
    }

    strcpy(fileSD,".SD");
    fnjoint(filepath,filename,fileSD);
    if ((fp_SD = fopen(fileSD,"a")) != NULL){

#ifdef xt3
      setvbuf(fp_SD,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      if (iter0==1){

        fprintf(fp_SD,"\n");
        fprintf(fp_SD,"***********************************************************\n");
        fprintf(fp_SD,"***********************************************************\n");
        fprintf(fp_SD,"              History of geometry optimization             \n");
        fprintf(fp_SD,"***********************************************************\n");
        fprintf(fp_SD,"***********************************************************\n\n");


        fprintf(fp_SD,"  MD_iter   SD_scaling     |Maximum force|   Maximum step        Utot\n");
        fprintf(fp_SD,"                           (Hartree/Bohr)        (Ang)         (Hartree)\n\n");
      }

      fprintf(fp_SD,"  %3d  %15.8f  %15.8f  %15.8f  %15.8f\n",
              iter0,SD_scaling,Max_Force,Max_Step*BohrR,Utot);
      fclose(fp_SD);
    }
    else{
      printf("Could not open a file in MD_pac.!\n");
    }

    if (MD_Opt_OK==1 || iter0==MD_IterNumber){

      strcpy(fileCoord,".crd");
      fnjoint(filepath,filename,fileCoord);
      if ((fp_crd = fopen(fileCoord,"w")) != NULL){

#ifdef xt3
        setvbuf(fp_crd,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

        fprintf(fp_crd,"\n");
        fprintf(fp_crd,"***********************************************************\n");
        fprintf(fp_crd,"***********************************************************\n");
        fprintf(fp_crd,"       xyz-coordinates (Ang) and forces (Hartree/Bohr)  \n");
        fprintf(fp_crd,"***********************************************************\n");
        fprintf(fp_crd,"***********************************************************\n\n");

        fprintf(fp_crd,"<coordinates.forces\n");
        fprintf(fp_crd,"  %i\n",atomnum);
        for (k=1; k<=atomnum; k++){
          i = WhatSpecies[k];
          j = Spe_WhatAtom[i];
          fprintf(fp_crd," %4d  %4s   %9.5f %9.5f %9.5f  %15.12f %15.12f %15.12f\n",
                  k,Atom_Symbol[j],
	          Gxyz[k][1]*BohrR,Gxyz[k][2]*BohrR,Gxyz[k][3]*BohrR,
	    	  -Gxyz[k][17],-Gxyz[k][18],-Gxyz[k][19]);
        }
        fprintf(fp_crd,"coordinates.forces>\n");
        fclose(fp_crd);
      }
      else
        printf("error(1) in MD_pac.c\n");
    }

  } /* if (myid==Host_ID) */

}




/*hmweng start RF method*/
void RF(int iter, int iter0)
{
  /* 1au=2.4189*10^-2 fs, 1fs=41.341105 au */

  int k1,k2,m1,m2;
  char *func_name="RF";
  char *JOBB="L";
  double dt,diff,Max_Step;
  double *A,*B,sumB,max_A, RR,dRi[4],dRj[4];
  double *work, *work2,**ahes;
/*hmweng c0, c1 are not used and lamda is newly defined */  
/*  double sum1,tmp1,tmp2,c0,c1; */
  double sum1,tmp1,tmp2,lamda;
  double mixing,force_Max;
  static double sMD_TimeStep,dx_max=0.05; 
  double diff_dx;
  int *ipiv;
  INTEGER i,j,k,iatom,N,LDA,LWORK,info;
  int diis_iter;
  char fileCoord[YOUSO10];
  char fileSD[YOUSO10];
  FILE *fp_crd,*fp_SD;
  char buf[fp_bsize];          /* setvbuf */
  char fileE[YOUSO10];

  /* variables for MPI */
  int Gc_AN;
  int numprocs,myid,ID;  

  /* MPI myid */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  /* share Gxyz */
  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    ID = G2ID[Gc_AN];
    MPI_Bcast(&Gxyz[Gc_AN][1],  3, MPI_DOUBLE, ID, mpi_comm_level1);
    MPI_Bcast(&Gxyz[Gc_AN][17], 3, MPI_DOUBLE, ID, mpi_comm_level1);
  }

  if (iter<M_GDIIS_HISTORY)
    diis_iter = iter;
  else   
    diis_iter = M_GDIIS_HISTORY;

  /* shift one */

  for (i=(diis_iter-2); 0<=i; i--) {
    for (iatom=1; iatom<=atomnum; iatom++) {
      for (k=1; k<=3; k++) {
	GxyzHistoryIn[i+1][iatom][k] = GxyzHistoryIn[i][iatom][k];
	GxyzHistoryR[i+1][iatom][k]  = GxyzHistoryR[i][iatom][k];
      }
    }
  }

  /* add GxyzHistoryIn and GxyzHisotryR */

  for (iatom=1; iatom<=atomnum; iatom++)   {

    for (k=1;k<=3;k++) {
      if (atom_Fixed_XYZ[iatom][k]==0){
	GxyzHistoryIn[0][iatom][k] =  Gxyz[iatom][k];
	GxyzHistoryR[0][iatom][k]  =  Gxyz[iatom][16+k];
      }
      else{
	GxyzHistoryIn[0][iatom][k] = Gxyz[iatom][k];
	GxyzHistoryR[0][iatom][k]  = 0.0;
      }
    }
  }

  if (myid!=Host_ID) goto Last_Bcast; 

  /*********************** for myid==Host_ID **************************/

  /* allocation of arrays */

  A = (double*)malloc(sizeof(double)*(diis_iter+2)*(diis_iter+2));
  for (i=0; i<(diis_iter+2)*(diis_iter+2); i++) A[i] = 0.0;
  B = (double*)malloc(sizeof(double)*(diis_iter+2));

  /* Max of force */

  force_Max=0.0;
  for (iatom=1; iatom<=atomnum; iatom++)   {
    for (k=1;k<=3;k++) {
      if (atom_Fixed_XYZ[iatom][k]==0){
	if (force_Max< fabs(Gxyz[iatom][16+k]) ) force_Max = fabs(Gxyz[iatom][16+k]);
      }
    }
  }

  sMD_TimeStep = 0.05/(0.01*41.341105);

  if (2<=level_stdout){
    printf("<%s>  |Maximum force| (Hartree/Bohr) = %15.12f tuned_dt= %f\n",func_name,force_Max, sMD_TimeStep);
    printf("<%s>  Criterion      (Hartree/Bohr) = %15.12f\n", func_name, MD_Opt_criterion);
  }

  for (i=0; i<diis_iter; i++) {
    for (j=0; j<diis_iter; j++) {

      RR = 0.0;

      for (iatom=1; iatom<=atomnum; iatom++)   {
	for (k=1; k<=3; k++) {
	  dRi[k] = GxyzHistoryR[i][iatom][k];
	  dRj[k] = GxyzHistoryR[j][iatom][k];
	}

	RR += dRi[1]*dRj[1] + dRi[2]*dRj[2] + dRi[3]*dRj[3];
      }

      A[ i*(diis_iter+1)+j ]= RR;
    }
  }

  /* find max of A */

  max_A = 0.0;

  for (i=0;i<diis_iter;i++) {
    for (j=0;j<diis_iter;j++) {
      RR = fabs(A[i*(diis_iter+1)+j]) ;
      if (max_A< RR ) max_A = RR;
    }
  }

  max_A = 1.0/max_A;

  for (i=0; i<diis_iter; i++) {
    for (j=0; j<diis_iter; j++) {
      A[ i*(diis_iter+1)+j ] *= max_A;
    }
  }

  for (i=0; i<diis_iter; i++) {
    A[ i*(diis_iter+1)+diis_iter ] = 1.0;
    A[ diis_iter*(diis_iter+1)+i ] = 1.0;
  }

  A[diis_iter*(diis_iter+1)+diis_iter] = 0.0;

  for (i=0; i<diis_iter; i++) B[i] = 0.0;
  B[diis_iter] = 1.0;

  if (2<=level_stdout){
    printf("<%s>  DIIS matrix\n",func_name);
    for (i=0; i<(diis_iter+1); i++) {
      printf("<%s> ",func_name);
      for (j=0; j<(diis_iter+1); j++) {
        printf("%6.3f ",A[i*(diis_iter+1)+j]);
      }
      printf("\n");
    }
  }

  /* lapack routine */

  N = diis_iter+1;
  LDA = diis_iter+1;
  LWORK = diis_iter+1;
  work = (double*)malloc(sizeof(double)*LWORK);
  ipiv = (int*)malloc(sizeof(int)*(diis_iter+1));

  i = 1; 

  if (2<=level_stdout){
    printf("M_GDIIS_HISTORY=%2d diis_iter=%2d\n",M_GDIIS_HISTORY,diis_iter);
  }

  F77_NAME(dsysv,DSYSV)( JOBB, &N, &i, A, &LDA,  ipiv, B, &LDA, work, &LWORK, &info);

  if (info!=0) {
    printf("<%s> dsysv_, info=%d\n",func_name,info);
    printf("<%s> \n",func_name);
    printf("<%s> ERROR, aborting\n",func_name);
    printf("<%s> \n",func_name);

    MD_Opt_OK =1; 
    /* no change */

    goto Last_Bcast ;
  }

  if (2<=level_stdout){
    printf("<%s> diis alpha=",func_name);
    sumB = 0;
    for (i=0; i<diis_iter; i++) {
      printf("%f ",B[i]);
      sumB += B[i];
    }
    printf("%lf\n",B[diis_iter]);
  }

  if (force_Max<MD_Opt_criterion ){
    MD_Opt_OK = 1;
    Max_Force = force_Max;
  }

  /****************************************************
   write informatins to *.ene
  ****************************************************/

  if (myid==Host_ID){  

    sprintf(fileE,"%s%s.ene",filepath,filename);
    iterout_md(iter,MD_TimeStep*(iter-1),fileE);
  }

  if (MD_Opt_OK!=1 && iter!=MD_IterNumber){

    /* initialize */

    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1; j<=3; j++) {
	Gxyz[iatom][j] = 0.0;
      }
    }

    /* add tilde{R} */

    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1; j<=3; j++) {
	if (atom_Fixed_XYZ[iatom][j]==0){
	  for (i=0; i<diis_iter; i++) {
	    Gxyz[iatom][j] += GxyzHistoryR[i][iatom][j]*B[i];
	  }
	}
      }
    }

    sumB = 0.0;
    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1; j<=3; j++) {
	sumB += Gxyz[iatom][j]*Gxyz[iatom][j] ;
      }
    }

    sumB = sqrt(sumB)/(double)atomnum;

    /************************************************************
   update an approximate Hessian matrix 
   
   H_k = H_{k-1} + y*y^t/(y^t*y) - H_{k-1}*s * s^t * H_{k-1}^t /(s^t * H_{k-1} * s)
   y = GxyzHistoryR[0][iatom][k]  - GxyzHistoryR[1][iatom][k]
   s = GxyzHistoryIn[0][iatom][k] - GxyzHistoryIn[1][iatom][k]
    ************************************************************/

    if (iter==1){
      for (i=1; i<=3*atomnum; i++){     
	for (j=1; j<=3*atomnum; j++){     
	  Hessian[i][j] = 0.0;
	}
	Hessian[i][i] = 1.0;
      }
    }

    else {
   
      /* H*s  */

      for (i=1; i<=3*atomnum; i++){     

	sum1 = 0.0;  
	for (k=1; k<=atomnum; k++){     

	  sum1 += ( Hessian[i][3*k-2]*(GxyzHistoryIn[0][k][1] - GxyzHistoryIn[1][k][1])
		   +Hessian[i][3*k-1]*(GxyzHistoryIn[0][k][2] - GxyzHistoryIn[1][k][2])
		   +Hessian[i][3*k]*(GxyzHistoryIn[0][k][3]   - GxyzHistoryIn[1][k][3]));
	}
 
	/* store H*s */

	Hessian[0][i] = sum1;
      }        

      /* tmp1 = y^t*s, tmp2 = s^t*H*s */ 
    
      tmp1 = 0.0;
      tmp2 = 0.0;
    
      for (k=1; k<=atomnum; k++){     

	tmp1 += ((GxyzHistoryIn[0][k][1] - GxyzHistoryIn[1][k][1])*(GxyzHistoryR[0][k][1] - GxyzHistoryR[1][k][1])
		 +(GxyzHistoryIn[0][k][2] - GxyzHistoryIn[1][k][2])*(GxyzHistoryR[0][k][2] - GxyzHistoryR[1][k][2])
		 +(GxyzHistoryIn[0][k][3] - GxyzHistoryIn[1][k][3])*(GxyzHistoryR[0][k][3] - GxyzHistoryR[1][k][3]));

	tmp2 += ((GxyzHistoryIn[0][k][1] - GxyzHistoryIn[1][k][1])*Hessian[0][3*k-2]
		 +(GxyzHistoryIn[0][k][2] - GxyzHistoryIn[1][k][2])*Hessian[0][3*k-1] 
		 +(GxyzHistoryIn[0][k][3] - GxyzHistoryIn[1][k][3])*Hessian[0][3*k  ]); 
      }

      /* c0=(tmp1+tmp2)/(tmp1*tmp1), c1=-1.0/tmp1 */

      /*
	c0 = (tmp1 + tmp2)/(tmp1*tmp1);
	c1 =-1.0/tmp1;
      */
      /* update the approximate Hessian by the BFGS method */

      m1 = 0;
      for (i=1; i<=atomnum; i++){   
	for (k1=1; k1<=3; k1++){     

	  m1++;

	  m2 = 0;
	  for (j=1; j<=atomnum; j++){     
	    for (k2=1; k2<=3; k2++){     

	      m2++;

	      Hessian[m1][m2] += 

		1.0/tmp1*(GxyzHistoryR[0][i][k1]-GxyzHistoryR[1][i][k1])*(GxyzHistoryR[0][j][k2]-GxyzHistoryR[1][j][k2])
		- 1.0/tmp2*Hessian[0][m1]*Hessian[0][m2];

	    }
	  }
	}
      }
    }

    /************************************************************
            Construct augmented Hessian matrix  
           | H    g | ( s )         ( s ) 
           |        |       = lamda             
           | g    0 | ( 1 )         ( 1 )
    ************************************************************/
  
    Hessian[3*atomnum+1][3*atomnum+1]=0.0;
  
    m2 = 0;
    for (i=1; i<=atomnum; i++){   
      for (k=1; k<=3; k++){   
	m2++;
	Hessian[3*atomnum+1][m2]=Gxyz[i][k];
	Hessian[m2][3*atomnum+1]=Gxyz[i][k];  
      }	
    }
  
    /************************************************************
     find the lowest eigenvalue and corresponding eigenvector
           of the augmented Hessian matrix
    ************************************************************/  

    /*  lamda = find_Emin(3*atomnum+1, Hessian); */

    work2 = (double*)malloc(sizeof(double)*(3*atomnum+3));

    ahes = (double**)malloc(sizeof(double*)*(3*atomnum+3));
    for (i=0; i<3*atomnum+3; i++){
      ahes[i] = (double*)malloc(sizeof(double)*(3*atomnum+3));
    }

    for(i=0; i<3*atomnum+2; i++){
      for(j=0;j<3*atomnum+2; j++){
	ahes[i][j] = Hessian[i][j];
      }
    }
    Eigen_lapack(ahes, work2, 3*atomnum+1, 1); 
 
    for(i=1; i<=1; i++){
      printf("Eigenvalue=%15.8f\n",work2[i]);
    }

    printf("EigenVector is\n");

    for(i=1;i<=3*atomnum+1; i++){
      printf("%15.8f\n",ahes[i][1]);
      Hessian[0][i]=ahes[i][1]/ahes[3*atomnum+1][1];
    }
 
    if (2<=level_stdout){
      printf("<%s> |tilde{R}|=%E\n",func_name, sumB);
    }

    /* initialize */

    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1;j<=3;j++) {
	Gxyz[iatom][j] = 0.0;
      }
    }

    /* calculate the DIIS coordinates tilde{x} */

    m1 = 0;
    for (iatom=1;iatom<=atomnum;iatom++) {
      for (j=1;j<=3;j++) {

	for (i=0; i<diis_iter; i++) {
	  Gxyz[iatom][j] += GxyzHistoryIn[i][iatom][j]*B[i]; 
	}

	/* a quasi Newton method */ 

	m1++;
	Gxyz[iatom][j] += Hessian[0][m1]; 

      }
    }

    /************************************************************
        In case of a too large updating, do a modest updating
    ************************************************************/  

    Max_Step = 0.0;
    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1; j<=3; j++) {

        diff = fabs(Gxyz[iatom][j] - GxyzHistoryIn[0][iatom][j]);
        if (Max_Step<diff) Max_Step = diff;
      }
    }

    if (Criterion_Max_Step<Max_Step){

      for (iatom=1; iatom<=atomnum; iatom++) {
	for (j=1; j<=3; j++) {

	  diff = Gxyz[iatom][j] - GxyzHistoryIn[0][iatom][j];
	  Gxyz[iatom][j] = GxyzHistoryIn[0][iatom][j] + diff/Max_Step*Criterion_Max_Step;
	}
      }
    }

    /* find Max_Step */
    Max_Step = 0.0;
    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1; j<=3; j++) {

        diff = fabs(Gxyz[iatom][j] - GxyzHistoryIn[0][iatom][j]);
        if (Max_Step<diff) Max_Step = diff;
      }
    }

    /************************************************************
                   show cooridinates and gradients
    ************************************************************/  

    if (2<=level_stdout){

      printf("<%s> diff_x= %f , dE= %f\n",func_name, Max_Step, fabs(Utot-Past_Utot[1]) );

      /* print atomic positions */
      printf("<%s> atomnum= %d\n",func_name,atomnum);
      for (i=1; i<=atomnum; i++){
	j = Spe_WhatAtom[WhatSpecies[i]];
	printf("  %3d %s XYZ(ang) Fxyz(a.u.)= %9.4f %9.4f %9.4f  %9.4f %9.4f %9.4f\n",
	       i,Atom_Symbol[j],
	       BohrR*Gxyz[i][1],BohrR*Gxyz[i][2],BohrR*Gxyz[i][3],
	       Gxyz[i][17],Gxyz[i][18],Gxyz[i][19] ); 
      }   
    }

    Past_Utot[1]=Utot;

    Max_Force = force_Max;
    SD_scaling_user = SD_scaling*Max_Force*0.2;

    /* free arrays */

    free(A);
    free(B);
    free(ipiv);
    free(work);
    free(work2);

    for (i=0; i<(3*atomnum+3); i++){
      free(ahes[i]);
    }
    free(ahes);

  } /* if (MD_Opt_OK!=1 && iter!=MD_IterNumber) */

  /*********************** end of "myid==Host_ID" **************************/

 Last_Bcast: 

  MPI_Bcast(&MD_Opt_OK,1,MPI_INT, Host_ID, mpi_comm_level1);

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    /* ID = G2ID[Gc_AN]; */
    MPI_Bcast(&Gxyz[Gc_AN][1],  3, MPI_DOUBLE, Host_ID, mpi_comm_level1);
    /*    MPI_Bcast(&Gxyz[Gc_AN][17], 3, MPI_DOUBLE, Host_ID, mpi_comm_level1); */
  }

  if (myid==Host_ID){ 

    if (0<level_stdout){ 

      printf("<%s>  |Maximum force| (Hartree/Bohr) =%15.12f\n",
	     func_name,Max_Force);fflush(stdout);
      printf("<%s>  Criterion       (Hartree/Bohr) =%15.12f\n",
	     func_name,MD_Opt_criterion);fflush(stdout);

      printf("\n");
      for (i=1; i<=atomnum; i++){
	printf("     atom=%4d, XYZ(ang) Fxyz(a.u.)=%9.4f %9.4f %9.4f  %9.4f %9.4f %9.4f\n",
	       i,BohrR*Gxyz[i][1],BohrR*Gxyz[i][2],BohrR*Gxyz[i][3],
	       Gxyz[i][17],Gxyz[i][18],Gxyz[i][19] ); fflush(stdout);
      }   
    }

    strcpy(fileSD,".SD");
    fnjoint(filepath,filename,fileSD);

    if ((fp_SD = fopen(fileSD,"a")) != NULL){

#ifdef xt3
      setvbuf(fp_SD,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      if (iter0==1){

        fprintf(fp_SD,"\n");
        fprintf(fp_SD,"***********************************************************\n");
        fprintf(fp_SD,"***********************************************************\n");
        fprintf(fp_SD,"              History of geometry optimization             \n");
        fprintf(fp_SD,"***********************************************************\n");
        fprintf(fp_SD,"***********************************************************\n\n");


        fprintf(fp_SD,"  MD_iter   SD_scaling     |Maximum force|   Maximum step        Utot\n");
        fprintf(fp_SD,"                           (Hartree/Bohr)        (Ang)         (Hartree)\n\n");
      }

      fprintf(fp_SD,"  %3d  %15.8f  %15.8f  %15.8f  %15.8f\n",
              iter0,SD_scaling,Max_Force,Max_Step*BohrR,Utot);
      fclose(fp_SD);
    }
    else{
      printf("Could not open a file in MD_pac.!\n");
    }

    if (MD_Opt_OK==1 || iter0==MD_IterNumber){

      strcpy(fileCoord,".crd");
      fnjoint(filepath,filename,fileCoord);
      if ((fp_crd = fopen(fileCoord,"w")) != NULL){

#ifdef xt3
        setvbuf(fp_crd,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

        fprintf(fp_crd,"\n");
        fprintf(fp_crd,"***********************************************************\n");
        fprintf(fp_crd,"***********************************************************\n");
        fprintf(fp_crd,"       xyz-coordinates (Ang) and forces (Hartree/Bohr)  \n");
        fprintf(fp_crd,"***********************************************************\n");
        fprintf(fp_crd,"***********************************************************\n\n");

        fprintf(fp_crd,"<coordinates.forces\n");
        fprintf(fp_crd,"  %i\n",atomnum);
        for (k=1; k<=atomnum; k++){
          i = WhatSpecies[k];
          j = Spe_WhatAtom[i];
          fprintf(fp_crd," %4d  %4s   %9.5f %9.5f %9.5f  %15.12f %15.12f %15.12f\n",
                  k,Atom_Symbol[j],
	          Gxyz[k][1]*BohrR,Gxyz[k][2]*BohrR,Gxyz[k][3]*BohrR,
	    	  -Gxyz[k][17],-Gxyz[k][18],-Gxyz[k][19]);
        }
        fprintf(fp_crd,"coordinates.forces>\n");
        fclose(fp_crd);
      }
      else
        printf("error(1) in MD_pac.c\n");
    }

  } /* if (myid==Host_ID) */

}



void GDIIS_BFGS(int iter, int iter0)
{
  /* 1au=2.4189*10^-2 fs, 1fs=41.341105 au */

  int k1,k2,m1,m2;
  char *func_name="BFGS";
  char *JOBB="L";
  double dt;
  double *A,*B,sumB,max_A, RR,dRi[4],dRj[4];
  double *work;
  double sum1,tmp1,tmp2,c0,c1;
  double mixing,force_Max;
  static double sMD_TimeStep,dx_max=0.05; 
  double diff_dx,diff,Max_Step;
  int *ipiv;
  INTEGER i,j,k,iatom,N,LDA,LWORK,info;
  int diis_iter;
  char fileCoord[YOUSO10];
  char fileSD[YOUSO10];
  FILE *fp_crd,*fp_SD;
  char buf[fp_bsize];          /* setvbuf */
  char fileE[YOUSO10];

  /* variables for MPI */
  int Gc_AN;
  int numprocs,myid,ID;  

  /* MPI myid */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  /* share Gxyz */
  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    ID = G2ID[Gc_AN];
    MPI_Bcast(&Gxyz[Gc_AN][1],  3, MPI_DOUBLE, ID, mpi_comm_level1);
    MPI_Bcast(&Gxyz[Gc_AN][17], 3, MPI_DOUBLE, ID, mpi_comm_level1);
  }

  if (iter<M_GDIIS_HISTORY)
    diis_iter = iter;
  else   
    diis_iter = M_GDIIS_HISTORY;

  /* shift one */

  for (i=(diis_iter-2); 0<=i; i--) {
    for (iatom=1; iatom<=atomnum; iatom++) {
      for (k=1; k<=3; k++) {
	GxyzHistoryIn[i+1][iatom][k] = GxyzHistoryIn[i][iatom][k];
	GxyzHistoryR[i+1][iatom][k]  = GxyzHistoryR[i][iatom][k];
      }
    }
  }

  /* Max of force */

  force_Max=0.0;
  for (iatom=1; iatom<=atomnum; iatom++)   {
    for (k=1;k<=3;k++) {
      if (atom_Fixed_XYZ[iatom][k]==0){
	if (force_Max< fabs(Gxyz[iatom][16+k]) ) force_Max = fabs(Gxyz[iatom][16+k]);
      }
    }
  }

  sMD_TimeStep = 0.05/(0.01*41.341105);

  if (2<=level_stdout && myid==Host_ID){
    printf("<%s>  |Maximum force| (Hartree/Bohr) = %15.12f tuned_dt= %f\n",func_name,force_Max, sMD_TimeStep);
    printf("<%s>  Criterion      (Hartree/Bohr) = %15.12f\n", func_name, MD_Opt_criterion);
  }

  /* add GxyzHistoryIn and GxyzHisotryR */

  for (iatom=1; iatom<=atomnum; iatom++)   {

    for (k=1;k<=3;k++) {
      if (atom_Fixed_XYZ[iatom][k]==0){
	GxyzHistoryIn[0][iatom][k] = Gxyz[iatom][k];
	GxyzHistoryR[0][iatom][k]  = Gxyz[iatom][16+k];
      }
      else{
	GxyzHistoryIn[0][iatom][k] = Gxyz[iatom][k];
	GxyzHistoryR[0][iatom][k]  = 0.0;
      }
    }
  }

  if (myid!=Host_ID)  goto Last_Bcast; 

  /*********************** for myid==Host_ID **************************/

  /* allocation of arrays */

  A = (double*)malloc(sizeof(double)*(diis_iter+1)*(diis_iter+1));
  for (i=0; i<(diis_iter+1)*(diis_iter+1); i++) A[i] = 0.0;
  B = (double*)malloc(sizeof(double)*(diis_iter+1));

  for (i=0; i<diis_iter; i++) {
    for (j=0; j<diis_iter; j++) {

      RR = 0.0;

      for (iatom=1; iatom<=atomnum; iatom++)   {
	for (k=1; k<=3; k++) {
	  dRi[k] = GxyzHistoryR[i][iatom][k];
	  dRj[k] = GxyzHistoryR[j][iatom][k];
	}

	RR += dRi[1]*dRj[1] + dRi[2]*dRj[2] + dRi[3]*dRj[3];
      }

      A[ i*(diis_iter+1)+j ]= RR;
    }
  }

  /* find max of A */

  max_A = 0.0;

  for (i=0;i<diis_iter;i++) {
    for (j=0;j<diis_iter;j++) {
      RR = fabs(A[i*(diis_iter+1)+j]) ;
      if (max_A< RR ) max_A = RR;
    }
  }

  max_A = 1.0/max_A;

  for (i=0; i<diis_iter; i++) {
    for (j=0; j<diis_iter; j++) {
      A[ i*(diis_iter+1)+j ] *= max_A;
    }
  }

  for (i=0; i<diis_iter; i++) {
    A[ i*(diis_iter+1)+diis_iter ] = 1.0;
    A[ diis_iter*(diis_iter+1)+i ] = 1.0;
  }

  A[diis_iter*(diis_iter+1)+diis_iter] = 0.0;

  for (i=0; i<diis_iter; i++) B[i] = 0.0;
  B[diis_iter] = 1.0;

  if (2<=level_stdout){
    printf("<%s>  DIIS matrix\n",func_name);
    for (i=0; i<(diis_iter+1); i++) {
      printf("<%s> ",func_name);
      for (j=0; j<(diis_iter+1); j++) {
        printf("%10.5f ",A[i*(diis_iter+1)+j]);
      }
      printf("\n");
    }
  }

  /* lapack routine */

  N=diis_iter+1;
  LDA=diis_iter+1;
  LWORK=diis_iter+1;
  work=(double*)malloc(sizeof(double)*LWORK);
  ipiv = (int*)malloc(sizeof(int)*(diis_iter+1));

  i = 1; 

  if (2<=level_stdout){
    printf("M_GDIIS_HISTORY=%2d diis_iter=%2d\n",M_GDIIS_HISTORY,diis_iter);
  }

  F77_NAME(dsysv,DSYSV)( JOBB, &N, &i, A, &LDA,  ipiv, B, &LDA, work, &LWORK, &info);

  if (info!=0) {
    printf("<%s> dsysv_, info=%d\n",func_name,info);
    printf("<%s> \n",func_name);
    printf("<%s> ERROR, aborting\n",func_name);
    printf("<%s> \n",func_name);

    MD_Opt_OK =1; 
    /* no change */

    goto Last_Bcast ;
  }

  if (2<=level_stdout){
    printf("<%s> diis alpha=",func_name);
    sumB = 0;
    for (i=0; i<diis_iter; i++) {
      printf("%10.5f ",B[i]);
      sumB += B[i];
    }
    printf("%lf\n",B[diis_iter]);
  }

  if (force_Max<MD_Opt_criterion )  MD_Opt_OK = 1;

  /****************************************************
   write informatins to *.ene
  ****************************************************/

  if (myid==Host_ID){  

    sprintf(fileE,"%s%s.ene",filepath,filename);
    iterout_md(iter,MD_TimeStep*(iter-1),fileE);
  }

  if (MD_Opt_OK!=1 && iter!=MD_IterNumber){

    /* initialize */

    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1; j<=3; j++) {
	Gxyz[iatom][j] = 0.0;
      }
    }

    /* add tilde{R} */

    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1; j<=3; j++) {
	if (atom_Fixed_XYZ[iatom][j]==0){
	  for (i=0; i<diis_iter; i++) {
	    Gxyz[iatom][j] += GxyzHistoryR[i][iatom][j]*B[i];
	  }
	}
      }
    }

    sumB = 0.0;
    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1; j<=3; j++) {
	sumB += Gxyz[iatom][j]*Gxyz[iatom][j] ;
      }
    }

    sumB = sqrt(sumB)/(double)atomnum;

    /************************************************************
     update an approximate inverse of Hessian matrix 

     y = GxyzHistoryR[0][iatom][k]  - GxyzHistoryR[1][iatom][k]
     s = GxyzHistoryIn[0][iatom][k] - GxyzHistoryIn[1][iatom][k]
    ************************************************************/

    if (iter==1){
      for (i=1; i<=3*atomnum; i++){     
	for (j=1; j<=3*atomnum; j++){     
	  InvHessian[i][j] = 0.0;
	}
	InvHessian[i][i] = 1.0;
      }
    }

    else {
   
      /* invH*y  */

      for (i=1; i<=3*atomnum; i++){     

	sum1 = 0.0;  
	for (k=1; k<=atomnum; k++){     

	  sum1 += ( InvHessian[i][3*k-2]*(GxyzHistoryR[0][k][1] - GxyzHistoryR[1][k][1])
		   +InvHessian[i][3*k-1]*(GxyzHistoryR[0][k][2] - GxyzHistoryR[1][k][2])
		   +InvHessian[i][3*k  ]*(GxyzHistoryR[0][k][3] - GxyzHistoryR[1][k][3]));
	}
 
	/* store invH*y */

	InvHessian[0][i] = sum1;
      }        

      /* tmp1 = s^t*y, tmp2 = y^t*H*y */ 

      tmp1 = 0.0;
      tmp2 = 0.0;

      for (k=1; k<=atomnum; k++){     

	tmp1 += ( (GxyzHistoryIn[0][k][1] - GxyzHistoryIn[1][k][1])*(GxyzHistoryR[0][k][1] - GxyzHistoryR[1][k][1])
		 +(GxyzHistoryIn[0][k][2] - GxyzHistoryIn[1][k][2])*(GxyzHistoryR[0][k][2] - GxyzHistoryR[1][k][2])
		 +(GxyzHistoryIn[0][k][3] - GxyzHistoryIn[1][k][3])*(GxyzHistoryR[0][k][3] - GxyzHistoryR[1][k][3]));

	tmp2 += ( (GxyzHistoryR[0][k][1] - GxyzHistoryR[1][k][1])*InvHessian[0][3*k-2]
		 +(GxyzHistoryR[0][k][2] - GxyzHistoryR[1][k][2])*InvHessian[0][3*k-1] 
		 +(GxyzHistoryR[0][k][3] - GxyzHistoryR[1][k][3])*InvHessian[0][3*k  ]); 
      }

      /* c0=(tmp1+tmp2)/(tmp1*tmp1), c1=-1.0/tmp1 */

      c0 = (tmp1 + tmp2)/(tmp1*tmp1);
      c1 =-1.0/tmp1;

      /* update the approximate Hessian by the BFGS method */

      m1 = 0;
      for (i=1; i<=atomnum; i++){   
	for (k1=1; k1<=3; k1++){     

	  m1++;

	  m2 = 0;
	  for (j=1; j<=atomnum; j++){     
	    for (k2=1; k2<=3; k2++){     

	      m2++;

	      InvHessian[m1][m2] += 

		c0*(GxyzHistoryIn[0][i][k1]-GxyzHistoryIn[1][i][k1])*(GxyzHistoryIn[0][j][k2]-GxyzHistoryIn[1][j][k2])
		+ c1*InvHessian[0][m1]*(GxyzHistoryIn[0][j][k2]-GxyzHistoryIn[1][j][k2])
		+ c1*(GxyzHistoryIn[0][i][k1]-GxyzHistoryIn[1][i][k1])*InvHessian[0][m2];

	    }
	  }
	}
      }
    }

    /************************************************************
            perform a quasi Newton-Raphson method  
    ************************************************************/

    /* H^-1*g */

    m1 = 0;
    for (i=1; i<=atomnum; i++){   
      for (k1=1; k1<=3; k1++){     

	m1++;
  
	m2 = 0;
	sum1 = 0.0;
	for (j=1; j<=atomnum; j++){     
	  for (k2=1; k2<=3; k2++){     

	    m2++;

	    sum1 += InvHessian[m1][m2]*Gxyz[j][k2];
	  }
	}

        InvHessian[0][m1] = sum1;

      }
    }

    if (2<=level_stdout){
      printf("<%s> |tilde{R}|=%E\n",func_name, sumB);
    }

    /* initialize */

    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1; j<=3; j++) {
	Gxyz[iatom][j] = 0.0;
      }
    }

    /* calculate the DIIS coordinates tilde{x} */

    m1 = 0;
    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1; j<=3; j++) {

	for (i=0; i<diis_iter; i++) {
	  Gxyz[iatom][j] += GxyzHistoryIn[i][iatom][j]*B[i]; 
	}

	/* a quasi Newton method */ 

	m1++;
	Gxyz[iatom][j] -= InvHessian[0][m1]; 

      }
    }

    /************************************************************
        In case of a too large updating, do a modest updating
    ************************************************************/  

    Max_Step = 0.0;
    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1; j<=3; j++) {

        diff = fabs(Gxyz[iatom][j] - GxyzHistoryIn[0][iatom][j]);
        if (Max_Step<diff) Max_Step = diff;
      }
    }

    if (Criterion_Max_Step<Max_Step){

      for (iatom=1; iatom<=atomnum; iatom++) {
	for (j=1; j<=3; j++) {

	  diff = Gxyz[iatom][j] - GxyzHistoryIn[0][iatom][j];
	  Gxyz[iatom][j] = GxyzHistoryIn[0][iatom][j] + diff/Max_Step*Criterion_Max_Step;
	}
      }
    }

    /* find Max_Step */
    Max_Step = 0.0;
    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1; j<=3; j++) {

        diff = fabs(Gxyz[iatom][j] - GxyzHistoryIn[0][iatom][j]);
        if (Max_Step<diff) Max_Step = diff;
      }
    }

    /************************************************************
                   show cooridinates and gradients
    ************************************************************/  

    if (2<=level_stdout){

      printf("<%s> diff_x= %f , dE= %f\n",func_name, Max_Step, fabs(Utot-Past_Utot[1]) );

      /* print atomic positions */
      printf("<%s> atomnum= %d\n",func_name,atomnum);
      for (i=1; i<=atomnum; i++){
	j = Spe_WhatAtom[WhatSpecies[i]];
	printf("  %3d %s XYZ(ang) Fxyz(a.u.)= %9.4f %9.4f %9.4f  %9.4f %9.4f %9.4f\n",
	       i,Atom_Symbol[j],
	       BohrR*Gxyz[i][1],BohrR*Gxyz[i][2],BohrR*Gxyz[i][3],
	       Gxyz[i][17],Gxyz[i][18],Gxyz[i][19] ); 
      }   
    }

  } /* if (MD_Opt_OK!=1 && iter!=MD_IterNumber) */

  Past_Utot[1]=Utot;

  Max_Force = force_Max;
  SD_scaling_user = SD_scaling*Max_Force*0.2;

  /* free arrays */

  free(A);
  free(B);
  free(ipiv);
  free(work);

  /*********************** end of "myid==Host_ID" **************************/

 Last_Bcast: 

  MPI_Bcast(&MD_Opt_OK,1,MPI_INT, Host_ID, mpi_comm_level1);

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    /* ID = G2ID[Gc_AN]; */
    MPI_Bcast(&Gxyz[Gc_AN][1],  3, MPI_DOUBLE, Host_ID, mpi_comm_level1);
    /*    MPI_Bcast(&Gxyz[Gc_AN][17], 3, MPI_DOUBLE, Host_ID, mpi_comm_level1); */
  }

  if (myid==Host_ID){ 

    if (0<level_stdout){ 

      printf("<%s>  |Maximum force| (Hartree/Bohr) =%15.12f\n",
	     func_name,Max_Force);fflush(stdout);
      printf("<%s>  Criterion       (Hartree/Bohr) =%15.12f\n",
	     func_name,MD_Opt_criterion);fflush(stdout);

      printf("\n");
      for (i=1; i<=atomnum; i++){
	printf("     atom=%4d, XYZ(ang) Fxyz(a.u.)=%9.4f %9.4f %9.4f  %9.4f %9.4f %9.4f\n",
	       i,BohrR*Gxyz[i][1],BohrR*Gxyz[i][2],BohrR*Gxyz[i][3],
	       Gxyz[i][17],Gxyz[i][18],Gxyz[i][19] ); fflush(stdout);
      }   
    }

    strcpy(fileSD,".SD");
    fnjoint(filepath,filename,fileSD);
    if ((fp_SD = fopen(fileSD,"a")) != NULL){

#ifdef xt3
      setvbuf(fp_SD,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      if (iter0==1){

        fprintf(fp_SD,"\n");
        fprintf(fp_SD,"***********************************************************\n");
        fprintf(fp_SD,"***********************************************************\n");
        fprintf(fp_SD,"              History of geometry optimization             \n");
        fprintf(fp_SD,"***********************************************************\n");
        fprintf(fp_SD,"***********************************************************\n\n");


        fprintf(fp_SD,"  MD_iter   SD_scaling     |Maximum force|   Maximum step        Utot\n");
        fprintf(fp_SD,"                           (Hartree/Bohr)        (Ang)         (Hartree)\n\n");
      }

      fprintf(fp_SD,"  %3d  %15.8f  %15.8f  %15.8f  %15.8f\n",
              iter0,SD_scaling,Max_Force,Max_Step*BohrR,Utot);
      fclose(fp_SD);
    }
    else{
      printf("Could not open a file in MD_pac.!\n");
    }

    if (MD_Opt_OK==1 || iter0==MD_IterNumber){

      strcpy(fileCoord,".crd");
      fnjoint(filepath,filename,fileCoord);
      if ((fp_crd = fopen(fileCoord,"w")) != NULL){

#ifdef xt3
        setvbuf(fp_crd,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

        fprintf(fp_crd,"\n");
        fprintf(fp_crd,"***********************************************************\n");
        fprintf(fp_crd,"***********************************************************\n");
        fprintf(fp_crd,"       xyz-coordinates (Ang) and forces (Hartree/Bohr)  \n");
        fprintf(fp_crd,"***********************************************************\n");
        fprintf(fp_crd,"***********************************************************\n\n");

        fprintf(fp_crd,"<coordinates.forces\n");
        fprintf(fp_crd,"  %i\n",atomnum);
        for (k=1; k<=atomnum; k++){
          i = WhatSpecies[k];
          j = Spe_WhatAtom[i];
          fprintf(fp_crd," %4d  %4s   %9.5f %9.5f %9.5f  %15.12f %15.12f %15.12f\n",
                  k,Atom_Symbol[j],
	          Gxyz[k][1]*BohrR,Gxyz[k][2]*BohrR,Gxyz[k][3]*BohrR,
	    	  -Gxyz[k][17],-Gxyz[k][18],-Gxyz[k][19]);
        }
        fprintf(fp_crd,"coordinates.forces>\n");
        fclose(fp_crd);
      }
      else
        printf("error(1) in MD_pac.c\n");
    }

  } /* if (myid==Host_ID) */

}









void GDIIS_EF(int iter, int iter0)
{
  /* 1au=2.4189*10^-2 fs, 1fs=41.341105 au */

  int k1,k2,m1,m2,isp,wan,wsp;
  static double Utot0,scaling_factor;
  char *func_name="EF";
  char *JOBB="L";
  double dt,MinKo;
  double *A,*B,sumB,max_A, RR,dRi[4],dRj[4];
  double *ko,**U;
  double *work;
/*hmweng c0, c1 are not used and lamda is newly defined */  
/*  double sum1,tmp1,tmp2,c0,c1; */
  double sum,sum1,tmp1,tmp2,lamda,c0,c1;
  double mixing,force_Max;
  static double sMD_TimeStep,dx_max=0.05; 
  double diff_dx,diff,Max_Step;
  int *ipiv;
  INTEGER i,j,k,iatom,N,LDA,LWORK,info;
  int diis_iter;
  char fileCoord[YOUSO10];
  char fileSD[YOUSO10];
  FILE *fp_crd,*fp_SD;
  char buf[fp_bsize];          /* setvbuf */
  char fileE[YOUSO10];

  /* variables for MPI */
  int Gc_AN;
  int numprocs,myid,ID;  

  /* MPI myid */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  /* share Gxyz */
  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    ID = G2ID[Gc_AN];
    MPI_Bcast(&Gxyz[Gc_AN][1],  3, MPI_DOUBLE, ID, mpi_comm_level1);
    MPI_Bcast(&Gxyz[Gc_AN][17], 3, MPI_DOUBLE, ID, mpi_comm_level1);
  }

  if (iter<M_GDIIS_HISTORY)
    diis_iter = iter;
  else   
    diis_iter = M_GDIIS_HISTORY;

  /* shift one */

  for (i=(diis_iter-2); 0<=i; i--) {
    for (iatom=1; iatom<=atomnum; iatom++) {
      for (k=1; k<=3; k++) {
	GxyzHistoryIn[i+1][iatom][k] = GxyzHistoryIn[i][iatom][k];
	GxyzHistoryR[i+1][iatom][k]  = GxyzHistoryR[i][iatom][k];
      }
    }
  }

  /* add GxyzHistoryIn and GxyzHisotryR */

  for (iatom=1; iatom<=atomnum; iatom++)   {

    for (k=1; k<=3; k++) {
      if (atom_Fixed_XYZ[iatom][k]==0){
	GxyzHistoryIn[0][iatom][k] = Gxyz[iatom][k];
	GxyzHistoryR[0][iatom][k]  = Gxyz[iatom][16+k];
      }
      else{
	GxyzHistoryIn[0][iatom][k] = Gxyz[iatom][k];
	GxyzHistoryR[0][iatom][k]  = 0.0;
      }
    }
  }
  
  if (myid!=Host_ID)   goto Last_Bcast; 

  /*********************** for myid==Host_ID **************************/

  /* allocation of arrays */

  A = (double*)malloc(sizeof(double)*(diis_iter+1)*(diis_iter+1));
  for (i=0; i<(diis_iter+1)*(diis_iter+1); i++) A[i] = 0.0;
  B = (double*)malloc(sizeof(double)*(diis_iter+1));
  ko = (double*)malloc(sizeof(double)*(3*atomnum+2));
  U = (double**)malloc(sizeof(double*)*(3*atomnum+2));
  for (i=0; i<(3*atomnum+2); i++){
    U[i] = (double*)malloc(sizeof(double)*(3*atomnum+2));
  }

  /* Max of force */

  force_Max=0.0;
  for (iatom=1; iatom<=atomnum; iatom++)   {

    wsp = WhatSpecies[iatom];
    wan = Spe_WhatAtom[wsp];

    for (k=1;k<=3;k++) {

      /*
      if (atom_Fixed_XYZ[iatom][k]==0){
      */

      if (atom_Fixed_XYZ[iatom][k]==0 && wan!=0){

	if (force_Max< fabs(Gxyz[iatom][16+k]) ) force_Max = fabs(Gxyz[iatom][16+k]);
      }

    }
  }

  sMD_TimeStep = 0.05/(0.01*41.341105);

  if (2<=level_stdout){
    printf("<%s>  |Maximum force| (Hartree/Bohr) = %15.12f tuned_dt= %f\n",func_name,force_Max, sMD_TimeStep);
    printf("<%s>  Criterion      (Hartree/Bohr) = %15.12f\n", func_name, MD_Opt_criterion);
  }

  for (i=0; i<diis_iter; i++) {
    for (j=0; j<diis_iter; j++) {

      RR = 0.0;

      for (iatom=1; iatom<=atomnum; iatom++)   {
	for (k=1; k<=3; k++) {
	  dRi[k] = GxyzHistoryR[i][iatom][k];
	  dRj[k] = GxyzHistoryR[j][iatom][k];
	}

	RR += dRi[1]*dRj[1] + dRi[2]*dRj[2] + dRi[3]*dRj[3];
      }

      A[ i*(diis_iter+1)+j ]= RR;
    }
  }

  /* find max of A */

  max_A = 0.0;

  for (i=0;i<diis_iter;i++) {
    for (j=0;j<diis_iter;j++) {
      RR = fabs(A[i*(diis_iter+1)+j]) ;
      if (max_A< RR ) max_A = RR;
    }
  }

  max_A = 1.0/max_A;

  for (i=0; i<diis_iter; i++) {
    for (j=0; j<diis_iter; j++) {
      A[ i*(diis_iter+1)+j ] *= max_A;
    }
  }

  for (i=0; i<diis_iter; i++) {
    A[ i*(diis_iter+1)+diis_iter ] = 1.0;
    A[ diis_iter*(diis_iter+1)+i ] = 1.0;
  }

  A[diis_iter*(diis_iter+1)+diis_iter] = 0.0;

  for (i=0; i<diis_iter; i++) B[i] = 0.0;
  B[diis_iter] = 1.0;
  
  if (2<=level_stdout){
    printf("<%s>  DIIS matrix\n",func_name);
    for (i=0; i<(diis_iter+1); i++) {
      printf("<%s> ",func_name);
      for (j=0; j<(diis_iter+1); j++) {
        printf("%6.3f ",A[i*(diis_iter+1)+j]);
      }
      printf("\n");
    }
  }

  /* lapack routine */

  N=diis_iter+1;
  LDA=diis_iter+1;
  LWORK=diis_iter+1;
  work=(double*)malloc(sizeof(double)*LWORK);
  ipiv = (int*)malloc(sizeof(int)*(diis_iter+1));

  i = 1; 

  if (2<=level_stdout){
    printf("M_GDIIS_HISTORY=%2d diis_iter=%2d\n",M_GDIIS_HISTORY,diis_iter);
  }

  F77_NAME(dsysv,DSYSV)( JOBB, &N, &i, A, &LDA,  ipiv, B, &LDA, work, &LWORK, &info);

  if (info!=0) {
    printf("<%s> dsysv_, info=%d\n",func_name,info);
    printf("<%s> \n",func_name);
    printf("<%s> ERROR, aborting\n",func_name);
    printf("<%s> \n",func_name);

    MD_Opt_OK =1; 
    /* no change */

    goto Last_Bcast ;
  }

  if (2<=level_stdout){
    printf("<%s> diis alpha=",func_name);
    sumB = 0;
    for (i=0; i<diis_iter; i++) {
      printf("%f ",B[i]);
      sumB += B[i];
    }
    printf("%lf\n",B[diis_iter]);
  }

  if (force_Max<MD_Opt_criterion )  MD_Opt_OK = 1;
    
  /****************************************************
   write informatins to *.ene
  ****************************************************/
    
  if (myid==Host_ID){  
    
    sprintf(fileE,"%s%s.ene",filepath,filename);
    iterout_md(iter,MD_TimeStep*(iter-1),fileE);
  } 

  if (MD_Opt_OK!=1 && iter!=MD_IterNumber){

    /* initialize */

    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1; j<=3; j++) {
	Gxyz[iatom][j] = 0.0;
      }
    }

    /* add tilde{R} */

    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1; j<=3; j++) {
	if (atom_Fixed_XYZ[iatom][j]==0){
	  for (i=0; i<diis_iter; i++) {
	    Gxyz[iatom][j] += GxyzHistoryR[i][iatom][j]*B[i];
	  }
	}
      }
    }

    /* store tilde{R} into Hessian[][0] */

    m1 = 0;
    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1; j<=3; j++) {
	m1++;
	Hessian[m1][0] = Gxyz[iatom][j];
      }
    }

    sumB = 0.0;
    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1; j<=3; j++) {
	sumB += Gxyz[iatom][j]*Gxyz[iatom][j] ;
      }
    }

    sumB = sqrt(sumB)/(double)atomnum;

    if (2<=level_stdout){
      printf("<%s> |tilde{R}|=%E\n",func_name, sumB);
    }

    /************************************************************
   update an approximate Hessian matrix 
   
   H_k = H_{k-1} + y*y^t/(y^t*y) - H_{k-1}*s * s^t * H_{k-1}^t /(s^t * H_{k-1} * s)
   y = GxyzHistoryR[0][iatom][k]  - GxyzHistoryR[1][iatom][k]
   s = GxyzHistoryIn[0][iatom][k] - GxyzHistoryIn[1][iatom][k]
    ************************************************************/

    if (iter==1){

      scaling_factor = 1.0;

      for (i=1; i<=3*atomnum; i++){     
	for (j=1; j<=3*atomnum; j++){     
	  Hessian[i][j] = 0.0;
	}
	Hessian[i][i] = 1.0;
      }
    }

    else {
   
      /* H*s  */

      for (i=1; i<=3*atomnum; i++){     

	sum1 = 0.0;  
	for (k=1; k<=atomnum; k++){     

	  sum1 += (Hessian[i][3*k-2]*(GxyzHistoryIn[0][k][1] - GxyzHistoryIn[1][k][1])
		   +Hessian[i][3*k-1]*(GxyzHistoryIn[0][k][2] - GxyzHistoryIn[1][k][2])
		   +Hessian[i][3*k  ]*(GxyzHistoryIn[0][k][3] - GxyzHistoryIn[1][k][3]));
	}
 
	/* store H*s */

	Hessian[0][i] = sum1;
      }        

      /* tmp1 = y^t*s, tmp2 = s^t*H*s */ 
    
      tmp1 = 0.0;
      tmp2 = 0.0;
    
      for (k=1; k<=atomnum; k++){     

	tmp1 += ((GxyzHistoryIn[0][k][1] - GxyzHistoryIn[1][k][1])*(GxyzHistoryR[0][k][1] - GxyzHistoryR[1][k][1])
		 +(GxyzHistoryIn[0][k][2] - GxyzHistoryIn[1][k][2])*(GxyzHistoryR[0][k][2] - GxyzHistoryR[1][k][2])
		 +(GxyzHistoryIn[0][k][3] - GxyzHistoryIn[1][k][3])*(GxyzHistoryR[0][k][3] - GxyzHistoryR[1][k][3]));

	tmp2 += ((GxyzHistoryIn[0][k][1] - GxyzHistoryIn[1][k][1])*Hessian[0][3*k-2]
		 +(GxyzHistoryIn[0][k][2] - GxyzHistoryIn[1][k][2])*Hessian[0][3*k-1] 
		 +(GxyzHistoryIn[0][k][3] - GxyzHistoryIn[1][k][3])*Hessian[0][3*k  ]); 
      }

      /* c0=1.0/tmp1, c1=1.0/tmp2 */

      /*
	if (myid==Host_ID){
	printf("tmp1=%15.12f tmp2=%15.12f\n",1.0/tmp1,1.0/tmp2);
	}
      */

      c0 = 1.0/tmp1;
      c1 = 1.0/tmp2;

      /* update the approximate Hessian by the BFGS method if 0.0<c0 */
    
      if (0.0<c0){
	m1 = 0;
	for (i=1; i<=atomnum; i++){   
	  for (k1=1; k1<=3; k1++){     

	    m1++;

	    m2 = 0;
	    for (j=1; j<=atomnum; j++){     
	      for (k2=1; k2<=3; k2++){     

		m2++;

		Hessian[m1][m2] += 

		  c0*(GxyzHistoryR[0][i][k1]-GxyzHistoryR[1][i][k1])*(GxyzHistoryR[0][j][k2]-GxyzHistoryR[1][j][k2])
		  -c1*Hessian[0][m1]*Hessian[0][m2];

	      }
	    }
	  }
	}

      }
    } 

    /************************************************************
             diagonalize the approximate Hessian
    ************************************************************/

    for (i=1; i<=3*atomnum; i++){
      for (j=1; j<=3*atomnum; j++){
	U[i][j] = Hessian[i][j];
      }
    }

    Eigen_lapack(U,ko,3*atomnum,3*atomnum);


    /*
      if (myid==Host_ID){
      for (i=1; i<=3*atomnum; i++){
      printf("i=%3d ko=%15.12f\n",i,ko[i]);
      }

      printf("EV\n"); 
      for (i=1; i<=3*atomnum; i++){
      for (j=1; j<=3*atomnum; j++){
      printf("%8.4f ",U[i][j]); 
      }
      printf("\n");
      }
      }
    */


    isp = 0;

    /*
      for (i=1; i<=3*atomnum; i++){
      if (ko[i]<1.0e-1) isp = i;
      }
    */

    if (atomnum<=4) MinKo = 0.10;
    else            MinKo = 0.02;

    for (i=1; i<=3*atomnum; i++){
      if (ko[i]<MinKo) ko[i] = MinKo;
    }

    if (isp!=0 && myid==Host_ID && 0<level_stdout){
      printf("Hessian is ill-conditioned.\n");
    } 

    /************************************************************
             U*lambda^{-1} U^{dag} * g
    ************************************************************/

    for (i=(isp+1); i<=3*atomnum; i++){
      sum = 0.0;
      for (j=1; j<=3*atomnum; j++){
	sum += U[j][i]*Hessian[j][0];
      }
      U[i][0] = sum;
    }  
  
    for (i=1; i<=3*atomnum; i++){
      sum = 0.0;
      for (j=(isp+1); j<=3*atomnum; j++){
	sum += U[i][j]*U[j][0]/ko[j];
      }
      U[0][i] = sum;
    }

    /************************************************************
            calculate the DIIS coordinates tilde{x}
              and perform a quasi Newton method
    ************************************************************/

    /* update scaling_factor */

    if (Utot0<Utot && iter!=1) scaling_factor = 0.95*scaling_factor;

    /* initialize */

    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1;j<=3;j++) {
	Gxyz[iatom][j] = 0.0;
      }
    }

    m1 = 0;
    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1; j<=3; j++) {

	/* DIIS coordinate */

	for (i=0; i<diis_iter; i++) {
	  Gxyz[iatom][j] += GxyzHistoryIn[i][iatom][j]*B[i]; 
	}

	/* with a quasi Newton method */ 

	m1++;

	Gxyz[iatom][j] -= scaling_factor*U[0][m1];

      }
    }

    /************************************************************
        In case of a too large updating, do a modest updating
    ************************************************************/  

    Max_Step = 0.0;
    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1; j<=3; j++) {

        diff = fabs(Gxyz[iatom][j] - GxyzHistoryIn[0][iatom][j]);
        if (Max_Step<diff) Max_Step = diff;
      }
    }

    if (Criterion_Max_Step<Max_Step){

      for (iatom=1; iatom<=atomnum; iatom++) {
	for (j=1; j<=3; j++) {

	  diff = Gxyz[iatom][j] - GxyzHistoryIn[0][iatom][j];
	  Gxyz[iatom][j] = GxyzHistoryIn[0][iatom][j] + diff/Max_Step*Criterion_Max_Step;
	}
      }
    }

    /* find Max_Step */
    Max_Step = 0.0;
    for (iatom=1; iatom<=atomnum; iatom++) {
      for (j=1; j<=3; j++) {

        diff = fabs(Gxyz[iatom][j] - GxyzHistoryIn[0][iatom][j]);
        if (Max_Step<diff) Max_Step = diff;
      }
    }

    /************************************************************
                    show coordinates and gradients
    ************************************************************/  

    if (2<=level_stdout){

      printf("<%s> diff_x= %f , dE= %f\n",func_name, Max_Step, fabs(Utot-Past_Utot[1]) );

      /* print atomic positions */
      printf("<%s> atomnum= %d\n",func_name,atomnum);
      for (i=1; i<=atomnum; i++){
	j = Spe_WhatAtom[WhatSpecies[i]];
	printf("  %3d %s XYZ(ang) Fxyz(a.u.)= %9.4f %9.4f %9.4f  %9.4f %9.4f %9.4f\n",
	       i,Atom_Symbol[j],
	       BohrR*Gxyz[i][1],BohrR*Gxyz[i][2],BohrR*Gxyz[i][3],
	       Gxyz[i][17],Gxyz[i][18],Gxyz[i][19] ); 
      }   
    }

  } /* if (MD_Opt_OK!=1 && iter!=MD_IterNumber) */

  Past_Utot[1]=Utot;

  Max_Force = force_Max;
  SD_scaling_user = SD_scaling*Max_Force*0.2;

  /* save Utot */

  Utot0 = Utot;

  /* free arrays */

  free(A);
  free(B);
  free(ko);
  for (i=0; i<(3*atomnum+2); i++){
    free(U[i]);
  }
  free(U);

  free(ipiv);
  free(work);

  /*********************** end of "myid==Host_ID" **************************/

 Last_Bcast: 

  MPI_Bcast(&MD_Opt_OK,1,MPI_INT, Host_ID, mpi_comm_level1);

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    /* ID = G2ID[Gc_AN]; */
    MPI_Bcast(&Gxyz[Gc_AN][1],  3, MPI_DOUBLE, Host_ID, mpi_comm_level1);
    /*    MPI_Bcast(&Gxyz[Gc_AN][17], 3, MPI_DOUBLE, Host_ID, mpi_comm_level1); */
  }


  if (myid==Host_ID){ 

    if (0<level_stdout){ 

      printf("<%s>  |Maximum force| (Hartree/Bohr) =%15.12f\n",
	     func_name,Max_Force);fflush(stdout);
      printf("<%s>  Criterion       (Hartree/Bohr) =%15.12f\n",
	     func_name,MD_Opt_criterion);fflush(stdout);

      printf("\n");
      for (i=1; i<=atomnum; i++){
	printf("     atom=%4d, XYZ(ang) Fxyz(a.u.)=%9.4f %9.4f %9.4f  %9.4f %9.4f %9.4f\n",
	       i,BohrR*Gxyz[i][1],BohrR*Gxyz[i][2],BohrR*Gxyz[i][3],
	       Gxyz[i][17],Gxyz[i][18],Gxyz[i][19] ); fflush(stdout);
      }   
    }

    strcpy(fileSD,".SD");
    fnjoint(filepath,filename,fileSD);
    if ((fp_SD = fopen(fileSD,"a")) != NULL){

#ifdef xt3
      setvbuf(fp_SD,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      if (iter0==1){

        fprintf(fp_SD,"\n");
        fprintf(fp_SD,"***********************************************************\n");
        fprintf(fp_SD,"***********************************************************\n");
        fprintf(fp_SD,"              History of geometry optimization             \n");
        fprintf(fp_SD,"***********************************************************\n");
        fprintf(fp_SD,"***********************************************************\n\n");


        fprintf(fp_SD,"  MD_iter   SD_scaling     |Maximum force|   Maximum step        Utot\n");
        fprintf(fp_SD,"                           (Hartree/Bohr)        (Ang)         (Hartree)\n\n");
      }

      fprintf(fp_SD,"  %3d  %15.8f  %15.8f  %15.8f  %15.8f\n",
              iter0,SD_scaling,Max_Force,Max_Step*BohrR,Utot);
      fclose(fp_SD);
    }
    else{
      printf("Could not open a file in MD_pac.!\n");
    }

    if (MD_Opt_OK==1 || iter0==MD_IterNumber){

      strcpy(fileCoord,".crd");
      fnjoint(filepath,filename,fileCoord);
      if ((fp_crd = fopen(fileCoord,"w")) != NULL){

#ifdef xt3
        setvbuf(fp_crd,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

        fprintf(fp_crd,"\n");
        fprintf(fp_crd,"***********************************************************\n");
        fprintf(fp_crd,"***********************************************************\n");
        fprintf(fp_crd,"       xyz-coordinates (Ang) and forces (Hartree/Bohr)  \n");
        fprintf(fp_crd,"***********************************************************\n");
        fprintf(fp_crd,"***********************************************************\n\n");

        fprintf(fp_crd,"<coordinates.forces\n");
        fprintf(fp_crd,"  %i\n",atomnum);
        for (k=1; k<=atomnum; k++){
          i = WhatSpecies[k];
          j = Spe_WhatAtom[i];
          fprintf(fp_crd," %4d  %4s   %9.5f %9.5f %9.5f  %15.12f %15.12f %15.12f\n",
                  k,Atom_Symbol[j],
	          Gxyz[k][1]*BohrR,Gxyz[k][2]*BohrR,Gxyz[k][3]*BohrR,
	    	  -Gxyz[k][17],-Gxyz[k][18],-Gxyz[k][19]);
        }
        fprintf(fp_crd,"coordinates.forces>\n");
        fclose(fp_crd);
      }
      else
        printf("error(1) in MD_pac.c\n");
    }

  } /* if (myid==Host_ID) */

}



















void NVT_VS(int iter)
{
  /* added by mari */
  /********************************************************
   This routine is added by Mari Ohfuti (May 20004).                

   a constant temperature molecular dynamics by a velocity
   scaling method with velocity-Verlet integrator
  ********************************************************/
  /******************************************************* 
   1 a.u.=2.4189*10^-2 fs, 1fs=41.341105 a.u. 
   Atom weight trasformation: proton = 1836.1526 a.u 
  ********************************************************/

  /****************************************************
    Gxyz[][1] = x-coordinate at current step
    Gxyz[][2] = y-coordinate at current step
    Gxyz[][3] = z-coordinate at current step

    Gxyz[][14] = dEtot/dx at previous step
    Gxyz[][15] = dEtot/dy at previous step
    Gxyz[][16] = dEtot/dz at previous step

    Gxyz[][17] = dEtot/dx at current step
    Gxyz[][18] = dEtot/dy at current step
    Gxyz[][19] = dEtot/dz at current step

    Gxyz[][20] = atomic mass

    Gxyz[][21] = x-coordinate at previous step
    Gxyz[][22] = y-coordinate at previous step
    Gxyz[][23] = z-coordinate at previous step

    Gxyz[][24] = x-component of velocity at current step
    Gxyz[][25] = y-component of velocity at current step
    Gxyz[][26] = z-component of velocity at current step

    Gxyz[][27] = x-component of velocity at t+dt/2
    Gxyz[][28] = y-component of velocity at t+dt/2
    Gxyz[][29] = z-component of velocity at t+dt/2

    Gxyz[][30] = hx
    Gxyz[][31] = hy
    Gxyz[][32] = hz

  ****************************************************/

  double dt,dt2,sum,My_Ukc,x,t,xyz0[4],xyz0_l[4];
  double Wscale;
  int Mc_AN,Gc_AN,i,j,k,l;
  int numprocs,myid,ID;
  char fileE[YOUSO10];

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  MD_Opt_OK = 0;
  dt = 41.3411*MD_TimeStep;
  dt2 = dt*dt;
  Wscale = 1836.1526;

  /****************************************************
                Velocity Verlet algorithm
  ****************************************************/

  if (iter==1){
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];

      for (j=1; j<=3; j++){

        if (atom_Fixed_XYZ[Gc_AN][j]==0){

	  Gxyz[Gc_AN][j] = Gxyz[Gc_AN][j]+dt*Gxyz[Gc_AN][23+j]
                          -Gxyz[Gc_AN][16+j]/(Gxyz[Gc_AN][20]*Wscale)*dt2*0.50;
          Gxyz[Gc_AN][13+j] = Gxyz[Gc_AN][16+j];

	}
      }
    }
  }
  else{
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];

      for (j=1; j<=3; j++){

        if (atom_Fixed_XYZ[Gc_AN][j]==0){
 	  Gxyz[Gc_AN][23+j] = Gxyz[Gc_AN][23+j]-(Gxyz[Gc_AN][16+j]
                             +Gxyz[Gc_AN][13+j])/(Gxyz[Gc_AN][20]*Wscale)*dt*0.50;

	}
      }
    }
  }

  /****************************************************
                     Kinetic Energy 
  ****************************************************/

  Ukc=0.0;
  My_Ukc = 0.0;

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];

    sum = 0.0;
    for (j=1; j<=3; j++){
      sum = sum + Gxyz[Gc_AN][j+23]*Gxyz[Gc_AN][j+23];
    }
    My_Ukc = My_Ukc + 0.5*Gxyz[Gc_AN][20]*Wscale*sum;
  }

  /****************************************************
   MPI: Ukc 
  ****************************************************/

  MPI_Allreduce(&My_Ukc, &Ukc, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

  /* calculation of temperature (K) */
  Temp = Ukc/(1.5*kB*(double)atomnum)*eV2Hartree;

  /* calculation of a given temperature (K) */
  for (i=1; i<=TempNum; i++) {
    if( (iter>NumScale[i-1]) && (iter<=NumScale[i]) ) {
      GivenTemp = TempPara[i-1][2] + (TempPara[i][2] - TempPara[i-1][2])*
          ((double)iter-(double)TempPara[i-1][1])/((double)TempPara[i][1]-(double)TempPara[i-1][1]);
    }
  }

  /****************************************************
   write informatins to *.ene
  ****************************************************/

  if (myid==Host_ID){  

    sprintf(fileE,"%s%s.ene",filepath,filename);
    iterout_md(iter,MD_TimeStep*(iter-1),fileE);
  }

  /****************************************************
   velocity scaling 
  ****************************************************/

  if(iter!=1) {

    x = 1.0;
    for (i=1; i<=TempNum; i++) {

      if( (iter>NumScale[i-1]) && (iter<=NumScale[i]) ) {

        /**************************************************
         find a scaling parameter, x, when MD step matches
         at the step where the temperature scaling is made.
         Otherwise, x = 1.0.
        **************************************************/

        if((iter-NumScale[i-1])%IntScale[i]==0) {

          GivenTemp = TempPara[i-1][2] + (TempPara[i][2] - TempPara[i-1][2])*
               ((double)iter-(double)TempPara[i-1][1])/((double)TempPara[i][1]-(double)TempPara[i-1][1]);
 
          x = GivenTemp + (Temp-GivenTemp)*RatScale[i];
          x = sqrt(1.5*kB*x/(Ukc*eV2Hartree)*(double)atomnum);
        }
      }
    }

    /* do scaling */

    if (MD_Opt_OK!=1 && iter!=MD_IterNumber){

      for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
	Gc_AN = M2G[Mc_AN];
	for (j=1; j<=3; j++){

	  if (atom_Fixed_XYZ[Gc_AN][j]==0){

	    Gxyz[Gc_AN][23+j] = Gxyz[Gc_AN][23+j]*x;
	    Gxyz[Gc_AN][j] = Gxyz[Gc_AN][j]+dt*Gxyz[Gc_AN][23+j]
	      -Gxyz[Gc_AN][16+j]/(Gxyz[Gc_AN][20]*Wscale)*dt2*0.50;
	    Gxyz[Gc_AN][13+j] = Gxyz[Gc_AN][16+j];

	  }
	}
      }

    }

  }

  /****************************************************
   MPI: Gxyz
  ****************************************************/

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    ID = G2ID[Gc_AN];
    MPI_Bcast(&Gxyz[Gc_AN][1],  3, MPI_DOUBLE, ID, mpi_comm_level1);
    MPI_Bcast(&Gxyz[Gc_AN][14], 3, MPI_DOUBLE, ID, mpi_comm_level1);
    MPI_Bcast(&Gxyz[Gc_AN][24], 3, MPI_DOUBLE, ID, mpi_comm_level1);
  }

}




void NVT_VS2(int iter)
{
  /* added by T.Ohwaki */
  /********************************************************
   This routine is added by T.Ohwaki (May 20004).                

   a constant temperature molecular dynamics by a velocity
   scaling method for each atom with velocity-Verlet integrator
  ********************************************************/
  /******************************************************* 
   1 a.u.=2.4189*10^-2 fs, 1fs=41.341105 a.u. 
   Atom weight trasformation: proton = 1836.1526 a.u 
  ********************************************************/

  /****************************************************
    Gxyz[][1] = x-coordinate at current step
    Gxyz[][2] = y-coordinate at current step
    Gxyz[][3] = z-coordinate at current step

    Gxyz[][14] = dEtot/dx at previous step
    Gxyz[][15] = dEtot/dy at previous step
    Gxyz[][16] = dEtot/dz at previous step

    Gxyz[][17] = dEtot/dx at current step
    Gxyz[][18] = dEtot/dy at current step
    Gxyz[][19] = dEtot/dz at current step

    Gxyz[][20] = atomic mass

    Gxyz[][21] = x-coordinate at previous step
    Gxyz[][22] = y-coordinate at previous step
    Gxyz[][23] = z-coordinate at previous step

    Gxyz[][24] = x-component of velocity at current step
    Gxyz[][25] = y-component of velocity at current step
    Gxyz[][26] = z-component of velocity at current step
  ****************************************************/

  double dt,dt2,sum,My_Ukc,x,t,xyz0[4],xyz0_l[4];
  double Wscale;

  double *Atomic_Temp,*Atomic_Ukc,*Atomic_Scale;

  int Mc_AN,Gc_AN,i,j,k,l;
  int numprocs,myid,ID;

  char fileE[YOUSO10];

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  MD_Opt_OK = 0;
  dt = 41.3411*MD_TimeStep;
  dt2 = dt*dt;
  Wscale = 1836.1526;

  Atomic_Temp  = (double*)malloc(sizeof(double)*(atomnum+1));
  Atomic_Ukc   = (double*)malloc(sizeof(double)*(atomnum+1));
  Atomic_Scale = (double*)malloc(sizeof(double)*(atomnum+1));

  /****************************************************
                Velocity Verlet algorithm
  ****************************************************/

  if (iter==1){
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];

      for (j=1; j<=3; j++){

        if (atom_Fixed_XYZ[Gc_AN][j]==0){

	  Gxyz[Gc_AN][j] = Gxyz[Gc_AN][j]+dt*Gxyz[Gc_AN][23+j]
                          -Gxyz[Gc_AN][16+j]/(Gxyz[Gc_AN][20]*Wscale)*dt2*0.50;
          Gxyz[Gc_AN][13+j] = Gxyz[Gc_AN][16+j];

	}
      }
    }
  }
  else{
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];

      for (j=1; j<=3; j++){

        if (atom_Fixed_XYZ[Gc_AN][j]==0){
 	  Gxyz[Gc_AN][23+j] = Gxyz[Gc_AN][23+j]-(Gxyz[Gc_AN][16+j]
                             +Gxyz[Gc_AN][13+j])/(Gxyz[Gc_AN][20]*Wscale)*dt*0.50;

	}
      }
    }
  }

  /****************************************************
      Kinetic Energy & Temperature (for each atom)
  ****************************************************/

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];

    sum = 0.0;
    for (j=1; j<=3; j++){
      sum = sum + Gxyz[Gc_AN][j+23]*Gxyz[Gc_AN][j+23];
    }

    Atomic_Ukc[Gc_AN]  = 0.5*Gxyz[Gc_AN][20]*Wscale*sum;
    Atomic_Temp[Gc_AN] = Atomic_Ukc[Gc_AN]/(1.5*kB)*eV2Hartree;
  }


  /****************************************************
            Kinetic Energy (for whole system) 
  ****************************************************/

  Ukc=0.0;
  My_Ukc = 0.0;

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];

    sum = 0.0;
    for (j=1; j<=3; j++){
      sum = sum + Gxyz[Gc_AN][j+23]*Gxyz[Gc_AN][j+23];
    }
    My_Ukc = My_Ukc + 0.5*Gxyz[Gc_AN][20]*Wscale*sum;
  }

  /****************************************************
   MPI: Ukc 
  ****************************************************/

  MPI_Allreduce(&My_Ukc, &Ukc, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

  /* calculation of temperature (K) */
  Temp = Ukc/(1.5*kB*(double)atomnum)*eV2Hartree;

  /****************************************************
   write informatins to *.ene
  ****************************************************/

  if (myid==Host_ID){  

    sprintf(fileE,"%s%s.ene",filepath,filename);
    iterout_md(iter,MD_TimeStep*(iter-1),fileE);
  }

  /****************************************************
   velocity scaling for each atom
  ****************************************************/

  if(iter!=1) {

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];
      Atomic_Scale[Gc_AN] = 1.0;
    }

    for (i=1; i<=TempNum; i++) {

      if( (iter>NumScale[i-1]) && (iter<=NumScale[i]) ) {

        /*************************************************************
         find a scaling parameter, Atomic_Scale, when MD step matches
         at the step where the temperature scaling is made.
         Otherwise, Atomic_Scale = 1.0.
        *************************************************************/

        if((iter-NumScale[i-1])%IntScale[i]==0) {

          GivenTemp = TempPara[i-1][2] + (TempPara[i][2] - TempPara[i-1][2])*
               ((double)iter-(double)TempPara[i-1][1])/((double)TempPara[i][1]-(double)TempPara[i-1][1]);
 
          for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
            Gc_AN = M2G[Mc_AN];

            Atomic_Scale[Gc_AN] = sqrt((GivenTemp + (Atomic_Temp[Gc_AN]-GivenTemp)*RatScale[i])/Atomic_Temp[Gc_AN]);

          }

        }
      }
    }

    /* do scaling */

    if (MD_Opt_OK!=1 && iter!=MD_IterNumber){

      for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
	Gc_AN = M2G[Mc_AN];
	for (j=1; j<=3; j++){

	  if (atom_Fixed_XYZ[Gc_AN][j]==0){

	    Gxyz[Gc_AN][23+j] = Gxyz[Gc_AN][23+j]*Atomic_Scale[Gc_AN];
	    Gxyz[Gc_AN][j] = Gxyz[Gc_AN][j]+dt*Gxyz[Gc_AN][23+j]
	                    -Gxyz[Gc_AN][16+j]/(Gxyz[Gc_AN][20]*Wscale)*dt2*0.50;
	    Gxyz[Gc_AN][13+j] = Gxyz[Gc_AN][16+j];

	  }
	}
      }

    }

  }

  /****************************************************
   MPI: Gxyz
  ****************************************************/

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    ID = G2ID[Gc_AN];
    MPI_Bcast(&Gxyz[Gc_AN][1],  3, MPI_DOUBLE, ID, mpi_comm_level1);
    MPI_Bcast(&Gxyz[Gc_AN][14], 3, MPI_DOUBLE, ID, mpi_comm_level1);
    MPI_Bcast(&Gxyz[Gc_AN][24], 3, MPI_DOUBLE, ID, mpi_comm_level1);
  }

  free(Atomic_Temp);
  free(Atomic_Ukc);
  free(Atomic_Scale);

}



void NVT_NH(int iter)
{
  /***********************************************************
   a constant temperature molecular dynamics by a Nose-Hoover
   method with velocity-Verlet integrator
  ***********************************************************/
  /*********************************************************** 
   1 a.u.=2.4189*10^-2 fs, 1fs=41.341105 a.u. 
   Atom weight trasformation: proton = 1836.1526 a.u 
  ***********************************************************/
  /****************************************************
    Gxyz[][1] = x-coordinate at current step
    Gxyz[][2] = y-coordinate at current step
    Gxyz[][3] = z-coordinate at current step

    Gxyz[][14] = dEtot/dx at previous step
    Gxyz[][15] = dEtot/dy at previous step
    Gxyz[][16] = dEtot/dz at previous step

    Gxyz[][17] = dEtot/dx at current step
    Gxyz[][18] = dEtot/dy at current step
    Gxyz[][19] = dEtot/dz at current step

    Gxyz[][20] = atomic mass

    Gxyz[][21] = x-coordinate at previous step
    Gxyz[][22] = y-coordinate at previous step
    Gxyz[][23] = z-coordinate at previous step

    Gxyz[][24] = x-component of velocity at current step
    Gxyz[][25] = y-component of velocity at current step
    Gxyz[][26] = z-component of velocity at current step

    Gxyz[][27] = x-component of velocity at t+dt/2
    Gxyz[][28] = y-component of velocity at t+dt/2
    Gxyz[][29] = z-component of velocity at t+dt/2

    Gxyz[][30] = hx
    Gxyz[][31] = hy
    Gxyz[][32] = hz

  ****************************************************/

  int Mc_AN,Gc_AN,i,j,k,l,po,num,NH_switch;
  int numprocs,myid,ID;

  double dt,dt2,sum,My_sum,My_Ukc,x,t,xyz0[4],xyz0_l[4];
  double scaled_force,Wscale,back;
  double dzeta,dv,h_zeta;
  double My_sum1,sum1,My_sum2,sum2;
  char fileE[YOUSO10];

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  MD_Opt_OK = 0;
  dt = 41.3411*MD_TimeStep;
  dt2 = dt*dt;
  Wscale = 1836.1526;

  /* find a given temperature by a linear interpolation */

  NH_switch = 0;
  i = 1;
  do {

    if ( TempPara[i][1]<=iter && iter<TempPara[i+1][1] ){

      GivenTemp = TempPara[i][2] + (TempPara[i+1][2] - TempPara[i][2])*
            ((double)iter-(double)TempPara[i][1])/((double)TempPara[i+1][1]-(double)TempPara[i][1]);

      NH_switch = 1; 
    }

    i++;
  } while (NH_switch==0 && i<=(TempNum-1));  

  /****************************************************
                Velocity Verlet algorithm
  ****************************************************/

  if (iter==1){

    /****************************************************
                       Kinetic Energy 
    ****************************************************/

    My_Ukc = 0.0;

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];
      sum = 0.0;
      for (j=1; j<=3; j++){
        if (atom_Fixed_XYZ[Gc_AN][j]==0){
          sum += Gxyz[Gc_AN][j+23]*Gxyz[Gc_AN][j+23];
	}
      }
      My_Ukc = My_Ukc + 0.5*Gxyz[Gc_AN][20]*Wscale*sum;
    }

    /****************************************************
     MPI, Ukc
    ****************************************************/

    MPI_Allreduce(&My_Ukc, &Ukc, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

    /* calculation of temperature (K) */
    Temp = Ukc/(1.5*kB*(double)atomnum)*eV2Hartree;

    /****************************************************
     write informatins to *.ene
    ****************************************************/

    if (myid==Host_ID){  

      sprintf(fileE,"%s%s.ene",filepath,filename);
      iterout_md(iter,MD_TimeStep*(iter-1),fileE);
    }

    /****************************************************
      first step in velocity Verlet 
    ****************************************************/

    NH_czeta = 0.0;
    NH_R = 0.0;

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];
      for (j=1; j<=3; j++){

        if (atom_Fixed_XYZ[Gc_AN][j]==0){
        
          scaled_force = -Gxyz[Gc_AN][16+j]/(Gxyz[Gc_AN][20]*Wscale);  

          /* v( r+0.5*dt ) */
          Gxyz[Gc_AN][26+j] = Gxyz[Gc_AN][23+j] + (scaled_force - NH_czeta*Gxyz[Gc_AN][23+j])*0.5*dt;
          Gxyz[Gc_AN][23+j] = Gxyz[Gc_AN][26+j];

          /* r( r+dt ) */
          Gxyz[Gc_AN][20+j] = Gxyz[Gc_AN][j];
 	  Gxyz[Gc_AN][j] =  Gxyz[Gc_AN][j] + Gxyz[Gc_AN][26+j]*dt;
	}

      }
    }

    /* zeta( t+0.5*dt ) */

    NH_nzeta = NH_czeta + (Ukc - 1.5*kB*(double)atomnum*GivenTemp/eV2Hartree)*dt/(TempQ*Wscale);
    NH_czeta = NH_nzeta;

    /* R( r+dt ) */
    NH_R = NH_R + NH_nzeta*dt;
  }

  else if (NH_switch==1) {

    /*****************************************************
     second step:

     refinement of v and zeta by a Newton-Raphson method
    *****************************************************/

    po = 0;
    num = 0;

    do {

      /* Ukc */

      Ukc = 0.0;
      My_Ukc = 0.0;

      for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
        Gc_AN = M2G[Mc_AN];

        sum = 0.0;
        for (j=1; j<=3; j++){
          if (atom_Fixed_XYZ[Gc_AN][j]==0){
            sum = sum + Gxyz[Gc_AN][j+23]*Gxyz[Gc_AN][j+23];
	  }
        }
        My_Ukc += 0.5*Gxyz[Gc_AN][20]*Wscale*sum;
      }

      /* MPI: Ukc */

      MPI_Allreduce(&My_Ukc, &Ukc, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

      /* calculation of h */ 

      h_zeta = NH_nzeta
              + (Ukc - 1.5*kB*(double)atomnum*GivenTemp/eV2Hartree)*dt/(TempQ*Wscale) - NH_czeta;

      for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
        Gc_AN = M2G[Mc_AN];
        for (j=1; j<=3; j++){

          if (atom_Fixed_XYZ[Gc_AN][j]==0){

            scaled_force = -Gxyz[Gc_AN][16+j]/(Gxyz[Gc_AN][20]*Wscale);  
            Gxyz[Gc_AN][j+29] = Gxyz[Gc_AN][26+j] + (scaled_force - NH_czeta*Gxyz[Gc_AN][23+j])*0.5*dt
                                -Gxyz[Gc_AN][23+j];
	  }
	}  
      }

      /* sum1 */
     
      sum1=0.0;
      My_sum1 = 0.0;

      for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
        Gc_AN = M2G[Mc_AN];

        sum = 0.0;
        for (j=1; j<=3; j++){
          if (atom_Fixed_XYZ[Gc_AN][j]==0){
            sum += Gxyz[Gc_AN][j+29]*Gxyz[Gc_AN][j+23];
	  }
        }
        My_sum1 += Gxyz[Gc_AN][20]*Wscale*sum*dt/(TempQ*Wscale);
      }

      /* MPI: sum1 */

      MPI_Allreduce(&My_sum1, &sum1, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

      /* sum2 */
     
      sum2=0.0;
      My_sum2 = 0.0;

      for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
        Gc_AN = M2G[Mc_AN];

        sum = 0.0;
        for (j=1; j<=3; j++){
          if (atom_Fixed_XYZ[Gc_AN][j]==0){
            sum += Gxyz[Gc_AN][j+23]*Gxyz[Gc_AN][j+23];
	  }
        }
        My_sum2 -= 0.5*Gxyz[Gc_AN][20]*Wscale*sum*dt*dt/(TempQ*Wscale);
      }

      /* MPI: sum2 */

      MPI_Allreduce(&My_sum2, &sum2, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

      /* new NH_czeta and new v */

      dzeta = (-h_zeta*(NH_czeta*0.5*dt+1.0)-sum1)/(-(NH_czeta*0.5*dt+1.0)+sum2);

      My_sum = 0.0;
      for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
        Gc_AN = M2G[Mc_AN];
        for (j=1; j<=3; j++){

          if (atom_Fixed_XYZ[Gc_AN][j]==0){
            dv = (Gxyz[Gc_AN][j+29] - 0.5*Gxyz[Gc_AN][j+23]*dt*dzeta)/(NH_czeta*0.5*dt + 1.0); 
            Gxyz[Gc_AN][j+23] += dv;
            My_sum += dv*dv; 
	  }
        }
      }

      NH_czeta += dzeta; 

      MPI_Allreduce(&My_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

      sum += dzeta*dzeta;

      if (sum<1.0e-12) po = 1;

      num++; 

      if (20<num) po = 1;

    } while(po==0);

    /****************************************************
                       Kinetic Energy 
    ****************************************************/

    Ukc = 0.0;
    My_Ukc = 0.0;

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];

      sum = 0.0;
      for (j=1; j<=3; j++){
        if (atom_Fixed_XYZ[Gc_AN][j]==0){
          sum = sum + Gxyz[Gc_AN][j+23]*Gxyz[Gc_AN][j+23];
	}
      }
      My_Ukc += 0.5*Gxyz[Gc_AN][20]*Wscale*sum;
    }

    /****************************************************
     MPI: Ukc 
    ****************************************************/

    MPI_Allreduce(&My_Ukc, &Ukc, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

    /* calculation of temperature (K) */

    Temp = Ukc/(1.5*kB*(double)atomnum)*eV2Hartree;

    /****************************************************
     write informatins to *.ene
    ****************************************************/

    if (myid==Host_ID){  

      sprintf(fileE,"%s%s.ene",filepath,filename);
      iterout_md(iter,MD_TimeStep*(iter-1),fileE);
    }

    /****************************************************
     Nose-Hoover Hamiltonian which is a conserved quantity
    ****************************************************/

    NH_Ham = Utot + Ukc + 0.5*NH_czeta*NH_czeta*TempQ*Wscale
                        + 3.0*kB*(double)atomnum*GivenTemp*NH_R/eV2Hartree; 

    /*****************************************************
     first step:

       v(t)    ->  v(t+0.5*dt) 
       r(t)    ->  r(t+dt) 
       zeta(t) ->  zeta(t+0.5*dt)
       R(t)    ->  R(t+dt) 
    *****************************************************/

    if (MD_Opt_OK!=1 && iter!=MD_IterNumber){

      for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
	Gc_AN = M2G[Mc_AN];
	for (j=1; j<=3; j++){

	  if (atom_Fixed_XYZ[Gc_AN][j]==0){
        
	    scaled_force = -Gxyz[Gc_AN][16+j]/(Gxyz[Gc_AN][20]*Wscale);  

	    /* v( r+0.5*dt ) */
	    Gxyz[Gc_AN][26+j] = Gxyz[Gc_AN][23+j] + (scaled_force - NH_czeta*Gxyz[Gc_AN][23+j])*0.5*dt;
	    Gxyz[Gc_AN][23+j] = Gxyz[Gc_AN][26+j];

	    /* r( r+dt ) */
	    Gxyz[Gc_AN][20+j] = Gxyz[Gc_AN][j];
	    Gxyz[Gc_AN][j] =  Gxyz[Gc_AN][j] + Gxyz[Gc_AN][26+j]*dt;
	  }
	}
      }

      /* zeta( t+0.5*dt ) */

      NH_nzeta = NH_czeta + (Ukc - 1.5*kB*(double)atomnum*GivenTemp/eV2Hartree)*dt/(TempQ*Wscale);
      NH_czeta = NH_nzeta;

      /* R( r+dt ) */
      NH_R = NH_R + NH_nzeta*dt;

    }

  }

  else {

    /****************************************************
      second step in velocity Verlet 
    ****************************************************/

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];
      for (j=1; j<=3; j++){

        if (atom_Fixed_XYZ[Gc_AN][j]==0){
          scaled_force = -Gxyz[Gc_AN][16+j]/(Gxyz[Gc_AN][20]*Wscale);  
          Gxyz[Gc_AN][23+j] = Gxyz[Gc_AN][26+j] + scaled_force*0.5*dt;
	}
      }
    }

    /****************************************************
                       Kinetic Energy 
    ****************************************************/

    My_Ukc = 0.0;

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];

      sum = 0.0;
      for (j=1; j<=3; j++){
        if (atom_Fixed_XYZ[Gc_AN][j]==0){
          sum = sum + Gxyz[Gc_AN][j+23]*Gxyz[Gc_AN][j+23];
	}
      }
      My_Ukc = My_Ukc + 0.5*Gxyz[Gc_AN][20]*Wscale*sum;
    }

    /****************************************************
     MPI: Ukc 
    ****************************************************/

    MPI_Allreduce(&My_Ukc, &Ukc, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

    /* calculation of temperature (K) */
    Temp = Ukc/(1.5*kB*(double)atomnum)*eV2Hartree;

    /****************************************************
     write informatins to *.ene
    ****************************************************/

    if (myid==Host_ID){  

      sprintf(fileE,"%s%s.ene",filepath,filename);
      iterout_md(iter,MD_TimeStep*(iter-1),fileE);
    }

    /****************************************************
      first step in velocity Verlet 
    ****************************************************/

    if (MD_Opt_OK!=1 && iter!=MD_IterNumber){

      for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
	Gc_AN = M2G[Mc_AN];
	for (j=1; j<=3; j++){

	  if (atom_Fixed_XYZ[Gc_AN][j]==0){

	    scaled_force = -Gxyz[Gc_AN][16+j]/(Gxyz[Gc_AN][20]*Wscale);  
	    /* v( r+0.5*dt ) */
	    Gxyz[Gc_AN][26+j] = Gxyz[Gc_AN][23+j] + scaled_force*0.5*dt;
	    Gxyz[Gc_AN][23+j] = Gxyz[Gc_AN][26+j];

	    /* r( r+dt ) */
	    Gxyz[Gc_AN][20+j] = Gxyz[Gc_AN][j];
	    Gxyz[Gc_AN][j] =  Gxyz[Gc_AN][j] + Gxyz[Gc_AN][26+j]*dt;
	  }
	}
      }
    }

  }

  /****************************************************
   MPI: Gxyz
  ****************************************************/

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    ID = G2ID[Gc_AN];
    MPI_Bcast(&Gxyz[Gc_AN][1],   3, MPI_DOUBLE, ID, mpi_comm_level1);
    MPI_Bcast(&Gxyz[Gc_AN][14], 19, MPI_DOUBLE, ID, mpi_comm_level1);
  }
}





static void EvsLC(int iter)
{
  int i,j,Gc_AN,myid;
  char fileE[YOUSO10];
  static double tv0[4][4];
  static int firsttime=1;

  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
   store tv into tv0 
  ****************************************************/

  if (firsttime){

    for (i=1; i<=3; i++){
      for (j=1; j<=3; j++){
	tv0[i][j] = tv[i][j];
      }
    }

    firsttime = 0;
  }

  /****************************************************
   write informatins to *.ene
  ****************************************************/

  if (myid==Host_ID){  
    sprintf(fileE,"%s%s.ene",filepath,filename);
    iterout_md(iter,MD_TimeStep*(iter-1),fileE);
  }

  /****************************************************
    scale tv 
  ****************************************************/

  for (i=1; i<=3; i++){
    for (j=1; j<=3; j++){
      tv[i][j] += tv0[i][j]*MD_EvsLattice_Step/100.0;
    }
  }

  /****************************************************
    update cartesian coordinates
  ****************************************************/

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){

    Gxyz[Gc_AN][1] =  Cell_Gxyz[Gc_AN][1]*tv[1][1]
                    + Cell_Gxyz[Gc_AN][2]*tv[2][1]
                    + Cell_Gxyz[Gc_AN][3]*tv[3][1] + Grid_Origin[1];

    Gxyz[Gc_AN][2] =  Cell_Gxyz[Gc_AN][1]*tv[1][2]
                    + Cell_Gxyz[Gc_AN][2]*tv[2][2]
                    + Cell_Gxyz[Gc_AN][3]*tv[3][2] + Grid_Origin[2];

    Gxyz[Gc_AN][3] =  Cell_Gxyz[Gc_AN][1]*tv[1][3]
                    + Cell_Gxyz[Gc_AN][2]*tv[2][3]
                    + Cell_Gxyz[Gc_AN][3]*tv[3][3] + Grid_Origin[3];
  }

}


