/**********************************************************************
  MD_pac.c:

     MD_pac.c is a subroutine to perform molecular dynamics
     simulations and geometry optimization.

  Log of MD_pac.c:

     22/Nov/2001  Released by T. Ozaki
     15/Dec/2003  DIISs are added by H. Kino
     14/May/2004  NVT_VS is added by M. Ohfuti
     25/May/2004  Modified by T. Ozaki
     14/Jul/2007  RF added by H.M. Weng
     08/Jan/2010  NVT_VS2 added by T. Ohwaki 
     23/Dec/2012  RestartFiles4GeoOpt added by T. Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "openmx_common.h"
#include "lapack_prototypes.h"
#include "mpi.h"


#define Criterion_Max_Step            0.06

 
static void NoMD(int iter);
static void VerletXYZ(int iter);
static void NVT_VS(int iter);  /* added by mari */
static void NVT_VS2(int iter); /* added by Ohwaki */
static void NVT_VS4(int iter); /* added by Ohwaki */
static void NVT_NH(int iter); 
static void NVT_Langevin(int iter); /* added by Ohwaki */
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
static void Delta_Factor(int iter);
static void Correct_Force();
static int RestartFiles4GeoOpt(char *mode);

int SuccessReadingfiles;




double MD_pac(int iter, char *fname_input)
{
  double time0;
  double TStime,TEtime;
  int numprocs,myid;

  dtime(&TStime);
 
  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (myid==Host_ID && 0<level_stdout){
    printf("\n*******************************************************\n"); 
    printf("             MD or geometry opt. at MD =%2d              \n",iter);
    printf("*******************************************************\n\n"); 
  }

  /*********************************************** 
    read restart files for geometry optimization
    MD_switch==3: Steepest_Descent
    MD_switch==4: Geometry_Opt_DIIS_EF
    MD_switch==5: Geometry_Opt_DIIS_BFGS
    MD_switch==6: Geometry_Opt_RF
    MD_switch==7: Geometry_Opt_DIIS
  ***********************************************/

  if (iter==1 && GeoOpt_RestartFromFile){
    if ( MD_switch==3 || MD_switch==4 || MD_switch==5 || MD_switch==6 || MD_switch==7 ){
      if (RestartFiles4GeoOpt("read")){
        SuccessReadingfiles = 1;
        if (myid==Host_ID){
          printf("Found restart files for geometry optimization\n");
        }
      }
    }
  }
  else {
    SuccessReadingfiles = 0;
  }

  /* correct forces for MD simulations so that the sum of forces can be zero. */

  if (  MD_switch==1 || MD_switch==2 || MD_switch==9 
     || MD_switch==11 || MD_switch==14 || MD_switch==15 ) Correct_Force();

  /* Call a subroutine based on MD_switch */

  switch (MD_switch) {
    case  0: NoMD(iter);                    break;
    case  1: VerletXYZ(iter);               break;
    case  2: NVT_VS(iter);                  break;  /* added by mari */
    case  3: Steepest_Descent(iter+MD_Current_Iter,1);      break;
    case  4: Geometry_Opt_DIIS_EF(iter+MD_Current_Iter);    break;
    case  5: Geometry_Opt_DIIS_BFGS(iter+MD_Current_Iter);  break;
    case  6: Geometry_Opt_RF(iter+MD_Current_Iter);         break;  /* added by hmweng */
    case  7: Geometry_Opt_DIIS(iter+MD_Current_Iter);       break;
    case  8:                                break;  /* not used */
    case  9: NVT_NH(iter);                  break;
    case 10:                                break;  /* not used */
    case 11: NVT_VS2(iter);                 break;  /* added by Ohwaki */
    case 12: EvsLC(iter);                   break;
    case 14: NVT_VS4(iter);                 break;  /* added by Ohwaki */
    case 15: NVT_Langevin(iter);            break;  /* added by Ohwaki */
    case 16: Delta_Factor(iter);            break;  /* delta-factor */
  }

  /***************************************************************
    correct atoms which are out of the first unit cell during 
    molecular dynamics simulations. The correction is not applied 
    for geometry optimization.
  ***************************************************************/

  if (   MD_switch==1 ||
         MD_switch==2 ||
         MD_switch==9 ||
         MD_switch==11||
         MD_switch==14||
         MD_switch==15 ){

    Correct_Position_In_First_Cell();
  }

  /* making of an input file with the final structure */
  if (Runtest_flag==0){
    Make_InputFile_with_FinalCoord(fname_input,iter);
  }

  /*********************************************** 
    save restart files for geometry optimization
    MD_switch==4: Geometry_Opt_DIIS_EF
    MD_switch==5: Geometry_Opt_DIIS_BFGS
    MD_switch==6: Geometry_Opt_RF
    MD_switch==7: Geometry_Opt_DIIS
  ***********************************************/

  if ( MD_switch==3 || MD_switch==4 || MD_switch==5 || MD_switch==6 || MD_switch==7 ){
    RestartFiles4GeoOpt("write");
  }


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
    iterout_md(iter+MD_Current_Iter,MD_TimeStep*(iter+MD_Current_Iter-1),fileE);
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
  Wscale = unified_atomic_mass_unit/electron_mass;

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
      iterout_md(iter+MD_Current_Iter,MD_TimeStep*(iter+MD_Current_Iter-1),fileE);
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
      iterout_md(iter+MD_Current_Iter,MD_TimeStep*(iter+MD_Current_Iter-1),fileE);
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
  double My_Max_Force,Wscale;
  char fileCoord[YOUSO10];
  char fileSD[YOUSO10];
  FILE *fp_crd,*fp_SD;
  int numprocs,myid,ID;
  double tmp1,MaxStep;
  char buf[fp_bsize];          /* setvbuf */
  char fileE[YOUSO10];
  char file_name[YOUSO10];
  FILE *fp;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  MD_Opt_OK = 0;
  Wscale = unified_atomic_mass_unit/electron_mass;

  /******************************************
              read *.rst4gopt.SD1
  ******************************************/

  if (SuccessReadingfiles){

    sprintf(file_name,"%s%s_rst/%s.rst4gopt.SD1",filepath,filename,filename);

    if ((fp = fopen(file_name,"rb")) != NULL){

      fread(&Correct_Position_flag, sizeof(int), 1, fp);
      fread(&SD_scaling,      sizeof(double), 1, fp);
      fread(&SD_scaling_user, sizeof(double), 1, fp);
      fread(&Past_Utot,       sizeof(double), 10, fp);

      fclose(fp);
    }    
    else{
      printf("Failure of saving %s\n",file_name);
    }
  }

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
  SD_init = dt*dt/Wscale;
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
    iterout_md(iter+MD_Current_Iter,MD_TimeStep*(iter+MD_Current_Iter-1),fileE);
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

  /******************************************
              save *.rst4gopt.SD1 
  ******************************************/

  if (myid==Host_ID){

    sprintf(file_name,"%s%s_rst/%s.rst4gopt.SD1",filepath,filename,filename);

    if ((fp = fopen(file_name,"wb")) != NULL){

      fwrite(&Correct_Position_flag, sizeof(int), 1, fp);
      fwrite(&SD_scaling,      sizeof(double), 1, fp);
      fwrite(&SD_scaling_user, sizeof(double), 1, fp);
      fwrite(&Past_Utot,       sizeof(double), 10, fp);

      fclose(fp);
    }    
    else{
      printf("Failure of saving %s\n",file_name);
    }
  }

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
  char file_name[YOUSO10];
  FILE *fp;
  int tmp_array[10];
  int everyiter,buf_iter;
  int numprocs,myid,ID;  

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  /******************************************
              read *.rst4gopt.RF1
  ******************************************/

  if (SuccessReadingfiles){

    sprintf(file_name,"%s%s_rst/%s.rst4gopt.RF1",filepath,filename,filename);

    if ((fp = fopen(file_name,"rb")) != NULL){

      fread(tmp_array, sizeof(int), 6, fp);

      local_iter            = tmp_array[0];
      SD_iter               = tmp_array[1];
      GDIIS_iter            = tmp_array[2];
      flag                  = tmp_array[3];
      Every_iter            = tmp_array[4];
      Correct_Position_flag = tmp_array[5];
 
      fread(&SD_scaling,      sizeof(double), 1, fp);
      fread(&SD_scaling_user, sizeof(double), 1, fp);
      fread(&Past_Utot,       sizeof(double), 10, fp);

      fclose(fp);
    }    
    else{
      printf("Failure of saving %s\n",file_name);
    }
  }

  /* set parameters */

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

  /******************************************
              save *.rst4gopt.RF1 
  ******************************************/

  if (myid==Host_ID){

    sprintf(file_name,"%s%s_rst/%s.rst4gopt.RF1",filepath,filename,filename);

    if ((fp = fopen(file_name,"wb")) != NULL){

      tmp_array[0] = local_iter;
      tmp_array[1] = SD_iter;
      tmp_array[2] = GDIIS_iter;
      tmp_array[3] = flag;
      tmp_array[4] = Every_iter;
      tmp_array[5] = Correct_Position_flag;
 
      fwrite(tmp_array, sizeof(int), 6, fp);

      fwrite(&SD_scaling,      sizeof(double), 1, fp);
      fwrite(&SD_scaling_user, sizeof(double), 1, fp);
      fwrite(&Past_Utot,       sizeof(double), 10, fp);

      fclose(fp);
    }    
    else{
      printf("Failure of saving %s\n",file_name);
    }
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
  char file_name[YOUSO10];
  FILE *fp;
  int tmp_array[10];
  int everyiter,buf_iter;
  int numprocs,myid,ID;  

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  /******************************************
              read *.rst4gopt.DIIS1
  ******************************************/

  if (SuccessReadingfiles){

    sprintf(file_name,"%s%s_rst/%s.rst4gopt.DIIS1",filepath,filename,filename);

    if ((fp = fopen(file_name,"rb")) != NULL){

      fread(tmp_array, sizeof(int), 6, fp);

      local_iter            = tmp_array[0];
      SD_iter               = tmp_array[1];
      GDIIS_iter            = tmp_array[2];
      flag                  = tmp_array[3];
      Every_iter            = tmp_array[4];
      Correct_Position_flag = tmp_array[5];
 
      fread(&SD_scaling,      sizeof(double), 1, fp);
      fread(&SD_scaling_user, sizeof(double), 1, fp);
      fread(&Past_Utot,       sizeof(double), 10, fp);

      fclose(fp);
    }    
    else{
      printf("Failure of saving %s\n",file_name);
    }
  }

  /* set parameters */

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

  /******************************************
              save *.rst4gopt.DIIS1 
  ******************************************/

  if (myid==Host_ID){

    sprintf(file_name,"%s%s_rst/%s.rst4gopt.DIIS1",filepath,filename,filename);

    if ((fp = fopen(file_name,"wb")) != NULL){

      tmp_array[0] = local_iter;
      tmp_array[1] = SD_iter;
      tmp_array[2] = GDIIS_iter;
      tmp_array[3] = flag;
      tmp_array[4] = Every_iter;
      tmp_array[5] = Correct_Position_flag;
 
      fwrite(tmp_array, sizeof(int), 6, fp);

      fwrite(&SD_scaling,      sizeof(double), 1, fp);
      fwrite(&SD_scaling_user, sizeof(double), 1, fp);
      fwrite(&Past_Utot,       sizeof(double), 10, fp);

      fclose(fp);
    }    
    else{
      printf("Failure of saving %s\n",file_name);
    }
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
  char file_name[YOUSO10];
  FILE *fp;
  int tmp_array[10];
  int everyiter,buf_iter;
  int numprocs,myid,ID;  

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  /******************************************
              read *.rst4gopt.BFGS1
  ******************************************/

  if (SuccessReadingfiles){

    sprintf(file_name,"%s%s_rst/%s.rst4gopt.BFGS1",filepath,filename,filename);

    if ((fp = fopen(file_name,"rb")) != NULL){

      fread(tmp_array, sizeof(int), 6, fp);

      local_iter            = tmp_array[0];
      SD_iter               = tmp_array[1];
      GDIIS_iter            = tmp_array[2];
      flag                  = tmp_array[3];
      Every_iter            = tmp_array[4];
      Correct_Position_flag = tmp_array[5];
 
      fread(&SD_scaling,      sizeof(double), 1, fp);
      fread(&SD_scaling_user, sizeof(double), 1, fp);
      fread(&Past_Utot,       sizeof(double), 10, fp);

      fclose(fp);
    }    
    else{
      printf("Failure of saving %s\n",file_name);
    }
  }

  /* set parameters */

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

  /******************************************
              save *.rst4gopt.BFGS1 
  ******************************************/

  if (myid==Host_ID){

    sprintf(file_name,"%s%s_rst/%s.rst4gopt.BFGS1",filepath,filename,filename);

    if ((fp = fopen(file_name,"wb")) != NULL){

      tmp_array[0] = local_iter;
      tmp_array[1] = SD_iter;
      tmp_array[2] = GDIIS_iter;
      tmp_array[3] = flag;
      tmp_array[4] = Every_iter;
      tmp_array[5] = Correct_Position_flag;
 
      fwrite(tmp_array, sizeof(int), 6, fp);

      fwrite(&SD_scaling,      sizeof(double), 1, fp);
      fwrite(&SD_scaling_user, sizeof(double), 1, fp);
      fwrite(&Past_Utot,       sizeof(double), 10, fp);

      fclose(fp);
    }    
    else{
      printf("Failure of saving %s\n",file_name);
    }
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
  char file_name[YOUSO10];
  FILE *fp;
  int tmp_array[10];
  int everyiter,buf_iter;
  int numprocs,myid,ID;  

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  /******************************************
              read *.rst4gopt.EF1
  ******************************************/

  if (SuccessReadingfiles){

    sprintf(file_name,"%s%s_rst/%s.rst4gopt.EF1",filepath,filename,filename);

    if ((fp = fopen(file_name,"rb")) != NULL){

      fread(tmp_array, sizeof(int), 6, fp);

      local_iter            = tmp_array[0];
      SD_iter               = tmp_array[1];
      GDIIS_iter            = tmp_array[2];
      flag                  = tmp_array[3];
      Every_iter            = tmp_array[4];
      Correct_Position_flag = tmp_array[5];
 
      fread(&SD_scaling,      sizeof(double), 1, fp);
      fread(&SD_scaling_user, sizeof(double), 1, fp);
      fread(&Past_Utot,       sizeof(double), 10, fp);

      fclose(fp);
    }    
    else{
      printf("Failure of saving %s\n",file_name);
    }
  }

  /* set parameters */

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

  /******************************************
              save *.rst4gopt.EF1 
  ******************************************/

  if (myid==Host_ID){

    sprintf(file_name,"%s%s_rst/%s.rst4gopt.EF1",filepath,filename,filename);

    if ((fp = fopen(file_name,"wb")) != NULL){

      tmp_array[0] = local_iter;
      tmp_array[1] = SD_iter;
      tmp_array[2] = GDIIS_iter;
      tmp_array[3] = flag;
      tmp_array[4] = Every_iter;
      tmp_array[5] = Correct_Position_flag;
 
      fwrite(tmp_array, sizeof(int), 6, fp);

      fwrite(&SD_scaling,      sizeof(double), 1, fp);
      fwrite(&SD_scaling_user, sizeof(double), 1, fp);
      fwrite(&Past_Utot,       sizeof(double), 10, fp);

      fclose(fp);
    }    
    else{
      printf("Failure of saving %s\n",file_name);
    }
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
    iterout_md(iter+MD_Current_Iter,MD_TimeStep*(iter+MD_Current_Iter-1),fileE);
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
    iterout_md(iter+MD_Current_Iter,MD_TimeStep*(iter+MD_Current_Iter-1),fileE);
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

    /*********************************************************************************
     update an approximate Hessian matrix 
   
     H_k = H_{k-1} + y*y^t/(y^t*y) - H_{k-1}*s * s^t * H_{k-1}^t /(s^t * H_{k-1} * s)
     y = GxyzHistoryR[0][iatom][k]  - GxyzHistoryR[1][iatom][k]
     s = GxyzHistoryIn[0][iatom][k] - GxyzHistoryIn[1][iatom][k]
    *********************************************************************************/

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

    /*
    for(i=1; i<=1; i++){
      printf("Eigenvalue=%15.8f\n",work2[i]);
    }
    printf("EigenVector is\n");
    */

    for(i=1;i<=3*atomnum+1; i++){

      /*
      printf("%15.8f\n",ahes[i][1]);
      */

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
    iterout_md(iter+MD_Current_Iter,MD_TimeStep*(iter+MD_Current_Iter-1),fileE);
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
  double tmp_array[10];
  char file_name[YOUSO10];
  FILE *fp;

  /* variables for MPI */
  int Gc_AN;
  int numprocs,myid,ID;  

  /* MPI myid */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  /******************************************
              read *.rst4gopt.EF2 
  ******************************************/

  if (SuccessReadingfiles){

    sprintf(file_name,"%s%s_rst/%s.rst4gopt.EF2",filepath,filename,filename);

    if ((fp = fopen(file_name,"rb")) != NULL){

      fread(tmp_array, sizeof(double), 2, fp);

      Utot0          = tmp_array[0];
      scaling_factor = tmp_array[1];
 
      fclose(fp);
    }    
    else{
      printf("Failure of saving %s\n",file_name);
    }
  }

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
    iterout_md(iter+MD_Current_Iter,MD_TimeStep*(iter+MD_Current_Iter-1),fileE);
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

          printf("ABC1 iatom=%2d j=%2d i=%2d %15.12f %15.12f %15.12f\n",
                 iatom,j,i,GxyzHistoryIn[i][iatom][j],B[i],Gxyz[iatom][j]);
	}


	/* with a quasi Newton method */ 

	m1++;

	Gxyz[iatom][j] -= scaling_factor*U[0][m1];

        printf("ABC2 iatom=%2d j=%2d %15.12f %15.12f %15.12f\n",
	       iatom,j,scaling_factor,U[0][m1],Gxyz[iatom][j]);

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

  /******************************************
              save *.rst4gopt.EF2 
  ******************************************/

  if (myid==Host_ID){

    sprintf(file_name,"%s%s_rst/%s.rst4gopt.EF2",filepath,filename,filename);

    if ((fp = fopen(file_name,"wb")) != NULL){

      tmp_array[0] = Utot0;
      tmp_array[1] = scaling_factor;
 
      fwrite(tmp_array, sizeof(double), 2, fp);
      fclose(fp);
    }    
    else{
      printf("Failure of saving %s\n",file_name);
    }
  }

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
  Wscale = unified_atomic_mass_unit/electron_mass;

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
    iterout_md(iter+MD_Current_Iter,MD_TimeStep*(iter+MD_Current_Iter-1),fileE);
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




void NVT_Langevin(int iter)
{
  /* added by T.Ohwaki */
  /********************************************************
   This routine is added by Tsuruku Ohwaki (May 2012).

   a constant temperature molecular dynamics by Langevin
   heat-bath method with velocity-Verlet integrator
  ********************************************************/
  /*******************************************************
   1 a.u.=2.4189*10^-2 fs, 1fs=41.341105 a.u.
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

    Gxyz[][24] = x-component of velocity at t
    Gxyz[][25] = y-component of velocity at t
    Gxyz[][26] = z-component of velocity at t

    Gxyz[][27] = x-component of tilde velocity at t+dt
    Gxyz[][28] = y-component of tilde velocity at t+dt
    Gxyz[][29] = z-component of tilde velocity at t+dt

    Gxyz[][30] = hx
    Gxyz[][31] = hy
    Gxyz[][32] = hz

  ****************************************************/

  double dt,dt2,sum,My_Ukc,x,t;
  double Wscale,rt,rt_mdt,rt_pdt,vt,vt_mdt,vt_pdt;
  double rtmp1,rtmp2,tmp,Lang_sig,RandomF;
  double vtt,ft_mdt,ftt,ft,vtt_pdt;
  double tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
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
  Wscale = unified_atomic_mass_unit/electron_mass;

  /* calculation of a given temperature (K) */

  i = 1;
  do {
    if ( TempPara[i][1]<=iter && iter<TempPara[i+1][1] ){

      GivenTemp = TempPara[i][2] + (TempPara[i+1][2] - TempPara[i][2])*
                 ((double)iter - (double)TempPara[i][1])
                /((double)TempPara[i+1][1] - (double)TempPara[i][1]);

    }
    i++;
  } while (i<=(TempNum-1));

  /****************************************************
   add random forces on Gxyz[17-19]
  ****************************************************/

  /*
  srand((unsigned)time(NULL));
  */

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];
    Lang_sig = sqrt(2.0*Gxyz[Gc_AN][20]*Wscale*FricFac*kB*GivenTemp/dt/eV2Hartree);

    for (j=1; j<=3; j++){

      if (atom_Fixed_XYZ[Gc_AN][j]==0){

	rtmp1 = (double)rand()/RAND_MAX;
	rtmp2 = (double)rand()/RAND_MAX;

	RandomF = Lang_sig * sqrt(-2.0*log(rtmp1)) * cos(2.0*PI*rtmp2);
	Gxyz[Gc_AN][16+j] -= RandomF;

      }
    }
  }

  /****************************************************
                Velocity Verlet algorithm
  ****************************************************/

  if (iter==1){

    /****************************************************
      first step in velocity Verlet 
    ****************************************************/

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];

      for (j=1; j<=3; j++){

        if (atom_Fixed_XYZ[Gc_AN][j]==0){

          /* update coordinates */
          Gxyz[Gc_AN][20+j] = Gxyz[Gc_AN][j];

          vt = Gxyz[Gc_AN][23+j];
          vt_mdt = Gxyz[Gc_AN][26+j];

          Gxyz[Gc_AN][j] = Gxyz[Gc_AN][j] + dt*vt
                        -0.5*dt2*(Gxyz[Gc_AN][16+j]/(Gxyz[Gc_AN][20]*Wscale)+FricFac*vt);

          /* store current gradient */
          Gxyz[Gc_AN][13+j] = Gxyz[Gc_AN][16+j];
        }
      }
    }
  }

  else{

    /****************************************************
      for the second step and then onward 
      in velocity Verlet 
    ****************************************************/

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

      Gc_AN = M2G[Mc_AN];

      tmp = 1.0/(1.0 + 0.50*dt*FricFac);
      tmp2 = 1.0 - 0.50*dt*FricFac;

      for (j=1; j<=3; j++){

        if (atom_Fixed_XYZ[Gc_AN][j]==0){

          /* modified Euler */

	  if (0){
          rt = Gxyz[Gc_AN][j];
          vtt = Gxyz[Gc_AN][26+j];
          vt_mdt = Gxyz[Gc_AN][23+j];
          ft_mdt = Gxyz[Gc_AN][13+j];

          ftt = -Gxyz[Gc_AN][16+j]/(Gxyz[Gc_AN][20]*Wscale) - FricFac*vtt;
          vt = vt_mdt + 0.5*dt*(ft_mdt + ftt);

          ft = -Gxyz[Gc_AN][16+j]/(Gxyz[Gc_AN][20]*Wscale) - FricFac*vt;

          vtt_pdt = vt + dt*ft;
          rt_pdt = rt + 0.5*dt*(vt + vtt_pdt); 

          Gxyz[Gc_AN][j] = rt_pdt;
          Gxyz[Gc_AN][26+j] = vtt_pdt;
          Gxyz[Gc_AN][23+j] = vt;
          Gxyz[Gc_AN][13+j] = ft;
	  }

          /* Taylor expansion */

          if (0){

          vt = Gxyz[Gc_AN][23+j];
          vt_mdt = Gxyz[Gc_AN][26+j];

          Gxyz[Gc_AN][j] = Gxyz[Gc_AN][j] + dt*vt
                        -0.5*dt2*(Gxyz[Gc_AN][16+j]/(Gxyz[Gc_AN][20]*Wscale)+FricFac*vt);

          Gxyz[Gc_AN][23+j] = vt_mdt - 2.0*dt*Gxyz[Gc_AN][16+j]/(Gxyz[Gc_AN][20]*Wscale) - 2.0*dt*FricFac*vt;
          Gxyz[Gc_AN][26+j] = vt;

	  /*
          printf("ABC iter=%2d Gc_AN=%2d j=%2d Gxyz[Gc_AN][23+j]=%15.12f\n",iter,Gc_AN,j,Gxyz[Gc_AN][23+j]); 
	  */
	  }

          /* simple Taylor expansion */
	  if (0){
          vt = Gxyz[Gc_AN][23+j];

          Gxyz[Gc_AN][j] = Gxyz[Gc_AN][j] + dt*vt
                        -0.5*dt2*(Gxyz[Gc_AN][16+j]/(Gxyz[Gc_AN][20]*Wscale)+FricFac*vt);

          Gxyz[Gc_AN][23+j] = vt - dt*Gxyz[Gc_AN][16+j]/(Gxyz[Gc_AN][20]*Wscale) - dt*FricFac*vt;

	  /*
          printf("ABC iter=%2d Gc_AN=%2d j=%2d Gxyz[Gc_AN][23+j]=%15.12f\n",iter,Gc_AN,j,Gxyz[Gc_AN][23+j]); 
	  */
	  }


          /* Brunger-Brooks-Karplus integrator */ 
	  if (0){

          /* update coordinates */

          rt_mdt = Gxyz[Gc_AN][20+j];
          rt = Gxyz[Gc_AN][j];

          rt_pdt = 2.0*rt - tmp2*rt_mdt - dt2/(Gxyz[Gc_AN][20]*Wscale)*Gxyz[Gc_AN][16+j]; 
          rt_pdt *= tmp; 

          /* update velocity */

          vt = 0.5*(rt_pdt - rt_mdt)/dt;

          /* store data */

          Gxyz[Gc_AN][20+j] = Gxyz[Gc_AN][j];  /* store coordinate at t    */
          Gxyz[Gc_AN][j] = rt_pdt;             /* store coordinate at t+dt */
          Gxyz[Gc_AN][23+j] = vt;              /* store velocity at t      */
	  }


          /* Velocity Verlet formulation of the BBK integrator */
	  if (0){

          /* velocity at the current step */
          Gxyz[Gc_AN][23+j] = tmp2*Gxyz[Gc_AN][23+j] 
                            -(Gxyz[Gc_AN][16+j] + Gxyz[Gc_AN][13+j])
                            /(Gxyz[Gc_AN][20]*Wscale)*dt*0.50;

          Gxyz[Gc_AN][23+j] *= tmp;

          /* update coordinates */
          Gxyz[Gc_AN][20+j] = Gxyz[Gc_AN][j];
          Gxyz[Gc_AN][j] = Gxyz[Gc_AN][j] + dt*Gxyz[Gc_AN][23+j]*tmp2
                          -Gxyz[Gc_AN][16+j]/(Gxyz[Gc_AN][20]*Wscale)*dt2*0.50;

          /* store current gradient */
          Gxyz[Gc_AN][13+j] = Gxyz[Gc_AN][16+j];
	  }


          /* Velocity Verlet formulation of the BBK integrator with velocity correction by Dr. Ohwaki */
	  if (0){

          /* velocity at the current step */

          tmp1 = Gxyz[Gc_AN][23+j];

          Gxyz[Gc_AN][23+j] = Gxyz[Gc_AN][23+j] 
                            -(Gxyz[Gc_AN][16+j] + Gxyz[Gc_AN][13+j])
                            /(Gxyz[Gc_AN][20]*Wscale)*dt*0.50;

          Gxyz[Gc_AN][23+j] *= tmp;

          /* add contribution of friction force */
          Gxyz[Gc_AN][16+j] += Gxyz[Gc_AN][20]*Wscale*FricFac*Gxyz[Gc_AN][23+j];


          /* update coordinates */
          Gxyz[Gc_AN][20+j] = Gxyz[Gc_AN][j];
          Gxyz[Gc_AN][j] = Gxyz[Gc_AN][j] + dt*Gxyz[Gc_AN][23+j]
                          -Gxyz[Gc_AN][16+j]/(Gxyz[Gc_AN][20]*Wscale)*dt2*0.50;

          /* store current gradient */
          Gxyz[Gc_AN][13+j] = Gxyz[Gc_AN][16+j]
                             -Gxyz[Gc_AN][20]*Wscale*FricFac*Gxyz[Gc_AN][23+j]
                             +Gxyz[Gc_AN][20]*Wscale*FricFac*tmp1;
	  }

          /* Velocity Verlet formulation by Ozaki's integrator */

	  if (1){

          /* velocity at the current step */
          Gxyz[Gc_AN][23+j] = Gxyz[Gc_AN][23+j] 
                            -(Gxyz[Gc_AN][16+j] + Gxyz[Gc_AN][13+j])/(Gxyz[Gc_AN][20]*Wscale)*dt*0.50 
                            -FricFac*(Gxyz[Gc_AN][j] - Gxyz[Gc_AN][20+j]);

          /* update coordinates */
          Gxyz[Gc_AN][20+j] = Gxyz[Gc_AN][j];
          Gxyz[Gc_AN][j] = Gxyz[Gc_AN][j]
                        + tmp*(dt*Gxyz[Gc_AN][23+j]-Gxyz[Gc_AN][16+j]/(Gxyz[Gc_AN][20]*Wscale)*dt2*0.50);

          /* store current gradient */
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
    MPI_Bcast(&Gxyz[Gc_AN][0], 40, MPI_DOUBLE, ID, mpi_comm_level1);
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
    iterout_md(iter+MD_Current_Iter,MD_TimeStep*(iter+MD_Current_Iter-1),fileE);
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
  Wscale = unified_atomic_mass_unit/electron_mass;

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
    iterout_md(iter+MD_Current_Iter,MD_TimeStep*(iter+MD_Current_Iter-1),fileE);
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


void NVT_VS4(int iter)
{
  /* added by T.Ohwaki */
  /********************************************************
   This routine is added by T.Ohwaki (2011/11/11).

   a constant temperature molecular dynamics
   by a velocity scaling method for each atom group
   with velocity-Verlet integrator
  ********************************************************/
  /*******************************************************
   1 a.u.=2.4189*10^-2 fs, 1fs=41.341105 a.u.
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

  double dt,dt2,sum,My_Ukc,x,t;
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
  Wscale = unified_atomic_mass_unit/electron_mass;

  double sum_AtGr[num_AtGr+1],My_Ukc_AtGr[num_AtGr+1],Ukc_AtGr[num_AtGr+1];
  double AtomGr_Scale[num_AtGr+1];

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

  /******************************************************
      Kinetic Energy & Temperature (for each atom group)
  ******************************************************/

  for (k=1; k<=num_AtGr; k++){ 
    sum_AtGr[k] = 0.0; 
  }

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    k = AtomGr[Gc_AN];

    for (j=1; j<=3; j++){
      sum_AtGr[k] += Gxyz[Gc_AN][j+23]*Gxyz[Gc_AN][j+23]*Gxyz[Gc_AN][20];
    }
  }

  for (k=1; k<=num_AtGr; k++){
    My_Ukc_AtGr[k]  = 0.5*Wscale*sum_AtGr[k];
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
   MPI: Ukc_AtGr
  ****************************************************/

  for (k=1; k<=num_AtGr; k++){
    MPI_Allreduce(&My_Ukc_AtGr[k], &Ukc_AtGr[k], 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

    /* calculation of temperature (K) */
    Temp_AtGr[k] = Ukc_AtGr[k]/(1.5*kB*(double)atnum_AtGr[k])*eV2Hartree;
  }

  /****************************************************
   write informatins to *.ene
  ****************************************************/

  if (myid==Host_ID){
    sprintf(fileE,"%s%s.ene",filepath,filename);
    iterout_md(iter+MD_Current_Iter,MD_TimeStep*(iter+MD_Current_Iter-1),fileE);
  }

  /****************************************************
   velocity scaling for each atom group
  ****************************************************/

  if(iter!=1) {

  for (k=1; k<=num_AtGr; k++){
      AtomGr_Scale[k] = 1.0;
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
                     ((double)iter - (double)TempPara[i-1][1])
                    /((double)TempPara[i][1] - (double)TempPara[i-1][1]);

          for (k=1; k<=num_AtGr; k++){

            AtomGr_Scale[k] = sqrt((GivenTemp + (Temp_AtGr[k]-GivenTemp)*RatScale[i])/Temp_AtGr[k]);

          }

        }
      }
    }

    /* do scaling */

    if (MD_Opt_OK!=1 && iter!=MD_IterNumber){

      for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
        Gc_AN = M2G[Mc_AN];
        k = AtomGr[Gc_AN];
        for (j=1; j<=3; j++){

          if (atom_Fixed_XYZ[Gc_AN][j]==0){

            Gxyz[Gc_AN][23+j] = Gxyz[Gc_AN][23+j]*AtomGr_Scale[k];
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




void NVT_NH(int iter)
{
  /***********************************************************
   a constant temperature molecular dynamics by a Nose-Hoover
   method with velocity-Verlet integrator
  ***********************************************************/
  /*********************************************************** 
   1 a.u.=2.4189*10^-2 fs, 1fs=41.341105 a.u. 
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
  Wscale = unified_atomic_mass_unit/electron_mass;

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
      iterout_md(iter+MD_Current_Iter,MD_TimeStep*(iter+MD_Current_Iter-1),fileE);
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
      iterout_md(iter+MD_Current_Iter,MD_TimeStep*(iter+MD_Current_Iter-1),fileE);
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
      iterout_md(iter+MD_Current_Iter,MD_TimeStep*(iter+MD_Current_Iter-1),fileE);
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


static void Correct_Force()
{
  int Mc_AN,Gc_AN;
  double my_sumx,my_sumy,my_sumz;
  double sumx,sumy,sumz;

  my_sumx = 0.0;
  my_sumy = 0.0;
  my_sumz = 0.0;

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];

    if (atom_Fixed_XYZ[Gc_AN][1]==0) my_sumx += Gxyz[Gc_AN][17];  
    if (atom_Fixed_XYZ[Gc_AN][2]==0) my_sumy += Gxyz[Gc_AN][18];  
    if (atom_Fixed_XYZ[Gc_AN][3]==0) my_sumz += Gxyz[Gc_AN][19];
  }  

  MPI_Allreduce(&my_sumx, &sumx, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
  MPI_Allreduce(&my_sumy, &sumy, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
  MPI_Allreduce(&my_sumz, &sumz, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

  sumx /= (double)atomnum;
  sumy /= (double)atomnum;
  sumz /= (double)atomnum;

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];

    if (atom_Fixed_XYZ[Gc_AN][1]==0) Gxyz[Gc_AN][17] -= sumx;  
    if (atom_Fixed_XYZ[Gc_AN][2]==0) Gxyz[Gc_AN][18] -= sumy;  
    if (atom_Fixed_XYZ[Gc_AN][3]==0) Gxyz[Gc_AN][19] -= sumz;
  }
}





static void Delta_Factor(int iter)
{
  int i,j,Gc_AN,myid;
  char fileE[YOUSO10];
  double scaling_factor;
  static double tv0[4][4];
  static double Cell_Volume0;
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

    Cell_Volume0 = Cell_Volume;

    firsttime = 0;
  }

  /****************************************************
   write informatins to *.ene
  ****************************************************/

  if (myid==Host_ID){  
    sprintf(fileE,"%s%s.ene",filepath,filename);
    iterout_md(iter+MD_Current_Iter,MD_TimeStep*(iter+MD_Current_Iter-1),fileE);
  }

  /****************************************************
    scale tv 
  ****************************************************/

  if      (iter==1) scaling_factor = pow(0.940,1.0/3.0);
  else if (iter==2) scaling_factor = pow(0.960,1.0/3.0);
  else if (iter==3) scaling_factor = pow(0.980,1.0/3.0);
  else if (iter==4) scaling_factor = pow(1.020,1.0/3.0);
  else if (iter==5) scaling_factor = pow(1.040,1.0/3.0);
  else if (iter==6) scaling_factor = pow(1.060,1.0/3.0);

  if (scaling_factor<7){

    for (i=1; i<=3; i++){
      for (j=1; j<=3; j++){
	tv[i][j] = tv0[i][j]*scaling_factor;
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
    iterout_md(iter+MD_Current_Iter,MD_TimeStep*(iter+MD_Current_Iter-1),fileE);
  }

  /****************************************************
    scale tv 
  ****************************************************/

  for (i=1; i<=3; i++){
    if (MD_EvsLattice_flag[i-1]==1){
      for (j=1; j<=3; j++){
        tv[i][j] += tv0[i][j]*MD_EvsLattice_Step/100.0;
      }
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



void Calc_Temp_Atoms(int iter)
{
  int i,Mc_AN,Gc_AN,j,My_Nr,Nr;
  double My_Ukc,sum,Wscale,dt,v;

  /* calculation of temperature (K) of the atomic system */

  dt = 41.3411*MD_TimeStep;
  Wscale = unified_atomic_mass_unit/electron_mass;
  Ukc = 0.0;
  My_Ukc = 0.0;
  My_Nr = 0;
  Nr = 0;

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    sum = 0.0;
    for (j=1; j<=3; j++){
      if (atom_Fixed_XYZ[Gc_AN][j]==0){
	sum += Gxyz[Gc_AN][j+23]*Gxyz[Gc_AN][j+23];
        My_Nr++;
      }
    }
    My_Ukc += 0.5*Gxyz[Gc_AN][20]*Wscale*sum;
  }

  MPI_Allreduce(&My_Ukc, &Ukc, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
  MPI_Allreduce(&My_Nr, &Nr, 1, MPI_INT, MPI_SUM, mpi_comm_level1);

  Temp = Ukc/(0.5*kB*(double)Nr)*eV2Hartree;

  /* for old version */
  /*
  Temp = Ukc/(1.5*kB*(double)atomnum)*eV2Hartree;
  */

  /* calculation of given temperature (K) of the atomic system */

  for (i=1; i<=TempNum; i++) {
    if( (iter>TempPara[i-1][1]) && (iter<=TempPara[i][1]) ) {
      GivenTemp = TempPara[i-1][2] + (TempPara[i][2] - TempPara[i-1][2])*
          ((double)iter-(double)TempPara[i-1][1])/((double)TempPara[i][1]-(double)TempPara[i-1][1]);
    }
  }
}




int RestartFiles4GeoOpt(char *mode)
{
  int i,j,num;
  int success_flag;
  int numprocs,myid;
  char file_name[YOUSO10];
  FILE *fp;
  
  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
   
  /* write */
  if ( strcasecmp(mode,"write")==0) {

    sprintf(file_name,"%s%s_rst/%s.rst4gopt",filepath,filename,filename);

    if (myid==Host_ID){

      if ( (fp = fopen(file_name,"wb")) != NULL ){

	num = M_GDIIS_HISTORY + 1;
      
	for(i=0; i<num; i++) {
	  for(j=0; j<(atomnum+1); j++) {
	    fwrite(GxyzHistoryIn[i][j], sizeof(double), 4, fp);
	  }
	}

	for(i=0; i<num; i++) {
	  for(j=0; j<(atomnum+1); j++) {
	    fwrite(GxyzHistoryR[i][j], sizeof(double), 4, fp);
	  }
	}

	for(i=0; i<Extrapolated_Charge_History; i++) {
	  fwrite(His_Gxyz[i], sizeof(double), atomnum*3, fp);
	}

	/* Geometry_Opt_DIIS_EF */
	if (MD_switch==4){
	  for (i=0; i<(3*atomnum+2); i++){
	    fwrite(Hessian[i], sizeof(double), 3*atomnum+2, fp);
	  }
	}

	/* BFGS */
	if (MD_switch==5){
	  for (i=0; i<(3*atomnum+2); i++){
	    fwrite(InvHessian[i], sizeof(double), 3*atomnum+2, fp);
	  }
	}

	/* RF */
	if (MD_switch==6){
	  for (i=0; i<(3*atomnum+2); i++){
	    fwrite(Hessian[i], sizeof(double), 3*atomnum+2, fp);
	  }
	}

	success_flag = 1;
	fclose(fp);
      }
      else{
        printf("Failure of saving %s\n",file_name);
        success_flag = 0;
      }
    }

    else{
      success_flag = 1;
    }

    /* MPI_Bcast */    
    MPI_Bcast(&success_flag, 1, MPI_INT, Host_ID, mpi_comm_level1);
  }

  /* read */
  else if (strcasecmp(mode,"read")==0) {

    sprintf(file_name,"%s%s_rst/%s.rst4gopt",filepath,filename,filename);

    if ((fp = fopen(file_name,"rb")) != NULL){

      num = M_GDIIS_HISTORY + 1;
      
      for(i=0; i<num; i++) {
	for(j=0; j<(atomnum+1); j++) {
	  fread(GxyzHistoryIn[i][j], sizeof(double), 4, fp);
	}
      }

      for(i=0; i<num; i++) {
	for(j=0; j<(atomnum+1); j++) {
	  fread(GxyzHistoryR[i][j], sizeof(double), 4, fp);
	}
      }

      for(i=0; i<Extrapolated_Charge_History; i++) {
	fread(His_Gxyz[i], sizeof(double), atomnum*3, fp);
      }

      /* Geometry_Opt_DIIS_EF */
      if (MD_switch==4){
	for (i=0; i<(3*atomnum+2); i++){
	  fread(Hessian[i], sizeof(double), 3*atomnum+2, fp);
	}
      }

      /* BFGS */
      if (MD_switch==5){
	for (i=0; i<(3*atomnum+2); i++){
	  fread(InvHessian[i], sizeof(double), 3*atomnum+2, fp);
	}
      }

      /* RF */
      if (MD_switch==6){
	for (i=0; i<(3*atomnum+2); i++){
	  fread(Hessian[i], sizeof(double), 3*atomnum+2, fp);
	}
      }

      success_flag = 1;
      fclose(fp);

    }
    else{
      printf("Failure of reading %s\n",file_name);
      success_flag = 0;
    }
  }

  return success_flag;

}

