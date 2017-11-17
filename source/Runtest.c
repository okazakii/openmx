/**********************************************************************
  Runtest.c:

     Runtest.c is a subroutine to check whether OpenMX runs normally 
     on many platforms or not by comparing the stored *.out and generated
     *.out on your machine.

  Log of Runtest.c:

     25/Oct/2004  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
/*  stat section */
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
/*  end stat section */
#include "openmx_common.h"
#include "Inputtools.h"
#include "mpi.h"
#include "tran_prototypes.h"
#include "tran_variables.h"


static int run_main(int argc, char *argv[], int numprocs0, int myid0);
int stringcomp( const void *a, const void *b);


typedef struct {
  char fn[YOUSO10];
} fname_type;

 

void Runtest(char *mode, int argc, char *argv[]) 
{
  FILE *fp,*fp0,*fp1,*fp2,*fp3;
  int Num_DatFiles,i,j,k,fp_OK;
  int Num_Atoms;
  int NGrid1_1,NGrid1_2,NGrid1_3;
  int NGrid2_1,NGrid2_2,NGrid2_3;
  double Utot1,Utot2,dU,dF;
  double Spread1,Spread2;
  double Omega1,Omega2;
  double gx,gy,gz,fx,fy,fz;
  double sum1,sum2;
  double time1,TotalTime;
  char fname[YOUSO10];
  char fname0[YOUSO10];
  char fname1[YOUSO10];
  char fname2[YOUSO10];
  char fname_dat[YOUSO10];
  char fname_dat2[YOUSO10];
  char fname_out1[YOUSO10];
  char fname_out2[YOUSO10];
  char ftmp[YOUSO10];
  fname_type *fndat;
  char operate[800];
  int numprocs,myid;

  char *dir;
  char *input_dir;
  char *output_file;
  DIR *dp;
  struct dirent *entry;

  MPI_Request request;
  MPI_Status  status;

  /* set up MPI */

  MPI_Comm_size(MPI_COMM_WORLD1, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD1, &myid);

  if (strcasecmp(mode,"S")==0){  
    input_dir = "input_example";
    output_file = "runtest.result";
  }
  else if (strcasecmp(mode,"L")==0){  
    input_dir = "large_example";
    output_file = "runtestL.result";
  }
  else if (strcasecmp(mode,"L2")==0){  
    input_dir = "large2_example";
    output_file = "runtestL2.result";
  }
  else if (strcasecmp(mode,"G")==0){  
    input_dir = "geoopt_example";
    output_file = "runtestG.result";
  }
  else if (strcasecmp(mode,"WF")==0){  
    input_dir = "wf_example";
    output_file = "runtestWF.result";
  }
  else if (strcasecmp(mode,"NEGF")==0){  
    input_dir = "negf_example";
    output_file = "runtestNEGF.result";
  }

  /* set Runtest_flag */

  Runtest_flag = 1;

  /* initialize TotalTime */

  TotalTime = 0.0;

  /* print the header */

  if (myid==Host_ID){

    printf("\n*******************************************************\n");  fflush(stdout);
    printf("*******************************************************\n");    fflush(stdout);
    printf(" Welcome to OpenMX   Ver. %s                           \n",Version_OpenMX); fflush(stdout);
    printf(" Copyright (C), 2002-2013, T.Ozaki                     \n");    fflush(stdout);
    printf(" OpenMX comes with ABSOLUTELY NO WARRANTY.             \n");    fflush(stdout);
    printf(" This is free software, and you are welcome to         \n");    fflush(stdout);
    printf(" redistribute it under the constitution of the GNU-GPL.\n");    fflush(stdout);
    printf("*******************************************************\n");    fflush(stdout);
    printf("*******************************************************\n\n\n");fflush(stdout);  

    printf("\n");
    printf(" OpenMX is now in the mode to check whether OpenMX runs normally\n"); fflush(stdout);
    printf(" on your machine or not by comparing the stored *.out and\n");        fflush(stdout);
    printf(" generated *.out \n"); fflush(stdout);
    printf("\n");fflush(stdout);

    /* set dir */

    dir = input_dir;

    /* count the number of dat files */

    if(( dp = opendir(dir) ) == NULL ){
      printf("could not find the directry '%s'\n",input_dir);
      MPI_Finalize();
      exit(0);
    }

    Num_DatFiles = 0;
    while((entry = readdir(dp)) != NULL){

      if ( strstr(entry->d_name,".dat")!=NULL ){ 
          
        Num_DatFiles++;
      }
    }
    closedir(dp);

    fndat = (fname_type*)malloc(sizeof(fname_type)*Num_DatFiles);

    /* store the name of dat files */

    if(( dp = opendir(dir) ) == NULL ){
      printf("could not find the directry '%s'\n",input_dir);
      MPI_Finalize();
      exit(0);
    }

    Num_DatFiles = 0;
    while((entry = readdir(dp)) != NULL){
 
      if ( strstr(entry->d_name,".dat")!=NULL ){ 

        sprintf(fndat[Num_DatFiles].fn,"%s/%s",input_dir,entry->d_name);  
        Num_DatFiles++;
      }
    }
    closedir(dp);

    /* sorting fndat */

    qsort(fndat, Num_DatFiles, sizeof(fname_type), stringcomp);  

    /*
    for (i=0; i<Num_DatFiles; i++){
      printf("i=%2d %s\n",i,fndat[i].fn);
    } 
    */

  } /* if (myid==Host_ID) */

  sprintf(fname2,"%s",output_file);

  if (myid==Host_ID){
    fp = fopen(fname2, "r");   
    if (fp!=NULL){
      fclose(fp); 
      sprintf(operate,"%s",fname2);
      remove(operate);
    }
  }

  if (myid==Host_ID){
    printf(" %2d dat files are found in the directory '%s'.\n\n\n",Num_DatFiles,input_dir);
  }

  MPI_Bcast(&Num_DatFiles, 1, MPI_INT, Host_ID, MPI_COMM_WORLD1);

  /***********************************************************
                      start calculations
  ***********************************************************/

  for (i=0; i<Num_DatFiles; i++){

    if (myid==Host_ID){
      sprintf(fname_dat,"%s",fndat[i].fn);
    }  

    MPI_Bcast(&fname_dat, YOUSO10, MPI_CHAR, Host_ID, MPI_COMM_WORLD1);

    /* run openmx */

    argv[1] = fname_dat;
    run_main(argc, argv, numprocs, myid); 

    /***********************************************************
          comparison between two files and save the result               
    ***********************************************************/

    if (myid==Host_ID){

      input_open(fname_dat);
      input_string("System.Name",fname_dat2,"default");
      input_close();

      /* compare two out files */

      sprintf(fname_out1,"%s.out",fname_dat2);
      sprintf(fname_out2,"%s/%s.out",input_dir,fname_dat2);

      /* generated file */

      input_open(fname_out1);
      input_double("Utot.",&Utot1,(double)0.0);

      /* for Wannier functions */

      if (strcasecmp(mode,"WF")==0){  
        input_double("Sum.of.Spreads.",&Spread1,(double)0.0);
        input_double("Total.Omega.=",&Omega1,(double)0.0);
      }

      input_int("Num.Grid1.",&NGrid1_1,(int)0);
      input_int("Num.Grid2.",&NGrid1_2,(int)0);
      input_int("Num.Grid3.",&NGrid1_3,(int)0);

      input_double("Elapsed.Time.",&time1,(double)0.0);

      TotalTime += time1;

      if (fp3=input_find("<coordinates.forces")) {
          
	fscanf(fp3,"%d",&Num_Atoms);

	sum1 = 0.0;
	for (j=1; j<=Num_Atoms; j++){  
	  fscanf(fp3,"%d %s %lf %lf %lf %lf %lf %lf",
		 &k,ftmp,&gx,&gy,&gz,&fx,&fy,&fz);
	  sum1 += fx + fy + fz;
	}

	if ( ! input_last("coordinates.forces>") ) {
	  printf("Format error for coordinates.forces\n");
	}
      }
      else {
	sum1 = 1000.0;
      }

      input_close();

      /* stored file */

      input_open(fname_out2);

      /* Utot */

      input_double("Utot.",&Utot2,(double)0.0);

      /* for Wannier functions */

      if (strcasecmp(mode,"WF")==0){  
        input_double("Sum.of.Spreads.",&Spread2,(double)0.0);
        input_double("Total.Omega.=",&Omega2,(double)0.0);
      }

      /* grids */

      input_int("Num.Grid1.",&NGrid2_1,(int)0);
      input_int("Num.Grid2.",&NGrid2_2,(int)0);
      input_int("Num.Grid3.",&NGrid2_3,(int)0);

      /* coordinates and forces */

      if (fp3=input_find("<coordinates.forces")) {
          
	fscanf(fp3,"%d",&Num_Atoms);

	sum2 = 0.0;
	for (j=1; j<=Num_Atoms; j++){  
	  fscanf(fp3,"%d %s %lf %lf %lf %lf %lf %lf",
		 &k,ftmp,&gx,&gy,&gz,&fx,&fy,&fz);
	  sum2 += fx + fy + fz;
	}

	if ( ! input_last("coordinates.forces>") ) {
	  /* format error */
	  printf("Format error for coordinates.forces\n");
	}
      }
      else {
	sum2 = 100.0;
      }

      input_close();

      dU = fabs(Utot1 - Utot2);
      dF = fabs(sum1 - sum2);

      /* write the result to a file, runtest.result */

      if ( (fp2 = fopen(fname2,"a")) != NULL ){

	if (  (NGrid1_1!=NGrid2_1)
	      || (NGrid1_2!=NGrid2_2)
	      || (NGrid1_3!=NGrid2_3) )
	  {
	    fprintf(fp2,"  Invalid comparison due to the different number of grids.\n");
	    fprintf(fp2,"  You may use a different radix for FFT.\n");
	  }

        if (strcasecmp(mode,"WF")==0){  

	  fprintf(fp2,"%4d  %-32.30s Elapsed time(s)=%8.2f  diff spread=%15.12f  diff Omega=%15.12f\n",
		  i+1,fname_dat,time1,fabs(Spread1-Spread2),fabs(Omega1-Omega2));
	}
        else{

	  fprintf(fp2,"%4d  %-32.30s Elapsed time(s)=%8.2f  diff Utot=%15.12f  diff Force=%15.12f\n",
		i+1,fname_dat,time1,dU,dF);
	}

	if (i==(Num_DatFiles-1)){
	  fprintf(fp2,"\n\nTotal elapsed time (s) %11.2f\n",TotalTime);
	}

	fclose(fp2);
      }
    }

  }

  /* tell us the end of calculation */
  if (myid==Host_ID){
    printf("\n\n\n\n");
    printf("The comparison can be found in a file '%s'.\n\n\n",output_file);
  }

  if (myid==Host_ID){
    free(fndat);
  }


  MPI_Barrier(MPI_COMM_WORLD1);
  MPI_Finalize();
  exit(0);
}



int stringcomp( const void *a, const void *b)
{
  return strcmp( (char*)a, (char*)b);
}



int run_main(int argc, char *argv[], int numprocs0, int myid0) 
{
  int MD_iter,i,j,po,ip;
  char fileE[YOUSO10] = ".ene"; 
  char fileDRC[YOUSO10] = ".md";
  char fileMemory[YOUSO10]; 
  char fileRestart[YOUSO10];
  char operate[200];
  double TStime,TEtime;

  /* for idle CPUs */
  int tag;
  int complete;
  MPI_Request request;
  MPI_Status  status;

  /* for measuring elapsed time */

  dtime(&TStime);

  /* allocation of CompTime */
  CompTime = (double**)malloc(sizeof(double*)*numprocs0); 
  for (i=0; i<numprocs0; i++){
    CompTime[i] = (double*)malloc(sizeof(double)*30); 
    for (j=0; j<30; j++) CompTime[i][j] = 0.0;
  }

  if (myid0==Host_ID){  
    printf("\n*******************************************************\n"); 
    printf("*******************************************************\n"); 
    printf(" Welcome to OpenMX   Ver. %s                           \n",Version_OpenMX); 
    printf(" Copyright (C), 2002-2009, T.Ozaki                     \n"); 
    printf(" OpenMX comes with ABSOLUTELY NO WARRANTY.             \n"); 
    printf(" This is free software, and you are welcome to         \n"); 
    printf(" redistribute it under the constitution of the GNU-GPL.\n");
    printf("*******************************************************\n"); 
    printf("*******************************************************\n\n"); 
  } 

  Init_List_YOUSO();
  remake_headfile = 0;
  ScaleSize = 1.2; 

  /****************************************************
                   Read the input file
  ****************************************************/

  init_alloc_first();
  CompTime[myid0][1] = readfile(argv);
  MPI_Barrier(MPI_COMM_WORLD1);

  /* initialize PrintMemory routine */

  sprintf(fileMemory,"%s%s.memory%i",filepath,filename,myid0);
  PrintMemory(fileMemory,0,"init"); 
  PrintMemory_Fix();
 
  /* initialize */
  
  init();
  fnjoint(filepath,filename,fileE);
  fnjoint(filepath,filename,fileDRC);

  /****************************************************
      SCF-DFT calculations and MD and geometrical
      optimization.
  ****************************************************/

  MD_iter = 1;

  do {

    CompTime[myid0][2] += truncation(MD_iter,1);

    if (ML_flag==1 && myid0==Host_ID) Get_VSZ(MD_iter);  

    if (Solver==4) {
      TRAN_Calc_GridBound( mpi_comm_level1, atomnum, WhatSpecies, Spe_Atom_Cut1,
                           Ngrid1, Grid_Origin, Gxyz, tv, gtv, rgtv, Left_tv, Right_tv );

      /* output: TRAN_region[], TRAN_grid_bound */
    }

    CompTime[myid0][3] += DFT(MD_iter,(MD_iter-1)%orbitalOpt_per_MDIter+1);
    if (myid0==Host_ID) iterout(MD_iter,MD_TimeStep*MD_iter,fileE,fileDRC);

    if (ML_flag==0) CompTime[myid0][4] += MD_pac(MD_iter,argv[1]);

    MD_iter++;

  } while(MD_Opt_OK==0 && MD_iter<=MD_IterNumber);

  if ( TRAN_output_hks ) {
    /* left is dummy */
    TRAN_RestartFile(mpi_comm_level1, "write","left",filepath,TRAN_hksoutfilename);
  }

  /****************************************************
               calculate Voronoi charge
  ****************************************************/
 
  if (Voronoi_Charge_flag==1) Voronoi_Charge();

  /****************************************************
  making of a file *.frac for the fractional coordinates
  ****************************************************/

  Make_FracCoord(argv[1]);

  /****************************************************
   generate Wannier functions added by Hongming Weng
  ****************************************************/

  /* hmweng */
  if(Wannier_Func_Calc){
    if (myid0==Host_ID) printf("Calling Generate_Wannier...\n");fflush(0);

    Generate_Wannier(argv[1]);
  }

  /****************************************************
                  Making of output files
  ****************************************************/

  CompTime[myid0][20] = OutData(argv[1]);

  /****************************************************
    write connectivity, Hamiltonian, overlap, density
    matrices, and etc. to a file, filename.scfout 
  ****************************************************/

  if (HS_fileout==1) SCF2File("write",argv[1]);

  /* elapsed time */

  dtime(&TEtime);
  CompTime[myid0][0] = TEtime - TStime;
  Output_CompTime();
  for (i=0; i<numprocs0; i++){
    free(CompTime[i]);
  }
  free(CompTime);

  /* merge log files */

  Merge_LogFile(argv[1]);

  /* free arrays */

  Free_Arrays(0);

  /* print memory */

  PrintMemory("total",0,"sum");

  /****************************************************
         reconstruct the original MPI group
  ****************************************************/

  {
    int *new_ranks; 
    MPI_Group  new_group,old_group; 

    new_ranks = (int*)malloc(sizeof(int)*numprocs0);
    for (i=0; i<numprocs0; i++) {
      new_ranks[i]=i; /* a new group is made of original rank=0:Pnum[k]-1 */
    }

    MPI_Comm_group(MPI_COMM_WORLD1, &old_group);

    /* define a new group */
    MPI_Group_incl(old_group,numprocs0,new_ranks,&new_group);
    MPI_Comm_create(MPI_COMM_WORLD1,new_group,&mpi_comm_level1);

    MPI_Group_free(&new_group);
    free(new_ranks); /* never forget cleaning! */
  }

  MPI_Barrier(MPI_COMM_WORLD1);
  if (myid0==Host_ID){
    printf("\nThe calculation was normally finished.\n");
  }

  return 0;
}





