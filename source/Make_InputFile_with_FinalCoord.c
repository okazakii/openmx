/**********************************************************************
  Make_InputFile_with_FinalCoord.c:

     Make_InputFile_with_FinalCoord.c is a subrutine to make an input file
     with the final coordinate of the system.

  Log of Make_InputFile_with_FinalCoord.c:

     19/Sep./2007  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <string.h>

/*  stat section */
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
/*  end stat section */
#include "openmx_common.h"
#include "tran_variables.h"
#include "mpi.h"


void Make_InputFile_with_FinalCoord_Normal(char *file, int MD_iter);
void Make_InputFile_with_FinalCoord_NEGF(char *file, int MD_iter);



void Make_InputFile_with_FinalCoord(char *file, int MD_iter)
{

  if (Solver!=4)
    Make_InputFile_with_FinalCoord_Normal(file, MD_iter);
  else 
    Make_InputFile_with_FinalCoord_NEGF(file, MD_iter);

}




void Make_InputFile_with_FinalCoord_Normal(char *file, int MD_iter)
{
  int i,Gc_AN,c,n1,k;
  int restart_flag,fixed_flag;
  int geoopt_restart_flag;
  int velocity_flag;
  int rstfile_num;
  double c1,c2,c3;
  double vx,vy,vz;
  double tmpxyz[4];
  char st[800];
  char st1[800];
  char rm_operate[YOUSO10];
  char fname[YOUSO10];
  char fname1[YOUSO10];
  char fname2[YOUSO10];
  FILE *fp1,*fp2;
  char *tp;
  int numprocs,myid;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (myid==Host_ID){

    /* initialize */

    restart_flag = 0;
    geoopt_restart_flag = 0;
    fixed_flag = 0;
    velocity_flag = 0; 

    /* the new input file */    

    sprintf(fname1,"%s#",file);
    fp1 = fopen(fname1,"w");
    fseek(fp1,0,SEEK_END);

    /* the original input file */    

    fp2 = fopen(file,"r");

    if (fp2!=NULL){

      while (fgets(st,800,fp2)!=NULL){

        string_tolower(st,st1); 

        /* find the specification of <atoms.speciesandcoordinates */

        if (strncmp(st1,"<atoms.speciesandcoordinates",28)==0){

          fprintf(fp1,"%s",st);

          /* replace the atomic coordinates */

          for (i=1; i<=atomnum; i++){

            fgets(st,800,fp2);
            string_tolower(st,st1);

            /* serial number */
	    tp = strtok(st, " ");
	    if (tp!=NULL) fprintf(fp1,"%4s",tp);

            /* name of species */
            tp =strtok(NULL, " ");  
	    if (tp!=NULL) fprintf(fp1," %4s",tp);

            /* "Ang" */ 
            if (coordinates_unit==0){

	      /* x-coordinate */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",Gxyz[i][1]*BohrR);

	      /* y-coordinate */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",Gxyz[i][2]*BohrR);

	      /* z-coordinate */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",Gxyz[i][3]*BohrR);
            }

            /* AU */
            else if (coordinates_unit==1){

	      /* x-coordinate */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",Gxyz[i][1]);

	      /* y-coordinate */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",Gxyz[i][2]);

	      /* z-coordinate */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",Gxyz[i][3]);
            }

            /* FRAC */
            else if (coordinates_unit==2){

              /* The zero is taken as the origin of the unit cell. */

              tmpxyz[1] = Gxyz[i][1] - Grid_Origin[1];
              tmpxyz[2] = Gxyz[i][2] - Grid_Origin[2];
              tmpxyz[3] = Gxyz[i][3] - Grid_Origin[3];

   	      c1 = Dot_Product(tmpxyz,rtv[1])*0.5/PI;
              c2 = Dot_Product(tmpxyz,rtv[2])*0.5/PI;
              c3 = Dot_Product(tmpxyz,rtv[3])*0.5/PI;

	      /* a-axis */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",c1);

	      /* b-axis */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",c2);

	      /* c-axis */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",c3);
            }

	    while (tp!=NULL){
	      tp =strtok(NULL, " ");  
	      if (tp!=NULL) fprintf(fp1,"     %s",tp);
	    }

          } 

        }

        /* scf.restart */

        else if (strncmp(st1,"scf.restart",11)==0){
          fprintf(fp1,"scf.restart    on\n");
          restart_flag = 1;
	}

        /* geoopt.restart */

        else if (strncmp(st1,"geoopt.restart",11)==0){
          fprintf(fp1,"geoopt.restart    on\n");
          geoopt_restart_flag = 1;
          fprintf(fp1,"MD.Current.Iter  %2d\n",MD_iter+MD_Current_Iter);
	}

        /* scf.fixed.grid */

        else if (strncmp(st1,"scf.fixed.grid",14)==0){
          fprintf(fp1,"scf.fixed.grid   %18.14f  %18.14f  %18.14f\n",
                  Grid_Origin[1],Grid_Origin[2],Grid_Origin[3]);
          fixed_flag = 1;
        }  

        /* MD.Init.Velocity && VerletXYZ (NEV) */

        else if (strncmp(st1,"<md.init.velocity",16)==0
		 && (MD_switch==1 || MD_switch==2 || MD_switch==9 || MD_switch==11
		     || MD_switch==14 || MD_switch==15)){

          for (i=1; i<=(atomnum+1); i++){
            fgets(st,800,fp2);
	  }

          fprintf(fp1,"<MD.Init.Velocity\n");

          for (i=1; i<=atomnum; i++){

            vx = Gxyz[i][24]/(0.4571028*0.000001);
            vy = Gxyz[i][25]/(0.4571028*0.000001);
            vz = Gxyz[i][26]/(0.4571028*0.000001);
            fprintf(fp1," %4d    %20.14f  %20.14f  %20.14f\n", i, vx, vy, vz);
	  }

          fprintf(fp1,"MD.Init.Velocity>\n");

          velocity_flag = 1;

          /* MD.Current.Iter */

          fprintf(fp1,"\n\n");
          fprintf(fp1,"MD.Current.Iter  %2d\n",MD_iter+MD_Current_Iter);
        }  

        else{
          fprintf(fp1,"%s",st);
	}
      }

      fclose(fp2); 
    }

    /* add the restart flag if it was not found. */

    if (restart_flag==0){
      fprintf(fp1,"\n\nscf.restart    on\n");
    }

    /* add the restart flag for geometry optimization if it was not found. */

    if (geoopt_restart_flag==0 &&
       (MD_switch==3 || MD_switch==4 || MD_switch==5 || MD_switch==6 || MD_switch==7) ){
      fprintf(fp1,"\n\ngeoopt.restart    on\n");
      fprintf(fp1,"MD.Current.Iter  %2d\n",MD_iter+MD_Current_Iter);
    }

    /* add scf.fixed.grid if it was not found. */

    if (fixed_flag==0){
      fprintf(fp1,"\n\nscf.fixed.grid   %18.14f  %18.14f  %18.14f\n",
                   Grid_Origin[1],Grid_Origin[2],Grid_Origin[3]);
    }

    /* add velocity frag if it was not found. */

    if (velocity_flag==0 && (MD_switch==1 || MD_switch==2 || MD_switch==9 || MD_switch==11
			     || MD_switch==14 || MD_switch==15)){

      fprintf(fp1,"\n\n<MD.Init.Velocity\n");

      for (i=1; i<=atomnum; i++){
	fprintf(fp1," %4d    %20.14f  %20.14f  %20.14f\n",
		i,
		Gxyz[i][24]/(0.4571028*0.000001),
		Gxyz[i][25]/(0.4571028*0.000001),
		Gxyz[i][26]/(0.4571028*0.000001));
      }

      fprintf(fp1,"MD.Init.Velocity>\n");

      /* MD.Current.Iter */

      fprintf(fp1,"\n\n");
      fprintf(fp1,"MD.Current.Iter  %2d\n",MD_iter+MD_Current_Iter);
    }

    /* fclose */
    fclose(fp1); 

  } /* if (myid==Host_ID) */

}



void Make_InputFile_with_FinalCoord_NEGF(char *file, int MD_iter)
{
  int i,Gc_AN,c,n1,k;
  int restart_flag,fixed_flag;
  int velocity_flag;
  int rstfile_num;
  double c1,c2,c3;
  char st[800];
  char st1[800];
  char rm_operate[YOUSO10];
  char fname[YOUSO10];
  char fname1[YOUSO10];
  char fname2[YOUSO10];
  FILE *fp1,*fp2;
  char *tp;
  int numprocs,myid;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (myid==Host_ID){

    /* initialize */

    restart_flag = 0;
    fixed_flag = 0;
    velocity_flag = 0; 

    /* the new input file */    

    sprintf(fname1,"%s#",file);
    fp1 = fopen(fname1,"w");
    fseek(fp1,0,SEEK_END);

    /* the original input file */    

    fp2 = fopen(file,"r");

    if (fp2!=NULL){

      while (fgets(st,800,fp2)!=NULL){

        string_tolower(st,st1); 

        /* find the specification of <leftleadatoms.speciesandcoordinates */

        if (strncmp(st1,"<leftleadatoms.speciesandcoordinates",35)==0){

          fprintf(fp1,"%s",st);

          /* replace the atomic coordinates */

          for (i=1; i<=Latomnum; i++){

            fgets(st,800,fp2);
            string_tolower(st,st1);

            /* serial number */
	    tp = strtok(st, " ");
	    if (tp!=NULL) fprintf(fp1,"%4s",tp);

            /* name of species */
            tp =strtok(NULL, " ");  
	    if (tp!=NULL) fprintf(fp1," %4s",tp);

            /* "Ang" */ 
            if (coordinates_unit==0){

	      /* x-coordinate */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",Gxyz[i][1]*BohrR);

	      /* y-coordinate */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",Gxyz[i][2]*BohrR);

	      /* z-coordinate */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",Gxyz[i][3]*BohrR);
            }

            /* AU */
            else if (coordinates_unit==1){

	      /* x-coordinate */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",Gxyz[i][1]);

	      /* y-coordinate */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",Gxyz[i][2]);

	      /* z-coordinate */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",Gxyz[i][3]);
            }

	    while (tp!=NULL){
	      tp =strtok(NULL, " ");  
	      if (tp!=NULL) fprintf(fp1,"     %s",tp);
	    }

          } 
	}

        /* find the specification of <atoms.speciesandcoordinates */

        else if (strncmp(st1,"<atoms.speciesandcoordinates",28)==0){

          fprintf(fp1,"%s",st);

          /* replace the atomic coordinates */

          for (i=1; i<=Catomnum; i++){

            fgets(st,800,fp2);
            string_tolower(st,st1);

            /* serial number */
	    tp = strtok(st, " ");
	    if (tp!=NULL) fprintf(fp1,"%4s",tp);

            /* name of species */
            tp =strtok(NULL, " ");  
	    if (tp!=NULL) fprintf(fp1," %4s",tp);

            /* "Ang" */ 
            if (coordinates_unit==0){

	      /* x-coordinate */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",Gxyz[i+Latomnum][1]*BohrR);

	      /* y-coordinate */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",Gxyz[i+Latomnum][2]*BohrR);

	      /* z-coordinate */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",Gxyz[i+Latomnum][3]*BohrR);
            }

            /* AU */
            else if (coordinates_unit==1){

	      /* x-coordinate */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",Gxyz[i+Latomnum][1]);

	      /* y-coordinate */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",Gxyz[i+Latomnum][2]);

	      /* z-coordinate */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",Gxyz[i+Latomnum][3]);
            }

	    while (tp!=NULL){
	      tp =strtok(NULL, " ");  
	      if (tp!=NULL) fprintf(fp1,"     %s",tp);
	    }
          }
        }

        /* find the specification of <rightleadatoms.speciesandcoordinates */

        else if (strncmp(st1,"<rightleadatoms.speciesandcoordinates",36)==0){

          fprintf(fp1,"%s",st);

          /* replace the atomic coordinates */

          for (i=1; i<=Ratomnum; i++){

            fgets(st,800,fp2);
            string_tolower(st,st1);

            /* serial number */
	    tp = strtok(st, " ");
	    if (tp!=NULL) fprintf(fp1,"%4s",tp);

            /* name of species */
            tp =strtok(NULL, " ");  
	    if (tp!=NULL) fprintf(fp1," %4s",tp);

            /* "Ang" */ 
            if (coordinates_unit==0){

	      /* x-coordinate */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",Gxyz[i+Latomnum+Catomnum][1]*BohrR);

	      /* y-coordinate */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",Gxyz[i+Latomnum+Catomnum][2]*BohrR);

	      /* z-coordinate */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",Gxyz[i+Latomnum+Catomnum][3]*BohrR);
            }

            /* AU */
            else if (coordinates_unit==1){

	      /* x-coordinate */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",Gxyz[i+Latomnum+Catomnum][1]);

	      /* y-coordinate */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",Gxyz[i+Latomnum+Catomnum][2]);

	      /* z-coordinate */
	      tp =strtok(NULL, " ");  
	      fprintf(fp1,"  %18.14f",Gxyz[i+Latomnum+Catomnum][3]);
            }

	    while (tp!=NULL){
	      tp =strtok(NULL, " ");  
	      if (tp!=NULL) fprintf(fp1,"     %s",tp);
	    }

          } 

        }

        /* scf.restart */

        else if (strncmp(st1,"scf.restart",11)==0){
          fprintf(fp1,"scf.restart    on\n");
          restart_flag = 1;
	}

        /* scf.fixed.grid */

        else if (strncmp(st1,"scf.fixed.grid",14)==0){

          fprintf(fp1,"scf.fixed.grid   %18.14f  %18.14f  %18.14f\n",
                       Grid_Origin[1],Grid_Origin[2],Grid_Origin[3]);
          fixed_flag = 1;
        }  

        /* MD.Init.Velocity && VerletXYZ (NEV) */

        else if (strncmp(st1,"<md.init.velocity",16)==0
		 && (MD_switch==1 || MD_switch==2 || MD_switch==9 || MD_switch==11
		     || MD_switch==14 || MD_switch==15)){

          for (i=1; i<=(atomnum+1); i++){
            fgets(st,800,fp2);
	  }

          fprintf(fp1,"<MD.Init.Velocity\n");

          for (i=1; i<=atomnum; i++){
            fprintf(fp1," %4d    %20.14f  %20.14f  %20.14f\n",
                    i,
                    Gxyz[i][24]/(0.4571028*0.000001),
                    Gxyz[i][25]/(0.4571028*0.000001),
                    Gxyz[i][26]/(0.4571028*0.000001));
	  }

          fprintf(fp1,"MD.Init.Velocity>\n");

          velocity_flag = 1;
        }  

        else{
          fprintf(fp1,"%s",st);
	}
      }

      fclose(fp2); 
    }

    /* add the restart flag if it was not found. */

    if (restart_flag==0){
      fprintf(fp1,"\n\nscf.restart    on\n");
    }

    /* add scf.fixed.grid if it was not found. */

    if (fixed_flag==0){
      fprintf(fp1,"\n\nscf.fixed.grid   %18.14f  %18.14f  %18.14f\n",
                   Grid_Origin[1],Grid_Origin[2],Grid_Origin[3]);
    }

    /* add velocity frag if it was not found. */

    if (velocity_flag==0 && (MD_switch==1 || MD_switch==2 || MD_switch==9 || MD_switch==11)){

      fprintf(fp1,"\n\n<MD.Init.Velocity\n");

      for (i=1; i<=atomnum; i++){
	fprintf(fp1," %4d    %20.14f  %20.14f  %20.14f\n",
		i,
		Gxyz[i][24]/(0.4571028*0.000001),
		Gxyz[i][25]/(0.4571028*0.000001),
		Gxyz[i][26]/(0.4571028*0.000001));
      }

      fprintf(fp1,"MD.Init.Velocity>\n");
    }

    /* fclose */
    fclose(fp1); 

  } /* if (myid==Host_ID) */

}










