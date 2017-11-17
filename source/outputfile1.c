#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "openmx_common.h"
#include "mpi.h"

void outputfile1(int f_switch, int MD_iter, int orbitalOpt_iter,
                 int Cnt_Now, int SCF_iter, char fname[YOUSO10], 
                 double ChemP_e0[2])
{
  FILE *fp;
  int myid;
  char buf[fp_bsize];          /* setvbuf */

  /* MPI */
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (myid==Host_ID){

    if (f_switch==1){
      fp = fopen(fname, "a");   
      if (fp==NULL){
        printf("error in outputfile1\n");
        exit(0); 
      }
      else{

        if (SCF_iter==1){
          fprintf(fp,"\n***********************************************************\n");
          fprintf(fp,"***********************************************************\n");
          fprintf(fp,"                  SCF history at MD=%2d                    \n",
                  MD_iter);
          fprintf(fp,"***********************************************************\n");
          fprintf(fp,"***********************************************************\n\n");
        }

        if ( Cnt_switch==0 || Cnt_Now!=1 ){
          fprintf(fp,"   SCF= %3d  NormRD= %15.12f  Uele= %15.12f\n",
                  SCF_iter,sqrt(fabs(NormRD[0])),Uele);
	}
        else if ( Cnt_switch==1 && Cnt_Now==1 ){
          fprintf(fp,"  OrbOpt=%3d  SCF= %3d  NormRD= %15.12f  Uele= %15.12f\n",
                  orbitalOpt_iter,SCF_iter,sqrt(fabs(NormRD[0])),Uele);
	}

        fclose(fp); 
      } 
    }

    else if (f_switch==2){
      fp = fopen(fname, "a");   

#ifdef xt3
      setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      if (fp==NULL){
        printf("error in outputfile1\n");
        exit(0); 
      }
      else{
        fprintf(fp,"\n*******************************************************\n"); 
        fprintf(fp,"        Total energy (Hartree) at MD =%2d        \n",MD_iter);
        fprintf(fp,"*******************************************************\n\n"); 

	fprintf(fp,"  Uele.    %20.12f\n\n",Uele);
	fprintf(fp,"  Ukin.    %20.12f\n",Ukin);
	fprintf(fp,"  UH0.     %20.12f\n",UH0);
	fprintf(fp,"  UH1.     %20.12f\n",UH1);
	fprintf(fp,"  Una.     %20.12f\n",Una);
	fprintf(fp,"  Unl.     %20.12f\n",Unl);
	fprintf(fp,"  Uxc0.    %20.12f\n",Uxc0);
	fprintf(fp,"  Uxc1.    %20.12f\n",Uxc1);
	fprintf(fp,"  Ucore.   %20.12f\n",Ucore);
	fprintf(fp,"  Uhub.    %20.12f\n",Uhub);
	fprintf(fp,"  Ucs.     %20.12f\n",Ucs);
        fprintf(fp,"  Uzs.     %20.12f\n",Uzs);
        fprintf(fp,"  Uzo.     %20.12f\n",Uzo);
        fprintf(fp,"  Uef.     %20.12f\n",Uef);
        fprintf(fp,"  UvdW     %20.12f\n",UvdW);
	fprintf(fp,"  Utot.    %20.12f\n\n",Utot);

	fprintf(fp,"  Note:\n\n");
	fprintf(fp,"  Utot = Ukin+UH0+UH1+Una+Unl+Uxc0+Uxc1+Ucore+Uhub+Ucs+Uzs+Uzo+Uef+UvdW\n\n");
        fprintf(fp,"  Uene:   band energy\n");
	fprintf(fp,"  Ukin:   kinetic energy\n");
	fprintf(fp,"  UH0:    electric part of screened Coulomb energy\n");
	fprintf(fp,"  UH1:    difference electron-electron Coulomb energy\n");
	fprintf(fp,"  Una:    neutral atom potential energy\n");
	fprintf(fp,"  Unl:    non-local potential energy\n");
	fprintf(fp,"  Uxc0:   exchange-correlation energy for alpha spin\n");
	fprintf(fp,"  Uxc1:   exchange-correlation energy for beta spin\n");
	fprintf(fp,"  Ucore:  core-core Coulomb energy\n");
	fprintf(fp,"  Uhub:   LDA+U energy\n");
	fprintf(fp,"  Ucs:    constraint energy for spin orientation\n");
        fprintf(fp,"  Uzs:    Zeeman term for spin magnetic moment\n");
        fprintf(fp,"  Uzo:    Zeeman term for orbital magnetic moment\n");
        fprintf(fp,"  Uef:    electric energy by electric field\n");
        fprintf(fp,"  UvdW:   semi-empirical vdW energy \n\n"); /* okuno */
        fprintf(fp,"  (see also PRB 72, 045121(2005) for the energy contributions)\n\n");

        fprintf(fp,"\n\n");

        if (Solver==4){ /* NEGF */
          fprintf(fp,"  Chemical potential of left lead  (Hartree) %20.12f\n",ChemP_e0[0]);
          fprintf(fp,"  Chemical potential of right lead (Hartree) %20.12f\n",ChemP_e0[1]);
	}
        else{ 
          fprintf(fp,"  Chemical potential (Hartree) %20.12f\n",ChemP);
	}

        fclose(fp); 
      } 
    }

    else if (f_switch==3){

      fp = fopen(fname, "a");   

#ifdef xt3
      setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      if (fp==NULL){
        printf("error in outputfile1\n");
        exit(0); 
      }

      else{

        if (orbitalOpt_iter==1){
          fprintf(fp,"\n***********************************************************\n");
          fprintf(fp,"***********************************************************\n");
          fprintf(fp,"         History of orbital optimization   MD=%2d          \n",
                  MD_iter);
          fprintf(fp,"*********     Gradient Norm ((Hartree/borh)^2)     ********\n");
          fprintf(fp,"              Required criterion= %15.12f                  \n",
                  orbitalOpt_criterion);
          fprintf(fp,"***********************************************************\n\n");
        }

        fprintf(fp,"   iter= %3d  Gradient Norm= %15.12f  Uele= %15.12f\n",
                orbitalOpt_iter,Oopt_NormD[1],Uele);
        fclose(fp); 
      } 
    }
  }

}

