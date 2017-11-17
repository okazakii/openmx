/**********************************************************************
  iterout.c:

     iterout.c is a subroutine to output xyz-coordinates
     at each MD step to filename.md and filename.md2.

  Log of iterout.c:

     22/Nov/2001  Released by T.Ozaki
     14/May/2004  Modified by M.Ohfuti
     22/Nov/2011  Modified by T.Ohwaki

***********************************************************************/

#include <stdio.h>
#include "openmx_common.h"

void iterout(int iter,double drctime,char fileSE[YOUSO10],char fileSDRC[YOUSO10])
{
  int i,j,k,kk,myid;
  double dt,itermax,aa,dx,dy,dz,xido,yido,zido;
  char fileXYZ[YOUSO10];
  FILE *fp;
  char buf[fp_bsize];          /* setvbuf */

  if(MD_switch==1 || MD_switch==2 || MD_switch==9 || MD_switch==11 || MD_switch==15){
    Calc_Temp_Atoms(iter);
  }

  if(MD_switch==14){
  }

  MPI_Comm_rank(mpi_comm_level1,&myid);
  if (myid!=Host_ID) return;

  /****************************************************
     cartesian coordinates for MD or geometry opt.
  ****************************************************/

  if ( ((iter-1) % 1)==0 ){

    if ((fp = fopen(fileSDRC,"a")) != NULL){

#ifdef xt3
      setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      fprintf(fp,"%i \n",atomnum);

      if(MD_switch==1 || MD_switch==2 || MD_switch==9 || MD_switch==11 || MD_switch==15) {
        fprintf(fp,"time= %8.3f (fs) Energy= %8.5f (Hartree) Temperature= %8.3f (Given Temp.= %8.3f)\n",drctime,Utot,Temp,GivenTemp);
      }

      else if (MD_switch==14) {
        fprintf(fp,"time= %8.3f (fs) Energy= %8.5f (Hartree) ",drctime,Utot);
        for (kk=1; kk<=num_AtGr; kk++){ 
          fprintf(fp,"Temperature[%d]= %8.3f ",kk,Temp_AtGr[kk]);
        }
        fprintf(fp,"\n");
      }

      else {
        fprintf(fp,"   time= %8.3f (fs)  Energy= %8.5f (Hartree)\n",drctime,Utot);
      } 

      for (k=1; k<=atomnum; k++){
        i = WhatSpecies[k];
        j = Spe_WhatAtom[i];

        if ( MD_switch==1 || 
             MD_switch==2 ||  
             MD_switch==9 ||
             MD_switch==11||
             MD_switch==14||
             MD_switch==15
            ) {

          /***********************************************
             5.291772083*10^{-11} m / 2.418884*10^{-17} s 
             = 2.1876917*10^6 m/s                         
             = 1 a.u. for velocity 

             1 m/s = 0.4571028 * 10^{-6} a.u.
          ***********************************************/

          fprintf(fp,"%4s   %8.5f  %8.5f  %8.5f  %14.6f %14.6f %14.6f  %8.5f  %8.5f  %8.5f  %8.5f\n",
                  Atom_Symbol[j],                
	   	  Gxyz[k][1]*BohrR,Gxyz[k][2]*BohrR,Gxyz[k][3]*BohrR,
		  Gxyz[k][24]/(0.4571028*0.000001),
                  Gxyz[k][25]/(0.4571028*0.000001),
                  Gxyz[k][26]/(0.4571028*0.000001),
                  InitN_USpin[k],
                  InitN_DSpin[k],
                  InitN_USpin[k]+InitN_DSpin[k],
                  InitN_USpin[k]-InitN_DSpin[k]);

        }
        else {

          fprintf(fp,"%4s   %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n",
                  Atom_Symbol[j],                
                  Gxyz[k][1]*BohrR,Gxyz[k][2]*BohrR,Gxyz[k][3]*BohrR,
		  Gxyz[k][17],Gxyz[k][18],Gxyz[k][19]);

        }
      }
      fclose(fp);
    }
    else{
      printf("failure of saving md file\n");
      fclose(fp);
    }
  }

  /****************************************************
      cartesian coordinates of the final structure
  ****************************************************/

  sprintf(fileXYZ,"%s2",fileSDRC);
  if ((fp = fopen(fileXYZ,"w")) != NULL){

#ifdef xt3
    setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

    fprintf(fp,"%i \n",atomnum);

    if (MD_switch==2 || MD_switch==11 || MD_switch==14 || MD_switch==15) {
      fprintf(fp,"time= %8.3f (fs) Energy= %8.5f (Hartree) Temperature= %8.3f\n",drctime,Utot,Temp);
    }
    else {
      fprintf(fp,"   time= %8.3f (fs)  Energy= %8.5f (Hartree)\n",drctime,Utot);
    } 

    for (k=1; k<=atomnum; k++){

      i = WhatSpecies[k];

      fprintf(fp,"%6d   %4s  %12.7f  %12.7f  %12.7f   %8.5f  %8.5f\n",
                k,
                SpeName[i],
                Gxyz[k][1]*BohrR,Gxyz[k][2]*BohrR,Gxyz[k][3]*BohrR,
                0.5*Spe_Core_Charge[i],0.5*Spe_Core_Charge[i]);                      
    }
    fclose(fp);
  }
  else{
    printf("failure of saving md file\n");
    fclose(fp);
  }

}
