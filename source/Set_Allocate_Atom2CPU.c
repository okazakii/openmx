/**********************************************************************
  Set_Allocate_Atom2CPU.c:

    Set_Allocate_Atom2CPU.c is a subroutine to allocate atoms to processors
    for the MPI parallel computation.

  Log of Set_Allocate_Atom2CPU.c:

     22/Nov/2001  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "openmx_common.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

static void Output_Atom2CPU();
static void Estimate_NL();

static int largest_bin_more_than_1datum(double *data,int n,int *bin,int numproc,int *odata); 
static int split_bin(double *data,int n,int *bin,int ilarge,int *odata); 
static double sum_in_bin_odata(double *data,int n,int istart,int iend, int *odata); 
static void Conventional_Allocation(int MD_iter, int isw, int NL_switch);
static void GDC_Allocation(int MD_iter, int isw, int NL_switch);





int Set_Allocate_Atom2CPU(int MD_iter, int isw, int NL_switch)
{
  double time0;
  time_t TStime,TEtime;

  time(&TStime);

  if (atomnum<=MYID_MPI_COMM_WORLD){
    Matomnum = 0; 
    MSpeciesNum = 0;
    return 0;
  }

  if (isw==1){
    Conventional_Allocation(MD_iter, isw, NL_switch);
  }

  /* for the generalized divide and conquer method */
  else if (NL_switch==2 && Solver==6) { 
    GDC_Allocation(MD_iter, isw, NL_switch);
  }

  else { 
    Conventional_Allocation(MD_iter, isw, NL_switch);
  }

  time(&TEtime);
  time0 = TEtime - TStime;
  return time0;
}



void GDC_Allocation(int MD_iter, int isw, int NL_switch)
{
  int time0,Num,ID,ID1,np;
  int i,j,k,num,num1,num2,ct_AN,Mct_AN;
  int Gh_AN,h_AN,Ghj_AN,hj_AN,wanA,wanB,tno1;
  int n,m,m2,n3,l1,l2,l3,numtop;
  int Rhj,Rh,R1,R2,h_AN0,h_AN_End;
  int pos[3],po,po0,po1,po_end;
  int numprocs,myid,Gc_AN,Mc_AN;
  int Matomnum2,ia,find_AN,Gc_AN2;
  int *MatomN,*Matom_GDCN;
  int Nc,*NAinCore,*CAinCore;
  int Mc_AN_GDC,Av_NAinCore,Optimum_NCluster;
  int Real_NCluster;
  double RAv_NAinCore,SD_RAv_NAinCore,Av_NAinCore0,Av_NA;
  double rcut,rcutA,rcutB;
  double sumw,d,tmp0,av_num_atoms;
  double av_weight,r2,dx,dy,dz;
  double eachw,diff,diff0,Natom;
  double Cxyz[4],tmp[4],CellV,WMatomnum;
  double *WMatomN,*bound,*weight,*num_atom;
  double a_max,a_min,da,av,ao,sum,coef1,coef3,alpha;
  double dis_top10[20],rcut_buffer,rcut_core;
  double *Cell_Gxyz2,*Gatom_num,*div_weight;
  int *found_flag,*found_flag2,*div_atom;
  int FNAN_top10[20],GAN_top10[20];
  double *Jun_FNAN,*Gatom_num2;
  int N_odata;
  int *odata, *bin;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************
            allocation of arrays
  *****************************************/

  Cell_Gxyz2 = (double*)malloc(sizeof(double)*(atomnum+1)); 
  Gatom_num  = (double*)malloc(sizeof(double)*(atomnum+1)); 
  Gatom_num2 = (double*)malloc(sizeof(double)*(atomnum+1)); 
  Jun_FNAN   = (double*)malloc(sizeof(double)*(atomnum+1)); 
  found_flag = (int*)malloc(sizeof(int)*(atomnum+1)); 
  found_flag2= (int*)malloc(sizeof(int)*(atomnum+1)); 
  weight     = (double*)malloc(sizeof(double)*(atomnum+1));
  div_atom   = (int*)malloc(sizeof(int)*(numprocs+3));
  div_weight = (double*)malloc(sizeof(double)*(numprocs+3));

  /****************************************
          calculate cell coordinates
  *****************************************/

  Cross_Product(tv[2],tv[3],tmp);
  CellV = Dot_Product(tv[1],tmp); 
  Cell_Volume = fabs(CellV); 
        
  Cross_Product(tv[2],tv[3],tmp);
  rtv[1][1] = 2.0*PI*tmp[1]/CellV;
  rtv[1][2] = 2.0*PI*tmp[2]/CellV;
  rtv[1][3] = 2.0*PI*tmp[3]/CellV;

  Cross_Product(tv[3],tv[1],tmp);
  rtv[2][1] = 2.0*PI*tmp[1]/CellV;
  rtv[2][2] = 2.0*PI*tmp[2]/CellV;
  rtv[2][3] = 2.0*PI*tmp[3]/CellV;

  Cross_Product(tv[1],tv[2],tmp);
  rtv[3][1] = 2.0*PI*tmp[1]/CellV;
  rtv[3][2] = 2.0*PI*tmp[2]/CellV;
  rtv[3][3] = 2.0*PI*tmp[3]/CellV;

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    Cxyz[1] = Gxyz[ct_AN][1];
    Cxyz[2] = Gxyz[ct_AN][2];
    Cxyz[3] = Gxyz[ct_AN][3];

    Cell_Gxyz[ct_AN][1] = Dot_Product(Cxyz,rtv[1])*0.5/PI;
    Cell_Gxyz[ct_AN][2] = Dot_Product(Cxyz,rtv[2])*0.5/PI;
    Cell_Gxyz[ct_AN][3] = Dot_Product(Cxyz,rtv[3])*0.5/PI;
  }
    
  a_max = -1.0e+10;
  a_min =  1.0e+10;
    
  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    if (a_max<Cell_Gxyz[ct_AN][1]) a_max = Cell_Gxyz[ct_AN][1];
    if (a_min>Cell_Gxyz[ct_AN][1]) a_min = Cell_Gxyz[ct_AN][1];
  }

  ao = a_min - 0.1;
  av = (a_max + 0.1) - (a_min - 0.1);

  /****************************************
          sorting of Cell_Gxyz[][1]
  *****************************************/

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    Cell_Gxyz2[ct_AN] = Cell_Gxyz[ct_AN][1];
    Gatom_num[ct_AN]  = (double)ct_AN;
    found_flag[ct_AN] = -1;
  }

  /* Sorting */
  qsort_double((long)atomnum,Cell_Gxyz2,Gatom_num); 

  /****************************************
           clustering of atoms              
  *****************************************/

  /***************************
   0) using FNAN + SNAN
  **************************/

  /*******************************************************
   find average core size and optimum number of clusters
  *******************************************************/

  tmp0 = 0.0;
  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    wanA = WhatSpecies[ct_AN];
    rcutA = Spe_Atom_Cut1[wanA];
    tmp0 += rcutA; 
  }  
  tmp0 /= (double)atomnum;

  if (tmp0<BCR)  rcut_buffer = BCR;
  else           rcut_buffer = tmp0;

  rcut_core = 0.5*rcut_buffer;

  n = 0;
  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    n++;
    for (h_AN=1; h_AN<=(FNAN[ct_AN]+True_SNAN[ct_AN]); h_AN++){
      if ( Dis[ct_AN][h_AN]<rcut_core )  n++;
    }
  }

  Av_NAinCore = n/atomnum;
  Av_NAinCore0 = (double)n/(double)atomnum;
  Optimum_NCluster = (int)((double)atomnum/Av_NAinCore0 + 0.5); 

  /*******************************************************
   reordering atoms
  *******************************************************/

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    Jun_FNAN[ct_AN] = (double)FNAN[ct_AN];
    Gatom_num2[ct_AN]  = (double)ct_AN;
  }
  qsort_double((long)atomnum,Jun_FNAN,Gatom_num2); 

  /*******************************************************
   the main loop for constructing flagments
  *******************************************************/

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){

    Gc_AN = Gatom_num2[ct_AN];
  
    do{ 

      wanA = WhatSpecies[Gc_AN];
      rcutA = Spe_Atom_Cut1[wanA];

      if (found_flag[Gc_AN]==-1){

        /********************************************
                constructing of core region
        *********************************************/

	found_flag[Gc_AN] = 0;

	for (h_AN=1; h_AN<=(FNAN[Gc_AN]+True_SNAN[Gc_AN]); h_AN++){
	  Gh_AN = natn[Gc_AN][h_AN];
	  wanB = WhatSpecies[Gh_AN];
	  rcutB = Spe_Atom_Cut1[wanB];
	  rcut = rcutA + rcutB;

	  if ( (Dis[Gc_AN][h_AN]<rcut_core) && found_flag[Gh_AN]==-1){
	    found_flag[Gh_AN]  = Gc_AN;
	    found_flag2[Gh_AN] = ncn[Gc_AN][h_AN];
	  } 
	}

	/*******************************************************
         find top 10 atoms with respect to distance from Gc_AN
        *******************************************************/

        numtop = 10;

        for (m=0; m<numtop; m++){
	  dis_top10[m]  = 0.0;
          FNAN_top10[m] = 0;
          GAN_top10[m]  = 0;
	}

        for (i=1; i<=atomnum; i++){

          if (found_flag[i]==Gc_AN){

            R1 = found_flag2[i];

            if (True_SNAN[i]==0){
              h_AN0    = 1;
              h_AN_End = FNAN[i];
            } 
            else{
              h_AN0    = FNAN[i] + 1;
              h_AN_End = FNAN[i] + True_SNAN[i];
            } 

 	    for (h_AN=h_AN0; h_AN<=h_AN_End; h_AN++){

   	      Gh_AN = natn[i][h_AN];
              R2 = ncn[i][h_AN];
 	      wanB = WhatSpecies[Gh_AN];

              if (found_flag[Gh_AN]==-1){ 

                po0 = 0;
                m = 0;
                do{
		  if (Gh_AN==GAN_top10[m]) po0 = 1;
                   m++;
                } while (po0==0 && m<numtop);

                if (po0==0){

		  l1 = atv_ijk[R2][1] - atv_ijk[R1][1];
		  l2 = atv_ijk[R2][2] - atv_ijk[R1][2];
		  l3 = atv_ijk[R2][3] - atv_ijk[R1][3];

		  dx = Gxyz[Gh_AN][1] + (double)l1*tv[1][1]
		                      + (double)l2*tv[2][1]
 	 	                      + (double)l3*tv[3][1] - Gxyz[Gc_AN][1];

		  dy = Gxyz[Gh_AN][2] + (double)l1*tv[1][2]
		                      + (double)l2*tv[2][2]
	                              + (double)l3*tv[3][2] - Gxyz[Gc_AN][2];

		  dz = Gxyz[Gh_AN][3] + (double)l1*tv[1][3]
		                      + (double)l2*tv[2][3]
		                      + (double)l3*tv[3][3] - Gxyz[Gc_AN][3];

		  r2 = dx*dx + dy*dy + dz*dz;

		  if (dis_top10[numtop-1]<r2){

		    po = 0;
		    m = 0;
		    do{
		      if (dis_top10[m]<r2)  po = 1;                    
		      else                  m++;
		    } while (po==0);

		    for (n=numtop-1; m<n; n--){
		      dis_top10[n]  =  dis_top10[n-1];
		      FNAN_top10[n] = FNAN_top10[n-1];
		      GAN_top10[n]  =  GAN_top10[n-1];
		    }

		    dis_top10[m]  = r2; 
		    FNAN_top10[m] = FNAN[Gh_AN];
		    GAN_top10[m]  = Gh_AN;

		  }
		}
	      }
	    }
          }
        }

	/********************************************************
         find atom having the largest FNAN among the top 10 atoms
        ********************************************************/

        po = 0;
        Max_FNAN = 0;
        
        for (n=0; n<numtop; n++){
          if (Max_FNAN<FNAN_top10[n]){
            Max_FNAN = FNAN[Gh_AN];
            Ghj_AN = Gh_AN; 
            po = 1;
          }                      
	}

        if (po==1) Gc_AN = Ghj_AN;

      }
      else {
        po = 0;
      }

    } while (po==1);
  }

  /*******************************************************
   reconstrucing the clusters so that the optimum number 
   of clusters can be acieved as much as possible
  *******************************************************/

  Nc = 0;

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    if (found_flag[Gc_AN]==-1){
      printf("error(2) in Set_Allocate_Atom2CPU\n"); 
      MPI_Finalize();
      exit(0); 
    }
    if (found_flag[Gc_AN]==0) Nc++;
  }

  CAinCore = (int*)malloc(sizeof(int)*(Nc+1));
  NAinCore = (int*)malloc(sizeof(int)*(Nc+1));

  Nc = 0;

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){

    if (found_flag[Gc_AN]==0){

      Nc++; 

      n = 1;
      for (i=1; i<=atomnum; i++){
        if (found_flag[i]==Gc_AN) n++;
      }

      CAinCore[Nc] = Gc_AN;
      NAinCore[Nc] = n;
    }
  }

  qsort_int((long)Nc,NAinCore,CAinCore); 

  for (n=1; n<=Nc/2; n++){
    i   = CAinCore[n];
    j   = NAinCore[n];

    CAinCore[n] = CAinCore[Nc-n+1];
    NAinCore[n] = NAinCore[Nc-n+1];
    CAinCore[Nc-n+1] = i;
    NAinCore[Nc-n+1] = j;
  }

  n = Nc;
  po_end = 0; 
  n3 = 0;

  do {

    if (Optimum_NCluster<n){

      Gc_AN = CAinCore[n];

      for (i=1; i<=atomnum; i++){

	if (found_flag[i]==Gc_AN || (i==Gc_AN && found_flag[i]==0) ){

	  tmp0 = 1.0e+100;
	  po = 0;

          /* find a neighbouring cluster */
          
	  for (j=1; j<=(FNAN[i]+True_SNAN[i]); j++){

	    k = natn[i][j];

	    if ( found_flag[k]!=Gc_AN && k!=Gc_AN && Dis[i][j]<tmp0){

	      m = k; 
	      if (found_flag[m]!=0)  m = found_flag[m];
             
	      po0 = 0;
	      l1 = 1;
	      do {
		if (CAinCore[l1]==m){
		  po0 = 1;
		}
		else l1++;
	      } while(po0==0);

	      if (NAinCore[l1]<Av_NAinCore){
                m2 = m; 
		tmp0 = Dis[i][j];
		po = 1;
	      }
	    }
	  }

          /* subtracting and adding */

	  if (po==1 && NAinCore[n]!=0){

	    found_flag[i] = m2;
	    NAinCore[n]--;

	    po0 = 0;
	    l1 = 1;
	    do {
	      if (CAinCore[l1]==m2){
		po0 = 1;
	      }
	      else l1++;
	    }while(po0==0);

	    NAinCore[l1]++; 
	  }
	}
      }

      np = n;

      n = 0;
      for (m=1; m<=Nc; m++)  if (NAinCore[m]!=0) n++;

      if (n==np){
        Av_NAinCore++;
        n3++;
      }
      else{
        n3 = 0;
      }  

      if (n3==4)  po_end = 1; 

      qsort_int((long)Nc,NAinCore,CAinCore); 

      for (m=1; m<=Nc/2; m++){
	i   = CAinCore[m];
	j   = NAinCore[m];

	CAinCore[m] = CAinCore[Nc-m+1];
	NAinCore[m] = NAinCore[Nc-m+1];
	CAinCore[Nc-m+1] = i;
	NAinCore[Nc-m+1] = j;
      }
    }
    else{
      po_end = 1; 
    }

  } while(po_end==0);


  Real_NCluster = n;

  RAv_NAinCore = 0.0;
  SD_RAv_NAinCore = 0.0;
  for (n=1; n<=Real_NCluster; n++){
    RAv_NAinCore += NAinCore[n];
    tmp0 = (double)NAinCore[n] - Av_NAinCore0;
    SD_RAv_NAinCore += tmp0*tmp0;
  }    
  RAv_NAinCore /= (double)Real_NCluster;
  SD_RAv_NAinCore = sqrt(SD_RAv_NAinCore/(double)Real_NCluster);

  /***************************
   1) using FNAN + SNAN
  **************************/

  GDC_opt = 0.6;

  /*

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    Jun_FNAN[ct_AN] = (double)FNAN[ct_AN];
    Gatom_num2[ct_AN]  = (double)ct_AN;
  }
  qsort_double((long)atomnum,Jun_FNAN,Gatom_num2); 

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){

    Gc_AN = Gatom_num2[ct_AN];

    do{ 

      wanA = WhatSpecies[Gc_AN];
      rcutA = Spe_Atom_Cut1[wanA];

      if (found_flag[Gc_AN]==-1){

	found_flag[Gc_AN] = 0;

	for (h_AN=1; h_AN<=(FNAN[Gc_AN]+True_SNAN[Gc_AN]); h_AN++){
	  Gh_AN = natn[Gc_AN][h_AN];
	  wanB = WhatSpecies[Gh_AN];
	  rcutB = Spe_Atom_Cut1[wanB];
	  rcut = rcutA + rcutB;

	  if ( (Dis[Gc_AN][h_AN]<(GDC_opt*rcut)) && found_flag[Gh_AN]==-1){
	    found_flag[Gh_AN] = Gc_AN;
	  } 
	}

        po = 0;
        Max_FNAN = -1;

        for (i=1; i<=atomnum; i++){
          if (found_flag[i]==Gc_AN){

 	    for (h_AN=(FNAN[i]+1); h_AN<=True_SNAN[i]; h_AN++){
   	      Gh_AN = natn[i][h_AN];
 	      wanB = WhatSpecies[Gh_AN];

              if (Max_FNAN<FNAN[Gh_AN] && found_flag[Gh_AN]==-1){ 
                Max_FNAN = FNAN[Gh_AN];
                Ghj_AN = Gh_AN; 
                po = 1;
              }
	    }

          }
        }

        if (po==1) Gc_AN = Ghj_AN;

      }
      else {
        po = 0;
      }

    } while (po==1);
  }
  */

  /***************************
   2) simple version
  ***************************/

  /*
  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){

    Gc_AN = Gatom_num[ct_AN];

    wanA = WhatSpecies[Gc_AN];
    rcutA = Spe_Atom_Cut1[wanA];
    
    if (found_flag[Gc_AN]==-1){

      found_flag[Gc_AN] = 0;

      for (h_AN=1; h_AN<=(FNAN[Gc_AN]+True_SNAN[Gc_AN]); h_AN++){
	Gh_AN = natn[Gc_AN][h_AN];
	wanB = WhatSpecies[Gh_AN];
	rcutB = Spe_Atom_Cut1[wanB];
	rcut = rcutA + rcutB;

	if ( (Dis[Gc_AN][h_AN]<(GDC_opt*rcut)) && found_flag[Gh_AN]==-1){
          found_flag[Gh_AN] = Gc_AN;
	} 
      }
    }
  }
  */

  /***************************
   3) simple version + alpha
  ***************************/

  /*
  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){

    Gc_AN = Gatom_num[ct_AN];

    do{ 

      wanA = WhatSpecies[Gc_AN];
      rcutA = Spe_Atom_Cut1[wanA];
    
      if (found_flag[Gc_AN]==-1){

        found_flag[Gc_AN] = 0;

        for (h_AN=1; h_AN<=(FNAN[Gc_AN]+True_SNAN[Gc_AN]); h_AN++){
	  Gh_AN = natn[Gc_AN][h_AN];
	  wanB = WhatSpecies[Gh_AN];
	  rcutB = Spe_Atom_Cut1[wanB];
	  rcut = rcutA + rcutB;

	  if ( (Dis[Gc_AN][h_AN]<(GDC_opt*rcut)) && found_flag[Gh_AN]==-1){
            found_flag[Gh_AN] = Gc_AN;
	  } 
        }
      

        po = 0;
        Max_FNAN = -1;

        for (i=1; i<=atomnum; i++){
          if (found_flag[i]==Gc_AN){
 	    for (h_AN=(FNAN[Gc_AN]+1); h_AN<=True_SNAN[Gc_AN]; h_AN++){
   	      Gh_AN = natn[Gc_AN][h_AN];
 	      wanB = WhatSpecies[Gh_AN];

              if (Max_FNAN<FNAN[Gh_AN] && found_flag[Gh_AN]==-1){ 
                Max_FNAN = FNAN[Gh_AN];
                Ghj_AN = Gh_AN; 
                po = 1;
              }

	    }
          }
        }

        if (po==1) Gc_AN = Ghj_AN;

      }
      else {
        po = 0;
      }

    } while (po==1);
  }
  */

  /* error check */
  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    if (found_flag[ct_AN]==-1){
      printf("error in clustering\n");
      MPI_Finalize();
      exit(0);
    }
  }

  if (myid==Host_ID && 2<=level_stdout){
    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      printf(" ct_AN=%2d found_flag=%4d\n",ct_AN,found_flag[ct_AN]);
    }
  }

  /* count atomnum_GDC */
  atomnum_GDC = 0;
  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    if (found_flag[ct_AN]==0) atomnum_GDC++;
  }

  if (atomnum_GDC<numprocs){
    printf("Failure of division in the GDC method\n"); 
    MPI_Finalize();
    exit(0);  
  }

  /****************************************
    freeing of arrays:
  *****************************************/
  
  if (alloc_first[19]==0){

    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){
      free(natn_GDC[ct_AN]);
    }
    free(natn_GDC);

    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){
      free(ncn_GDC[ct_AN]);
    }
    free(ncn_GDC);

  }  

  /****************************************
    find SNAN_GDC
  *****************************************/

  SNAN_GDC[0] = 0;

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){

    SNAN_GDC[ct_AN] = 0;

    if (found_flag[ct_AN]==0){

      for (j=1; j<=atomnum; j++){
        if (found_flag[j]==ct_AN){

	  for (hj_AN=0; hj_AN<=(FNAN[j]+True_SNAN[j]); hj_AN++){

            Ghj_AN = natn[j][hj_AN];
            Rhj    =  ncn[j][hj_AN]; 

	    po = 0;
	    h_AN = 0;

	    do{

	      Gh_AN = natn[ct_AN][h_AN];
              Rh    =  ncn[ct_AN][h_AN];

	      if (Ghj_AN==Gh_AN && Rhj==Rh) po = 1;

	      h_AN++;

	    } while(po==0 && h_AN<=(FNAN[ct_AN]+True_SNAN[ct_AN]));

	    if (po==0){
	      SNAN_GDC[ct_AN]++;
	    }

          }

        }
      }
    }

  } /* ct_AN */

  /****************************************
    allocation of arrays:
  *****************************************/
  
  natn_GDC = (int**)malloc(sizeof(int*)*(atomnum+1));
  for (ct_AN=0; ct_AN<=atomnum; ct_AN++){
    natn_GDC[ct_AN] = (int*)malloc(sizeof(int)*(SNAN_GDC[ct_AN]+2));
  }

  ncn_GDC = (int**)malloc(sizeof(int*)*(atomnum+1));
  for (ct_AN=0; ct_AN<=atomnum; ct_AN++){
    ncn_GDC[ct_AN] = (int*)malloc(sizeof(int)*(SNAN_GDC[ct_AN]+2));
  }

  alloc_first[19] = 0;

  /****************************************
      refind SNAN_GDC
      find natn_GDC and ncn_GDC
  *****************************************/

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){

    SNAN_GDC[ct_AN] = 0;

    if (found_flag[ct_AN]==0){

      for (j=1; j<=atomnum; j++){
        if (found_flag[j]==ct_AN){

	  for (hj_AN=0; hj_AN<=(FNAN[j]+True_SNAN[j]); hj_AN++){

            Ghj_AN = natn[j][hj_AN];
            Rhj    =  ncn[j][hj_AN]; 

            po = 0;
            h_AN = 0;

            do{

              Gh_AN = natn[ct_AN][h_AN];
              Rh    =  ncn[ct_AN][h_AN];

	      if (Ghj_AN==Gh_AN && Rhj==Rh) po = 1;

              h_AN++;

  	    } while(po==0 && h_AN<=(FNAN[ct_AN]+True_SNAN[ct_AN]));

            if (po==0){

              po0 = 0; 
              for (k=0; k<SNAN_GDC[ct_AN]; k++){
                if (natn_GDC[ct_AN][k]==Ghj_AN && ncn_GDC[ct_AN][k]==Rhj) po0 = 1;
              }

              if (po0==0){
                natn_GDC[ct_AN][SNAN_GDC[ct_AN]] = Ghj_AN;
                ncn_GDC[ct_AN][SNAN_GDC[ct_AN]]  = Rhj;
                SNAN_GDC[ct_AN]++;
	      }

            }
	  }

        } 
      }
    }
  } /* ct_AN */

  /****************************************
       calculate weight and av_weight
  *****************************************/

  sumw = 0.0;

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){

    weight[ct_AN] = 0.0;

    if (found_flag[ct_AN]==0){

      tmp0 = 0.0; 

      for (h_AN=0; h_AN<=(FNAN[ct_AN]+True_SNAN[ct_AN]+SNAN_GDC[ct_AN]); h_AN++){

        k = FNAN[ct_AN] + True_SNAN[ct_AN];

	if (h_AN<=(FNAN[ct_AN]+True_SNAN[ct_AN])) Gh_AN =     natn[ct_AN][h_AN];
        else                                      Gh_AN = natn_GDC[ct_AN][h_AN-k-1];

        wanA = WhatSpecies[Gh_AN];
        tno1 = Spe_Total_CNO[wanA];  

        tmp0 += (double)tno1;
      }

      weight[ct_AN] = tmp0*tmp0*tmp0;
    }

    sumw += weight[ct_AN]; 
  }

  av_weight = sumw/(double)numprocs;

  /****************************************
    division of system   
  *****************************************/

  div_atom[0] = 0; 
  div_weight[0] = 0.0; 

  for (ID=0; ID<numprocs; ID++){

    diff = 1.0e+200;
    sumw = 0.0;

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      Gc_AN = (int)Gatom_num[ct_AN];
    
      sumw += weight[Gc_AN];
    
      if (found_flag[Gc_AN]==0){

        diff0 = fabs((double)(ID+1)*av_weight-sumw);

        for (ID1=0; ID1<ID; ID1++){
          if (div_atom[ID1+1]==ct_AN) diff0 += 1.0e+201;
	}

        if ( diff0<diff ){
          diff = diff0;
          find_AN = ct_AN;
          div_weight[ID+1] = sumw;
        }
      }
    }

    div_atom[ID+1] = find_AN;
  }
  div_atom[numprocs] = atomnum;

  for (ID=1; ID<=numprocs; ID++){
    if ( (div_atom[ID]-div_atom[ID-1])==0 ) {
      printf("Failure of division in the GDC method\n"); 
      MPI_Finalize();
      exit(0);  
    }
  }

  /****************************************
       set Matomnum, M2G, and G2ID  
  *****************************************/

  Matomnum = 0;
  for (ct_AN=(div_atom[myid]+1); ct_AN<=div_atom[myid+1]; ct_AN++){

    Gc_AN = Gatom_num[ct_AN];

    if (found_flag[Gc_AN]==0){

      Matomnum++;

      for (j=1; j<=atomnum; j++){
        if (found_flag[j]==Gc_AN){
          Matomnum++;
	}
      }
    }
  }

  if (alloc_first[10]==0) free(M2G);
  M2G = (int*)malloc(sizeof(int)*(Matomnum+2));
  alloc_first[10] = 0;

  Matomnum = 0;
  for (ct_AN=(div_atom[myid]+1); ct_AN<=div_atom[myid+1]; ct_AN++){
    Gc_AN = Gatom_num[ct_AN];
    if (found_flag[Gc_AN]==0){

      Matomnum++;
      M2G[Matomnum] = Gc_AN;

      for (j=1; j<=atomnum; j++){
        if (found_flag[j]==Gc_AN){
          Matomnum++;
          M2G[Matomnum] = j;
	}
      }
    }
  }

  /* set G2ID */
  for (ID=0; ID<numprocs; ID++){
    for (ct_AN=(div_atom[ID]+1); ct_AN<=div_atom[ID+1]; ct_AN++){
      Gc_AN = Gatom_num[ct_AN];
      if (found_flag[Gc_AN]==0){

        G2ID[Gc_AN] = ID;

        for (j=1; j<=atomnum; j++){
          if (found_flag[j]==Gc_AN){
            G2ID[j] = ID;
	  }
        }
      }
    }
  }


  /****************************************
      set Matomnum_GDC and Mnatn_GDC
  *****************************************/

  if (alloc_first[20]==0){

    free(MFNAN_GDC);
    for (Mc_AN=0; Mc_AN<=Matomnum_GDC; Mc_AN++){
      free(Mnatn_GDC[Mc_AN]);
    }
    free(Mnatn_GDC);
  }

  Matomnum_GDC = 0;
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];   
    if (found_flag[Gc_AN]==0) Matomnum_GDC++;
  }

  MFNAN_GDC = (int*)malloc(sizeof(int)*(Matomnum_GDC+1));
  Mnatn_GDC = (int**)malloc(sizeof(int*)*(Matomnum_GDC+1));
  Mnatn_GDC[0] = (int*)malloc(sizeof(int)*1);

  Matomnum_GDC = 0;
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];   
    if (found_flag[Gc_AN]==0){

      Matomnum_GDC++;

      k = 0;
      for (j=1; j<=Matomnum; j++){
        Ghj_AN = M2G[j];   
        if (found_flag[Ghj_AN]==Gc_AN){
          k++;
	} 
      }

      MFNAN_GDC[Matomnum_GDC] = k;   
      Mnatn_GDC[Matomnum_GDC] = (int*)malloc(sizeof(int)*(k+1));
      Mnatn_GDC[Matomnum_GDC][0] = Mc_AN;
      
      k = 0;
      for (j=1; j<=Matomnum; j++){
        Ghj_AN = M2G[j];   
        if (found_flag[Ghj_AN]==Gc_AN){
          k++;
          Mnatn_GDC[Matomnum_GDC][k] = j;
	} 
      }
    }
  }

  alloc_first[20] = 0;

  /****************************************
                  stdout
  *****************************************/

  MatomN = (int*)malloc(sizeof(int)*numprocs);
  Matom_GDCN = (int*)malloc(sizeof(int)*numprocs);
  WMatomN = (double*)malloc(sizeof(double)*numprocs);

  MatomN[myid] = Matomnum;
  Matom_GDCN[myid] = Matomnum_GDC;
  WMatomN[myid] = div_weight[myid+1] - div_weight[myid];

  for (ID=0; ID<numprocs; ID++){
    MPI_Bcast(&MatomN[ID],      1, MPI_INT,    ID, mpi_comm_level1);
    MPI_Bcast(&Matom_GDCN[ID],  1, MPI_INT,    ID, mpi_comm_level1);
    MPI_Bcast(&WMatomN[ID],     1, MPI_DOUBLE, ID, mpi_comm_level1);
  } 

  if (myid==Host_ID && 1<=MD_iter){

    Av_NA = 0.0;
    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      if (found_flag[ct_AN]==0){
        Av_NA += 1.0 + (double)FNAN[ct_AN] + (double)True_SNAN[ct_AN] + (double)SNAN_GDC[ct_AN];
      }
    }
    Av_NA /= (double)Real_NCluster;

    printf("\n\n\n");fflush(stdout);

    printf("\n");
    printf("*********************************************************************\n"); 
    printf(" Clustering in generalized divide-conquer method at MD_iter=%5d \n", MD_iter );
    printf("*********************************************************************\n\n"); 
    printf("   Optimum number of clusters              %3d\n",Optimum_NCluster );fflush(stdout);
    printf("      Real number of clusters              %3d\n\n",Real_NCluster  );fflush(stdout);

    printf("   Optimum number of core atoms in each cluster %5.2f\n",     Av_NAinCore0);fflush(stdout);
    printf("   Average number of core atoms in each cluster %5.2f\n",     RAv_NAinCore);fflush(stdout);
    printf("   Standard deviation                           %5.2f\n\n",SD_RAv_NAinCore);fflush(stdout);

    printf("   Average number of atoms in each cluster      %5.2f\n\n",          Av_NA);fflush(stdout);

    printf("   The details\n\n");fflush(stdout);

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      if (found_flag[ct_AN]==0){

        n = 1;
        po = 0;
        do {
          if (ct_AN==CAinCore[n]) po = 1; 
          else                    n++;   
        } while (po==0);

        printf("   Center Atom %5d   NAinCore %4d  FNAN %4d  SNAN %4d\n",
                  ct_AN,NAinCore[n],FNAN[ct_AN],True_SNAN[ct_AN]+SNAN_GDC[ct_AN]); 
      }
    }

    if (0<level_stdout){
      printf("\n");
      printf("*******************************************************\n"); 
      printf("  Allocation of atoms to proccesors at MD_iter=%5d     \n", MD_iter );
      printf("*******************************************************\n\n"); 
    }
  }

  if (myid==Host_ID && 1<=MD_iter && 0<level_stdout){
    printf(" Number of divided clusters =%4d\n\n",atomnum_GDC);
  }

  for (ID=0; ID<numprocs; ID++){

    if (myid==Host_ID && 1<=MD_iter && 0<level_stdout){
      printf(" proc = %3d  # of clusters and atoms = %4d %4d  estimated weight=%16.5f\n",
              ID,Matom_GDCN[ID],MatomN[ID],WMatomN[ID]*1.0e-7);
    }

    if (MatomN[ID]==0){

      if (myid==Host_ID){
	printf("Failure of Atoms to CPUs\n");
	printf("Reduce the number of CPUs\n");
      }

      MPI_Finalize();
      exit(1);
    }
  }     

  if (myid==Host_ID && 1<=MD_iter && 0<level_stdout) printf("\n\n\n");

  /*
  if (myid==Host_ID){
    n = 0; 
    for (Mc_AN_GDC=1; Mc_AN_GDC<=Matomnum_GDC; Mc_AN_GDC++){
      Mc_AN = Mnatn_GDC[Mc_AN_GDC][0];
      Gc_AN = M2G[Mc_AN];
      printf("%3d Gc_AN = %3d   # of atoms in core = %3d\n",
              Mc_AN_GDC,Gc_AN,MFNAN_GDC[Mc_AN_GDC]+1);

      n = n + MFNAN_GDC[Mc_AN_GDC] + 1;
    }
    printf("Total Number of Atoms = %4d\n",n);
    printf("\n\n");
  }
  */

  /****************************************
            freeing of arrays
  *****************************************/

  free(Cell_Gxyz2); 
  free(Gatom_num);
  free(found_flag);
  free(found_flag2);
  free(weight);
  free(div_atom);
  free(div_weight);
  free(MatomN);
  free(Matom_GDCN);
  free(WMatomN);
  free(Gatom_num2);
  free(Jun_FNAN);

  free(CAinCore);
  free(NAinCore);

}




void Conventional_Allocation(int MD_iter, int isw, int NL_switch)
{
  int time0,Num,ID;
  int i,j,k,num,num1,num2,ct_AN,Mct_AN;
  int pos[3],po,po0,po1;
  int numprocs,myid,Gc_AN,Mc_AN;
  int Matomnum2,ia;
  int *MatomN,*pos_div;
  double sumw,d,tmp0,av_num_atoms;
  double eachw,diff,diff0,dx,Natom;
  double Cxyz[4],tmp[4],CellV,WMatomnum;
  double *WMatomN,*bound,*weight,*num_atom;
  double a_max,a_min,da,av,ao,sum,coef1,coef3,alpha;
  double rndx,rndy,rndz; 
  int N_odata;
  int *odata, *bin;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
     partition of the system by the one-dimentional
                 domain decomposition
  ****************************************************/

  if (isw == 0){

    if (2<=level_stdout){
      if ( 1<=NL_switch ){
        WMatomN = (double*)malloc(sizeof(double)*numprocs);
	WMatomN[myid] = 0.0;
	for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
	  Gc_AN = M2G[Mc_AN];
	  WMatomN[myid] += time_per_atom[Gc_AN];
	}
	for (ID=0; ID<numprocs; ID++){
	  MPI_Bcast(&WMatomN[ID], 1, MPI_DOUBLE, ID, mpi_comm_level1);
	} 

	if (myid==Host_ID && 0<level_stdout){
          for (ID=0; ID<numprocs; ID++){
            printf(" proc = %3d elapsed time = %16.5f\n",ID,WMatomN[ID]);
	  }
	}
        free(WMatomN);
      }
    }

    Cross_Product(tv[2],tv[3],tmp);
    CellV = Dot_Product(tv[1],tmp); 
    Cell_Volume = fabs(CellV); 
        
    Cross_Product(tv[2],tv[3],tmp);
    rtv[1][1] = 2.0*PI*tmp[1]/CellV;
    rtv[1][2] = 2.0*PI*tmp[2]/CellV;
    rtv[1][3] = 2.0*PI*tmp[3]/CellV;

    Cross_Product(tv[3],tv[1],tmp);
    rtv[2][1] = 2.0*PI*tmp[1]/CellV;
    rtv[2][2] = 2.0*PI*tmp[2]/CellV;
    rtv[2][3] = 2.0*PI*tmp[3]/CellV;

    Cross_Product(tv[1],tv[2],tmp);
    rtv[3][1] = 2.0*PI*tmp[1]/CellV;
    rtv[3][2] = 2.0*PI*tmp[2]/CellV;
    rtv[3][3] = 2.0*PI*tmp[3]/CellV;

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      Cxyz[1] = Gxyz[ct_AN][1];
      Cxyz[2] = Gxyz[ct_AN][2];
      Cxyz[3] = Gxyz[ct_AN][3];

      rndx = rnd(0.05); 
      rndy = rnd(0.05); 
      rndz = rnd(0.05); 

      MPI_Bcast(&rndx, 1, MPI_DOUBLE, Host_ID, mpi_comm_level1);
      MPI_Bcast(&rndy, 1, MPI_DOUBLE, Host_ID, mpi_comm_level1);
      MPI_Bcast(&rndz, 1, MPI_DOUBLE, Host_ID, mpi_comm_level1);

      Cell_Gxyz[ct_AN][1] = Dot_Product(Cxyz,rtv[1])*0.5/PI + rndx;
      Cell_Gxyz[ct_AN][2] = Dot_Product(Cxyz,rtv[2])*0.5/PI + rndy;
      Cell_Gxyz[ct_AN][3] = Dot_Product(Cxyz,rtv[3])*0.5/PI + rndz;
    }
    
    a_max = -1.0e+10;
    a_min =  1.0e+10;
    
    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      if (a_max<Cell_Gxyz[ct_AN][1]) a_max = Cell_Gxyz[ct_AN][1];
      if (a_min>Cell_Gxyz[ct_AN][1]) a_min = Cell_Gxyz[ct_AN][1];
    }

    ao = a_min - 0.1;
    av = (a_max + 0.1) - (a_min - 0.1);

    /****************************************
       estimate the nearest neighbor list
    *****************************************/

    if (NL_switch==1)  Estimate_NL();

    /****************************************
              allocation of arrays:
    *****************************************/

    num = 300;   /* the number of a division factor */
    bound = (double*)malloc(sizeof(double)*(num*numprocs+1));
    pos_div = (int*)malloc(sizeof(int)*(numprocs+1));
    num_atom = (double*)malloc(sizeof(double)*(num*numprocs+1));
    weight = (double*)malloc(sizeof(double)*(atomnum+1));

    /****************************************
           count the number of atoms
             with a weight factor
    *****************************************/

    alpha = 0.0001;
    coef1 = 0.001;
    coef3 = alpha*coef1;

    /* the number of atoms */
    if ( NL_switch==0 ){

      for (ct_AN=1; ct_AN<=atomnum; ct_AN++) weight[ct_AN] = 1.0;
      av_num_atoms = (double)atomnum/(double)numprocs;
    }

    /* the number of atoms weighted by (# of orbitals) */
    else if ( NL_switch==1 && (Solver==2 || Solver==3 || Solver==7 || Solver==9) ){

      sum = 0.0;
      for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
        weight[ct_AN] = coef1*(double)FNAN[ct_AN];
        sum += weight[ct_AN];
      }
      av_num_atoms = sum/(double)numprocs;

    }

    /* the number of atoms weighted by (# of orbitals)^3 */
    else if ( NL_switch==1 && Solver==5 ){

      sum = 0.0;
      for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
        num1 = SNAN[ct_AN]*SNAN[ct_AN]*SNAN[ct_AN];
        weight[ct_AN] = coef3*(double)num1 + coef1*(double)FNAN[ct_AN];
        sum += weight[ct_AN];
      }
      av_num_atoms = sum/(double)numprocs;

    }

    /* the number of atoms weighted by (# of orbitals)^2 */
    else if ( NL_switch==1 && Solver==8 ){

      sum = 0.0;
      for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
        weight[ct_AN] = 1.0e-4*(double)FNAN[ct_AN];
        sum += weight[ct_AN];
      }
      av_num_atoms = sum/(double)numprocs;
    }

    /* the number of atoms weighted by elapsed times at the previous MD step */
    else if ( NL_switch==2 ){

      sum = 0.0;
      for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
        weight[ct_AN] = time_per_atom[ct_AN];
        sum += weight[ct_AN];
      }
      av_num_atoms = sum/(double)numprocs;

    }

    /* the number of atoms */
    else {

      for (ct_AN=1; ct_AN<=atomnum; ct_AN++) weight[ct_AN] = 1.0;
      av_num_atoms = (double)atomnum/(double)numprocs;
    }

    /* set num_atom */

    da = fabs(av)/(double)(numprocs*num);
    bound[0] = ao;
    num_atom[0] = 0.0;
    for (i=1; i<=numprocs*num; i++){
      bound[i] = ao + da*(double)i;
      num_atom[i] = 0.0;
      for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
        if (Cell_Gxyz[ct_AN][1]<bound[i]){
          num_atom[i] += weight[ct_AN];
	}
      }
    }

    /****************************************
                divide the cell
    *****************************************/

    /* int N_odata;
       int *odata, *bin; */

    N_odata=0; 
    odata=(int*)malloc(sizeof(int)*(atomnum+1));
    bin = (int*)malloc(sizeof(int)*(numprocs+1)) ; 

    for (k=0; k<=numprocs; k++){
      bin[k]=0;
    }

    bin[0]=1;

    /*  set position of division */
    pos_div[0] = 0;
    for (k=1; k<=numprocs; k++){
      Natom = av_num_atoms*(double)k; 
      diff = 10000000.0;
      for (i=1; i<=numprocs*num; i++){
        diff0 = fabs(num_atom[i]-Natom);
        if (diff0<diff){
          diff = diff0;
          pos_div[k] = i;
        }
      }

      po0 = pos_div[k-1];
      po1 = pos_div[k];

      /* get Matomnum2 */
      Matomnum2 = 1; 
      for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
        if (   bound[po0]<=Cell_Gxyz[ct_AN][1] 
	       && Cell_Gxyz[ct_AN][1]<bound[po1] ){

          Matomnum2++; 
          odata[++N_odata] = ct_AN;

        }
      }

      bin[k] = N_odata+1; 

    } /* for (k=... */

    /* shift if some of bin[k]-bin[k-1] are zero */

    for (k=1; k<=numprocs; k++) {

      while (bin[k]<=bin[k-1] && bin[k]) {

	/* shift data after k */
	for (i=k; i<numprocs; i++) {
	  bin[i] = bin[i+1];
	} 
	bin[numprocs] = 0;  

      }
    }   

    /* if some bin[] is zero,  split other bin[] and increase the number of bin[] */

    while (bin[numprocs]==0) {

      int ilarge=-1,isplit=-1;  /* initial values are for error check */

      ilarge = largest_bin_more_than_1datum(weight, atomnum,bin,numprocs,odata);
      if (ilarge) {
	isplit = split_bin(weight, atomnum,bin, ilarge,odata);
      }

      /* error check, but this must not occur. */
      if ( ilarge==0 || isplit==0 ) {
        if (Host_ID==myid) {
          printf("error in bin[]: ");
          for (i=0; i<=numprocs; i++)  printf("%d ", bin[i]);
          printf("\nilarge=%d isplit=%d\n",ilarge, isplit);
        }
        MPI_Finalize();
        exit(10);
      }

      if (myid==Host_ID) {
	printf("correct: ilarge=%d isplit=%d\n",ilarge,isplit);
      }

      /*
	bin[ilarge+1]=bin[ilarge];  ...
	bin[ilarge]=isplit;
      */

      for (i=numprocs;i>ilarge;i--) {
	bin[i]=bin[i-1];
      }
      bin[ilarge] = isplit; 

    }

    for (k=1; k<=numprocs; k++){

      Matomnum2=bin[k]-bin[k-1]; /* bin[k-1]: bin[k]-1  is the range */
      po0 = pos_div[k-1];
      po1 = pos_div[k];

      if (myid==(k-1)){
	if (alloc_first[10]==0) free(M2G);
	M2G = (int*)malloc(sizeof(int)*(Matomnum2+1));
	alloc_first[10] = 0;
      }

      /* store M2G */

      for (i=bin[k-1];i<bin[k]; i++) {
	ct_AN = odata[i];
	G2ID[ct_AN] = k-1; 
	if (myid==(k-1)) M2G[i-bin[k-1]+1] = ct_AN;
      }

      if (myid==(k-1)){
	Matomnum = Matomnum2;
	WMatomnum = num_atom[po1] - num_atom[po0];
      }

    } /* for (k=... */

    if (myid==Host_ID && 2<=level_stdout){
      for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
	printf(" Gc_AN=%i G2ID=%i\n",ct_AN,G2ID[ct_AN]); 
      }
    }

    MatomN = (int*)malloc(sizeof(int)*numprocs);
    WMatomN = (double*)malloc(sizeof(double)*numprocs);

    MatomN[myid]  = Matomnum;
    WMatomN[myid] = WMatomnum;

    for (ID=0; ID<numprocs; ID++){
      MPI_Bcast(&MatomN[ID],  1, MPI_INT, ID, mpi_comm_level1);
      MPI_Bcast(&WMatomN[ID], 1, MPI_DOUBLE, ID, mpi_comm_level1);
    } 

    /* find Max Matomnum */

    Max_Matomnum = 0;
    for (ID=0; ID<numprocs; ID++){
      if (Max_Matomnum<MatomN[ID]) Max_Matomnum = MatomN[ID];
    }     

    /* print the information */

    if (myid==Host_ID && 1<=MD_iter && 0<level_stdout){
      printf("\n");
      printf("*******************************************************\n"); 
      printf("  Allocation of atoms to proccesors at MD_iter=%5d     \n", MD_iter );
      printf("*******************************************************\n\n"); 
    }

    for (ID=0; ID<numprocs; ID++){

      if (myid==Host_ID && 1<=MD_iter && 0<level_stdout){
	printf(" proc = %3d  # of atoms=%4d  estimated weight=%16.5f\n",
	       ID,MatomN[ID],WMatomN[ID]);
      }

      if (MatomN[ID]==0){

	if (myid==Host_ID){
	  printf("Failure of Atoms to CPUs\n");
	  printf("Reduce the number of CPUs\n");
	}

	MPI_Finalize();
	exit(1);
      }
    }     

    if (myid==Host_ID && 1<=MD_iter && 0<level_stdout) printf("\n\n\n");

    /****************************************
                freeing of arrays:
    *****************************************/

    free(bin);
    free(odata);

    free(MatomN);
    free(WMatomN);
    free(bound);
    free(pos_div);
    free(num_atom);
    free(weight);

    MPI_Barrier(mpi_comm_level1);

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++) time_per_atom[ct_AN] = 0.0;

  }

  /****************************************************
                  partition of species
  ****************************************************/

  else if (isw == 1){

    for (i=0; i<numprocs; i++){
      Species_Top[i] = 0;
      Species_End[i] = 0;
    }

    if (SpeciesNum<numprocs) Num_Procs2 = SpeciesNum;
    else                     Num_Procs2 = numprocs;

    num1 = SpeciesNum/Num_Procs2;
    num2 = SpeciesNum%Num_Procs2;
    for (i=0; i<Num_Procs2; i++){
      Species_Top[i] = num1*i;
      Species_End[i] = num1*(i + 1) - 1;
    }
    if (num2!=0){
      for (i=0; i<num2; i++){
	Species_Top[i] = Species_Top[i] + i;
	Species_End[i] = Species_End[i] + i + 1;
      }
      for (i=num2; i<Num_Procs2; i++){
	Species_Top[i] = Species_Top[i] + num2;
	Species_End[i] = Species_End[i] + num2;
      }
    }

    if (myid==Host_ID && 2<=level_stdout){
      for (i=0; i<Num_Procs2; i++){
	printf("proc=%4d  Species_Top=%4d  Species_End=%4d\n",
	       i,Species_Top[i],Species_End[i]);
      }
    }

    if (myid<Num_Procs2)
      MSpeciesNum = Species_End[myid] - Species_Top[myid] + 1;
    else 
      MSpeciesNum = 0;

    if (myid==Host_ID && 2<=level_stdout){
      printf("myid=%i  MSpeciesNum=%i\n",myid,MSpeciesNum);
    }
 
  }

  /*
  if(isw == 1){
    Set_Inf_SndRcv();
    if (myid==Host_ID) Output_Atom2CPU();
  }
  */
}




void Estimate_NL()
{
  /****************************************************
          FNAN, SNAN by the physical truncation
  ****************************************************/

  int i,j,k,ct_AN,wanA,wanB,Rn;
  int N,TN_Rn,tno1,tno2;
  int numprocs,myid,tag=999,ID;
  double rcutA,rcutB,rcut,rcut2;
  double r,r2,dx,dy,dz,BCR2;
  double **cell_atv;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
                    set cell vectors
  ****************************************************/

  if (CellNN_flag==0)  N = 2;
  else                 N = 1;

  TN_Rn = (2*N+1)*(2*N+1)*(2*N+1);
 
  /* allocation of cell_atv */
  cell_atv = (double**)malloc(sizeof(double*)*TN_Rn);
  for (Rn=0; Rn<TN_Rn; Rn++){
    cell_atv[Rn] = (double*)malloc(sizeof(double)*4);
  }

  Rn = 1;
 
  for (i=-N; i<=N; i++){
    for (j=-N; j<=N; j++){
      for (k=-N; k<=N; k++){
	if (i==0 && j==0 && k==0){
	  cell_atv[0][1] = 0.0;
	  cell_atv[0][2] = 0.0;
	  cell_atv[0][3] = 0.0;
	}
	else{
	  cell_atv[Rn][1] = (double)i*tv[1][1] + (double)j*tv[2][1] + (double)k*tv[3][1];
	  cell_atv[Rn][2] = (double)i*tv[1][2] + (double)j*tv[2][2] + (double)k*tv[3][2];
	  cell_atv[Rn][3] = (double)i*tv[1][3] + (double)j*tv[2][3] + (double)k*tv[3][3];
	  Rn++;
	}
      }
    }
  }

  /****************************************************
             find FNAN, SNAN, tno1, and tno2 
  ****************************************************/

  BCR2 = BCR*BCR;

#pragma omp parallel shared(level_stdout,NOHS_C,cell_atv,Gxyz,TN_Rn,atomnum,SNAN,FNAN,Spe_Total_CNO,Spe_Atom_Cut1,WhatSpecies,M2G,BCR2,Matomnum) private(i,ct_AN,wanA,rcutA,tno1,j,wanB,rcutB,rcut,rcut2,tno2,Rn,dx,dy,dz,r2)
  {

    for (i=1; i<=Matomnum; i++){
        
      ct_AN = M2G[i];
      wanA = WhatSpecies[ct_AN];
      rcutA = Spe_Atom_Cut1[wanA];
      tno1 = Spe_Total_CNO[wanA];

      FNAN[ct_AN] = 0;
      SNAN[ct_AN] = 0;

      for (j=1; j<=atomnum; j++){

	wanB = WhatSpecies[j];
	rcutB = Spe_Atom_Cut1[wanB];
	rcut = rcutA + rcutB;
	rcut2 = rcut*rcut;
	tno2 = Spe_Total_CNO[wanB];

	for (Rn=0; Rn<TN_Rn; Rn++){

	  if ((ct_AN==j) && Rn==0){

	    FNAN[ct_AN] += tno2;  
	    SNAN[ct_AN] += tno2;
	  }
 
	  else{

	    if (((Rn==0)&&(ct_AN!=j)) || Rn!=0){
	      dx = fabs(Gxyz[ct_AN][1] - Gxyz[j][1] - cell_atv[Rn][1]);
	      dy = fabs(Gxyz[ct_AN][2] - Gxyz[j][2] - cell_atv[Rn][2]);
	      dz = fabs(Gxyz[ct_AN][3] - Gxyz[j][3] - cell_atv[Rn][3]);
	      r2 = dx*dx + dy*dy + dz*dz;
	    }

	    if (0.01<r2 && r2<=rcut2)    FNAN[ct_AN] += tno2;
	    if (NOHS_C==1){
	      if (0.01<r2 && r2<=rcut2)  SNAN[ct_AN] += tno2;
	    }
	    else{
	      if (0.01<r2 && r2<=BCR2)   SNAN[ct_AN] += tno2;
	    }

	  }
	}
      }

      if (2<=level_stdout){
	printf("<Estimate_NL>  ct_AN=%2d FNAN SNAN %2d %2d\n",ct_AN,FNAN[ct_AN],SNAN[ct_AN]);
      }

    } /* i */

  } /* #pragma omp parallel */

  /* freeing of cell_atv */
  for (Rn=0; Rn<TN_Rn; Rn++){
    free(cell_atv[Rn]);
  }
  free(cell_atv);

  /* MPI, FNAN and SNAN */

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    ID = G2ID[ct_AN];
    MPI_Bcast(&FNAN[ct_AN], 1, MPI_INT, ID, mpi_comm_level1);
    MPI_Bcast(&SNAN[ct_AN], 1, MPI_INT, ID, mpi_comm_level1);


  }

}







/*
void Output_Atom2CPU()
{

  int i,numprocs;
  char file_A2C[YOUSO10] = ".A2C";
  FILE *fp;

  MPI_Comm_size(mpi_comm_level1,&numprocs);

  fnjoint(filepath,filename,file_A2C);

  if ((fp = fopen(file_A2C,"w")) != NULL){

    fprintf(fp,"\n");
    fprintf(fp,"***********************************************************\n");
    fprintf(fp,"***********************************************************\n");
    fprintf(fp,"               Allocation of atoms to CPUs                 \n");
    fprintf(fp,"***********************************************************\n");
    fprintf(fp,"***********************************************************\n");
    fprintf(fp,"\n");

    fprintf(fp,"   Average weight = %5.4f\n\n",eachw);
    fprintf(fp,"   CPU    Atoms    Top     End     CPU_Weight\n");
    for (i=0; i<numprocs; i++){
      fprintf(fp,"  %4d   %4d    %4d    %4d    %4d\n",
              i,
              Gatom_End[i]-Gatom_Top[i]+1,
              Gatom_Top[i], 
              Gatom_End[i],
              CPU_Weight[i]);
    } 
    fclose(fp);
  }
  else
    printf("Failure of saving the Set_Allocate_Atom2CPU.\n");
}
*/



int  largest_bin_more_than_1datum(double *data,int n,int *bin,int numproc,
                                  int *odata)
{
  int i;
  int ilarge; double vlarge,val;

  ilarge=0;
  vlarge=-1000000.0;

  for (i=1; i<=numproc; i++) {

    val = sum_in_bin_odata(data,n,bin[i-1],bin[i],odata);

    if (val>vlarge && bin[i]-bin[i-1]>1 ) {
      ilarge = i;
      vlarge = val;
    }
  } 

/*
  printf("largest bin=%d\n",ilarge);
*/

  return ilarge;

}


int split_bin(double *data,int n,int *bin,int ilarge,int *odata)
{
  int i,ismall;
  double val1,val2;
  double diff,smallestdiff;

  ismall=0;
  smallestdiff=1e+20;

  for (i=bin[ilarge-1]+1; i<bin[ilarge];i++) {
    val1 = sum_in_bin_odata(data,n,bin[ilarge-1],i,odata); 
    val2 = sum_in_bin_odata(data,n,i,bin[ilarge],odata); 
    if ( fabs(val1-val2) <smallestdiff ) {
      smallestdiff = fabs(val1-val2);
      ismall = i;
    } 
  } 

/*
  printf("isplit=%d\n",ismall);
*/

  return ismall;
}



double  sum_in_bin_odata(double *data,int n,int istart,int iend, int *odata)
{
  double sum;
  int i;

  for (sum=0.0, i=istart; i<iend;i++) {
    sum += data[odata[i]];
  }

  /*
  printf("[%d:%d] sum=%lf\n",istart,iend-1,sum);
  */

  return sum;
}

