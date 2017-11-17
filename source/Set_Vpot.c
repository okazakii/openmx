/**********************************************************************
  Set_Vpot.c:

    Set_Vpot.c is a subroutine to calculate the value of local potential
    on each grid point.

  Log of Set_Vpot.c:

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

#ifdef noomp
#include "mimic_omp.h"
#else
#include <omp.h>
#endif



void Set_Vpot(int SCF_iter, int XC_P_switch, double *****CDM)
{
  /****************************************************
        XC_P_switch:
            0  \epsilon_XC (XC energy density)  
            1  \mu_XC      (XC potential)  
            2  \epsilon_XC - \mu_XC
  ****************************************************/

  int i,j,k,n,nmax,ri,Mc_AN,Gc_AN,Rn,GNc,GRc;
  int Nc,Nd,n1,n2,n3,Cwan,spin,MN,ct_AN;
  int h_AN,Gh_AN,Hwan,Rnh,size1,size2;
  int My_Max,Max_Size,top_num;
  int hNgrid1,hNgrid2,hNgrid3;
  double Gx,Gy,Gz,sum,x,y,z,dx,dy,dz,r,xc,yc,zc;
  double Cxyz[4];
  double *tmp_array;
  double *tmp_array2;
  double **AtomVNA_Grid,*AtomVNA2_Grid;
  int *Snd_Size,*Rcv_Size;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  double Stime_atom, Etime_atom;

  /* for OpenMP */
  int OMPID,Nthrds,Nprocs;

  MPI_Status stat;
  MPI_Request request;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
    allocation of array:

    int Snd_Size[numprocs]
    int Rcv_Size[numprocs]

    double AtomVNA_Grid[Matomnum+MatomnumF+1]
                       [GridN_Atom[Gc_AN]]
  ****************************************************/

  if ( SCF_iter<=2 && ProExpn_VNA==0 ){
    Snd_Size = (int*)malloc(sizeof(int)*numprocs);
    Rcv_Size = (int*)malloc(sizeof(int)*numprocs); 

    AtomVNA_Grid = (double**)malloc(sizeof(double*)*(Matomnum+MatomnumF+1)); 
    AtomVNA_Grid[0] = (double*)malloc(sizeof(double)*1); 
    for (Mc_AN=1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
      Gc_AN = F_M2G[Mc_AN];
      AtomVNA_Grid[Mc_AN] = (double*)malloc(sizeof(double)*GridN_Atom[Gc_AN]);
    }
  }

  /****************************************************
                       Vxc on grid
  ****************************************************/

  Set_XC_Grid(XC_P_switch,XC_switch,
               Density_Grid[0],Density_Grid[1],
               Density_Grid[2],Density_Grid[3],
               Vxc_Grid[0], Vxc_Grid[1],
               Vxc_Grid[2], Vxc_Grid[3] );

  /****************************************************
          The neutral atom potential on grids
  ****************************************************/

  if (SCF_iter<=2 && ProExpn_VNA==0){

    for (MN=0; MN<My_NumGrid1; MN++) VNA_Grid[MN] = 0.0;

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

      dtime(&Stime_atom);

      Gc_AN = M2G[Mc_AN];    
      Cwan = WhatSpecies[Gc_AN];

#pragma omp parallel shared(AtomVNA_Grid,GridN_Atom,atv,Gxyz,Gc_AN,Cwan,Mc_AN,GridListAtom,CellListAtom) private(OMPID,Nthrds,Nprocs,Nc,GNc,GRc,Cxyz,dx,dy,dz,r)
      {

	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	for (Nc=OMPID*GridN_Atom[Gc_AN]/Nthrds; Nc<(OMPID+1)*GridN_Atom[Gc_AN]/Nthrds; Nc++){

	  GNc = GridListAtom[Mc_AN][Nc];
	  GRc = CellListAtom[Mc_AN][Nc];

	  Get_Grid_XYZ(GNc,Cxyz);
	  dx = Cxyz[1] + atv[GRc][1] - Gxyz[Gc_AN][1];
	  dy = Cxyz[2] + atv[GRc][2] - Gxyz[Gc_AN][2];
	  dz = Cxyz[3] + atv[GRc][3] - Gxyz[Gc_AN][3];

	  r = sqrt(dx*dx + dy*dy + dz*dz);
	  AtomVNA_Grid[Mc_AN][Nc] = VNAF(Cwan,r);
	}

#pragma omp flush(AtomVNA_Grid)

      } /* #pragma omp parallel */

      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
    }

    for (Mc_AN=Matomnum+1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
      Gc_AN = F_M2G[Mc_AN];    
      for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
        AtomVNA_Grid[Mc_AN][Nc] = 0.0;
      }
    }

    /******************************************************
     MPI:
          AtomVNA_Grid
    ******************************************************/

    /* find data size for sending and receiving */

    tag = 999;
    My_Max = -10000;
    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      if (ID!=0){

        /*  sending size */
        if (F_Snd_Num[IDS]!=0){
          /* find data size */ 
          size1 = 0; 
          for (n=0; n<F_Snd_Num[IDS]; n++){
            Gc_AN = Snd_GAN[IDS][n];
            size1 += GridN_Atom[Gc_AN];
          }
          Snd_Size[IDS] = size1;
          MPI_Isend(&size1, 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
        }
        else{
          Snd_Size[IDS] = 0;
        }

        /*  receiving size */
        if (F_Rcv_Num[IDR]!=0){
          MPI_Recv(&size2, 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
          Rcv_Size[IDR] = size2;
        }
        else{
          Rcv_Size[IDR] = 0;
        }
        if (F_Snd_Num[IDS]!=0) MPI_Wait(&request,&stat);
      } 
      else{
        Snd_Size[IDS] = 0;
        Rcv_Size[IDR] = 0;
      }
      if (My_Max<Snd_Size[IDS]) My_Max = Snd_Size[IDS];
      if (My_Max<Rcv_Size[IDR]) My_Max = Rcv_Size[IDR];
    }  

    MPI_Allreduce(&My_Max, &Max_Size, 1, MPI_INT, MPI_MAX, mpi_comm_level1);
    tmp_array  = (double*)malloc(sizeof(double)*Max_Size);
    tmp_array2 = (double*)malloc(sizeof(double)*Max_Size);

    /* send and receive AtomVNA_Grid */

    tag = 999;
    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      if (ID!=0){

        /*  sending of data  */

        if (F_Snd_Num[IDS]!=0){

          /* find data size */
          size1 = Snd_Size[IDS];

          /* multidimentional array to vector array */
          k = -1; 
          for (n=0; n<F_Snd_Num[IDS]; n++){
            Mc_AN = Snd_MAN[IDS][n];
            Gc_AN = Snd_GAN[IDS][n];
            for (i=0; i<GridN_Atom[Gc_AN]; i++){
              k++;
              tmp_array[k] = AtomVNA_Grid[Mc_AN][i];
            }          
  	  } 
          /* MPI_Isend */
          MPI_Isend(&tmp_array[0], size1, MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);
        }

        /* receiving of block data */
        if (F_Rcv_Num[IDR]!=0){
  
	  /* find data size */
          size2 = Rcv_Size[IDR]; 
          for (i=0; i<size2; i++) tmp_array2[i] = 0.0;

          /* MPI_Recv */
          MPI_Recv(&tmp_array2[0], size2, MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

          k = -1;
          Mc_AN = F_TopMAN[IDR] - 1;
          for (n=0; n<F_Rcv_Num[IDR]; n++){
            Mc_AN++;
            Gc_AN = Rcv_GAN[IDR][n];
            for (i=0; i<GridN_Atom[Gc_AN]; i++){
              k++;
              AtomVNA_Grid[Mc_AN][i] = tmp_array2[k];
            }          
          }
        }

        if (F_Snd_Num[IDS]!=0) MPI_Wait(&request,&stat);

      }
    }  

    /******************************************************
                superposition of AtomVNA
    ******************************************************/

    for (Mc_AN=1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){

      dtime(&Stime_atom);

      Gc_AN = F_M2G[Mc_AN];
      Cwan = WhatSpecies[Gc_AN];
      for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
        MN = MGridListAtom[Mc_AN][Nc];
        if (0<=MN){ 
          VNA_Grid[MN] += AtomVNA_Grid[Mc_AN][Nc];
        }
      }

      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
    }

    /* freeing */
    free(tmp_array);
    free(tmp_array2);

  } /*  if (SCF_iter<=2){ */

  /**********************************************************
    neutral atom potential (AtomVNA2_Grid) in terms of FNAN2   
  **********************************************************/

  if (SCF_iter<=2 && ProExpn_VNA==0){

    AtomVNA2_Grid = (double*)malloc(sizeof(double)*FNAN2_Grid);

    /* MPI */
    tag = 999;
    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      if (ID!=0){

	/*****************************
               sending of data 
	*****************************/

	if (Num_Snd_FNAN2_Grid[IDS]!=0){

	  tmp_array = (double*)malloc(sizeof(double)*Num_Snd_FNAN2_Grid[IDS]);

	  /* vector array */
	  for (i=0; i<Num_Snd_FNAN2_Grid[IDS]; i++){
	    Gc_AN = Snd_FNAN2_At[IDS][i];
	    Mc_AN = F_G2M[Gc_AN];
	    Nc    = Snd_FNAN2_Nc[IDS][i];
	    tmp_array[i] = AtomVNA_Grid[Mc_AN][Nc];
	  }

	  /* MPI_Isend */
	  MPI_Isend(&tmp_array[0], Num_Snd_FNAN2_Grid[IDS], MPI_DOUBLE,
		    IDS, tag, mpi_comm_level1, &request);
	}

	/*****************************
              receiving of data
	*****************************/

	if (Num_Rcv_FNAN2_Grid[IDR]!=0){
	  top_num = TopMAN2_Grid[IDR];
	  /* MPI_Recv */
	  MPI_Recv(&AtomVNA2_Grid[top_num], Num_Rcv_FNAN2_Grid[IDR], MPI_DOUBLE,
		   IDR, tag, mpi_comm_level1, &stat);
	}

	if (Num_Snd_FNAN2_Grid[IDS]!=0){
	  MPI_Wait(&request,&stat);
	  free(tmp_array);
	}

      }
    }

    /* VNA_Grid += AtomVNA2_Grid */
    for (i=0; i<FNAN2_Grid; i++){
      MN = Rcv_FNAN2_MN[i];
      VNA_Grid[MN] += AtomVNA2_Grid[i];
    }

    free(AtomVNA2_Grid);
  }

  /****************************************************
                 external electric field
  ****************************************************/

  if ( SCF_iter<=2 && E_Field_switch==1 ){

    /* the center of the system */

    xc = 0.0;
    yc = 0.0;
    zc = 0.0;

    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      xc += Gxyz[ct_AN][1];
      yc += Gxyz[ct_AN][2];
      zc += Gxyz[ct_AN][3];
    }

    xc = xc/(double)atomnum;
    yc = yc/(double)atomnum;
    zc = zc/(double)atomnum;

    hNgrid1 = Ngrid1/2;
    hNgrid2 = Ngrid2/2;
    hNgrid3 = Ngrid3/2;

#pragma omp parallel shared(E_Field,My_Cell1,Num_Cells0,Ngrid2,Ngrid3,VEF_Grid,length_gtv,hNgrid1,hNgrid2,hNgrid3,xc,yc,zc) private(OMPID,Nthrds,Nprocs,nmax,n,i,j,k,ri,MN,dx,dy,dz)
    {

      OMPID = omp_get_thread_num();
      Nthrds = omp_get_num_threads();
      Nprocs = omp_get_num_procs();
      nmax = Num_Cells0*Ngrid2*Ngrid3; 

      for (n=OMPID*nmax/Nthrds; n<(OMPID+1)*nmax/Nthrds; n++){

	i = n/(Ngrid2*Ngrid3);
	j = (n-i*Ngrid2*Ngrid3)/Ngrid3;
	k = n - i*Ngrid2*Ngrid3 - j*Ngrid3; 
	ri = My_Cell1[i];

	MN = i*Ngrid2*Ngrid3 + j*Ngrid3 + k; 

        dx = (double)ri*length_gtv[1] + Grid_Origin[1] - xc;
        dy =  (double)j*length_gtv[2] + Grid_Origin[2] - xc;
        dz =  (double)k*length_gtv[3] + Grid_Origin[3] - xc;
	VEF_Grid[MN] = dx*E_Field[0] + dy*E_Field[1] + dz*E_Field[2];

	/*
	VEF_Grid[MN] = (double)(ri-hNgrid1)*E_Field[0]*length_gtv[1]
	              + (double)(j-hNgrid2)*E_Field[1]*length_gtv[2]
	              + (double)(k-hNgrid3)*E_Field[2]*length_gtv[3];
	*/
      }

#pragma omp flush(VEF_Grid)

    } /* #pragma omp parallel */
  }

  /****************************************************
                         Sum
  ****************************************************/

  /* spin non-collinear */
  if (SpinP_switch==3){  

    /******************************
     diagonal part 
     spin=0:  11
     spin=1:  22
    ******************************/

    if ( E_Field_switch==1 ){

      if (ProExpn_VNA==0){
        for (spin=0; spin<=1; spin++){
          for (MN=0; MN<My_NumGrid1; MN++){
            Vpot_Grid[spin][MN] = F_dVHart_flag*dVHart_Grid[MN]
                                + F_Vxc_flag*Vxc_Grid[spin][MN]
                                + F_VNA_flag*VNA_Grid[MN]
                                + F_VEF_flag*VEF_Grid[MN];
          }
        }
      }
      else{
        for (spin=0; spin<=1; spin++){
          for (MN=0; MN<My_NumGrid1; MN++){
            Vpot_Grid[spin][MN] = F_dVHart_flag*dVHart_Grid[MN]
                                + F_Vxc_flag*Vxc_Grid[spin][MN]
                                + F_VEF_flag*VEF_Grid[MN];
          }
        }
      }

    }

    else{  

      if (ProExpn_VNA==0){
        for (spin=0; spin<=1; spin++){
          for (MN=0; MN<My_NumGrid1; MN++){
            Vpot_Grid[spin][MN] = F_dVHart_flag*dVHart_Grid[MN] 
                                + F_Vxc_flag*Vxc_Grid[spin][MN]
                                + F_VNA_flag*VNA_Grid[MN];
          }
        }
      }
      else{

        for (spin=0; spin<=1; spin++){
          for (MN=0; MN<My_NumGrid1; MN++){
            Vpot_Grid[spin][MN] = F_dVHart_flag*dVHart_Grid[MN]
                                + F_Vxc_flag*Vxc_Grid[spin][MN];
          }
        }
      }
    }

    /******************************
     off-diagonal part 
     spin=2:  real 12
     spin=3:  imaginary 12
    ******************************/

    for (spin=2; spin<=3; spin++){
      for (MN=0; MN<My_NumGrid1; MN++){
        Vpot_Grid[spin][MN] = F_Vxc_flag*Vxc_Grid[spin][MN];
      }
    }
  }

  /* spin collinear */
  else{

    if ( E_Field_switch==1 ){

      if (ProExpn_VNA==0){
        for (spin=0; spin<=SpinP_switch; spin++){
          for (MN=0; MN<My_NumGrid1; MN++){
            Vpot_Grid[spin][MN] = F_dVHart_flag*dVHart_Grid[MN]
                                + F_Vxc_flag*Vxc_Grid[spin][MN]
                                + F_VNA_flag*VNA_Grid[MN]
                                + F_VEF_flag*VEF_Grid[MN];
          }
        }
      }
      else{
        for (spin=0; spin<=SpinP_switch; spin++){
          for (MN=0; MN<My_NumGrid1; MN++){
            Vpot_Grid[spin][MN] = F_dVHart_flag*dVHart_Grid[MN]
                                + F_Vxc_flag*Vxc_Grid[spin][MN]
                                + F_VEF_flag*VEF_Grid[MN];
          }
        }
      }

    }
    else{  
     
      if (ProExpn_VNA==0){
        for (spin=0; spin<=SpinP_switch; spin++){
          for (MN=0; MN<My_NumGrid1; MN++){
            Vpot_Grid[spin][MN] = F_dVHart_flag*dVHart_Grid[MN]
                                + F_Vxc_flag*Vxc_Grid[spin][MN]
                                + F_VNA_flag*VNA_Grid[MN];
          }
        }
      }
      else{
        for (spin=0; spin<=SpinP_switch; spin++){
          for (MN=0; MN<My_NumGrid1; MN++){
            Vpot_Grid[spin][MN] = F_dVHart_flag*dVHart_Grid[MN]
                                + F_Vxc_flag*Vxc_Grid[spin][MN];
          }
        }
      }

    }
  }

  /*
  for (spin=0; spin<=SpinP_switch; spin++){
    for (MN=0; MN<My_NumGrid1; MN++){
      printf("spin=%2d MN=%7d Vxc=%10.5f d %7.3f %7.3f %7.3f %7.3f\n",
              spin,MN,Vxc_Grid[spin][MN],
              Density_Grid[0][MN],Density_Grid[1][MN],Density_Grid[2][MN],Density_Grid[3][MN]);
    }
  }
  */

  /*
  {

  int nn1,nn0,MN1,MN2;

  n2 = Ngrid2/2;
  n3 = Ngrid3/2;

  for (n1=Start_Grid1[myid]; n1<=End_Grid1[myid]; n1++){
    nn1 = My_Cell0[n1];
    nn0 = n1 - Start_Grid1[myid]; 
    MN1 = nn1*Ngrid2*Ngrid3;
    MN2 = n2*Ngrid3;
    MN = MN1 + MN2 + n3;
    printf("%2d %18.15f\n",n1,Vpot_Grid[0][MN]);
  }

  }

  MPI_Finalize();
  exit(0);
  */

  /****************************************************
    freeing of array:

    double AtomVNA_Grid[Matomnum+MatomnumF+1]
                       [GridN_Atom[Gc_AN]]
  ****************************************************/

  if (SCF_iter<=2 && ProExpn_VNA==0){
    free(Snd_Size);
    free(Rcv_Size);

    for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
      free(AtomVNA_Grid[Mc_AN]);
    }
    free(AtomVNA_Grid);
  }

}
