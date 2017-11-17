/**********************************************************************
  TRAN_Poisson.c:

     TRAN_Poisson.c is a subroutine to solve Poisson's equation with 
     the boundary conditions using FFT and a finite difference method.

  Log of TRAN_Poisson.c:

     09/July/2008  Released by T.Ozaki

***********************************************************************/

#define  measure_time   0

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "openmx_common.h"


#ifdef nompi
#include "mimic_mpi.h"
#else
#include <mpi.h>
#endif

#include "tran_prototypes.h"
#include "tran_variables.h"
#include "lapack_prototypes.h"


#ifdef fftw2 
#include <fftw.h>
#else
#include <fftw3.h> 
#endif


static void FFT2D_Poisson(int NG1, int NG2,
			  int Ng1, int Ng2,  int Ng3,
			  int My_Ng1, int My_Ng2,
			  int sgn2, int sgn3,
			  double ***ReF1, double ***ImF1,
			  double ***ReF2, double ***ImF2);

static void MPI_Communication_F1_F2(int NG1, int NG2,
				    int Ng1, int Ng2,  int Ng3,
				    int My_Ng1, int My_Ng2,
				    double ***ReF1, double ***ImF1,
				    double ***ReF2, double ***ImF2);

double TRAN_Poisson_FD(double ***ReV1, double ***ImV1,
		       double ***ReV2, double ***ImV2);

double TRAN_Poisson_FD_TCN(double ***ReV1, double ***ImV1,
		           double ***ReV2, double ***ImV2);

double TRAN_Poisson_FD_LCN(double ***ReV1, double ***ImV1,
		           double ***ReV2, double ***ImV2);

double TRAN_Poisson_FFT(double ***ReV1, double ***ImV1,
		        double ***ReV2, double ***ImV2);

double TRAN_Poisson_FFT_Extended(double ***ReV1, double ***ImV1,
				 double ***ReV2, double ***ImV2);
	


double TRAN_Poisson(double ***ReV1, double ***ImV1,
		    double ***ReV2, double ***ImV2)
{ 
  double time;

  if (atomnum<=MYID_MPI_COMM_WORLD) return 0.0;

  switch (TRAN_Poisson_flag){

    case 1:
      time = TRAN_Poisson_FD(ReV1, ImV1, ReV2, ImV2);
    break;

    case 2:
      time = TRAN_Poisson_FD_TCN(ReV1, ImV1, ReV2, ImV2);
    break;

    case 3:
      time = TRAN_Poisson_FD_LCN(ReV1, ImV1, ReV2, ImV2);
    break;

    case 4:
      time = TRAN_Poisson_FFT(ReV1, ImV1, ReV2, ImV2);
    break;

    case 5:
      time = TRAN_Poisson_FFT_Extended(ReV1, ImV1, ReV2, ImV2);
    break;
  }

  return time;
}


double TRAN_Poisson_FFT_Extended(double ***ReV1, double ***ImV1,
				 double ***ReV2, double ***ImV2)
{ 
  int i,k1,k2,k3,kk2,ip,ri,k,j;
  int MN0,MN1,MN2,MN,nnn1,n1,n2,n3,nn1;
  int Ngrid1_extended;
  double time0,G2,x0,x1;
  double tmp0,sk1,sk2,sk3,tot_den,sum;
  double x,y,z,Gx,Gy,Gz,Eden[2];
  double ReTmp,ImTmp,c_coe,p_coe;
  double da2,Gpara2,a,cc,cv,gateV;
  double TStime,TEtime;
  double v0,v1,b,v;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  double *tmp_array0;
  double *tmp_array1;
  double tv_ex[4][4];
  double rtv_ex[4][4],tmp[4];
  double CellV_ex,Cell_Volume_ex;
  dcomplex ac,bc,vc,vf0,vf1,vb0,vb1;

  fftw_complex *in, *out;
  fftw_plan pm,pp;

  MPI_Status stat;
  MPI_Request request;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (myid==Host_ID) printf("<TRAN_Poisson_FFT_Extended>  Solving Poisson's equation...\n");

  MPI_Barrier(mpi_comm_level1);
  dtime(&TStime);

  /****************************************************
          FFT of charge density on the b-c plane
  ****************************************************/
  
  FFT2D_Density(0,ReV1,ImV1,ReV2,ImV2);

  /****************************************************
          
  ****************************************************/

  /* total number of grids along the a-axis */

  Ngrid1_extended = IntNgrid1_e[0] + Ngrid1 + IntNgrid1_e[1];

  /*
  printf("ABC %2d %2d %2d\n",IntNgrid1_e[0],Ngrid1,IntNgrid1_e[1]);fflush(stdout);
  MPI_Barrier(mpi_comm_level1);
  */

  /* set tv_ex and rtv_ex */
 
  if      (1.0<fabs(tv[1][1])) ip = 1;
  else if (1.0<fabs(tv[1][2])) ip = 2;
  else if (1.0<fabs(tv[1][3])) ip = 3;

  tv_ex[1][1] = 0.0;
  tv_ex[1][2] = 0.0;
  tv_ex[1][3] = 0.0;
  tv_ex[1][ip] = (double)Ngrid1_extended*length_gtv[1];

  tv_ex[2][1] = tv[2][1];
  tv_ex[2][2] = tv[2][2];
  tv_ex[2][3] = tv[2][3];

  tv_ex[3][1] = tv[3][1];
  tv_ex[3][2] = tv[3][2];
  tv_ex[3][3] = tv[3][3];

  Cross_Product(tv_ex[2],tv_ex[3],tmp);
  CellV_ex = Dot_Product(tv_ex[1],tmp); 
  Cell_Volume_ex = fabs(CellV_ex);
  
  Cross_Product(tv_ex[2],tv_ex[3],tmp);
  rtv_ex[1][1] = 2.0*PI*tmp[1]/CellV_ex;
  rtv_ex[1][2] = 2.0*PI*tmp[2]/CellV_ex;
  rtv_ex[1][3] = 2.0*PI*tmp[3]/CellV_ex;
  
  Cross_Product(tv_ex[3],tv_ex[1],tmp);
  rtv_ex[2][1] = 2.0*PI*tmp[1]/CellV_ex;
  rtv_ex[2][2] = 2.0*PI*tmp[2]/CellV_ex;
  rtv_ex[2][3] = 2.0*PI*tmp[3]/CellV_ex;
  
  Cross_Product(tv_ex[1],tv_ex[2],tmp);
  rtv_ex[3][1] = 2.0*PI*tmp[1]/CellV_ex;
  rtv_ex[3][2] = 2.0*PI*tmp[2]/CellV_ex;
  rtv_ex[3][3] = 2.0*PI*tmp[3]/CellV_ex;

  /* allocate in, out, pm, and pp */

  in  = fftw_malloc(sizeof(fftw_complex)*(Ngrid1_extended+2)); 
  out = fftw_malloc(sizeof(fftw_complex)*(Ngrid1_extended+2)); 

  pm = fftw_plan_dft_1d(Ngrid1_extended, in,out,-1,FFTW_ESTIMATE);
  pp = fftw_plan_dft_1d(Ngrid1_extended, in,out, 1,FFTW_ESTIMATE);

  /* 1D FFT along a-axis */

  tmp0 = 4.0*PI/(double)Ngrid1_extended;

  for (k2=0; k2<My_NGrid2_Poisson; k2++){

    kk2 = k2 + Start_Grid2[myid];

    if (kk2<Ngrid2/2) sk2 = (double)kk2;
    else              sk2 = (double)(kk2 - Ngrid2);

    for (k3=0; k3<Ngrid3; k3++){

      if (k3<Ngrid3/2) sk3 = (double)k3;
      else             sk3 = (double)(k3 - Ngrid3);

      for (n1=0; n1<IntNgrid1_e[0]; n1++){
        in[n1][0] = dDen_IntBoundary[0][n1][kk2][k3].r;
        in[n1][1] = dDen_IntBoundary[0][n1][kk2][k3].i;
      }

      for (n1=0; n1<Ngrid1; n1++){
	in[n1+IntNgrid1_e[0]][0] = ReV2[k2][n1][k3];
	in[n1+IntNgrid1_e[0]][1] = ImV2[k2][n1][k3];
      }

      for (n1=0; n1<IntNgrid1_e[1]; n1++){
        in[IntNgrid1_e[0]+Ngrid1+n1][0] = dDen_IntBoundary[1][n1][kk2][k3].r;
        in[IntNgrid1_e[0]+Ngrid1+n1][1] = dDen_IntBoundary[1][n1][kk2][k3].i;
      }

      fftw_execute(pm);

      for (k1=0; k1<Ngrid1_extended; k1++){

        if (k1<Ngrid1_extended/2) sk1 = (double)k1;
        else                      sk1 = (double)(k1 - Ngrid1_extended);

        Gx = sk1*rtv_ex[1][1] + sk2*rtv_ex[2][1] + sk3*rtv_ex[3][1];
        Gy = sk1*rtv_ex[1][2] + sk2*rtv_ex[2][2] + sk3*rtv_ex[3][2]; 
        Gz = sk1*rtv_ex[1][3] + sk2*rtv_ex[2][3] + sk3*rtv_ex[3][3];

        G2 = Gx*Gx + Gy*Gy + Gz*Gz;

        if (k1==0 && kk2==0 && k3==0){
          in[k1][0] = 0.0;
          in[k1][1] = 0.0;
	}
        else{
          in[k1][0] = tmp0*out[k1][0]/G2;          
          in[k1][1] = tmp0*out[k1][1]/G2;
	}
      }

      fftw_execute(pp);

      /* get dVhart */

      for (n1=0; n1<Ngrid1; n1++){

        ReV2[k2][n1][k3] = out[n1+IntNgrid1_e[0]][0];
        ImV2[k2][n1][k3] = out[n1+IntNgrid1_e[0]][1];
      }

    }
  }

  fftw_destroy_plan(pm);  
  fftw_destroy_plan(pp);  

  free(in);
  free(out);

  /****************************************************
   MPI communication: ReV2 and ImV2 into ReV1 and ImV1
  ****************************************************/
 
  MPI_Communication_F1_F2(2, 1, Ngrid2, Ngrid1, Ngrid3,
                          My_NGrid2_Poisson, My_NGrid1_Poisson,
                          ReV2, ImV2, ReV1, ImV1);

  /****************************************************
      Inverse FFT of ReV1 and ImV1 on the b-c plane
  ****************************************************/

  FFT2D_Poisson(1, 2, Ngrid1, Ngrid2, Ngrid3,
                My_NGrid1_Poisson, My_NGrid2_Poisson,
                1, 1, ReV1, ImV1, ReV2, ImV2);

  /****************************************************
   find the difference Hartree potential in real space
  ****************************************************/

  /* initialize */

  for (MN=0; MN<My_NumGrid1; MN++){
    dVHart_Grid[MN] = 0.0;
  }

  /* use their potential using MPI */

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    /* Isend */
    if (Num_ISnd_Grid1[IDS]!=0){
      tmp_array0 = (double*)malloc(sizeof(double)*Num_ISnd_Grid1[IDS]*Ngrid2*Ngrid3);

      for (i=0; i<Num_ISnd_Grid1[IDS]; i++){ 
	n1 = ISnd_Grid1[IDS][i];
	nn1 = n1 - Start_Grid1[myid];
	MN1 = i*Ngrid2*Ngrid3;
	for (n2=0; n2<Ngrid2; n2++){
	  MN2 = n2*Ngrid3;

          for (n3=0; n3<Ngrid3; n3++){
	    MN = MN1 + MN2 + n3;
	    tmp_array0[MN] = ReV1[nn1][n2][n3];
	  }
	}
      } 

      tag = 999;
      MPI_Isend(&tmp_array0[0], Num_ISnd_Grid1[IDS]*Ngrid2*Ngrid3,
		MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);
    }

    /* Recv */
    if (Num_IRcv_Grid1[IDR]!=0){

      tmp_array1 = (double*)malloc(sizeof(double)*Num_IRcv_Grid1[IDR]*Ngrid2*Ngrid3); 

      tag = 999;
      MPI_Recv(&tmp_array1[0], Num_IRcv_Grid1[IDR]*Ngrid2*Ngrid3,
	       MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

      for (i=0; i<Num_IRcv_Grid1[IDR]; i++){ 
	n1 = IRcv_Grid1[IDR][i];
	nn1 = My_Cell0[n1];
	MN1 = nn1*Ngrid2*Ngrid3;
	MN0 = i*Ngrid2*Ngrid3;

	for (n2=0; n2<Ngrid2; n2++){
	  MN2 = n2*Ngrid3;
	  for (n3=0; n3<Ngrid3; n3++){
	    dVHart_Grid[MN1 + MN2 + n3] = tmp_array1[MN0 + MN2 + n3];
	  }
	}
      }
      free(tmp_array1);
    }

    if (Num_ISnd_Grid1[IDS]!=0){
      MPI_Wait(&request,&stat);
      free(tmp_array0);
    }
  }

  /* use own potential */
  for (n1=0; n1<My_NGrid1_Poisson; n1++){
    nn1 = Start_Grid1[myid] + n1;
    nnn1 = My_Cell0[nn1];
    if (My_Cell0[nn1]!=-1){
      MN1 = nnn1*Ngrid2*Ngrid3; 
      for (n2=0; n2<Ngrid2; n2++){
	MN2 = n2*Ngrid3;
	for (n3=0; n3<Ngrid3; n3++){
	  MN = MN1 + MN2 + n3;
	  dVHart_Grid[MN] = ReV1[n1][n2][n3];
	}
      }
    }
  }

  /****************************************************
                add the boundary condition
  ****************************************************/

  if      (1.0<fabs(tv[1][1])) ip = 1;
  else if (1.0<fabs(tv[1][2])) ip = 2;
  else if (1.0<fabs(tv[1][3])) ip = 3;

  x0 = Grid_Origin[ip];
  x1 = (double)Ngrid1*length_gtv[1] + Grid_Origin[ip];

  for (i=0; i<Num_Cells0; i++){

    ri = My_Cell1[i];
    x = (double)ri*length_gtv[1] + Grid_Origin[ip];

    for (j=0; j<Ngrid2; j++){
      for (k=0; k<Ngrid3; k++){

        v0 = dVHart_Grid_e[0][j*Ngrid3+k];
        v1 = dVHart_Grid_e[1][j*Ngrid3+k];

        a = (v1 - v0)/(x1 - x0);
        b = v0 - a*x0;

	/*
        v = a*x + b;
	*/

        v = a*x;

	MN = i*Ngrid2*Ngrid3 + j*Ngrid3 + k; 
	dVHart_Grid[MN] += v;

      }
    }
  }

  /****************************************************
               apply the gate voltage
  ****************************************************/

  if      (1.0<fabs(tv[1][1])) ip = 1;
  else if (1.0<fabs(tv[1][2])) ip = 2;
  else if (1.0<fabs(tv[1][3])) ip = 3;

  cc = 0.0;
  for (i=1; i<=Catomnum; i++){
    cc += Gxyz[Latomnum+i][ip];
  }  
  cc /= (double)Catomnum;

  a = 1.0*( fabs(tv[1][ip]) - fabs(Left_tv[1][ip]) - fabs(Right_tv[1][ip]) );

  for (i=0; i<Num_Cells0; i++){

    ri = My_Cell1[i];
    cv = (double)ri*length_gtv[1] + Grid_Origin[ip];

    gateV = tran_gate_voltage*exp( -pow( (cv-cc)/a, 8.0) );

    for (j=0; j<Ngrid2; j++){
      for (k=0; k<Ngrid3; k++){

	MN = i*Ngrid2*Ngrid3 + j*Ngrid3 + k; 
	dVHart_Grid[MN] += gateV;

      }
    }
  }

  /* for time */
  MPI_Barrier(mpi_comm_level1);
  dtime(&TEtime);
  time0 = TEtime - TStime;
  return time0;
}




double TRAN_Poisson_FFT(double ***ReV1, double ***ImV1,
		        double ***ReV2, double ***ImV2)
{ 
  int k1,k2,k3,kk2;
  int nn1,MN1,MN2,ID0,ID1;
  int ri,i,j,k,ip,MN;
  double x0,x1,v0,v1,v;
  double cc,cv,a,b,gateV;
  double time0;
  double tmp0,sk1,sk2,sk3,tot_den;
  double x,y,z,Gx,Gy,Gz,Eden[2],DenA,G2;
  double ReTmp,ImTmp,c_coe,p_coe;
  double *dVH_FFT0;
  double *dVH_FFT1;
  double Av_dVH_FFT0,Av_dVH_FFT1;
  double Av_dVH0,Av_dVH1;
  double TStime,TEtime;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  
  MPI_Status stat;
  MPI_Request request;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (myid==Host_ID) printf("<TRAN_Poisson_FFT>  Solving Poisson's equation...\n");

  MPI_Barrier(mpi_comm_level1);
  dtime(&TStime);

  /****************************************************
                 FFT of charge density 
  ****************************************************/

  FFT_Density(0,ReV1,ImV1,ReV2,ImV2);

  /****************************************************
                       4*PI/G2/N^3
  ****************************************************/

  tmp0 = 4.0*PI/(double)Ngrid1/(double)Ngrid2/(double)Ngrid3;

  for (k2=0; k2<My_NGrid2_Poisson; k2++){

    kk2 = k2 + Start_Grid2[myid];

    if (kk2<Ngrid2/2) sk2 = (double)kk2;
    else              sk2 = (double)(kk2 - Ngrid2);

    for (k1=0; k1<Ngrid1; k1++){

      if (k1<Ngrid1/2) sk1 = (double)k1;
      else             sk1 = (double)(k1 - Ngrid1);

      for (k3=0; k3<Ngrid3; k3++){

        if (k3<Ngrid3/2) sk3 = (double)k3;
        else             sk3 = (double)(k3 - Ngrid3);

        Gx = sk1*rtv[1][1] + sk2*rtv[2][1] + sk3*rtv[3][1];
        Gy = sk1*rtv[1][2] + sk2*rtv[2][2] + sk3*rtv[3][2]; 
        Gz = sk1*rtv[1][3] + sk2*rtv[2][3] + sk3*rtv[3][3];
        G2 = Gx*Gx + Gy*Gy + Gz*Gz;

        if (k1==0 && kk2==0 && k3==0){
          ReV2[k2][k1][k3] = 0.0;
          ImV2[k2][k1][k3] = 0.0;
        }
        else{
          ReV2[k2][k1][k3] = tmp0*ReV2[k2][k1][k3]/G2;
          ImV2[k2][k1][k3] = tmp0*ImV2[k2][k1][k3]/G2; 
        }
      }
    }
  }

  /****************************************************
        find the Hartree potential in real space
  ****************************************************/
  
  Get_Value_inReal(0,ReV2,ImV2,ReV1,ImV1,dVHart_Grid,dVHart_Grid);

  /****************************************************
    MPI communication of the difference Hartree 
    potential at the boudaries.
  ****************************************************/

  ID0 = Cell_ID0[0]; 
  ID1 = Cell_ID0[Ngrid1-1]; 

  dVH_FFT0 = (double*)malloc(sizeof(double)*Ngrid2*Ngrid3);  
  dVH_FFT1 = (double*)malloc(sizeof(double)*Ngrid2*Ngrid3);  

  if (myid==ID0){

    nn1 = My_Cell0[0];  
    MN = nn1*Ngrid2*Ngrid3;

    for (j=0; j<Ngrid2; j++){
      for (k=0; k<Ngrid3; k++){
        MN1 = j*Ngrid3 + k;
        MN2 = MN + MN1;
        dVH_FFT0[MN1] = dVHart_Grid[MN2]; 
      }
    }    
  }

  MPI_Bcast(&dVH_FFT0[0], Ngrid2*Ngrid3, MPI_DOUBLE, ID0, mpi_comm_level1);

  if (myid==ID1){

    nn1 = My_Cell0[Ngrid1-1];  
    MN = nn1*Ngrid2*Ngrid3;

    for (j=0; j<Ngrid2; j++){
      for (k=0; k<Ngrid3; k++){
        MN1 = j*Ngrid3 + k;
        MN2 = MN + MN1;
        dVH_FFT1[MN1] = dVHart_Grid[MN2]; 
      }
    }    
  }

  MPI_Bcast(&dVH_FFT1[0], Ngrid2*Ngrid3, MPI_DOUBLE, ID1, mpi_comm_level1);

  Av_dVH_FFT0 = 0.0;
  Av_dVH_FFT1 = 0.0;
  Av_dVH0 = 0.0;
  Av_dVH1 = 0.0;

  for (j=0; j<Ngrid2; j++){
    for (k=0; k<Ngrid3; k++){

      MN1 = j*Ngrid3 + k;

      Av_dVH_FFT0 += dVH_FFT0[MN1];
      Av_dVH_FFT1 += dVH_FFT1[MN1];

      Av_dVH0 += dVHart_Grid_e[0][MN1];
      Av_dVH1 += dVHart_Grid_e[1][MN1];
    }
  } 

  Av_dVH_FFT0 = Av_dVH_FFT0/(double)(Ngrid2*Ngrid3);
  Av_dVH_FFT1 = Av_dVH_FFT1/(double)(Ngrid2*Ngrid3);
  Av_dVH0 = Av_dVH0/(double)(Ngrid2*Ngrid3);
  Av_dVH1 = Av_dVH1/(double)(Ngrid2*Ngrid3);

  /****************************************************
                add the boundary condition
  ****************************************************/

  if      (1.0e-4<fabs(tv[1][1])) ip = 1;
  else if (1.0e-4<fabs(tv[1][2])) ip = 2;
  else if (1.0e-4<fabs(tv[1][3])) ip = 3;

  x0 = Grid_Origin[ip];
  x1 = (double)Ngrid1*length_gtv[1] + Grid_Origin[ip];

  for (i=0; i<Num_Cells0; i++){

    ri = My_Cell1[i];
    x = (double)ri*length_gtv[1] + Grid_Origin[ip];

    for (j=0; j<Ngrid2; j++){
      for (k=0; k<Ngrid3; k++){

        v0 = Av_dVH0 - Av_dVH_FFT0;
        v1 = Av_dVH1 - Av_dVH_FFT1;

        a = (v1 - v0)/(x1 - x0);
        b = v0 - a*x0;
        v = a*x + b;

	MN = i*Ngrid2*Ngrid3 + j*Ngrid3 + k; 
	dVHart_Grid[MN] += v;

      }
    }
  }

  /****************************************************
                 apply the gate voltage
  ****************************************************/

  if      (1.0<fabs(tv[1][1])) ip = 1;
  else if (1.0<fabs(tv[1][2])) ip = 2;
  else if (1.0<fabs(tv[1][3])) ip = 3;

  cc = 0.0;
  for (i=1; i<=Catomnum; i++){
    cc += Gxyz[Latomnum+i][ip];
  }  
  cc /= (double)Catomnum;

  a = 1.0*( fabs(tv[1][ip]) - fabs(Left_tv[1][ip]) - fabs(Right_tv[1][ip]) );

  for (i=0; i<Num_Cells0; i++){

    ri = My_Cell1[i];
    cv = (double)ri*length_gtv[1] + Grid_Origin[ip];

    gateV = tran_gate_voltage*exp( -pow( (cv-cc)/a, 8.0) );

    for (j=0; j<Ngrid2; j++){
      for (k=0; k<Ngrid3; k++){

	MN = i*Ngrid2*Ngrid3 + j*Ngrid3 + k; 
	dVHart_Grid[MN] += gateV;

      }
    }
  }

  /*  
  if (myid==0 && 0){
  
  for (i=0; i<Num_Cells0; i++){
  
    ri = My_Cell1[i];
    j = 0;
    k = 0;
  
    MN = i*Ngrid2*Ngrid3 + j*Ngrid3 + k; 
  
    printf("%2d %2d %15.12f\n",i,ri,dVHart_Grid[MN]); fflush(stdout);
  }
  }
  */

  /* freeing of arrays */

  free(dVH_FFT0);
  free(dVH_FFT1);

  /* for time */
  MPI_Barrier(mpi_comm_level1);
  dtime(&TEtime);
  time0 = TEtime - TStime;
  return time0;
}






double TRAN_Poisson_FD(double ***ReV1, double ***ImV1,
		       double ***ReV2, double ***ImV2)
{ 
  int i,k1,k2,k3,kk2,ip,ri,k,j;
  int MN0,MN1,MN2,MN,nnn1,n1,n2,n3,nn1;
  double time0;
  double tmp0,sk1,sk2,sk3,tot_den,sum;
  double x,y,z,Gx,Gy,Gz,Eden[2];
  double ReTmp,ImTmp,c_coe,p_coe;
  double da2,Gpara2,a,cc,cv,gateV;
  double TStime,TEtime;
  dcomplex *DL,*D,*DU,*B;
  INTEGER n,nrhs,ldb,info;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  double *tmp_array0;
  double *tmp_array1;

  MPI_Status stat;
  MPI_Request request;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (myid==Host_ID) printf("<TRAN_Poisson_FD>  Solving Poisson's equation...\n");

  MPI_Barrier(mpi_comm_level1);
  dtime(&TStime);

  /****************************************************
                 allocation of arrays
  ****************************************************/

  DL = (dcomplex*)malloc(sizeof(dcomplex)*Ngrid1);
  D  = (dcomplex*)malloc(sizeof(dcomplex)*Ngrid1);
  DU = (dcomplex*)malloc(sizeof(dcomplex)*Ngrid1);
  B  = (dcomplex*)malloc(sizeof(dcomplex)*Ngrid1);

  /****************************************************
          FFT of charge density on the b-c plane
  ****************************************************/
  
  FFT2D_Density(0,ReV1,ImV1,ReV2,ImV2);

  /****************************************************
          solve finite difference equations
  ****************************************************/
        
  da2 = Dot_Product(gtv[1],gtv[1]);

  for (k2=0; k2<My_NGrid2_Poisson; k2++){

    kk2 = k2 + Start_Grid2[myid];

    if (kk2<Ngrid2/2) sk2 = (double)kk2;
    else              sk2 = (double)(kk2 - Ngrid2);

    for (k3=0; k3<Ngrid3; k3++){

      if (k3<Ngrid3/2) sk3 = (double)k3;
      else             sk3 = (double)(k3 - Ngrid3);

      Gx = sk2*rtv[2][1] + sk3*rtv[3][1];
      Gy = sk2*rtv[2][2] + sk3*rtv[3][2]; 
      Gz = sk2*rtv[2][3] + sk3*rtv[3][3];
      Gpara2 = Gx*Gx + Gy*Gy + Gz*Gz;

      for (k1=0; k1<(Ngrid1-1); k1++){
        DL[k1].r = 1.0;
        DL[k1].i = 0.0;
        DU[k1].r = 1.0;
        DU[k1].i = 0.0;
      }         

      for (k1=0; k1<(Ngrid1-0); k1++){
        D[k1].r = -2.0 - da2*Gpara2;
        D[k1].i = 0.0;
      }

      /* scale terms with Gpara=0 */ 

      if (kk2==0 && k3==0){

        for (k1=0; k1<Ngrid1; k1++){
          ReV2[k2][k1][k3] *= TRAN_Poisson_Gpara_Scaling;
        }
      }

      /* set B */

      for (k1=0; k1<(Ngrid1-0); k1++){
        B[k1].r = -4.0*PI*da2*ReV2[k2][k1][k3];
        B[k1].i = -4.0*PI*da2*ImV2[k2][k1][k3];
      }

      /* add the boundary condition */

      B[0       ].r -= VHart_Boundary[0][Ngrid1_e[0]-1][kk2][k3].r;
      B[0       ].i -= VHart_Boundary[0][Ngrid1_e[0]-1][kk2][k3].i;
      B[Ngrid1-1].r -= VHart_Boundary[1][0            ][kk2][k3].r;
      B[Ngrid1-1].i -= VHart_Boundary[1][0            ][kk2][k3].i;

      /* solve the linear equation */ 

      n = Ngrid1-0;
      nrhs = 1;
      ldb = Ngrid1-0;
     
      F77_NAME(zgtsv,ZGTSV)(&n, &nrhs, DL, D, DU, B, &ldb, &info);

      /* store B to ReV2 and ImV2 */

      for (k1=0; k1<(Ngrid1-0); k1++){
        ReV2[k2][k1][k3] = B[k1].r;
        ImV2[k2][k1][k3] = B[k1].i;
      }

      /* set the boundary for the left part */

      /*
      ReV2[k2][0][k3] = VHart_Boundary[0][kk2][k3].r;
      ImV2[k2][0][k3] = VHart_Boundary[0][kk2][k3].i;
      */

    } /* k3 */
  } /* k2 */

  /****************************************************
   MPI communication: ReV2 and ImV2 into ReV1 and ImV1
  ****************************************************/
 
  MPI_Communication_F1_F2(2, 1, Ngrid2, Ngrid1, Ngrid3,
                          My_NGrid2_Poisson, My_NGrid1_Poisson,
                          ReV2, ImV2, ReV1, ImV1);

  /****************************************************
      Inverse FFT of ReV1 and ImV1 on the b-c plane
  ****************************************************/

  FFT2D_Poisson(1, 2, Ngrid1, Ngrid2, Ngrid3,
                My_NGrid1_Poisson, My_NGrid2_Poisson,
                1, 1, ReV1, ImV1, ReV2, ImV2);

  /****************************************************
   find the difference Hartree potential in real space
  ****************************************************/

  /* initialize */

  for (MN=0; MN<My_NumGrid1; MN++){
    dVHart_Grid[MN] = 0.0;
  }

  /* use their potential using MPI */

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    /* Isend */
    if (Num_ISnd_Grid1[IDS]!=0){
      tmp_array0 = (double*)malloc(sizeof(double)*Num_ISnd_Grid1[IDS]*Ngrid2*Ngrid3);

      for (i=0; i<Num_ISnd_Grid1[IDS]; i++){ 
	n1 = ISnd_Grid1[IDS][i];
	nn1 = n1 - Start_Grid1[myid];
	MN1 = i*Ngrid2*Ngrid3;
	for (n2=0; n2<Ngrid2; n2++){
	  MN2 = n2*Ngrid3;

          for (n3=0; n3<Ngrid3; n3++){
	    MN = MN1 + MN2 + n3;
	    tmp_array0[MN] = ReV1[nn1][n2][n3];
	  }
	}
      } 

      tag = 999;
      MPI_Isend(&tmp_array0[0], Num_ISnd_Grid1[IDS]*Ngrid2*Ngrid3,
		MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);
    }

    /* Recv */
    if (Num_IRcv_Grid1[IDR]!=0){

      tmp_array1 = (double*)malloc(sizeof(double)*Num_IRcv_Grid1[IDR]*Ngrid2*Ngrid3); 

      tag = 999;
      MPI_Recv(&tmp_array1[0], Num_IRcv_Grid1[IDR]*Ngrid2*Ngrid3,
	       MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

      for (i=0; i<Num_IRcv_Grid1[IDR]; i++){ 
	n1 = IRcv_Grid1[IDR][i];
	nn1 = My_Cell0[n1];
	MN1 = nn1*Ngrid2*Ngrid3;
	MN0 = i*Ngrid2*Ngrid3;

	for (n2=0; n2<Ngrid2; n2++){
	  MN2 = n2*Ngrid3;
	  for (n3=0; n3<Ngrid3; n3++){
	    dVHart_Grid[MN1 + MN2 + n3] = tmp_array1[MN0 + MN2 + n3];
	  }
	}
      }
      free(tmp_array1);
    }

    if (Num_ISnd_Grid1[IDS]!=0){
      MPI_Wait(&request,&stat);
      free(tmp_array0);
    }
  }

  /* use own potential */
  for (n1=0; n1<My_NGrid1_Poisson; n1++){
    nn1 = Start_Grid1[myid] + n1;
    nnn1 = My_Cell0[nn1];
    if (My_Cell0[nn1]!=-1){
      MN1 = nnn1*Ngrid2*Ngrid3; 
      for (n2=0; n2<Ngrid2; n2++){
	MN2 = n2*Ngrid3;
	for (n3=0; n3<Ngrid3; n3++){
	  MN = MN1 + MN2 + n3;
	  dVHart_Grid[MN] = ReV1[n1][n2][n3];
	}
      }
    }
  }

  /****************************************************
               apply the gate voltage
  ****************************************************/

  if      (1.0<fabs(tv[1][1])) ip = 1;
  else if (1.0<fabs(tv[1][2])) ip = 2;
  else if (1.0<fabs(tv[1][3])) ip = 3;

  cc = 0.0;
  for (i=1; i<=Catomnum; i++){
    cc += Gxyz[Latomnum+i][ip];
  }  
  cc /= (double)Catomnum;

  a = 1.0*( fabs(tv[1][ip]) - fabs(Left_tv[1][ip]) - fabs(Right_tv[1][ip]) );

  for (i=0; i<Num_Cells0; i++){

    ri = My_Cell1[i];
    cv = (double)ri*length_gtv[1] + Grid_Origin[ip];

    gateV = tran_gate_voltage*exp( -pow( (cv-cc)/a, 8.0) );

    for (j=0; j<Ngrid2; j++){
      for (k=0; k<Ngrid3; k++){

	MN = i*Ngrid2*Ngrid3 + j*Ngrid3 + k; 
	dVHart_Grid[MN] += gateV;

      }
    }
  }

  /****************************************************
                   freeing of arrays
  ****************************************************/

  free(DL);
  free(D);
  free(DU);
  free(B);

  /* for time */
  MPI_Barrier(mpi_comm_level1);
  dtime(&TEtime);
  time0 = TEtime - TStime;
  return time0;
}




double TRAN_Poisson_FD_TCN(double ***ReV1, double ***ImV1,
		           double ***ReV2, double ***ImV2)
{ 
  int i,k1,k2,k3,kk2,ip,ri,k,j;
  int MN0,MN1,MN2,MN,nnn1,n1,n2,n3,nn1;
  double time0;
  double tmp0,sk1,sk2,sk3,tot_den,sum;
  double x,y,z,Gx,Gy,Gz,Eden[2];
  double ReTmp,ImTmp,c_coe,p_coe;
  double da2,Gpara2,a,cc,cv,gateV;
  double TStime,TEtime;
  dcomplex *DL,*D,*DU,*B;
  INTEGER n,nrhs,ldb,info;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  double *tmp_array0;
  double *tmp_array1;

  MPI_Status stat;
  MPI_Request request;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (myid==Host_ID) printf("<TRAN_Poisson_FD_TCN>  Solving Poisson's equation...\n");

  MPI_Barrier(mpi_comm_level1);
  dtime(&TStime);

  /****************************************************
                 allocation of arrays
  ****************************************************/

  DL = (dcomplex*)malloc(sizeof(dcomplex)*Ngrid1);
  D  = (dcomplex*)malloc(sizeof(dcomplex)*Ngrid1);
  DU = (dcomplex*)malloc(sizeof(dcomplex)*Ngrid1);
  B  = (dcomplex*)malloc(sizeof(dcomplex)*Ngrid1);

  /****************************************************
          FFT of charge density on the b-c plane
  ****************************************************/
  
  FFT2D_Density(0,ReV1,ImV1,ReV2,ImV2);

  /****************************************************
          solve finite difference equations
  ****************************************************/
        
  da2 = Dot_Product(gtv[1],gtv[1]);

  for (k2=0; k2<My_NGrid2_Poisson; k2++){

    kk2 = k2 + Start_Grid2[myid];

    if (kk2<Ngrid2/2) sk2 = (double)kk2;
    else              sk2 = (double)(kk2 - Ngrid2);

    for (k3=0; k3<Ngrid3; k3++){

      if (k3<Ngrid3/2) sk3 = (double)k3;
      else             sk3 = (double)(k3 - Ngrid3);

      Gx = sk2*rtv[2][1] + sk3*rtv[3][1];
      Gy = sk2*rtv[2][2] + sk3*rtv[3][2]; 
      Gz = sk2*rtv[2][3] + sk3*rtv[3][3];
      Gpara2 = Gx*Gx + Gy*Gy + Gz*Gz;

      for (k1=0; k1<(Ngrid1-1); k1++){
        DL[k1].r = 1.0;
        DL[k1].i = 0.0;
        DU[k1].r = 1.0;
        DU[k1].i = 0.0;
      }         

      for (k1=0; k1<(Ngrid1-0); k1++){
        D[k1].r = -2.0 - da2*Gpara2;
        D[k1].i = 0.0;
      }

      /* satisfy a charge neutrality condition */ 

      if (kk2==0 && k3==0){

        sum = 0.0;
        for (k1=0; k1<Ngrid1; k1++){
          sum += ReV2[k2][k1][k3];
        }
        
        for (k1=0; k1<Ngrid1; k1++){
          ReV2[k2][k1][k3] -= sum/(double)Ngrid1;
        }
      }

      /* set B */

      for (k1=0; k1<(Ngrid1-0); k1++){
        B[k1].r = -4.0*PI*da2*ReV2[k2][k1][k3];
        B[k1].i = -4.0*PI*da2*ImV2[k2][k1][k3];
      }

      /* add the boundary condition */

      B[0       ].r -= VHart_Boundary[0][Ngrid1_e[0]-1][kk2][k3].r;
      B[0       ].i -= VHart_Boundary[0][Ngrid1_e[0]-1][kk2][k3].i;
      B[Ngrid1-1].r -= VHart_Boundary[1][0            ][kk2][k3].r;
      B[Ngrid1-1].i -= VHart_Boundary[1][0            ][kk2][k3].i;

      /* solve the linear equation */ 

      n = Ngrid1-0;
      nrhs = 1;
      ldb = Ngrid1-0;
     
      F77_NAME(zgtsv,ZGTSV)(&n, &nrhs, DL, D, DU, B, &ldb, &info);

      /* store B to ReV2 and ImV2 */

      for (k1=0; k1<(Ngrid1-0); k1++){
        ReV2[k2][k1][k3] = B[k1].r;
        ImV2[k2][k1][k3] = B[k1].i;
      }

      /* set the boundary for the left part */

      /*
      ReV2[k2][0][k3] = VHart_Boundary[0][kk2][k3].r;
      ImV2[k2][0][k3] = VHart_Boundary[0][kk2][k3].i;
      */

    } /* k3 */
  } /* k2 */

  /****************************************************
   MPI communication: ReV2 and ImV2 into ReV1 and ImV1
  ****************************************************/
 
  MPI_Communication_F1_F2(2, 1, Ngrid2, Ngrid1, Ngrid3,
                          My_NGrid2_Poisson, My_NGrid1_Poisson,
                          ReV2, ImV2, ReV1, ImV1);

  /****************************************************
      Inverse FFT of ReV1 and ImV1 on the b-c plane
  ****************************************************/

  FFT2D_Poisson(1, 2, Ngrid1, Ngrid2, Ngrid3,
                My_NGrid1_Poisson, My_NGrid2_Poisson,
                1, 1, ReV1, ImV1, ReV2, ImV2);

  /****************************************************
   find the difference Hartree potential in real space
  ****************************************************/

  /* initialize */

  for (MN=0; MN<My_NumGrid1; MN++){
    dVHart_Grid[MN] = 0.0;
  }

  /* use their potential using MPI */

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    /* Isend */
    if (Num_ISnd_Grid1[IDS]!=0){
      tmp_array0 = (double*)malloc(sizeof(double)*Num_ISnd_Grid1[IDS]*Ngrid2*Ngrid3);

      for (i=0; i<Num_ISnd_Grid1[IDS]; i++){ 
	n1 = ISnd_Grid1[IDS][i];
	nn1 = n1 - Start_Grid1[myid];
	MN1 = i*Ngrid2*Ngrid3;
	for (n2=0; n2<Ngrid2; n2++){
	  MN2 = n2*Ngrid3;

          for (n3=0; n3<Ngrid3; n3++){
	    MN = MN1 + MN2 + n3;
	    tmp_array0[MN] = ReV1[nn1][n2][n3];
	  }
	}
      } 

      tag = 999;
      MPI_Isend(&tmp_array0[0], Num_ISnd_Grid1[IDS]*Ngrid2*Ngrid3,
		MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);
    }

    /* Recv */
    if (Num_IRcv_Grid1[IDR]!=0){

      tmp_array1 = (double*)malloc(sizeof(double)*Num_IRcv_Grid1[IDR]*Ngrid2*Ngrid3); 

      tag = 999;
      MPI_Recv(&tmp_array1[0], Num_IRcv_Grid1[IDR]*Ngrid2*Ngrid3,
	       MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

      for (i=0; i<Num_IRcv_Grid1[IDR]; i++){ 
	n1 = IRcv_Grid1[IDR][i];
	nn1 = My_Cell0[n1];
	MN1 = nn1*Ngrid2*Ngrid3;
	MN0 = i*Ngrid2*Ngrid3;

	for (n2=0; n2<Ngrid2; n2++){
	  MN2 = n2*Ngrid3;
	  for (n3=0; n3<Ngrid3; n3++){
	    dVHart_Grid[MN1 + MN2 + n3] = tmp_array1[MN0 + MN2 + n3];
	  }
	}
      }
      free(tmp_array1);
    }

    if (Num_ISnd_Grid1[IDS]!=0){
      MPI_Wait(&request,&stat);
      free(tmp_array0);
    }
  }

  /* use own potential */
  for (n1=0; n1<My_NGrid1_Poisson; n1++){
    nn1 = Start_Grid1[myid] + n1;
    nnn1 = My_Cell0[nn1];
    if (My_Cell0[nn1]!=-1){
      MN1 = nnn1*Ngrid2*Ngrid3; 
      for (n2=0; n2<Ngrid2; n2++){
	MN2 = n2*Ngrid3;
	for (n3=0; n3<Ngrid3; n3++){
	  MN = MN1 + MN2 + n3;
	  dVHart_Grid[MN] = ReV1[n1][n2][n3];
	}
      }
    }
  }

  /****************************************************
               apply the gate voltage
  ****************************************************/

  if      (1.0<fabs(tv[1][1])) ip = 1;
  else if (1.0<fabs(tv[1][2])) ip = 2;
  else if (1.0<fabs(tv[1][3])) ip = 3;

  cc = 0.0;
  for (i=1; i<=Catomnum; i++){
    cc += Gxyz[Latomnum+i][ip];
  }  
  cc /= (double)Catomnum;

  a = 1.0*( fabs(tv[1][ip]) - fabs(Left_tv[1][ip]) - fabs(Right_tv[1][ip]) );

  for (i=0; i<Num_Cells0; i++){

    ri = My_Cell1[i];
    cv = (double)ri*length_gtv[1] + Grid_Origin[ip];

    gateV = tran_gate_voltage*exp( -pow( (cv-cc)/a, 8.0) );

    for (j=0; j<Ngrid2; j++){
      for (k=0; k<Ngrid3; k++){

	MN = i*Ngrid2*Ngrid3 + j*Ngrid3 + k; 
	dVHart_Grid[MN] += gateV;

      }
    }
  }

  /****************************************************
                   freeing of arrays
  ****************************************************/

  free(DL);
  free(D);
  free(DU);
  free(B);

  /* for time */
  MPI_Barrier(mpi_comm_level1);
  dtime(&TEtime);
  time0 = TEtime - TStime;
  return time0;
}



double TRAN_Poisson_FD_LCN(double ***ReV1, double ***ImV1,
		           double ***ReV2, double ***ImV2)
{ 
  int i,k1,k2,k3,kk2,ip,ri,k,j;
  int MN0,MN1,MN2,MN,nnn1,n1,n2,n3,nn1;
  double time0;
  double tmp0,sk1,sk2,sk3,tot_den,sum;
  double x,y,z,Gx,Gy,Gz,Eden[2];
  double ReTmp,ImTmp,c_coe,p_coe;
  double da2,Gpara2,a,cc,cv,gateV;
  double TStime,TEtime;
  dcomplex *DL,*D,*DU,*B;
  INTEGER n,nrhs,ldb,info;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  double *tmp_array0;
  double *tmp_array1;

  MPI_Status stat;
  MPI_Request request;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (myid==Host_ID) printf("<TRAN_Poisson_FD_LCN>  Solving Poisson's equation...\n");

  MPI_Barrier(mpi_comm_level1);
  dtime(&TStime);

  /****************************************************
                 allocation of arrays
  ****************************************************/

  DL = (dcomplex*)malloc(sizeof(dcomplex)*Ngrid1);
  D  = (dcomplex*)malloc(sizeof(dcomplex)*Ngrid1);
  DU = (dcomplex*)malloc(sizeof(dcomplex)*Ngrid1);
  B  = (dcomplex*)malloc(sizeof(dcomplex)*Ngrid1);

  /****************************************************
          FFT of charge density on the b-c plane
  ****************************************************/
  
  FFT2D_Density(0,ReV1,ImV1,ReV2,ImV2);

  /****************************************************
          solve finite difference equations
  ****************************************************/
        
  da2 = Dot_Product(gtv[1],gtv[1]);

  for (k2=0; k2<My_NGrid2_Poisson; k2++){

    kk2 = k2 + Start_Grid2[myid];

    if (kk2<Ngrid2/2) sk2 = (double)kk2;
    else              sk2 = (double)(kk2 - Ngrid2);

    for (k3=0; k3<Ngrid3; k3++){

      if (k3<Ngrid3/2) sk3 = (double)k3;
      else             sk3 = (double)(k3 - Ngrid3);

      Gx = sk2*rtv[2][1] + sk3*rtv[3][1];
      Gy = sk2*rtv[2][2] + sk3*rtv[3][2]; 
      Gz = sk2*rtv[2][3] + sk3*rtv[3][3];
      Gpara2 = Gx*Gx + Gy*Gy + Gz*Gz;

      for (k1=0; k1<(Ngrid1-1); k1++){
        DL[k1].r = 1.0;
        DL[k1].i = 0.0;
        DU[k1].r = 1.0;
        DU[k1].i = 0.0;
      }         

      for (k1=0; k1<(Ngrid1-0); k1++){
        D[k1].r = -2.0 - da2*Gpara2;
        D[k1].i = 0.0;
      }

      /* satisfy a charge neutrality condition */ 

      if (kk2==0 && k3==0){

        for (k1=0; k1<Ngrid1; k1++){
          ReV2[k2][k1][k3] = 0.0;
        }
      }

      /* set B */

      for (k1=0; k1<(Ngrid1-0); k1++){
        B[k1].r = -4.0*PI*da2*ReV2[k2][k1][k3];
        B[k1].i = -4.0*PI*da2*ImV2[k2][k1][k3];
      }

      /* add the boundary condition */

      B[0       ].r -= VHart_Boundary[0][Ngrid1_e[0]-1][kk2][k3].r;
      B[0       ].i -= VHart_Boundary[0][Ngrid1_e[0]-1][kk2][k3].i;
      B[Ngrid1-1].r -= VHart_Boundary[1][0            ][kk2][k3].r;
      B[Ngrid1-1].i -= VHart_Boundary[1][0            ][kk2][k3].i;

      /* solve the linear equation */ 

      n = Ngrid1-0;
      nrhs = 1;
      ldb = Ngrid1-0;
     
      F77_NAME(zgtsv,ZGTSV)(&n, &nrhs, DL, D, DU, B, &ldb, &info);

      /* store B to ReV2 and ImV2 */

      for (k1=0; k1<(Ngrid1-0); k1++){
        ReV2[k2][k1][k3] = B[k1].r;
        ImV2[k2][k1][k3] = B[k1].i;
      }

      /* set the boundary for the left part */

      /*
      ReV2[k2][0][k3] = VHart_Boundary[0][kk2][k3].r;
      ImV2[k2][0][k3] = VHart_Boundary[0][kk2][k3].i;
      */

    } /* k3 */
  } /* k2 */

  /****************************************************
   MPI communication: ReV2 and ImV2 into ReV1 and ImV1
  ****************************************************/
 
  MPI_Communication_F1_F2(2, 1, Ngrid2, Ngrid1, Ngrid3,
                          My_NGrid2_Poisson, My_NGrid1_Poisson,
                          ReV2, ImV2, ReV1, ImV1);

  /****************************************************
      Inverse FFT of ReV1 and ImV1 on the b-c plane
  ****************************************************/

  FFT2D_Poisson(1, 2, Ngrid1, Ngrid2, Ngrid3,
                My_NGrid1_Poisson, My_NGrid2_Poisson,
                1, 1, ReV1, ImV1, ReV2, ImV2);

  /****************************************************
   find the difference Hartree potential in real space
  ****************************************************/

  /* initialize */

  for (MN=0; MN<My_NumGrid1; MN++){
    dVHart_Grid[MN] = 0.0;
  }

  /* use their potential using MPI */

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    /* Isend */
    if (Num_ISnd_Grid1[IDS]!=0){
      tmp_array0 = (double*)malloc(sizeof(double)*Num_ISnd_Grid1[IDS]*Ngrid2*Ngrid3);

      for (i=0; i<Num_ISnd_Grid1[IDS]; i++){ 
	n1 = ISnd_Grid1[IDS][i];
	nn1 = n1 - Start_Grid1[myid];
	MN1 = i*Ngrid2*Ngrid3;
	for (n2=0; n2<Ngrid2; n2++){
	  MN2 = n2*Ngrid3;

          for (n3=0; n3<Ngrid3; n3++){
	    MN = MN1 + MN2 + n3;
	    tmp_array0[MN] = ReV1[nn1][n2][n3];
	  }
	}
      } 

      tag = 999;
      MPI_Isend(&tmp_array0[0], Num_ISnd_Grid1[IDS]*Ngrid2*Ngrid3,
		MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);
    }

    /* Recv */
    if (Num_IRcv_Grid1[IDR]!=0){

      tmp_array1 = (double*)malloc(sizeof(double)*Num_IRcv_Grid1[IDR]*Ngrid2*Ngrid3); 

      tag = 999;
      MPI_Recv(&tmp_array1[0], Num_IRcv_Grid1[IDR]*Ngrid2*Ngrid3,
	       MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

      for (i=0; i<Num_IRcv_Grid1[IDR]; i++){ 
	n1 = IRcv_Grid1[IDR][i];
	nn1 = My_Cell0[n1];
	MN1 = nn1*Ngrid2*Ngrid3;
	MN0 = i*Ngrid2*Ngrid3;

	for (n2=0; n2<Ngrid2; n2++){
	  MN2 = n2*Ngrid3;
	  for (n3=0; n3<Ngrid3; n3++){
	    dVHart_Grid[MN1 + MN2 + n3] = tmp_array1[MN0 + MN2 + n3];
	  }
	}
      }
      free(tmp_array1);
    }

    if (Num_ISnd_Grid1[IDS]!=0){
      MPI_Wait(&request,&stat);
      free(tmp_array0);
    }
  }

  /* use own potential */
  for (n1=0; n1<My_NGrid1_Poisson; n1++){
    nn1 = Start_Grid1[myid] + n1;
    nnn1 = My_Cell0[nn1];
    if (My_Cell0[nn1]!=-1){
      MN1 = nnn1*Ngrid2*Ngrid3; 
      for (n2=0; n2<Ngrid2; n2++){
	MN2 = n2*Ngrid3;
	for (n3=0; n3<Ngrid3; n3++){
	  MN = MN1 + MN2 + n3;
	  dVHart_Grid[MN] = ReV1[n1][n2][n3];
	}
      }
    }
  }

  /****************************************************
               apply the gate voltage
  ****************************************************/

  if      (1.0<fabs(tv[1][1])) ip = 1;
  else if (1.0<fabs(tv[1][2])) ip = 2;
  else if (1.0<fabs(tv[1][3])) ip = 3;

  cc = 0.0;
  for (i=1; i<=Catomnum; i++){
    cc += Gxyz[Latomnum+i][ip];
  }  
  cc /= (double)Catomnum;

  a = 1.0*( fabs(tv[1][ip]) - fabs(Left_tv[1][ip]) - fabs(Right_tv[1][ip]) );

  for (i=0; i<Num_Cells0; i++){

    ri = My_Cell1[i];
    cv = (double)ri*length_gtv[1] + Grid_Origin[ip];

    gateV = tran_gate_voltage*exp( -pow( (cv-cc)/a, 8.0) );

    for (j=0; j<Ngrid2; j++){
      for (k=0; k<Ngrid3; k++){

	MN = i*Ngrid2*Ngrid3 + j*Ngrid3 + k; 
	dVHart_Grid[MN] += gateV;

      }
    }
  }

  /****************************************************
                   freeing of arrays
  ****************************************************/

  free(DL);
  free(D);
  free(DU);
  free(B);

  /* for time */
  MPI_Barrier(mpi_comm_level1);
  dtime(&TEtime);
  time0 = TEtime - TStime;
  return time0;
}



void FFT2D_Poisson(int NG1, int NG2,
                   int Ng1, int Ng2,  int Ng3,
		   int My_Ng1, int My_Ng2,
		   int sgn2, int sgn3,
		   double ***ReF1, double ***ImF1,
		   double ***ReF2, double ***ImF2)
{
  int k1,k2,k3;
  int n1,n2,n3;
  int nn1,nn2,MN,MN1,MN2,MN0;
  int Num1,RNum2,Num2,NsizeS,NsizeR;
  double *tmp_array0,*tmp_array1;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  double Stime_proc, Etime_proc;

  MPI_Status stat;
  MPI_Request request;
  fftw_complex *in, *out;
  fftw_plan p;
  
  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
    allocation of arrays:

    fftw_complex  in[List_YOUSO[17]];
    fftw_complex out[List_YOUSO[17]];
  ****************************************************/

#ifdef fftw2
  in  = (fftw_complex*)malloc(sizeof(fftw_complex)*List_YOUSO[17]); 
  out = (fftw_complex*)malloc(sizeof(fftw_complex)*List_YOUSO[17]); 
#else
  in  = fftw_malloc(sizeof(fftw_complex)*List_YOUSO[17]); 
  out = fftw_malloc(sizeof(fftw_complex)*List_YOUSO[17]); 
#endif

  /*------------------ S n1,n2,k3 ------------------*/

  if (measure_time==1) dtime(&Stime_proc);

#ifdef fftw2
  p = fftw_create_plan(Ng3, sgn3, FFTW_ESTIMATE);
#else
  p = fftw_plan_dft_1d(Ng3,in,out,sgn3,FFTW_ESTIMATE);
#endif

  for (n1=0; n1<My_Ng1; n1++){
    for (n2=0; n2<Ng2; n2++){
      
      for (n3=0; n3<Ng3; n3++){

#ifdef fftw2
        c_re(in[n3]) = ReF1[n1][n2][n3];
        c_im(in[n3]) = ImF1[n1][n2][n3];
#else
        in[n3][0] = ReF1[n1][n2][n3];
        in[n3][1] = ImF1[n1][n2][n3];
#endif

      }

#ifdef fftw2
      fftw_one(p, in, out);
#else
      fftw_execute(p);
#endif
      
      for (k3=0; k3<Ng3; k3++){

#ifdef fftw2
        ReF1[n1][n2][k3] = c_re(out[k3]);
        ImF1[n1][n2][k3] = c_im(out[k3]);
#else
        ReF1[n1][n2][k3] = out[k3][0];
        ImF1[n1][n2][k3] = out[k3][1];
#endif

      }
      
    } 
  }   
  fftw_destroy_plan(p);  

  /*------------------ T n1,k2,k3 ------------------*/

#ifdef fftw2
  p = fftw_create_plan(Ng2, sgn2, FFTW_ESTIMATE);
#else
  p = fftw_plan_dft_1d(Ng2,in,out,sgn2,FFTW_ESTIMATE);
#endif

  for (n1=0; n1<My_Ng1; n1++){
    for (k3=0; k3<Ng3; k3++){

      for (n2=0; n2<Ng2; n2++){

#ifdef fftw2
        c_re(in[n2]) = ReF1[n1][n2][k3];
        c_im(in[n2]) = ImF1[n1][n2][k3];
#else
        in[n2][0] = ReF1[n1][n2][k3];
        in[n2][1] = ImF1[n1][n2][k3];
#endif

      }

#ifdef fftw2
      fftw_one(p, in, out);
#else
      fftw_execute(p);
#endif

      for (k2=0; k2<Ng2; k2++){

#ifdef fftw2
        ReF1[n1][k2][k3] = c_re(out[k2]);
        ImF1[n1][k2][k3] = c_im(out[k2]);
#else
        ReF1[n1][k2][k3] = out[k2][0];
        ImF1[n1][k2][k3] = out[k2][1];
#endif

      }

    }
  }

  fftw_destroy_plan(p);  

  if (measure_time==1){
    dtime(&Etime_proc);
    printf("myid=%2d  Time yz FFT    = %15.12f\n",myid,Etime_proc-Stime_proc);
  }

  /****************************************************
    freeing of arrays:

    fftw_complex  in[List_YOUSO[17]];
    fftw_complex out[List_YOUSO[17]];
  ****************************************************/

#ifdef fftw2
  free(in);
  free(out);
#else
  fftw_free(in);
  fftw_free(out);
#endif

}




void MPI_Communication_F1_F2(int NG1, int NG2,
			     int Ng1, int Ng2,  int Ng3,
			     int My_Ng1, int My_Ng2,
			     double ***ReF1, double ***ImF1,
			     double ***ReF2, double ***ImF2)
{  
  int k1,k2,k3;
  int n1,n2,n3;
  int nn1,nn2,MN,MN1,MN2,MN0;
  int Num1,RNum2,Num2,NsizeS,NsizeR;
  double *tmp_array0,*tmp_array1;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  double Stime_proc, Etime_proc;

  MPI_Status stat;
  MPI_Request request;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  
  /*********************************************
      MPI:   ReF1 and ImF1 -> ReF2 and ImF2
  *********************************************/

  if (measure_time==1){
    MPI_Barrier(mpi_comm_level1);
    dtime(&Stime_proc);
  }

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){

      /* sending */ 

      if (NG2==1)
        Num2 = End_Grid1[IDS] - Start_Grid1[IDS] + 1;
      else if (NG2==2)
        Num2 = End_Grid2[IDS] - Start_Grid2[IDS] + 1;

      NsizeS = My_Ng1*Num2*Ng3;
      tmp_array0 = (double*)malloc(sizeof(double)*2*NsizeS); 

      for (n1=0; n1<My_Ng1; n1++){
        MN1 = n1*Num2*Ng3;
        for (n2=0; n2<Num2; n2++){
          MN2 = n2*Ng3;

          if (NG2==1)
            nn2 = n2 + Start_Grid1[IDS];
          else if (NG2==2)
            nn2 = n2 + Start_Grid2[IDS];

          MN0 = MN1 + MN2; 
          for (n3=0; n3<Ng3; n3++){
            MN = MN0 + n3;
            tmp_array0[MN] = ReF1[n1][nn2][n3];            
 	  }

          MN0 = MN1 + MN2 + NsizeS;
          for (n3=0; n3<Ng3; n3++){
            MN = MN0 + n3;
            tmp_array0[MN] = ImF1[n1][nn2][n3];
 	  }

	}
      }

      tag = 999;
      MPI_Isend(&tmp_array0[0], 2*NsizeS,
                MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);

      /* receiving */ 

      if (NG1==1)
        Num1 = End_Grid1[IDR] - Start_Grid1[IDR] + 1;
      else if (NG1==2)
        Num1 = End_Grid2[IDR] - Start_Grid2[IDR] + 1;

      if (NG2==1)
        RNum2 = End_Grid1[myid] - Start_Grid1[myid] + 1;
      else if (NG2==2)
        RNum2 = End_Grid2[myid] - Start_Grid2[myid] + 1;

      NsizeR = Num1*RNum2*Ng3; 
      tmp_array1 = (double*)malloc(sizeof(double)*2*NsizeR);

      tag = 999;
      MPI_Recv(&tmp_array1[0], 2*NsizeR,
                MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

      for (n1=0; n1<Num1; n1++){
        MN1 = n1*RNum2*Ng3;

        if (NG1==1)
          nn1 = n1 + Start_Grid1[IDR];
        else if (NG1==2)
          nn1 = n1 + Start_Grid2[IDR];

        for (n2=0; n2<RNum2; n2++){
          MN2 = n2*Ng3;

          MN0 = MN1 + MN2;
          for (n3=0; n3<Ng3; n3++){
            MN = MN0 + n3;
            ReF2[n2][nn1][n3] = tmp_array1[MN];
 	  }

          MN0 = MN1 + MN2 + NsizeR;
          for (n3=0; n3<Ng3; n3++){
            MN = MN0 + n3;
            ImF2[n2][nn1][n3] = tmp_array1[MN];
 	  }

	}
      }      

      MPI_Wait(&request,&stat);
      free(tmp_array0);
      free(tmp_array1);
  
    }

    else{

      if (NG2==1)
        Num2 = End_Grid1[myid] - Start_Grid1[myid] + 1;
      else if (NG2==2)
        Num2 = End_Grid2[myid] - Start_Grid2[myid] + 1;

      for (n1=0; n1<My_Ng1; n1++){

        if (NG1==1)
          nn1 = n1 + Start_Grid1[myid];
        else if (NG1==2)
          nn1 = n1 + Start_Grid2[myid];

        for (n2=0; n2<Num2; n2++){

          if (NG2==1)
            nn2 = n2 + Start_Grid1[myid];
          else if (NG2==2)
            nn2 = n2 + Start_Grid2[myid];

          for (n3=0; n3<Ng3; n3++){
            ReF2[n2][nn1][n3] = ReF1[n1][nn2][n3];            
            ImF2[n2][nn1][n3] = ImF1[n1][nn2][n3];
 	  }
	}
      }
    }

  }

  if (measure_time==1){
    MPI_Barrier(mpi_comm_level1);
    dtime(&Etime_proc);
    printf("myid=%2d  Time transpose = %15.12f\n",myid,Etime_proc-Stime_proc);
  }
  
} 





void FFT2D_Density(int den_flag, 
                   double ***ReV1, double ***ImV1,
                   double ***ReV2, double ***ImV2)
{
  int ct_AN,n1,n2,n3,k1,k2,k3,kk2,N3[4];
  int N,i,j;
  int nn0,nn1,nnn1,MN,MN0,MN1,MN2;
  double time0; 
  int h_AN,Gh_AN,spin,Nc,GNc,GN,Nh,Nog;
  double tmp0,tmp1,sk1,sk2,sk3,tot_den;
  double x,y,z,Gx,Gy,Gz,Eden[2],DenA,G2;
  double ReTmp,ImTmp,c_coe,p_coe;
  double *tmp_array0;
  double *tmp_array1;
  double TStime,TEtime;
  int numprocs,myid,tag=999,ID,IDS,IDR;

  MPI_Status stat;
  MPI_Request request;

  if (atomnum<=MYID_MPI_COMM_WORLD) return;
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /* initialize */
  for (n1=0; n1<My_NGrid1_Poisson; n1++){
    for (n2=0; n2<Ngrid2; n2++){
      for (n3=0; n3<Ngrid3; n3++){
        ReV1[n1][n2][n3] = 0.0;
        ImV1[n1][n2][n3] = 0.0;
      }
    }
  }

  /* use their densities using MPI */
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    /* Isend */
    if (Num_Snd_Grid1[IDS]!=0){

      tmp_array0 = (double*)malloc(sizeof(double)*Num_Snd_Grid1[IDS]*Ngrid2*Ngrid3); 

      for (i=0; i<Num_Snd_Grid1[IDS]; i++){ 
	n1 = Snd_Grid1[IDS][i];
        nn1 = My_Cell0[n1];
        MN1 = nn1*Ngrid2*Ngrid3;
        MN0 = i*Ngrid2*Ngrid3;

        for (n2=0; n2<Ngrid2; n2++){
          MN2 = n2*Ngrid3;

          switch(den_flag) {

          case 0:

	    for (n3=0; n3<Ngrid3; n3++){
	      MN = MN1 + MN2 + n3;
	      DenA = 2.0*ADensity_Grid[MN];
	      tmp_array0[MN0+MN2+n3] = Density_Grid[0][MN] + Density_Grid[1][MN] - DenA;
	    }

          break;

          case 1:

	    for (n3=0; n3<Ngrid3; n3++){
	      MN = MN1 + MN2 + n3;
	      tmp_array0[MN0+MN2+n3] = Density_Grid[0][MN];
	    }

          break;

          case 2:

	    for (n3=0; n3<Ngrid3; n3++){
	      MN = MN1 + MN2 + n3;
	      tmp_array0[MN0+MN2+n3] = Density_Grid[1][MN];
	    }

          break;

          case 3:

	    for (n3=0; n3<Ngrid3; n3++){
	      MN = MN1 + MN2 + n3;
	      tmp_array0[MN0+MN2+n3] = 2.0*ADensity_Grid[MN];
	    }

          break;

	  }

	}
      }

      tag = 999;
      MPI_Isend(&tmp_array0[0], Num_Snd_Grid1[IDS]*Ngrid2*Ngrid3,
                MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);
    }

    /* Recv */

    if (Num_Rcv_Grid1[IDR]!=0){

      tmp_array1 = (double*)malloc(sizeof(double)*Num_Rcv_Grid1[IDR]*Ngrid2*Ngrid3); 

      tag = 999;
      MPI_Recv(&tmp_array1[0], Num_Rcv_Grid1[IDR]*Ngrid2*Ngrid3,
                MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

      for (i=0; i<Num_Rcv_Grid1[IDR]; i++){ 
	n1 = Rcv_Grid1[IDR][i];
	nn1 = My_Cell0[n1];
	nn0 = n1 - Start_Grid1[myid];
	MN0 = i*Ngrid2*Ngrid3;
	for (n2=0; n2<Ngrid2; n2++){
	  MN2 = n2*Ngrid3;
	  for (n3=0; n3<Ngrid3; n3++){
	    MN = MN0 + MN2 + n3;
	    ReV1[nn0][n2][n3] = tmp_array1[MN];
	    ImV1[nn0][n2][n3] = 0.0;
	  }
	}
      }

      free(tmp_array1);
    }

    if (Num_Snd_Grid1[IDS]!=0){
      MPI_Wait(&request,&stat);
      free(tmp_array0);
    }

  }

  /* use own densities */
  for (n1=Start_Grid1[myid]; n1<=End_Grid1[myid]; n1++){
    nn1 = My_Cell0[n1];
    nn0 = n1 - Start_Grid1[myid]; 
    if (nn1!=-1){
      MN1 = nn1*Ngrid2*Ngrid3;
      for (n2=0; n2<Ngrid2; n2++){
        MN2 = n2*Ngrid3;

        switch(den_flag) {

	case 0:

	  for (n3=0; n3<Ngrid3; n3++){
	    MN = MN1 + MN2 + n3;
	    DenA = 2.0*ADensity_Grid[MN];
	    ReV1[nn0][n2][n3] = Density_Grid[0][MN] + Density_Grid[1][MN] - DenA;
	    ImV1[nn0][n2][n3] = 0.0;
	  }

        break;

	case 1:

	  for (n3=0; n3<Ngrid3; n3++){
	    MN = MN1 + MN2 + n3;
	    ReV1[nn0][n2][n3] = Density_Grid[0][MN];
	    ImV1[nn0][n2][n3] = 0.0;
	  }

        break;

	case 2:

	  for (n3=0; n3<Ngrid3; n3++){
	    MN = MN1 + MN2 + n3;
	    ReV1[nn0][n2][n3] = Density_Grid[1][MN];
	    ImV1[nn0][n2][n3] = 0.0;
	  }

        break;

	case 3:

	  for (n3=0; n3<Ngrid3; n3++){
	    MN = MN1 + MN2 + n3;
	    ReV1[nn0][n2][n3] = 2.0*ADensity_Grid[MN];
	    ImV1[nn0][n2][n3] = 0.0;
	  }

        break;

	}

      }    
    }
  }

  /****************************************************
                       FFT of Dens
  ****************************************************/

  FFT2D_Poisson(1, 2, Ngrid1, Ngrid2, Ngrid3,
                My_NGrid1_Poisson, My_NGrid2_Poisson,
                -1, -1, ReV1, ImV1, ReV2, ImV2);

  MPI_Communication_F1_F2(1, 2, Ngrid1, Ngrid2, Ngrid3,
                          My_NGrid1_Poisson, My_NGrid2_Poisson,
                          ReV1, ImV1, ReV2, ImV2);

  tmp0 = 1.0/(double)(Ngrid2*Ngrid3);
  for (k2=0; k2<My_NGrid2_Poisson; k2++){
    for (k1=0; k1<Ngrid1; k1++){
      for (k3=0; k3<Ngrid3; k3++){
        ReV2[k2][k1][k3] *= tmp0;
        ImV2[k2][k1][k3] *= tmp0;
      }
    }
  }

}



void Get_Value_inReal2D(double ***ReV2, double ***ImV2, double ***ReV1, double ***ImV1, double *Value_Grid)
{
  int i,MN,MN0,MN1,MN2,nn1,nnn1,n1,n2,n3;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  double *tmp_array0;
  double *tmp_array1;

  MPI_Status stat;
  MPI_Request request;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
   MPI communication: ReV2 and ImV2 into ReV1 and ImV1
  ****************************************************/

  MPI_Communication_F1_F2(2, 1, Ngrid2, Ngrid1, Ngrid3,
                          My_NGrid2_Poisson, My_NGrid1_Poisson,
                          ReV2, ImV2, ReV1, ImV1);

  /****************************************************
      Inverse FFT of ReV1 and ImV1 on the b-c plane
  ****************************************************/

  FFT2D_Poisson(1, 2, Ngrid1, Ngrid2, Ngrid3,
                My_NGrid1_Poisson, My_NGrid2_Poisson,
                1, 1, ReV1, ImV1, ReV2, ImV2);

  /****************************************************
   find the difference Hartree potential in real space
  ****************************************************/

  /* initialize */

  for (MN=0; MN<My_NumGrid1; MN++){
    dVHart_Grid[MN] = 0.0;
  }

  /* use their potential using MPI */

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    /* Isend */
    if (Num_ISnd_Grid1[IDS]!=0){
      tmp_array0 = (double*)malloc(sizeof(double)*Num_ISnd_Grid1[IDS]*Ngrid2*Ngrid3);

      for (i=0; i<Num_ISnd_Grid1[IDS]; i++){ 
	n1 = ISnd_Grid1[IDS][i];
	nn1 = n1 - Start_Grid1[myid];
	MN1 = i*Ngrid2*Ngrid3;
	for (n2=0; n2<Ngrid2; n2++){
	  MN2 = n2*Ngrid3;

          for (n3=0; n3<Ngrid3; n3++){
	    MN = MN1 + MN2 + n3;
	    tmp_array0[MN] = ReV1[nn1][n2][n3];
	  }
	}
      } 

      tag = 999;
      MPI_Isend(&tmp_array0[0], Num_ISnd_Grid1[IDS]*Ngrid2*Ngrid3,
		MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);
    }

    /* Recv */
    if (Num_IRcv_Grid1[IDR]!=0){

      tmp_array1 = (double*)malloc(sizeof(double)*Num_IRcv_Grid1[IDR]*Ngrid2*Ngrid3); 

      tag = 999;
      MPI_Recv(&tmp_array1[0], Num_IRcv_Grid1[IDR]*Ngrid2*Ngrid3,
	       MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

      for (i=0; i<Num_IRcv_Grid1[IDR]; i++){ 
	n1 = IRcv_Grid1[IDR][i];
	nn1 = My_Cell0[n1];
	MN1 = nn1*Ngrid2*Ngrid3;
	MN0 = i*Ngrid2*Ngrid3;

	for (n2=0; n2<Ngrid2; n2++){
	  MN2 = n2*Ngrid3;
	  for (n3=0; n3<Ngrid3; n3++){
	    Value_Grid[MN1 + MN2 + n3] = tmp_array1[MN0 + MN2 + n3];
	  }
	}
      }
      free(tmp_array1);
    }

    if (Num_ISnd_Grid1[IDS]!=0){
      MPI_Wait(&request,&stat);
      free(tmp_array0);
    }
  }

  /* use own potential */
  for (n1=0; n1<My_NGrid1_Poisson; n1++){
    nn1 = Start_Grid1[myid] + n1;
    nnn1 = My_Cell0[nn1];
    if (My_Cell0[nn1]!=-1){
      MN1 = nnn1*Ngrid2*Ngrid3; 
      for (n2=0; n2<Ngrid2; n2++){
	MN2 = n2*Ngrid3;
	for (n3=0; n3<Ngrid3; n3++){
	  MN = MN1 + MN2 + n3;
	  Value_Grid[MN] = ReV1[n1][n2][n3];
	}
      }
    }
  }

}
