/**********************************************************************
  Poisson.c:
  
     Poisson.c is a subrutine to solve Poisson's equation using
     fast Fourier transformation.

  Log of Poisson.c:

     22/Nov/2001  Released by T.Ozaki

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
#include "mpi.h"
#endif

#ifdef fftw2 
#include <fftw.h>
#else
#include <fftw3.h> 
#endif


double Poisson(int fft_charge_flag,
               double ***ReV1, double ***ImV1,
               double ***ReV2, double ***ImV2)
{ 
  int k1,k2,k3,kk2;
  double time0;
  double tmp0,sk1,sk2,sk3,tot_den;
  double x,y,z,Gx,Gy,Gz,Eden[2],DenA,G2;
  double ReTmp,ImTmp,c_coe,p_coe;
  double *tmp_array0;
  double *tmp_array1;
  double TStime,TEtime;
  int numprocs,myid,tag=999,ID,IDS,IDR;

  MPI_Status stat;
  MPI_Request request;

  /* MPI */
  if (atomnum<=MYID_MPI_COMM_WORLD) return 0.0;
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (myid==Host_ID && 0<level_stdout){
    printf("<Poisson>  Poisson's equation using FFT...\n");
  }

  MPI_Barrier(mpi_comm_level1);
  dtime(&TStime);

  /****************************************************
                 FFT of charge density 
  ****************************************************/

  if (fft_charge_flag==1) FFT_Density(0,ReV1,ImV1,ReV2,ImV2);

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

  /* for time */
  MPI_Barrier(mpi_comm_level1);
  dtime(&TEtime);
  time0 = TEtime - TStime;
  return time0;
}




void FFT_Poisson(int NG1, int NG2,
                 int Ng1, int Ng2,  int Ng3,
                 int My_Ng1, int My_Ng2,
                 int sgn1, int sgn2, int sgn3,
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

  /*********************************************
    MPI    ReF1 and ImF1 -> ReF2 and ImF2
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

  /*------------------ U k1,k2,k3 ------------------*/

  if (measure_time==1){
    dtime(&Stime_proc);
  }

#ifdef fftw2
  p = fftw_create_plan(Ng1, sgn1, FFTW_ESTIMATE);
#else
  p = fftw_plan_dft_1d(Ng1,in,out,sgn1,FFTW_ESTIMATE);
#endif

  for (k2=0; k2<My_Ng2; k2++){
    for (k3=0; k3<Ng3; k3++){

      for (n1=0; n1<Ng1; n1++){

#ifdef fftw2
        c_re(in[n1]) = ReF2[k2][n1][k3];
        c_im(in[n1]) = ImF2[k2][n1][k3];
#else
	in[n1][0] = ReF2[k2][n1][k3];
        in[n1][1] = ImF2[k2][n1][k3];
#endif

      }

#ifdef fftw2
      fftw_one(p, in, out);
#else
      fftw_execute(p);
#endif
      
      for (k1=0; k1<Ng1; k1++){

#ifdef fftw2
        ReF2[k2][k1][k3] = c_re(out[k1]);
        ImF2[k2][k1][k3] = c_im(out[k1]);
#else
        ReF2[k2][k1][k3] = out[k1][0];
        ImF2[k2][k1][k3] = out[k1][1];
#endif

      }
    }
  }
  fftw_destroy_plan(p);  
  fftw_cleanup();

  if (measure_time==1){
    dtime(&Etime_proc);
    printf("myid=%2d  Time x FFT     = %15.12f\n",myid,Etime_proc-Stime_proc);
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


void FFT_Density(int den_flag,
                 double ***ReV1, double ***ImV1,
                 double ***ReV2, double ***ImV2)
{
  int ct_AN,n1,n2,n3,k1,k2,k3,kk2,N3[4];
  int N,i,j;
  int nn0,nn1,nnn1,MN,MN0,MN1,MN2;
  double time0;
  int h_AN,Gh_AN,spin,Nc,GNc,GN,Nh,Nog,factor;
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

  if (den_flag!=4)      factor = 1;
  else if (den_flag==4) factor = 2;
  else {
    printf("invalid den_flag\n");
    MPI_Finalize();
    exit(0);
  }

  /* initialize */
  for (n1=0; n1<My_NGrid1_Poisson; n1++){
    for (n2=0; n2<Ngrid2; n2++){
      for (n3=0; n3<Ngrid3; n3++){
        ReV1[n1][n2][n3] = 0.0;
        ImV1[n1][n2][n3] = 0.0;
      }
    }
  }

  /* use their densities using MPI  */
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    /* Isend */
    if (Num_Snd_Grid1[IDS]!=0){

      tmp_array0 = (double*)malloc(sizeof(double)*Num_Snd_Grid1[IDS]*Ngrid2*Ngrid3*factor); 

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

          case 4:

            for (n3=0; n3<Ngrid3; n3++){
              MN = MN1 + MN2 + n3;
              tmp_array0[2*MN0+2*MN2+2*n3  ] = Density_Grid[2][MN];
              tmp_array0[2*MN0+2*MN2+2*n3+1] = Density_Grid[3][MN];
	    }

          break;

	  } /* switch(den_flag) */

	}
      }

      tag = 999;
      MPI_Isend(&tmp_array0[0], Num_Snd_Grid1[IDS]*Ngrid2*Ngrid3*factor,
                MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);
    }

    /* Recv */

    if (Num_Rcv_Grid1[IDR]!=0){

      tmp_array1 = (double*)malloc(sizeof(double)*Num_Rcv_Grid1[IDR]*Ngrid2*Ngrid3*factor); 

      tag = 999;
      MPI_Recv(&tmp_array1[0], Num_Rcv_Grid1[IDR]*Ngrid2*Ngrid3*factor,
                MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

      if (den_flag!=4){

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
      }

      else if (den_flag==4){

	for (i=0; i<Num_Rcv_Grid1[IDR]; i++){ 
	  n1 = Rcv_Grid1[IDR][i];
	  nn1 = My_Cell0[n1];
	  nn0 = n1 - Start_Grid1[myid];
	  MN0 = i*Ngrid2*Ngrid3;
	  for (n2=0; n2<Ngrid2; n2++){
	    MN2 = n2*Ngrid3;
	    for (n3=0; n3<Ngrid3; n3++){
	      MN = MN0 + MN2 + n3;
	      ReV1[nn0][n2][n3] = tmp_array1[2*MN0 + 2*MN2 + 2*n3  ];
	      ImV1[nn0][n2][n3] = tmp_array1[2*MN0 + 2*MN2 + 2*n3+1];
	    }
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

        case 4:

          for (n3=0; n3<Ngrid3; n3++){
            MN = MN1 + MN2 + n3;
            ReV1[nn0][n2][n3] = Density_Grid[2][MN];
            ImV1[nn0][n2][n3] = Density_Grid[3][MN];
          }     

        break;

        } /* switch(den_flag) */
      }    
    }
  }

  /****************************************************
                       FFT of Dens
  ****************************************************/

  /*
  n2 = Ngrid2/2;
  n3 = Ngrid3/2;

  for (n1=0; n1<Ngrid1; n1++){
    printf("%2d %15.12f %15.12f\n",n1,ReV1[n1][n2][n3],ImV2[n1][n2][n3]);
  }
  */

  FFT_Poisson(1, 2, Ngrid1, Ngrid2, Ngrid3,
              My_NGrid1_Poisson, My_NGrid2_Poisson,
              -1, -1, -1, ReV1, ImV1, ReV2, ImV2);
}






void Get_Value_inReal(int complex_flag,
                      double ***ReV2, double ***ImV2,
                      double ***ReV1, double ***ImV1,
                      double *Value_Grid, double *iValue_Grid)
{
  int ct_AN,n1,n2,n3,k1,k2,k3,kk2,N3[4];
  int N,i,j;
  int nn0,nn1,nnn1,MN,MN0,MN1,MN2;
  double time0;
  int h_AN,Gh_AN,spin,spinmax;
  int Nc,GNc,GN,Nh,Nog,factor;
  double tmp0,tmp1,sk1,sk2,sk3,tot_den;
  double x,y,z,Gx,Gy,Gz,Eden[2],DenA,G2;
  double ReTmp,ImTmp,c_coe,p_coe;
  double *tmp_array0;
  double *tmp_array1;
  double TStime,TEtime;
  int numprocs,myid,tag=999,ID,IDS,IDR;

  MPI_Status stat;
  MPI_Request request;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if      (complex_flag==0) factor = 1;
  else if (complex_flag==1) factor = 2; 

  /****************************************************
                  Inverse FFT of Rhok
  ****************************************************/

  FFT_Poisson(2, 1, Ngrid2, Ngrid1, Ngrid3,
              My_NGrid2_Poisson, My_NGrid1_Poisson,
              1, 1, 1, ReV2, ImV2, ReV1, ImV1);

  /****************************************************
                   ReV1 -> Density_Grid
  ****************************************************/

  /* initialize */
  if   (complex_flag==0){
    for (MN=0; MN<My_NumGrid1; MN++){
      Value_Grid[MN] = 0.0;
    }
  }
  else{
    for (MN=0; MN<My_NumGrid1; MN++){
      Value_Grid[MN] = 0.0;
     iValue_Grid[MN] = 0.0;
    }
  } 

  /* use their potential using MPI */
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    /* Isend */
    if (Num_ISnd_Grid1[IDS]!=0){
      tmp_array0 = (double*)malloc(sizeof(double)*Num_ISnd_Grid1[IDS]*Ngrid2*Ngrid3*factor);

      for (i=0; i<Num_ISnd_Grid1[IDS]; i++){ 
	n1 = ISnd_Grid1[IDS][i];
	nn1 = n1 - Start_Grid1[myid];
	MN1 = i*Ngrid2*Ngrid3;
	for (n2=0; n2<Ngrid2; n2++){
	  MN2 = n2*Ngrid3;

	  if   (complex_flag==0){
	    for (n3=0; n3<Ngrid3; n3++){
	      MN = MN1 + MN2 + n3;
	      tmp_array0[MN] = ReV1[nn1][n2][n3];
	    }
	  }
	  else{
	    for (n3=0; n3<Ngrid3; n3++){
	      tmp_array0[2*MN1 + 2*MN2 + 2*n3  ] = ReV1[nn1][n2][n3];
	      tmp_array0[2*MN1 + 2*MN2 + 2*n3+1] = ImV1[nn1][n2][n3];
	    }
	  }
	}
      } 

      tag = 999;
      MPI_Isend(&tmp_array0[0], Num_ISnd_Grid1[IDS]*Ngrid2*Ngrid3*factor,
		MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);
    }

    /* Recv */
    if (Num_IRcv_Grid1[IDR]!=0){

      tmp_array1 = (double*)malloc(sizeof(double)*Num_IRcv_Grid1[IDR]*Ngrid2*Ngrid3*factor); 

      tag = 999;
      MPI_Recv(&tmp_array1[0], Num_IRcv_Grid1[IDR]*Ngrid2*Ngrid3*factor,
	       MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

      for (i=0; i<Num_IRcv_Grid1[IDR]; i++){ 
	n1 = IRcv_Grid1[IDR][i];
	nn1 = My_Cell0[n1];
	MN1 = nn1*Ngrid2*Ngrid3;
	MN0 = i*Ngrid2*Ngrid3;

	for (n2=0; n2<Ngrid2; n2++){
	  MN2 = n2*Ngrid3;

	  if   (complex_flag==0){
	    for (n3=0; n3<Ngrid3; n3++){
	      Value_Grid[MN1 + MN2 + n3] = tmp_array1[MN0 + MN2 + n3];
	    }
	  }
	  else{
	    for (n3=0; n3<Ngrid3; n3++){
	      Value_Grid[MN1 + MN2 + n3] = tmp_array1[2*MN0 + 2*MN2 + 2*n3  ];
	     iValue_Grid[MN1 + MN2 + n3] = tmp_array1[2*MN0 + 2*MN2 + 2*n3+1];
	    }
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

        if   (complex_flag==0){
  	  for (n3=0; n3<Ngrid3; n3++){
	    MN = MN1 + MN2 + n3;
	    Value_Grid[MN] = ReV1[n1][n2][n3];
	  }
	}
        else{
  	  for (n3=0; n3<Ngrid3; n3++){
	    MN = MN1 + MN2 + n3;
	    Value_Grid[MN] = ReV1[n1][n2][n3];
	   iValue_Grid[MN] = ImV1[n1][n2][n3];
	  }
        }
      }
    }
  }

  /*
  for (n1=0; n1<Ngrid1; n1++){
    printf("%2d %18.15f %18.15f\n",n1,ReV1[n1][Ngrid2/2][Ngrid3/2],ImV1[n1][Ngrid2/2][Ngrid3/2]);
  }
  */

}
