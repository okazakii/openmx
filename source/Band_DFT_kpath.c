/**********************************************************************
  Band_DFT_kpath.c:

     Band_DFT_kpath.c is a subroutine to calculate eigenenergies
     at given k-points for writing the band dispersion.

     input: information of k points
     output: None
             file.Band

  Log of Band_DFT_kpath.c:

     12/May/2003  Released by H.Kino
     25/Dec/2003  a non-collinear part (added by T.Ozaki)

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "openmx_common.h"

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

static void Band_DFT_kpath_Col(
                      int nkpath, int *n_perk,
                      double ***kpath, char ***kname, 
                      int  SpinP_switch, 
                      double *****nh,
                      double *****ImNL,
                      double ****CntOLP);

static void Band_DFT_kpath_NonCol(
                      int nkpath, int *n_perk,
                      double ***kpath, char ***kname, 
                      int  SpinP_switch, 
                      double *****nh,
                      double *****ImNL,
                      double ****CntOLP);

void Band_DFT_kpath( int nkpath, int *n_perk,
                     double ***kpath, char ***kname, 
                     int  SpinP_switch, 
                     double *****nh,
                     double *****ImNL,
                     double ****CntOLP)
{

  if (SpinP_switch==0 || SpinP_switch==1){
    Band_DFT_kpath_Col(nkpath, n_perk, kpath, kname, SpinP_switch, nh, ImNL, CntOLP);
  }
  else if (SpinP_switch==3){
    Band_DFT_kpath_NonCol(nkpath, n_perk, kpath, kname, SpinP_switch, nh, ImNL, CntOLP);
  }

}


void Band_DFT_kpath_Col( int nkpath, int *n_perk,
                         double ***kpath, char ***kname, 
                         int  SpinP_switch, 
                         double *****nh,
                         double *****ImNL,
                         double ****CntOLP)
{
  static int firsttime=1;
  int i,j,k,l,n,wan,SOMO,HOMO,ITZ;
  int *MP,*arpo,num_kloop0,T_knum;
  int i1,j1,po,spin,spinsize,n1;
  int num2,RnB,l1,l2,l3,S_knum,E_knum;
  int ct_AN,h_AN,wanA,tnoA,wanB,tnoB;
  int GA_AN,Anum,ID,ID1,s1,e1,kloop,kloop0;
  double time0;
  int LB_AN,GB_AN,Bnum;
  double snum_i,snum_j,snum_k,k1,k2,k3,sum,sumi,Num_State,FermiF;
  double x,Dnum,Dnum2,AcP,ChemP_MAX,ChemP_MIN,EV_cut0;
  double **ko,*M1;
  double *koS;
  double ****Dummy_ImNL;
  dcomplex **H,**S,**C,**TmpM;
  dcomplex Ctmp1,Ctmp2;
  float ****EigenVal;
  int ii,ij,ik;
  double u2,v2,uv,vu,tmp;
  double TZ,Eele_Zero,dum,sumE,kRn,si,co;
  double Resum,ResumE,Redum,Redum2,Imdum;
  double TStime,TEtime, SiloopTime,EiloopTime;
  char file_EV[YOUSO10];
  FILE *fp_EV;

  double FermiEps = 1.0e-14;
  double x_cut = 30.0;
  int numprocs,myid;
  int *ik_list,*i_perk_list;
  int i_perk,id;

  char file_Band[YOUSO10];
  FILE *fp_Band;
  char buf[fp_bsize];          /* setvbuf */

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  /*d*/
  if (myid==Host_ID && 0<level_stdout){
    printf("Band_DFT_kpath start\n");fflush(stdout);
  }  

  dtime(&TStime);

  /****************************************************
   Allocation

   int       MP[List_YOUSO[1]]
   double    ko[List_YOUSO[23]][n+1]
   double    koS[n+1];
   dcomplex  H[n+1][n+1]
   dcomplex  S[n+1][n+1]
   double    M1[n+1]
   dcomplex  C[n+1][n+1]
  ****************************************************/

  MP = (int*)malloc(sizeof(int)*List_YOUSO[1]);
  
  n = 0;
  for (i=1; i<=atomnum; i++){
    wanA = WhatSpecies[i];
    n += Spe_Total_CNO[wanA];
  }

  ko = (double**)malloc(sizeof(double*)*List_YOUSO[23]);
  for (i=0; i<List_YOUSO[23]; i++){
    ko[i] = (double*)malloc(sizeof(double)*(n+1));
  }

  koS = (double*)malloc(sizeof(double)*(n+1));

  if (firsttime) {
  PrintMemory("Band_DFT: ko",sizeof(double)*List_YOUSO[23]*(n+1),NULL);
  }

  H = (dcomplex**)malloc(sizeof(dcomplex*)*(n+1));
  for (j=0; j<n+1; j++){
    H[j] = (dcomplex*)malloc(sizeof(dcomplex)*(n+1));
  }

  TmpM = (dcomplex**)malloc(sizeof(dcomplex*)*(n+1));
  for (j=0; j<n+1; j++){
    TmpM[j] = (dcomplex*)malloc(sizeof(dcomplex)*(n+1));
  }

  S = (dcomplex**)malloc(sizeof(dcomplex*)*(n+1));
  for (i=0; i<n+1; i++){
    S[i] = (dcomplex*)malloc(sizeof(dcomplex)*(n+1));
  }

  M1 = (double*)malloc(sizeof(double)*(n+1));

  C = (dcomplex**)malloc(sizeof(dcomplex*)*(n+1));
  for (j=0; j<n+1; j++){
    C[j] = (dcomplex*)malloc(sizeof(dcomplex)*(n+1));
  }

  if      (SpinP_switch==0){ spinsize=1; }
  else if (SpinP_switch==1){ spinsize=2; }
  else if (SpinP_switch==3){ spinsize=1; } 

  EigenVal = (float****)malloc(sizeof(float***)*(nkpath+1));
  for (ik=0; ik<=nkpath; ik++) {
    EigenVal[ik] = (float***)malloc(sizeof(float**)*(n_perk[ik]+1));
    for (i_perk=0; i_perk<=n_perk[ik]; i_perk++) {
      EigenVal[ik][i_perk] = (float**)malloc(sizeof(float*)*spinsize);
      for (spin=0; spin<spinsize; spin++){
        EigenVal[ik][i_perk][spin] = (float*)malloc(sizeof(float)*(n+2));
      }
    }
  }

  /* no spin-orbit coupling */
  if (SO_switch==0){
    Dummy_ImNL = (double****)malloc(sizeof(double***)*1);
    Dummy_ImNL[0] = (double***)malloc(sizeof(double**)*1);
    Dummy_ImNL[0][0] = (double**)malloc(sizeof(double*)*1);
    Dummy_ImNL[0][0][0] = (double*)malloc(sizeof(double)*1);
  }

  if (firsttime)
  PrintMemory("Band_DFT: D",
        sizeof(double)*List_YOUSO[23]*List_YOUSO[27]*
                       List_YOUSO[28]*List_YOUSO[29]*(n+1),NULL);

  /* for PrintMemory */
  firsttime=0;

  TZ = 0.0;
  for (i=1; i<=atomnum; i++){
    wan = WhatSpecies[i];
    TZ = TZ + Spe_Core_Charge[wan];
  }
  ITZ = TZ;
  SOMO = ITZ % 2;
  HOMO = (TZ-SOMO)/2;

  dtime(&SiloopTime);
  Eele_Zero = 0.0;

  if (myid==Host_ID && 0<level_stdout){
    printf("kpath\n");
    for (ik=1;ik<=nkpath; ik++) {
      printf("%d (%lf %lf %lf)->(%lf %lf %lf)\n", ik,
              kpath[ik][1][1],kpath[ik][1][2],kpath[ik][1][3],
              kpath[ik][2][1],kpath[ik][2][2],kpath[ik][2][3]);fflush(stdout);
    }
  }

  /*****************************************************
        distribute node for the parallel calculation
  *****************************************************/

  T_knum = 0;
  for (ik=1; ik<=nkpath; ik++) {
    T_knum += n_perk[ik];
  }

  ik_list     = (int*)malloc(sizeof(int)*(T_knum+3));
  i_perk_list = (int*)malloc(sizeof(int)*(T_knum+3));
  arpo = (int*)malloc(sizeof(int)*numprocs);
  
  T_knum = 0;
  for (ik=1; ik<=nkpath; ik++) {
    for (i_perk=1; i_perk<=n_perk[ik]; i_perk++) {
      ik_list[T_knum] = ik;
      i_perk_list[T_knum] = i_perk; 
      T_knum++;  
    }
  }

  /* allocate k-points into proccessors */

  if (T_knum<=myid){
    S_knum = -10;
    E_knum = -100;
    num_kloop0 = 1;
  }
  else if (T_knum<numprocs) {
    S_knum = myid;
    E_knum = myid;
    num_kloop0 = 1;
  }
  else {
    tmp = (double)T_knum/(double)numprocs; 
    num_kloop0 = (int)tmp + 1;
    S_knum = (int)((double)myid*(tmp+0.0001)); 
    E_knum = (int)((double)(myid+1)*(tmp+0.0001)) - 1;
    if (myid==(numprocs-1)) E_knum = T_knum - 1;
    if (E_knum<0)           E_knum = 0;
  }

  /*****************************************************
           calculate eigenvalues along k-paths
  *****************************************************/

  for (kloop0=0; kloop0<num_kloop0; kloop0++){

    kloop = kloop0 + S_knum;
    arpo[myid] = -1;
    if (S_knum<=kloop && kloop<=E_knum) arpo[myid] = kloop;
    for (ID=0; ID<numprocs; ID++){
      MPI_Bcast(&arpo[ID], 1, MPI_INT, ID, mpi_comm_level1);
    }

    /*
    if (myid==Host_ID){
      printf("%d (%lf %lf %lf)->(%lf %lf %lf)\n", ik,
              kpath[ik][1][1],kpath[ik][1][2],kpath[ik][1][3],
              kpath[ik][2][1],kpath[ik][2][2],kpath[ik][2][3]);fflush(stdout);
    }
    */

    /* set S and diagonalize it */
    
    for (ID=0; ID<numprocs; ID++){

      kloop = arpo[ID];

      if (0<=kloop){

        ik     = ik_list[kloop];
        i_perk = i_perk_list[kloop];

        id=1;
        k1 = kpath[ik][1][id]+
          (kpath[ik][2][id]-kpath[ik][1][id])*(i_perk-1)/(n_perk[ik]-1) + Shift_K_Point;
        id=2;
        k2 = kpath[ik][1][id]+
          (kpath[ik][2][id]-kpath[ik][1][id])*(i_perk-1)/(n_perk[ik]-1) - Shift_K_Point;
        id=3;
        k3 = kpath[ik][1][id]+
          (kpath[ik][2][id]-kpath[ik][1][id])*(i_perk-1)/(n_perk[ik]-1) + 2.0*Shift_K_Point;

        Overlap_Band(ID,CntOLP,TmpM,MP,k1,k2,k3);
        n = TmpM[0][0].r;

        if (myid==ID){
	  for (i1=1; i1<=n; i1++){
	    for (j1=1; j1<=n; j1++){
	      S[i1][j1].r = TmpM[i1][j1].r;
	      S[i1][j1].i = TmpM[i1][j1].i;
	    } 
	  } 
        } 
      }
    }

    kloop = arpo[myid];

    if (0<=kloop) {

      EigenBand_lapack(S,ko[0],n,1);

      if (3<=level_stdout){
	printf(" myid=%2d kloop %2d  k1 k2 k3 %10.6f %10.6f %10.6f\n",
 	       myid,kloop,k1,k2,k3);
	for (i1=1; i1<=n; i1++){
	  printf("  Eigenvalues of OLP  %2d  %15.12f\n",i1,ko[0][i1]);
	}
      }

      /* minus eigenvalues to 1.0e-14 */

      for (l=1; l<=n; l++){
        if (ko[0][l]<0.0) ko[0][l] = 1.0e-14;
        koS[l] = ko[0][l];
      }

      /* calculate S*1/sqrt(ko) */

      for (l=1; l<=n; l++) M1[l] = 1.0/sqrt(ko[0][l]);

      /* S * M1 */

      for (i1=1; i1<=n; i1++){
	for (j1=1; j1<=n; j1++){
	  S[i1][j1].r = S[i1][j1].r*M1[j1];
	  S[i1][j1].i = S[i1][j1].i*M1[j1];
	} 
      } 

    } /* if (0<=kloop) */
    
    /* construct H' and diagonalize it */

    for (spin=0; spin<=SpinP_switch; spin++){

      for (ID=0; ID<numprocs; ID++){

        kloop = arpo[ID];

        if (0<=kloop) {

	  ik     = ik_list[kloop];
	  i_perk = i_perk_list[kloop];

	  id=1;
	  k1 = kpath[ik][1][id]+
	    (kpath[ik][2][id]-kpath[ik][1][id])*(i_perk-1)/(n_perk[ik]-1) + Shift_K_Point;
	  id=2;
	  k2 = kpath[ik][1][id]+
	    (kpath[ik][2][id]-kpath[ik][1][id])*(i_perk-1)/(n_perk[ik]-1) - Shift_K_Point;
	  id=3;
	  k3 = kpath[ik][1][id]+
	    (kpath[ik][2][id]-kpath[ik][1][id])*(i_perk-1)/(n_perk[ik]-1) + 2.0*Shift_K_Point;

          Hamiltonian_Band(ID, nh[spin], TmpM, MP, k1, k2, k3);

          if (myid==ID){
	    for (i1=1; i1<=n; i1++){
	      for (j1=1; j1<=n; j1++){
	        H[i1][j1].r = TmpM[i1][j1].r;
	        H[i1][j1].i = TmpM[i1][j1].i;
	      } 
	    } 
          } 
	} /* if (0<=kloop) */
      } /* ID */

      kloop = arpo[myid];

      if (0<=kloop){

        ik     = ik_list[kloop];
	i_perk = i_perk_list[kloop];

        /****************************************************
 	                 M1 * U^t * H * U * M1
        ****************************************************/

        /* transpose S */

        if (spin==0){
	  for (i1=1; i1<=n; i1++){
	    for (j1=i1+1; j1<=n; j1++){
	      Ctmp1 = S[i1][j1];
	      Ctmp2 = S[j1][i1];
	      S[i1][j1] = Ctmp2;
	      S[j1][i1] = Ctmp1;
	    }
	  }
	}

        /* H * U * M1 */

        for (j1=1; j1<=n; j1++){
 	  for (i1=1; i1<=n; i1++){

	    sum  = 0.0;
	    sumi = 0.0;

	    for (l=1; l<=n; l++){
	      sum  += H[i1][l].r*S[j1][l].r - H[i1][l].i*S[j1][l].i;
	      sumi += H[i1][l].r*S[j1][l].i + H[i1][l].i*S[j1][l].r;
	    }

	    C[j1][i1].r = sum;
	    C[j1][i1].i = sumi;
	  }
	}     

        /* M1 * U^+ H * U * M1 */

        for (i1=1; i1<=n; i1++){
          for (j1=1; j1<=n; j1++){
	    sum  = 0.0;
            sumi = 0.0;
	    for (l=1; l<=n; l++){
	      sum  +=  S[i1][l].r*C[j1][l].r + S[i1][l].i*C[j1][l].i;
	      sumi +=  S[i1][l].r*C[j1][l].i - S[i1][l].i*C[j1][l].r;
	    }
	    H[i1][j1].r = sum;
	    H[i1][j1].i = sumi;
	  }
	} 

        /* H to C */

	for (i1=1; i1<=n; i1++){
	  for (j1=1; j1<=n; j1++){
	    C[i1][j1] = H[i1][j1];
	  }
	}

        /* penalty for ill-conditioning states */

	EV_cut0 = Threshold_OLP_Eigen;

	for (i1=1; i1<=n; i1++){

	  if (koS[i1]<EV_cut0){
	    C[i1][i1].r += pow((koS[i1]/EV_cut0),-2.0) - 1.0;
	  }
 
	  /* cutoff the interaction between the ill-conditioned state */
 
	  if (1.0e+3<C[i1][i1].r){
	    for (j1=1; j1<=n; j1++){
	      C[i1][j1] = Complex(0.0,0.0);
	      C[j1][i1] = Complex(0.0,0.0);
	    }
	    C[i1][i1].r = 1.0e+4;
	  }
	}

        /* solve eigenvalue problem */

	n1 = n;
        EigenBand_lapack(C,ko[spin],n1,0);

	for (l=1; l<=n1; l++){
	  EigenVal[ik][i_perk][spin][l-1] = (float)ko[spin][l];
	}

      } /* if (0<=kloop) */
    } /* spin */
  } /* kloop0 */

  /****************************************************
     MPI:

     EigenVal
  ****************************************************/

  tmp = (double)T_knum/(double)numprocs; 

  for (spin=0; spin<=SpinP_switch; spin++){
    for (kloop=0; kloop<T_knum; kloop++){

      ik     = ik_list[kloop];
      i_perk = i_perk_list[kloop];

      for (ID=0; ID<numprocs; ID++){

	if (T_knum<=ID){
	  s1 = -10;
	  e1 = -100;
	}
	else if (T_knum<numprocs) {
	  s1 = ID;
	  e1 = ID;
	}
	else {
	  s1 = (int)((double)ID*(tmp+0.0001));
          e1 = (int)((double)(ID+1)*(tmp+0.0001)) - 1;
          if (ID==(numprocs-1)) e1 = T_knum - 1;
          if (e1<0)             e1 = 0;
	}
        if (s1<=kloop && kloop<=e1)  ID1 = ID;
      }

      MPI_Bcast(&EigenVal[ik][i_perk][spin][0], n+1, MPI_FLOAT, ID1, mpi_comm_level1);

    } /* kloop */
  } /* spin */

  /****************************************************
                    write a file 
  ****************************************************/

  if (myid==Host_ID) {

    strcpy(file_Band,".Band");
    fnjoint(filepath,filename,file_Band);  

    if ((fp_Band = fopen(file_Band,"w"))==NULL) {

#ifdef xt3
      setvbuf(fp_Band,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      printf("<Band_DFT_kpath> can not open a file (%s)\n",file_Band);
      return;
    }

    fprintf(fp_Band," %d  %d  %lf\n",n,SpinP_switch,ChemP);
    for (i=1;i<=3;i++) 
      for (j=1;j<=3;j++) {
	fprintf(fp_Band,"%lf ", rtv[i][j]);
      }
    fprintf(fp_Band,"\n");
    fprintf(fp_Band,"%d\n",nkpath);
    for (i=1;i<=nkpath;i++) {
      fprintf(fp_Band,"%d %lf %lf %lf  %lf %lf %lf  %s %s\n",
	      n_perk[i],
	      kpath[i][1][1], kpath[i][1][2], kpath[i][1][3],
	      kpath[i][2][1], kpath[i][2][2], kpath[i][2][3],
	      kname[i][1],kname[i][2]);
    }

    for (ik=1; ik<=nkpath; ik++) {
      for (i_perk=1; i_perk<=n_perk[ik]; i_perk++) {

        id=1;
        k1 = kpath[ik][1][id]+
          (kpath[ik][2][id]-kpath[ik][1][id])*(i_perk-1)/(n_perk[ik]-1) + Shift_K_Point;
        id=2;
        k2 = kpath[ik][1][id]+
          (kpath[ik][2][id]-kpath[ik][1][id])*(i_perk-1)/(n_perk[ik]-1) - Shift_K_Point;
        id=3;
        k3 = kpath[ik][1][id]+
          (kpath[ik][2][id]-kpath[ik][1][id])*(i_perk-1)/(n_perk[ik]-1) + 2.0*Shift_K_Point;

        for (spin=0; spin<=SpinP_switch; spin++){

	  fprintf(fp_Band,"%d %lf %lf %lf\n", n,k1,k2,k3);

	  for (l=0; l<n; l++) {
	    fprintf(fp_Band,"%lf ",EigenVal[ik][i_perk][spin][l]);
	  }
	  fprintf(fp_Band  ,"\n");
	}
      }
    }

    fclose(fp_Band);

  } /* if (myid==Host_ID) */
  
  /****************************************************
                       free arrays
  ****************************************************/

  free(MP);

  for (i=0; i<List_YOUSO[23]; i++){
    free(ko[i]);
  }
  free(ko);

  free(koS);

  for (j=0; j<n+1; j++){
    free(H[j]);
  }
  free(H);

  for (j=0; j<n+1; j++){
    free(TmpM[j]);
  }
  free(TmpM);

  for (i=0; i<n+1; i++){
    free(S[i]);
  }
  free(S);

  free(M1);

  for (j=0; j<n+1; j++){
    free(C[j]);
  }
  free(C);

  for (ik=0; ik<=nkpath; ik++) {
    for (i_perk=0; i_perk<=n_perk[ik]; i_perk++) {
      for (spin=0; spin<spinsize; spin++){
        free(EigenVal[ik][i_perk][spin]);
      }
      free(EigenVal[ik][i_perk]);
    }
    free(EigenVal[ik]);
  }
  free(EigenVal);

  /* no spin-orbit coupling */
  if (SO_switch==0){
    free(Dummy_ImNL[0][0][0]);
    free(Dummy_ImNL[0][0]);
    free(Dummy_ImNL[0]);
    free(Dummy_ImNL);
  }

  free(ik_list);
  free(i_perk_list);
  free(arpo);

  /* elapsed time */
  dtime(&TEtime);
}







void Band_DFT_kpath_NonCol( int nkpath, int *n_perk,
                            double ***kpath, char ***kname, 
                            int  SpinP_switch, 
                            double *****nh,
                            double *****ImNL,
                            double ****CntOLP)
{
  int i,j,k,l,n,wan,SOMO,HOMO,ITZ;
  int *MP;
  int i1,j1,po,spin,n1;
  int num2,RnB,l1,l2,l3,ii1,jj1,m,n2;
  int ct_AN,h_AN,wanA,tnoA,wanB,tnoB;
  int GA_AN,Anum;
  double time0;
  int LB_AN,GB_AN,Bnum;
  double EV_cut0;
  double snum_i,snum_j,snum_k,k1,k2,k3,sum,sumi,Num_State,FermiF;
  double x,Dnum,Dnum2,AcP,ChemP_MAX,ChemP_MIN;
  double *ko,*M1;
  double *koS;
  double *****Dummy_ImNL;
  dcomplex **H,**S,**C,**TmpM;
  dcomplex Ctmp1,Ctmp2;
  float ***EigenVal;
  int ii,ij,ik;
  int *ik_list,*i_perk_list,*arpo;
  double u2,v2,uv,vu,tmp;
  double TZ,Eele_Zero,dum,sumE,kRn,si,co;
  double Resum,ResumE,Redum,Redum2,Imdum;
  double TStime,TEtime, SiloopTime,EiloopTime;
  char file_EV[YOUSO10];
  FILE *fp_EV;

  double FermiEps = 1.0e-14;
  double x_cut = 30.0;
  int numprocs,myid,ID,ID1;
  int e1,s1,kloop,T_knum,mn,S_knum,E_knum;
  int num_kloop0,kloop0,spinsize; 

  int i_perk,id;
  char buf[fp_bsize];          /* setvbuf */
  char file_Band[YOUSO10];
  FILE *fp_Band;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  /*d*/

  if (myid==Host_ID && 0<level_stdout){
    printf("Band_DFT_kpath start\n");fflush(stdout);
  }

  dtime(&TStime);

  /****************************************************
             calculation of the array size
  ****************************************************/

  n = 0;
  for (i=1; i<=atomnum; i++){
    wanA  = WhatSpecies[i];
    n  = n + Spe_Total_CNO[wanA];
  }
  n2 = 2*n + 2;

  /****************************************************
   Allocation

   int       MP[List_YOUSO[1]]
   double    ko[n2]
   double    koS[n+1];
   dcomplex  H[n2][n2]
   dcomplex  S[n2][n2]
   double    M1[n2]
   dcomplex  C[n2][n2]
  ****************************************************/

  MP = (int*)malloc(sizeof(int)*List_YOUSO[1]);

  ko = (double*)malloc(sizeof(double)*n2);
  koS = (double*)malloc(sizeof(double)*(n+1));

  H = (dcomplex**)malloc(sizeof(dcomplex*)*n2);
  for (j=0; j<n2; j++){
    H[j] = (dcomplex*)malloc(sizeof(dcomplex)*n2);
  }

  TmpM = (dcomplex**)malloc(sizeof(dcomplex*)*n2);
  for (j=0; j<n2; j++){
    TmpM[j] = (dcomplex*)malloc(sizeof(dcomplex)*n2);
  }

  S = (dcomplex**)malloc(sizeof(dcomplex*)*n2);
  for (i=0; i<n2; i++){
    S[i] = (dcomplex*)malloc(sizeof(dcomplex)*n2);
  }

  M1 = (double*)malloc(sizeof(double)*n2);

  C = (dcomplex**)malloc(sizeof(dcomplex*)*n2);
  for (j=0; j<n2; j++){
    C[j] = (dcomplex*)malloc(sizeof(dcomplex)*n2);
  }

  if      (SpinP_switch==0){ spinsize=1; }
  else if (SpinP_switch==1){ spinsize=2; }
  else if (SpinP_switch==3){ spinsize=1; } 

  /*  
  printf("myid=%2d nkpath+1=%2d\n",myid,nkpath+1);
  for (ik=0; ik<=nkpath; ik++) {
    printf("myid=%2d ik=%2d n_perk=%2d\n",myid,ik,n_perk[ik]+1); 
  }  
  printf("myid=%2d n2=%2d\n",myid,n2);
  */

  EigenVal = (float***)malloc(sizeof(float**)*(nkpath+1));
  for (ik=0; ik<=nkpath; ik++) {
    EigenVal[ik] = (float**)malloc(sizeof(float*)*(n_perk[ik]+1));
    for (i_perk=0; i_perk<=n_perk[ik]; i_perk++) {
      EigenVal[ik][i_perk] = (float*)malloc(sizeof(float)*n2);
    }
  }

  /* non-spin-orbit coupling and non-LDA+U */
  if (SO_switch==0 && Hub_U_switch==0 && Constraint_NCS_switch==0 
      && Zeeman_NCS_switch==0 && Zeeman_NCO_switch==0){
    Dummy_ImNL = (double*****)malloc(sizeof(double****)*1);
    Dummy_ImNL[0] = (double****)malloc(sizeof(double***)*1);
    Dummy_ImNL[0][0] = (double***)malloc(sizeof(double**)*1);
    Dummy_ImNL[0][0][0] = (double**)malloc(sizeof(double*)*1);
    Dummy_ImNL[0][0][0][0] = (double*)malloc(sizeof(double)*1);
  }

  /*****************************************************
    TZ, ITZ, and HOMO 
  *****************************************************/

  TZ = 0.0;
  for (i=1; i<=atomnum; i++){
    wan = WhatSpecies[i];
    TZ = TZ + Spe_Core_Charge[wan];
  }
  ITZ = TZ;
  HOMO = TZ;

  dtime(&SiloopTime);
  Eele_Zero = 0.0;

  if (myid==Host_ID && 0<level_stdout){
    printf("kpath\n");
    for (ik=1; ik<=nkpath; ik++) {
      printf("%d (%lf %lf %lf)->(%lf %lf %lf)\n", ik,
              kpath[ik][1][1],kpath[ik][1][2],kpath[ik][1][3],
              kpath[ik][2][1],kpath[ik][2][2],kpath[ik][2][3]);fflush(stdout);
    }
  }

  /*****************************************************
        distribute node for the parallel calculation
  *****************************************************/

  T_knum = 0;
  for (ik=1; ik<=nkpath; ik++) {
    T_knum += n_perk[ik];
  }

  ik_list     = (int*)malloc(sizeof(int)*(T_knum+3));
  i_perk_list = (int*)malloc(sizeof(int)*(T_knum+3));
  arpo = (int*)malloc(sizeof(int)*numprocs);
  
  T_knum = 0;
  for (ik=1; ik<=nkpath; ik++) {
    for (i_perk=1; i_perk<=n_perk[ik]; i_perk++) {
      ik_list[T_knum] = ik;
      i_perk_list[T_knum] = i_perk; 
      T_knum++;  
    }
  }

  /* allocate k-points into proccessors */

  if (T_knum<=myid){
    S_knum = -10;
    E_knum = -100;
    num_kloop0 = 1;
  }
  else if (T_knum<numprocs) {
    S_knum = myid;
    E_knum = myid;
    num_kloop0 = 1;
  }
  else {
    tmp = (double)T_knum/(double)numprocs; 
    num_kloop0 = (int)tmp + 1;
    S_knum = (int)((double)myid*(tmp+0.0001)); 
    E_knum = (int)((double)(myid+1)*(tmp+0.0001)) - 1;
    if (myid==(numprocs-1)) E_knum = T_knum - 1;
    if (E_knum<0)           E_knum = 0;
  }

  /*****************************************************
           calculate eigenvalues along k-paths
  *****************************************************/

  for (kloop0=0; kloop0<num_kloop0; kloop0++){

    kloop = kloop0 + S_knum;
    arpo[myid] = -1;
    if (S_knum<=kloop && kloop<=E_knum) arpo[myid] = kloop;
    for (ID=0; ID<numprocs; ID++){
      MPI_Bcast(&arpo[ID], 1, MPI_INT, ID, mpi_comm_level1);
    }

    /* set S and diagonalize it */
    
    for (ID=0; ID<numprocs; ID++){

      kloop = arpo[ID];

      if (0<=kloop){

        ik     = ik_list[kloop];
        i_perk = i_perk_list[kloop];

        id=1;
        k1 = kpath[ik][1][id]+
          (kpath[ik][2][id]-kpath[ik][1][id])*(i_perk-1)/(n_perk[ik]-1) + Shift_K_Point;
        id=2;
        k2 = kpath[ik][1][id]+
          (kpath[ik][2][id]-kpath[ik][1][id])*(i_perk-1)/(n_perk[ik]-1) - Shift_K_Point;
        id=3;
        k3 = kpath[ik][1][id]+
          (kpath[ik][2][id]-kpath[ik][1][id])*(i_perk-1)/(n_perk[ik]-1) + 2.0*Shift_K_Point;

        Overlap_Band(ID,CntOLP,TmpM,MP,k1,k2,k3);
        n = TmpM[0][0].r;

        if (myid==ID){
	  for (i1=1; i1<=n; i1++){
	    for (j1=1; j1<=n; j1++){
	      S[i1][j1].r = TmpM[i1][j1].r;
	      S[i1][j1].i = TmpM[i1][j1].i;
	    } 
	  } 
        } 
      }
    } /* ID */

    kloop = arpo[myid];

    if (0<=kloop){

      EigenBand_lapack(S,ko,n,1);

      /* minus eigenvalues to 1.0e-14 */
      for (l=1; l<=n; l++){
        if (ko[l]<0.0){
          ko[l] = 1.0e-14;

          if (2<=level_stdout){
  	    printf("found an eigenvalue smaller than %10.8f of OLP kloop=%2d\n",
                           Threshold_OLP_Eigen,kloop);
	  }
	}

        koS[l] = ko[l];
      }

      /* calculate S*1/sqrt(ko) */

      for (l=1; l<=n; l++) M1[l] = 1.0/sqrt(ko[l]);

      /* S * M1  */

      for (i1=1; i1<=n; i1++){
	for (j1=1; j1<=n; j1++){
	  S[i1][j1].r = S[i1][j1].r*M1[j1];
	  S[i1][j1].i = S[i1][j1].i*M1[j1];
	} 
      } 

    } /* if (0<=kloop) */

    /* set H' and diagonalize it */

    for (ID=0; ID<numprocs; ID++){

      kloop = arpo[ID];

      if (0<=kloop){

	ik     = ik_list[kloop];
	i_perk = i_perk_list[kloop];

	id=1;
	k1 = kpath[ik][1][id]+
	  (kpath[ik][2][id]-kpath[ik][1][id])*(i_perk-1)/(n_perk[ik]-1) + Shift_K_Point;
	id=2;
	k2 = kpath[ik][1][id]+
	  (kpath[ik][2][id]-kpath[ik][1][id])*(i_perk-1)/(n_perk[ik]-1) - Shift_K_Point;
	id=3;
	k3 = kpath[ik][1][id]+
	  (kpath[ik][2][id]-kpath[ik][1][id])*(i_perk-1)/(n_perk[ik]-1) + 2.0*Shift_K_Point;

        /****************************************************
                   make a full Hamiltonian matrix
        ****************************************************/

        /* non-spin-orbit coupling and non-LDA+U */
        if (SO_switch==0 && Hub_U_switch==0 && Constraint_NCS_switch==0 
            && Zeeman_NCS_switch==0 && Zeeman_NCO_switch==0){
          Hamiltonian_Band_NC(ID, nh, Dummy_ImNL, TmpM, MP, k1, k2, k3);
	}
        /* spin-orbit coupling or LDA+U */
        else{
          Hamiltonian_Band_NC(ID, nh,       ImNL, TmpM, MP, k1, k2, k3);
	}

	if (myid==ID){
	  for (i1=1; i1<=2*n; i1++){
	    for (j1=1; j1<=2*n; j1++){
	      H[i1][j1].r = TmpM[i1][j1].r;
	      H[i1][j1].i = TmpM[i1][j1].i;
	    } 
	  } 
	} 

      } /* if (0<=kloop) */
    } /* ID */

    kloop = arpo[myid];

    if (0<=kloop){

      ik     = ik_list[kloop];
      i_perk = i_perk_list[kloop];

      /****************************************************
                     M1 * U^+ * H * U * M1
      ****************************************************/

      /* transpose S */

      for (i1=1; i1<=n; i1++){
	for (j1=i1+1; j1<=n; j1++){
	  Ctmp1 = S[i1][j1];
	  Ctmp2 = S[j1][i1];
	  S[i1][j1] = Ctmp2;
	  S[j1][i1] = Ctmp1;
	}
      }

      /* H * U * M1 */

      for (j1=1; j1<=n; j1++){
        for (i1=1; i1<=2*n; i1++){
	  for (m=0; m<=1; m++){

	    sum  = 0.0;
	    sumi = 0.0;
            mn = m*n;
	    for (l=1; l<=n; l++){
	      sum  += H[i1][l+mn].r*S[j1][l].r - H[i1][l+mn].i*S[j1][l].i;
	      sumi += H[i1][l+mn].r*S[j1][l].i + H[i1][l+mn].i*S[j1][l].r;
	    }

	    jj1 = 2*j1 - 1 + m;

	    C[jj1][i1].r = sum;
	    C[jj1][i1].i = sumi;
	  }
	}
      }     

      /* M1 * U^+ H * U * M1 */

      for (i1=1; i1<=n; i1++){

	for (m=0; m<=1; m++){

	  ii1 = 2*i1 - 1 + m;

	  for (j1=1; j1<=2*n; j1++){
	    sum  = 0.0;
	    sumi = 0.0;
            mn = m*n;
	    for (l=1; l<=n; l++){
	      sum  +=  S[i1][l].r*C[j1][l+mn].r + S[i1][l].i*C[j1][l+mn].i;
	      sumi +=  S[i1][l].r*C[j1][l+mn].i - S[i1][l].i*C[j1][l+mn].r;
	    }
	    H[ii1][j1].r = sum;
	    H[ii1][j1].i = sumi;
	  }
	}
      }     

      /* H to C */

      for (i1=1; i1<=2*n; i1++){
	for (j1=1; j1<=2*n; j1++){
	  C[i1][j1].r = H[i1][j1].r;
	  C[i1][j1].i = H[i1][j1].i;
	}
      }

      /* penalty for ill-conditioning states */

      EV_cut0 = Threshold_OLP_Eigen;

      for (i1=1; i1<=n; i1++){

	if (koS[i1]<EV_cut0){
	  C[2*i1-1][2*i1-1].r += pow((koS[i1]/EV_cut0),-2.0) - 1.0;
	  C[2*i1  ][2*i1  ].r += pow((koS[i1]/EV_cut0),-2.0) - 1.0;
	}

	/* cutoff the interaction between the ill-conditioned state */

	if (1.0e+3<C[2*i1-1][2*i1-1].r){
	  for (j1=1; j1<=2*n; j1++){
	    C[2*i1-1][j1    ] = Complex(0.0,0.0);
	    C[j1    ][2*i1-1] = Complex(0.0,0.0);
	    C[2*i1  ][j1    ] = Complex(0.0,0.0);
	    C[j1    ][2*i1  ] = Complex(0.0,0.0);
	  }
	  C[2*i1-1][2*i1-1] = Complex(1.0e+4,0.0);
	  C[2*i1  ][2*i1  ] = Complex(1.0e+4,0.0);
	}
      }

      /* solve eigenvalue problem */

      n1 = 2*n;
      EigenBand_lapack(C, ko, n1,0);

      for (l=1; l<=n1; l++){
        EigenVal[ik][i_perk][l-1] = (float)ko[l];
      }

    } /* if (0<=kloop) */
  } /* kloop0 */

  /****************************************************
     MPI:

     EigenVal
  ****************************************************/

  tmp = (double)T_knum/(double)numprocs; 

  for (kloop=0; kloop<T_knum; kloop++){

    ik     = ik_list[kloop];
    i_perk = i_perk_list[kloop];

    for (ID=0; ID<numprocs; ID++){

      if (T_knum<=ID){
	s1 = -10;
	e1 = -100;
      }
      else if (T_knum<numprocs) {
	s1 = ID;
	e1 = ID;
      }
      else {
	s1 = (int)((double)ID*(tmp+0.0001));
	e1 = (int)((double)(ID+1)*(tmp+0.0001)) - 1;
	if (ID==(numprocs-1)) e1 = T_knum - 1;
	if (e1<0)             e1 = 0;
      }
      if (s1<=kloop && kloop<=e1)  ID1 = ID;
    }

    MPI_Bcast(&EigenVal[ik][i_perk][0], 2*n+1, MPI_FLOAT, ID1, mpi_comm_level1);

  } /* kloop */

  /****************************************************
                    write a file 
  ****************************************************/

  if (myid==Host_ID) {

    strcpy(file_Band,".Band");
    fnjoint(filepath,filename,file_Band);  

    if ((fp_Band = fopen(file_Band,"w"))==NULL) {

#ifdef xt3
      setvbuf(fp_Band,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      printf("<Band_DFT_kpath> can not open a file (%s)\n",file_Band);
      return;
    }

    fprintf(fp_Band," %d  %d  %lf\n",2*n,0,ChemP);  /* set SpinP_switch==0 */
    for (i=1;i<=3;i++) 
      for (j=1;j<=3;j++) {
	fprintf(fp_Band,"%lf ", rtv[i][j]);
      }
    fprintf(fp_Band,"\n");
    fprintf(fp_Band,"%d\n",nkpath);
    for (i=1;i<=nkpath;i++) {
      fprintf(fp_Band,"%d %lf %lf %lf  %lf %lf %lf  %s %s\n",
	      n_perk[i],
	      kpath[i][1][1], kpath[i][1][2], kpath[i][1][3],
	      kpath[i][2][1], kpath[i][2][2], kpath[i][2][3],
	      kname[i][1],kname[i][2]);
    }

    for (ik=1; ik<=nkpath; ik++) {
      for (i_perk=1; i_perk<=n_perk[ik]; i_perk++) {

        id=1;
        k1 = kpath[ik][1][id]+
          (kpath[ik][2][id]-kpath[ik][1][id])*(i_perk-1)/(n_perk[ik]-1) + Shift_K_Point;
        id=2;
        k2 = kpath[ik][1][id]+
          (kpath[ik][2][id]-kpath[ik][1][id])*(i_perk-1)/(n_perk[ik]-1) - Shift_K_Point;
        id=3;
        k3 = kpath[ik][1][id]+
          (kpath[ik][2][id]-kpath[ik][1][id])*(i_perk-1)/(n_perk[ik]-1) + 2.0*Shift_K_Point;

	fprintf(fp_Band,"%d %lf %lf %lf\n", 2*n,k1,k2,k3);

	for (l=0; l<2*n; l++) {
	  fprintf(fp_Band,"%lf ",EigenVal[ik][i_perk][l]);
	}
	fprintf(fp_Band  ,"\n");
      }
    }

    fclose(fp_Band);

  } /* if (myid==Host_ID) */

  /****************************************************
                       free arrays
  ****************************************************/

  free(MP);
  free(ko);
  free(koS);

  for (j=0; j<n2; j++){
    free(H[j]);
  }
  free(H);

  for (j=0; j<n2; j++){
    free(TmpM[j]);
  }
  free(TmpM);

  for (i=0; i<n2; i++){
    free(S[i]);
  }
  free(S);

  free(M1);

  for (j=0; j<n2; j++){
    free(C[j]);
  }
  free(C);

  /* non-spin-orbit coupling and non-LDA+U */  
  if (SO_switch==0 && Hub_U_switch==0 && Constraint_NCS_switch==0 
      && Zeeman_NCS_switch==0 && Zeeman_NCO_switch==0){
    free(Dummy_ImNL[0][0][0][0]);
    free(Dummy_ImNL[0][0][0]);
    free(Dummy_ImNL[0][0]);
    free(Dummy_ImNL[0]);
    free(Dummy_ImNL);
  }

  for (ik=0; ik<=nkpath; ik++) {
    for (i_perk=0; i_perk<=n_perk[ik]; i_perk++) {
      free(EigenVal[ik][i_perk]);
    }
    free(EigenVal[ik]);
  }
  free(EigenVal);

  free(ik_list);
  free(i_perk_list);
  free(arpo);

  /* elapsed time */
  dtime(&TEtime);
}