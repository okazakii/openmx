void Force_HNL(double *****CDM0)
{
  /****************************************************
                  Force arising from HNL
  ****************************************************/

  int Mc_AN,Gc_AN,Cwan,i,j,h_AN,q_AN,Mq_AN;
  int jan,kl,Qwan,Gq_AN,Gh_AN,Mh_AN,Hwan,ian;
  int l1,l2,l3,l,LL,Mul1,tno0,ncp;
  int tno1,tno2,size1,size2,n,kk,num,po,po1,po2;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  int **S_array,**R_array;
  int S_comm_flag,R_comm_flag;
  int SA_num,q,Sc_AN,GSc_AN;
  int Sc_wan,Sh_AN,GSh_AN,Sh_wan;
  int Sh_AN2,fan,jg,j0,jg0,Mj_AN0;
  int Original_Mc_AN;

  double rcutA,rcutB,rcut;
  double dEx,dEy,dEz,ene;
  double Stime_atom, Etime_atom;
  double **HNLx,**HNLy,**HNLz;
  int *Snd_DS_NL_Size,*Rcv_DS_NL_Size;  
  int *Indicator;
  Type_DS_NL *tmp_array;
  Type_DS_NL *tmp_array2;

  /* for OpenMP */
  int OMPID,Nthrds,Nthrds0,Nprocs,Nloop,ODNloop;
  int *OneD2h_AN,*OneD2q_AN;
  double *dEx_threads;
  double stime,etime;
  double stime1,etime1;

  MPI_Status stat;
  MPI_Request request;

  /* MPI */

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  dtime(&stime);

  /****************************
       allocation of arrays 
  *****************************/

  Indicator = (int*)malloc(sizeof(int)*numprocs);

  S_array = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    S_array[ID] = (int*)malloc(sizeof(int)*3);
  }

  R_array = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    R_array[ID] = (int*)malloc(sizeof(int)*3);
  }

  Snd_DS_NL_Size = (int*)malloc(sizeof(int)*numprocs);
  Rcv_DS_NL_Size = (int*)malloc(sizeof(int)*numprocs);

  /* initialize the temporal array storing the force contribution */

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = F_M2G[Mc_AN];
    Gxyz[Gc_AN][41] = 0.0;
    Gxyz[Gc_AN][42] = 0.0;
    Gxyz[Gc_AN][43] = 0.0;
  }

  /*************************************************************
                 contraction of DS_NL
  *************************************************************/

  if (Cnt_switch==1){
    for (so=0; so<(SO_switch+1); so++){
      Cont_Matrix1(DS_NL[so][0],CntDS_NL[so][0]);
      Cont_Matrix1(DS_NL[so][1],CntDS_NL[so][1]);
      Cont_Matrix1(DS_NL[so][2],CntDS_NL[so][2]);
      Cont_Matrix1(DS_NL[so][3],CntDS_NL[so][3]);
    }
  }

  /*****************************************}********************** 
      THE FIRST CASE:
      In case of I=i or I=j 
      for d [ \sum_k <i|k>ek<k|j> ]/dRI  
  ****************************************************************/

  /*******************************************************
   *******************************************************
       multiplying overlap integrals WITH COMMUNICATION
   *******************************************************
   *******************************************************/

  MPI_Barrier(mpi_comm_level1);
  dtime(&stime);

  HNLx = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (j=0; j<List_YOUSO[7]; j++){
    HNLx[j] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
  }

  HNLy = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (j=0; j<List_YOUSO[7]; j++){
    HNLy[j] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
  }

  HNLz = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
  for (j=0; j<List_YOUSO[7]; j++){
    HNLz[j] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
  }

  for (ID=0; ID<numprocs; ID++){
    F_Snd_Num_WK[ID] = 0;
    F_Rcv_Num_WK[ID] = 0;
  }

  do {

    /***********************************                                                            
            set the size of data                                                                      
    ************************************/

    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      /* find the data size to send the block data */

      if ( 0<(F_Snd_Num[IDS]-F_Snd_Num_WK[IDS]) ){

	size1 = 0;
	n = F_Snd_Num_WK[IDS];

	Mc_AN = Snd_MAN[IDS][n];
	Gc_AN = Snd_GAN[IDS][n];
	Cwan = WhatSpecies[Gc_AN];
	tno1 = Spe_Total_NO[Cwan];

	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	  Gh_AN = natn[Gc_AN][h_AN];
	  Hwan = WhatSpecies[Gh_AN];
	  tno2 = Spe_Total_VPS_Pro[Hwan];
	  size1 += (VPS_j_dependency[Hwan]+1)*tno1*tno2;
	}

	Snd_DS_NL_Size[IDS] = size1;
	MPI_Isend(&size1, 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
      }
      else{
	Snd_DS_NL_Size[IDS] = 0;
      }

      /* receiving of the size of the data */

      if ( 0<(F_Rcv_Num[IDR]-F_Rcv_Num_WK[IDR]) ){
	MPI_Recv(&size2, 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
	Rcv_DS_NL_Size[IDR] = size2;
      }
      else{
	Rcv_DS_NL_Size[IDR] = 0;
      }

      if ( 0<(F_Snd_Num[IDS]-F_Snd_Num_WK[IDS]) )  MPI_Wait(&request,&stat);

    } /* ID */

      /***********************************
                 data transfer
      ************************************/

    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      /******************************
            sending of the data 
      ******************************/

      if ( 0<(F_Snd_Num[IDS]-F_Snd_Num_WK[IDS]) ){

	size1 = Snd_DS_NL_Size[IDS];

	/* allocation of the array */

	tmp_array = (double*)malloc(sizeof(double)*size1);

	/* multidimentional array to the vector array */

	num = 0;
	n = F_Snd_Num_WK[IDS];

	Mc_AN = Snd_MAN[IDS][n];
	Gc_AN = Snd_GAN[IDS][n];
	Cwan = WhatSpecies[Gc_AN]; 
	tno1 = Spe_Total_NO[Cwan];

	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	  Gh_AN = natn[Gc_AN][h_AN];        
	  Hwan = WhatSpecies[Gh_AN];
	  tno2 = Spe_Total_VPS_Pro[Hwan];

	  for (so=0; so<=VPS_j_dependency[Hwan]; so++){
	    for (i=0; i<tno1; i++){
	      for (j=0; j<tno2; j++){
		tmp_array[num] = DS_NL[so][0][Mc_AN][h_AN][i][j];
		num++;
	      } 
	    } 
	  }
	}

	MPI_Isend(&tmp_array[0], size1, MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request);
      }

      /******************************
          receiving of the block data
      ******************************/

      if ( 0<(F_Rcv_Num[IDR]-F_Rcv_Num_WK[IDR]) ){
        
	size2 = Rcv_DS_NL_Size[IDR];
	tmp_array2 = (double*)malloc(sizeof(double)*size2);
	MPI_Recv(&tmp_array2[0], size2, MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

	/* store */

	num = 0;
	n = F_Rcv_Num_WK[IDR];
	Original_Mc_AN = F_TopMAN[IDR] + n;

	Gc_AN = Rcv_GAN[IDR][n];
	Cwan = WhatSpecies[Gc_AN];
	tno1 = Spe_Total_NO[Cwan];
	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	  Gh_AN = natn[Gc_AN][h_AN];
	  Hwan = WhatSpecies[Gh_AN];
	  tno2 = Spe_Total_VPS_Pro[Hwan];

	  for (so=0; so<=VPS_j_dependency[Hwan]; so++){
	    for (i=0; i<tno1; i++){
	      for (j=0; j<tno2; j++){
		DS_NL[so][0][Matomnum+1][h_AN][i][j] = tmp_array2[num];
		num++;
	      }
	    }
	  }
	}

	/* free tmp_array2 */
	free(tmp_array2);

	/*****************************************************************
                             multiplying overlap integrals
	*****************************************************************/

	for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

	  dtime(&Stime_atom);

	  dEx = 0.0; 
	  dEy = 0.0; 
	  dEz = 0.0; 

	  Gc_AN = M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  fan = FNAN[Gc_AN];

	  h_AN = 0;
	  Gh_AN = natn[Gc_AN][h_AN];
	  Mh_AN = F_G2M[Gh_AN];
	  Hwan = WhatSpecies[Gh_AN];
	  ian = Spe_Total_CNO[Hwan];

	  n = F_Rcv_Num_WK[IDR];
	  jg = Rcv_GAN[IDR][n];

	  for (j0=0; j0<=fan; j0++){

	    jg0 = natn[Gc_AN][j0];
	    Mj_AN0 = F_G2M[jg0];

	    po2 = 0;
	    if (Original_Mc_AN==Mj_AN0){
	      po2 = 1;
	      q_AN = j0;
	    }

	    if (po2==1){

	      Gq_AN = natn[Gc_AN][q_AN];
	      Mq_AN = F_G2M[Gq_AN];
	      Qwan = WhatSpecies[Gq_AN];
	      jan = Spe_Total_CNO[Qwan];
	      kl = RMI1[Mc_AN][h_AN][q_AN];

	      if (Cnt_switch==0) {
		dHNL(0,Mc_AN,h_AN,q_AN,DS_NL,HNLx,HNLy,HNLz);
	      }
	      else { 
		dHNL(0,Mc_AN,h_AN,q_AN,CntDS_NL,HNLx,HNLy,HNLz);
	      }

	      /* contribution of force = Trace(CDM0*dH) */
	      /* spin non-polarization */

	      if (SpinP_switch==0){

		for (i=0; i<ian; i++){
		  for (j=0; j<jan; j++){
		    if (q_AN==h_AN){

		      dEx += 2.0*CDM0[0][Mh_AN][kl][i][j]*HNLx[0][i][j].r;
		      dEy += 2.0*CDM0[0][Mh_AN][kl][i][j]*HNLy[0][i][j].r;
		      dEz += 2.0*CDM0[0][Mh_AN][kl][i][j]*HNLz[0][i][j].r;
		    }
		    else{

		      dEx += 4.0*CDM0[0][Mh_AN][kl][i][j]*HNLx[0][i][j].r;
		      dEy += 4.0*CDM0[0][Mh_AN][kl][i][j]*HNLy[0][i][j].r;
		      dEz += 4.0*CDM0[0][Mh_AN][kl][i][j]*HNLz[0][i][j].r;
		    }
		  }
		}
	      }

	      /* else */

	      else{

		for (i=0; i<ian; i++){
		  for (j=0; j<jan; j++){
		    if (q_AN==h_AN){
		      dEx += (  CDM0[0][Mh_AN][kl][i][j]
				+ CDM0[1][Mh_AN][kl][i][j] )*HVNAx[i][j];
		    }
		    else{
		      dEx += 2.0*(  CDM0[0][Mh_AN][kl][i][j]
				    + CDM0[1][Mh_AN][kl][i][j] )*HVNAx[i][j];
		    } 
		  }
		}
	      }

	    } /* if (po2==1) */
	  } /* j0 */             

	    /* force from #4B */

	  Gxyz[Gc_AN][40+kk] += dEx;

	  /* timing */
	  dtime(&Etime_atom);
	  time_per_atom[Gc_AN] += Etime_atom - Stime_atom;

	} /* Mc_AN */

	  /********************************************
            increment of F_Rcv_Num_WK[IDR] 
	  ********************************************/

	F_Rcv_Num_WK[IDR]++;

      } /* if ( 0<(F_Snd_Num[IDS]-F_Snd_Num_WK[IDS]) ) */

      if ( 0<(F_Snd_Num[IDS]-F_Snd_Num_WK[IDS]) ) {

	MPI_Wait(&request,&stat);
	free(tmp_array);  /* freeing of array */

	/********************************************
             increment of F_Snd_Num_WK[IDS]
	********************************************/

	F_Snd_Num_WK[IDS]++;
      } 

    } /* ID */

      /*****************************************************
        check whether all the communications have finished
      *****************************************************/

    po = 0;
    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      if ( 0<(F_Snd_Num[IDS]-F_Snd_Num_WK[IDS]) ) po += F_Snd_Num[IDS]-F_Snd_Num_WK[IDS];
      if ( 0<(F_Rcv_Num[IDR]-F_Rcv_Num_WK[IDR]) ) po += F_Rcv_Num[IDR]-F_Rcv_Num_WK[IDR];
    }

  } while (po!=0);

  for (j=0; j<List_YOUSO[7]; j++){
    free(HNLx[j]);
  }
  free(HNLx);

  for (j=0; j<List_YOUSO[7]; j++){
    free(HNLy[j]);
  }
  free(HNLy);

  for (j=0; j<List_YOUSO[7]; j++){
    free(HNLz[j]);
  }
  free(HNLz);

  dtime(&etime);
  if(myid==0 && measure_time){
    printf("Time for part3 of force#4=%18.5f\n",etime-stime);fflush(stdout);
  } 

  /*******************************************************
   *******************************************************
      THE FIRST CASE:
      multiplying overlap integrals WITHOUT COMMUNICATION
   *******************************************************
   *******************************************************/

  dtime(&stime);

#pragma omp parallel shared(time_per_atom,Gxyz,CDM0,SpinP_switch,CntHVNA2,HVNA2,DS_VNA,Cnt_switch,RMI1,FNAN,Spe_Total_CNO,WhatSpecies,F_G2M,natn,M2G,Matomnum,List_YOUSO) private(HVNAx,OMPID,Nthrds,Nprocs,Mc_AN,Stime_atom,Etime_atom,dEx,Gc_AN,h_AN,Gh_AN,Mh_AN,Hwan,ian,q_AN,Gq_AN,Mq_AN,Qwan,jan,kl,i,j)
  {

    /* allocation of array */

    HVNAx = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
    for (j=0; j<List_YOUSO[7]; j++){
      HVNAx[j] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    }

    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    for (Mc_AN=(OMPID*Matomnum/Nthrds+1); Mc_AN<((OMPID+1)*Matomnum/Nthrds+1); Mc_AN++){

      dtime(&Stime_atom);

      dEx = 0.0; 

      Gc_AN = M2G[Mc_AN];
      h_AN = 0;
      Gh_AN = natn[Gc_AN][h_AN];
      Mh_AN = F_G2M[Gh_AN];
      Hwan = WhatSpecies[Gh_AN];
      ian = Spe_Total_CNO[Hwan];

      for (q_AN=h_AN; q_AN<=FNAN[Gc_AN]; q_AN++){

	Gq_AN = natn[Gc_AN][q_AN];
	Mq_AN = F_G2M[Gq_AN];

	if (Mq_AN<=Matomnum){

	  Qwan = WhatSpecies[Gq_AN];
	  jan = Spe_Total_CNO[Qwan];
	  kl = RMI1[Mc_AN][h_AN][q_AN];

	  if (Cnt_switch==0) {
	    dHVNA(kk,0,Mc_AN,h_AN,q_AN,DS_VNA,HVNA2,HVNAx);
	  }
	  else { 
	    dHVNA(kk,0,Mc_AN,h_AN,q_AN,DS_VNA,CntHVNA2,HVNAx);
	  }

	  if (SpinP_switch==0){

	    for (i=0; i<ian; i++){
	      for (j=0; j<jan; j++){
		if (q_AN==h_AN){
		  dEx += 2.0*CDM0[0][Mh_AN][kl][i][j]*HVNAx[i][j];
		}
		else{
		  dEx += 4.0*CDM0[0][Mh_AN][kl][i][j]*HVNAx[i][j];
		}
	      }
	    }
	  }

	  /* else */

	  else{

	    for (i=0; i<ian; i++){
	      for (j=0; j<jan; j++){
		if (q_AN==h_AN){
		  dEx += (  CDM0[0][Mh_AN][kl][i][j]
			    + CDM0[1][Mh_AN][kl][i][j] )*HVNAx[i][j];
		}
		else{
		  dEx += 2.0*(  CDM0[0][Mh_AN][kl][i][j]
			        + CDM0[1][Mh_AN][kl][i][j] )*HVNAx[i][j];
		} 
	      }
	    }
	  }
	}
      }

      /* force from #4B */

      Gxyz[Gc_AN][40+kk] += dEx;

      /* timing */
      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;

    } /* Mc_AN */

      /* freeing of array */

    for (j=0; j<List_YOUSO[7]; j++){
      free(HVNAx[j]);
    }
    free(HVNAx);

  } /* #pragma omp parallel */

  dtime(&etime);
  if(myid==0 && measure_time){
    printf("Time for part4 of force#4=%18.5f\n",etime-stime);fflush(stdout);
  } 
         
  /*************************************************************
     THE SECOND CASE:
     In case of I=k with I!=i and I!=j
     d [ \sum_k <i|k>ek<k|j> ]/dRI  
  *************************************************************/

  /************************************************************ 
     MPI communication of DS_VNA whose basis part is not located 
     on own site but projector part is located on own site. 
  ************************************************************/

  MPI_Barrier(mpi_comm_level1);
  dtime(&stime);

  for (ID=0; ID<numprocs; ID++) Indicator[ID] = 0;

  for (Mc_AN=1; Mc_AN<=Max_Matomnum; Mc_AN++){

    dtime(&Stime_atom);

    dtime(&stime1);

    if (Mc_AN<=Matomnum){ 
      Gc_AN = M2G[Mc_AN];
      Cwan = WhatSpecies[Gc_AN];
    }

    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      i = Indicator[IDS]; 
      po = 0;

      Gh_AN = Pro_Snd_GAtom[IDS][i]; 

      if (Gh_AN!=0){

	/* find the range with the same global atomic number */

	do {

	  i++;
	  if (Gh_AN!=Pro_Snd_GAtom[IDS][i]) po = 1;
	} while(po==0);

	i--;
	SA_num = i - Indicator[IDS] + 1;

	/* find the data size to send the block data */

	size1 = 0;
	for (q=Indicator[IDS]; q<=(Indicator[IDS]+SA_num-1); q++){
	  Sc_AN = Pro_Snd_MAtom[IDS][q]; 
	  GSc_AN = F_M2G[Sc_AN];
	  Sc_wan = WhatSpecies[GSc_AN];
	  tno1 = Spe_Total_CNO[Sc_wan];
	  tno2 = (List_YOUSO[35]+1)*(List_YOUSO[35]+1)*List_YOUSO[34];
	  size1 += 2*tno1*tno2;
	  size1 += 3;
	}

      } /* if (Gh_AN!=0) */

      else {
	SA_num = 0;
	size1 = 0;
      }
        
      S_array[IDS][0] = Gh_AN;
      S_array[IDS][1] = SA_num;
      S_array[IDS][2] = size1;

      if (ID!=0){
	MPI_Isend(&S_array[IDS][0], 3, MPI_INT, IDS, tag, mpi_comm_level1, &request);
	MPI_Recv( &R_array[IDR][0], 3, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
	MPI_Wait(&request,&stat);
      }
      else {
	R_array[myid][0] = S_array[myid][0];
	R_array[myid][1] = S_array[myid][1];
	R_array[myid][2] = S_array[myid][2];
      }

      if (R_array[IDR][0]==Gc_AN) R_comm_flag = 1;
      else                        R_comm_flag = 0;

      if (ID!=0){
	MPI_Isend(&R_comm_flag, 1, MPI_INT, IDR, tag, mpi_comm_level1, &request);
	MPI_Recv( &S_comm_flag, 1, MPI_INT, IDS, tag, mpi_comm_level1, &stat);
	MPI_Wait(&request,&stat);
      }
      else{
	S_comm_flag = R_comm_flag;
      }

      /*****************************************
                       send the data
      *****************************************/
        
      /* if (S_comm_flag==1) then, send data to IDS */
         
      if (S_comm_flag==1){

	/* allocate tmp_array */

	tmp_array = (Type_DS_VNA*)malloc(sizeof(Type_DS_VNA)*size1);

	/* multidimentional array to vector array */

	num = 0;

	for (q=Indicator[IDS]; q<=(Indicator[IDS]+SA_num-1); q++){

	  Sc_AN = Pro_Snd_MAtom[IDS][q]; 
	  GSc_AN = F_M2G[Sc_AN];
	  Sc_wan = WhatSpecies[GSc_AN];
	  tno1 = Spe_Total_CNO[Sc_wan];

	  Sh_AN = Pro_Snd_LAtom[IDS][q]; 
	  GSh_AN = natn[GSc_AN][Sh_AN];
	  Sh_wan = WhatSpecies[GSh_AN];
	  tno2 = (List_YOUSO[35]+1)*(List_YOUSO[35]+1)*List_YOUSO[34];

	  Sh_AN2 = Pro_Snd_LAtom2[IDS][q]; 

	  tmp_array[num] = (Type_DS_VNA)Sc_AN;  num++;
	  tmp_array[num] = (Type_DS_VNA)Sh_AN;  num++; 
	  tmp_array[num] = (Type_DS_VNA)Sh_AN2; num++; 

	  for (i=0; i<tno1; i++){
	    for (j=0; j<tno2; j++){
	      tmp_array[num] = DS_VNA[0][Sc_AN][Sh_AN][i][j];
	      num++;
	    }
	  }

	  for (i=0; i<tno1; i++){
	    for (j=0; j<tno2; j++){
	      tmp_array[num] = DS_VNA[1][Sc_AN][Sh_AN][i][j];
	      num++;
	    }
	  }
         
	}

	if (ID!=0){
	  MPI_Isend(&tmp_array[0], size1, MPI_Type_DS_VNA, IDS, tag, mpi_comm_level1, &request);
	}

	/* update Indicator[IDS] */

	Indicator[IDS] += SA_num;

      } /* if (S_comm_flag==1) */ 

        /*****************************************
                     receive the data
	*****************************************/

        /* if (R_comm_flag==1) then, receive the data from IDR */

      if (R_comm_flag==1){

	size2 = R_array[IDR][2]; 
	tmp_array2 = (Type_DS_VNA*)malloc(sizeof(Type_DS_VNA)*size2);

	if (ID!=0){
	  MPI_Recv(&tmp_array2[0], size2, MPI_Type_DS_VNA, IDR, tag, mpi_comm_level1, &stat);
	}
	else{
	  for (i=0; i<size2; i++) tmp_array2[i] = tmp_array[i];
	}

	/* store */

	num = 0;

	for (n=0; n<R_array[IDR][1]; n++){
            
	  Sc_AN  = (int)tmp_array2[num]; num++;
	  Sh_AN  = (int)tmp_array2[num]; num++;
	  Sh_AN2 = (int)tmp_array2[num]; num++;

	  GSc_AN = natn[Gc_AN][Sh_AN2]; 
	  Sc_wan = WhatSpecies[GSc_AN]; 

	  tno1 = Spe_Total_CNO[Sc_wan];
	  tno2 = (List_YOUSO[35]+1)*(List_YOUSO[35]+1)*List_YOUSO[34];

	  for (i=0; i<tno1; i++){
	    for (j=0; j<tno2; j++){
	      DS_VNA[0][Matomnum+1][Sh_AN2][i][j] = tmp_array2[num];
	      num++;
	    }
	  }

	  for (i=0; i<tno1; i++){
	    for (j=0; j<tno2; j++){
	      DS_VNA[1][Matomnum+1][Sh_AN2][i][j] = tmp_array2[num];
	      num++;
	    }
	  }
	}         

	/* free tmp_array2 */
	free(tmp_array2);
 
      } /* if (R_comm_flag==1) */

      if (S_comm_flag==1){
	if (ID!=0) MPI_Wait(&request,&stat);
	free(tmp_array);  /* freeing of array */
      }
      
    } /* ID */

    dtime(&etime1);
    if(myid==0 && measure_time){
      printf("Time for part5A of force#4=%18.5f\n",etime1-stime1);fflush(stdout);
    } 

    dtime(&stime1);

    if (Mc_AN<=Matomnum){ 

      /* get Nthrds0 */  
#pragma omp parallel shared(Nthrds0)
      {
	Nthrds0 = omp_get_num_threads();
      }

      /* allocation of array */
      dEx_threads = (double*)malloc(sizeof(double)*Nthrds0);
      for (Nloop=0; Nloop<Nthrds0; Nloop++) dEx_threads[Nloop] = 0.0;

      /* one-dimensionalize the h_AN and q_AN loops */ 

      OneD2h_AN = (int*)malloc(sizeof(int)*(FNAN[Gc_AN]+1)*(FNAN[Gc_AN]+2)/2);
      OneD2q_AN = (int*)malloc(sizeof(int)*(FNAN[Gc_AN]+1)*(FNAN[Gc_AN]+2)/2);

      ODNloop = 0;
      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	for (q_AN=h_AN; q_AN<=FNAN[Gc_AN]; q_AN++){

	  kl = RMI1[Mc_AN][h_AN][q_AN];

	  if (0<=kl){
	    OneD2h_AN[ODNloop] = h_AN;
	    OneD2q_AN[ODNloop] = q_AN; 
	    ODNloop++;      
	  }
	}
      }

#pragma omp parallel shared(ODNloop,OneD2h_AN,OneD2q_AN,Mc_AN,Gc_AN,kk,dEx_threads,CDM0,SpinP_switch,CntHVNA2,HVNA2,DS_VNA,Cnt_switch,RMI1,Spe_Total_CNO,WhatSpecies,F_G2M,natn,FNAN,List_YOUSO) private(OMPID,Nthrds,Nprocs,HVNAx,i,j,h_AN,Gh_AN,Mh_AN,Hwan,ian,q_AN,Gq_AN,Mq_AN,Qwan,jan,kl,Nloop)
      {
          
	/* allocation of array */
           
	HVNAx = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
	for (j=0; j<List_YOUSO[7]; j++){
	  HVNAx[j] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
	}

	/* get info. on OpenMP */ 
          
	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	for (Nloop=OMPID*ODNloop/Nthrds; Nloop<(OMPID+1)*ODNloop/Nthrds; Nloop++){

	  /* get h_AN and q_AN */

	  h_AN = OneD2h_AN[Nloop];
	  q_AN = OneD2q_AN[Nloop];

	  /* set informations on h_AN */

	  Gh_AN = natn[Gc_AN][h_AN];
	  Mh_AN = F_G2M[Gh_AN];
	  Hwan = WhatSpecies[Gh_AN];
	  ian = Spe_Total_CNO[Hwan];

	  /* set informations on q_AN */

	  Gq_AN = natn[Gc_AN][q_AN];
	  Mq_AN = F_G2M[Gq_AN];
	  Qwan = WhatSpecies[Gq_AN];
	  jan = Spe_Total_CNO[Qwan];
	  kl = RMI1[Mc_AN][h_AN][q_AN];

	  if (0<=kl){

	    if (Cnt_switch==0)
	      dHVNA(kk,1,Mc_AN,h_AN,q_AN,DS_VNA,HVNA2,HVNAx);
	    else 
	      dHVNA(kk,1,Mc_AN,h_AN,q_AN,DS_VNA,CntHVNA2,HVNAx);

	    /* contribution of force = Trace(CDM0*dH) */

	    /* spin non-polarization */

	    if (SpinP_switch==0){

	      for (i=0; i<ian; i++){
		for (j=0; j<jan; j++){
		  if (q_AN==h_AN){
		    dEx_threads[OMPID] += 2.0*CDM0[0][Mh_AN][kl][i][j]*HVNAx[i][j];
		  }
		  else{
		    dEx_threads[OMPID] += 4.0*CDM0[0][Mh_AN][kl][i][j]*HVNAx[i][j];
		  }
		}
	      }
	    }

	    /* else */

	    else{

	      for (i=0; i<ian; i++){
		for (j=0; j<jan; j++){
		  if (q_AN==h_AN){
		    dEx_threads[OMPID] += (   CDM0[0][Mh_AN][kl][i][j]
				              + CDM0[1][Mh_AN][kl][i][j] )*HVNAx[i][j];
		  }
		  else{
		    dEx_threads[OMPID] += 2.0*(   CDM0[0][Mh_AN][kl][i][j]
				                  + CDM0[1][Mh_AN][kl][i][j] )*HVNAx[i][j];
		  } 
		}
	      }
	    }

	  } /* if (0<=kl) */

	} /* Nloop */

	  /* freeing of array */

	for (j=0; j<List_YOUSO[7]; j++){
	  free(HVNAx[j]);
	}
	free(HVNAx);

      } /* #pragma omp parallel */

	/* sum of dEx_threads */

      dEx = 0.0;
      for (Nloop=0; Nloop<Nthrds0; Nloop++){
	dEx += dEx_threads[Nloop];
      }

      /* force from #4B */

      Gxyz[Gc_AN][40+kk] += dEx;

      /* timing */
      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;

      /* freeing of array */
      free(OneD2q_AN);
      free(OneD2h_AN);
      free(dEx_threads);

    } /* if (Mc_AN<=Matomnum) */

    dtime(&etime1);
    if(myid==0 && measure_time){
      printf("Time for part5B of force#4=%18.5f\n",etime1-stime1);fflush(stdout);
    } 

  } /* Mc_AN */

  dtime(&etime);
  if(myid==0 && measure_time){
    printf("Time for part5 of force#4=%18.5f\n",etime-stime);fflush(stdout);
  } 


  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];

    if (2<=level_stdout){
      printf("<Force>  force(4B) myid=%2d  Mc_AN=%2d Gc_AN=%2d  %15.12f %15.12f %15.12f\n",
	     myid,Mc_AN,Gc_AN,Gxyz[Gc_AN][41],Gxyz[Gc_AN][42],Gxyz[Gc_AN][43]);fflush(stdout);
    }

    Gxyz[Gc_AN][17] += Gxyz[Gc_AN][41];
    Gxyz[Gc_AN][18] += Gxyz[Gc_AN][42];
    Gxyz[Gc_AN][19] += Gxyz[Gc_AN][43];
  }

  /***********************************
            freeing of arrays 
  ************************************/

  free(Indicator);

  for (ID=0; ID<numprocs; ID++){
    free(S_array[ID]);
  }
  free(S_array);

  for (ID=0; ID<numprocs; ID++){
    free(R_array[ID]);
  }
  free(R_array);

  free(Snd_DS_VNA_Size);
  free(Rcv_DS_VNA_Size);
}
