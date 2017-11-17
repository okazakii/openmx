/**********************************************************************
  neb_check.c:

     neb_check.c is a subroutine to check 
     whether the calculation is NEB or not.

  Log of neb_check.c:

    13/April/2011  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "openmx_common.h"
#include "Inputtools.h"

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
  

int neb_check(char *argv[]) 
{ 
  int i,j,flag;
  char *s_vec[40];
  int i_vec[40];

  if (input_open(argv[1])==0){
    MPI_Finalize(); 
    exit(0);
  }

  i=0;
  s_vec[i]="NOMD";                    i_vec[i]=0;  i++;
  s_vec[i]="NVE" ;                    i_vec[i]=1;  i++;
  s_vec[i]="NVT_VS";                  i_vec[i]=2;  i++; /* modified by mari */
  s_vec[i]="OPT";                     i_vec[i]=3;  i++;
  s_vec[i]="EF";                      i_vec[i]=4;  i++; 
  s_vec[i]="BFGS";                    i_vec[i]=5;  i++; 
  s_vec[i]="RF";                      i_vec[i]=6;  i++; /* RF method by hmweng */
  s_vec[i]="DIIS";                    i_vec[i]=7;  i++;
  s_vec[i]="Constraint_DIIS";         i_vec[i]=8;  i++; /* not used */
  s_vec[i]="NVT_NH";                  i_vec[i]=9;  i++; 
  s_vec[i]="Opt_LBFGS";               i_vec[i]=10; i++; 
  s_vec[i]="NVT_VS2";                 i_vec[i]=11; i++; /* modified by Ohwaki */
  s_vec[i]="EvsLC";                   i_vec[i]=12; i++; 
  s_vec[i]="NEB";                     i_vec[i]=13; i++; 

  j = input_string2int("MD.Type",&MD_switch, i, s_vec,i_vec);
  if (j==-1){
    MPI_Finalize();
    exit(0);
  }

  input_close();

  flag = 0;
  if (MD_switch==13) flag = 1;
   
  return flag;
}
