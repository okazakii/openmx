/**********************************************************************
  init_alloc_first.c:

     init_alloc_first.c is a subroutine to initialize an array 
     alloc_first[];

  Log of init_alloc_first.c:

     24/May/2003  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <math.h>
#include "openmx_common.h"

void init_alloc_first()
{

  /***********************************************
   truncation.c

    GListTAtoms1
    GListTAtoms2
  ***********************************************/
  alloc_first[0] = 1;

  /***********************************************
   truncation.c

  ***********************************************/
  alloc_first[1] = 1;

  /***********************************************
   truncation.c

    GridListAtom
    CellListAtom
  ***********************************************/
  alloc_first[2] = 1;

  /***********************************************
   truncation.c

    Density_Grid
    ADensity_Grid
    PCCDensity_Grid
    Vxc_Grid
    VNA_Grid
    dVHart_Grid
    Vpot_Grid
    Orbs_Grid
    COrbs_Grid
    dOrbs_Grid
    dCOrbs_Grid
  ***********************************************/
  alloc_first[3] = 1;

  /***********************************************
   truncation.c

     H0
     CntH0
     OLP
     CntOLP
     H
     CntH
     DS_NL
     CntDS_NL
     DM
     ResidualDM
     EDM
     PDM
     IOLP 
     CntCoes
  ***********************************************/
  alloc_first[4] = 1;

  /***********************************************
   truncation.c

     NumOLG
  ***********************************************/
  alloc_first[5] = 1;

  /***********************************************
   truncation.c

      RMI1
      RMI2
  ***********************************************/
  alloc_first[6] = 1;

  /***********************************************
   truncation.c

      ratv
      atv
      atv_ijk
  ***********************************************/
  alloc_first[7] = 1;

  /***********************************************
   Allocation_Arrays.c

      natn
      ncn
      Dis
  ***********************************************/
  alloc_first[8] = 1;

  /***********************************************
   Allocation_Arrays.c

      S
  ***********************************************/
  alloc_first[9] = 1;

  /***********************************************
   Set_Allocate_Atom2CPU.c

      M2G
  ***********************************************/
  alloc_first[10] = 1;

  /***********************************************
   in Set_Inf_SndRcv() of truncation.c

   Snd_MAN[numprocs][FS_Snd_Num[ID1]]
   Snd_GAN[numprocs][FS_Snd_Num[ID1]]
  ***********************************************/
  alloc_first[11] = 1;

  /***********************************************
   in Set_Inf_SndRcv() of truncation.c

   int Rcv_GAN[numprocs]
              [F_Rcv_Num[ID]+S_Rcv_Num[ID]]
  ***********************************************/
  alloc_first[12] = 1;

  /***********************************************
   Set_Allocate_Atom2CPU.c

      F_M2G
      S_M2G
  ***********************************************/
  alloc_first[13] = 1;

  /***********************************************
   allocate_grids2atoms() of truncation.c.

     My_Cell1
  ***********************************************/
  alloc_first[14] = 1;

  /***********************************************
   allocate_grids2atoms() of truncation.c.

     My_Cell0
     Cell_ID0
  ***********************************************/
  alloc_first[15] = 1;

  /***********************************************
   allocate_grids2atoms() of truncation.c.

     Num_Rcv_Grid1
     Num_Snd_Grid1
     Rcv_Grid1
     Snd_Grid1
  ***********************************************/
  alloc_first[16] = 1;

  /***********************************************
   allocate_grids2atoms() of truncation.c.

     Num_IRcv_Grid1
     Num_ISnd_Grid1
     IRcv_Grid1
     ISnd_Grid1
  ***********************************************/
  alloc_first[17] = 1;

  /***********************************************
   allocate_grids2atoms() of truncation.c.

     Rcv_FNAN2_MN
     Rcv_FNAN2_GA
     TopMAN2_Grid
     Num_Rcv_FNAN2_Grid
     Num_Snd_FNAN2_Grid
     Snd_FNAN2_At
     Snd_FNAN2_Nc
  ***********************************************/
  alloc_first[18] = 1;


  /***********************************************
   GDC_Allocation() of Set_Allocate_Atom2CPU.c.

  ***********************************************/
  alloc_first[19] = 1;

  /***********************************************
   GDC_Allocation() of Set_Allocate_Atom2CPU.c.

  ***********************************************/
  alloc_first[20] = 1;

  /***********************************************
   Setup_EC of truncation.c.

   NAtom_EC 
   MAtom_EC
   LAtom_EC
  ***********************************************/
  alloc_first[21] = 1;

  /***********************************************
   Set_Inf_SndRcv of truncation.c.

   Pro_Snd_GAtom
   Pro_Snd_MAtom
   Pro_Snd_LAtom
  ***********************************************/
  alloc_first[22] = 1;

  /***********************************************
   Generating_MP_Special_Kpt.c.

   NE_T_k_op
   NE_KGrids1
   NE_KGrids2
   NE_KGrids3
  ***********************************************/
  alloc_first[23] = 1;

  /***********************************************
   Generating_MP_Special_Kpt.c.

   Wannier_ProSpeName
   Wannier_ProName
   Wannier_Pos
   Wannier_X_Direction
   Wannier_Z_Direction
  ***********************************************/
  alloc_first[24] = 1;

  alloc_first[25] = 1;

}
