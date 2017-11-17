/*----------------------------------------------------------------------
  demo.c
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "eri.h"
#include "demo_sub.h"


extern void demo_test1(double I4[3][3], double sec[2]);
extern void demo_test2(double I4[5][3], double sec[2]);
extern void demo_test3(double err[2], double sec[2]);
extern void demo_test4(double I4[10][3], double sec[2]);



int main(void)
{
  int i;
  double I4_test1[3][3], I4_test2[5][3], deriv_err[2], I4_test4[10][3];
  double sec1[2], sec2[2], sec3[2], sec4[2];

  const char* card1[3] = {"(ss|ss)", "(ps|sp)", "(dd|dd)"};
  const char* card2[5] = { "(12|34)", "(13|42)", "(14|23)", 
                           "(21|34)", "(31|24)"};

  double err1, err2, ERRTHRESH;
  int flag;

  flag = 1;
  ERRTHRESH = 1e-5;

  { 
    printf("----------------------------------------------------------------\n");
    printf("  LIBERI DEMO PROGRAM\n");
    printf("----------------------------------------------------------------\n");
    printf("LIBERI VERSION %s\n", ERI_Version());
    printf("\n");
  }

  demo_test1(I4_test1, sec1);
  demo_test2(I4_test2, sec2);
  demo_test3(deriv_err, sec3); 
  demo_test4(I4_test4, sec4);

  {
    printf("\n");
    printf("----------------------------------------------------------------\n");
    printf("  SUMMARY\n");
    printf("----------------------------------------------------------------\n");
    printf("GTO:\n");
    printf("           EXACT       LINEAR-SBT (ERR)    LOG-SBT    (ERR)\n");
    for (i=0; i<3; i++) {
      err1 = fabs(I4_test1[i][0]-I4_test1[i][2]);
      err2 = fabs(I4_test1[i][1]-I4_test1[i][2]);
      if (err1 > ERRTHRESH || err2 > ERRTHRESH) flag=0;
      printf("  %s  %10.7f  %10.7f (%5.0e)  %10.7f (%5.0e)\n",
           card1[i], I4_test1[i][2], I4_test1[i][0], err1, 
           I4_test1[i][1], err2 );
    }
    printf("  TIME (LINEAR-SBT) = %8.4f SEC\n", sec1[0]);
    printf("  TIME (LOG-SBT)    = %8.4f SEC\n", sec1[1]);
    printf("\n");
    printf("PAO:\n");
    printf("           REFERENCE   LINEAR-SBT (ERR)    LOG-SBT    (ERR)\n");
    for (i=0; i<5; i++) {
      err1 = fabs(I4_test2[i][0]-I4_test2[i][2]);
      err2 = fabs(I4_test2[i][1]-I4_test2[i][2]);
      if (err1 > ERRTHRESH || err2 > ERRTHRESH) flag=0;
      printf("  %s  %10.7f  %10.7f (%5.0e)  %10.7f (%5.0e)\n",
           card2[i], I4_test2[i][2], 
           I4_test2[i][0], err1, I4_test2[i][1], err2);
    }
    printf("  TIME (LINEAR-SBT) = %8.4f SEC\n", sec2[0]);
    printf("  TIME (LOG-SBT)    = %8.4f SEC\n", sec2[1]);
    printf("\n");
    printf("DERIVATIVES:\n");
    if (deriv_err[0]>ERRTHRESH || deriv_err[1]>ERRTHRESH) flag=0;
    printf("  MAXERR (LINEAR-SBT) = %5.0e\n", deriv_err[0]);
    printf("  MAXERR (LOG-SBT)    = %5.0e\n", deriv_err[1]);
    printf("  TIME   (LINEAR-SBT) = %8.4f SEC\n", sec3[0]);
    printf("  TIME   (LOG-SBT)    = %8.4f SEC\n", sec3[1]);
    printf("\n");
    printf("SCREENING:\n");
    printf("  OMEGA  EXACT       LINEAR-SBT (ERR)   LOG-SBT    (ERR)\n");
    for (i=0; i<10; i++) {
      err1 = fabs(I4_test4[i][0]-I4_test4[i][2]);
      err2 = fabs(I4_test4[i][1]-I4_test4[i][2]);
      if (err1>ERRTHRESH || err2>ERRTHRESH) flag=0;
      printf("  %5.2f  %10.7f  %10.7f (%5.0e)  %10.7f (%5.0e)\n",
        0.05*(double)(i+1), I4_test4[i][2], 
        I4_test4[i][0], err1, I4_test4[i][1], err2);
    }
    printf("  TIME (LINEAR-SBT) = %8.4f SEC\n", sec4[0]);
    printf("  TIME (LOG-SBT)    = %8.4f SEC\n", sec4[1]);
    printf("\n");

  }

  /* result */
  if (flag) {
    printf("SELF-CHECK TEST HAS COMPLETED SUCCESSFULLY.\n");
  } else {
    printf("SELF-CHECK TEST HAS COMPLETED.\n");
    printf("POOR NUMERICAL ACCURACY WAS FOUND (SEE ABOVE).\n");
  }
  printf("\n");
   
  return 0;
}


/* EOF */
