#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h> 
  
struct timeval2 {
  long tv_sec;    /* second */
  long tv_usec;   /* microsecond */
};


void dtime(double *t)
{

  /* real time */
  struct timeval timev;
  gettimeofday(&timev, NULL);
  *t = timev.tv_sec + (double)timev.tv_usec*1e-6;

  /* user time + system time */
  /*
  float tarray[2];
  clock_t times(), wall;
  struct tms tbuf;
  wall = times(&tbuf);
  tarray[0] = (float) (tbuf.tms_utime / (float)CLK_TCK);
  tarray[1] = (float) (tbuf.tms_stime / (float)CLK_TCK);
  *t = (double) (tarray[0]+tarray[1]);
  printf("dtime: %lf\n",*t);
  */

}


