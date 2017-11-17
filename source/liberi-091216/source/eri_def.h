/*----------------------------------------------------------------------
  eri_def.h

  Coded by M. Toyoda, May 2009, JAIST/RCIS.
----------------------------------------------------------------------*/
#ifndef ERI_DEF_H_INCLUDED
#define ERI_DEF_H_INCLUDED


#define PI 3.1415926535897932384626

/*
 * 0 < ngl <= ERI_NGLMAX
 * 0 < ngird <= ERI_NGRIDMAX
 * 0 < lmax <= ERI_LMAXMAX
*/
#define ERI_NGLMAX 256
#define ERI_NGRIDMAX 8192
#define ERI_LMAXMAX 30
#define ERI_NQMAX   20 /* for logFSBT */

#ifdef DEBUG
  /* If DEBUG is defined, the following macros are available */
  #include <stdio.h>
  #define TRACE  ( msg )              fprintf(stderr, msg)
  #define TRACE1 ( msg , arg)         fprintf(stderr, msg, arg )
  #define TRACE2 ( msg , arg1, arg2 ) fprintf(stderr, msg, arg1, arg2)
  #define STEPTRACE( msg ) fprintf(stderr, "%s (%d): %s\n", __FILE__, __LINE__, msg )
#else
  /* DEBUG is not defined */
  #define TRACE( msg ) 
  #define TRACE1( msg , arg) 
  #define TRACE2( msg , arg1, arg2 ) 
  #define STEPTRACE_FILE( msg ) 
  #define STEPTRACE( msg ) 
#endif 


/* inline */
#if defined(__STDC_VERSION__) && __STDC_VERSION__>=199901L
#define ERI_INLINE inline
#else
#define ERI_INLINE static 
#endif


#endif /* ERI_DEF_H_INCLUDED */

