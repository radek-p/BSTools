
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"

/* /////////////////////////////////////////// */
/* degree elevation of Bezier curves */

boolean mbs_multiBCDegElevd ( int ncurves, int spdimen,
                              int inpitch, int indegree, const double *inctlpoints,
                              int deltadeg,
                              int outpitch, int *outdegree, double *outctlpoints )
{
  int   i, j, k, n, m;
  double *p, *q;

  if ( deltadeg < 0 ) {/* We can do the degree elevation, not reduction here. */
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_5, ERRMSG_5 );
    return false;
  }
    
/* If inctlpoints is NULL then we assume that data are already pointed by */
/* outctlpoints and we do not have to copy. But we may have to rearrange  */
/* the data if the final pitch is different than the initial pitch.       */
  if ( inctlpoints )
    pkv_Selectd ( ncurves, spdimen*(indegree+1), inpitch, outpitch,
                  inctlpoints, outctlpoints );
  else if ( outpitch != inpitch )
    pkv_Rearranged ( ncurves, spdimen*(indegree+1), inpitch, outpitch,
                     outctlpoints );

/* Now the proper degree elevation, by 1 with every loop execution. */
  for ( n = indegree; n < indegree+deltadeg; n++ ) {
    m = n+1;
    for ( k = 0, p = outctlpoints;  k < ncurves;  k++, p+= outpitch ) {
      q = p + n*spdimen;
      memcpy ( q+spdimen, q, spdimen*sizeof(double) );
      for ( i = n;  i > 0;  i--, q -= spdimen )
        for ( j = 0; j < spdimen; j++ )
          q[j] = ((double)i*q[j-spdimen]+(double)(m-i)*q[j])/(double)m;
    }
  }
  *outdegree = indegree+deltadeg;
  return true;
} /*mbs_multiBCDegElevd*/

