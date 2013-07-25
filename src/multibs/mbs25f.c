
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
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
/* Knot insertion for closed curves            */

int mbs_multiKnotInsClosedf ( int degree, int *lastknot, float *knots,
                              int ncurves, int spdimen, int inpitch, int outpitch,
                              float *ctlpoints, float t )
{
  float T;
  int   i, k, K, shift, lkn;

  k  = mbs_multiKnotInsf ( degree, lastknot, knots, ncurves, spdimen,
                           inpitch, outpitch, ctlpoints, t );
  lkn = *lastknot;
  K = lkn-2*degree;
  T  = knots[lkn-degree] - knots[degree];

  for ( i = k-degree+1; i <= k; i++ ) {
    if ( i+K < lkn-degree )
      shift = K*spdimen;
    else if ( i-K >= 0 )
      shift = -K*spdimen;
    else
      continue;
    pkv_Movef ( ncurves, spdimen, outpitch, shift, &ctlpoints[i*spdimen] );
  }
  for ( i = k; i+K <= lkn; i++ )
    knots[i+K] = knots[i]+T;
  for ( i = k+1; i-K >= 0; i-- )
    knots[i-K] = knots[i]-T;
  return k;
} /*mbs_multiKnotInsClosedf*/

