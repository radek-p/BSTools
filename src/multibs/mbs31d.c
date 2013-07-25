
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */
#include <stdlib.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"


/* testing polylines for monotonicity */

boolean mbs_MonotonicPolylined ( int spdimen, int npoints, int pitch,
                                 const double *points, const double *v )
{
  int         i, j;
  long double a, b;

  a = pkn_ScalarProductd ( spdimen, v, &points[0] );
  for ( i = 1, j = pitch;  i < npoints;  i++, j += pitch ) {
    b = pkn_ScalarProductd ( spdimen, v, &points[j] );
    if ( b < a )
      return false;
    a = b;
  }
  return true;
} /*mbs_MonotonicPolylined*/

boolean mbs_MonotonicPolylineRd ( int spdimen, int npoints, int pitch,
                                  const double *points, const double *v )
{
  int         i, j, s;
  long double a, b, c;

  s = spdimen-1;
  c = points[s];
  if ( c > 0.0 ) {
    a = pkn_ScalarProductd ( s, v, &points[0] ) / c;
    for ( i = 1, j = pitch;  i < npoints;  i++, j += pitch ) {
      c = points[j+s];
      if ( c <= 0.0 )
        return false;
      b = pkn_ScalarProductd ( s, v, &points[j] ) / c;
      if ( b < a )
        return false;
      a = b;
    }
    return true;
  }
  else if ( c < 0.0 ) {
    a = pkn_ScalarProductd ( s, v, &points[0] ) / c;
    for ( i = 1, j = pitch;  i < npoints;  i++, j += pitch ) {
      c = points[j+s];
      if ( c >= 0.0 )
        return false;
      b = pkn_ScalarProductd ( s, v, &points[j] ) / c;
      if ( b < a )
        return false;
      a = b;
    }
    return true;
  }
  else return false;
} /*mbs_MonotonicPolylineRd*/

