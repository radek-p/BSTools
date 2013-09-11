
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

/* ////////////////////////////////////////// */
/* change of n+1 leftmost knots in a B-spline representation of curves, */
/* using the Oslo algorithm */
boolean mbs_multiBSChangeLeftKnotsf ( int ncurves, int spdimen, int degree,
                                      float *knots, int pitch, float *ctlpoints,
                                      float *newknots )
{
  void   *sp;
  int    i, j, k;
  float  *d, *e, *t;
  double alpha;

  sp = pkv_GetScratchMemTop ();
  t = pkv_GetScratchMemf ( 2*degree );
  d = pkv_GetScratchMemf ( ncurves*spdimen*(degree+1) );
  e = pkv_GetScratchMemf ( ncurves*spdimen*(degree+1) );
  if ( !t || !d || !e ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  memcpy ( t, newknots, (degree+1)*sizeof(float) );
  memcpy ( &t[degree+1], &knots[degree+1], (degree-1)*sizeof(float) );
  pkv_Selectf ( ncurves, spdimen*(degree+1), pitch, spdimen*(degree+1),
                ctlpoints, e );
  for ( k = 0; k < degree; k++ ) {
    memcpy ( d, e, ncurves*spdimen*(degree+1)*sizeof(float) );
    for ( j = 1; j <= degree; j++ )
      for ( i = degree; i >= j; i-- ) {
        alpha = (t[k+j]-knots[i])/
                (knots[i+degree+1-j]-knots[i]);
        pkn_MatrixLinCombf ( ncurves, spdimen,
                             spdimen*(degree+1), &d[(i-1)*spdimen], 1.0-alpha,
                             spdimen*(degree+1), &d[i*spdimen], alpha,
                             spdimen*(degree+1), &d[i*spdimen] );
      }
    pkv_Selectf ( ncurves, spdimen, spdimen*(degree+1), pitch,
                  &d[degree*spdimen], &ctlpoints[k*spdimen] );
  }
  memcpy ( knots, newknots, (degree+1)*sizeof(float) );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_multiBSChangeLeftKnotsf*/

/* ////////////////////////////////////////// */
/* change of n+1 rightmost knots in a B-spline representation of curves */

boolean mbs_multiBSChangeRightKnotsf ( int ncurves, int spdimen, int degree,
                                       int lastknot, float *knots,
                                       int pitch, float *ctlpoints,
                                       float *newknots )
{
  void   *sp;
  int    i, j, k;
  float  *d, *e, *t;
  double alpha;

  sp = pkv_GetScratchMemTop ();
  t = pkv_GetScratchMemf ( 2*degree );
  d = pkv_GetScratchMemf ( ncurves*spdimen*(degree+1) );
  e = pkv_GetScratchMemf ( ncurves*spdimen*(degree+1) );
  if ( !t || !d || !e ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  memcpy ( t, &knots[lastknot-2*degree+1], (degree-1)*sizeof(float) );
  memcpy ( &t[degree-1], newknots, (degree+1)*sizeof(float) );
  pkv_Selectf ( ncurves, spdimen*(degree+1), pitch, spdimen*(degree+1),
                &ctlpoints[(lastknot-2*degree-1)*spdimen], e );
  for ( k = 0; k < degree; k++ ) {
    memcpy ( d, e, ncurves*spdimen*(degree+1)*sizeof(float) );
    for ( j = 0; j < degree; j++ )
      for ( i = degree-1; i >= j; i-- ) {
        alpha = (t[k+j]-knots[lastknot-2*degree+i])/
                (knots[lastknot-degree+i-j]-knots[lastknot-2*degree+i]);
        pkn_MatrixLinCombf ( ncurves, spdimen,
                             spdimen*(degree+1), &d[i*spdimen], 1.0-alpha,
                             spdimen*(degree+1), &d[(i+1)*spdimen], alpha,
                             spdimen*(degree+1), &d[(i+1)*spdimen] );
      }
    pkv_Selectf ( ncurves, spdimen, spdimen*(degree+1), pitch,
                  &d[degree*spdimen], &ctlpoints[(lastknot-2*degree+k)*spdimen] );
  }
  memcpy ( &knots[lastknot-degree], newknots, (degree+1)*sizeof(float) );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_multiBSChangeRightKnotsf*/

