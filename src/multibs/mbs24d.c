
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
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

#include "msgpool.h"

/* ////////////////////////////////////////// */
/* change of n+1 leftmost knots in a B-spline representation of curves, */
/* using the Oslo algorithm */

void mbs_multiBSChangeLeftKnotsd ( int ncurves, int spdimen, int degree,
                                   double *knots, int pitch, double *ctlpoints,
                                   double *newknots )
{
  void        *sp;
  int         i, j, k;
  double      *d, *e, *t;
  long double alpha;

  sp = pkv_GetScratchMemTop ();
  t = pkv_GetScratchMemd ( 2*degree );
  d = pkv_GetScratchMemd ( ncurves*spdimen*(degree+1) );
  e = pkv_GetScratchMemd ( ncurves*spdimen*(degree+1) );
  if ( !t || !d || !e )
    pkv_SignalError ( LIB_MULTIBS, 55, ERRMSG_0 );
  memcpy ( t, newknots, (degree+1)*sizeof(double) );
  memcpy ( &t[degree+1], &knots[degree+1], (degree-1)*sizeof(double) );
  pkv_Selectd ( ncurves, spdimen*(degree+1), pitch, spdimen*(degree+1),
                ctlpoints, e );
  for ( k = 0; k < degree; k++ ) {
    memcpy ( d, e, ncurves*spdimen*(degree+1)*sizeof(double) );
    for ( j = 1; j <= degree; j++ )
      for ( i = degree; i >= j; i-- ) {
        alpha = (t[k+j]-knots[i])/
                (knots[i+degree+1-j]-knots[i]);
        pkn_MatrixLinCombd ( ncurves, spdimen,
                   spdimen*(degree+1), &d[(i-1)*spdimen], (double)(1.0-alpha),
                   spdimen*(degree+1), &d[i*spdimen], (double)alpha,
                   spdimen*(degree+1), &d[i*spdimen] );
      }
    pkv_Selectd ( ncurves, spdimen, spdimen*(degree+1), pitch,
                  &d[degree*spdimen], &ctlpoints[k*spdimen] );
  }
  memcpy ( knots, newknots, (degree+1)*sizeof(double) );
  pkv_SetScratchMemTop ( sp );
} /*mbs_multiBSChangeLeftKnotsd*/

/* ////////////////////////////////////////// */
/* change of n+1 rightmost knots in a B-spline representation of curves */

void mbs_multiBSChangeRightKnotsd ( int ncurves, int spdimen, int degree,
                                    int lastknot, double *knots,
                                    int pitch, double *ctlpoints,
                                    double *newknots )
{
  void        *sp;
  int         i, j, k;
  double      *d, *e, *t;
  long double alpha;

  sp = pkv_GetScratchMemTop ();
  t = pkv_GetScratchMemd ( 2*degree );
  d = pkv_GetScratchMemd ( ncurves*spdimen*(degree+1) );
  e = pkv_GetScratchMemd ( ncurves*spdimen*(degree+1) );
  if ( !t || !d || !e )
    pkv_SignalError ( LIB_MULTIBS, 56, ERRMSG_0 );
  memcpy ( t, &knots[lastknot-2*degree+1], (degree-1)*sizeof(double) );
  memcpy ( &t[degree-1], newknots, (degree+1)*sizeof(double) );
  pkv_Selectd ( ncurves, spdimen*(degree+1), pitch, spdimen*(degree+1),
                &ctlpoints[(lastknot-2*degree-1)*spdimen], e );
  for ( k = 0; k < degree; k++ ) {
    memcpy ( d, e, ncurves*spdimen*(degree+1)*sizeof(double) );
    for ( j = 0; j < degree; j++ )
      for ( i = degree-1; i >= j; i-- ) {
        alpha = (t[k+j]-knots[lastknot-2*degree+i])/
                (knots[lastknot-degree+i-j]-knots[lastknot-2*degree+i]);
        pkn_MatrixLinCombd ( ncurves, spdimen,
                   spdimen*(degree+1), &d[i*spdimen], (double)(1.0-alpha),
                   spdimen*(degree+1), &d[(i+1)*spdimen], (double)alpha,
                   spdimen*(degree+1), &d[(i+1)*spdimen] );
      }
    pkv_Selectd ( ncurves, spdimen, spdimen*(degree+1), pitch,
                  &d[degree*spdimen], &ctlpoints[(lastknot-2*degree+k)*spdimen] );
  }
  memcpy ( &knots[lastknot-degree], newknots, (degree+1)*sizeof(double) );
  pkv_SetScratchMemTop ( sp );
} /*mbs_multiBSChangeRightKnotsd*/

