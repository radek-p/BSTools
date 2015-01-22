
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2015                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"


/* /////////////////////////////////////////// */
/* computing curvature and Frenet frames of Bezier curves */

boolean mbs_BCFrenetC2d ( int degree, const point2d *ctlpoints, double t,
                          point2d *cpoint, vector2d *fframe, double *curvature )
{
  point2d  c[3];
  vector2d dc0, dc1, dc; 
  double   num, denom;

  if ( degree < 2 ) {
    *curvature = 0.0;
    if ( degree == 1 ) {
      SubtractPoints2d ( &ctlpoints[1], &ctlpoints[0], &fframe[0] );
      NormalizeVector2d ( &fframe[0] );
      SetVector2d ( &fframe[1], -fframe[0].y, fframe[0].x );
      InterPoint2d ( &ctlpoints[0], &ctlpoints[1], t, cpoint );
    }
    else        /* for a curve of degree 0 the Frenet frame is undefined */
      *cpoint = ctlpoints[0];
  }
  else {
                /* deal with the curve of degree at least 2 */
    if ( !mbs_multiBCHornerd ( degree-2, 3, 2, 2, (double*)ctlpoints, t, (double*)c ) )
      return false;
                /* the last 2 steps of the de Casteljau algorithm */
                /* and computation of differences */
    SubtractPoints2d ( &c[1], &c[0], &dc0 );
    SubtractPoints2d ( &c[2], &c[1], &dc1 );
    InterPoint2d ( &c[0], &c[1], t, &c[0] );
    InterPoint2d ( &c[1], &c[2], t, &c[1] );
    SubtractPoints2d ( &c[1], &c[0], &dc );
    InterPoint2d ( &c[0], &c[1], t, cpoint );
                /* compute the Frenet frame vectors */
    fframe[0] = dc;
    NormalizeVector2d ( &fframe[0] );
    SetVector2d ( &fframe[1], -fframe[0].y, fframe[0].x );
                /* compute the curvature */
    num   = det2d ( &dc0, &dc1 );
    denom = DotProduct2d ( &dc, &dc );
    denom *= sqrt(denom);
    *curvature = (double)(degree-1)/(double)degree*num/denom;
  }
  return true;
} /*mbs_BCFrenetC2d*/

