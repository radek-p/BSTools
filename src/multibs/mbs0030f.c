
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

boolean mbs_BCFrenetC2f ( int degree, const point2f *ctlpoints, float t,
                          point2f *cpoint, vector2f *fframe, float *curvature )
{
  point2f  c[3];
  vector2f dc0, dc1, dc; 
  float   num, denom;

  if ( degree < 2 ) {
    *curvature = 0.0;
    if ( degree == 1 ) {
      SubtractPoints2f ( &ctlpoints[1], &ctlpoints[0], &fframe[0] );
      NormalizeVector2f ( &fframe[0] );
      SetVector2f ( &fframe[1], -fframe[0].y, fframe[0].x );
      InterPoint2f ( &ctlpoints[0], &ctlpoints[1], t, cpoint );
    }
    else        /* for a curve of degree 0 the Frenet frame is undefined */
      *cpoint = ctlpoints[0];
  }
  else {
                /* deal with the curve of degree at least 2 */
    if ( !mbs_multiBCHornerf ( degree-2, 3, 2, 2, (float*)ctlpoints, t, (float*)c ) )
      return false;
                /* the last 2 steps of the de Casteljau algorithm */
                /* and computation of differences */
    SubtractPoints2f ( &c[1], &c[0], &dc0 );
    SubtractPoints2f ( &c[2], &c[1], &dc1 );
    InterPoint2f ( &c[0], &c[1], t, &c[0] );
    InterPoint2f ( &c[1], &c[2], t, &c[1] );
    SubtractPoints2f ( &c[1], &c[0], &dc );
    InterPoint2f ( &c[0], &c[1], t, cpoint );
                /* compute the Frenet frame vectors */
    fframe[0] = dc;
    NormalizeVector2f ( &fframe[0] );
    SetVector2f ( &fframe[1], -fframe[0].y, fframe[0].x );
                /* compute the curvature */
    num   = (float)det2f ( &dc0, &dc1 );
    denom = (float)DotProduct2f ( &dc, &dc );
    denom *= (float)sqrt(denom);
    *curvature = (float)(degree-1)/(float)degree*num/denom;
  }
  return true;
} /*mbs_BCFrenetC2f*/

