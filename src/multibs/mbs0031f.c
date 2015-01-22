
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

boolean mbs_BCFrenetC2Rf ( int degree, const point3f *ctlpoints, float t,
                           point2f *cpoint, vector2f *fframe, float *curvature )
{
  vector3f c[3];
  vector2f dc;
  float   num, denom;

  if ( degree < 2 ) {
    *curvature = 0.0;
    if ( degree == 1 ) {
      SetVector2f ( &fframe[0],
          ctlpoints[1].x*ctlpoints[0].z-ctlpoints[0].x*ctlpoints[1].z,
          ctlpoints[1].y*ctlpoints[0].z-ctlpoints[0].y*ctlpoints[1].z );
      NormalizeVector2f ( &fframe[0] );
      SetVector2f ( &fframe[1], -fframe[0].y, fframe[0].x );
      InterPoint3f ( &ctlpoints[0], &ctlpoints[1], t, &c[0] );
      SetPoint2f ( cpoint, c[0].x/c[0].z, c[0].y/c[0].z );
    }
    else        /* for a curve of degree 0 the Frenet frame is undefined */
      SetPoint2f ( cpoint,
          ctlpoints[0].x/ctlpoints[0].z, ctlpoints[0].y/ctlpoints[0].z );
  }
  else {
                /* deal with the curve of degree at least 2 */
    if ( !mbs_multiBCHornerf ( degree-2, 3, 3, 3, (float*)ctlpoints, t, (float*)c ) )
      return false;
                /* the last 2 steps of the de Casteljau algorithm */
                /* and computation of the curvature */
    num = (float)det3f ( &c[0], &c[1], &c[2] );
    InterPoint3f ( &c[0], &c[1], t, &c[0] );
    InterPoint3f ( &c[1], &c[2], t, &c[1] );
    SetVector2f ( &dc,
        c[1].x*c[0].z-c[0].x*c[1].z, c[1].y*c[0].z-c[0].y*c[1].z );
    InterPoint3f ( &c[0], &c[1], t, &c[0] );
    SetPoint2f ( cpoint, c[0].x/c[0].z, c[0].y/c[0].z );
                /* compute the Frenet frame vectors */
    fframe[0] = dc;
    NormalizeVector2f ( &fframe[0] );
    SetVector2f ( &fframe[1], -fframe[0].y, fframe[0].x );
                /* compute the curvature */
    denom = (float)DotProduct2f ( &dc, &dc );
    denom *= (float)sqrt ( denom );
    *curvature = c[0].z*c[0].z*c[0].z*
                 (float)(degree-1)/(float)degree*num/denom;
  }
  return true;
} /*mbs_BCFrenetC2Rf*/

