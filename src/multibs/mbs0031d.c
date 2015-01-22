
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

boolean mbs_BCFrenetC2Rd ( int degree, const point3d *ctlpoints, double t,
                           point2d *cpoint, vector2d *fframe, double *curvature )
{
  vector3d c[3];
  vector2d dc;
  double   num, denom;

  if ( degree < 2 ) {
    *curvature = 0.0;
    if ( degree == 1 ) {
      SetVector2d ( &fframe[0],
          ctlpoints[1].x*ctlpoints[0].z-ctlpoints[0].x*ctlpoints[1].z,
          ctlpoints[1].y*ctlpoints[0].z-ctlpoints[0].y*ctlpoints[1].z );
      NormalizeVector2d ( &fframe[0] );
      SetVector2d ( &fframe[1], -fframe[0].y, fframe[0].x );
      InterPoint3d ( &ctlpoints[0], &ctlpoints[1], t, &c[0] );
      SetPoint2d ( cpoint, c[0].x/c[0].z, c[0].y/c[0].z );
    }
    else        /* for a curve of degree 0 the Frenet frame is undefined */
      SetPoint2d ( cpoint,
          ctlpoints[0].x/ctlpoints[0].z, ctlpoints[0].y/ctlpoints[0].z );
  }
  else {
                /* deal with the curve of degree at least 2 */
    if ( !mbs_multiBCHornerd ( degree-2, 3, 3, 3, (double*)ctlpoints, t, (double*)c ) )
      return false;
                /* the last 2 steps of the de Casteljau algorithm */
                /* and computation of the curvature */
    num = det3d ( &c[0], &c[1], &c[2] );
    InterPoint3d ( &c[0], &c[1], t, &c[0] );
    InterPoint3d ( &c[1], &c[2], t, &c[1] );
    SetVector2d ( &dc,
        c[1].x*c[0].z-c[0].x*c[1].z, c[1].y*c[0].z-c[0].y*c[1].z );
    InterPoint3d ( &c[0], &c[1], t, &c[0] );
    SetPoint2d ( cpoint, c[0].x/c[0].z, c[0].y/c[0].z );
                /* compute the Frenet frame vectors */
    fframe[0] = dc;
    NormalizeVector2d ( &fframe[0] );
    SetVector2d ( &fframe[1], -fframe[0].y, fframe[0].x );
                /* compute the curvature */
    denom = DotProduct2d ( &dc, &dc );
    denom *= sqrt ( denom );
    *curvature = c[0].z*c[0].z*c[0].z*
                 (double)(degree-1)/(double)degree*num/denom;
  }
  return true;
} /*mbs_BCFrenetC2Rd*/

