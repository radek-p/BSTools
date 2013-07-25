
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

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"


/* /////////////////////////////////////////// */
/* computing curvature and Frenet frames of Bezier curves */

void mbs_BCFrenetC2d ( int degree, const point2d *ctlpoints, double t,
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
    mbs_multiBCHornerd ( degree-2, 3, 2, 2, (double*)ctlpoints, t, (double*)c );
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
} /*mbs_BCFrenetC2d*/

void mbs_BCFrenetC2Rd ( int degree, const point3d *ctlpoints, double t,
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
    mbs_multiBCHornerd ( degree-2, 3, 3, 3, (double*)ctlpoints, t, (double*)c );
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
} /*mbs_BCFrenetC2Rd*/

void mbs_BCFrenetC3d ( int degree, const point3d *ctlpoints, double t,
                       point3d *cpoint, vector3d *fframe, double *curvatures )
{
  point3d  c[4];
  vector3d dc0, dc1, dc2;
  double   numt, denomt, denomc;
  int     i;

  switch ( degree ) {
case 0:
    *cpoint = ctlpoints[0];
    return;

case 1:
    InterPoint3d ( &ctlpoints[0], &ctlpoints[1], t, cpoint );
    SubtractPoints3d ( &ctlpoints[1], &ctlpoints[0], &fframe[0] );
    NormalizeVector3d ( &fframe[0] );
    curvatures[0] = curvatures[1] = 0.0;
    return;

case 2:
    memcpy ( c, ctlpoints, 3*sizeof(point3d) );
    numt = 0.0;
    goto tail;

default:
    mbs_multiBCHornerd ( degree-3, 4, 3, 3, (double*)ctlpoints, t, (double*)c );
    SubtractPoints3d ( &c[1], &c[0], &dc0 );
    SubtractPoints3d ( &c[2], &c[1], &dc1 );
    SubtractPoints3d ( &c[3], &c[2], &dc2 );
    numt = det3d ( &dc0, &dc1, &dc2 );
    for ( i = 0; i < 3; i++ )
      InterPoint3d ( &c[i], &c[i+1], t, &c[i] );

tail:
    SubtractPoints3d ( &c[1], &c[0], &dc0 );
    SubtractPoints3d ( &c[2], &c[1], &dc1 );
    CrossProduct3d ( &dc0, &dc1, &dc2 );
    denomt = DotProduct3d ( &dc2, &dc2 );
    for ( i = 0; i < 2; i++ )
      InterPoint3d ( &c[i], &c[i+1], t, &c[i] );
    SubtractPoints3d ( &c[1], &c[0], &dc1 );
    denomc = DotProduct3d ( &dc1, &dc1 );
    denomc *= sqrt ( denomc );
    InterPoint3d ( &c[0], &c[1], t, cpoint );

    curvatures[0] = (double)(degree-1)/(double)degree*sqrt(denomt)/denomc;
    curvatures[1] = (double)(degree-2)/(double)degree*numt/denomt;
    fframe[0] = dc1;
    NormalizeVector3d ( &fframe[0] );
    fframe[2] = dc2;
    NormalizeVector3d ( &fframe[2] );
    CrossProduct3d ( &fframe[2], &fframe[0], &fframe[1] );
    return;
  }
} /*mbs_BCFrenetC3d*/

void mbs_BCFrenetC3Rd ( int degree, const point4d *ctlpoints, double t,
                        point3d *cpoint, vector3d *fframe, double *curvatures )
{
  vector4d c[4];
  vector3d dc1, dc2;
  double    numt, denomt, denomc, ww, www;
  int      i;

  switch ( degree ) {
case 0:
    Point4to3d ( ctlpoints, cpoint );
    return;

case 1:
    InterPoint4d ( &ctlpoints[0], &ctlpoints[1], t, &c[0] );
    SetVector3d ( &fframe[0],
        ctlpoints[1].x*ctlpoints[0].w-ctlpoints[0].x*ctlpoints[1].w,
        ctlpoints[1].y*ctlpoints[0].w-ctlpoints[0].y*ctlpoints[1].w,
        ctlpoints[1].z*ctlpoints[0].w-ctlpoints[0].z*ctlpoints[1].w );
    NormalizeVector3d ( &fframe[0] );
    curvatures[0] = curvatures[1] = 0.0;
    return;

case 2:
    memcpy ( c, ctlpoints, 3*sizeof(point4d) );
    numt = 0.0;
    goto tail;

default:
    mbs_multiBCHornerd ( degree-3, 4, 4, 4, (double*)ctlpoints, t, (double*)c );
    numt = det4d ( &c[0], &c[1], &c[2], &c[3] );
    for ( i = 0; i < 3; i++ )
      InterPoint4d ( &c[i], &c[i+1], t, &c[i] );

tail:
    CrossProduct4P3d ( &c[2], &c[1], &c[0], &dc2 );
    denomt = DotProduct3d ( &dc2, &dc2 );
    for ( i = 0; i < 2; i++ )
      InterPoint4d ( &c[i], &c[i+1], t, &c[i] );
    OutProduct4P3d ( &c[1], &c[0], &dc1 );
    denomc = DotProduct3d ( &dc1, &dc1 );
    denomc *= sqrt ( denomc );
    InterPoint4d ( &c[0], &c[1], t, &c[0] );
    Point4to3d ( c, cpoint );
    www = c[0].w;  ww = www*www;  www *= ww;

    curvatures[0] = www*(double)(degree-1)/(double)degree*sqrt(denomt)/denomc;
    curvatures[1] = -ww*(double)(degree-2)/(double)degree*numt/denomt;
    fframe[0] = dc1;
    NormalizeVector3d ( &fframe[0] );
    fframe[2] = dc2;
    NormalizeVector3d ( &fframe[2] );
    CrossProduct3d ( &fframe[2], &fframe[0], &fframe[1] );
    return;
  }
} /*mbs_BCFrenetC3Rd*/

