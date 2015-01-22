
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

boolean mbs_BCFrenetC3Rd ( int degree, const point4d *ctlpoints, double t,
                           point3d *cpoint, vector3d *fframe, double *curvatures )
{
  vector4d c[4];
  vector3d dc1, dc2;
  double    numt, denomt, denomc, ww, www;
  int      i;

  switch ( degree ) {
case 0:
    Point4to3d ( ctlpoints, cpoint );
    return true;

case 1:
    InterPoint4d ( &ctlpoints[0], &ctlpoints[1], t, &c[0] );
    SetVector3d ( &fframe[0],
        ctlpoints[1].x*ctlpoints[0].w-ctlpoints[0].x*ctlpoints[1].w,
        ctlpoints[1].y*ctlpoints[0].w-ctlpoints[0].y*ctlpoints[1].w,
        ctlpoints[1].z*ctlpoints[0].w-ctlpoints[0].z*ctlpoints[1].w );
    NormalizeVector3d ( &fframe[0] );
    curvatures[0] = curvatures[1] = 0.0;
    return true;

case 2:
    memcpy ( c, ctlpoints, 3*sizeof(point4d) );
    numt = 0.0;
    goto tail;

default:
    if ( !mbs_multiBCHornerd ( degree-3, 4, 4, 4, (double*)ctlpoints, t, (double*)c ) )
      return false;
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
    return true;
  }
} /*mbs_BCFrenetC3Rd*/

