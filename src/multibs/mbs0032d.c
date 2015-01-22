
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

boolean mbs_BCFrenetC3d ( int degree, const point3d *ctlpoints, double t,
                          point3d *cpoint, vector3d *fframe, double *curvatures )
{
  point3d  c[4];
  vector3d dc0, dc1, dc2;
  double   numt, denomt, denomc;
  int     i;

  switch ( degree ) {
case 0:
    *cpoint = ctlpoints[0];
    return true;

case 1:
    InterPoint3d ( &ctlpoints[0], &ctlpoints[1], t, cpoint );
    SubtractPoints3d ( &ctlpoints[1], &ctlpoints[0], &fframe[0] );
    NormalizeVector3d ( &fframe[0] );
    curvatures[0] = curvatures[1] = 0.0;
    return true;

case 2:
    memcpy ( c, ctlpoints, 3*sizeof(point3d) );
    numt = 0.0;
    goto tail;

default:
    if ( !mbs_multiBCHornerd ( degree-3, 4, 3, 3, (double*)ctlpoints, t, (double*)c ) )
      return false;
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
    return true;
  }
} /*mbs_BCFrenetC3d*/

