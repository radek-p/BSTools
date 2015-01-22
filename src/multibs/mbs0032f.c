
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

boolean mbs_BCFrenetC3f ( int degree, const point3f *ctlpoints, float t,
                          point3f *cpoint, vector3f *fframe, float *curvatures )
{
  point3f  c[4];
  vector3f dc0, dc1, dc2;
  float   numt, denomt, denomc;
  int     i;

  switch ( degree ) {
case 0:
    *cpoint = ctlpoints[0];
    return true;

case 1:
    InterPoint3f ( &ctlpoints[0], &ctlpoints[1], t, cpoint );
    SubtractPoints3f ( &ctlpoints[1], &ctlpoints[0], &fframe[0] );
    NormalizeVector3f ( &fframe[0] );
    curvatures[0] = curvatures[1] = 0.0;
    return true;

case 2:
    memcpy ( c, ctlpoints, 3*sizeof(point3f) );
    numt = 0.0;
    goto tail;

default:
    if ( !mbs_multiBCHornerf ( degree-3, 4, 3, 3, (float*)ctlpoints, t, (float*)c ) )
      return false;
    SubtractPoints3f ( &c[1], &c[0], &dc0 );
    SubtractPoints3f ( &c[2], &c[1], &dc1 );
    SubtractPoints3f ( &c[3], &c[2], &dc2 );
    numt = (float)det3f ( &dc0, &dc1, &dc2 );
    for ( i = 0; i < 3; i++ )
      InterPoint3f ( &c[i], &c[i+1], t, &c[i] );

tail:
    SubtractPoints3f ( &c[1], &c[0], &dc0 );
    SubtractPoints3f ( &c[2], &c[1], &dc1 );
    CrossProduct3f ( &dc0, &dc1, &dc2 );
    denomt = (float)DotProduct3f ( &dc2, &dc2 );
    for ( i = 0; i < 2; i++ )
      InterPoint3f ( &c[i], &c[i+1], t, &c[i] );
    SubtractPoints3f ( &c[1], &c[0], &dc1 );
    denomc = (float)DotProduct3f ( &dc1, &dc1 );
    denomc *= (float)sqrt ( denomc );
    InterPoint3f ( &c[0], &c[1], t, cpoint );

    curvatures[0] = (float)((float)(degree-1)/(float)degree*sqrt(denomt)/denomc);
    curvatures[1] = (float)((float)(degree-2)/(float)degree*numt/denomt);
    fframe[0] = dc1;
    NormalizeVector3f ( &fframe[0] );
    fframe[2] = dc2;
    NormalizeVector3f ( &fframe[2] );
    CrossProduct3f ( &fframe[2], &fframe[0], &fframe[1] );
    return true;
  }
} /*mbs_BCFrenetC3f*/

