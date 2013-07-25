
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

void mbs_BCFrenetC2f ( int degree, const point2f *ctlpoints, float t,
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
    mbs_multiBCHornerf ( degree-2, 3, 2, 2, (float*)ctlpoints, t, (float*)c );
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
} /*mbs_BCFrenetC2f*/

void mbs_BCFrenetC2Rf ( int degree, const point3f *ctlpoints, float t,
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
    mbs_multiBCHornerf ( degree-2, 3, 3, 3, (float*)ctlpoints, t, (float*)c );
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
} /*mbs_BCFrenetC2Rf*/

void mbs_BCFrenetC3f ( int degree, const point3f *ctlpoints, float t,
                       point3f *cpoint, vector3f *fframe, float *curvatures )
{
  point3f  c[4];
  vector3f dc0, dc1, dc2;
  float   numt, denomt, denomc;
  int     i;

  switch ( degree ) {
case 0:
    *cpoint = ctlpoints[0];
    return;

case 1:
    InterPoint3f ( &ctlpoints[0], &ctlpoints[1], t, cpoint );
    SubtractPoints3f ( &ctlpoints[1], &ctlpoints[0], &fframe[0] );
    NormalizeVector3f ( &fframe[0] );
    curvatures[0] = curvatures[1] = 0.0;
    return;

case 2:
    memcpy ( c, ctlpoints, 3*sizeof(point3f) );
    numt = 0.0;
    goto tail;

default:
    mbs_multiBCHornerf ( degree-3, 4, 3, 3, (float*)ctlpoints, t, (float*)c );
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
    return;
  }
} /*mbs_BCFrenetC3f*/

void mbs_BCFrenetC3Rf ( int degree, const point4f *ctlpoints, float t,
                        point3f *cpoint, vector3f *fframe, float *curvatures )
{
  vector4f c[4];
  vector3f dc1, dc2;
  float    numt, denomt, denomc, ww, www;
  int      i;

  switch ( degree ) {
case 0:
    Point4to3f ( ctlpoints, cpoint );
    return;

case 1:
    InterPoint4f ( &ctlpoints[0], &ctlpoints[1], t, &c[0] );
    SetVector3f ( &fframe[0],
        ctlpoints[1].x*ctlpoints[0].w-ctlpoints[0].x*ctlpoints[1].w,
        ctlpoints[1].y*ctlpoints[0].w-ctlpoints[0].y*ctlpoints[1].w,
        ctlpoints[1].z*ctlpoints[0].w-ctlpoints[0].z*ctlpoints[1].w );
    NormalizeVector3f ( &fframe[0] );
    curvatures[0] = curvatures[1] = 0.0;
    return;

case 2:
    memcpy ( c, ctlpoints, 3*sizeof(point4f) );
    numt = 0.0;
    goto tail;

default:
    mbs_multiBCHornerf ( degree-3, 4, 4, 4, (float*)ctlpoints, t, (float*)c );
    numt = (float)det4f ( &c[0], &c[1], &c[2], &c[3] );
    for ( i = 0; i < 3; i++ )
      InterPoint4f ( &c[i], &c[i+1], t, &c[i] );

tail:
    CrossProduct4P3f ( &c[2], &c[1], &c[0], &dc2 );
    denomt = (float)DotProduct3f ( &dc2, &dc2 );
    for ( i = 0; i < 2; i++ )
      InterPoint4f ( &c[i], &c[i+1], t, &c[i] );
    OutProduct4P3f ( &c[1], &c[0], &dc1 );
    denomc = (float)DotProduct3f ( &dc1, &dc1 );
    denomc *= (float)sqrt ( denomc );
    InterPoint4f ( &c[0], &c[1], t, &c[0] );
    Point4to3f ( c, cpoint );
    www = c[0].w;  ww = www*www;  www *= ww;

    curvatures[0] = (float)((float)(degree-1)/(float)degree*sqrt(denomt)/denomc);
    curvatures[1] = (float)((float)(degree-2)/(float)degree*numt/denomt);
    fframe[0] = dc1;
    NormalizeVector3f ( &fframe[0] );
    fframe[2] = dc2;
    NormalizeVector3f ( &fframe[2] );
    CrossProduct3f ( &fframe[2], &fframe[0], &fframe[1] );
    return;
  }
} /*mbs_BCFrenetC3Rf*/

