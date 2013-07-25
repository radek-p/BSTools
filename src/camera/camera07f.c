
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2008                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <math.h>

#include "camera.h"

#include "msgpool.h"


/* ///////////////////////////////////////////////////////////////////////// */
/* Clipping points and lines. On input points are given in the global        */
/* coordinates. The clipping procedures output projections, in the same      */
/* form, that the CameraProjectPoint3 procedure does.                        */

boolean CameraClipPoint3f ( CameraRecf *CPos, point3f *p, point3f *q )
{
  int      i;
  vector4f *cpl;

  for ( i = 0; i < CPos->ncplanes; i++ ) {
    cpl = &CPos->cplane[i];
    if ( DotProduct3f ( (vector3f*)cpl, p ) + cpl->w < 0.0 )
      return false;
  }
  CameraProjectPoint3f ( CPos, p, q );
  return true;
} /*CameraClipPoint3f*/

/* ///////////////////////////////////////////////////////////////////////// */
/* An adaptation of the Liang & Barsky line clipping algorithm               */

static boolean _LBTestf ( float p, float q, float *t0, float *t1 )
{
  float r;

  if ( p < 0.0 ) {
    r = q/p;
    if ( r > *t1 ) return false;
    else if ( r > *t0 ) *t0 = r;
  }
  else if ( p > 0.0 ) {
    r = q/p;
    if ( r < *t0 ) return false;
    else if ( r < *t1 ) *t1 = r;
  }
  else if ( q < 0.0 ) return false;
  return true;
} /*_LBTestf*/

boolean CameraClipLine3f ( CameraRecf *CPos,
                           point3f *p0, float t0, point3f *p1, float t1,
                           point3f *q0, point3f *q1 )
{
  int      i;
  float    p, q;
  vector3f v;
  vector4f *cpl;

  SubtractPoints3f ( p1, p0, &v );
  for ( i = 0; i < CPos->ncplanes; i++ ) {
    cpl = &CPos->cplane[i];
    p = (float)(-DotProduct3f ( (vector3f*)cpl, &v ));
    q = (float)DotProduct3f ( (vector3f*)cpl, p0 ) + cpl->w;
    if ( !_LBTestf ( p, q, &t0, &t1 ) )
      return false;
  }
  InterPoint3f ( p0, p1, t0, &v );
  CameraProjectPoint3f ( CPos, &v, q0 );
  InterPoint3f ( p0, p1, t1, &v );
  CameraProjectPoint3f ( CPos, &v, q1 );
  return true;
} /*CameraClipLine3f*/

