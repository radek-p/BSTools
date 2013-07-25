
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

boolean CameraClipPoint3d ( CameraRecd *CPos, point3d *p, point3d *q )
{
  int      i;
  vector4d *cpl;

  for ( i = 0; i < CPos->ncplanes; i++ ) {
    cpl = &CPos->cplane[i];
    if ( DotProduct3d ( (vector3d*)cpl, p ) + cpl->w < 0.0 )
      return false;
  }
  CameraProjectPoint3d ( CPos, p, q );
  return true;
} /*CameraClipPoint3d*/

/* ///////////////////////////////////////////////////////////////////////// */
/* An adaptation of the Liang & Barsky line clipping algorithm               */

static boolean _LBTestd ( double p, double q, double *t0, double *t1 )
{
  double r;

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
} /*_LBTestd*/

boolean CameraClipLine3d ( CameraRecd *CPos,
                           point3d *p0, double t0, point3d *p1, double t1,
                           point3d *q0, point3d *q1 )
{
  int      i;
  double   p, q;
  vector3d v;
  vector4d *cpl;

  SubtractPoints3d ( p1, p0, &v );
  for ( i = 0; i < CPos->ncplanes; i++ ) {
    cpl = &CPos->cplane[i];
    p = -DotProduct3d ( (vector3d*)cpl, &v );
    q = DotProduct3d ( (vector3d*)cpl, p0 ) + cpl->w;
    if ( !_LBTestd ( p, q, &t0, &t1 ) )
      return false;
  }
  InterPoint3d ( p0, p1, t0, &v );
  CameraProjectPoint3d ( CPos, &v, q0 );
  InterPoint3d ( p0, p1, t1, &v );
  CameraProjectPoint3d ( CPos, &v, q1 );
  return true;
} /*CameraClipLine3d*/

