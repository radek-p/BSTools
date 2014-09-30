
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "pkgeomclip.h"

static boolean _LB_Testd ( double p, double q, double *t0, double *t1 )
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
  else if ( q < 0.0 )
    return false;
  return true;
} /*_LB_Testd*/

int LiangBarskyClip2d ( point2d *p0, point2d *p1, double t0, double t1,
                        Box2d *box,
                        point2d *q0, point2d *q1 )
{
  double d, u0, u1;

  u0 = t0;
  u1 = t1;
  d = p1->x-p0->x;
  if ( _LB_Testd ( -d, p0->x-box->x0, &u0, &u1 ) )
    if ( _LB_Testd ( d, box->x1-p0->x, &u0, &u1 ) ) {
      d = p1->y-p0->y;
      if ( _LB_Testd ( -d, p0->y-box->y0, &u0, &u1 ) )
        if ( _LB_Testd ( d, box->y1-p0->y, &u0, &u1 ) ) {
          if ( q0 )
            InterPoint2d ( p0, p1, u0, q0 );
          if ( q1 )
            InterPoint2d ( p0, p1, u1, q1 );
          if ( u0 > t0 || u1 < t1 )
            return PKGEOM_CLIP_PART;
          else
            return PKGEOM_CLIP_ENTIRE;
        }
    }
  return PKGEOM_CLIP_NONE;
} /*LiangBarskyClipd*/

