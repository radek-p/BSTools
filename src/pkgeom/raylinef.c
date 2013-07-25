
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"


/* find the pair of closest points on a ray and on a line segment */
float pkg_LineRayDist3f ( ray3f *ray, point3f *q0, point3f *q1,
                          float *rparam, float *lparam,
                          point3f *rp, point3f *lp )
{
  float    dist, a[6], x[2];
  vector3f w, qp;
  point3f  _p, _q;

  SubtractPoints3f ( q0, &ray->p, &qp );
  SubtractPoints3f ( q0, q1, &w );
  pkv_Selectf ( 3, 1, 1, 2, &ray->v.x, a );
  pkv_Selectf ( 3, 1, 1, 2, &w.x, &a[1] );
  if ( !rp ) rp = &_p;
  if ( !lp ) lp = &_q;
  if ( pkn_multiSolveRLSQf ( 3, 2, a, 1, 1, &qp.x, 1, x ) ) {
    AddVector3Mf ( &ray->p, &ray->v, x[0], rp );
    AddVector3Mf ( q0, &w, -x[1], lp );
  }
  else {  /* the line segment is parallel to the ray */
    x[0] = DotProduct3f ( &ray->v, &qp ) / DotProduct3f ( &ray->v, &ray->v );
    x[1] = 0.0;
    AddVector3Mf ( &ray->p, &ray->v, x[0], rp );
    *lp = *q0;
  } 
  if ( rparam ) *rparam = x[0];
  if ( lparam ) *lparam = x[1];
  SubtractPoints3f ( lp, rp, &w );
  dist = sqrt ( DotProduct3f ( &w, &w ) );
  return dist;
} /*pkg_LineRayDist3f*/

