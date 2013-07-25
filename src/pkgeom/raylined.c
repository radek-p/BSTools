
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
double pkg_LineRayDist3d ( ray3d *ray, point3d *q0, point3d *q1,
                           double *rparam, double *lparam,
                           point3d *rp, point3d *lp )
{
  double   dist, a[6], x[2];
  vector3d w, qp;
  point3d  _p, _q;

  SubtractPoints3d ( q0, &ray->p, &qp );
  SubtractPoints3d ( q0, q1, &w );
  pkv_Selectd ( 3, 1, 1, 2, &ray->v.x, a );
  pkv_Selectd ( 3, 1, 1, 2, &w.x, &a[1] );
  if ( !rp ) rp = &_p;
  if ( !lp ) lp = &_q;
  if ( pkn_multiSolveRLSQd ( 3, 2, a, 1, 1, &qp.x, 1, x ) ) {
    AddVector3Md ( &ray->p, &ray->v, x[0], rp );
    AddVector3Md ( q0, &w, -x[1], lp );
  }
  else {  /* the line segment is parallel to the ray */
    x[0] = DotProduct3d ( &ray->v, &qp ) / DotProduct3d ( &ray->v, &ray->v );
    x[1] = 0.0;
    AddVector3Md ( &ray->p, &ray->v, x[0], rp );
    *lp = *q0;
  } 
  if ( rparam ) *rparam = x[0];
  if ( lparam ) *lparam = x[1];
  SubtractPoints3d ( lp, rp, &w );
  dist = sqrt ( DotProduct3d ( &w, &w ) );
  return dist;
} /*pkg_LineRayDist3d*/

