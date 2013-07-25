
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>  
#include <memory.h>
#include <pthread.h>

#include "pkvaria.h"
#include "pkgeom.h" 
#include "multibs.h"
#include "raybez.h"

/* ///////////////////////////////////////////////////////////////////////// */
static char ClipTestd ( double p, double q, double *t0, double *t1 )
{
#define EPS 1.0e-10
  double r;

  if ( p < 0.0 ) {
    r = q/p;
    if ( r > *t1+EPS ) return 0;
    else if ( r > *t0 ) *t0 = r;
  }
  else if ( p > 0.0 ) {
    r = q/p;
    if ( r < *t0-EPS ) return 0;
    else if ( r < *t1 ) *t1 = r;
  }
  else if ( q < -EPS ) return 0;
  return 1;
#undef EPS
} /*ClipTestd*/

char rbez_TestRayBBoxd ( ray3d *ray, Box3d *box )
{
  double t0, t1;

  t0 = 0.0;  t1 = 1.0e308;
  if ( ClipTestd ( -ray->v.x, ray->p.x - box->x0, &t0, &t1 ) )
    if ( ClipTestd (  ray->v.x, box->x1 - ray->p.x, &t0, &t1 ) )
      if ( ClipTestd ( -ray->v.y, ray->p.y - box->y0, &t0, &t1 ) )
        if ( ClipTestd (  ray->v.y, box->y1 - ray->p.y, &t0, &t1 ) )
          if ( ClipTestd ( -ray->v.z, ray->p.z - box->z0, &t0, &t1 ) )
            if ( ClipTestd (  ray->v.z, box->z1 - ray->p.z, &t0, &t1 ) )
              return 1;
  return 0;
} /*rbez_TestRayBBoxd*/

