
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2014                            */
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

boolean rbez_TestRayBBoxd ( ray3d *ray, Box3d *box )
{
  double t0, t1;

  t0 = 0.0;  t1 = 1.0e308;
  if ( ClipTestd ( -ray->v.x, ray->p.x - box->x0, &t0, &t1 ) )
    if ( ClipTestd (  ray->v.x, box->x1 - ray->p.x, &t0, &t1 ) )
      if ( ClipTestd ( -ray->v.y, ray->p.y - box->y0, &t0, &t1 ) )
        if ( ClipTestd (  ray->v.y, box->y1 - ray->p.y, &t0, &t1 ) )
          if ( ClipTestd ( -ray->v.z, ray->p.z - box->z0, &t0, &t1 ) )
            if ( ClipTestd (  ray->v.z, box->z1 - ray->p.z, &t0, &t1 ) )
              return true;
  return false;
} /*rbez_TestRayBBoxd*/

void rbez_FindSumBBoxd ( Box3d *box1, Box3d *box2, Box3d *box )
{
  box->x0 = min ( box1->x0, box2->x0 );
  box->x1 = max ( box1->x1, box2->x1 );
  box->y0 = min ( box1->y0, box2->y0 );
  box->y1 = max ( box1->y1, box2->y1 );
  box->z0 = min ( box1->z0, box2->z0 );
  box->z1 = max ( box1->z1, box2->z1 );
} /*rbez_FindSumBBoxd*/ 

boolean rbez_NarrowBBoxSumd ( Box3d *box1, Box3d *box2, Box3d *box )
{
  Box3d   sum;
  boolean change;

  rbez_FindSumBBoxd ( box1, box2, &sum );
  change = false;
  if ( sum.x0 > box->x0 ) { box->x0 = sum.x0;  change = true; }
  if ( sum.x1 < box->x1 ) { box->x1 = sum.x1;  change = true; }
  if ( sum.y0 > box->y0 ) { box->y0 = sum.y0;  change = true; }
  if ( sum.y1 < box->y1 ) { box->y1 = sum.y1;  change = true; }
  if ( sum.z0 > box->z0 ) { box->z0 = sum.z0;  change = true; }
  if ( sum.z1 < box->z1 ) { box->z1 = sum.z1;  change = true; }
  return change;
} /*rbez_NarrowBBoxSumd*/

