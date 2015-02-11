
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014, 2015                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"
#include "raybez.h"

#include "raybezprivated.h"

/* ///////////////////////////////////////////////////////////////////////// */
void rbez_InitBBox3d ( Box3d *bbox, point3d *p )
{
  bbox->x0 = bbox->x1 = p->x;
  bbox->y0 = bbox->y1 = p->y;
  bbox->z0 = bbox->z1 = p->z;
} /*rbez_InitBBox3d*/

void rbez_ExtendBBox3d ( Box3d *bbox, point3d *p )
{
  if ( p->x < bbox->x0 )      bbox->x0 = p->x;
  else if ( p->x > bbox->x1 ) bbox->x1 = p->x;
  if ( p->y < bbox->y0 )      bbox->y0 = p->y;
  else if ( p->y > bbox->y1 ) bbox->y1 = p->y;
  if ( p->z < bbox->z0 )      bbox->z0 = p->z;
  else if ( p->z > bbox->z1 ) bbox->z1 = p->z;
} /*rbez_ExtendBBox3d*/

void rbez_Extend2BBox3d ( Box3d *bbox, point3d *p1, point3d *p2 )
{
  if ( p1->x < p2->x ) {
    if ( p1->x < bbox->x0 ) bbox->x0 = p1->x;
    if ( p2->x > bbox->x1 ) bbox->x1 = p2->x;
  }
  else {
    if ( p2->x < bbox->x0 ) bbox->x0 = p2->x;
    if ( p1->x > bbox->x1 ) bbox->x1 = p1->x;
  }
  if ( p1->y < p2->y ) {
    if ( p1->y < bbox->y0 ) bbox->y0 = p1->y;
    if ( p2->y > bbox->y1 ) bbox->y1 = p2->y;
  }
  else {
    if ( p2->y < bbox->y0 ) bbox->y0 = p2->y;
    if ( p1->y > bbox->y1 ) bbox->y1 = p1->y;
  }
  if ( p1->z < p2->z ) {
    if ( p1->z < bbox->z0 ) bbox->z0 = p1->z;
    if ( p2->z > bbox->z1 ) bbox->z1 = p2->z;
  }
  else {
    if ( p2->z < bbox->z0 ) bbox->z0 = p2->z;
    if ( p1->z > bbox->z1 ) bbox->z1 = p1->z;
  }
} /*rbez_Extend2BBox3d*/

/* there must be ncols > 1, nrows >= 1 */
void rbez_FindCPBoundingBox3d ( int nrows, int ncols, int pitch,
                                point3d *cp, double eps, Box3d *bbox )
{
  int     i, j, k, ncp;
  point3d *q, *a, *b;

  ncp = ncols*nrows;
  i = 0;
  q = cp;
  rbez_InitBBox3d ( bbox, &q[0] );
  if ( ncp & 0x01 )
    j = k = 1;
  else {
    rbez_ExtendBBox3d ( bbox, &q[1] );
    j = k = 2;
  }
  while ( k < ncp ) {
    if ( j >= ncols ) { i++;  q = &cp[i*pitch];  a = q;  j = 1; }
                 else a = &q[j++];
    if ( j >= ncols ) { i++;  q = &cp[i*pitch];  b = q;  j = 1; }
                 else b = &q[j++];
    rbez_Extend2BBox3d ( bbox, a, b );
    k += 2;
  }
  bbox->x0 -= eps;  bbox->x1 += eps;
  bbox->y0 -= eps;  bbox->y1 += eps;
  bbox->z0 -= eps;  bbox->z1 += eps;
} /*rbez_FindCPBoundingBox3d*/

