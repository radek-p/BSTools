
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
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
/* there must be ncols > 1, nrows >= 1 */
void rbez_FindCPBoundingBox3Rd ( int nrows, int ncols, int pitch,
                                 point4d *cp, double eps, Box3d *bbox )
{
  int     i, j, k, ncp;
  point4d *q;
  point3d a, b;

  ncp = ncols*nrows;
  i = 0;
  q = cp;
  Point4to3d ( &q[0], &a );
  if ( ncp & 0x01 ) {
    bbox->x0 = bbox->x1 = a.x;
    bbox->y0 = bbox->y1 = a.y;
    bbox->z0 = bbox->z1 = a.z;
    j = k = 1;
  }
  else {
    Point4to3d ( &q[1], &b );
    pkv_Sort2d ( &a.x, &b.x );  bbox->x0 = a.x;  bbox->x1 = b.x;
    pkv_Sort2d ( &a.y, &b.y );  bbox->y0 = a.y;  bbox->y1 = b.y;
    pkv_Sort2d ( &a.z, &b.z );  bbox->z0 = a.z;  bbox->z1 = b.z;
    j = k = 2;
  }
  while ( k < ncp ) {
    if ( j >= ncols ) { i++;  q = &cp[i*pitch];  Point4to3d ( &q[0], &a );  j = 1; }
                 else Point4to3d ( &q[j++], &a );
    if ( j >= ncols ) { i++;  q = &cp[i*pitch];  Point4to3d ( &q[0], &b );  j = 1; }
                 else Point4to3d ( &q[j++], &b );
    pkv_Sort2d ( &a.x, &b.x );
    bbox->x0 = min ( bbox->x0, a.x );  bbox->x1 = max ( bbox->x1, b.x );
    pkv_Sort2d ( &a.y, &b.y );
    bbox->y0 = min ( bbox->y0, a.y );  bbox->y1 = max ( bbox->y1, b.y );
    pkv_Sort2d ( &a.z, &b.z );
    bbox->z0 = min ( bbox->z0, a.z );  bbox->z1 = max ( bbox->z1, b.z );
    k += 2;
  }
  bbox->x0 -= eps;  bbox->x1 += eps;
  bbox->y0 -= eps;  bbox->y1 += eps;
  bbox->z0 -= eps;  bbox->z1 += eps;
} /*rbez_FindCPBoundingBox3d*/

