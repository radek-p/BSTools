
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2011, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <math.h>

#include "pknum.h"
#include "pkgeom.h"


void TransPoint3fd ( const trans3f *tr, const point3d *p, point3d *q )
{
  /* q := tr * p */
  point3d r;

  r.x = tr->U0.a11 * p->x + tr->U0.a12 * p->y + tr->U0.a13 * p->z + tr->U0.a14;
  r.y = tr->U0.a21 * p->x + tr->U0.a22 * p->y + tr->U0.a23 * p->z + tr->U0.a24;
  r.z = tr->U0.a31 * p->x + tr->U0.a32 * p->y + tr->U0.a33 * p->z + tr->U0.a34;
  *q = r;
} /*TransPoint3fd*/

void TransVector3fd ( const trans3f *tr, const vector3d *v, vector3d *w )
{
  /* w := tr * v */
  vector3d r;

  r.x = tr->U0.a11 * v->x + tr->U0.a12 * v->y + tr->U0.a13 * v->z;
  r.y = tr->U0.a21 * v->x + tr->U0.a22 * v->y + tr->U0.a23 * v->z;
  r.z = tr->U0.a31 * v->x + tr->U0.a32 * v->y + tr->U0.a33 * v->z;
  *w = r;
} /*TransVector3fd*/

void TransContra3fd ( const trans3f *tri, const vector3d *v, vector3d *w )
{
  /* w := v * tri; v is a contravariant vector,                */
  /* tri is an inversion of the transformation being evaluated */
  vector3d r;

  r.x = tri->U0.a11 * v->x + tri->U0.a21 * v->y + tri->U0.a31 * v->z;
  r.y = tri->U0.a12 * v->x + tri->U0.a22 * v->y + tri->U0.a32 * v->z;
  r.z = tri->U0.a13 * v->x + tri->U0.a23 * v->y + tri->U0.a33 * v->z;
  *w = r;
} /*TransContra3fd*/

void Trans3Point2fd ( const trans3f *tr, const point2d *p, point2d *q )
{
  point2d r;

  r.x = tr->U0.a11*p->x + tr->U0.a12*p->y + tr->U0.a14;
  r.y = tr->U0.a21*p->x + tr->U0.a22*p->y + tr->U0.a24;
  *q = r;
} /*Trans3Point2fd*/

void Trans3fTod ( const trans3f *trf, trans3d *trd )
{
  pkv_Selectfd ( 3, 4, &trf->U1.a[1][0]-&trf->U1.a[0][0],
                 &trd->U1.a[1][0]-&trd->U1.a[0][0],
                 &trf->U1.a[0][0], &trd->U1.a[0][0] );
  trd->U1.detsgn = trf->U1.detsgn;
} /*Trans3fTod*/

