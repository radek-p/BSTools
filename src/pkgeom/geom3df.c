
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2011                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <math.h>

#include "pknum.h"
#include "pkgeom.h"


void TransPoint3df ( const trans3d *tr, const point3f *p, point3f *q )
{
  /* q := tr * p */
  point3f r;

  r.x = tr->U0.a11 * p->x + tr->U0.a12 * p->y + tr->U0.a13 * p->z + tr->U0.a14;
  r.y = tr->U0.a21 * p->x + tr->U0.a22 * p->y + tr->U0.a23 * p->z + tr->U0.a24;
  r.z = tr->U0.a31 * p->x + tr->U0.a32 * p->y + tr->U0.a33 * p->z + tr->U0.a34;
  *q = r;
} /*TransPoint3df*/

void TransVector3df ( const trans3d *tr, const vector3f *v, vector3f *w )
{
  /* w := tr * v */
  vector3f r;

  r.x = tr->U0.a11 * v->x + tr->U0.a12 * v->y + tr->U0.a13 * v->z;
  r.y = tr->U0.a21 * v->x + tr->U0.a22 * v->y + tr->U0.a23 * v->z;
  r.z = tr->U0.a31 * v->x + tr->U0.a32 * v->y + tr->U0.a33 * v->z;
  *w = r;
} /*TransVector3df*/

void TransContra3df ( const trans3d *tri, const vector3f *v, vector3f *w )
{
  /* w := v * tri; v is a contravariant vector,                */
  /* tri is an inversion of the transformation being evaluated */
  vector3f r;

  r.x = tri->U0.a11 * v->x + tri->U0.a21 * v->y + tri->U0.a31 * v->z;
  r.y = tri->U0.a12 * v->x + tri->U0.a22 * v->y + tri->U0.a32 * v->z;
  r.z = tri->U0.a13 * v->x + tri->U0.a23 * v->y + tri->U0.a33 * v->z;
  *w = r;
} /*TransContra3df*/

void Trans3Point2df ( const trans3d *tr, const point2f *p, point2f *q )
{
  point2f r;

  r.x = tr->U0.a11*p->x + tr->U0.a12*p->y + tr->U0.a14;
  r.y = tr->U0.a21*p->x + tr->U0.a22*p->y + tr->U0.a24;
  *q = r;
} /*Trans3Point2df*/

void Trans3dTof ( const trans3d *trd, trans3f *trf )
{
  pkv_Selectdf ( 3, 4, &trd->U1.a[1][0]-&trd->U1.a[0][0],
                 &trf->U1.a[1][0]-&trf->U1.a[0][0],
                 &trd->U1.a[0][0], &trf->U1.a[0][0] );      
  trf->U1.detsgn = trd->U1.detsgn;
} /*Trans3fTod*/

