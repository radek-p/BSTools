
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2011                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <math.h>

#include "pkgeom.h"

void SetPoint2d ( point2d *p, double x, double y )
{
  p->x = x;
  p->y = y;
} /*SetPoint2d*/

void TransPoint2d ( const trans2d *tr, const point2d *p, point2d *q )
{
  point2d r;

  r.x = tr->U0.a11 * p->x + tr->U0.a12 * p->y + tr->U0.a13;
  r.y = tr->U0.a21 * p->x + tr->U0.a22 * p->y + tr->U0.a23;
  *q = r;
} /*TransPoint2d*/

void TransVector2d ( const trans2d *tr, const vector2d *v, vector2d *w )
{
  vector2d r;

  r.x = tr->U0.a11 * v->x + tr->U0.a12 * v->y;
  r.y = tr->U0.a21 * v->x + tr->U0.a22 * v->y;
  *w = r;
} /*TransVector2d*/

void Trans2Point3d ( const trans2d *tr, const point3d *p, point3d *q )
{
  point3d r;

  r.x = tr->U0.a11*p->x + tr->U0.a12*p->y + tr->U0.a13*p->z;
  r.y = tr->U0.a21*p->x + tr->U0.a22*p->y + tr->U0.a23*p->z;
  r.z = p->z;
  *q = r;
} /*Trans2Point3d*/

void IdentTrans2d ( trans2d *tr )
{
  short i, j;

  for (i = 0; i <= 1; i++) {
    for (j = 0; j <= 2; j++)
      tr->U1.a[i][j] = 0.0;
    tr->U1.a[i][i] = 1.0;
  }
  tr->U1.detsgn = 1;
} /*IdentTrans2d*/

void CompTrans2d ( trans2d *s, trans2d *t, trans2d *u )
{
  byte i, j, k;

  for (i = 0; i <= 1; i++) {
    for (j = 0; j <= 2; j++) {
      s->U1.a[i][j] = 0.0;
      for (k = 0; k <= 1; k++)
	s->U1.a[i][j] += t->U1.a[i][k] * u->U1.a[k][j];
    }
    s->U1.a[i][2] += t->U1.a[i][2];
  }
  s->U1.detsgn = (short)(t->U1.detsgn * u->U1.detsgn);
} /*CompTrans2d*/

void ShiftTrans2d ( trans2d *tr, double tx, double ty )
{
  tr->U0.a13 += tx;
  tr->U0.a23 += ty;
} /*ShiftTrans2d*/

void RotTrans2d ( trans2d *tr, double angle )
{
  double c, s, m, n;
  byte i;

  if ( !angle )
    return;
  c = cos ( angle );
  s = sin ( angle );
  for (i = 0; i <= 2; i++) {
    m = tr->U1.a[0][i];
    n = tr->U1.a[1][i];
    tr->U1.a[0][i] = c * m - s * n;
    tr->U1.a[1][i] = s * m + c * n;
  }
} /*RotTrans2d*/

void ScaleTrans2d ( trans2d *t, double sx, double sy )
{
  byte i;
  double d;

  for (i = 0; i <= 2; i++) {
    t->U1.a[0][i] = sx * t->U1.a[0][i];
    t->U1.a[1][i] = sy * t->U1.a[1][i];
  }
  d = sx * sy;
  if (d < 0.0)
    t->U1.detsgn = (short)(-t->U1.detsgn);
  else if ( !d )
    t->U1.detsgn = 0;
} /*ScaleTrans2d*/

void Trans2Shiftd ( trans2d *tr, double tx, double ty )
{
  tr->U0.a13 += tr->U0.a11*tx + tr->U0.a12*ty;
  tr->U0.a23 += tr->U0.a21*tx + tr->U0.a22*ty;
} /*Trans2Shiftd*/

void Trans2Rotd ( trans2d *tr, double angle )
{
  double c, s, m, n;
  byte i;

  if ( !angle )
    return;
  c = cos ( angle );
  s = sin ( angle );
  for ( i = 0; i < 2; i++ ) {
    m = tr->U1.a[i][0];
    n = tr->U1.a[i][1];
    tr->U1.a[i][0] = (c * m + s * n);
    tr->U1.a[i][1] = (s * m - c * n);
  }
} /*Trans2Rotd*/

void Trans2Scaled ( trans2d *tr, double sx, double sy )
{
  double d;

  tr->U0.a11 *= sx;  tr->U0.a12 *= sy;
  tr->U0.a21 *= sx;  tr->U0.a22 *= sy;
  d = sx * sy;
  if ( d < 0.0 )
    tr->U1.detsgn = (short)(-tr->U1.detsgn);
  else if ( !d )
    tr->U1.detsgn = 0;
} /*Trans2Scaled*/

boolean InvertTrans2d ( trans2d *tr )
{
  double d;
  trans2d b;

  if (tr->U1.detsgn == 0)
    return false;
  d = 1.0 / (tr->U0.a11 * tr->U0.a22 - tr->U0.a12 * tr->U0.a21);
  b.U0.a11 = tr->U0.a22 * d;
  b.U0.a12 = -tr->U0.a12 * d;
  b.U0.a13 = (tr->U0.a12 * tr->U0.a23 - tr->U0.a22 * tr->U0.a13) * d;
  b.U0.a21 = -tr->U0.a21 * d;
  b.U0.a22 = tr->U0.a11 * d;
  b.U0.a23 = (tr->U0.a13 * tr->U0.a21 - tr->U0.a11 * tr->U0.a23) * d;
  b.U1.detsgn = tr->U1.detsgn;
  *tr = b;
  return true;
} /*InvertTrans2d*/

void MultVector2d ( double a, const vector2d *v, vector2d *w )
{
  w->x = a * v->x;
  w->y = a * v->y;
} /*MultVector2d*/

void SubtractPoints2d ( const point2d *p1, const point2d * p2, vector2d *v )
{
  v->x = p1->x - p2->x;
  v->y = p1->y - p2->y;
} /*SubtractPoints2d*/

void AddVector2d ( const point2d *p, const vector2d *v, point2d *q )
{
  q->x = p->x + v->x;
  q->y = p->y + v->y;
} /*AddVector2d*/

void AddVector2Md ( const point2d *p, const vector2d *v, double t, point2d *q )
{
  q->x = p->x + t * v->x;
  q->y = p->y + t * v->y;
} /*AddVector2Md*/

void InterPoint2d ( const point2d *p1, const point2d *p2, double t, point2d *q )
{
  double s;

  s = 1.0 - t;
  q->x = s * p1->x + t * p2->x;
  q->y = s * p1->y + t * p2->y;
} /*InterPointd*/

void MidPoint2d ( const point2d *p1, const point2d *p2, point2d *q )
{
  q->x = 0.5 * (p1->x + p2->x);
  q->y = 0.5 * (p1->y + p2->y);
} /*MidPoint2d*/

void Interp3Vectors2d ( const vector2d *p0, const vector2d *p1, const vector2d *p2,
                        const double *coeff, vector2d *p )
{
  p->x = coeff[0]*p0->x + coeff[1]*p1->x + coeff[2]*p2->x;
  p->y = coeff[0]*p0->y + coeff[1]*p1->y + coeff[2]*p2->y;
} /*Interp3Vectors2d*/

void NormalizeVector2d ( vector2d *v )
{
  double l;

  if (v->x == 0.0 && v->y == 0.0)
    return;
  l = sqrt ( v->x * v->x + v->y * v->y );
  v->x /= l;
  v->y /= l;
} /*NormalizeVector2d*/

double DotProduct2d ( const vector2d *v1, const vector2d *v2 )
{
  return (v1->x * v2->x + v1->y * v2->y);
} /*DotProduct2d*/

double det2d ( const vector2d *v1, const vector2d *v2 )
{
  return (v1->x * v2->y - v1->y * v2->x);
} /*det2d*/

void OrtVector2d ( const vector2d *v1, const vector2d *v2, vector2d *v )
{
#define tol 1.0e-15
  double nv, t;

  nv = DotProduct2d ( v1, v1 );
  if (nv < tol)
    *v = *v2;   
  else {
    t = DotProduct2d ( v1, v2 ) / nv;
    SetVector2d ( v, v2->x - t * v1->x, v2->y - t * v1->y );
  }
#undef tol
} /*OrtVector2d*/

void ProjectPointOnLine2d ( const point2d *p0, const point2d *p1,
                            point2d *q )
{
  vector2d v, w;
  double   d, e;

  SubtractPoints2d ( p1, p0, &v );
  SubtractPoints2d ( q, p0, &w );
  d = DotProduct2d ( &v, &v );
  e = DotProduct2d ( &v, &w );
  AddVector2Md ( p0, &v, e/d, q );
} /*ProjectPointOnLine2d*/

double Point2Distanced ( point2d *p, point2d *q )
{
  vector2d v;

  SubtractPoints2d ( p, q, &v );
  return sqrt ( DotProduct2d ( &v, &v ) );
} /*Point2Distanced*/

