
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2011                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>

#include "pkgeom.h"

void SetPoint2f ( point2f *p, float x, float y )
{
  p->x = x;
  p->y = y;
} /*SetPoint2f*/

void TransPoint2f ( const trans2f *tr, const point2f *p, point2f *q )
{
  point2f r;

  r.x = tr->U0.a11 * p->x + tr->U0.a12 * p->y + tr->U0.a13;
  r.y = tr->U0.a21 * p->x + tr->U0.a22 * p->y + tr->U0.a23;
  *q = r;
} /*TransPoint2f*/

void TransVector2f ( const trans2f *tr, const vector2f *v, vector2f *w )
{
  vector2f r;

  r.x = tr->U0.a11 * v->x + tr->U0.a12 * v->y;
  r.y = tr->U0.a21 * v->x + tr->U0.a22 * v->y;
  *w = r;
} /*TransVector2f*/

void Trans2Point3f ( const trans2f *tr, const point3f *p, point3f *q )
{
  point3f r;

  r.x = tr->U0.a11*p->x + tr->U0.a12*p->y + tr->U0.a13*p->z;
  r.y = tr->U0.a21*p->x + tr->U0.a22*p->y + tr->U0.a23*p->z;
  r.z = p->z;
  *q = r;
} /*Trans2Point3f*/

void IdentTrans2f ( trans2f *tr )
{
  short i, j;

  for (i = 0; i <= 1; i++) {
    for (j = 0; j <= 2; j++)
      tr->U1.a[i][j] = 0.0;
    tr->U1.a[i][i] = 1.0;
  }
  tr->U1.detsgn = 1;
} /*IdentTrans2f*/

void CompTrans2f ( trans2f *s, trans2f *t, trans2f *u )
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
} /*CompTrans2f*/

void ShiftTrans2f ( trans2f *tr, float tx, float ty )
{
  tr->U0.a13 += tx;
  tr->U0.a23 += ty;
} /*ShiftTrans2f*/

void RotTrans2f ( trans2f *tr, float angle )
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
    tr->U1.a[0][i] = (float)(c * m - s * n);
    tr->U1.a[1][i] = (float)(s * m + c * n);
  }
} /*RotTrans2f*/

void ScaleTrans2f ( trans2f *t, float sx, float sy )
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
} /*ScaleTrans2f*/

void Trans2Shiftf ( trans2f *tr, float tx, float ty )
{
  tr->U0.a13 += tr->U0.a11*tx + tr->U0.a12*ty;
  tr->U0.a23 += tr->U0.a21*tx + tr->U0.a22*ty;
} /*Trans2Shiftf*/

void Trans2Rotf ( trans2f *tr, float angle )
{
  double c, s, m, n;
  byte i;

  if ( !angle )
    return;
  c = (float)cos ( angle );
  s = (float)sin ( angle );
  for ( i = 0; i < 2; i++ ) {
    m = tr->U1.a[i][0];
    n = tr->U1.a[i][1];
    tr->U1.a[i][0] = (float)(c * m + s * n);
    tr->U1.a[i][1] = (float)(s * m - c * n);
  }
} /*Trans2Rotf*/

void Trans2Scalef ( trans2f *tr, float sx, float sy )
{
  double d;

  tr->U0.a11 *= sx;  tr->U0.a12 *= sy;
  tr->U0.a21 *= sx;  tr->U0.a22 *= sy;
  d = sx * sy;
  if ( d < 0.0 )
    tr->U1.detsgn = (short)(-tr->U1.detsgn);
  else if ( !d )
    tr->U1.detsgn = 0;
} /*Trans2Scalef*/

boolean InvertTrans2f ( trans2f *tr )
{
  double  d;
  trans2f b;

  if (tr->U1.detsgn == 0)
    return false;
  d = 1.0 / (tr->U0.a11 * tr->U0.a22 - tr->U0.a12 * tr->U0.a21);
  b.U0.a11 = (float)(tr->U0.a22 * d);
  b.U0.a12 = (float)(-tr->U0.a12 * d);
  b.U0.a13 = (float)((tr->U0.a12 * tr->U0.a23 - tr->U0.a22 * tr->U0.a13) * d);
  b.U0.a21 = (float)(-tr->U0.a21 * d);
  b.U0.a22 = (float)(tr->U0.a11 * d);
  b.U0.a23 = (float)((tr->U0.a13 * tr->U0.a21 - tr->U0.a11 * tr->U0.a23) * d);
  b.U1.detsgn = tr->U1.detsgn;
  *tr = b;
  return true;
} /*InvertTrans2f*/

void MultVector2f ( double a, const vector2f *v, vector2f *w )
{
  w->x = (float)(a * v->x);
  w->y = (float)(a * v->y);
} /*MultVector2f*/

void SubtractPoints2f ( const point2f *p1, const point2f * p2, vector2f *v )
{
  v->x = p1->x - p2->x;
  v->y = p1->y - p2->y;
} /*SubtractPoints2f*/

void AddVector2f ( const point2f *p, const vector2f *v, point2f *q )
{
  q->x = p->x + v->x;
  q->y = p->y + v->y;
} /*AddVector2f*/

void AddVector2Mf ( const point2f *p, const vector2f *v, double t, point2f *q )
{
  q->x = (float)(p->x + t * v->x);
  q->y = (float)(p->y + t * v->y);
} /*AddVector2Mf*/

void InterPoint2f ( const point2f *p1, const point2f *p2, double t, point2f *q )
{
  double s;

  s = 1.0 - t;
  q->x = (float)(s * p1->x + t * p2->x);
  q->y = (float)(s * p1->y + t * p2->y);
} /*InterPointf*/

void MidPoint2f ( const point2f *p1, const point2f *p2, point2f *q )
{
  q->x = (float)0.5 * (p1->x + p2->x);
  q->y = (float)0.5 * (p1->y + p2->y);
} /*MidPoint2f*/

void Interp3Vectors2f ( const vector2f *p0, const vector2f *p1, const vector2f *p2,
                        const float *coeff, vector2f *p )
{
  p->x = coeff[0]*p0->x + coeff[1]*p1->x + coeff[2]*p2->x;
  p->y = coeff[0]*p0->y + coeff[1]*p1->y + coeff[2]*p2->y;
} /*Interp3Vectors2f*/

void NormalizeVector2f ( vector2f *v )
{
  double l;

  if (v->x == 0.0 && v->y == 0.0)
    return;
  l = sqrt ( v->x * v->x + v->y * v->y );
  v->x /= (float)l;
  v->y /= (float)l;
} /*NormalizeVector2f*/

double DotProduct2f ( const vector2f *v1, const vector2f *v2 )
{
  return (v1->x * v2->x + v1->y * v2->y);
} /*DotProduct2f*/

double det2f ( const vector2f *v1, const vector2f *v2 )
{
  return (v1->x * v2->y - v1->y * v2->x);
} /*det2f*/

void OrtVector2f ( const vector2f *v1, const vector2f *v2, vector2f *v )
{
#define tol 1.0e-15
  double nv, t;

  nv = DotProduct2f ( v1, v1 );
  if (nv < tol)
    *v = *v2;   
  else {
    t = DotProduct2f ( v1, v2 ) / nv;
    SetVector2f ( v, (float)(v2->x - t * v1->x), (float)(v2->y - t * v1->y) );
  }
#undef tol
} /*OrtVector2f*/

void ProjectPointOnLine2f ( const point2f *p0, const point2f *p1,
                            point2f *q )
{
  vector2f v, w;
  float    d, e;

  SubtractPoints2f ( p1, p0, &v );
  SubtractPoints2f ( q, p0, &w );
  d = DotProduct2f ( &v, &v );
  e = DotProduct2f ( &v, &w );
  AddVector2Mf ( p0, &v, e/d, q );
} /*ProjectPointOnLine2f*/

double Point2Distancef ( point2f *p, point2f *q )
{
  vector2f v;

  SubtractPoints2f ( p, q, &v );
  return sqrt ( DotProduct2f ( &v, &v ) );
} /*Point2Distancef*/

