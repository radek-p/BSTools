
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <math.h>

#include "pknum.h"
#include "pkgeom.h"


void SetPoint4f ( point4f *p, float x, float y, float z, float w )
{
  p->x = x;
  p->y = y;
  p->z = z;
  p->w = w;
} /*SetPoint4f*/

void Point4to3f ( const point4f *P, point3f *p )
{
  SetPoint3f ( p, P->x/P->w, P->y/P->w, P->z/P->w );
} /*Point4to3f*/

void Point3to4f ( const point3f *p, float w, point4f *P )
{
  SetPoint4f ( P, p->x*w, p->y*w, p->z*w, w );
} /*Point3to4f*/

void Point4to2f ( const point4f *P, point2f *p )
{
  SetPoint2f ( p, P->x/P->w, P->y/P->w );
} /*Point4to2f*/

void Point2to4f ( const point2f *p, float w, point4f *P )
{
  SetPoint4f ( P, p->x*w, p->y*w, 0.0, w );
} /*Point2to4f*/

void Trans3Point4f ( const trans3f *tr, const point4f *p, point4f *q )
{
  /* q := tr * p */
  point4f r;

  r.x = tr->U0.a11 * p->x + tr->U0.a12 * p->y + tr->U0.a13 * p->z + tr->U0.a14 * p->w;
  r.y = tr->U0.a21 * p->x + tr->U0.a22 * p->y + tr->U0.a23 * p->z + tr->U0.a24 * p->w;
  r.z = tr->U0.a31 * p->x + tr->U0.a32 * p->y + tr->U0.a33 * p->z + tr->U0.a34 * p->w;
  r.w = p->w;
  *q = r;
} /*Trans3Point4f*/

void MultVector4f ( double a, const vector4f *v, vector4f *w )
{
  w->x = (float)(a * v->x);
  w->y = (float)(a * v->y);
  w->z = (float)(a * v->z);
  w->w = (float)(a * v->w);
} /*MultVector4f*/

void SubtractPoints4f ( const point4f *p1, const point4f *p2, vector4f *v )
{
  v->x = p1->x - p2->x;
  v->y = p1->y - p2->y;
  v->z = p1->z - p2->z;
  v->w = p1->w - p2->w;
} /*SubtractPoints4f*/

void AddVector4f ( const point4f *p, const vector4f *v, point4f *q )
{
  /* q := p + v */
  q->x = p->x + v->x;
  q->y = p->y + v->y;
  q->z = p->z + v->z;
  q->w = p->w + v->w;
} /*AddVector4f*/

void AddVector4Mf ( const point4f *p, const vector4f *v, double t, point4f *q )
{
  /* q := p + t * v */
  q->x = (float)(p->x + t * v->x);
  q->y = (float)(p->y + t * v->y);
  q->z = (float)(p->z + t * v->z);
  q->w = (float)(p->w + t * v->w);
} /*AddVector4Mf*/

void InterPoint4f ( const point4f *p1, const point4f *p2, double t, point4f *q )
{
  double s;

  s = 1.0 - t;
  q->x = (float)(s * p1->x + t * p2->x);
  q->y = (float)(s * p1->y + t * p2->y);
  q->z = (float)(s * p1->z + t * p2->z);
  q->w = (float)(s * p1->w + t * p2->w);
} /*InterPoint4f*/

void MidPoint4f ( const point4f *p1, const point4f *p2, point4f *q )
{
  q->x = (float)0.5 * (p1->x + p2->x);
  q->y = (float)0.5 * (p1->y + p2->y);
  q->z = (float)0.5 * (p1->z + p2->z);
  q->w = (float)0.5 * (p1->w + p2->w);
} /*MidPoint4f*/

void Interp3Vectors4f ( const vector4f *p0, const vector4f *p1, const vector4f *p2,
                        const float *coeff, vector4f *p )
{
  p->x = coeff[0]*p0->x + coeff[1]*p1->x + coeff[2]*p2->x;
  p->y = coeff[0]*p0->y + coeff[1]*p1->y + coeff[2]*p2->y;
  p->z = coeff[0]*p0->z + coeff[1]*p1->z + coeff[2]*p2->z;
  p->w = coeff[0]*p0->w + coeff[1]*p1->w + coeff[2]*p2->w;
} /*Interp3Vectors4d*/

void NormalizeVector4f ( vector4f *v )
{
  /* if v <> 0 then v := v/length(v) */
#define tol 1.0e-15
  double r;

  r = sqrt ( v->x*v->x + v->y*v->y + v->z*v->z + v->w*v->w );
  if (r <= tol)
    return;
  r = 1.0 / r;
  v->x *= (float)r;
  v->y *= (float)r;
  v->z *= (float)r;
  v->w *= (float)r;
#undef tol
} /*NormalizeVector4f*/

double det4f ( const vector4f *v0, const vector4f *v1,
               const vector4f *v2, const vector4f *v3 )                          
{
  float a[16];

  memcpy ( &a[ 0], v0, sizeof(vector4f) );
  memcpy ( &a[ 4], v1, sizeof(vector4f) );
  memcpy ( &a[ 8], v2, sizeof(vector4f) );
  memcpy ( &a[12], v3, sizeof(vector4f) );
  return pkn_detf ( 4, a );
} /*det4f*/

void CrossProduct4f ( const vector4f *v0, const vector4f *v1,
                      const vector4f *v2, vector4f *v )
{
  double d12, d13, d14, d23, d24, d34;

  d12 = v0->x*v1->y - v0->y*v1->x;
  d13 = v0->x*v1->z - v0->z*v1->x;
  d14 = v0->x*v1->w - v0->w*v1->x;
  d23 = v0->y*v1->z - v0->z*v1->y;
  d24 = v0->y*v1->w - v0->w*v1->y;
  d34 = v0->z*v1->w - v0->w*v1->z;
  SetVector4f ( v,
               (float)(-d34*v2->y + d24*v2->z - d23*v2->w),
               (float)( d34*v2->x - d14*v2->z + d13*v2->w),
               (float)(-d24*v2->x + d14*v2->y - d12*v2->w),
               (float)( d23*v2->x - d13*v2->y + d12*v2->z) );
} /*CrossProduct4f*/

void CrossProduct4P3f ( const vector4f *v0, const vector4f *v1,
                        const vector4f *v2, vector3f *v )
{
  double d12, d13, d14, d23, d24, d34;

  d12 = v0->x*v1->y - v0->y*v1->x;
  d13 = v0->x*v1->z - v0->z*v1->x;
  d14 = v0->x*v1->w - v0->w*v1->x;
  d23 = v0->y*v1->z - v0->z*v1->y;
  d24 = v0->y*v1->w - v0->w*v1->y;
  d34 = v0->z*v1->w - v0->w*v1->z;
  SetVector3f ( v,
               (float)(-d34*v2->y + d24*v2->z - d23*v2->w),
               (float)( d34*v2->x - d14*v2->z + d13*v2->w),
               (float)(-d24*v2->x + d14*v2->y - d12*v2->w) );
} /*CrossProduct4P3f*/

double DotProduct4f ( const vector4f *v0, const vector4f *v1 )
{
  return v0->x*v1->x + v0->y*v1->y + v0->z*v1->z + v0->w*v1->w; 
} /*DotProduct4f*/

void OutProduct4P3f ( const vector4f *v0, const vector4f *v1, vector3f *v )
{
  SetVector3f ( v,
                v0->x*v1->w - v0->w*v1->x,
                v0->y*v1->w - v0->w*v1->y,
                v0->z*v1->w - v0->w*v1->z );
} /*OutProduct4P3f*/

void OrtVector4f ( const vector4f *v1, const vector4f *v2, vector4f *v )
{
#define tol 1.0e-15
  float nv, t;

  nv = DotProduct4f ( v1, v1 );
  if (nv < tol)
    *v = *v2;
  else {
    t = DotProduct4f ( v1, v2 ) / nv;
    SetVector4f ( v, v2->x - t * v1->x, v2->y - t * v1->y, v2->z - t * v1->z,
                  v2->w - t * v1->w );
  }
#undef tol
} /*OrtVector4f*/

void ProjectPointOnLine4f ( const point4f *p0, const point4f *p1,
                            point4f *q )
{
  vector4f v, w;
  float    d, e;

  SubtractPoints4f ( p1, p0, &v );
  SubtractPoints4f ( q, p0, &w );
  d = DotProduct4f ( &v, &v );
  e = DotProduct4f ( &v, &w );
  AddVector4Mf ( p0, &v, e/d, q );
} /*ProjectPointOnLine4f*/

void ProjectPointOnPlane4f ( const point4f *p0, const point4f *p1, const point4f *p2,
                             point4f *q )
{
  vector4f v1, v2, w;
  double   d, e;

  SubtractPoints4f ( p1, p0, &v1 );
  SubtractPoints4f ( p2, p0, &v2 );
  OrtVector4f ( &v1, &v2, &v2 );
  SubtractPoints4f ( q, p0, &w );
  d = DotProduct4f ( &v1, &v1 );
  e = DotProduct4f ( &v1, &w );
  AddVector4Mf ( p0, &v1, e/d, q );
  d = DotProduct4f ( &v2, &v2 );
  e = DotProduct4f ( &v2, &w );
  AddVector4Mf ( q, &v2, e/d, q );
} /*ProjectPointOnPlane4f*/

