
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


void SetPoint4d ( point4d *p, double x, double y, double z, double w )
{
  p->x = x;
  p->y = y;
  p->z = z;
  p->w = w;
} /*SetPoint4d*/

void Point4to3d ( const point4d *P, point3d *p )
{
  SetPoint3d ( p, P->x/P->w, P->y/P->w, P->z/P->w );
} /*Point4to3d*/

void Point3to4d ( const point3d *p, double w, point4d *P )
{
  SetPoint4d ( P, p->x*w, p->y*w, p->z*w, w );
} /*Point3to4d*/

void Point4to2d ( const point4d *P, point2d *p )
{
  SetPoint2d ( p, P->x/P->w, P->y/P->w );
} /*Point4to2d*/

void Point2to4d ( const point2d *p, float w, point4d *P )
{
  SetPoint4d ( P, p->x*w, p->y*w, 0.0, w );
} /*Point2to4d*/

void Trans3Point4d ( const trans3d *tr, const point4d *p, point4d *q )
{
  /* q := tr * p */
  point4d r;

  r.x = tr->U0.a11 * p->x + tr->U0.a12 * p->y + tr->U0.a13 * p->z + tr->U0.a14 * p->w;
  r.y = tr->U0.a21 * p->x + tr->U0.a22 * p->y + tr->U0.a23 * p->z + tr->U0.a24 * p->w;
  r.z = tr->U0.a31 * p->x + tr->U0.a32 * p->y + tr->U0.a33 * p->z + tr->U0.a34 * p->w;
  r.w = p->w;
  *q = r;
} /*Trans3Point4d*/

void MultVector4d ( double a, const vector4d *v, vector4d *w )
{
  w->x = a * v->x;
  w->y = a * v->y;
  w->z = a * v->z;
  w->w = a * v->w;
} /*MultVector4d*/

void SubtractPoints4d ( const point4d *p1, const point4d *p2, vector4d *v )
{
  v->x = p1->x - p2->x;
  v->y = p1->y - p2->y;
  v->z = p1->z - p2->z;
  v->w = p1->w - p2->w;
} /*SubtractPoints4d*/

void AddVector4d ( const point4d *p, const vector4d *v, point4d *q )
{
  /* q := p + v */
  q->x = p->x + v->x;
  q->y = p->y + v->y;
  q->z = p->z + v->z;
  q->w = p->w + v->w;
} /*AddVector4d*/

void AddVector4Md ( const point4d *p, const vector4d *v, double t, point4d *q )
{
  /* q := p + t * v */
  q->x = p->x + t * v->x;
  q->y = p->y + t * v->y;
  q->z = p->z + t * v->z;
  q->w = p->w + t * v->w;
} /*AddVector4Md*/

void InterPoint4d ( const point4d *p1, const point4d *p2, double t, point4d *q )
{
  double s;

  s = 1.0 - t;
  q->x = s * p1->x + t * p2->x;
  q->y = s * p1->y + t * p2->y;
  q->z = s * p1->z + t * p2->z;
  q->w = s * p1->w + t * p2->w;
} /*InterPoint4d*/

void MidPoint4d ( const point4d *p1, const point4d *p2, point4d *q )
{
  q->x = 0.5 * (p1->x + p2->x);
  q->y = 0.5 * (p1->y + p2->y);
  q->z = 0.5 * (p1->z + p2->z);
  q->w = 0.5 * (p1->w + p2->w);
} /*MidPoint4d*/

void Interp3Vectors4d ( const vector4d *p0, const vector4d *p1, const vector4d *p2,
                        const double *coeff, vector4d *p )
{
  p->x = coeff[0]*p0->x + coeff[1]*p1->x + coeff[2]*p2->x;
  p->y = coeff[0]*p0->y + coeff[1]*p1->y + coeff[2]*p2->y;
  p->z = coeff[0]*p0->z + coeff[1]*p1->z + coeff[2]*p2->z;
  p->w = coeff[0]*p0->w + coeff[1]*p1->w + coeff[2]*p2->w;
} /*Interp3Vectors4d*/

void NormalizeVector4d ( vector4d *v )
{
  /* if v <> 0 then v := v/length(v) */
#define tol 1.0e-15
  double r;

  r = sqrt ( v->x*v->x + v->y*v->y + v->z*v->z + v->w*v->w );
  if (r <= tol)
    return;
  r = 1.0 / r;
  v->x *= r;
  v->y *= r;
  v->z *= r;
  v->w *= r;
#undef tol
} /*NormalizeVector4d*/

double det4d ( const vector4d *v0, const vector4d *v1,
               const vector4d *v2, const vector4d *v3 )                          
{
  double a[16];

  memcpy ( &a[ 0], v0, sizeof(vector4d) );
  memcpy ( &a[ 4], v1, sizeof(vector4d) );
  memcpy ( &a[ 8], v2, sizeof(vector4d) );
  memcpy ( &a[12], v3, sizeof(vector4d) );
  return pkn_detd ( 4, a );
} /*det4d*/

void CrossProduct4d ( const vector4d *v0, const vector4d *v1,
                      const vector4d *v2, vector4d *v )
{
  double d12, d13, d14, d23, d24, d34;

  d12 = v0->x*v1->y - v0->y*v1->x;
  d13 = v0->x*v1->z - v0->z*v1->x;
  d14 = v0->x*v1->w - v0->w*v1->x;
  d23 = v0->y*v1->z - v0->z*v1->y;
  d24 = v0->y*v1->w - v0->w*v1->y;
  d34 = v0->z*v1->w - v0->w*v1->z;
  SetVector4d ( v,
               -d34*v2->y + d24*v2->z - d23*v2->w,
                d34*v2->x - d14*v2->z + d13*v2->w,
               -d24*v2->x + d14*v2->y - d12*v2->w,
                d23*v2->x - d13*v2->y + d12*v2->z );
} /*CrossProduct4d*/

void CrossProduct4P3d ( const vector4d *v0, const vector4d *v1,
                        const vector4d *v2, vector3d *v )
{
  double d12, d13, d14, d23, d24, d34;

  d12 = v0->x*v1->y - v0->y*v1->x;
  d13 = v0->x*v1->z - v0->z*v1->x;
  d14 = v0->x*v1->w - v0->w*v1->x;
  d23 = v0->y*v1->z - v0->z*v1->y;
  d24 = v0->y*v1->w - v0->w*v1->y;
  d34 = v0->z*v1->w - v0->w*v1->z;
  SetVector3d ( v,
               -d34*v2->y + d24*v2->z - d23*v2->w,
                d34*v2->x - d14*v2->z + d13*v2->w,
               -d24*v2->x + d14*v2->y - d12*v2->w );
} /*CrossProduct4P3d*/

double DotProduct4d ( const vector4d *v0, const vector4d *v1 )
{
  return v0->x*v1->x + v0->y*v1->y + v0->z*v1->z + v0->w*v1->w; 
} /*DotProduct4d*/

void OutProduct4P3d ( const vector4d *v0, const vector4d *v1, vector3d *v )
{
  SetVector3d ( v,
                v0->x*v1->w - v0->w*v1->x,
                v0->y*v1->w - v0->w*v1->y,
                v0->z*v1->w - v0->w*v1->z );
} /*OutProduct4P3d*/

void OrtVector4d ( const vector4d *v1, const vector4d *v2, vector4d *v )
{
#define tol 1.0e-15
  double nv, t;

  nv = DotProduct4d ( v1, v1 );
  if (nv < tol)
    *v = *v2;
  else {
    t = DotProduct4d ( v1, v2 ) / nv;
    SetVector4d ( v, v2->x - t * v1->x, v2->y - t * v1->y, v2->z - t * v1->z,
                  v2->w - t * v1->w );
  }
#undef tol
} /*OrtVector4d*/

void ProjectPointOnLine4d ( const point4d *p0, const point4d *p1,
                            point4d *q )
{
  vector4d v, w;
  double   d, e;

  SubtractPoints4d ( p1, p0, &v );
  SubtractPoints4d ( q, p0, &w );
  d = DotProduct4d ( &v, &v );
  e = DotProduct4d ( &v, &w );
  AddVector4Md ( p0, &v, e/d, q );
} /*ProjectPointOnLine4d*/

void ProjectPointOnPlane4d ( const point4d *p0, const point4d *p1, const point4d *p2,
                             point4d *q )
{
  vector4d v1, v2, w;
  double   d, e;

  SubtractPoints4d ( p1, p0, &v1 );
  SubtractPoints4d ( p2, p0, &v2 );
  OrtVector4d ( &v1, &v2, &v2 );
  SubtractPoints4d ( q, p0, &w );
  d = DotProduct4d ( &v1, &v1 );
  e = DotProduct4d ( &v1, &w );
  AddVector4Md ( p0, &v1, e/d, q );
  d = DotProduct4d ( &v2, &v2 );
  e = DotProduct4d ( &v2, &w );
  AddVector4Md ( q, &v2, e/d, q );
} /*ProjectPointOnPlane4d*/

