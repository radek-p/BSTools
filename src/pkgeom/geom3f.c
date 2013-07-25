
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <math.h>

#include "pknum.h"
#include "pkgeom.h"


void SetPoint3f ( point3f *p, float x, float y, float z )
{
  p->x = x;
  p->y = y;
  p->z = z;
} /*SetPoint3f*/

void Point3to2f ( const point3f *P, point2f *p )
{
  SetPoint2f ( p, P->x/P->z, P->y/P->z );
} /*Point3to2f*/

void Point2to3f ( const point2f *p, float w, point3f *P )
{
  SetPoint3f ( P, p->x*w, p->y*w, w );
} /*Point2to3f*/

void TransPoint3f ( const trans3f *tr, const point3f *p, point3f *q )
{
  /* q := tr * p */
  point3f r;

  r.x = tr->U0.a11 * p->x + tr->U0.a12 * p->y + tr->U0.a13 * p->z + tr->U0.a14;
  r.y = tr->U0.a21 * p->x + tr->U0.a22 * p->y + tr->U0.a23 * p->z + tr->U0.a24;
  r.z = tr->U0.a31 * p->x + tr->U0.a32 * p->y + tr->U0.a33 * p->z + tr->U0.a34;
  *q = r;
} /*TransPoint3f*/

void TransVector3f ( const trans3f *tr, const vector3f *v, vector3f *w )
{
  /* w := tr * v */
  vector3f r;

  r.x = tr->U0.a11 * v->x + tr->U0.a12 * v->y + tr->U0.a13 * v->z;
  r.y = tr->U0.a21 * v->x + tr->U0.a22 * v->y + tr->U0.a23 * v->z;
  r.z = tr->U0.a31 * v->x + tr->U0.a32 * v->y + tr->U0.a33 * v->z;
  *w = r;
} /*TransVector3f*/

void TransContra3f ( const trans3f *tri, const vector3f *v, vector3f *w )
{
  /* w := v * tri; v is a contravariant vector,                */
  /* tri is an inversion of the transformation being evaluated */
  vector3f r;

  r.x = tri->U0.a11 * v->x + tri->U0.a21 * v->y + tri->U0.a31 * v->z;
  r.y = tri->U0.a12 * v->x + tri->U0.a22 * v->y + tri->U0.a32 * v->z;
  r.z = tri->U0.a13 * v->x + tri->U0.a23 * v->y + tri->U0.a33 * v->z;
  *w = r;
} /*TransContra3f*/

void Trans3Point2f ( const trans3f *tr, const point2f *p, point2f *q )
{
  point2f r;

  r.x = tr->U0.a11*p->x + tr->U0.a12*p->y + tr->U0.a14;
  r.y = tr->U0.a21*p->x + tr->U0.a22*p->y + tr->U0.a24;
  *q = r;
} /*Trans3Point2f*/

/* ////////////////////////////////////////////////////////////////////////// */
void Trans3Shiftf ( trans3f *tr, float tx, float ty, float tz )
{
  tr->U0.a14 += tr->U0.a11*tx + tr->U0.a12*ty + tr->U0.a13*tz;
  tr->U0.a24 += tr->U0.a21*tx + tr->U0.a22*ty + tr->U0.a23*tz;
  tr->U0.a34 += tr->U0.a31*tx + tr->U0.a32*ty + tr->U0.a33*tz;
} /*Trans3Shiftf*/

void Trans3Mirrorf ( trans3f *tr, vector3f *n )
{
  /* tr := tr * M(n) */
  vector3f h; 
  float    s;
  byte     i;

  h = *n;  
  NormalizeVector3f ( &h );
  for ( i = 0; i < 3; i++ ) { 
    s = -2.0 * (h.x * tr->U1.a[i][0] + h.y * tr->U1.a[i][1] + h.z * tr->U1.a[i][2]);   
    tr->U1.a[i][0] += s * h.x;   
    tr->U1.a[i][1] += s * h.y;   
    tr->U1.a[i][2] += s * h.z;   
  }
  tr->U1.detsgn = (short)(-tr->U1.detsgn);
} /*Trans3Mirrorf*/

void Trans3Rotf ( trans3f *tr, byte j, byte k, float angle )
{
  float c, s;
  float m, n;
  byte i;

  if ( !angle )
    return;
  c = (float)cos ( angle );
  s = (float)sin ( angle );
  for ( i = 0; i < 3; i++ ) {
    m = tr->U1.a[i][j-1];
    n = tr->U1.a[i][k-1];
    tr->U1.a[i][j-1] = c * m + s * n;
    tr->U1.a[i][k-1] = c * n - s * m;
  }
} /*Trans3Rotf*/

void Trans3RotVf ( trans3f *tr, vector3f *v, float angle )
{
  /* rotate by the angle, around the axis parallel to the v vector */
  /* v must be nonzero.                                            */
  vector3f vv, e, f;
  float    c, s, t;
  byte j;  
  trans3f tt;

  vv = *v;
  NormalizeVector3f ( &vv );
  s = (float)sin ( angle );
  c = (float)cos ( angle );
  for ( j = 0; j < 3; j++ ) {
    e.x = tr->U1.a[j][0];   
    e.y = tr->U1.a[j][1];
    e.z = tr->U1.a[j][2];
    t = vv.x * e.x + vv.y * e.y + vv.z * e.z;
    CrossProduct3f ( &e, &vv, &f );
    t = (1.0 - c) * t;
    tt.U1.a[j][0] = t * vv.x + c * tr->U1.a[j][0] + s * f.x;   
    tt.U1.a[j][1] = t * vv.y + c * tr->U1.a[j][1] + s * f.y;   
    tt.U1.a[j][2] = t * vv.z + c * tr->U1.a[j][2] + s * f.z;   
  }
  memcpy ( tr->U1.a, tt.U1.a, 12*sizeof(float) );
} /*Trans3RotVf*/

void Trans3Scalef ( trans3f *tr, float sx, float sy, float sz )
{
  double d;

  tr->U0.a11 *= sx;  tr->U0.a12 *= sy;  tr->U0.a13 *= sz;
  tr->U0.a21 *= sx;  tr->U0.a22 *= sy;  tr->U0.a23 *= sz;
  tr->U0.a31 *= sx;  tr->U0.a32 *= sy;  tr->U0.a33 *= sz;
  d = sx*sy*sz;
  if ( d < 0.0 )
    tr->U1.detsgn = (short)(-tr->U1.detsgn);
  else if ( !d )
    tr->U1.detsgn = 0;
} /*Trans3Scalef*/

/* ////////////////////////////////////////////////////////////////////////// */
void IdentTrans3f ( trans3f *tr )
{
  /* tr := I */
  int i, j;

  for ( i = 0; i <= 2; i++ ) {
    for ( j = 0; j <= 3; j++ )
      tr->U1.a[i][j] = 0.0;
    tr->U1.a[i][i] = 1.0;
  }
  tr->U1.detsgn = 1;
} /*IdentTrans3f*/

void CompTrans3f ( trans3f *s, trans3f *t, trans3f *u )
{
  /* s := t * u */
  byte i, j, k;
  double d;

  for ( i = 0; i <= 2; i++ ) {
    for ( j = 0; j <= 3; j++ ) {
      d = 0.0;
      for ( k = 0; k <= 2; k++ )
	d += t->U1.a[i][k] * u->U1.a[k][j];
      s->U1.a[i][j] = (float)d;
    }
    s->U1.a[i][3] += t->U1.a[i][3];
  }
  s->U1.detsgn = (short)(t->U1.detsgn * u->U1.detsgn);
} /*CompTrans3f*/

void GeneralAffineTrans3f ( trans3f *tr, vector3f *v1, vector3f *v2, vector3f *v3 )
{
  trans3f s, t;
  double  d;

  t.U0.a11 = v1->x;  t.U0.a12 = v2->x;  t.U0.a13 = v3->x;
  t.U0.a21 = v1->y;  t.U0.a22 = v2->y;  t.U0.a23 = v3->y;
  t.U0.a31 = v1->z;  t.U0.a32 = v2->z;  t.U0.a33 = v3->z;

  t.U0.a14 = 0.0;  t.U0.a24 = 0.0;  t.U0.a34 = 0.0;

  d = det3f ( v1, v2, v3 );
  if ( d > 0.0 )
    t.U1.detsgn = 1;
  else if ( d < 0.0 )
    t.U1.detsgn = -1;
  else
    t.U1.detsgn = 0;
  CompTrans3f ( &s, &t, tr );
  *tr = s;
} /*GeneralAffineTrans3f*/

void ShiftTrans3f ( trans3f *tr, float tx, float ty, float tz )
{
  /* tr := T(tx,ty,tz) * tr */
  tr->U0.a14 += tx;
  tr->U0.a24 += ty;
  tr->U0.a34 += tz;
} /*ShiftTrans3f*/

void MirrorTrans3f ( trans3f *tr, vector3f *n )
{
  /* tr := M(n) * tr */
  vector3f h;
  double   s;
  byte     i;

  h = *n;
  NormalizeVector3f ( &h );
  for (i = 0; i <= 3; i++) {
    s = -2.0 * (h.x * tr->U1.a[0][i] + h.y * tr->U1.a[1][i] + h.z * tr->U1.a[2][i]);
    tr->U1.a[0][i] += (float)(s * h.x);
    tr->U1.a[1][i] += (float)(s * h.y);
    tr->U1.a[2][i] += (float)(s * h.z);
  }
  tr->U1.detsgn = (short)(-tr->U1.detsgn);
} /*MirrorTrans3f*/

void RotTrans3f ( trans3f *tr, byte j, byte k, float angle )
{
  float c, s;
  float m, n;
  byte i;

  if ( !angle )
    return;
  c = (float)cos ( angle );
  s = (float)sin ( angle );
  for (i = 0; i <= 3; i++) {
    m = tr->U1.a[j-1][i];
    n = tr->U1.a[k-1][i];
    tr->U1.a[j-1][i] = (float)(c * m - s * n);
    tr->U1.a[k-1][i] = (float)(s * m + c * n);
  }
} /*RotTrans3f*/

void EulerRotTrans3f ( trans3f *tr, float psi, float theta, float phi )
{
  /* tr := E(psi,theta,phi) * tr */
  RotTrans3f ( tr, 1, 2, phi );
  RotTrans3f ( tr, 2, 3, theta );
  RotTrans3f ( tr, 1, 2, psi );
} /*EulerRotTrans3f*/

void FindRotVEulerf ( const vector3f *v, float angle,
                      float *psi, float *theta, float *phi )
{
  /* find Euler angles of rotation by angle around the axis parallel */
  /* to the v vector. v must be nonzero.                             */
  float a, b;

  if ( fabs(v->x) + fabs(v->y) > 0.0 ) {
    a = (float)atan2(v->x, v->y);
    b = (float)atan2(sqrt(v->x * v->x + v->y * v->y), v->z);
    CompEulerRotf ( 0.0, -b, -a, a, b, angle, psi, theta, phi );
  }
  else {
    *psi = 0.0;
    *theta = 0.0;
    if (v->z > 0.0)      *phi =  angle;
    else if (v->z < 0.0) *phi = -angle;
    else                 *phi = 0.0;
  }
} /*FindRotVEulerf*/

void RotVTrans3f ( trans3f *tr, vector3f *v, float angle )
{
  /* rotate by the angle, around the axis parallel to the v vector */
  /* v must be nonzero.                                            */
  vector3f vv, e, f;
  float    c, s, t;
  byte j;
  trans3f tt;

  vv = *v;
  NormalizeVector3f ( &vv );
  s = (float)sin ( angle );
  c = (float)cos ( angle );
  for (j = 0; j <= 3; j++) {
    e.x = tr->U1.a[0][j];
    e.y = tr->U1.a[1][j];
    e.z = tr->U1.a[2][j];
    t = vv.x * e.x + vv.y * e.y + vv.z * e.z;
    CrossProduct3f ( &vv, &e, &f );
    t = (1.0 - c) * t;
    tt.U1.a[0][j] = (float)(t * vv.x + c * tr->U1.a[0][j] + s * f.x);
    tt.U1.a[1][j] = (float)(t * vv.y + c * tr->U1.a[1][j] + s * f.y);
    tt.U1.a[2][j] = (float)(t * vv.z + c * tr->U1.a[2][j] + s * f.z);
  }
  memcpy ( tr->U1.a, tt.U1.a, 12*sizeof(float) );
} /*RotVTrans3f*/

void ScaleTrans3f ( trans3f *tr, float sx, float sy, float sz )
{
  /* tr := S(sx,sy,sz) * tr */
  byte   i;
  double d;

  for (i = 0; i <= 3; i++) {
    tr->U1.a[0][i] = sx * tr->U1.a[0][i];
    tr->U1.a[1][i] = sy * tr->U1.a[1][i];
    tr->U1.a[2][i] = sz * tr->U1.a[2][i];
  }
  d = sx * sy * sz;
  if ( d < 0.0 )
    tr->U1.detsgn = (short)(-tr->U1.detsgn);
  else if ( !d )
    tr->U1.detsgn = 0;
} /*ScaleTrans3f*/

boolean InvertTrans3f ( trans3f *tr )
{
  /* if tr not singular then tr := inversion(tr) */
  float b[16];

  memcpy ( b, &tr->U0.a11, 12*sizeof(float) );
  b[12] = b[13] = b[14] = 0.0;  b[15] = 1.0;
  if ( !pkn_GaussInvertMatrixf ( 4, b ) )
    return false;
  memcpy ( &tr->U0.a11, b, 12*sizeof(float) );
  return true;
} /*InvertTrans3d*/

float TrimAnglef ( float angle )
{
  /* result := angle between -pi and pi */
  long a;

  if (angle > PI) {
    a = ((long)(angle / PI) + 1) / 2;
    angle -= (float)(a * 2.0 * PI);
    return angle;
  }
  if (angle < -PI) {
    a = ((long)(-(angle / PI)) + 1) / 2;
    angle += (float)(a * 2.0 * PI);
  }
  return angle;
} /*TrimAnglef*/

void CompEulerRotf ( float psi1, float theta1, float phi1,
                     float psi2, float theta2, float phi2,
                     float *psi, float *theta, float *phi )
{
  /* E(psi, theta, phi) = E(psi2, theta2, phi2) * E(psi1, theta1, phi1)    */
  /* is a transformation - rotation by Euler angles psi, theta, phi, which */
  /* is a superposition of two rotations by Euler angles. This procedure   */
  /* finds psi, theta, phi given psi1 .. phi2 - values of Euler angles of  */
  /* the two rotations being superposed.                                   */
#define eps 1.0e-6
  double s_alpha, c_alpha, s_beta, c_beta, s_gamma, c_gamma, a31, a32,
	 s_theta, c_theta, _psi, _theta, _phi;

  /* 1st check for particular cases */
  if (fabs(theta2) <= eps)
  {
    _psi = psi2 + phi2 + psi1;
    _theta = theta1;
    _phi = phi1;
  }
  else if (fabs(theta1) <= eps)
  {
    _psi = psi2;
    _theta = theta2;
    _phi = phi2 + psi1 + phi1;
  }
  else if (fabs(phi2 + psi1) <= eps)
  {
    _psi = psi2;
    _theta = theta2 + theta1;
    _phi = phi1;
  }
  else
  {
    s_alpha = sin ( theta2 );
    c_alpha = cos ( theta2 );
    s_beta = sin ( phi2 + psi1 );
    c_beta = cos ( phi2 + psi1 );
    s_gamma = sin ( theta1 );
    c_gamma = cos ( theta1 );

    a31 = s_alpha * s_beta;
    a32 = s_alpha * c_beta * c_gamma + c_alpha * s_gamma;

    s_theta = sqrt ( a31 * a31 + a32 * a32 );
    c_theta = c_alpha * c_gamma - s_alpha * c_beta * s_gamma;

    if (s_theta > eps)
    {
      _phi = atan2 ( a31, a32 ) + phi1;
      _psi = atan2 ( s_gamma * s_beta,
		     c_alpha * c_beta * s_gamma + s_alpha * c_gamma ) + psi2;
    }
    else
    {   /* only the sum or difference of psi and phi is */
      s_theta = 0.0;
      _psi = psi2;
      /* determined. One of them may be arbitrary     */
      _phi = atan2 ( c_alpha * s_beta / c_theta, c_beta ) + phi1;
    }

    _theta = atan2 ( s_theta, c_theta );
    /* the general case */
  }
  *psi = TrimAnglef ( (float)_psi );
  *theta = TrimAnglef ( (float)_theta );
  *phi = TrimAnglef ( (float)_phi );
#undef eps
} /*CompEulerRotf*/

void CompRotV3f ( const vector3f *v1, float a1, const vector3f *v2, float a2,
                  vector3f *v, float *a )
{
  double c1, s1, c2, s2, c, s;
  vector3f w1, w2, w;

  w1 = *v1;  NormalizeVector3f ( &w1 );
  c1 = cos ( 0.5*a1 );  s1 = sin ( 0.5*a1 );
  w2 = *v2;  NormalizeVector3f ( &w2 );
  c2 = cos ( 0.5*a2 );  s2 = sin ( 0.5*a2 );
  c = c1*c2 - s1*s2*DotProduct3f ( &w1, &w2 );
  CrossProduct3f ( &w2, &w1, &w );
  MultVector3f ( s1*s2, &w, v );
  AddVector3Mf ( v, &w2, s2*c1, v );
  AddVector3Mf ( v, &w1, s1*c2, v );
  s = sqrt ( DotProduct3f ( v, v ) );
  *a = (float)(2.0*atan2 ( s, c ));
  if ( s )
    MultVector3f ( 1.0/s, v, v );
} /*CompRotV3f*/

void MultVector3f ( double a, const vector3f *v, vector3f *w )
{
  w->x = (float)(a * v->x);
  w->y = (float)(a * v->y);
  w->z = (float)(a * v->z);
} /*MultVector3f*/

void SubtractPoints3f ( const point3f *p1, const point3f *p2, vector3f *v )
{
  /* v := p1 - p2 */
  v->x = p1->x - p2->x;
  v->y = p1->y - p2->y;
  v->z = p1->z - p2->z;
} /*SubtractPoints3f*/

void AddVector3f ( const point3f *p, const vector3f *v, point3f *q )
{
  /* q := p + v */
  q->x = p->x + v->x;
  q->y = p->y + v->y;
  q->z = p->z + v->z;
} /*AddVector3f*/

void AddVector3Mf ( const point3f *p, const vector3f *v, double t, point3f *q )
{
  /* q := p + t * v */
  q->x = (float)(p->x + t * v->x);
  q->y = (float)(p->y + t * v->y);
  q->z = (float)(p->z + t * v->z);
} /*AddVector3Mf*/

void InterPoint3f ( const point3f *p1, const point3f *p2, double t, point3f *q )
{
  double s;

  s = 1.0 - t;
  q->x = (float)(s * p1->x + t * p2->x);
  q->y = (float)(s * p1->y + t * p2->y);
  q->z = (float)(s * p1->z + t * p2->z);
} /*InterPoint3f*/

void MidPoint3f ( const point3f *p1, const point3f *p2, point3f *q )
{
  q->x = (float)0.5 * (p1->x + p2->x);
  q->y = (float)0.5 * (p1->y + p2->y);
  q->z = (float)0.5 * (p1->z + p2->z);
} /*MidPoint3f*/

void Interp3Vectors3f ( const vector3f *p0, const vector3f *p1, const vector3f *p2,
                        const float *coeff, vector3f *p )
{
  p->x = coeff[0]*p0->x + coeff[1]*p1->x + coeff[2]*p2->x;
  p->y = coeff[0]*p0->y + coeff[1]*p1->y + coeff[2]*p2->y;
  p->z = coeff[0]*p0->z + coeff[1]*p1->z + coeff[2]*p2->z;
} /*Interp3Vectors3f*/

void NormalizeVector3f ( vector3f *v )
{
  /* if v <> 0 then v := v/length(v) */
#define tol 1.0e-15
  double r;

  r = sqrt ( v->x * v->x + v->y * v->y + v->z * v->z );
  if (r <= tol)
    return;
  r = 1.0 / r;
  v->x *= (float)r;
  v->y *= (float)r;
  v->z *= (float)r;
#undef tol
} /*NormalizeVector3f*/

double DotProduct3f ( const vector3f *v1, const vector3f *v2 )
{
  return (v1->x * v2->x + v1->y * v2->y + v1->z * v2->z);
} /*DotProduct3f*/

void CrossProduct3f ( const vector3f *v1, const vector3f *v2, vector3f *v )
{
  vector3f vv;

  vv.x = v1->y * v2->z - v1->z * v2->y;
  vv.y = v1->z * v2->x - v1->x * v2->z;
  vv.z = v1->x * v2->y - v1->y * v2->x;
  *v = vv;
} /*CrossProduct3f*/

void AddCrossProduct3f ( const vector3f *w, const vector3f *v1, const vector3f *v2,
                         vector3f *v )
{
  vector3f vv;

  vv.x = v1->y * v2->z - v1->z * v2->y;
  vv.y = v1->z * v2->x - v1->x * v2->z;
  vv.z = v1->x * v2->y - v1->y * v2->x;
  AddVector3f ( w, &vv, v );
} /*AddCrossProduct3f*/

double det3f ( const vector3f *v1, const vector3f *v2, const vector3f *v3 )
{
  float a[9];

  memcpy ( &a[0], v1, sizeof(vector3f) );
  memcpy ( &a[3], v2, sizeof(vector3f) );
  memcpy ( &a[6], v3, sizeof(vector3f) );
  return pkn_detf ( 3, a );
}  /*det3f*/

void OrtVector3f ( const vector3f *v1, const vector3f *v2, vector3f *v )
{
#define tol 1.0e-15
  double nv, t;

  nv = DotProduct3f ( v1, v1 );
  if (nv < tol)
    *v = *v2;
  else {
    t = DotProduct3f ( v1, v2 ) / nv;
    SetVector3f ( v, (float)(v2->x - t * v1->x), (float)(v2->y - t * v1->y),
                  (float)(v2->z - t * v1->z) );
  }
#undef tol
} /*OrtVector3f*/

void ProjectPointOnLine3f ( const point3f *p0, const point3f *p1,
                            point3f *q )
{
  vector3f v, w;
  float    d, e;

  SubtractPoints3f ( p1, p0, &v );
  SubtractPoints3f ( q, p0, &w );
  d = DotProduct3f ( &v, &v );
  e = DotProduct3f ( &v, &w );
  AddVector3Mf ( p0, &v, e/d, q );
} /*ProjectPointOnLine3f*/

void ProjectPointOnPlane3f ( const point3f *p0, const point3f *p1, const point3f *p2,
                             point3f *q )
{
  vector3f v1, v2, w;
  float    d, e;

  SubtractPoints3f ( p1, p0, &v1 );
  SubtractPoints3f ( p2, p0, &v2 );
  OrtVector3f ( &v1, &v2, &v2 );
  SubtractPoints3f ( q, p0, &w );
  d = DotProduct3f ( &v1, &v1 );
  e = DotProduct3f ( &v1, &w );
  AddVector3Mf ( p0, &v1, e/d, q );
  d = DotProduct3f ( &v2, &v2 );
  e = DotProduct3f ( &v2, &w );
  AddVector3Mf ( q, &v2, e/d, q );
} /*ProjectPointOnPlane3f*/

double Point3Distancef ( point3f *p, point3f *q )
{
  vector3f v;

  SubtractPoints3f ( p, q, &v );
  return sqrt ( DotProduct3f ( &v, &v ) );
} /*Point3Distancef*/

