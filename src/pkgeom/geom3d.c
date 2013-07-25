
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


void SetPoint3d ( point3d *p, double x, double y, double z )
{
  p->x = x;
  p->y = y;
  p->z = z;
} /*SetPoint3d*/

void Point3to2d ( const point3d *P, point2d *p )
{
  SetPoint2d ( p, P->x/P->z, P->y/P->z );
} /*Point3to2d*/

void Point2to3d ( const point2d *p, double w, point3d *P )
{
  SetPoint3d ( P, p->x*w, p->y*w, w );
} /*Point2to3d*/

void TransPoint3d ( const trans3d *tr, const point3d *p, point3d *q )
{
  /* q := tr * p */
  point3d r;

  r.x = tr->U0.a11 * p->x + tr->U0.a12 * p->y + tr->U0.a13 * p->z + tr->U0.a14;
  r.y = tr->U0.a21 * p->x + tr->U0.a22 * p->y + tr->U0.a23 * p->z + tr->U0.a24;
  r.z = tr->U0.a31 * p->x + tr->U0.a32 * p->y + tr->U0.a33 * p->z + tr->U0.a34;
  *q = r;
} /*TransPoint3d*/

void TransVector3d ( const trans3d *tr, const vector3d *v, vector3d *w )
{
  /* w := tr * v */
  vector3d r;

  r.x = tr->U0.a11 * v->x + tr->U0.a12 * v->y + tr->U0.a13 * v->z;
  r.y = tr->U0.a21 * v->x + tr->U0.a22 * v->y + tr->U0.a23 * v->z;
  r.z = tr->U0.a31 * v->x + tr->U0.a32 * v->y + tr->U0.a33 * v->z;
  *w = r;
} /*TransVector3d*/

void TransContra3d ( const trans3d *tri, const vector3d *v, vector3d *w )
{
  /* w := v * tri; v is a contravariant vector,                */
  /* tri is an inversion of the transformation being evaluated */
  vector3d r;

  r.x = tri->U0.a11 * v->x + tri->U0.a21 * v->y + tri->U0.a31 * v->z;
  r.y = tri->U0.a12 * v->x + tri->U0.a22 * v->y + tri->U0.a32 * v->z;
  r.z = tri->U0.a13 * v->x + tri->U0.a23 * v->y + tri->U0.a33 * v->z;
  *w = r;
} /*TransContra3d*/

void Trans3Point2d ( const trans3d *tr, const point2d *p, point2d *q )
{
  point2d r;

  r.x = tr->U0.a11*p->x + tr->U0.a12*p->y + tr->U0.a14;
  r.y = tr->U0.a21*p->x + tr->U0.a22*p->y + tr->U0.a24;
  *q = r;
} /*Trans3Point2d*/

/* ////////////////////////////////////////////////////////////////////////// */
void Trans3Shiftd ( trans3d *tr, double tx, double ty, double tz )
{
  tr->U0.a14 += tr->U0.a11*tx + tr->U0.a12*ty + tr->U0.a13*tz;
  tr->U0.a24 += tr->U0.a21*tx + tr->U0.a22*ty + tr->U0.a23*tz;
  tr->U0.a34 += tr->U0.a31*tx + tr->U0.a32*ty + tr->U0.a33*tz;
} /*Trans3Shiftd*/

void Trans3Mirrord ( trans3d *tr, vector3d *n )
{
  /* tr := tr * M(n) */
  vector3d h;
  double   s;
  byte     i;

  h = *n;
  NormalizeVector3d ( &h );
  for ( i = 0; i < 3; i++ ) {
    s = -2.0 * (h.x * tr->U1.a[i][0] + h.y * tr->U1.a[i][1] + h.z * tr->U1.a[i][2]);
    tr->U1.a[i][0] += s * h.x;
    tr->U1.a[i][1] += s * h.y;
    tr->U1.a[i][2] += s * h.z;
  }
  tr->U1.detsgn = (short)(-tr->U1.detsgn);
} /*Trans3Mirrord*/

void Trans3Rotd ( trans3d *tr, byte j, byte k, double angle )
{
  double c, s;
  double m, n;
  byte i;

  if ( !angle )
    return;
  c = cos ( angle );
  s = sin ( angle );
  for ( i = 0; i < 3; i++ ) {
    m = tr->U1.a[i][j-1];
    n = tr->U1.a[i][k-1];
    tr->U1.a[i][j-1] = c * m + s * n;
    tr->U1.a[i][k-1] = c * n - s * m;
  }
} /*Trans3Rotd*/

void Trans3RotVd ( trans3d *tr, vector3d *v, double angle )
{
  /* rotate by the angle, around the axis parallel to the v vector */
  /* v must be nonzero.                                            */
  vector3d vv, e, f;
  double c, s, t;
  byte j;
  trans3d tt;

  vv = *v;
  NormalizeVector3d ( &vv );
  s = sin ( angle );
  c = cos ( angle );
  for ( j = 0; j < 3; j++ ) {
    e.x = tr->U1.a[j][0];
    e.y = tr->U1.a[j][1];
    e.z = tr->U1.a[j][2];
    t = vv.x * e.x + vv.y * e.y + vv.z * e.z;
    CrossProduct3d ( &e, &vv, &f );
    t = (1.0 - c) * t;
    tt.U1.a[j][0] = t * vv.x + c * tr->U1.a[j][0] + s * f.x;
    tt.U1.a[j][1] = t * vv.y + c * tr->U1.a[j][1] + s * f.y;
    tt.U1.a[j][2] = t * vv.z + c * tr->U1.a[j][2] + s * f.z;
  }
  memcpy ( tr->U1.a, tt.U1.a, 12*sizeof(double) );
} /*Trans3RotVd*/

void Trans3Scaled ( trans3d *tr, double sx, double sy, double sz )
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
} /*Trans3Scaled*/

/* ////////////////////////////////////////////////////////////////////////// */
void IdentTrans3d ( trans3d *tr )
{
  /* tr := I */
  int i, j;

  for ( i = 0; i <= 2; i++ ) {
    for ( j = 0; j <= 3; j++ )
      tr->U1.a[i][j] = 0.0;
    tr->U1.a[i][i] = 1.0;
  }
  tr->U1.detsgn = 1;
} /*IdentTrans3d*/

void CompTrans3d ( trans3d *s, trans3d *t, trans3d *u )
{
  /* s := t * u */
  byte i, j, k;
  double d;

  for ( i = 0; i <= 2; i++ ) {
    for ( j = 0; j <= 3; j++ ) {
      d = 0.0;
      for ( k = 0; k <= 2; k++ )
	d += t->U1.a[i][k] * u->U1.a[k][j];
      s->U1.a[i][j] = d;
    }
    s->U1.a[i][3] += t->U1.a[i][3];
  }
  s->U1.detsgn = (short)(t->U1.detsgn * u->U1.detsgn);
} /*CompTrans3d*/

void GeneralAffineTrans3d ( trans3d *tr, vector3d *v1, vector3d *v2, vector3d *v3 )
{
  trans3d s, t;
  double d;

  t.U0.a11 = v1->x;  t.U0.a12 = v2->x;  t.U0.a13 = v3->x;
  t.U0.a21 = v1->y;  t.U0.a22 = v2->y;  t.U0.a23 = v3->y;
  t.U0.a31 = v1->z;  t.U0.a32 = v2->z;  t.U0.a33 = v3->z;

  t.U0.a14 = 0.0;  t.U0.a24 = 0.0;  t.U0.a34 = 0.0;

  d = det3d ( v1, v2, v3 );
  if (d > 0.0)
    t.U1.detsgn = 1;
  else if (d < 0.0)
    t.U1.detsgn = -1;
  else
    t.U1.detsgn = 0;
  CompTrans3d ( &s, &t, tr );
  *tr = s;
} /*GeneralAffineTrans3d*/

void ShiftTrans3d ( trans3d *tr, double tx, double ty, double tz )
{
  /* tr := T(tx,ty,tz) * tr */
  tr->U0.a14 += tx;
  tr->U0.a24 += ty;
  tr->U0.a34 += tz;
} /*ShiftTrans3d*/

void MirrorTrans3d ( trans3d *tr, vector3d *n )
{
  /* tr := M(n) * tr */
  vector3d h;
  double s;
  byte i;

  h = *n;
  NormalizeVector3d ( &h );
  for (i = 0; i <= 3; i++) {
    s = -2.0 * (h.x * tr->U1.a[0][i] + h.y * tr->U1.a[1][i] + h.z * tr->U1.a[2][i]);
    tr->U1.a[0][i] += s * h.x;
    tr->U1.a[1][i] += s * h.y;
    tr->U1.a[2][i] += s * h.z;
  }
  tr->U1.detsgn = (short)(-tr->U1.detsgn);
} /*MirrorTrans3d*/

void RotTrans3d ( trans3d *tr, byte j, byte k, double angle )
{
  double c, s;
  double m, n;
  byte i;

  if ( !angle )
    return;
  c = cos ( angle );
  s = sin ( angle );
  for (i = 0; i <= 3; i++) {
    m = tr->U1.a[j-1][i];
    n = tr->U1.a[k-1][i];
    tr->U1.a[j-1][i] = c * m - s * n;
    tr->U1.a[k-1][i] = s * m + c * n;
  }
} /*RotTrans3d*/

void EulerRotTrans3d ( trans3d *tr, double psi, double theta, double phi )
{
  /* tr := E(psi,theta,phi) * tr */
  RotTrans3d ( tr, 1, 2, phi );
  RotTrans3d ( tr, 2, 3, theta );
  RotTrans3d ( tr, 1, 2, psi );
} /*EulerRotTrans3d*/

void FindRotVEulerd ( const vector3d *v, double angle,
                      double *psi, double *theta, double *phi )
{
  /* find Euler angles of rotation by angle around the axis parallel */
  /* to the v vector. v must be nonzero.                             */
  double a, b;

  if (fabs(v->x) + fabs(v->y) > 0.0) {
    a = atan2(v->x, v->y);
    b = atan2(sqrt(v->x * v->x + v->y * v->y), v->z);
    CompEulerRotd ( 0.0, -b, -a, a, b, angle, psi, theta, phi );
  }
  else {
    *psi = 0.0;
    *theta = 0.0;
    if (v->z > 0.0)      *phi =  angle;
    else if (v->z < 0.0) *phi = -angle;
    else                 *phi = 0.0;
  }
} /*FindRotVEulerd*/

void RotVTrans3d ( trans3d *tr, vector3d *v, double angle )
{
  /* rotate by the angle, around the axis parallel to the v vector */
  /* v must be nonzero.                                            */
  vector3d vv, e, f;
  double c, s, t;
  byte j;
  trans3d tt;

  vv = *v;
  NormalizeVector3d ( &vv );
  s = sin ( angle );
  c = cos ( angle );
  for (j = 0; j <= 3; j++) {
    e.x = tr->U1.a[0][j];
    e.y = tr->U1.a[1][j];
    e.z = tr->U1.a[2][j];
    t = vv.x * e.x + vv.y * e.y + vv.z * e.z;
    CrossProduct3d ( &vv, &e, &f );
    t = (1.0 - c) * t;
    tt.U1.a[0][j] = t * vv.x + c * tr->U1.a[0][j] + s * f.x;
    tt.U1.a[1][j] = t * vv.y + c * tr->U1.a[1][j] + s * f.y;
    tt.U1.a[2][j] = t * vv.z + c * tr->U1.a[2][j] + s * f.z;
  }
  memcpy ( tr->U1.a, tt.U1.a, 12*sizeof(double) );
} /*RotVTrans3d*/

void ScaleTrans3d ( trans3d *tr, double sx, double sy, double sz )
{
  /* tr := S(sx,sy,sz) * tr */
  byte i;
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
} /*ScaleTrans3d*/

boolean InvertTrans3d ( trans3d *tr )
{
  /* if tr not singular then tr := inversion(tr) */
  double b[16];

  memcpy ( b, &tr->U0.a11, 12*sizeof(double) );
  b[12] = b[13] = b[14] = 0.0;  b[15] = 1.0;
  if ( !pkn_GaussInvertMatrixd ( 4, b ) )
    return false;
  memcpy ( &tr->U0.a11, b, 12*sizeof(double) );
  return true;
} /*InvertTrans3d*/

double TrimAngled ( double angle )
{
  /* result := angle between -pi and pi */
  long a;

  if (angle > PI) {
    a = ((long)(angle / PI) + 1) / 2;
    angle -= a * 2.0 * PI;
    return angle;
  }
  if (angle < -PI) {
    a = ((long)(-(angle / PI)) + 1) / 2;
    angle += a * 2.0 * PI;
  }
  return angle;
} /*TrimAngled*/

void CompEulerRotd ( double psi1, double theta1, double phi1,
                     double psi2, double theta2, double phi2,
                     double *psi, double *theta, double *phi )
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
  *psi = TrimAngled ( _psi );
  *theta = TrimAngled ( _theta );
  *phi = TrimAngled ( _phi );
#undef eps
} /*CompEulerRotd*/

void CompRotV3d ( const vector3d *v1, double a1, const vector3d *v2, double a2,
                  vector3d *v, double *a )
{
  double c1, s1, c2, s2, c, s;
  vector3d w1, w2, w;

  w1 = *v1;  NormalizeVector3d ( &w1 );
  c1 = cos ( 0.5*a1 );  s1 = sin ( 0.5*a1 );
  w2 = *v2;  NormalizeVector3d ( &w2 );
  c2 = cos ( 0.5*a2 );  s2 = sin ( 0.5*a2 );
  c = c1*c2 - s1*s2*DotProduct3d ( &w1, &w2 );
  CrossProduct3d ( &w2, &w1, &w );
  MultVector3d ( s1*s2, &w, v );
  AddVector3Md ( v, &w2, s2*c1, v );
  AddVector3Md ( v, &w1, s1*c2, v );
  s = sqrt ( DotProduct3d ( v, v ) );
  *a = 2.0*atan2 ( s, c );
  if ( s )
    MultVector3d ( 1.0/s, v, v );
} /*CompRotV3d*/

void MultVector3d ( double a, const vector3d *v, vector3d *w )
{
  w->x = a * v->x;
  w->y = a * v->y;
  w->z = a * v->z;
} /*MultVector3d*/

void SubtractPoints3d ( const point3d *p1, const point3d *p2, vector3d *v )
{
  /* v := p1 - p2 */
  v->x = p1->x - p2->x;
  v->y = p1->y - p2->y;
  v->z = p1->z - p2->z;
} /*SubtractPoints3d*/

void AddVector3d ( const point3d *p, const vector3d *v, point3d *q )
{
  /* q := p + v */
  q->x = p->x + v->x;
  q->y = p->y + v->y;
  q->z = p->z + v->z;
} /*AddVector3d*/

void AddVector3Md ( const point3d *p, const vector3d *v, double t, point3d *q )
{
  /* q := p + t * v */
  q->x = p->x + t * v->x;
  q->y = p->y + t * v->y;
  q->z = p->z + t * v->z;
} /*AddVector3Md*/

void InterPoint3d ( const point3d *p1, const point3d *p2, double t, point3d *q )
{
  double s;

  s = 1.0 - t;
  q->x = s * p1->x + t * p2->x;
  q->y = s * p1->y + t * p2->y;
  q->z = s * p1->z + t * p2->z;
} /*InterPoint3d*/

void MidPoint3d ( const point3d *p1, const point3d *p2, point3d *q )
{
  q->x = 0.5 * (p1->x + p2->x);
  q->y = 0.5 * (p1->y + p2->y);
  q->z = 0.5 * (p1->z + p2->z);
} /*MidPoint3d*/

void Interp3Vectors3d ( const vector3d *p0, const vector3d *p1, const vector3d *p2,
                        const double *coeff, vector3d *p )
{
  p->x = coeff[0]*p0->x + coeff[1]*p1->x + coeff[2]*p2->x;
  p->y = coeff[0]*p0->y + coeff[1]*p1->y + coeff[2]*p2->y;
  p->z = coeff[0]*p0->z + coeff[1]*p1->z + coeff[2]*p2->z;
} /*Interp3Vectors3f*/

void NormalizeVector3d ( vector3d *v )
{
  /* if v <> 0 then v := v/length(v) */
#define tol 1.0e-15
  double r;

  r = sqrt ( v->x * v->x + v->y * v->y + v->z * v->z );
  if (r <= tol)
    return;
  r = 1.0 / r;
  v->x *= r;
  v->y *= r;
  v->z *= r;
#undef tol
} /*NormalizeVector3d*/

double DotProduct3d ( const vector3d *v1, const vector3d *v2 )
{
  return (v1->x * v2->x + v1->y * v2->y + v1->z * v2->z);
} /*DotProduct3d*/

void CrossProduct3d ( const vector3d *v1, const vector3d *v2, vector3d *v )
{
  vector3d vv;

  vv.x = v1->y * v2->z - v1->z * v2->y;
  vv.y = v1->z * v2->x - v1->x * v2->z;
  vv.z = v1->x * v2->y - v1->y * v2->x;
  *v = vv;
} /*CrossProduct3d*/

void AddCrossProduct3d ( const vector3d *w, const vector3d *v1, const vector3d *v2,
                         vector3d *v )
{
  vector3d vv;

  vv.x = v1->y * v2->z - v1->z * v2->y;
  vv.y = v1->z * v2->x - v1->x * v2->z;
  vv.z = v1->x * v2->y - v1->y * v2->x;
  AddVector3d ( w, &vv, v );
} /*AddCrossProduct3d*/

double det3d ( const vector3d *v1, const vector3d *v2, const vector3d *v3 )
{
  double a[9];

  memcpy ( &a[0], v1, sizeof(vector3d) );
  memcpy ( &a[3], v2, sizeof(vector3d) );
  memcpy ( &a[6], v3, sizeof(vector3d) );
  return pkn_detd ( 3, a );
}  /*det3d*/

void OrtVector3d ( const vector3d *v1, const vector3d *v2, vector3d *v )
{
#define tol 1.0e-15
  double nv, t;

  nv = DotProduct3d ( v1, v1 );
  if (nv < tol)
    *v = *v2;
  else {
    t = DotProduct3d ( v1, v2 ) / nv;
    SetVector3d ( v, v2->x - t * v1->x, v2->y - t * v1->y,
                  v2->z - t * v1->z );
  }
#undef tol
} /*OrtVector3d*/

void ProjectPointOnLine3d ( const point3d *p0, const point3d *p1,
                            point3d *q )
{
  vector3d v, w;
  double   d, e;

  SubtractPoints3d ( p1, p0, &v );
  SubtractPoints3d ( q, p0, &w );
  d = DotProduct3d ( &v, &v );
  e = DotProduct3d ( &v, &w );
  AddVector3Md ( p0, &v, e/d, q );
} /*ProjectPointOnLine3d*/

void ProjectPointOnPlane3d ( const point3d *p0, const point3d *p1, const point3d *p2,
                             point3d *q )
{
  vector3d v1, v2, w;
  double   d, e;

  SubtractPoints3d ( p1, p0, &v1 );
  SubtractPoints3d ( p2, p0, &v2 );
  OrtVector3d ( &v1, &v2, &v2 );
  SubtractPoints3d ( q, p0, &w );
  d = DotProduct3d ( &v1, &v1 );
  e = DotProduct3d ( &v1, &w );
  AddVector3Md ( p0, &v1, e/d, q );
  d = DotProduct3d ( &v2, &v2 );
  e = DotProduct3d ( &v2, &w );
  AddVector3Md ( q, &v2, e/d, q );
} /*ProjectPointOnPlane3d*/

double Point3Distanced ( point3d *p, point3d *q )
{
  vector3d v;

  SubtractPoints3d ( p, q, &v );
  return sqrt ( DotProduct3d ( &v, &v ) );
} /*Point3Distanced*/

