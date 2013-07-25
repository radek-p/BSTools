
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* Header file for the libgeom library of C procedures - */
/* basic 2D, 3D and 4D affine geometry                   */ 

#ifndef PKGEOMD_H
#define PKGEOMD_H

#ifndef PKVARIA_H
#include "pkvaria.h"
#endif

#ifdef __cplusplus   
extern "C" {
#endif

typedef union trans2d {
  struct {
    /* affine transformation of a 2D space  */
    double a11, a12, a13,   /* coefficients of the transformation   */
           a21, a22, a23;
  } U0;
  struct {
    double a[2][3];
    short detsgn;   /* sign of the determinant : -1, 0 or 1 */
  } U1;
} trans2d;

typedef struct Box2d {
    double x0, x1, y0, y1;
  } Box2d;

typedef union trans3d {
  struct {
    /* affine transformation of a 3D space  */
    double a11, a12, a13, a14,   /* coefficients of the transformation   */
           a21, a22, a23, a24,
           a31, a32, a33, a34;
  } U0;
  struct {
    double a[3][4];
    short detsgn;   /* sign of the determinant : -1, 0 or 1 */
  } U1;
} trans3d;

typedef struct Box3d {
    double x0, x1, y0, y1, z0, z1;
  } Box3d;

typedef struct ray3d {
  /* ray, ie. a halfline in 3D space */
  point3d  p;  /* origin point of the ray */
  vector3d v;  /* direction               */
} ray3d;

typedef vector4d quaterniond;

typedef struct {
    double amin, amax;
    short  code;
  } ConvexCone2d;

/* /////////////////////////////////////////////////////////////////////////// */

void SetPoint2d ( point2d *p, double x, double y );
#define SetVector2d(v,x,y) SetPoint2d ( v, x, y )

void TransPoint2d ( const trans2d *tr, const point2d *p, point2d *q );
void TransVector2d ( const trans2d *tr, const vector2d *p, vector2d *q );
void IdentTrans2d ( trans2d *tr );
void CompTrans2d ( trans2d *s, trans2d *t, trans2d *u );
void ShiftTrans2d ( trans2d *tr, double tx, double ty );
void RotTrans2d ( trans2d *tr, double angle );
void ScaleTrans2d ( trans2d *tr, double sx, double sy );
void Trans2Shiftd ( trans2d *tr, double tx, double ty );
void Trans2Rotd ( trans2d *tr, double angle );
void Trans2Scaled ( trans2d *tr, double sx, double sy );
boolean InvertTrans2d ( trans2d *tr );

void MultVector2d ( double a, const vector2d *v, vector2d *w );
void SubtractPoints2d ( const point2d *p1, const point2d *p2, vector2d *v );
void AddVector2d ( const point2d *p, const vector2d *v, point2d *q );
void AddVector2Md ( const point2d *p, const vector2d *v, double t, point2d *q );
void InterPoint2d ( const point2d *p1, const point2d *p2, double t, point2d *q );
void MidPoint2d ( const point2d *p1, const point2d *p2, point2d *q );
void Interp3Vectors2d ( const vector2d *p0, const vector2d *p1, const vector2d *p2,
                        const double *coeff, vector2d *p );
void NormalizeVector2d ( vector2d *v );

double DotProduct2d ( const vector2d *v1, const vector2d *v2 );
double det2d ( const vector2d *v1, const vector2d *v2 );
void OrtVector2d ( const vector2d *v1, const vector2d *v2, vector2d *v );

void ProjectPointOnLine2d ( const point2d *p0, const point2d *p1,
                            point2d *q );

double Point2Distanced ( point2d *p, point2d *q );

short ExtendConvexCone2d ( ConvexCone2d *cone, double a );
short ExtendConvexCone2dv ( ConvexCone2d *cone, vector2d *v );
boolean InsideConvexCone2d ( ConvexCone2d *cone, double a );

/* /////////////////////////////////////////////////////////////////////////// */
void Trans2Point3d ( const trans2d *tr, const point3d *p, point3d *q );

/* /////////////////////////////////////////////////////////////////////////// */

void SetPoint3d ( point3d *p, double x, double y, double z );
#define SetVector3d(v,x,y,z) SetPoint3d ( v, x, y, z )

void Point3to2d ( const point3d *P, point2d *p );
void Point2to3d ( const point2d *p, double w, point3d *P );
void Point4to2d ( const point4d *P, point2d *p );
void Point2to4d ( const point2d *p, float w, point4d *P );

void TransPoint3d ( const trans3d *tr, const point3d *p, point3d *q );
void TransVector3d ( const trans3d *tr, const vector3d *v, vector3d *w );
void TransContra3d ( const trans3d *tri, const vector3d *v, vector3d *w );
void Trans3Point2d ( const trans3d *tr, const point2d *p, point2d *q );
void Trans3Point4d ( const trans3d *tr, const point4d *p, point4d *q );

void Trans3Shiftd ( trans3d *tr, double tx, double ty, double tz );
void Trans3Mirrord ( trans3d *tr, vector3d *n );
void Trans3Rotd ( trans3d *tr, byte j, byte k, double angle );
#define Trans3RotXd(tr,angle) Trans3Rotd ( tr, 2, 3, angle )
#define Trans3RotYd(tr,angle) Trans3Rotd ( tr, 3, 1, angle )
#define Trans3RotZd(tr,angle) Trans3Rotd ( tr, 1, 2, angle )
void Trans3RotVd ( trans3d *tr, vector3d *v, double angle );
void Trans3Scaled ( trans3d *tr, double sx, double sy, double sz );

void IdentTrans3d ( trans3d *tr );
void CompTrans3d ( trans3d *s, trans3d *t, trans3d *u );
void GeneralAffineTrans3d ( trans3d *tr, vector3d *v1, vector3d *v2,
                            vector3d *v3 );
void ShiftTrans3d ( trans3d *tr, double tx, double ty, double tz );
void MirrorTrans3d ( trans3d *tr, vector3d *n );

void RotTrans3d ( trans3d *tr, byte j, byte k, double angle );
void EulerRotTrans3d ( trans3d *tr, double psi, double theta, double phi );
#define RotXTrans3d(tr,angle) RotTrans3d ( tr, 2, 3, angle )
#define RotYTrans3d(tr,angle) RotTrans3d ( tr, 3, 1, angle )
#define RotZTrans3d(tr,angle) RotTrans3d ( tr, 1, 2, angle )

void FindRotVEulerd ( const vector3d *v, double angle,
                      double *psi, double *theta, double *phi );
void RotVTrans3d ( trans3d *tr, vector3d *v, double angle );
void ScaleTrans3d ( trans3d *tr, double sx, double sy, double sz );
boolean InvertTrans3d ( trans3d *tr );
double TrimAngled ( double angle );
void CompEulerRotd ( double psi1, double theta1, double phi1, double psi2,
                     double theta2, double phi2, double *psi, double *theta,
                     double *phi );
void CompRotV3d ( const vector3d *v1, double a1, const vector3d *v2, double a2,
                  vector3d *v, double *a );

void MultVector3d ( double a, const vector3d *v, vector3d *w );
void SubtractPoints3d ( const point3d *p1, const point3d *p2, vector3d *v );
void AddVector3d ( const point3d *p, const vector3d *v, point3d *q );
void AddVector3Md ( const point3d *p, const vector3d *v, double t,
                    point3d *q );
void InterPoint3d ( const point3d *p1, const point3d *p2, double t,
                    point3d *q );
void MidPoint3d ( const point3d *p1, const point3d *p2, point3d *q );
void Interp3Vectors3d ( const vector3d *p0, const vector3d *p1, const vector3d *p2,
                        const double *coeff, vector3d *p );
void NormalizeVector3d ( vector3d *v );

double DotProduct3d ( const vector3d *v1, const vector3d *v2 );
void CrossProduct3d ( const vector3d *v1, const vector3d *v2, vector3d *v );
double det3d ( const vector3d *v1, const vector3d *v2, const vector3d *v3 );
void OrtVector3d ( const vector3d *v1, const vector3d *v2, vector3d *v );

void ProjectPointOnLine3d ( const point3d *p0, const point3d *p1,
                            point3d *q );
void ProjectPointOnPlane3d ( const point3d *p0, const point3d *p1, const point3d *p2,
                             point3d *q );

double Point3Distanced ( point3d *p, point3d *q );

/* /////////////////////////////////////////////////////////////////////////// */
void SetPoint4d ( point4d *p, double X, double Y, double Z, double W );
#define SetVector4d(v,X,Y,Z,W) SetPoint4d ( v, X, Y, Z, W )

void Point4to3d ( const point4d *P, point3d *p );
void Point3to4d ( const point3d *p, double w, point4d *P );

void SubtractPoints4d ( const point4d *p1, const point4d *p2, vector4d *v );
void AddVector4d ( const point4d *p, const vector4d *v, point4d *q );
void AddVector4Md ( const point4d *p, const vector4d *v, double t,
                    point4d *q );
void MultVector4d ( double a, const vector4d *v, vector4d *w );

void InterPoint4d ( const point4d *p1, const point4d *p2, double t,
                    point4d *q );
void MidPoint4d ( const point4d *p1, const point4d *p2, point4d *q );
void Interp3Vectors4d ( const vector4d *p0, const vector4d *p1, const vector4d *p2,
                        const double *coeff, vector4d *p );
void NormalizeVector4d ( vector4d *v );

double det4d ( const vector4d *v0, const vector4d *v1,
               const vector4d *v2, const vector4d *v3 );
void CrossProduct4d (const vector4d *v0, const vector4d *v1,
                     const vector4d *v2, vector4d *v );
void CrossProduct4P3d ( const vector4d *v0, const vector4d *v1,
                        const vector4d *v2, vector3d *v );
double DotProduct4d ( const vector4d *v0, const vector4d *v1 );
void OutProduct4P3d ( const vector4d *v0, const vector4d *v1, vector3d *v );

void OrtVector4d ( const vector4d *v1, const vector4d *v2, vector4d *v );

void ProjectPointOnLine4d ( const point4d *p0, const point4d *p1,
                            point4d *q );
void ProjectPointOnPlane4d ( const point4d *p0, const point4d *p1, const point4d *p2,
                             point4d *q );

double pkg_LineRayDist3d ( ray3d *ray, point3d *q0, point3d *q1,
                           double *rparam, double *lparam,
                           point3d *rp, point3d *lp );

/* /////////////////////////////////////////////////////////////////////////// */
#define SetQuaterniond(q,a,x,y,z) SetVector4d ( q, x, y, z, a )
#define AddQuaternionsd AddVectors4d

void QuaternionMultd ( quaterniond *q1, quaterniond *q2, quaterniond *q );
void QuaternionMultInvd ( quaterniond *q1, quaterniond *q2, quaterniond *q );
void QuaternionInvMultd ( quaterniond *q1, quaterniond *q2, quaterniond *q );

void QuaternionRotTrans3d ( trans3d *tr, quaterniond *q );
void Trans3QuaternionRotd ( trans3d *tr, quaterniond *q );

#ifdef __cplusplus
}
#endif

#endif /*PKGEOMD_H*/

