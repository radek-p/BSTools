
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* Header file for the libgeom library of C procedures - */
/* basic 2D, 3D and 4D affine geometry                   */ 

#ifndef PKGEOMF_H
#define PKGEOMF_H

#ifndef PKVARIA_H
#include "pkvaria.h"
#endif

#ifdef __cplusplus   
extern "C" {
#endif

typedef union trans2f {
  struct {
    /* affine transformation of a 2D space  */
    float a11, a12, a13,   /* coefficients of the transformation   */
          a21, a22, a23;
  } U0;
  struct {
    float a[2][3];
    short detsgn;   /* sign of the determinant : -1, 0 or 1 */
  } U1;
} trans2f;

typedef struct Box2f {
    float x0, x1, y0, y1;
  } Box2f;

typedef union trans3f {
  struct {
    /* affine transformation of a 3D space  */
    float a11, a12, a13, a14,   /* coefficients of the transformation   */
          a21, a22, a23, a24,
          a31, a32, a33, a34;
  } U0;
  struct {
    float a[3][4];
    short detsgn;   /* sign of the determinant : -1, 0 or 1 */
  } U1;
} trans3f;

typedef struct Box3f {
    float x0, x1, y0, y1, z0, z1;
  } Box3f;

typedef struct ray3f {
  /* ray, ie. a halfline in 3D space */
  point3f  p;  /* origin point of the ray */
  vector3f v;  /* direction               */
} ray3f;

typedef vector4f quaternionf;

typedef struct {
    float amin, amax;
    short code;  
  } ConvexCone2f;

/* /////////////////////////////////////////////////////////////////////////// */

void SetPoint2f ( point2f *p, float x, float y );
#define SetVector2f(v,x,y) SetPoint2f ( v, x, y )

void TransPoint2f ( const trans2f *tr, const point2f *p, point2f *q );
void TransVector2f ( const trans2f *tr, const vector2f *p, vector2f *q );
void IdentTrans2f ( trans2f *tr );
void CompTrans2f ( trans2f *s, trans2f *t, trans2f *u );
void ShiftTrans2f ( trans2f *tr, float tx, float ty );
void RotTrans2f ( trans2f *tr, float angle );
void ScaleTrans2f ( trans2f *tr, float sx, float sy );
void Trans2Shiftf ( trans2f *tr, float tx, float ty );
void Trans2Rotf ( trans2f *tr, float angle );
void Trans2Scalef ( trans2f *tr, float sx, float sy );
boolean InvertTrans2f ( trans2f *tr );

void MultVector2f ( double a, const vector2f *v, vector2f *w );
void SubtractPoints2f ( const point2f *p1, const point2f *p2, vector2f *v );
void AddVector2f ( const point2f *p, const vector2f *v, point2f *q );
void AddVector2Mf ( const point2f *p, const vector2f *v, double t, point2f *q );
void InterPoint2f ( const point2f *p1, const point2f *p2, double t, point2f *q );
void Interp3Vectors2f ( const vector2f *p0, const vector2f *p1, const vector2f *p2,
                        const float *coeff, vector2f *p );
void NormalizeVector2f ( vector2f *v );
void MidPoint2f ( const point2f *p1, const point2f *p2, point2f *q );

double DotProduct2f ( const vector2f *v1, const vector2f *v2 );
double det2f ( const vector2f *v1, const vector2f *v2 );
void OrtVector2f ( const vector2f *v1, const vector2f *v2, vector2f *v );

void ProjectPointOnLine2f ( const point2f *p0, const point2f *p1,
                            point2f *q );

double Point2Distancef ( point2f *p, point2f *q );

short ExtendConvexCone2f ( ConvexCone2f *cone, float a );
short ExtendConvexCone2fv ( ConvexCone2f *cone, vector2f *v );
boolean InsideConvexCone2f ( ConvexCone2f *cone, float a );

/* /////////////////////////////////////////////////////////////////////////// */
void Trans2Point3f ( const trans2f *tr, const point3f *p, point3f *q );

/* /////////////////////////////////////////////////////////////////////////// */

void SetPoint3f ( point3f *p, float x, float y, float z );
#define SetVector3f(v,x,y,z) SetPoint3f ( v, x, y, z )

void Point3to2f ( const point3f *P, point2f *p );
void Point2to3f ( const point2f *p, float w, point3f *P );
void Point4to2f ( const point4f *P, point2f *p );
void Point2to4f ( const point2f *p, float w, point4f *P );

void TransPoint3f ( const trans3f *tr, const point3f *p, point3f *q );
void TransVector3f ( const trans3f *tr, const vector3f *v, vector3f *w );
void TransContra3f ( const trans3f *tri, const vector3f *v, vector3f *w );
void Trans3Point2f ( const trans3f *tr, const point2f *p, point2f *q );
void Trans3Point4f ( const trans3f *tr, const point4f *p, point4f *q );

void Trans3Shiftf ( trans3f *tr, float tx, float ty, float tz );
void Trans3Mirrorf ( trans3f *tr, vector3f *n );
void Trans3Rotf ( trans3f *tr, byte j, byte k, float angle );
#define Trans3RotXf(tr,angle) Trans3Rotf ( tr, 2, 3, angle )
#define Trans3RotYTf(tr,angle) Trans3Rotf ( tr, 3, 1, angle )
#define Trans3RotZTf(tr,angle) Trans3Rotf ( tr, 1, 2, angle )
void Trans3RotVf ( trans3f *tr, vector3f *v, float angle );
void Trans3Scalef ( trans3f *tr, float sx, float sy, float sz );

void IdentTrans3f ( trans3f *tr );
void CompTrans3f ( trans3f *s, trans3f *t, trans3f *u );
void GeneralAffineTrans3f ( trans3f *tr, vector3f *v1, vector3f *v2,
                            vector3f *v3 );
void ShiftTrans3f ( trans3f *tr, float tx, float ty, float tz );
void MirrorTrans3f ( trans3f *tr, vector3f *n );

void RotTrans3f ( trans3f *tr, byte j, byte k, float angle );
void EulerRotTrans3f ( trans3f *tr, float psi, float theta, float phi );
#define RotXTrans3f(tr,angle) RotTrans3f ( tr, 2, 3, angle )
#define RotYTrans3f(tr,angle) RotTrans3f ( tr, 3, 1, angle )
#define RotZTrans3f(tr,angle) RotTrans3f ( tr, 1, 2, angle )

void FindRotVEulerf ( const vector3f *v, float angle,
                      float *psi, float *theta, float *phi );
void RotVTrans3f ( trans3f *tr, vector3f *v, float angle );
void ScaleTrans3f ( trans3f *tr, float sx, float sy, float sz );
boolean InvertTrans3f ( trans3f *tr );
float TrimAnglef ( float angle );
void CompEulerRotf ( float psi1, float theta1, float phi1, float psi2,
                     float theta2, float phi2, float *psi, float *theta,
                     float *phi );
void CompRotV3f ( const vector3f *v1, float a1, const vector3f *v2, float a2,
                  vector3f *v, float *a );

void MultVector3f ( double a, const vector3f *v, vector3f *w );
void SubtractPoints3f ( const point3f *p1, const point3f *p2, vector3f *v );
void AddVector3f ( const point3f *p, const vector3f *v, point3f *q );
void AddVector3Mf ( const point3f *p, const vector3f *v, double t,
                    point3f *q );
void InterPoint3f ( const point3f *p1, const point3f *p2, double t,
                    point3f *q );
void MidPoint3f ( const point3f *p1, const point3f *p2, point3f *q );
void Interp3Vectors3f ( const vector3f *p0, const vector3f *p1, const vector3f *p2,
                        const float *coeff, vector3f *p );
void NormalizeVector3f ( vector3f *v );

double DotProduct3f ( const vector3f *v1, const vector3f *v2 );
void CrossProduct3f ( const vector3f *v1, const vector3f *v2, vector3f *v );
double det3f ( const vector3f *v1, const vector3f *v2, const vector3f *v3 );
void OrtVector3f ( const vector3f *v1, const vector3f *v2, vector3f *v );

void ProjectPointOnLine3f ( const point3f *p0, const point3f *p1,
                            point3f *q );
void ProjectPointOnPlane3f ( const point3f *p0, const point3f *p1, const point3f *p2,
                             point3f *q );

double Point3Distancef ( point3f *p, point3f *q );

/* /////////////////////////////////////////////////////////////////////////// */
void SetPoint4f ( point4f *p, float X, float Y, float Z, float W );
#define SetVector4f(v,X,Y,Z,W) SetPoint4f ( v, X, Y, Z, W )

void Point4to3f ( const point4f *P, point3f *p );
void Point3to4f ( const point3f *p, float w, point4f *P );

void SubtractPoints4f ( const point4f *p1, const point4f *p2, vector4f *v );
void AddVector4f ( const point4f *p, const vector4f *v, point4f *q );
void AddVector4Mf ( const point4f *p, const vector4f *v, double t,
                    point4f *q );
void MultVector4f ( double a, const vector4f *v, vector4f *w );

void InterPoint4f ( const point4f *p1, const point4f *p2, double t,
                    point4f *q );
void MidPoint4f ( const point4f *p1, const point4f *p2, point4f *q );
void Interp3Vectors4f ( const vector4f *p0, const vector4f *p1, const vector4f *p2,
                        const float *coeff, vector4f *p );
void NormalizeVector4f ( vector4f *v );

double det4f ( const vector4f *v0, const vector4f *v1,
               const vector4f *v2, const vector4f *v3 );
void CrossProduct4f (const vector4f *v0, const vector4f *v1,
                     const vector4f *v2, vector4f *v );
void CrossProduct4P3f ( const vector4f *v0, const vector4f *v1,
                        const vector4f *v2, vector3f *v );
double DotProduct4f ( const vector4f *v0, const vector4f *v1 );
void OutProduct4P3f ( const vector4f *v0, const vector4f *v1, vector3f *v );

void OrtVector4f ( const vector4f *v1, const vector4f *v2, vector4f *v );

void ProjectPointOnLine4f ( const point4f *p0, const point4f *p1,
                            point4f *q );
void ProjectPointOnPlane4f ( const point4f *p0, const point4f *p1, const point4f *p2,
                             point4f *q );

float pkg_LineRayDist3f ( ray3f *ray, point3f *q0, point3f *q1,
                          float *rparam, float *lparam,
                          point3f *rp, point3f *lp );

/* /////////////////////////////////////////////////////////////////////////// */
#define SetQuaternionf(q,a,x,y,z) SetVector4f ( q, x, y, z, a )
#define AddQuaternionsf AddVectors4f

void QuaternionMultf ( quaternionf *q1, quaternionf *q2, quaternionf *q );
void QuaternionMultInvf ( quaternionf *q1, quaternionf *q2, quaternionf *q );
void QuaternionInvMultf ( quaternionf *q1, quaternionf *q2, quaternionf *q );

void QuaternionRotTrans3f ( trans3f *tr, quaternionf *q );
void Trans3QuaternionRotf ( trans3f *tr, quaternionf *q );

#ifdef __cplusplus
}
#endif

#endif /*PKGEOMF_H*/

