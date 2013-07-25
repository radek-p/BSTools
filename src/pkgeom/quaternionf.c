
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

#include "pknum.h"
#include "pkgeom.h"


void QuaternionMultf ( quaternionf *q1, quaternionf *q2, quaternionf *q )
{
/* q = q1*q2 */
  SetQuaternionf ( q,
                   q1->w*q2->w - q1->x*q2->x - q1->y*q2->y - q1->z*q2->z,
                   q1->x*q2->w + q1->w*q2->x - q1->z*q2->y + q1->y*q2->z,
                   q1->y*q2->w + q1->z*q2->x + q1->w*q2->y - q1->x*q2->z,
                   q1->z*q2->w - q1->y*q2->x + q1->x*q2->y + q1->w*q2->z );
} /*QuaternionMultf*/

void QuaternionMultInvf ( quaternionf *q1, quaternionf *q2, quaternionf *q )
{
/* q = q1/q2 ( = q1*q2^{-1} ) */
  quaternionf q2i;
  double      d;

  d = DotProduct4f ( q2, q2 );
  if ( d ) {
    SetQuaternionf ( &q2i, q2->w/d, -q2->x/d, -q2->y/d, -q2->z/d );
    QuaternionMultf ( q1, &q2i, q );
  }
} /*QuaternionMultInvf*/

void QuaternionInvMultf ( quaternionf *q1, quaternionf *q2, quaternionf *q )
{
/* q = q1\q2 ( = q1^{-1}*q2 ) */
  quaternionf q1i;
  double     d;

  d = DotProduct4f ( q1, q1 );
  if ( d ) {
    SetQuaternionf ( &q1i, q1->w/d, -q1->x/d, -q1->y/d, -q1->z/d );
    QuaternionMultf ( &q1i, q2, q );
  }
} /*QuaternionInvMultf*/

void QuaternionRotTrans3f ( trans3f *tr, quaternionf *q )
{
  vector3f v, e, f;
  double   d, c, s, c2, s2, t;
  byte     j;
  trans3f  tt;

  d = DotProduct4f ( q, q );
  if ( d > 1.0e-20 ) {
    d = sqrt ( d );
    memcpy ( &v, q, sizeof(vector3f) );
    c = q->w/d;
    s = sqrt ( DotProduct3f ( &v, &v ) )/d;
    if ( fabs ( s ) > 1.0e-9 ) {
      c2 = c*c - s*s;
      s2 = 2.0*c*s;
      NormalizeVector3f ( &v );
      for ( j = 0; j <= 3; j++ ) {
        e.x = tr->U1.a[0][j];   
        e.y = tr->U1.a[1][j];   
        e.z = tr->U1.a[2][j];   
        t = DotProduct3f ( &v, &e );
        CrossProduct3f ( &v, &e, &f );
        t = (1.0 - c2) * t;
        tt.U1.a[0][j] = t * v.x + c2 * tr->U1.a[0][j] + s2 * f.x;
        tt.U1.a[1][j] = t * v.y + c2 * tr->U1.a[1][j] + s2 * f.y;
        tt.U1.a[2][j] = t * v.z + c2 * tr->U1.a[2][j] + s2 * f.z;
      }
      memcpy ( tr->U1.a, tt.U1.a, 12*sizeof(float) );
    }
  }
} /*QuaternionRotTrans3f*/

void Trans3QuaternionRotf ( trans3f *tr, quaternionf *q )
{
  vector3f v, e, f;
  double   d, c, s, c2, s2, t;
  byte     i;
  trans3f  tt;

  d = DotProduct4f ( q, q );
  if ( d > 1.0e-20 ) {
    d = sqrt ( d );
    memcpy ( &v, q, sizeof(vector3f) );
    c = q->w/d;
    s = sqrt ( DotProduct3f ( &v, &v ) )/d;
    if ( fabs ( s ) > 1.0e-9 ) {
      c2 = c*c - s*s;
      s2 = 2.0*c*s;
      NormalizeVector3f ( &v );
      for ( i = 0; i < 3; i++ ) {
        e.x = tr->U1.a[i][0];   
        e.y = tr->U1.a[i][1];   
        e.z = tr->U1.a[i][2];   
        t = DotProduct3f ( &v, &e );
        CrossProduct3f ( &e, &v, &f );
        t = (1.0 - c2) * t;
        tt.U1.a[i][0] = t * v.x + c2 * tr->U1.a[i][0] + s2 * f.x;
        tt.U1.a[i][1] = t * v.y + c2 * tr->U1.a[i][1] + s2 * f.y;
        tt.U1.a[i][2] = t * v.z + c2 * tr->U1.a[i][2] + s2 * f.z;
      }
      memcpy ( tr->U1.a, tt.U1.a, 12*sizeof(float) );
    }
  }
} /*Trans3QuaternionRotf*/

