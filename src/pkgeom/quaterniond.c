
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


void QuaternionMultd ( quaterniond *q1, quaterniond *q2, quaterniond *q )
{
/* q = q1*q2 */
  SetQuaterniond ( q,
                   q1->w*q2->w - q1->x*q2->x - q1->y*q2->y - q1->z*q2->z,
                   q1->x*q2->w + q1->w*q2->x - q1->z*q2->y + q1->y*q2->z,
                   q1->y*q2->w + q1->z*q2->x + q1->w*q2->y - q1->x*q2->z,
                   q1->z*q2->w - q1->y*q2->x + q1->x*q2->y + q1->w*q2->z );
} /*QuaternionMultd*/

void QuaternionMultInvd ( quaterniond *q1, quaterniond *q2, quaterniond *q )
{
/* q = q1/q2 ( = q1*q2^{-1} ) */
  quaterniond q2i;
  double      d;

  d = DotProduct4d ( q2, q2 );
  if ( d ) {
    SetQuaterniond ( &q2i, q2->w/d, -q2->x/d, -q2->y/d, -q2->z/d );
    QuaternionMultd ( q1, &q2i, q );
  }
} /*QuaternionMultInvd*/

void QuaternionInvMultd ( quaterniond *q1, quaterniond *q2, quaterniond *q )
{
/* q = q1\q2 ( = q1^{-1}*q2 ) */
  quaterniond q1i;
  double     d;

  d = DotProduct4d ( q1, q1 );
  if ( d ) {
    SetQuaterniond ( &q1i, q1->w/d, -q1->x/d, -q1->y/d, -q1->z/d );
    QuaternionMultd ( &q1i, q2, q );
  }
} /*QuaternionInvMultd*/

void QuaternionRotTrans3d ( trans3d *tr, quaterniond *q )
{
  vector3d v, e, f;
  double   d, c, s, c2, s2, t;
  byte     j;
  trans3d  tt;

  d = DotProduct4d ( q, q );
  if ( d > 1.0e-20 ) {
    d = sqrt ( d );
    memcpy ( &v, q, sizeof(vector3d) );
    c = q->w/d;
    s = sqrt ( DotProduct3d ( &v, &v ) )/d;
    if ( fabs ( s ) > 1.0e-9 ) {
      c2 = c*c - s*s;
      s2 = 2.0*c*s;
      NormalizeVector3d ( &v );
      for ( j = 0; j <= 3; j++ ) {
        e.x = tr->U1.a[0][j];   
        e.y = tr->U1.a[1][j];   
        e.z = tr->U1.a[2][j];   
        t = DotProduct3d ( &v, &e );
        CrossProduct3d ( &v, &e, &f );
        t = (1.0 - c2) * t;
        tt.U1.a[0][j] = t * v.x + c2 * tr->U1.a[0][j] + s2 * f.x;
        tt.U1.a[1][j] = t * v.y + c2 * tr->U1.a[1][j] + s2 * f.y;
        tt.U1.a[2][j] = t * v.z + c2 * tr->U1.a[2][j] + s2 * f.z;
      }
      memcpy ( tr->U1.a, tt.U1.a, 12*sizeof(double) );
    }
  }
} /*QuaternionRotTrans3d*/

void Trans3QuaternionRotd ( trans3d *tr, quaterniond *q )
{
  vector3d v, e, f;
  double   d, c, s, c2, s2, t;
  byte     i;
  trans3d  tt;

  d = DotProduct4d ( q, q );
  if ( d > 1.0e-20 ) {
    d = sqrt ( d );
    memcpy ( &v, q, sizeof(vector3d) );
    c = q->w/d;
    s = sqrt ( DotProduct3d ( &v, &v ) )/d;
    if ( fabs ( s ) > 1.0e-9 ) {
      c2 = c*c - s*s;
      s2 = 2.0*c*s;
      NormalizeVector3d ( &v );
      for ( i = 0; i < 3; i++ ) {
        e.x = tr->U1.a[i][0];   
        e.y = tr->U1.a[i][1];   
        e.z = tr->U1.a[i][2];   
        t = DotProduct3d ( &v, &e );
        CrossProduct3d ( &e, &v, &f );
        t = (1.0 - c2) * t;
        tt.U1.a[i][0] = t * v.x + c2 * tr->U1.a[i][0] + s2 * f.x;
        tt.U1.a[i][1] = t * v.y + c2 * tr->U1.a[i][1] + s2 * f.y;
        tt.U1.a[i][2] = t * v.z + c2 * tr->U1.a[i][2] + s2 * f.z;
      }
      memcpy ( tr->U1.a, tt.U1.a, 12*sizeof(double) );
    }
  }
} /*Trans3QuaternionRotd*/

