
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2015                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <pthread.h>

#define CONST_

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"
#include "raybez.h"

#include "raybezprivatef.h"

/* ///////////////////////////////////////////////////////////////////////// */
boolean _rbez_SubdividePatch2f ( int n, int m, point2f *p, point2f *q )
{
  int      i, j, k;
  vector2f v;
  float    a, lu, lv;

        /* choose the division direction */
  lu = lv = 0.0;
  for ( i = k = 0;  i <= n;  i++, k++ )
    for ( j = 0;  j < m;  j++, k++ ) {
      SubtractPoints2f ( &p[k+1], &p[k], &v );
      a = DotProduct2f ( &v, &v );
      if ( a > lv )
        lv = a;
    }
  for ( i = k = 0;  i < n;  i++ )
    for ( j = 0;  j <= m;  j++, k++ ) {
      SubtractPoints2f ( &p[k+m+1], &p[k], &v );
      a = DotProduct2f ( &v, &v );
      if ( a > lu )
        lu = a;
    }
        /* divide */
  if ( (float)(n*n)*lu > (float)(m*m)*lv ) {
    mbs_BisectBP2uf ( n, m, p, q );
    return true;
  }
  else {
    mbs_BisectBP2vf ( n, m, p, q );
    return false;
  }
} /*_rbez_SubdividePatch2f*/

boolean _rbez_ConvexHullTest2f ( int ncp, point2f *mcp )
{
  int          i;
  ConvexCone2f cone;

  cone.code = 0;
  for ( i = 0; i < ncp; i++ )
    if ( ExtendConvexCone2fv ( &cone, &mcp[i] ) == 2 )
      return true;
  return false;
} /*_rbez_ConvexHullTest2f*/

boolean _rbez_UniquenessTest2f ( int n, int m, int ncp, point2f *mcp,
                                 point2f *p, vector2f *du, vector2f *dv,
                                 float *K1, float *K2,
                                 float *workspace )
{
#define eps 1.0e-10
  vector2f v;
  float    c00, c01, c10, c11;
  float    d, e, f;
  point2f  *g, *gp, *gq;
  int      i, j;
  float    k1, k2, k;

  _mbs_BCHornerDerP2f ( n, m, mcp, 0.5, 0.5, p, du, dv, workspace );
                              /* compute inversion of differential matrix */
  d = du->x*dv->y - du->y*dv->x;
  e = fabs(du->x);
  f = fabs(du->y);  e = max ( e, f );
  f = fabs(dv->x);  e = max ( e, f );
  f = fabs(dv->y);  e = max ( e, f );
  if ( fabs(d) < eps*e )
    return false;
  c00 =  dv->y/d;  c01 = -dv->x/d;
  c10 = -du->y/d;  c11 =  du->x/d;
                              /* compose the mappings */
  g = &mcp[ncp];
  for ( i = 0; i < ncp; i++ ) {
    g[i].x = c00*mcp[i].x + c01*mcp[i].y;
    g[i].y = c10*mcp[i].x + c11*mcp[i].y;
  }
                              /* compute the constants K1, K2 */
  k1 = 0.0;
  for ( i = 0, gp = g, gq = gp+(m+1);
        i < n*(m+1);
        i++, gp++, gq++ ) {
    v.x = (float)n*(gq->x - gp->x) - 1.0;  v.y = (float)n*(gq->y - gp->y);
    k = v.x*v.x+v.y*v.y;
    if ( k > k1 ) {
      if ( k >= 1.0 ) return false;
      k1 = k;
    }
  }
  k2 = 0.0;
  for ( i = 0, gp = g;
        i <= n;
        i++, gp++ )
    for ( j = 0; j < m; j++, gp++ ) {
      v.x = (float)m*(gp[1].x - gp->x);  v.y = (float)m*(gp[1].y - gp->y) - 1.0;
      k = v.x*v.x+v.y*v.y;
      if ( k > k2 ) {
        if ( k+k1 >= 1.0 ) return false;
        k2 = k;
      }
    }
  *K1 = k1;  *K2 = k2;
  return true;
#undef eps
} /*_rbez_UniquenessTest2f*/

char _rbez_NewtonMethod2f ( int n, int m, point2f *mcp,
                            point2f *p, vector2f *pu, vector2f *pv,
                            point2f *z, float *workspace )
{
#define MAXITER 7
#define EPS     1.0e-6
#define DELTA   1.0e-6
  float    s, a[4];
  int      i, P[2], Q[2];

  z->x = z->y = 0.5;
  for ( i = 0; i < MAXITER; i++ ) {
    if ( i ) {
      if ( !_mbs_BCHornerDerP2f ( n, m, mcp, z->x, z->y,
                                  p, pu, pv, workspace ) )
        return RBEZ_NEWTON_ERROR;
    }
    if ( p->x*p->x+p->y*p->y < EPS*EPS )
      return RBEZ_NEWTON_YES;
    a[0] = pu->x;  a[1] = pv->x;  a[2] = pu->y;  a[3] = pv->y;
    if ( !pkn_GaussDecomposePLUQf ( 2, a, P, Q ) )
      return RBEZ_NEWTON_NO;
    pkn_multiSolvePLUQf ( 2, a, P, Q, 1, 1, &p->x );
    s = p->x*p->x+p->y*p->y;
    if ( s > 1.0 )
      return RBEZ_NEWTON_NO;
    z->x -= p->x;
    z->y -= p->y;
    if ( s < DELTA*DELTA )
      return RBEZ_NEWTON_YES;
  }
  return RBEZ_NEWTON_NO;
#undef EPS
#undef MAXITER
} /*_rbez_NewtonMethod2f*/

boolean _rbez_SecondTest2f ( point2f *z, int n, int m, float K1, float K2 )
{
  float r1, r2, R, R1, R2;

  if ( z->x < 0.0 ) r1 = 1.0 - 2.0*z->x;
  else if ( z->x > 1.0 ) r1 = 2.0*z->x - 1.0;
  else r1 = 1.0;
  if ( z->y < 0.0 ) r2 = 1.0 - 2.0*z->y;
  else if ( z->y > 1.0 ) r2 = 2.0*z->y - 1.0;
  else r2 = 1.0;

  R = pkv_rpower ( r1, n-1 ) * pkv_rpower ( r2, m-1 );
  R1 = R*r2;  R2 = R*r1;
  return K1*R1*R1 + K2*R2*R2 >= 1.0;
} /*_rbez_SecondTest2f*/

