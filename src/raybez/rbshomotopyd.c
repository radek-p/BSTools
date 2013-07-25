
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2013                                  */
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

/* ///////////////////////////////////////////////////////////////////////// */
#define MAXLEVEL 30

typedef struct {
    double  u0, u1, v0, v1, t0, t1;
    int     level;
  } stack_el;

static boolean _rbez_ConvexHullTestd ( int degree,
                                       point3d *bcp0, point3d *bcp1,
                                       point3d *bcq0, point3d *bcq1 )
{
  vector3d v[4] = {{1.0,1.0,1.0},{1.0,1.0,-1.0},{1.0,-1.0,1.0},{1.0,-1.0,-1.0}};
  vector3d w, z;
  int      i, j, k;

  SubtractPoints3d ( &bcp0[0], &bcq0[0], &z );
  if ( z.x <= 0.0 ) {
    for ( i = 0; i <= degree; i++ )
      for ( j = 0; j <= degree; j++ ) {
        if ( bcp0[i].x-bcq0[j].x >= 0.0 ) goto next1;
        if ( bcp1[i].x-bcq0[j].x >= 0.0 ) goto next1;
        if ( bcp0[i].x-bcq1[j].x >= 0.0 ) goto next1;
        if ( bcp1[i].x-bcq1[j].x >= 0.0 ) goto next1;
      }
  }
  else {
    for ( i = 0; i <= degree; i++ )
      for ( j = 0; j <= degree; j++ ) {
        if ( bcp0[i].x-bcq0[j].x <= 0.0 ) goto next1;
        if ( bcp1[i].x-bcq0[j].x <= 0.0 ) goto next1;
        if ( bcp0[i].x-bcq1[j].x <= 0.0 ) goto next1;
        if ( bcp1[i].x-bcq1[j].x <= 0.0 ) goto next1;
      }
  }
  return false;
next1:
  if ( z.y <= 0.0 ) {
    for ( i = 0; i <= degree; i++ )
      for ( j = 0; j <= degree; j++ ) {
        if ( bcp0[i].y-bcq0[j].y >= 0.0 ) goto next2;
        if ( bcp1[i].y-bcq0[j].y >= 0.0 ) goto next2;
        if ( bcp0[i].y-bcq1[j].y >= 0.0 ) goto next2;
        if ( bcp1[i].y-bcq1[j].y >= 0.0 ) goto next2;
      }
  }
  else {
    for ( i = 0; i <= degree; i++ )
      for ( j = 0; j <= degree; j++ ) {
        if ( bcp0[i].y-bcq0[j].y <= 0.0 ) goto next2;
        if ( bcp1[i].y-bcq0[j].y <= 0.0 ) goto next2;
        if ( bcp0[i].y-bcq1[j].y <= 0.0 ) goto next2;
        if ( bcp1[i].y-bcq1[j].y <= 0.0 ) goto next2;
      }
  }
  return false;
next2:
  if ( z.z <= 0.0 ) {
    for ( i = 0; i <= degree; i++ )
      for ( j = 0; j <= degree; j++ ) {
        if ( bcp0[i].z-bcq0[j].z >= 0.0 ) goto next3;
        if ( bcp1[i].z-bcq0[j].z >= 0.0 ) goto next3;
        if ( bcp0[i].z-bcq1[j].z >= 0.0 ) goto next3;
        if ( bcp1[i].z-bcq1[j].z >= 0.0 ) goto next3;
      }
  }
  else {
    for ( i = 0; i <= degree; i++ )
      for ( j = 0; j <= degree; j++ ) {
        if ( bcp0[i].z-bcq0[j].z <= 0.0 ) goto next3;
        if ( bcp1[i].z-bcq0[j].z <= 0.0 ) goto next3;
        if ( bcp0[i].z-bcq1[j].z <= 0.0 ) goto next3;
        if ( bcp1[i].z-bcq1[j].z <= 0.0 ) goto next3;
      }
  }
  return false;
next3:
  for ( k = 0; k < 4; k++ ) {
    if ( DotProduct3d ( &v[k], &z ) <= 0.0 ) {
      for ( i = 0; i <= degree; i++ )
        for ( j = 0; j <= degree; j++ ) {
          SubtractPoints3d ( &bcp0[i], &bcq0[j], &w );
          if ( DotProduct3d ( &v[k], &w ) >= 0.0 ) goto next4;
          SubtractPoints3d ( &bcp1[i], &bcq0[j], &w );
          if ( DotProduct3d ( &v[k], &w ) >= 0.0 ) goto next4;
          SubtractPoints3d ( &bcp0[i], &bcq1[j], &w );
          if ( DotProduct3d ( &v[k], &w ) >= 0.0 ) goto next4;
          SubtractPoints3d ( &bcp1[i], &bcq1[j], &w );
          if ( DotProduct3d ( &v[k], &w ) >= 0.0 ) goto next4;
        }
    }
    else {
      for ( i = 0; i <= degree; i++ )
        for ( j = 0; j <= degree; j++ ) {
          SubtractPoints3d ( &bcp0[i], &bcq0[j], &w );
          if ( DotProduct3d ( &v[k], &w ) <= 0.0 ) goto next4;
          SubtractPoints3d ( &bcp1[i], &bcq0[j], &w );
          if ( DotProduct3d ( &v[k], &w ) <= 0.0 ) goto next4;
          SubtractPoints3d ( &bcp0[i], &bcq1[j], &w );
          if ( DotProduct3d ( &v[k], &w ) <= 0.0 ) goto next4;
          SubtractPoints3d ( &bcp1[i], &bcq1[j], &w );
          if ( DotProduct3d ( &v[k], &w ) <= 0.0 ) goto next4;
        }
    }
    return false;
next4:
    ;
  }
  return true;
} /*_rbez_ConvexHullTestd*/

static boolean _rbez_UniquenessTestd ( int degree,
                                       point3d *bcp0, point3d *bcp1,
                                       point3d *bcq0, point3d *bcq1,
                                       vector3d *f, vector3d *fu,
                                       vector3d *fv, vector3d *ft,
                                       double *K, boolean *error )
{
  void     *sp;
  vector3d p0, p0u, p1, p1u, q0, q0v, q1, q1v, f0, f1, ff;
  vector3d *acp0, *acp1, *acq0, *acq1;
  double   a[9], KK, SK;
  int      i, j, P[3], Q[3];

  sp = pkv_GetScratchMemTop ();
  *error = false;
  acp0 = pkv_GetScratchMem ( 4*(degree+1)*sizeof(vector3d) );
  if ( !acp0 ) {
    *error = true;
    goto test_failed;
  }
  acp1 = &acp0[degree+1];
  acq0 = &acp1[degree+1];
  acq1 = &acq0[degree+1];
        /* find derivatives of f at the cube centre */
  mbs_BCHornerDerC3d ( degree, bcp0, 0.5, &p0, &p0u );
  mbs_BCHornerDerC3d ( degree, bcp1, 0.5, &p1, &p1u );
  mbs_BCHornerDerC3d ( degree, bcq0, 0.5, &q0, &q0v );
  mbs_BCHornerDerC3d ( degree, bcq1, 0.5, &q1, &q1v );
  SubtractPoints3d ( &p0, &q0, &f0 );
  SubtractPoints3d ( &p1, &q1, &f1 );
  MidPoint3d ( &f0, &f1, f );
  MidPoint3d ( &p0u, &p1u, fu );
  SetVector3d ( fv, -0.5*(q0v.x+q1v.x), -0.5*(q0v.y+q1v.y), -0.5*(q0v.z+q1v.z) );
  SubtractPoints3d ( &p1, &p0, ft );
  AddVector3d ( ft, &q0, ft );
  SubtractPoints3d ( ft, &q1, ft );
        /* compose the mapping */
  a[0] = fu->x;  a[1] = fv->x;  a[2] = ft->x;
  a[3] = fu->y;  a[4] = fv->y;  a[5] = ft->y;
  a[6] = fu->z;  a[7] = fv->z;  a[8] = ft->z;
  if ( !pkn_GaussDecomposePLUQd ( 3, a, P, Q ) )
    goto test_failed;
        /* transform the control points */
  memcpy ( acp0, bcp0, (degree+1)*sizeof(vector3d) );
  memcpy ( acp1, bcp1, (degree+1)*sizeof(vector3d) );
  memcpy ( acq0, bcq0, (degree+1)*sizeof(vector3d) );
  memcpy ( acq1, bcq1, (degree+1)*sizeof(vector3d) );
  for ( i = 0; i <= degree; i++ ) {
    pkn_multiSolvePLUQd ( 3, a, P, Q, 1, 1, &acp0[i].x );
    pkn_multiSolvePLUQd ( 3, a, P, Q, 1, 1, &acp1[i].x );
    pkn_multiSolvePLUQd ( 3, a, P, Q, 1, 1, &acq0[i].x );
    pkn_multiSolvePLUQd ( 3, a, P, Q, 1, 1, &acq1[i].x );
  }
        /* find the constants K_0, K_1, K_2 */
  memset ( K, 0, 3*sizeof(double) );
  for ( i = 0; i < degree; i++ ) {
    SubtractPoints3d ( &acp0[i+1], &acp0[i], &p0u );
    MultVector3d ( (double)degree, &p0u, &p0u );
    p0u.x -= 1.0;
    KK = DotProduct3d ( &p0u, &p0u );
    if ( KK > K[0] ) {
      if ( KK > 1.0 ) goto test_failed;
      K[0] = KK;
    }
    SubtractPoints3d ( &acp1[i+1], &acp1[i], &p1u );
    MultVector3d ( (double)degree, &p1u, &p1u );
    p1u.x -= 1.0;
    KK = DotProduct3d ( &p1u, &p1u );
    if ( KK > K[0] ) {
      if ( KK > 1.0 ) goto test_failed;
      K[0] = KK;
    }
  }
  SK = K[0];
  for ( i = 0; i < degree; i++ ) {
    SubtractPoints3d ( &acq0[i+1], &acq0[i], &q0v );
    MultVector3d ( (double)degree, &q0v, &q0v );
    q0v.y += 1.0;
    KK = DotProduct3d ( &q0v, &q0v );
    if ( KK > K[1] ) {
      if ( SK+KK > 1.0 ) goto test_failed;
      K[1] = KK;
    }
    SubtractPoints3d ( &acq1[i+1], &acq1[i], &q1v );
    MultVector3d ( (double)degree, &q1v, &q1v );
    q1v.y += 1.0;
    KK = DotProduct3d ( &q1v, &q1v );
    if ( KK > K[1] ) {
      if ( SK+KK > 1.0 ) goto test_failed;
      K[1] = KK;
    }
  }
  SK += K[1];
  for ( i = 0; i <= degree; i++ )
    for ( j = 0; j <= degree; j++ ) {
      SubtractPoints3d ( &acp0[i], &acq0[j], &f0 );
      SubtractPoints3d ( &acp1[i], &acq1[j], &f1 );
      SubtractPoints3d ( &f1, &f0, &ff );
      ff.z -= 1.0;
      KK = DotProduct3d ( &ff, &ff );
      if ( KK > K[2] ) {
        if ( SK+KK > 1.0 ) goto test_failed;
        K[2] = KK;
      }
    }
  pkv_SetScratchMemTop ( sp );
  return true;

test_failed:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_rbez_UniquenessTestd*/

static boolean _rbez_NewtonMethodd ( int degree,
                                     point3d *bcp0, point3d *bcp1,
                                     point3d *bcq0, point3d *bcq1,
                                     vector3d *f, vector3d *fu,
                                     vector3d *fv, vector3d *ft, double *z )
{
#define MAXITER 7
#define EPS     1.0e-8
#define DELTA   1.0e-8
  double s, a[9];
  int    i, j;
  vector3d p0, p0u, p1, p1u, q0, q0v, q1, q1v, f0, f1;

  z[0] = z[1] = z[2] = 0.5;
  for ( i = 0; i < MAXITER; i++ ) {
    if ( i ) {
      mbs_BCHornerDerC3d ( degree, bcp0, z[0], &p0, &p0u );
      mbs_BCHornerDerC3d ( degree, bcp1, z[0], &p1, &p1u );
      mbs_BCHornerDerC3d ( degree, bcq0, z[1], &q0, &q0v );
      mbs_BCHornerDerC3d ( degree, bcq1, z[1], &q1, &q1v );
      SubtractPoints3d ( &p0, &q0, &f0 );
      SubtractPoints3d ( &p1, &q1, &f1 );
      InterPoint3d ( &f0, &f1, z[2], f );
      s = DotProduct3d ( f, f );
      if ( s < EPS*EPS )
        return true;
      InterPoint3d ( &p0u, &p1u, z[2], fu );
      InterPoint3d ( &q0v, &q1v, z[2], fv );
      SetVector3d ( fv, -fv->x, -fv->y, -fv->z );
      SubtractPoints3d ( &p1, &q1, ft );
      SubtractPoints3d ( ft, &p0, ft );
      AddVector3d ( ft, &q0, ft );
    }
    else {
      s = DotProduct3d ( f, f );
      if ( s < EPS*EPS )
        return true;
    }
    a[0] = fu->x;  a[1] = fv->x;  a[2] = ft->x;
    a[3] = fu->y;  a[4] = fv->y;  a[5] = ft->y;
    a[6] = fu->z;  a[7] = fv->z;  a[8] = ft->z;
    if ( !pkn_multiGaussSolveLinEqd ( 3, a, 1, 1, &f->x ) )
      return false;
    z[0] -= f->x;  z[1] -= f->y;  z[2] -= f->z;
    for ( j = 0; j < 3; j++ )
      if ( z[j] < -0.25 || z[j] > 1.25 )
        return false;
    s = DotProduct3d ( f, f );
    if ( s < DELTA*DELTA )
      return true;
  }
  return false;
#undef DELTA
#undef EPS
#undef MAXITER
} /*_rbez_NewtonMethodd*/

static boolean _rbez_SolutionOKd ( double *z, stack_el *el, double *tfh )
{
  int i;

  for ( i = 0; i < 3; i++ )
    if ( z[i] < 0.0 || z[i] > 1.0 )
      return false;
  tfh[0] = (1.0-z[0])*el->u0 + z[0]*el->u1;
  tfh[1] = (1.0-z[1])*el->v0 + z[1]*el->v1;
  tfh[2] = (1.0-z[2])*el->t0 + z[2]*el->t1;
  return true;
}/*_rbez_SolutionOKd*/

static boolean _rbez_SecondTestd ( int degree, double *K, double *z )
{
  double r[3], RR, KK;
  int    i;

  for ( i = 0; i < 3; i++ )
    if ( z[i] < 0.0 )      r[i] = 1.0-2.0*z[i];
    else if ( z[i] > 1.0 ) r[i] = 2.0*z[i] - 1.0;
    else                   r[i] = 1.0;
  RR = pkv_rpower ( r[0]*r[1], degree-1 );
  for ( i = 0; i < 3; i++ )
    r[i] *= r[i];
  KK = RR*RR*(K[0]*r[1]*r[2] + K[1]*r[0]*r[2] + K[2]*r[0]*r[1]);
  return KK >= 1.0;
} /*_rbez_SecondTestd*/

static void _rbez_Subdivided ( stack_el *el, int degree,
                               point3d *bcp0, point3d *bcp1,
                               point3d *bcq0, point3d *bcq1,
                               point3d *acp0, point3d *acp1,
                               point3d *acq0, point3d *acq1 )
{
  vector3d v, v0, v1;
  double   lgt[3], vn;
  int      i, j, size;
  char     divdir;

  if ( el->u0 >= el->v0 ) {  /* a special case - just a single arc */
    divdir = ( el->u1-el->u0 >= el->v1-el->v0 ) ? 0 : 1;
    goto subdiv;
  }
        /* determine the division direction */
  memset ( lgt, 0, 3*sizeof(double) );
  for ( i = 0; i < degree; i++ ) {
    SubtractPoints3d ( &bcp0[i+1], &bcp0[i], &v );
    vn = DotProduct3d ( &v, &v );
    if ( vn > lgt[0] ) lgt[0] = vn;
    SubtractPoints3d ( &bcp1[i+1], &bcp1[i], &v );
    vn = DotProduct3d ( &v, &v );
    if ( vn > lgt[0] ) lgt[0] = vn;
  }
  for ( j = 0; j < degree; j++ ) {
    SubtractPoints3d ( &bcq0[j+1], &bcq0[j], &v );
    vn = DotProduct3d ( &v, &v );
    if ( vn > lgt[1] ) lgt[1] = vn;
    SubtractPoints3d ( &bcq1[j+1], &bcq1[j], &v );
    vn = DotProduct3d ( &v, &v );
    if ( vn > lgt[1] ) lgt[1] = vn;
  }
  for ( i = 0; i <= degree; i++ )
    for ( j = 0; j <= degree; j++ ) {
      SubtractPoints3d ( &bcp0[i], &bcq0[j], &v0 );
      SubtractPoints3d ( &bcp1[i], &bcq1[j], &v1 );
      SubtractPoints3d ( &v1, &v0, &v );
      vn = DotProduct3d ( &v, &v );
      if ( vn > lgt[2] ) lgt[2] = vn;
    }

  divdir = 0;
  if ( lgt[1] > lgt[0] )
    { divdir = 1;  lgt[0] = lgt[1]; }
  if ( lgt[2] > ((double)(degree*degree))*lgt[0] )
    divdir = 2;
        /* subdivide */
subdiv:
  size = (degree+1)*sizeof(point3d);
  el[0].level ++;
  el[1] = el[0];
  switch ( divdir ) {
case 0:
    el[0].u0 = el[1].u1 = 0.5*(el[0].u0+el[0].u1);
    mbs_BisectBC3d ( degree, bcp0, acp0 );
    mbs_BisectBC3d ( degree, bcp1, acp1 );
    memcpy ( acq0, bcq0, size );
    memcpy ( acq1, bcq1, size );
    break;
case 1:
    el[0].v0 = el[1].v1 = 0.5*(el[0].v0+el[0].v1);
    memcpy ( acp0, bcp0, size );
    memcpy ( acp1, bcp1, size );
    mbs_BisectBC3d ( degree, bcq0, acq0 );
    mbs_BisectBC3d ( degree, bcq1, acq1 );
    break;
case 2:
    el[0].t0 = el[1].t1 = 0.5*(el[0].t0+el[0].t1);
    memcpy ( acp0, bcp0, size );
    memcpy ( acq0, bcq0, size );
    for ( i = 0; i <= degree; i++ ) {
      MidPoint3d ( &bcp0[i], &bcp1[i], &acp1[i] );
      MidPoint3d ( &bcq0[i], &bcq1[i], &acq1[i] );
    }
    memcpy ( bcp0, acp1, size );
    memcpy ( bcq0, acq1, size );
    break;
  }
} /*_rbez_Subdivided*/

static boolean _rbez_HomotopyNeighboursDisjointd ( int degree,
                                                   point3d *bcp0, point3d *bcp1,
                                                   point3d *bcq0, point3d *bcq1,
                                                   boolean same )
{
  point3d  p0, q1;
  vector3d v;
  double   a, b;
  int      i;

  MidPoint3d ( &bcp0[0], &bcp1[0], &p0 );
  MidPoint3d ( &bcq0[degree], &bcq1[degree], &q1 );
  SubtractPoints3d ( &q1, &p0, &v );
        /* it is assumed that bcp?[degree] == bcq?[0]            */
        /* the procedure returns true if both polylines made of  */
        /* bcp?[0],...,bcp?[degree],bcq?[1],...,bcq?[degree] are */
        /* monotonic with respect to the line of direction of v  */
  a = DotProduct3d ( &v, &bcp0[0] );
  for ( i = 1; i <= degree; i++ ) {
    b = DotProduct3d ( &v, &bcp0[i] );
    if ( b <= a )
      return false;
    a = b;
  }
  if ( !same ) {
    for ( i = 1; i <= degree; i++ ) {
      b = DotProduct3d ( &v, &bcq0[i] );
      if ( b <= a )
        return false;
      a = b;
    }
  }
  a = DotProduct3d ( &v, &bcp1[0] );
  for ( i = 1; i <= degree; i++ ) {
    b = DotProduct3d ( &v, &bcp1[i] );
    if ( b <= a )
      return false;
    a = b;
  }
  if ( !same ) {
    for ( i = 1; i <= degree; i++ ) {
      b = DotProduct3d ( &v, &bcq1[i] );
      if ( b <= a )
        return false;
      a = b;
    }
  }
  return true;
} /*_rbez_HomotopyNeighboursDisjointd*/

static boolean _rbez_HomotopyBCDisjointd ( int degree,
                        double u0, double u1, point3d *bcp0, point3d *bcp1,
                        double v0, double v1, point3d *bcq0, point3d *bcq1,
                        double *tfh, boolean *error )
{
  void     *sp;
  stack_el *st;
  point3d  *pst, *acp0, *acp1, *acq0, *acq1;
  int      stp, size;
  vector3d p, pu, pv, pt;
  double   K[3], z[3];

  sp = pkv_GetScratchMemTop ();
  *error = false;
        /* allocate the stack */
  st = pkv_GetScratchMem ( (MAXLEVEL+1)*sizeof(stack_el) );
  pst = pkv_GetScratchMem ( (MAXLEVEL+1)*4*(degree+1)*sizeof(point3d) );
  if ( !st || !pst ) {
    *error = true;
    goto test_failed;
  }
        /* put the mapping on the stack */
  st[0].u0 = u0;  st[0].v0 = v0;  st[0].t0 = 0.0;
  st[0].u1 = u1;  st[0].v1 = v1;  st[0].t1 = 1.0;
  st[0].level = 0;
  size = (degree+1)*sizeof(point3d);
  memcpy ( pst, bcp0, size );
  memcpy ( &pst[degree+1], bcp1, size );
  memcpy ( &pst[2*(degree+1)], bcq0, size );
  memcpy ( &pst[3*(degree+1)], bcq1, size );
  stp = 1;
        /* now try to find a zero of the mapping */
  do {
    stp --;
    if ( st[stp].u0 >= st[stp].v1 )
      goto skip;
    acp0 = &pst[stp*4*(degree+1)];
    acp1 = &acp0[degree+1];
    acq0 = &acp1[degree+1];
    acq1 = &acq0[degree+1];
    if ( st[stp].u1 >= st[stp].v0 ) {
      if ( _rbez_HomotopyNeighboursDisjointd ( degree, acp0, acp1, acq0, acq1,
                                              st[stp].u0 >= st[stp].v0 ) )
        goto skip;
    }
    if ( _rbez_ConvexHullTestd ( degree, acp0, acp1, acq0, acq1 ) ) {
      if ( st[stp].u0 >= st[stp].v0 && st[stp].u1 >= st[stp].v1 )
        goto subdivide;
      if ( _rbez_UniquenessTestd ( degree, acp0, acp1, acq0, acq1,
                                  &p, &pu, &pv, &pt, K, error ) ) {
        if ( _rbez_NewtonMethodd ( degree, acp0, acp1, acq0, acq1,
                                  &p, &pu, &pv, &pt, z ) ) {
          if ( _rbez_SolutionOKd ( z, &st[stp], tfh ) )
            goto test_failed;
          else if ( _rbez_SecondTestd ( degree, K, z ) )
            goto subdivide;
        }
        else goto subdivide;
      }
      else {
subdivide:
        if ( st[stp].level < MAXLEVEL ) {
          _rbez_Subdivided ( &st[stp], degree, acp0, acp1, acq0, acq1,
                            &acp0[4*(degree+1)], &acp1[4*(degree+1)],
                            &acq0[4*(degree+1)], &acq1[4*(degree+1)] );
          stp += 2;
        }
        else {  /* a singularity */
          *tfh = 0.5*(st[stp].t0+st[stp].t1);
          goto test_failed;
        }
      }
    }
skip:
    ;
  } while ( stp > 0 );

  pkv_SetScratchMemTop ( sp );
  return true;

test_failed:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_rbez_HomotopyBCDisjointd*/

boolean rbez_HomotopicClosedBSC3d ( int degree, int lastknot, double *knots,
                                    point3d *cpoints0, point3d *cpoints1,
                                    double *tfh, boolean *error )
{
  void    *sp;
  int     ku, lkn;
  double  *kn, u0, u1, v0, v1;
  point3d *bcp0, *bcp1;
  int     i, j;

  sp = pkv_GetScratchMemTop ();
  *error = false;
        /* convert the curves to Bezier arcs */
  ku = mbs_NumKnotIntervalsd ( degree, lastknot, knots );
  lkn = mbs_LastknotMaxInsd ( degree, lastknot, knots, &ku );
  bcp0 = pkv_GetScratchMem ( 2*ku*(degree+1)*sizeof(point3d) );
  kn = pkv_GetScratchMemd ( lkn+1 );
  if ( !bcp0 || !kn ) {
    *error = true;
    goto test_failed;
  }
  bcp1 = &bcp0[ku*(degree+1)];
  mbs_BSToBezC3d ( degree, lastknot, knots, cpoints0, &ku, &lkn, kn, bcp0 );
  mbs_BSToBezC3d ( degree, lastknot, knots, cpoints1, &ku, &lkn, kn, bcp1 );
        /* test the homotopy violation */
  for ( i = 0; i < ku; i++ ) {
    u0 = kn[i*(degree+1)+degree];
    u1 = kn[(i+1)*(degree+1)];
    for ( j = i; j < ku; j++ ) {
      if ( i == 0 && j == ku-1 ) {
        v0 = u0 - (kn[(j+1)*(degree+1)]-kn[j*(degree+1)+degree]);
        v1 = u0;
        if ( !_rbez_HomotopyBCDisjointd ( degree,
                      v0, v1, &bcp0[j*(degree+1)], &bcp1[j*(degree+1)],
                      u0, u1, &bcp0[i*(degree+1)], &bcp1[i*(degree+1)],
                      tfh, error ) ) {
          tfh[1] += knots[lastknot-degree]-knots[degree];
          goto test_failed;
        }
      }
      else {
        v0 = kn[j*(degree+1)+degree];
        v1 = kn[(j+1)*(degree+1)];
        if ( !_rbez_HomotopyBCDisjointd ( degree,
                      u0, u1, &bcp0[i*(degree+1)], &bcp1[i*(degree+1)],
                      v0, v1, &bcp0[j*(degree+1)], &bcp1[j*(degree+1)],
                      tfh, error ) )
          goto test_failed;
      }
    }
  }
  pkv_SetScratchMemTop ( sp );
  return true;

test_failed:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*rbez_HomotopicClosedBSC3d*/

