
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2013                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#define CONST_

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "raybez.h"

#include "raybezprivatef.h"

/* ////////////////////////////////////////////////////////////////////////// */
static boolean SubdividePatch ( int n, int m, point2f *p, point2f *q )
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
} /*SubdividePatch*/

static boolean EnterSolution ( int degree, point4f *cp, int object_id, ray3f *ray,
                               float t0, float t1, float u0, float u1,
                               point2f *z,
                               RayObjectIntersf *inters )
{
  point3f p;
  float   t, u;

  if ( z->x >= 0.0 && z->x <= 1.0 && z->y >= 0.0 && z->y <= 1.0 ) {
    t = t0 + (t1-t0)*z->x;
    u = fabs ( u0 + (u1-u0)*z->y );
    mbs_BCHornerC3Rf ( degree, cp, t, &p );
    AddVector3Mf ( &ray->p, &ray->v, u, &inters->p );
    SubtractPoints3f ( &inters->p, &p, &inters->nv );
    NormalizeVector3f ( &inters->nv );
    inters->object_id = object_id;
    inters->u = t;
    inters->v = 0.0;
    inters->t = u;
    return true;
  }
  else
    return false;
} /*EnterSolution*/

int rbez_FindRayRBezcOffsetIntersf ( RBezCurveTreefp tree, ray3f *ray,
                                     int maxlevel, int maxinters,
                                     int *ninters, RayObjectIntersf *inters )
{
  typedef  struct {
      float    u0, u1, t0, t1;
      vector2f *mcp;
      int      level;
    } stack_el;

  void     *sp;
  int      _ninters;
  point4f  *cp, *tcp;
  vector4f *dtcp, *stcp, *sdtcp;
  int      i, j, b;
  vector3f rpn, *pdp;
  boolean  plus_e3;
  float    s, gamma;
  float    s0, s1, s01, s0s0, s0s1, s1s1, rad, rad2;
  float    *pp, *ww, *zw, *pdpp, *zdzw;
  int      degree, degreem1, deg2, deg2m1, deg3, mapcpsize, ncp;
  vector2f *mapcp, *mcp, p, pu, pv, z;
  float    K1, K2;
  int      stp;
  stack_el *stack;
  float    u0, u1, t0, t1;

  sp = pkv_GetScratchMemTop ();
  _ninters = 0;
  if ( !rbez_TestRayBBoxf ( ray, &tree->root->bbox ) )
    goto way_out;
        /* allocate workspace */
  degree = tree->degree;
  deg2 = 2*degree;
  deg2m1 = deg2-1;
  deg3 = 3*degree-1;
  degreem1 = degree-1;
  cp = tree->root->ctlpoints;
  mapcpsize = 3*(deg3+1);
  tcp = pkv_GetScratchMem ( (4*degree+2)*sizeof(point4f) );
  if ( !tcp )
    goto failure;
  dtcp = &tcp[degree+1];
  stcp = &dtcp[degree];
  sdtcp = &stcp[degree+1];
  pp = pkv_GetScratchMemf ( 3*(deg3+1) );
  if ( !pp )
    goto failure;
  ww = &pp[deg3+1];
  zw = &ww[deg3+1];
  pdp = pkv_GetScratchMem ( deg2*sizeof(vector3f) );
  pdpp = pkv_GetScratchMemf ( 2*(deg3+1) );
  if ( !pdp || !pdpp )
    goto failure;
  zdzw = &pdpp[deg3+1];

  mapcp = pkv_GetScratchMem ( (maxlevel+1)*mapcpsize*sizeof(vector2f) );
  stack = pkv_GetScratchMem ( (maxlevel+1)*sizeof(stack_el) );
  if ( !mapcp || !stack )
    goto failure;
  rad = tree->ext;
  rad2 = rad*rad;
        /* construct the Householder reflection */
  rpn = ray->v;
  if ( rpn.z >= 0.0 )
    { rpn.z += 1.0;  plus_e3 = false; }
  else
    { rpn.z -= 1.0;  plus_e3 = true; }
  gamma = -2.0/DotProduct3f ( &rpn, &rpn );
        /* transform the control points */
  for ( i = 0; i <= degree; i++ ) {
    AddVector3Mf ( (vector3f*)&cp[i], &ray->p, -cp[i].w, (vector3f*)&tcp[i] );
    tcp[i].w = cp[i].w;
    s = gamma*DotProduct3f ( (vector3f*)&tcp[i], &rpn );
    AddVector3Mf ( (vector3f*)&tcp[i], &rpn, s, (vector3f*)&tcp[i] );
  }
        /* compute the control points of the derivative */
  for ( i = 0; i < degree; i++ ) {
    SubtractPoints4f ( &tcp[i+1], &tcp[i], &dtcp[i] );
    MultVector4f ( (float)degree, &dtcp[i], &dtcp[i] );
  }
        /* find the range for the ray parameter */
  s0 = s1 = tcp[0].z/tcp[0].w;
  for ( i = 1; i <= degree; i++ ) {
    s = tcp[i].z/tcp[i].w;
    if ( s < s0 )
      s0 = s;
    else if ( s > s1 )
      s1 = s;
  }
  s0 -= rad;
  s1 += rad;
  if ( plus_e3 ) {
    if ( s1 <= 0.0 ) goto way_out;
    if ( s0 < 0.0 ) s0 = 0.0;
  }
  else {
    if ( s1 >= 0.0 ) goto way_out;
    if ( s0 > 0.0 ) s0 = 0.0;
  }
        /* construct the Bezier representation of the mapping, */
        /* whose zeros are to be found */
  memcpy ( stcp, tcp, (degree+1)*sizeof(point4f) );
  b = degree;
  for ( i = 1; i < degree; i++ ) {
    stcp[i].x *= (float)b;
    stcp[i].y *= (float)b;
    stcp[i].z *= (float)b;
    stcp[i].w *= (float)b;
    b = (b*(degree-i)) / (i+1);
  }
  memcpy ( sdtcp, dtcp, degree*sizeof(point4f) );
  b = degreem1;
  for ( i = 1; i < degreem1; i++ ) {
    sdtcp[i].x *= (float)b;
    sdtcp[i].y *= (float)b;
    sdtcp[i].z *= (float)b;
    sdtcp[i].w *= (float)b;
    b = (b*(degreem1-i)) / (i+1);
  }
  memset ( pp, 0, 3*(deg3+1)*sizeof(float) );
  for ( i = 0; i <= degree; i++ )
    for ( j = 0; j <= degree; j++ ) {
      pp[i+j] += DotProduct3f ( (vector3f*)&stcp[i], (vector3f*)&stcp[j] );
      ww[i+j] += stcp[i].w*stcp[j].w;
      zw[i+j] += stcp[i].z*stcp[j].w;
    }
  memset ( pdp, 0, deg2*sizeof(vector3f) );
  for ( i = 0; i <= degree; i++ )
    for ( j = 0; j <= degreem1; j++ ) {
      pdp[i+j].x += stcp[i].x*sdtcp[j].w-stcp[i].w*sdtcp[j].x;
      pdp[i+j].y += stcp[i].y*sdtcp[j].w-stcp[i].w*sdtcp[j].y;
      pdp[i+j].z += stcp[i].z*sdtcp[j].w-stcp[i].w*sdtcp[j].z;
    }
  memset ( pdpp, 0, 2*(deg3+1)*sizeof(float) );
  for ( i = 0; i <= deg2m1; i++ )
    for ( j = 0; j <= degree; j++ ) {
      pdpp[i+j] += DotProduct3f ( &pdp[i], (vector3f*)&stcp[j] );
      zdzw[i+j] += pdp[i].z*stcp[j].w;
    }
  for ( i = deg2; i < deg3; i++ )
    for ( j = i; j >= 0; j-- ) {
      pp[j+1] += pp[j];
      ww[j+1] += ww[j];
      zw[j+1] += zw[j];
    }
  b = deg3;
  for ( i = 1; i < deg3; i++ ) {
    pp[i] /= (float)b;
    ww[i] /= (float)b;
    zw[i] /= (float)b;
    pdpp[i] /= (float)b;
    zdzw[i] /= (float)b;
    b = (b*(deg3-i)) / (i+1);
  }
          /* actual mapping construction */
  s0s0 = s0*s0;
  s0s1 = s0*s1;
  s1s1 = s1*s1;
  s01 = 0.5*(s0+s1);
  for ( i = j = 0;  i <= deg3;  i++, j += 3 ) {
    mapcp[j].x = (s0s0-rad2)*ww[i] - 2.0*s0*zw[i] + pp[i];
    mapcp[j].y = s0*zdzw[i] - pdpp[i];
    mapcp[j+1].x = (s0s1-rad2)*ww[i] - 2.0*s01*zw[i] + pp[i];
    mapcp[j+1].y = s01*zdzw[i] - pdpp[i];
    mapcp[j+2].x = (s1s1-rad2)*ww[i] - 2.0*s1*zw[i] + pp[i];
    mapcp[j+2].y = s1*zdzw[i] - pdpp[i];
  }
  ncp = 3*(deg3+1);
        /* solving the equations using recursive subdivision */
  stack[0].u0 = s0;      stack[0].u1 = s1;
  stack[0].t0 = 0.0;     stack[0].t1 = 1.0;
  stack[0].mcp = mapcp;  stack[0].level = 0;
  stp = 1;
  do {
    stp --;
    u0 = stack[stp].u0;  u1 = stack[stp].u1;
    t0 = stack[stp].t0;  t1 = stack[stp].t1;
    mcp = stack[stp].mcp;
    if ( _rbez_ConvexHullTest2f ( ncp, mcp ) ) {
      if ( _rbez_UniquenessTest2f ( deg3, 2, ncp, mcp, &p, &pu, &pv, &K1, &K2 ) ) {
        switch ( _rbez_NewtonMethod2f ( deg3, 2, mcp, &p, &pu, &pv, &z ) ) {
      case RBEZ_NEWTON_YES:
          /* regular solution found */
          if ( EnterSolution ( degree, cp, tree->object_id, ray,
                               t0, t1, u0, u1, &z, &inters[_ninters] ) )
            _ninters ++;
          else if ( _rbez_SecondTest2f ( &z, deg3, 2, K1, K2 ) )
            goto subdivide;
          break;
      case RBEZ_NEWTON_NO:
          goto subdivide;
      case RBEZ_NEWTON_ERROR:
          goto failure;
        }
      }
      else {
subdivide:
        if ( stack[stp].level < maxlevel ) {
          stack[stp+1].mcp = &stack[stp].mcp[mapcpsize];
          stack[stp+1].level = stack[stp].level = stack[stp].level+1;
          stack[stp+1].u0 = u0;
          stack[stp+1].t0 = t0;
          if ( SubdividePatch ( deg3, 2, stack[stp].mcp, stack[stp+1].mcp ) ) {
            stack[stp+1].u1 = u1;
            stack[stp].t0 = stack[stp+1].t1 = 0.5*(t0+t1);
          }
          else {
            stack[stp].u0 = stack[stp+1].u1 = 0.5*(u0+u1);
            stack[stp+1].t1 = t1;
          }
          stp += 2;
        }
        else {
          /* solution at singularity found */
          z.x = z.y = 0.5;
          if ( EnterSolution ( degree, cp, tree->object_id, ray,
                               t0, t1, u0, u1, &z, &inters[_ninters] ) )
            _ninters ++;
        }
      }
    }
  } while ( stp > 0 && _ninters < maxinters );

way_out:
  pkv_SetScratchMemTop ( sp );
  *ninters = _ninters;
  return _ninters;

failure:
  pkv_SetScratchMemTop ( sp );
  return -1;
} /*rbez_FindRayRBezcOffsetIntersf*/

