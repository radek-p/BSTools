
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2013, 2015                            */
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
static boolean EnterSolution ( int degree, point3f *cp, int object_id, ray3f *ray,
                               float t0, float t1, float u0, float u1,
                               point2f *z, BezCurveTreeVertexfp vertex,
                               RayObjectIntersf *inters )
{
  point3f p;
  double  t, u;

  if ( z->x >= 0.0 && z->x <= 1.0 && z->y >= 0.0 && z->y <= 1.0 ) {
    t = t0 + (t1-t0)*z->x;
    u = fabs ( u0 + (u1-u0)*z->y );
    if ( !mbs_BCHornerC3f ( degree, cp, t, &p ) )
      return false;
    AddVector3Mf ( &ray->p, &ray->v, u, &inters->p );
    SubtractPoints3f ( &inters->p, &p, &inters->nv );
    NormalizeVector3f ( &inters->nv );
    inters->object_id = object_id;
    inters->u = vertex->t0 +(vertex->t1-vertex->t0)*t;
    inters->v = 0.0;
    inters->t = u;
    inters->extra_info = vertex;
    return true;
  }
  else
    return false;
} /*EnterSolution*/

typedef  struct {
    float    u0, u1, t0, t1;
    vector2f *mcp;
    int      level;
  } stack_el;

static int FindRayBezcOffsetIntersf ( BezCurveTreefp tree,
                 BezCurveTreeVertexfp vertex, ray3f *ray,
                 int maxlevel, int maxinters,
                 int *ninters, RayObjectIntersf *inters,
                 void *workspace )
{
  int      _ninters;
  point3f  *cp, *tcp;
  vector3f *dtcp, *stcp, *sdtcp;
  int      i, j, b;
  vector3f rpn;
  boolean  plus_e3;
  float    s, gamma;
  float    s0, s1, s01, s0s0, s0s1, s1s1, rad, rad2;
  float    *pp, *pdp, *ez, *edz;
  int      degree, deg, degreem1, mapcpsize, ncp;
  vector2f *mapcp, *mcp, p, pu, pv, z;
  float   K1, K2;
  int      stp;
  stack_el *stack;
  float   u0, u1, t0, t1;
  float   *wsp;

  _ninters = *ninters;
        /* allocate workspace */
  degree = tree->degree;
  deg = 2*degree;
  degreem1 = degree-1;
  mapcpsize = 3*(deg+1);
  tcp = workspace;
  dtcp = &tcp[degree+1];
  stcp = &dtcp[degree];
  sdtcp = &stcp[degree+1];
  pp = (float*)&sdtcp[degree];
  pdp = &pp[deg+1];
  ez = &pdp[deg+1];
  edz = &ez[deg+1];
  mapcp = (vector2f*)&edz[deg+1];
  stack = (stack_el*)&mapcp[(maxlevel+1)*mapcpsize];
  wsp = (float*)&stack[maxlevel+1];

  cp = vertex->ctlpoints;
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
    SubtractPoints3f ( &cp[i], &ray->p, &tcp[i] );
    s = gamma*DotProduct3f ( &tcp[i], &rpn );
    AddVector3Mf ( &tcp[i], &rpn, s, &tcp[i] );
  }
        /* compute the control points of the derivative */
  for ( i = 0; i < degree; i++ ) {
    SubtractPoints3f ( &tcp[i+1], &tcp[i], &dtcp[i] );
    MultVector3f ( (float)degree, &dtcp[i], &dtcp[i] );
  }
        /* find the range for the ray parameter */
  s0 = s1 = tcp[0].z;
  for ( i = 1; i <= degree; i++ )
    if ( tcp[i].z < s0 )
      s0 = tcp[i].z;
    else if ( tcp[i].z > s1 )
      s1 = tcp[i].z;
  s0 -= rad;
  s1 += rad;
  if ( plus_e3 ) {
    if ( s1 <= 0.0 ) goto way_out;
    if ( s0 < 0.0 ) s0 = 0.0;
  }
  else {
    if ( s0 >= 0.0 ) goto way_out;
    if ( s1 > 0.0 ) s1 = 0.0;
  }
        /* construct the Bezier representation of the mapping, */
        /* whose zeros are to be found */
  memcpy ( stcp, tcp, (degree+1)*sizeof(point3f) );
  b = degree;
  for ( i = 1; i < degree; i++ ) {
    stcp[i].x *= (float)b;
    stcp[i].y *= (float)b;
    stcp[i].z *= (float)b;
    b = (b*(degree-i)) / (i+1);
  }
  memcpy ( sdtcp, dtcp, degree*sizeof(point3f) );
  b = degreem1;
  for ( i = 1; i < degreem1; i++ ) {
    sdtcp[i].x *= (float)b;
    sdtcp[i].y *= (float)b;
    sdtcp[i].z *= (float)b;
    b = (b*(degreem1-i)) / (i+1);
  }
  memset ( pp, 0, (deg+1)*sizeof(float) );
  memset ( ez, 0, (deg+1)*sizeof(float) );
  for ( i = 0; i <= degree; i++ )
    for ( j = 0; j <= degree; j++ )
      pp[i+j] += DotProduct3f ( &stcp[i], &stcp[j] );
  memset ( pdp, 0, (deg+1)*sizeof(float) );
  for ( i = 0; i <= degree; i++ )
    for ( j = 0; j <= degreem1; j++ )
      pdp[i+j] += DotProduct3f ( &stcp[i], &sdtcp[j] );
  for ( i = deg; i > 0; i-- )
    pdp[i] += pdp[i-1];
  for ( i = 0; i <= degree; i++ )
    ez[i] = stcp[i].z;
  memset ( &ez[degree+1], 0, (deg-degree)*sizeof(float) );
  for ( i = degree; i < deg; i++ )
    for ( j = i; j >= 0; j-- )
      ez[j+1] += ez[j];
  for ( i = 0; i <= degreem1; i++ )
    edz[i] = sdtcp[i].z;
  memset ( &edz[degree], 0, (deg+1-degree)*sizeof(float) );
  for ( i = degreem1; i < deg; i++ )
    for ( j = i; j >= 0; j-- )
      edz[j+1] += edz[j];
  b = deg;
  for ( i = 1; i < deg; i++ ) {
    pp[i]  /= (float)b;
    pdp[i] /= (float)b;
    ez[i]  /= (float)b;
    edz[i] /= (float)b;
    b = (b*(deg-i)) / (i+1);
  }
          /* actual mapping construction */
  s0s0 = s0*s0;
  s0s1 = s0*s1;
  s1s1 = s1*s1;
  s01 = 0.5*(s0+s1);
  for ( i = j = 0;  i <= deg;  i++, j += 3 ) {
    s = pp[i] - rad2;
    mapcp[j].x = s + s0s0 - 2.0*s0*ez[i];
    mapcp[j].y = s0*edz[i] - pdp[i];
    mapcp[j+1].x = s + s0s1 - 2.0*s01*ez[i];
    mapcp[j+1].y = s01*edz[i] - pdp[i];
    mapcp[j+2].x = s + s1s1 - 2.0*s1*ez[i];
    mapcp[j+2].y = s1*edz[i] - pdp[i];
  }
  ncp = 3*(deg+1);
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
      if ( _rbez_UniquenessTest2f ( deg, 2, ncp, mcp,
                                    &p, &pu, &pv, &K1, &K2, wsp ) ) {
        switch ( _rbez_NewtonMethod2f ( deg, 2, mcp,
                                        &p, &pu, &pv, &z, wsp ) ) {
      case RBEZ_NEWTON_YES:
          /* regular solution found */
          if ( EnterSolution ( degree, cp, tree->object_id, ray,
                               t0, t1, u0, u1, &z, vertex, &inters[_ninters] ) )
            _ninters ++;
          else if ( _rbez_SecondTest2f ( &z, deg, 2, K1, K2 ) )
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
          if ( _rbez_SubdividePatch2f ( deg, 2, stack[stp].mcp,
                                        stack[stp+1].mcp ) ) {
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
                               t0, t1, u0, u1, &z, vertex, &inters[_ninters] ) )
            _ninters ++;
        }
      }
    }
  } while ( stp > 0 && _ninters < maxinters );

way_out:
  *ninters = _ninters;
  return _ninters;

failure:
  return -1;
} /*FindRayBezcOffsetIntersf*/

static int r_FindRayBezcOffsetIntersf ( BezCurveTreefp tree,
                 BezCurveTreeVertexfp vertex, ray3f *ray,
                 int maxlevel, int maxinters,
                 int *ninters, RayObjectIntersf *inters,
                 void *workspace )
{
  int result;

  if ( rbez_TestRayBBoxf ( ray, &vertex->bbox ) ) {
    if ( vertex->ctlpoints )  /* offset to a polynomial arc */
      result = FindRayBezcOffsetIntersf ( tree, vertex, ray,
                   maxlevel, maxinters, ninters, inters, workspace );
    else {  /* a spline - need to deal with pieces separately */
      if ( vertex->left )
        result = r_FindRayBezcOffsetIntersf ( tree, vertex->left, ray,
                       maxlevel, maxinters, ninters, inters, workspace );
      else
        result = 0;
      if ( vertex->right && result >= 0 )
        result += r_FindRayBezcOffsetIntersf ( tree, vertex->right, ray,
                       maxlevel, maxinters, ninters, inters, workspace );

    }
  }
  else
    result = 0;
  return result;
} /*r_FindRayBezcOffsetIntersf*/

int rbez_RayBezcOffsetWspSizef ( int degree, int maxlevel )
{
  int mapcpsize;

  mapcpsize = 3*(2*degree+1);
  return (4*degree+2)*sizeof(point3f) + 4*(2*degree+1)*sizeof(float) +
         (maxlevel+1)*mapcpsize*sizeof(vector2f) +
         (maxlevel+1)*sizeof(stack_el) + 24*sizeof(float);
} /*rbez_RayBezcOffsetWspSizef*/

int _rbez_FindRayBezcOffsetIntersf ( BezCurveTreefp tree, ray3f *ray,
                                     int maxlevel, int maxinters,
                                     int *ninters, RayObjectIntersf *inters,
                                     void *workspace )
{
  *ninters = 0;
  return r_FindRayBezcOffsetIntersf ( tree, tree->root, ray, maxlevel,
                                      maxinters, ninters, inters, workspace );
} /*_rbez_FindRayBezcOffsetIntersf*/

int rbez_FindRayBezcOffsetIntersf ( BezCurveTreefp tree, ray3f *ray,
                                    int maxlevel, int maxinters,
                                    int *ninters, RayObjectIntersf *inters )
{
  void *workspace;
  int  result;

  workspace = pkv_GetScratchMem (
                  rbez_RayBezcOffsetWspSizef ( tree->degree, maxlevel ) );
  if ( !workspace )
    return -1;
  *ninters = 0;
  result = r_FindRayBezcOffsetIntersf ( tree, tree->root, ray, maxlevel,
                                        maxinters, ninters, inters, workspace );
  pkv_SetScratchMemTop ( workspace );
  return result;
} /*rbez_FindRayBezcOffsetIntersf*/

