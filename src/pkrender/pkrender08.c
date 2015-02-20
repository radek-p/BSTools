
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2015                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "pkvaria.h"
#include "pkvthreads.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "camera.h"
#include "raybez.h"
#include "pkrender.h"

#include "pkrenderprivate.h"

static void _rend_FindRayTriangleIntersd ( triangle_data *trd, ray3d *ray,
                                           int *nint, RayObjectIntersd *inters )
{
#define TOL 1.0e-10
  double   a, b, t;
  vector3d pq;
  int      i;

  SubtractPoints3d ( &trd->p0, &ray->p, &pq );
  a = DotProduct3d ( &ray->v, &trd->n );
  if ( fabs(a) < TOL )
    return;
  b = DotProduct3d ( &pq, &trd->n );
  t = b/a;
  if ( t <= 0.0 )
    return;
  AddVector3Md ( &pq, &ray->v, -t, &pq );
  a = DotProduct3d ( &trd->a1, &pq );
  if ( a >= 0.0 )
    return;
  b = DotProduct3d ( &trd->a2, &pq );
  if ( b >= 0.0 || (a+b) <= -1.0 )
    return;
  i = *nint;
  inters[i].object_id = 0;
  AddVector3Md ( &ray->p, &ray->v, t, &inters[i].p );
  inters[i].nv = trd->n;
  inters[i].u = a;  inters[i].v = b;  inters[i].t = t;
  *nint = i+1;
#undef TOL
} /*_rend_FindRayTriangleIntersd*/

static boolean rFindRayTrInters ( renderobj *obj_tab, rendertnode *node, ray3d *ray,
                                  RayObjectIntersd *inters, int *k )
{
  void *sp;
  int  _k, nint;

  sp = pkv_GetScratchMemTop ();
  *k = -1;
  nint = 0;
  _k = node->k;
  if ( _k >= 0 ) {
    RayObjectIntersd *_inters, *_intp;
    int              i;
           /* it is a leaf, i.e. a Bezier patch or a Bezier curve; no ray/box */
           /* intersection is tested here, as the ray/patch or ray/curve offset */
           /* intersection procedures have their own tests */

    _inters = pkv_GetScratchMem ( MAXINTERS*sizeof(RayObjectIntersd) );
    if ( !_inters ) {
      pkv_SetScratchMemTop ( sp );
      return false;
    }
    switch ( obj_tab[_k].type ) {
  case obj_TRIANGLE:
      _rend_FindRayTriangleIntersd ( obj_tab[_k].triang.trdata, ray,
                                     &nint, _inters );
      break;
  case obj_BSPATCH:
      rbez_FindRayBezPatchIntersd ( obj_tab[_k].bsp.ptree, ray,
                    MAXPLEVEL, MAXINTERS, &nint, _inters );
      break;
  case obj_RBSPATCH:
      rbez_FindRayRBezPatchIntersd ( obj_tab[_k].rbsp.ptree, ray,
                    MAXPLEVEL, MAXINTERS, &nint, _inters );
      break;
  case obj_BEZCURVE:
      rbez_FindRayBezcOffsetIntersd ( obj_tab[_k].bezc.ctree, ray,
                    MAXCLEVEL, MAXINTERS, &nint, _inters );
      break;
  case obj_RBEZCURVE:
      rbez_FindRayRBezcOffsetIntersd ( obj_tab[_k].rbezc.ctree, ray,
                    MAXCLEVEL, MAXINTERS, &nint, _inters );
      break;
  default:
      pkv_SetScratchMemTop ( sp );
      return false;
    }
    if ( nint > 0 ) {
           /* ignore all intersections too close to the ray origin */
      for ( i = 0; i < nint; )
        if ( _inters[i].t < MIN_T )
          _inters[i] = _inters[--nint];
        else
          i ++;
      if ( nint ) {
           /* something was left */
        _intp = _inters;
        for ( i = 1; i < nint; i++ )
          if ( _inters[i].t < _intp->t )
            _intp = &_inters[i];
        *inters = *_intp;
        *k = _k;
      }
    }
  }
  else if ( rbez_TestRayBBoxd ( ray, &node->bbox ) ) {
    int              k0, k1;
    boolean          r0, r1;
    RayObjectIntersd inters0, inters1;

    r0 = rFindRayTrInters ( obj_tab, node->left, ray, &inters0, &k0 );
    r1 = rFindRayTrInters ( obj_tab, node->right, ray, &inters1, &k1 );
    if ( r0 ) {
      if ( r1 ) {
        if ( inters0.t < inters1.t )
          { *inters = inters0;  *k = k0; }
        else
          { *inters = inters1;  *k = k1; }
      }
      else
        { *inters = inters0;  *k = k0; }
      nint = 1;
    }
    else {
      if ( r1 )
        { *inters = inters1;  *k = k1;  nint = 1; }
      else nint = 0;
    }
  }
  pkv_SetScratchMemTop ( sp );
  return nint > 0;
} /*rFindRayTrInters*/

boolean FindRayTrInters ( pkRenderer *rend,
                          ray3d *ray, RayObjectIntersd *inters, int *k )
{
  if ( !rend->root )
    return false;
  return rFindRayTrInters ( rend->obj_tab, rend->root, ray, inters, k );
} /*FindRayTrInters*/

static boolean rFindShadowRayInters ( renderobj *obj_tab, ray3d *ray,
                                      rendertnode *node )
{
  int              _k, nint;
  RayObjectIntersd inters;

  _k = node->k;
  nint = 0;
  if ( _k >= 0 ) {
           /* it is a leaf, i.e. a Bezier patch or a Bezier curve; no ray/box */
           /* intersection is tested here, as the ray/patch or ray/curve offset */
           /* intersection procedures have their own tests */
    switch ( obj_tab[_k].type ) {
  case obj_TRIANGLE:
      _rend_FindRayTriangleIntersd ( obj_tab[_k].triang.trdata, ray,
                                     &nint, &inters );
      break;
  case obj_BSPATCH:
      rbez_FindRayBezPatchIntersd ( obj_tab[_k].bsp.ptree, ray,
                    MAXPLEVEL, 1, &nint, &inters );
      break;
  case obj_RBSPATCH:
      rbez_FindRayRBezPatchIntersd ( obj_tab[_k].rbsp.ptree, ray,
                    MAXPLEVEL, 1, &nint, &inters );
      break;
  case obj_BEZCURVE:
      rbez_FindRayBezcOffsetIntersd ( obj_tab[_k].bezc.ctree, ray,
                    MAXCLEVEL, 1, &nint, &inters );
      break;
  case obj_RBEZCURVE:
      rbez_FindRayRBezcOffsetIntersd ( obj_tab[_k].rbezc.ctree, ray,
                    MAXCLEVEL, 1, &nint, &inters );
      break;
  default:
      return false;
    }
    return nint > 0;
  }
  else if ( rbez_TestRayBBoxd ( ray, &node->bbox ) ) {
    if ( rFindShadowRayInters ( obj_tab, ray, node->left ) )
      return true;
    return rFindShadowRayInters ( obj_tab, ray, node->right );
  }
  else
    return false;
} /*rFindShadowRayInters*/

boolean IsInShadow ( pkRenderer *rend, point3d *p, vector3d *light )
{
  ray3d ray;

  if ( !rend->root )
    return false;
  ray.p = *p;
  ray.v = *light;
  AddVector3Md ( p, light, MIN_T, &ray.p );
  return rFindShadowRayInters ( rend->obj_tab, &ray, rend->root );
} /*IsInShadow*/

