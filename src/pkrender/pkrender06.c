
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

double ShapeFunc ( pkRenderer *rend, int func, ray3d *ray, renderobj *obj,
                   RayObjectIntersd *roint )
{
  int                  n, m;
  point3d              *cp;
  point4d              *rcp;
  vector3d             r, dp;
  BezPatchTreeVertexd  *vertex;
  RBezPatchTreeVertexd *rvertex;
  double               u, v, s;
  double               gaussian, mean;

  if ( func == shapefunc_NONE )
    return 0.0;
  switch ( obj->type ) {
case obj_BSPATCH:
    switch ( func ) {
  case shapefunc_GAUSSIAN:
  case shapefunc_MEAN:
      n = obj->bsp.ptree->n;
      m = obj->bsp.ptree->m;
      vertex = (BezPatchTreeVertexd*)roint->extra_info;
      while ( vertex->up ) {
        if ( !vertex->up->ctlpoints )
          break;
        else
          vertex = vertex->up;
      }
      cp = vertex->ctlpoints;
      u = (roint->u-vertex->u0)/(vertex->u1-vertex->u0);
      v = (roint->v-vertex->v0)/(vertex->v1-vertex->v0);
      mbs_GMCurvaturesBP3d ( n, m, cp, u, v, &gaussian, &mean );
      return func == shapefunc_GAUSSIAN ? gaussian : mean;
  case shapefunc_REFLECTION:
      goto reflection;
  case shapefunc_HIGHLIGHT:
      goto highlight;
  case shapefunc_LAMBERTIAN:
      goto lambertian;
  case shapefunc_CROSSECTIONS:
      goto crossections;
      return DotProduct3d ( &rend->sectiondir, &roint->p );
  default:
      return 0.0;
    }

case obj_RBSPATCH:
    switch ( func ) {
  case shapefunc_GAUSSIAN:
  case shapefunc_MEAN:
      n = obj->rbsp.ptree->n;
      m = obj->rbsp.ptree->m;
      rvertex = (RBezPatchTreeVertexd*)roint->extra_info;
      while ( rvertex->up ) {
        if ( !rvertex->up->ctlpoints )
          break;
        else
          rvertex = rvertex->up;
      }
      rcp = rvertex->ctlpoints;
      u = (roint->u-rvertex->u0)/(rvertex->u1-rvertex->u0);
      v = (roint->v-rvertex->v0)/(rvertex->v1-rvertex->v0);
      mbs_GMCurvaturesBP3Rd ( n, m, rcp, u, v, &gaussian, &mean );
      return func == shapefunc_GAUSSIAN ? gaussian : mean;
  case shapefunc_REFLECTION:
reflection:
      s = 2.0*DotProduct3d( &roint->nv, &ray->v )/
              DotProduct3d ( &roint->nv, &roint->nv );
      SetVector3d ( &r, ray->v.x - s*roint->nv.x,
                        ray->v.y - s*roint->nv.y,
                        ray->v.z - s*roint->nv.z );
      SubtractPoints3d ( &roint->p, &rend->rp0, &dp );
      s = det3d ( &r, &dp, &rend->rv2 ) / det3d ( &r, &rend->rv1, &rend->rv2 );
      return s;
  case shapefunc_HIGHLIGHT:
highlight:
      SubtractPoints3d ( &roint->p, &rend->rp0, &dp );
      s = det3d ( &roint->nv, &dp, &rend->rv2 ) /
            det3d ( &roint->nv, &rend->rv1, &rend->rv2 );
      return s;
  case shapefunc_LAMBERTIAN:
lambertian:
      return DotProduct3d ( &rend->lightdir[0], &roint->p );
  case shapefunc_CROSSECTIONS:
crossections:
      return DotProduct3d ( &rend->sectiondir, &roint->p );
  default:
      return 0.0;
    }

default:
    return 0.0;
  }
} /*ShapeFunc*/

double cShapeFunc ( pkRenderer *rend, ray3d *ray, renderobj *obj,
                    RayObjectIntersd *roint )
{
  return ShapeFunc ( rend, rend->c_shape_func, ray, obj, roint );
} /*cShapeFunc*/

boolean dShapeFunc ( pkRenderer *rend, ray3d *ray, renderobj *obj,
                     RayObjectIntersd *roint )
{
  double f;

  f = rend->dfsf*ShapeFunc ( rend, rend->d_shape_func, ray, obj, roint );
  if ( f < 0.0 )
    f = 1.0-f;
  return !(((int)f) & 0x01);
} /*dShapeFunc*/

