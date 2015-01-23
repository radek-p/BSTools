
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
void _rbez_ReflectCPointsRf ( int ncp, const point4f *ctlpoints,
                              const point3f *rayp, const vector3f *nh, float sh,
                              point2f *auxcp )
{
  int     i;
  vector3f v;
  float   w, r;

  for ( i = 0; i < ncp; i++ ) {
    w = ctlpoints[i].w;
    v.x = ctlpoints[i].x - rayp->x*w;
    v.y = ctlpoints[i].y - rayp->y*w;
    v.z = ctlpoints[i].z - rayp->z*w;
    r = (float)(sh * DotProduct3f ( &v, nh ));
    auxcp[i].x = v.x - r*nh->x;
    auxcp[i].y = v.y - r*nh->y;
  }
} /*_rbez_ReflectCPointsRf*/

boolean _rbez_RBezPSolutionOKf ( int object_id, ray3f *ray, point2f *z,
                                 int n, int m,
                                 RBezPatchTreeVertexf *vertex,
                                 int *ninters, RayObjectIntersf *inters,
                                 float *workspace )
{
  vector3f pv;
  float    t;

  if ( z->x >= 0.0 && z->x <= 1.0 && z->y >= 0.0 && z->y <= 1.0 ) {
    inters += *ninters;
    _mbs_BCHornerNvP3Rf ( n, m, vertex->ctlpoints,
                          z->x, z->y, &inters->p, &inters->nv, workspace );
    SubtractPoints3f ( &inters->p, &ray->p, &pv );
    t = DotProduct3f ( &pv, &ray->v );
    if ( t > 0.0 ) {
      inters->t = t;
      inters->u = vertex->u0 + z->x*(vertex->u1 - vertex->u0);
      inters->v = vertex->v0 + z->y*(vertex->v1 -vertex-> v0);
      inters->object_id = object_id;
      inters->extra_info = vertex;
      (*ninters)++;
    }
    return true;
  }
  else
    return false;
} /*_rbez_RBezPSolutionOKf*/

void _rbez_ROutputSingularSolutionf ( int object_id, ray3f *ray,
                                      int n, int m,
                                      RBezPatchTreeVertexf *vertex,
                                      int *ninters, RayObjectIntersf *inters,
                                      float *workspace )
{
  vector3f pv;
  float    t;

  inters += *ninters;
  _mbs_BCHornerNvP3Rf ( n, m, vertex->ctlpoints, 0.5, 0.5,
                        &inters->p, &inters->nv, workspace );
  SubtractPoints3f ( &inters->p, &ray->p, &pv );
  t = DotProduct3f ( &pv, &ray->v );
  if ( t > 0.0 ) {
    inters->t = t;
    inters->u = vertex->u0 + 0.5*(vertex->u1 - vertex->u0);
    inters->v = vertex->v0 + 0.5*(vertex->v1 -vertex-> v0);
    inters->object_id = object_id;
    inters->extra_info = vertex;
    (*ninters)++;
  }
} /*OutputSingularSolutionf*/

/* ///////////////////////////////////////////////////////////////////////// */
int rbez_RayRBezpWspSizef ( int n, int m, int maxlevel )
{
  return 2*(n+1)*(m+1)*sizeof(point2f) +
         (maxlevel+1)*sizeof(RBezPatchTreeVertexfp) +
         (6+2*(m+1))*4*sizeof(float);
} /*rbez_RayRBezpWspSizef*/

int _rbez_FindRayRBezPatchIntersf ( RBezPatchTreef *tree, ray3f *ray,
                                    int maxlevel, int maxinters,
                                    int *ninters, RayObjectIntersf *inters,
                                    void *workspace )
{
  int                   n, m, ncp;
  point2f               *auxcp;
  float                 *wsp;
  int                   stp;
  RBezPatchTreeVertexfp *stk, vertex, left, right;
  vector3f              nh;
  float                 sh, K1, K2;
  point2f               p, z;
  vector2f              du, dv;

  n = tree->n;  m = tree->m;
  ncp = (n+1)*(m+1);
  auxcp = (point2f*)workspace;
  stk = (RBezPatchTreeVertexfp*)&auxcp[2*ncp];
  wsp = (float*)&stk[maxlevel+1];
        /* construct the Householder reflection */
  nh = ray->v;
  if ( nh.z > 0.0 ) nh.z += 1.0;
               else nh.z -= 1.0;
  sh = 2.0/DotProduct3f ( &nh, &nh );

  *ninters = 0;
  stp = 0;
  stk[stp++] = tree->root;
  do {
    vertex = stk[--stp];
    if ( rbez_TestRayBBoxf ( ray, &vertex->bbox ) ) {
      if ( vertex->ctlpoints == NULL )
        goto DIVIDE;
      _rbez_ReflectCPointsRf ( ncp, vertex->ctlpoints, &ray->p, &nh, sh, auxcp );
      if ( _rbez_ConvexHullTest2f ( ncp, auxcp ) ) {
        if ( _rbez_UniquenessTest2f ( n, m, ncp, auxcp,
                                      &p, &du, &dv, &K1, &K2, wsp ) ) {
          switch ( _rbez_NewtonMethod2f ( n, m, auxcp,
                                          &p, &du, &dv, &z, wsp ) ) {
        case RBEZ_NEWTON_YES:
            if ( !_rbez_RBezPSolutionOKf ( tree->object_id, ray, &z,
                      n, m, vertex, ninters, inters, wsp ) ) {
              if ( _rbez_SecondTest2f ( &z, n, m, K1, K2 ) )
                goto DIVIDE;
            }
            break;
        case RBEZ_NEWTON_NO:
            goto DIVIDE;
        case RBEZ_NEWTON_ERROR:
            goto failure;
          }
        }
        else {
DIVIDE:
          if ( vertex->level < maxlevel ) {
                        /* divide and push pieces */
            left  = rbez_GetRBezLeftVertexf  ( tree, vertex );
            right = rbez_GetRBezRightVertexf ( tree, vertex );
            stk[stp++] = left;
            stk[stp++] = right;
          }
          else
            _rbez_ROutputSingularSolutionf ( tree->object_id, ray,
                          n, m, vertex, ninters, inters, wsp );
        }
      }
    }
  } while ( stp && *ninters < maxinters );
  return *ninters;

failure:
  return -1;
} /*_rbez_FindRayRBezPatchIntersf*/

int rbez_FindRayRBezPatchIntersf ( RBezPatchTreef *tree, ray3f *ray,
                                   int maxlevel, int maxinters,
                                   int *ninters, RayObjectIntersf *inters )
{
  void *sp;
  int  result;

  sp = pkv_GetScratchMem ( rbez_RayRBezpWspSizef ( tree->n, tree->m, maxlevel ) );
  if ( !sp )
    return -1;
  result = _rbez_FindRayRBezPatchIntersf ( tree, ray, maxlevel, maxinters,
                                           ninters, inters, sp );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*rbez_FindRayRBezPatchIntersf*/

