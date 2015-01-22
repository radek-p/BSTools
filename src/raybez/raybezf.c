
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

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "raybez.h"

#include "raybezprivatef.h"

/* ///////////////////////////////////////////////////////////////////////// */
void _rbez_ReflectCPointsf ( int ncp, const point3f *ctlpoints,
                             const point3f *rayp, const vector3f *nh, float sh,
                             point2f *auxcp )
{
  int      i;
  vector3f v;
  float    r;

  for ( i = 0; i < ncp; i++ ) {
    SubtractPoints3f ( &ctlpoints[i], rayp, &v );
    r = sh * DotProduct3f ( &v, nh );
    auxcp[i].x = v.x - r*nh->x;
    auxcp[i].y = v.y - r*nh->y;
  }
} /*_rbez_ReflectCPointsf*/

boolean _rbez_BezPSolutionOKf ( int object_id, ray3f *ray, point2f *z,
                                int n, int m,
                                BezPatchTreeVertexf *vertex,
                                int *ninters, RayObjectIntersf *inters )
{
  vector3f pv;
  float    t;

  if ( z->x >= 0.0 && z->x <= 1.0 && z->y >= 0.0 && z->y <= 1.0 ) {
    inters += *ninters;
    mbs_BCHornerNvP3f ( n, m, vertex->ctlpoints,
                        z->x, z->y, &inters->p, &inters->nv );
    SubtractPoints3f ( &inters->p, &ray->p, &pv );
    t = DotProduct3f ( &pv, &ray->v );
    if ( t > 0.0 ) {
      inters->t = t;
      inters->u = vertex->u0 + z->x*(vertex->u1 - vertex->u0);
      inters->v = vertex->v0 + z->y*(vertex->v1 - vertex->v0);
      inters->object_id = object_id;
      inters->extra_info = vertex;
      (*ninters)++;
    }
    return true;
  }
  else
    return false;
} /*_rbez_BezPSolutionOKf*/

void _rbez_OutputSingularSolutionf ( int object_id, ray3f *ray,
                                     int n, int m,
                                     BezPatchTreeVertexf *vertex,
                                     int *ninters, RayObjectIntersf *inters )
{
  vector3f pv;
  float    t;

  inters += *ninters;
  mbs_BCHornerNvP3f ( n, m, vertex->ctlpoints, 0.5, 0.5, &inters->p, &inters->nv );
  SubtractPoints3f ( &inters->p, &ray->p, &pv );
  t = DotProduct3f ( &pv, &ray->v );
  if ( t > 0.0 ) {
    inters->t = t;
    inters->u = 0.5*(vertex->u0 + vertex->u1);
    inters->v = 0.5*(vertex->v0 + vertex->v1);
    inters->object_id = object_id;
    inters->extra_info = vertex;
    (*ninters)++;
  }
} /*OutputSingularSolutionf*/

/* ///////////////////////////////////////////////////////////////////////// */
int rbez_FindRayBezPatchIntersf ( BezPatchTreef *tree, ray3f *ray,
                                  int maxlevel, int maxinters,
                                  int *ninters, RayObjectIntersf *inters )
{
  void                 *sp;
  int                  n, m, ncp;
  point2f              *auxcp;
  int                  stp;
  BezPatchTreeVertexfp *stk, vertex, left, right;
  vector3f             nh;
  float                sh, K1, K2;
  point2f              p, z;
  vector2f             du, dv;

  sp = pkv_GetScratchMemTop ();
  n = tree->n;  m = tree->m;
  ncp = (n+1)*(m+1);
  auxcp = pkv_GetScratchMem ( 2*ncp*sizeof(point2f) );
  stk = pkv_GetScratchMem ( (maxlevel+1)*sizeof(BezPatchTreeVertexfp) );
  if ( !auxcp || !stk )
    goto failure;
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
      _rbez_ReflectCPointsf ( ncp, vertex->ctlpoints, &ray->p, &nh, sh, auxcp );
      if ( _rbez_ConvexHullTest2f ( ncp, auxcp ) ) {
        if ( _rbez_UniquenessTest2f ( n, m, ncp, auxcp, &p, &du, &dv, &K1, &K2 ) ) {
          switch ( _rbez_NewtonMethod2f ( n, m, auxcp, &p, &du, &dv, &z ) ) {
        case RBEZ_NEWTON_YES:
            if ( !_rbez_BezPSolutionOKf ( tree->object_id, ray, &z,
                      n, m, vertex, ninters, inters ) ) {
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
            left  = rbez_GetBezLeftVertexf  ( tree, vertex );
            right = rbez_GetBezRightVertexf ( tree, vertex );
            stk[stp++] = left;
            stk[stp++] = right;
          }
          else
            _rbez_OutputSingularSolutionf ( tree->object_id, ray,
                                            n, m, vertex, ninters, inters );
        }
      }
    }
  } while ( stp && *ninters < maxinters );

  pkv_SetScratchMemTop ( sp );
  return *ninters;

failure:
  pkv_SetScratchMemTop ( sp );
  return -1;
} /*rbez_FindRayBezPatchIntersf*/

