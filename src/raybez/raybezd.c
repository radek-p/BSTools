
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2014                            */ 
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

#include "raybezprivated.h"

/* ///////////////////////////////////////////////////////////////////////// */
void _rbez_ReflectCPointsd ( int ncp, const point3d *ctlpoints,
                             const point3d *rayp, const vector3d *nh, double sh,
                             point2d *auxcp )
{
  int      i;
  vector3d v;
  double   r;

  for ( i = 0; i < ncp; i++ ) {
    SubtractPoints3d ( &ctlpoints[i], rayp, &v );
    r = sh * DotProduct3d ( &v, nh );
    auxcp[i].x = v.x - r*nh->x;
    auxcp[i].y = v.y - r*nh->y;
  }
} /*_rbez_ReflectCPointsd*/

boolean _rbez_BezPSolutionOKd ( int object_id, ray3d *ray, point2d *z,
                                int n, int m, point3d *cpoints,
                                BezPatchTreeVertexd *vertex,
                                int *ninters, RayObjectIntersd *inters )
{
  vector3d pv;
  double   t;

  if ( z->x >= 0.0 && z->x <= 1.0 && z->y >= 0.0 && z->y <= 1.0 ) {
    inters += *ninters;
    mbs_BCHornerNvP3d ( n, m, cpoints,
                        z->x, z->y, &inters->p, &inters->nv );
    SubtractPoints3d ( &inters->p, &ray->p, &pv );
    t = DotProduct3d ( &pv, &ray->v );
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
} /*_rbez_BezPSolutionOKd*/

void _rbez_OutputSingularSolutiond ( int object_id, ray3d *ray,
                                     int n, int m, point3d *cpoints,
                                     BezPatchTreeVertexd *vertex,
                                     int *ninters, RayObjectIntersd *inters )
{
  vector3d pv;
  double   t;

  inters += *ninters;
  mbs_BCHornerNvP3d ( n, m, cpoints, 0.5, 0.5, &inters->p, &inters->nv );
  SubtractPoints3d ( &inters->p, &ray->p, &pv );
  t = DotProduct3d ( &pv, &ray->v );
  if ( t > 0.0 ) {
    inters->t = t;
    inters->u = 0.5*(vertex->u0 + vertex->u1);
    inters->v = 0.5*(vertex->v0 + vertex->v1);
    inters->object_id = object_id;
    inters->extra_info = vertex;
    (*ninters)++;
  }
} /*OutputSingularSolutiond*/

/* ///////////////////////////////////////////////////////////////////////// */
int rbez_FindRayBezPatchIntersd ( BezPatchTreed *tree, ray3d *ray,
                                  int maxlevel, int maxinters,
                                  int *ninters, RayObjectIntersd *inters )
{
  void                 *sp;
  int                  n, m, ncp;
  point2d              *auxcp;
  int                  stp;
  BezPatchTreeVertexdp *stk, vertex, left, right;
  vector3d             nh;
  double               sh, K1, K2;
  point2d              p, z;
  vector2d             du, dv;

  sp = pkv_GetScratchMemTop ();
  n = tree->n;  m = tree->m;
  ncp = (n+1)*(m+1);
  auxcp = pkv_GetScratchMem ( 2*ncp*sizeof(point2d) );
  stk = pkv_GetScratchMem ( (maxlevel+1)*sizeof(BezPatchTreeVertexdp) );
  if ( !auxcp || !stk )
    goto failure;
        /* construct the Householder reflection */
  nh = ray->v;
  if ( nh.z > 0.0 ) nh.z += 1.0;
               else nh.z -= 1.0;
  sh = 2.0/DotProduct3d ( &nh, &nh );

  *ninters = 0;
  stp = 0;
  stk[stp++] = tree->root;
  do {
    vertex = stk[--stp];
    if ( rbez_TestRayBBoxd ( ray, &vertex->bbox ) ) {
      if ( vertex->ctlpoints == NULL )
        goto DIVIDE;
      _rbez_ReflectCPointsd ( ncp, vertex->ctlpoints, &ray->p, &nh, sh, auxcp );
      if ( _rbez_ConvexHullTest2d ( ncp, auxcp ) ) {
        if ( _rbez_UniquenessTest2d ( n, m, ncp, auxcp, &p, &du, &dv, &K1, &K2 ) ) {
          switch ( _rbez_NewtonMethod2d ( n, m, auxcp, &p, &du, &dv, &z ) ) {
        case RBEZ_NEWTON_YES:
            if ( !_rbez_BezPSolutionOKd ( tree->object_id, ray, &z,
                      n, m, vertex->ctlpoints, vertex,
                      ninters, inters ) ) {
              if ( _rbez_SecondTest2d ( &z, n, m, K1, K2 ) )
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
            left  = rbez_GetBezLeftVertexd  ( tree, vertex );
            right = rbez_GetBezRightVertexd ( tree, vertex );
            stk[stp++] = left;
            stk[stp++] = right;
          }
          else
            _rbez_OutputSingularSolutiond ( tree->object_id, ray,
                          n, m, vertex->ctlpoints, vertex, ninters, inters );
        }
      }
    }
  } while ( stp && *ninters < maxinters );

  pkv_SetScratchMemTop ( sp );
  return *ninters;

failure:
  pkv_SetScratchMemTop ( sp );
  return -1;
} /*rbez_FindRayBezPatchIntersd*/

