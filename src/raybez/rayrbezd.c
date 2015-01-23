
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

#include "raybezprivated.h"

/* ///////////////////////////////////////////////////////////////////////// */
void _rbez_ReflectCPointsRd ( int ncp, const point4d *ctlpoints,
                              const point3d *rayp, const vector3d *nh, double sh,
                              point2d *auxcp )
{
  int     i;
  vector3d v;
  double   w, r;

  for ( i = 0; i < ncp; i++ ) {
    w = ctlpoints[i].w;
    v.x = ctlpoints[i].x - rayp->x*w;
    v.y = ctlpoints[i].y - rayp->y*w;
    v.z = ctlpoints[i].z - rayp->z*w;
    r = (double)(sh * DotProduct3d ( &v, nh ));
    auxcp[i].x = v.x - r*nh->x;
    auxcp[i].y = v.y - r*nh->y;
  }
} /*_rbez_ReflectCPointsRd*/

boolean _rbez_RBezPSolutionOKd ( int object_id, ray3d *ray, point2d *z,
                                 int n, int m,
                                 RBezPatchTreeVertexd *vertex,
                                 int *ninters, RayObjectIntersd *inters,
                                 double *workspace )
{
  vector3d pv;
  double    t;

  if ( z->x >= 0.0 && z->x <= 1.0 && z->y >= 0.0 && z->y <= 1.0 ) {
    inters += *ninters;
    _mbs_BCHornerNvP3Rd ( n, m, vertex->ctlpoints,
                          z->x, z->y, &inters->p, &inters->nv, workspace );
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
} /*_rbez_RBezPSolutionOKd*/

void _rbez_ROutputSingularSolutiond ( int object_id, ray3d *ray,
                                      int n, int m,
                                      RBezPatchTreeVertexd *vertex,
                                      int *ninters, RayObjectIntersd *inters,
                                      double *workspace )
{
  vector3d pv;
  double    t;

  inters += *ninters;
  _mbs_BCHornerNvP3Rd ( n, m, vertex->ctlpoints, 0.5, 0.5,
                        &inters->p, &inters->nv, workspace );
  SubtractPoints3d ( &inters->p, &ray->p, &pv );
  t = DotProduct3d ( &pv, &ray->v );
  if ( t > 0.0 ) {
    inters->t = t;
    inters->u = vertex->u0 + 0.5*(vertex->u1 - vertex->u0);
    inters->v = vertex->v0 + 0.5*(vertex->v1 -vertex-> v0);
    inters->object_id = object_id;
    inters->extra_info = vertex;
    (*ninters)++;
  }
} /*OutputSingularSolutiond*/

/* ///////////////////////////////////////////////////////////////////////// */
int rbez_RayRBezpWspSized ( int n, int m, int maxlevel )
{
  return 2*(n+1)*(m+1)*sizeof(point2d) +
         (maxlevel+1)*sizeof(RBezPatchTreeVertexdp) +
         (6+2*(m+1))*4*sizeof(double);
} /*rbez_RayRBezpWspSized*/

int _rbez_FindRayRBezPatchIntersd ( RBezPatchTreed *tree, ray3d *ray,
                                    int maxlevel, int maxinters,
                                    int *ninters, RayObjectIntersd *inters,
                                    void *workspace )
{
  int                   n, m, ncp;
  point2d               *auxcp;
  double                *wsp;
  int                   stp;
  RBezPatchTreeVertexdp *stk, vertex, left, right;
  vector3d              nh;
  double                sh, K1, K2;
  point2d               p, z;
  vector2d              du, dv;

  n = tree->n;  m = tree->m;
  ncp = (n+1)*(m+1);
  auxcp = (point2d*)workspace;
  stk = (RBezPatchTreeVertexdp*)&auxcp[2*ncp];
  wsp = (double*)&stk[maxlevel+1];
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
      _rbez_ReflectCPointsRd ( ncp, vertex->ctlpoints, &ray->p, &nh, sh, auxcp );
      if ( _rbez_ConvexHullTest2d ( ncp, auxcp ) ) {
        if ( _rbez_UniquenessTest2d ( n, m, ncp, auxcp,
                                      &p, &du, &dv, &K1, &K2, wsp ) ) {
          switch ( _rbez_NewtonMethod2d ( n, m, auxcp,
                                          &p, &du, &dv, &z, wsp ) ) {
        case RBEZ_NEWTON_YES:
            if ( !_rbez_RBezPSolutionOKd ( tree->object_id, ray, &z,
                      n, m, vertex, ninters, inters, wsp ) ) {
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
            left  = rbez_GetRBezLeftVertexd  ( tree, vertex );
            right = rbez_GetRBezRightVertexd ( tree, vertex );
            stk[stp++] = left;
            stk[stp++] = right;
          }
          else
            _rbez_ROutputSingularSolutiond ( tree->object_id, ray,
                          n, m, vertex, ninters, inters, wsp );
        }
      }
    }
  } while ( stp && *ninters < maxinters );
  return *ninters;

failure:
  return -1;
} /*_rbez_FindRayRBezPatchIntersd*/

int rbez_FindRayRBezPatchIntersd ( RBezPatchTreed *tree, ray3d *ray,
                                   int maxlevel, int maxinters,
                                   int *ninters, RayObjectIntersd *inters )
{
  void *sp;
  int  result;

  sp = pkv_GetScratchMem ( rbez_RayRBezpWspSized ( tree->n, tree->m, maxlevel ) );
  if ( !sp )
    return -1;
  result = _rbez_FindRayRBezPatchIntersd ( tree, ray, maxlevel, maxinters,
                                           ninters, inters, sp );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*rbez_FindRayRBezPatchIntersd*/

