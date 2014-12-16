
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2014                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* Changes: */
/* 26.06.2005, P. Kiciak - added some tolerance EPS to the ClipTestd function */
/* 15.09.2008, P. Kiciak - added some tolerance EPS to the box dimensions */
/* 22.10.2010, P. Kiciak - added the object_id attribute */
/* 27.03.2012, P. Kiciak - moved out the ray/box intersection procedures */
/*  8.07.2012, P. Kiciak - added locking/unlocking a mutex for pthreads */
/* 31.08.2012, P. Kiciak - added the leaf attribute, to avoid */
/*                         unnecessary mutex locking */
/* 26.09.2012, P. Kiciak - added the normalvect attribute of the vertex */
/* 28.01.2013, P. Kiciak - moved auxiliary procedures (solution existence
                           and uniqueness tests, Newton method)
                           to a separate file */

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
static RBezPatchTreeVertexdp
  AllocRTreeVertexd ( RBezPatchTreedp tree,
                     RBezPatchTreeVertexdp up,
                     double u0, double u1, double v0, double v1 )
{
  RBezPatchTreeVertexdp vertex;

  PKV_MALLOC ( vertex, sizeof(RBezPatchTreeVertexd) + tree->cpsize );
  if ( vertex ) {
    vertex->u0 = u0;  vertex->u1 = u1;
    vertex->v0 = v0;  vertex->v1 = v1;
    vertex->left = vertex->right = NULL;
    vertex->up = up;
    vertex->leaf = true;
    if ( up )
      vertex->level = (short)(up->level + 1);
    else
      vertex->level = 0;
    vertex->ctlpoints = (point4d*)((char*)vertex+sizeof(RBezPatchTreeVertexd));
    vertex->normalvect = NULL;
  }
  return vertex;
} /*AllocRTreeVertexd*/

static void FindRBoundingBoxd ( RBezPatchTreedp tree,
                                RBezPatchTreeVertexdp vertex )
{
#define EPS 1.0e-10
  int      n, m, ncp;
  point4d  *cp, *cq;
  vector4d v;
  int      i, j;
  double    du, dv, d;

  n = tree->n;  m = tree->m;
  ncp = (n+1)*(m+1);
                                     /* find the bounding box */
  rbez_FindCPBoundingBox3Rd ( 1, ncp, 0, vertex->ctlpoints, EPS, &vertex->bbox );
                                     /* determine division direction */
  du = 0.0;
  for ( i = 0, cp = vertex->ctlpoints, cq = cp+(m+1);
        i < n*(m+1);
        i++, cp++, cq++ ) {
    SubtractPoints4d ( cq, cp, &v );
    d = DotProduct4d ( &v, &v );
    du = max ( du, d );
  }
  du = (double)n * sqrt(du);
  dv = 0.0;
  for ( i = 0, cp = vertex->ctlpoints;
        i <= n;
        i++, cp++ )
    for ( j = 0; j < m; j++, cp++ ) {
      SubtractPoints4d ( cp+1, cp, &v );
      d = DotProduct4d ( &v, &v );
      dv = max ( dv, d );
    }
  dv = (double)m * sqrt(dv);

  if ( du > dv ) {
    vertex->divdir = 0;
    vertex->maxder = du;
  }
  else {
    vertex->divdir = 1;
    vertex->maxder = dv;
  }
#undef EPS
} /*FindRBoundingBoxd*/

static void UpdateRBoundingBoxesd ( RBezPatchTreeVertexdp vertex )
{
  while ( vertex ) {
    if ( rbez_NarrowBBoxSumd ( &vertex->left->bbox, &vertex->right->bbox,
                               &vertex->bbox ) )
      vertex = vertex->up;
    else
      return;
  }
} /*UpdateRBoundingBoxesd*/

RBezPatchTreedp
  rbez_NewRBezPatchTreed ( int object_id,
                           unsigned char n, unsigned char m,
                           double u0, double u1, double v0, double v1,
                           CONST_ point4d *ctlpoints )
{
  RBezPatchTreedp       tree;

  PKV_MALLOC ( tree, sizeof(RBezPatchTreed) );
  if ( tree ) {
    tree->object_id = object_id;
    tree->n = n;    tree->m = m;
    tree->cpsize = (n+1)*(m+1)*sizeof(point4d);
    tree->root = AllocRTreeVertexd ( tree, NULL, u0, u1, v0, v1 );
    if ( tree->root ) {
      memcpy ( tree->root->ctlpoints, ctlpoints, tree->cpsize );
      mbs_BCHornerNvP3Rd ( n, m, ctlpoints, 0.5, 0.5,
                           &tree->root->pcent, &tree->root->nvcent );
      FindRBoundingBoxd ( tree, tree->root );
    }
    else {
      PKV_FREE ( tree );
      return NULL;
    }
  }
  return tree;
} /*rbez_NewRBezPatchTreed*/

static void r_DestroyRTreed ( RBezPatchTreeVertexdp vertex )
{
  if ( vertex->right )
    r_DestroyRTreed ( vertex->right );
  if ( vertex->left )
    r_DestroyRTreed ( vertex->left );
  if ( vertex->normalvect )
    PKV_FREE ( vertex->normalvect );
  PKV_FREE ( vertex );
} /*r_DestroyRTreed*/

void rbez_DestroyRBezPatchTreed ( RBezPatchTreedp tree )
{
  r_DestroyRTreed ( tree->root );
  PKV_FREE ( tree );
} /*rbez_DestroyRBezPatchTreed*/

static void DivideRVertexd ( RBezPatchTreedp tree, RBezPatchTreeVertexdp vertex )
{
  double h;

        /* allocate subtree roots */
  if ( !vertex->divdir ) {
    h = 0.5*(vertex->u0 + vertex->u1);
    vertex->left  = AllocRTreeVertexd ( tree, vertex,
                                        vertex->u0, h, vertex->v0, vertex->v1 );
    vertex->right = AllocRTreeVertexd ( tree, vertex,
                                        h, vertex->u1, vertex->v0, vertex->v1 );
    if ( !vertex->left || !vertex->right )
      goto dealloc;
    memcpy ( vertex->right->ctlpoints, vertex->ctlpoints, tree->cpsize );
    mbs_BisectBP4ud ( tree->n, tree->m,
                      vertex->right->ctlpoints, vertex->left->ctlpoints );
  }
  else {
    h = 0.5*(vertex->v0 + vertex->v1);
    vertex->left  = AllocRTreeVertexd ( tree, vertex,
                                        vertex->u0, vertex->u1, vertex->v0, h );
    vertex->right = AllocRTreeVertexd ( tree, vertex,
                                        vertex->u0, vertex->u1, h, vertex->v1 );
    if ( !vertex->left || !vertex->right ) {
dealloc:
      if ( vertex->left )  PKV_FREE ( vertex->left );
      if ( vertex->right ) PKV_FREE ( vertex->right );
      return;
    }
    memcpy ( vertex->right->ctlpoints, vertex->ctlpoints, tree->cpsize );
    mbs_BisectBP4vd ( tree->n, tree->m,
                      vertex->right->ctlpoints, vertex->left->ctlpoints );
  }
  mbs_BCHornerNvP3Rd ( tree->n, tree->m, vertex->left->ctlpoints,
                       0.5, 0.5, &vertex->left->pcent, &vertex->left->nvcent );
  mbs_BCHornerNvP3Rd ( tree->n, tree->m, vertex->right->ctlpoints,
                       0.5, 0.5, &vertex->right->pcent, &vertex->right->nvcent );
  FindRBoundingBoxd ( tree, vertex->left );
  FindRBoundingBoxd ( tree, vertex->right );
  UpdateRBoundingBoxesd ( vertex );
  vertex->leaf = false;
} /*DivideRVertexd*/

RBezPatchTreeVertexdp
  rbez_GetRBezLeftVertexd ( RBezPatchTreedp tree,
                            RBezPatchTreeVertexdp vertex )
{
  if ( vertex->leaf ) {
    if ( raybez_use_mutex )
      pthread_mutex_lock ( &raybez_mutex );
    if ( !vertex->left )
      DivideRVertexd ( tree, vertex );
    if ( raybez_use_mutex )
      pthread_mutex_unlock ( &raybez_mutex );
  }
  return vertex->left;
} /*rbez_GetRBezLeftVertexd*/

RBezPatchTreeVertexdp
  rbez_GetRBezRightVertexd ( RBezPatchTreedp tree,
                         RBezPatchTreeVertexdp vertex )
{
  if ( vertex->leaf ) {
    if ( raybez_use_mutex )
      pthread_mutex_lock ( &raybez_mutex );
    if ( !vertex->right )
      DivideRVertexd ( tree, vertex );
    if ( raybez_use_mutex )
      pthread_mutex_unlock ( &raybez_mutex );
  }
  return vertex->right;
} /*rbez_GetRBezRightVertexd*/

/* ///////////////////////////////////////////////////////////////////////// */
static void ConvertPatchd ( int ncp, const point4d *ctlpoints,
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
    r = sh * DotProduct3d ( &v, nh );
    auxcp[i].x = v.x - r*nh->x;
    auxcp[i].y = v.y - r*nh->y;
  }
} /*ConvertPatchd*/

static char SolutionOKd ( ray3d *ray, RBezPatchTreed *tree,
                         RBezPatchTreeVertexd *vertex, point2d *z,
                         int *ninters, RayObjectIntersd *inters )
{
  double   u, v;
  vector3d pv;

  if ( z->x >= 0.0 && z->x <= 1.0 && z->y >= 0.0 && z->y <= 1.0 ) {
    inters += *ninters;
    u = inters->u = vertex->u0 + z->x*(vertex->u1 - vertex->u0);
    u = (u - tree->root->u0)/(tree->root->u1 - tree->root->u0);
    v = inters->v = vertex->v0 + z->y*(vertex->v1 - vertex->v0);
    v = (v - tree->root->v0)/(tree->root->v1 - tree->root->v0);

    mbs_BCHornerNvP3Rd ( tree->n, tree->m, tree->root->ctlpoints,
                         u, v, &inters->p, &inters->nv );
    SubtractPoints3d ( &inters->p, &ray->p, &pv );
    inters->t = DotProduct3d ( &pv, &ray->v );
    if ( inters->t > 0.0 ) {
      inters->object_id = tree->object_id;
      (*ninters)++;
    }
    return true;
  }
  else
    return false;
} /*SolutionOKd*/

static void OutputSingularSolutiond ( ray3d *ray, RBezPatchTreed *tree,
                                      RBezPatchTreeVertexd *vertex,
                                      int *ninters, RayObjectIntersd *inters )
{
  double    u, v;
  vector3d pv;

  inters += *ninters;
  inters->u = 0.5*(vertex->u0 + vertex->u1);
  u = (inters->u - tree->root->u0)/(tree->root->u1 - tree->root->u0);
  inters->v = 0.5*(vertex->v0 + vertex->v1);
  v = (inters->v - tree->root->v0)/(tree->root->v1 - tree->root->v0);
  mbs_BCHornerNvP3Rd ( tree->n, tree->m, tree->root->ctlpoints,
                       u, v, &inters->p, &inters->nv );
  SubtractPoints3d ( &inters->p, &ray->p, &pv );
  inters->t = DotProduct3d ( &pv, &ray->v );
  if ( inters->t > 0.0 ) {
    inters->object_id = tree->object_id;
    (*ninters)++;
  }
} /*OutputSingularSolutiond*/

int rbez_FindRayRBezPatchIntersd ( RBezPatchTreed *tree, ray3d *ray,
                                   int maxlevel, int maxinters,
                                   int *ninters, RayObjectIntersd *inters )
{
  int      n, m, ncp;
  int      size_auxcp, size_stack;
  int      sp;           /* stack pointer */
  RBezPatchTreeVertexdp *stack;
  RBezPatchTreeVertexdp vertex, left, right;
  point2d  *auxcp;
  vector3d nh;           /* Householder reflection vector */
  double    sh;           /* scaling coefficient */
  double    K1, K2;       /* uniqueness test numbers */
  point2d  p, z;
  vector2d du, dv;

  n = tree->n;  m = tree->m;
  ncp = (n+1)*(m+1);
  stack = (RBezPatchTreeVertexdp*)pkv_GetScratchMem
              ( size_stack = maxlevel*sizeof(RBezPatchTreeVertexdp) );
                        /* the auxcp array is twice longer for the */
                        /* needs of solution uniqueness test */
  auxcp = (point2d*)pkv_GetScratchMem ( size_auxcp = 2*ncp*sizeof(point2d) );

  if ( !auxcp || !stack )
    goto failure;

                        /* construct the Householder reflection */
  nh = ray->v;
  if ( nh.z > 0 ) nh.z += 1.0;
    else          nh.z -= 1.0;
  sh = 2.0/DotProduct3d ( &nh, &nh );

                        /* initialize stack, push the patch and go around */
  *ninters = 0;
  sp = 0;
  stack[sp++] = tree->root;
  do {
    vertex = stack[--sp];
    if ( rbez_TestRayBBoxd ( ray, &vertex->bbox ) ) {
      ConvertPatchd ( ncp, vertex->ctlpoints, &ray->p, &nh, sh, auxcp );
      if ( _rbez_ConvexHullTest2d ( ncp, auxcp ) ) {
        if ( _rbez_UniquenessTest2d ( n, m, ncp, auxcp, &p, &du, &dv, &K1, &K2 ) ) {
          switch ( _rbez_NewtonMethod2d ( n, m, auxcp, &p, &du, &dv, &z ) ) {
        case RBEZ_NEWTON_YES:
            if ( !SolutionOKd ( ray, tree, vertex, &z, ninters, inters ) ) {
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
            stack[sp++] = left;
            stack[sp++] = right;
          }
          else
            OutputSingularSolutiond ( ray, tree, vertex, ninters, inters );
        }
      }
    }
  } while ( sp && *ninters < maxinters );

  pkv_FreeScratchMem ( size_auxcp+size_stack );
  return *ninters;

failure:
  pkv_FreeScratchMem ( size_auxcp+size_stack );
  return -1;
} /*FindRayRBezPatchIntersd*/

