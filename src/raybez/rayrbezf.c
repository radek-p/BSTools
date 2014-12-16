
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2014                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* Changes: */
/* 26.06.2005, P. Kiciak - added some tolerance EPS to the ClipTestf function */
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

#include "raybezprivatef.h"

/* ///////////////////////////////////////////////////////////////////////// */
static RBezPatchTreeVertexfp
  AllocRTreeVertexf ( RBezPatchTreefp tree,
                     RBezPatchTreeVertexfp up,
                     float u0, float u1, float v0, float v1 )
{
  RBezPatchTreeVertexfp vertex;

  PKV_MALLOC ( vertex, sizeof(RBezPatchTreeVertexf) + tree->cpsize );
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
    vertex->ctlpoints = (point4f*)((char*)vertex+sizeof(RBezPatchTreeVertexf));
    vertex->normalvect = NULL;
  }
  return vertex;
} /*AllocRTreeVertexf*/

static void FindRBoundingBoxf ( RBezPatchTreefp tree,
                                RBezPatchTreeVertexfp vertex )
{
#define EPS 5.0e-6
  int      n, m, ncp;
  point4f  *cp, *cq;
  vector4f v;
  int      i, j;
  float    du, dv, d;

  n = tree->n;  m = tree->m;
  ncp = (n+1)*(m+1);
                                     /* find the bounding box */
  rbez_FindCPBoundingBox3Rf ( 1, ncp, 0, vertex->ctlpoints, EPS, &vertex->bbox );
                                     /* determine division direction */
  du = 0.0;
  for ( i = 0, cp = vertex->ctlpoints, cq = cp+(m+1);
        i < n*(m+1);
        i++, cp++, cq++ ) {
    SubtractPoints4f ( cq, cp, &v );
    d = (float)DotProduct4f ( &v, &v );
    du = max ( du, d );
  }
  du = (float)(n * sqrt(du));
  dv = 0.0;
  for ( i = 0, cp = vertex->ctlpoints;
        i <= n;
        i++, cp++ )
    for ( j = 0; j < m; j++, cp++ ) {
      SubtractPoints4f ( cp+1, cp, &v );
      d = (float)DotProduct4f ( &v, &v );
      dv = max ( dv, d );
    }
  dv = (float)(m * sqrt(dv));

  if ( du > dv ) {
    vertex->divdir = 0;
    vertex->maxder = du;
  }
  else {
    vertex->divdir = 1;
    vertex->maxder = dv;
  }
#undef EPS
} /*FindRBoundingBoxf*/

static void UpdateRBoundingBoxesf ( RBezPatchTreeVertexfp vertex )
{
  while ( vertex ) {
    if ( rbez_NarrowBBoxSumf ( &vertex->left->bbox, &vertex->right->bbox,
                               &vertex->bbox ) )
      vertex = vertex->up;
    else
      return;
  }
} /*UpdateRBoundingBoxesf*/

RBezPatchTreefp
  rbez_NewRBezPatchTreef ( int object_id,
                           unsigned char n, unsigned char m,
                           float u0, float u1, float v0, float v1,
                           CONST_ point4f *ctlpoints )
{
  RBezPatchTreefp       tree;

  PKV_MALLOC ( tree, sizeof(RBezPatchTreef) );
  if ( tree ) {
    tree->object_id = object_id;
    tree->n = n;    tree->m = m;
    tree->cpsize = (n+1)*(m+1)*sizeof(point4f);
    tree->root = AllocRTreeVertexf ( tree, NULL, u0, u1, v0, v1 );
    if ( tree->root ) {
      memcpy ( tree->root->ctlpoints, ctlpoints, tree->cpsize );
      mbs_BCHornerNvP3Rf ( n, m, ctlpoints, 0.5, 0.5,
                           &tree->root->pcent, &tree->root->nvcent );
      FindRBoundingBoxf ( tree, tree->root );
    }
    else {
      PKV_FREE ( tree );
      return NULL;
    }
  }
  return tree;
} /*rbez_NewRBezPatchTreef*/

static void r_DestroyRTreef ( RBezPatchTreeVertexfp vertex )
{
  if ( vertex->right )
    r_DestroyRTreef ( vertex->right );
  if ( vertex->left )
    r_DestroyRTreef ( vertex->left );
  if ( vertex->normalvect )
    PKV_FREE ( vertex->normalvect );
  PKV_FREE ( vertex );
} /*r_DestroyRTreef*/

void rbez_DestroyRBezPatchTreef ( RBezPatchTreefp tree )
{
  r_DestroyRTreef ( tree->root );
  PKV_FREE ( tree );
} /*rbez_DestroyRBezPatchTreef*/

static void DivideRVertexf ( RBezPatchTreefp tree, RBezPatchTreeVertexfp vertex )
{
  float h;

        /* allocate subtree roots */
  if ( !vertex->divdir ) {
    h = (float)(0.5*(vertex->u0 + vertex->u1));
    vertex->left  = AllocRTreeVertexf ( tree, vertex,
                                        vertex->u0, h, vertex->v0, vertex->v1 );
    vertex->right = AllocRTreeVertexf ( tree, vertex,
                                        h, vertex->u1, vertex->v0, vertex->v1 );
    if ( !vertex->left || !vertex->right )
      goto dealloc;
    memcpy ( vertex->right->ctlpoints, vertex->ctlpoints, tree->cpsize );
    mbs_BisectBP4uf ( tree->n, tree->m,
                      vertex->right->ctlpoints, vertex->left->ctlpoints );
  }
  else {
    h = (float)(0.5*(vertex->v0 + vertex->v1));
    vertex->left  = AllocRTreeVertexf ( tree, vertex,
                                        vertex->u0, vertex->u1, vertex->v0, h );
    vertex->right = AllocRTreeVertexf ( tree, vertex,
                                        vertex->u0, vertex->u1, h, vertex->v1 );
    if ( !vertex->left || !vertex->right ) {
dealloc:
      if ( vertex->left )  PKV_FREE ( vertex->left );
      if ( vertex->right ) PKV_FREE ( vertex->right );
      return;
    }
    memcpy ( vertex->right->ctlpoints, vertex->ctlpoints, tree->cpsize );
    mbs_BisectBP4vf ( tree->n, tree->m,
                      vertex->right->ctlpoints, vertex->left->ctlpoints );
  }
  mbs_BCHornerNvP3Rf ( tree->n, tree->m, vertex->left->ctlpoints,
                       0.5, 0.5, &vertex->left->pcent, &vertex->left->nvcent );
  mbs_BCHornerNvP3Rf ( tree->n, tree->m, vertex->right->ctlpoints,
                       0.5, 0.5, &vertex->right->pcent, &vertex->right->nvcent );
  FindRBoundingBoxf ( tree, vertex->left );
  FindRBoundingBoxf ( tree, vertex->right );
  UpdateRBoundingBoxesf ( vertex );
  vertex->leaf = false;
} /*DivideRVertexf*/

RBezPatchTreeVertexfp
  rbez_GetRBezLeftVertexf ( RBezPatchTreefp tree,
                            RBezPatchTreeVertexfp vertex )
{
  if ( vertex->leaf ) {
    if ( raybez_use_mutex )
      pthread_mutex_lock ( &raybez_mutex );
    if ( !vertex->left )
      DivideRVertexf ( tree, vertex );
    if ( raybez_use_mutex )
      pthread_mutex_unlock ( &raybez_mutex );
  }
  return vertex->left;
} /*rbez_GetRBezLeftVertexf*/

RBezPatchTreeVertexfp
  rbez_GetRBezRightVertexf ( RBezPatchTreefp tree,
                             RBezPatchTreeVertexfp vertex )
{
  if ( vertex->leaf ) {
    if ( raybez_use_mutex )
      pthread_mutex_lock ( &raybez_mutex );
    if ( !vertex->right )
      DivideRVertexf ( tree, vertex );
    if ( raybez_use_mutex )
      pthread_mutex_unlock ( &raybez_mutex );
  }
  return vertex->right;
} /*rbez_GetRBezRightVertexf*/

/* ///////////////////////////////////////////////////////////////////////// */
static void ConvertPatchf ( int ncp, const point4f *ctlpoints,
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
} /*ConvertPatchf*/

static char SolutionOKf ( ray3f *ray, RBezPatchTreef *tree,
                         RBezPatchTreeVertexf *vertex, point2f *z,
                         int *ninters, RayObjectIntersf *inters )
{
  float   u, v;
  vector3f pv;

  if ( z->x >= 0.0 && z->x <= 1.0 && z->y >= 0.0 && z->y <= 1.0 ) {
    inters += *ninters;
    u = inters->u = vertex->u0 + z->x*(vertex->u1 - vertex->u0);
    u = (u - tree->root->u0)/(tree->root->u1 - tree->root->u0);
    v = inters->v = vertex->v0 + z->y*(vertex->v1 - vertex->v0);
    v = (v - tree->root->v0)/(tree->root->v1 - tree->root->v0);

    mbs_BCHornerNvP3Rf ( tree->n, tree->m, tree->root->ctlpoints,
                         u, v, &inters->p, &inters->nv );
    SubtractPoints3f ( &inters->p, &ray->p, &pv );
    inters->t = (float)DotProduct3f ( &pv, &ray->v );
    if ( inters->t > 0.0 ) {
      inters->object_id = tree->object_id;
      (*ninters)++;
    }
    return true;
  }
  else
    return false;
} /*SolutionOKf*/

static void OutputSingularSolutionf ( ray3f *ray, RBezPatchTreef *tree,
                                      RBezPatchTreeVertexf *vertex,
                                      int *ninters, RayObjectIntersf *inters )
{
  float    u, v;
  vector3f pv;

  inters += *ninters;
  inters->u = (float)(0.5*(vertex->u0 + vertex->u1));
  u = (inters->u - tree->root->u0)/(tree->root->u1 - tree->root->u0);
  inters->v = (float)(0.5*(vertex->v0 + vertex->v1));
  v = (inters->v - tree->root->v0)/(tree->root->v1 - tree->root->v0);
  mbs_BCHornerNvP3Rf ( tree->n, tree->m, tree->root->ctlpoints,
                       u, v, &inters->p, &inters->nv );
  SubtractPoints3f ( &inters->p, &ray->p, &pv );
  inters->t = (float)DotProduct3f ( &pv, &ray->v );
  if ( inters->t > 0.0 ) {
    inters->object_id = tree->object_id;
    (*ninters)++;
  }
} /*OutputSingularSolutionf*/

int rbez_FindRayRBezPatchIntersf ( RBezPatchTreef *tree, ray3f *ray,
                                   int maxlevel, int maxinters,
                                   int *ninters, RayObjectIntersf *inters )
{
  int      n, m, ncp;
  int      size_auxcp, size_stack;
  int      sp;           /* stack pointer */
  RBezPatchTreeVertexfp *stack;
  RBezPatchTreeVertexfp vertex, left, right;
  point2f  *auxcp;
  vector3f nh;           /* Householder reflection vector */
  float    sh;           /* scaling coefficient */
  float    K1, K2;       /* uniqueness test numbers */
  point2f  p, z;
  vector2f du, dv;

  n = tree->n;  m = tree->m;
  ncp = (n+1)*(m+1);
  stack = (RBezPatchTreeVertexfp*)pkv_GetScratchMem
              ( size_stack = maxlevel*sizeof(RBezPatchTreeVertexfp) );
                        /* the auxcp array is twice longer for the */
                        /* needs of solution uniqueness test */
  auxcp = (point2f*)pkv_GetScratchMem ( size_auxcp = 2*ncp*sizeof(point2f) );

  if ( !auxcp || !stack )
    goto failure;

                        /* construct the Householder reflection */
  nh = ray->v;
  if ( nh.z > 0 ) nh.z += 1.0;
    else          nh.z -= 1.0;
  sh = (float)(2.0/DotProduct3f ( &nh, &nh ));

                        /* initialize stack, push the patch and go around */
  *ninters = 0;
  sp = 0;
  stack[sp++] = tree->root;
  do {
    vertex = stack[--sp];
    if ( rbez_TestRayBBoxf ( ray, &vertex->bbox ) ) {
      ConvertPatchf ( ncp, vertex->ctlpoints, &ray->p, &nh, sh, auxcp );
      if ( _rbez_ConvexHullTest2f ( ncp, auxcp ) ) {
        if ( _rbez_UniquenessTest2f ( n, m, ncp, auxcp, &p, &du, &dv, &K1, &K2 ) ) {
          switch ( _rbez_NewtonMethod2f ( n, m, auxcp, &p, &du, &dv, &z ) ) {
        case RBEZ_NEWTON_YES:
            if ( !SolutionOKf ( ray, tree, vertex, &z, ninters, inters ) ) {
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
            stack[sp++] = left;
            stack[sp++] = right;
          }
          else
            OutputSingularSolutionf ( ray, tree, vertex, ninters, inters );
        }
      }
    }
  } while ( sp && *ninters < maxinters );

  pkv_FreeScratchMem ( size_auxcp+size_stack );
  return *ninters;

failure:
  pkv_FreeScratchMem ( size_auxcp+size_stack );
  return -1;
} /*FindRayRBezPatchIntersf*/

