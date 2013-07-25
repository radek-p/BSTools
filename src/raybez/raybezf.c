
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
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
static BezPatchTreeVertexfp
  AllocTreeVertexf ( BezPatchTreefp tree,
                     BezPatchTreeVertexfp up,
                     float u0, float u1, float v0, float v1 )
{
  BezPatchTreeVertexfp vertex;

  PKV_MALLOC ( vertex, sizeof(BezPatchTreeVertexf) + tree->cpsize );
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
    vertex->ctlpoints = (point3f*)((char*)vertex+sizeof(BezPatchTreeVertexf));
    vertex->normalvect = NULL;
  }
  return vertex;
} /*AllocTreeVertexf*/

static void FindBoundingBoxf ( BezPatchTreefp tree,
                                BezPatchTreeVertexfp vertex )
{
#define EPS 5.0e-6
  int      n, m, ncp;
  point3f  *cp, *cq, a, b;
  vector3f v;
  int      i, j;
  float    du, dv, d;

  n = tree->n;  m = tree->m;
  ncp = (n+1)*(m+1);
                                     /* find the bounding box */
  cp = vertex->ctlpoints;
  a = *cp;  cp ++;
  if ( ncp & 1 ) {
    vertex->bbox.x0 = vertex->bbox.x1 = a.x;
    vertex->bbox.y0 = vertex->bbox.y1 = a.y;
    vertex->bbox.z0 = vertex->bbox.z1 = a.z;
    i = 2;
  }
  else {
    b = *cp;  cp ++;
    pkv_Sort2f ( &a.x, &b.x );  vertex->bbox.x0 = a.x;  vertex->bbox.x1 = b.x;
    pkv_Sort2f ( &a.y, &b.y );  vertex->bbox.y0 = a.y;  vertex->bbox.y1 = b.y;
    pkv_Sort2f ( &a.z, &b.z );  vertex->bbox.z0 = a.z;  vertex->bbox.z1 = b.z;
    i = 3;
  }
  for ( ; i < ncp; i += 2 ) {
    a = *cp;  cp++;
    b = *cp;  cp++;
    pkv_Sort2f ( &a.x, &b.x );
    vertex->bbox.x0 = min ( vertex->bbox.x0, a.x );
    vertex->bbox.x1 = max ( vertex->bbox.x1, b.x );
    pkv_Sort2f ( &a.y, &b.y );
    vertex->bbox.y0 = min ( vertex->bbox.y0, a.y );
    vertex->bbox.y1 = max ( vertex->bbox.y1, b.y );
    pkv_Sort2f ( &a.z, &b.z );
    vertex->bbox.z0 = min ( vertex->bbox.z0, a.z );
    vertex->bbox.z1 = max ( vertex->bbox.z1, b.z );
  }
  vertex->bbox.x0 -= EPS;  vertex->bbox.x1 += EPS;
  vertex->bbox.y0 -= EPS;  vertex->bbox.y1 += EPS;
  vertex->bbox.z0 -= EPS;  vertex->bbox.z1 += EPS;
                                     /* determine division direction */
  du = 0.0;
  for ( i = 0, cp = vertex->ctlpoints, cq = cp+(m+1);
        i < n*(m+1);
        i++, cp++, cq++ ) {
    SubtractPoints3f ( cq, cp, &v );
    d = (float)DotProduct3f ( &v, &v );
    du = max ( du, d );
  }
  du = (float)(n * sqrt(du));
  dv = 0.0;
  for ( i = 0, cp = vertex->ctlpoints;
        i <= n;
        i++, cp++ )
    for ( j = 0; j < m; j++, cp++ ) {
      SubtractPoints3f ( cp+1, cp, &v );
      d = (float)DotProduct3f ( &v, &v );
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
} /*FindBoundingBoxf*/

static void UpdateBoundingBoxesf ( BezPatchTreeVertexfp vertex )
{
  float a;
  char  change;

  while ( vertex ) {
    change = 0;
    a = min ( vertex->left->bbox.x0, vertex->right->bbox.x0 );
    if ( vertex->bbox.x0 < a ) { vertex->bbox.x0 = a;  change = 1; }
    a = max ( vertex->left->bbox.x1, vertex->right->bbox.x1 );
    if ( vertex->bbox.x1 > a ) { vertex->bbox.x1 = a;  change = 1; }
    a = min ( vertex->left->bbox.y0, vertex->right->bbox.y0 );
    if ( vertex->bbox.y0 < a ) { vertex->bbox.y0 = a;  change = 1; }
    a = max ( vertex->left->bbox.y1, vertex->right->bbox.y1 );
    if ( vertex->bbox.y1 > a ) { vertex->bbox.y1 = a;  change = 1; }
    a = min ( vertex->left->bbox.z0, vertex->right->bbox.z0 );
    if ( vertex->bbox.z0 < a ) { vertex->bbox.z0 = a;  change = 1; }
    a = max ( vertex->left->bbox.z1, vertex->right->bbox.z1 );
    if ( vertex->bbox.z1 > a ) { vertex->bbox.z1 = a;  change = 1; }

    if ( change )
      vertex = vertex->up;
    else
      return;
  }
} /*UpdateBoundingBoxesf*/

BezPatchTreefp
  rbez_NewBezPatchTreef ( int object_id,
                          unsigned char n, unsigned char m,
                          float u0, float u1, float v0, float v1,
                          CONST_ point3f *ctlpoints )
{
  BezPatchTreefp       tree;

  PKV_MALLOC ( tree, sizeof(BezPatchTreef) );
  if ( tree ) {
    tree->object_id = object_id;
    tree->n = n;    tree->m = m;
    tree->cpsize = (n+1)*(m+1)*sizeof(point3f);
    tree->root = AllocTreeVertexf ( tree, NULL, u0, u1, v0, v1 );
    if ( tree->root ) {
      memcpy ( tree->root->ctlpoints, ctlpoints, tree->cpsize );
      mbs_BCHornerNvP3f ( n, m, ctlpoints, 0.5, 0.5,
                          &tree->root->pcent, &tree->root->nvcent );
      FindBoundingBoxf ( tree, tree->root );
    }
    else {
      PKV_FREE ( tree );
      return NULL;
    }
  }
  return tree;
} /*rbez_NewBezPatchTreef*/

static void r_DestroyTreef ( BezPatchTreeVertexfp vertex )
{
  if ( vertex->right )
    r_DestroyTreef ( vertex->right );
  if ( vertex->left )
    r_DestroyTreef ( vertex->left );
  if ( vertex->normalvect )
    PKV_FREE ( vertex->normalvect );
  PKV_FREE ( vertex );
} /*r_DestroyTreef*/

void rbez_DestroyBezPatchTreef ( BezPatchTreefp tree )
{
  r_DestroyTreef ( tree->root );
  PKV_FREE ( tree );
} /*rbez_DestroyBezPatchTreef*/

static void DivideVertexf ( BezPatchTreefp tree, BezPatchTreeVertexfp vertex )
{
  float h;

        /* allocate subtree roots */
  if ( !vertex->divdir ) {
    h = (float)(0.5*(vertex->u0 + vertex->u1));
    vertex->left  = AllocTreeVertexf ( tree, vertex,
                                        vertex->u0, h, vertex->v0, vertex->v1 );
    vertex->right = AllocTreeVertexf ( tree, vertex,
                                        h, vertex->u1, vertex->v0, vertex->v1 );
    if ( !vertex->left || !vertex->right )
      goto dealloc;
    memcpy ( vertex->right->ctlpoints, vertex->ctlpoints, tree->cpsize );
    mbs_BisectBP3uf ( tree->n, tree->m,
                      vertex->right->ctlpoints, vertex->left->ctlpoints );
  }
  else {
    h = (float)(0.5*(vertex->v0 + vertex->v1));
    vertex->left  = AllocTreeVertexf ( tree, vertex,
                                        vertex->u0, vertex->u1, vertex->v0, h );
    vertex->right = AllocTreeVertexf ( tree, vertex,
                                        vertex->u0, vertex->u1, h, vertex->v1 );
    if ( !vertex->left || !vertex->right ) {
dealloc:
      if ( vertex->left )  PKV_FREE ( vertex->left );
      if ( vertex->right ) PKV_FREE ( vertex->right );
      return;
    }
    memcpy ( vertex->right->ctlpoints, vertex->ctlpoints, tree->cpsize );
    mbs_BisectBP3vf ( tree->n, tree->m,
                      vertex->right->ctlpoints, vertex->left->ctlpoints );
  }
  mbs_BCHornerNvP3f ( tree->n, tree->m, vertex->left->ctlpoints,
                      0.5, 0.5, &vertex->left->pcent, &vertex->left->nvcent );
  mbs_BCHornerNvP3f ( tree->n, tree->m, vertex->right->ctlpoints,
                      0.5, 0.5, &vertex->right->pcent, &vertex->right->nvcent );
  FindBoundingBoxf ( tree, vertex->left );
  FindBoundingBoxf ( tree, vertex->right );
  UpdateBoundingBoxesf ( vertex );
  vertex->leaf = false;
} /*DivideVertexf*/

BezPatchTreeVertexfp
  rbez_GetBezLeftVertexf ( BezPatchTreefp tree,
                            BezPatchTreeVertexfp vertex )
{
  if ( vertex->leaf ) {
    if ( raybez_use_mutex )
      pthread_mutex_lock ( &raybez_mutex );
    if ( !vertex->left )
      DivideVertexf ( tree, vertex );
    if ( raybez_use_mutex )
      pthread_mutex_unlock ( &raybez_mutex );
  }
  return vertex->left;
} /*rbez_GetBezLeftVertexf*/

BezPatchTreeVertexfp
  rbez_GetBezRightVertexf ( BezPatchTreefp tree,
                             BezPatchTreeVertexfp vertex )
{
  if ( vertex->leaf ) {
    if ( raybez_use_mutex )
      pthread_mutex_lock ( &raybez_mutex );
    if ( !vertex->right )
      DivideVertexf ( tree, vertex );
    if ( raybez_use_mutex )
      pthread_mutex_unlock ( &raybez_mutex );
  }
  return vertex->right;
} /*rbez_GetBezRightVertexf*/

/* ///////////////////////////////////////////////////////////////////////// */
static void ConvertPatchf ( int ncp, const point3f *ctlpoints,
                            const point3f *rayp, const vector3f *nh, float sh,
                            point2f *auxcp )
{
  int      i;
  vector3f v;
  float    r;

  for ( i = 0; i < ncp; i++ ) {
    SubtractPoints3f ( &ctlpoints[i], rayp, &v );
    r = (float)(sh * DotProduct3f ( &v, nh ));
    auxcp[i].x = v.x - r*nh->x;
    auxcp[i].y = v.y - r*nh->y;
  }
} /*ConvertPatchf*/

static char SolutionOKf ( ray3f *ray, BezPatchTreef *tree,
                          BezPatchTreeVertexf *vertex, point2f *z,
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

    mbs_BCHornerNvP3f ( tree->n, tree->m, tree->root->ctlpoints,
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

static void OutputSingularSolutionf ( ray3f *ray, BezPatchTreef *tree,
                                      BezPatchTreeVertexf *vertex,
                                      int *ninters, RayObjectIntersf *inters )
{
  float    u, v;
  vector3f pv;

  inters += *ninters;
  inters->u = (float)(0.5*(vertex->u0 + vertex->u1));
  u = (inters->u - tree->root->u0)/(tree->root->u1 - tree->root->u0);
  inters->v = (float)(0.5*(vertex->v0 + vertex->v1));
  v = (inters->v - tree->root->v0)/(tree->root->v1 - tree->root->v0);
  mbs_BCHornerNvP3f ( tree->n, tree->m, tree->root->ctlpoints,
                      u, v, &inters->p, &inters->nv );
  SubtractPoints3f ( &inters->p, &ray->p, &pv );
  inters->t = (float)DotProduct3f ( &pv, &ray->v );
  if ( inters->t > 0.0 ) {
    inters->object_id = tree->object_id;
    (*ninters)++;
  }
} /*OutputSingularSolutionf*/

int rbez_FindRayBezPatchIntersf ( BezPatchTreef *tree, ray3f *ray,
                                   int maxlevel, int maxinters,
                                   int *ninters, RayObjectIntersf *inters )
{
  int      n, m, ncp;
  int      size_auxcp, size_stack;
  int      sp;           /* stack pointer */
  BezPatchTreeVertexfp *stack;
  BezPatchTreeVertexfp vertex, left, right;
  point2f  *auxcp;
  vector3f nh;           /* Householder reflection vector */
  float    sh;           /* scaling coefficient */
  float    K1, K2;       /* uniqueness test numbers */
  point2f  p, z;
  vector2f du, dv;

  n = tree->n;  m = tree->m;
  ncp = (n+1)*(m+1);
  stack = (BezPatchTreeVertexfp*)pkv_GetScratchMem
              ( size_stack = maxlevel*sizeof(BezPatchTreeVertexfp) );
                        /* the auxcp array is twice longer for the */
                        /* needs of solution uniqueness test */
  auxcp = (point2f*)pkv_GetScratchMem ( size_auxcp = 2*ncp*sizeof(point2f) );

  if ( !auxcp || !stack )
    exit ( 0 );

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
          if ( _rbez_NewtonMethod2f ( n, m, auxcp, &p, &du, &dv, &z ) ) {
            if ( !SolutionOKf ( ray, tree, vertex, &z, ninters, inters ) )
              if ( _rbez_SecondTest2f ( &z, n, m, K1, K2 ) )
                goto DIVIDE;
          }
          else goto DIVIDE;
        }
        else {
DIVIDE:
          if ( vertex->level < maxlevel ) {
                        /* divide and push pieces */
            left  = rbez_GetBezLeftVertexf  ( tree, vertex );
            right = rbez_GetBezRightVertexf ( tree, vertex );
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
} /*rbez_FindRayBezPatchIntersf*/

