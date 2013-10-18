
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
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
static BezPatchTreeVertexdp
  AllocTreeVertexd ( BezPatchTreedp tree,
                     BezPatchTreeVertexdp up,
                     double u0, double u1, double v0, double v1 )
{
  BezPatchTreeVertexdp vertex;

  PKV_MALLOC ( vertex, sizeof(BezPatchTreeVertexd) + tree->cpsize );
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
    vertex->ctlpoints = (point3d*)((char*)vertex+sizeof(BezPatchTreeVertexd));
    vertex->normalvect = NULL;
  }
  return vertex;
} /*AllocTreeVertexd*/

static void FindBoundingBoxd ( BezPatchTreedp tree,
                                BezPatchTreeVertexdp vertex )
{
#define EPS 1.0e-10
  int      n, m, ncp;
  point3d  *cp, *cq, a, b;
  vector3d v;
  int      i, j;
  double    du, dv, d;

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
    b = *cp;  cp++;
    pkv_Sort2d ( &a.x, &b.x );  vertex->bbox.x0 = a.x;  vertex->bbox.x1 = b.x;
    pkv_Sort2d ( &a.y, &b.y );  vertex->bbox.y0 = a.y;  vertex->bbox.y1 = b.y;
    pkv_Sort2d ( &a.z, &b.z );  vertex->bbox.z0 = a.z;  vertex->bbox.z1 = b.z;
    i = 3;
  }
  for ( ; i < ncp; i += 2 ) {
    a = *cp;  cp++;
    b = *cp;  cp++;
    pkv_Sort2d ( &a.x, &b.x );
    vertex->bbox.x0 = min ( vertex->bbox.x0, a.x );
    vertex->bbox.x1 = max ( vertex->bbox.x1, b.x );
    pkv_Sort2d ( &a.y, &b.y );
    vertex->bbox.y0 = min ( vertex->bbox.y0, a.y );
    vertex->bbox.y1 = max ( vertex->bbox.y1, b.y );
    pkv_Sort2d ( &a.z, &b.z );
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
    SubtractPoints3d ( cq, cp, &v );
    d = DotProduct3d ( &v, &v );
    du = max ( du, d );
  }
  du = (double)n * sqrt(du);
  dv = 0.0;
  for ( i = 0, cp = vertex->ctlpoints;
        i <= n;
        i++, cp++ )
    for ( j = 0; j < m; j++, cp++ ) {
      SubtractPoints3d ( cp+1, cp, &v );
      d = DotProduct3d ( &v, &v );
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
} /*FindBoundingBoxd*/

static void UpdateBoundingBoxesd ( BezPatchTreeVertexdp vertex )
{
  double a;
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
} /*UpdateBoundingBoxesd*/

BezPatchTreedp
  rbez_NewBezPatchTreed ( int object_id,
                          unsigned char n, unsigned char m,
                          double u0, double u1, double v0, double v1,
                          CONST_ point3d *ctlpoints )
{
  BezPatchTreedp       tree;

  PKV_MALLOC ( tree, sizeof(BezPatchTreed) );
  if ( tree ) {
    tree->object_id = object_id;
    tree->n = n;    tree->m = m;
    tree->cpsize = (n+1)*(m+1)*sizeof(point3d);
    tree->root = AllocTreeVertexd ( tree, NULL, u0, u1, v0, v1 );
    if ( tree->root ) {
      memcpy ( tree->root->ctlpoints, ctlpoints, tree->cpsize );
      if ( !mbs_BCHornerNvP3d ( n, m, ctlpoints, 0.5, 0.5,
                                &tree->root->pcent, &tree->root->nvcent ) )
        goto failure;
      FindBoundingBoxd ( tree, tree->root );
    }
    else {
failure:
      PKV_FREE ( tree );
      return NULL;
    }
  }
  return tree;
} /*rbez_NewBezPatchTreed*/

static void r_DestroyTreed ( BezPatchTreeVertexdp vertex )
{
  if ( vertex->right )
    r_DestroyTreed ( vertex->right );
  if ( vertex->left )
    r_DestroyTreed ( vertex->left );
  if ( vertex->normalvect )
    PKV_FREE ( vertex->normalvect );
  PKV_FREE ( vertex );
} /*r_DestroyTreed*/

void rbez_DestroyBezPatchTreed ( BezPatchTreedp tree )
{
  r_DestroyTreed ( tree->root );
  PKV_FREE ( tree );
} /*rbez_DestroyBezPatchTreed*/

static void DivideVertexd ( BezPatchTreedp tree, BezPatchTreeVertexdp vertex )
{
  double h;

        /* allocate subtree roots */
  if ( !vertex->divdir ) {
    h = 0.5*(vertex->u0 + vertex->u1);
    vertex->left  = AllocTreeVertexd ( tree, vertex,
                                        vertex->u0, h, vertex->v0, vertex->v1 );
    vertex->right = AllocTreeVertexd ( tree, vertex,
                                        h, vertex->u1, vertex->v0, vertex->v1 );
    if ( !vertex->left || !vertex->right )
      goto dealloc;
    memcpy ( vertex->right->ctlpoints, vertex->ctlpoints, tree->cpsize );
    mbs_BisectBP3ud ( tree->n, tree->m,
                      vertex->right->ctlpoints, vertex->left->ctlpoints );
  }
  else {
    h = 0.5*(vertex->v0 + vertex->v1);
    vertex->left  = AllocTreeVertexd ( tree, vertex,
                                        vertex->u0, vertex->u1, vertex->v0, h );
    vertex->right = AllocTreeVertexd ( tree, vertex,
                                        vertex->u0, vertex->u1, h, vertex->v1 );
    if ( !vertex->left || !vertex->right )
      goto dealloc;
    memcpy ( vertex->right->ctlpoints, vertex->ctlpoints, tree->cpsize );
    mbs_BisectBP3vd ( tree->n, tree->m,
                      vertex->right->ctlpoints, vertex->left->ctlpoints );
  }
  if ( !mbs_BCHornerNvP3d ( tree->n, tree->m, vertex->left->ctlpoints,
                      0.5, 0.5, &vertex->left->pcent, &vertex->left->nvcent ) )
    goto dealloc;
  if ( !mbs_BCHornerNvP3d ( tree->n, tree->m, vertex->right->ctlpoints,
                      0.5, 0.5, &vertex->right->pcent, &vertex->right->nvcent ) )
    goto dealloc;
  FindBoundingBoxd ( tree, vertex->left );
  FindBoundingBoxd ( tree, vertex->right );
  UpdateBoundingBoxesd ( vertex );
  vertex->leaf = false;
  return;

dealloc:
  if ( vertex->left )  PKV_FREE ( vertex->left );
  if ( vertex->right ) PKV_FREE ( vertex->right );
} /*DivideVertexd*/

BezPatchTreeVertexdp
  rbez_GetBezLeftVertexd ( BezPatchTreedp tree,
                            BezPatchTreeVertexdp vertex )
{
  if ( vertex->leaf ) {
    if ( raybez_use_mutex )
      pthread_mutex_lock ( &raybez_mutex );
    if ( !vertex->left )
      DivideVertexd ( tree, vertex );
    if ( raybez_use_mutex )
      pthread_mutex_unlock ( &raybez_mutex );
  }
  return vertex->left;
} /*rbez_GetBezLeftVertexd*/

BezPatchTreeVertexdp
  rbez_GetBezRightVertexd ( BezPatchTreedp tree,
                         BezPatchTreeVertexdp vertex )
{
  if ( vertex->leaf ) {
    if ( raybez_use_mutex )
      pthread_mutex_lock ( &raybez_mutex );
    if ( !vertex->right )
      DivideVertexd ( tree, vertex );
    if ( raybez_use_mutex )
      pthread_mutex_unlock ( &raybez_mutex );
  }
  return vertex->right;
} /*rbez_GetBezRightVertexd*/

/* ///////////////////////////////////////////////////////////////////////// */
static void ConvertPatchd ( int ncp, const point3d *ctlpoints,
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
} /*ConvertPatchd*/

static boolean SolutionOKd ( ray3d *ray, BezPatchTreed *tree,
                             BezPatchTreeVertexd *vertex, point2d *z,
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

    mbs_BCHornerNvP3d ( tree->n, tree->m, tree->root->ctlpoints,
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

static void OutputSingularSolutiond ( ray3d *ray, BezPatchTreed *tree,
                                      BezPatchTreeVertexd *vertex,
                                      int *ninters, RayObjectIntersd *inters )
{
  double    u, v;
  vector3d pv;

  inters += *ninters;
  inters->u = 0.5*(vertex->u0 + vertex->u1);
  u = (inters->u - tree->root->u0)/(tree->root->u1 - tree->root->u0);
  inters->v = 0.5*(vertex->v0 + vertex->v1);
  v = (inters->v - tree->root->v0)/(tree->root->v1 - tree->root->v0);
  mbs_BCHornerNvP3d ( tree->n, tree->m, tree->root->ctlpoints,
                      u, v, &inters->p, &inters->nv );
  SubtractPoints3d ( &inters->p, &ray->p, &pv );
  inters->t = DotProduct3d ( &pv, &ray->v );
  if ( inters->t > 0.0 ) {
    inters->object_id = tree->object_id;
    (*ninters)++;
  }
} /*OutputSingularSolutiond*/

int rbez_FindRayBezPatchIntersd ( BezPatchTreed *tree, ray3d *ray,
                                   int maxlevel, int maxinters,
                                   int *ninters, RayObjectIntersd *inters )
{
  int      n, m, ncp;
  int      size_auxcp, size_stack;
  int      sp;           /* stack pointer */
  BezPatchTreeVertexdp *stack;
  BezPatchTreeVertexdp vertex, left, right;
  point2d  *auxcp;
  vector3d nh;           /* Householder reflection vector */
  double    sh;           /* scaling coefficient */
  double    K1, K2;       /* uniqueness test numbers */
  point2d  p, z;
  vector2d du, dv;

  n = tree->n;  m = tree->m;
  ncp = (n+1)*(m+1);
  stack = (BezPatchTreeVertexdp*)pkv_GetScratchMem
              ( size_stack = maxlevel*sizeof(BezPatchTreeVertexdp) );
                        /* the auxcp array is twice longer for the */
                        /* needs of solution uniqueness test */
  auxcp = (point2d*)pkv_GetScratchMem ( size_auxcp = 2*ncp*sizeof(point2d) );

  if ( !auxcp || !stack )
    exit ( 0 );

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
            left  = rbez_GetBezLeftVertexd  ( tree, vertex );
            right = rbez_GetBezRightVertexd ( tree, vertex );
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
} /*rbez_FindRayBezPatchIntersd*/

