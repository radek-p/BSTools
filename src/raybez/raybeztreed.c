
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012, 2014                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* Changes: */
/*  8.07.2012, P. Kiciak - added locking/unlocking a mutex for pthreads */
/* 31.08.2012, P. Kiciak - added the leaf attribute, to avoid */
/*                         unnecessary mutex locking */
/* 16.01.2013, P. Kiciak - separated the tree processing code from */
/*                         computing intersections of a ray with the tube */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>

#define CONST_

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "raybez.h"

/* ////////////////////////////////////////////////////////////////////////// */
static BezCurveTreeVertexdp
    AllocCTreeVertexd ( BezCurveTreedp tree,
                        BezCurveTreeVertexdp up,
                        double t0, double t1 )
{
  BezCurveTreeVertexdp vertex;

  PKV_MALLOC ( vertex, sizeof(BezCurveTreeVertexd) + tree->cpsize );
  if ( vertex ) {
    vertex->t0 = t0;
    vertex->t1 = t1;
    vertex->left = vertex->right = NULL;
    vertex->up = up;
    vertex->leaf = true;
    if ( up )
      vertex->level = (short)(up->level + 1);
    else
      vertex->level = 0;
    vertex->ctlpoints = (point3d*)((char*)vertex+sizeof(BezCurveTreeVertexd));
  }
  return vertex;
} /*AllocCTreeVertexd*/

static void FindCBoundingBoxd ( BezCurveTreedp tree,
                                BezCurveTreeVertexdp vertex )
{
#define EPS 1.0e-10
  int ncp;
  point3d *cp, a, b;
  int     i;
  double  ext, d, e;

  ncp = tree->degree + 1;
  cp = vertex->ctlpoints;
  a = *cp; cp ++;
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
    a = *cp;  cp ++;
    b = *cp;  cp ++;
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
  ext = tree->ext+EPS;
  vertex->bbox.x0 -= ext;  vertex->bbox.x1 += ext;
  vertex->bbox.y0 -= ext;  vertex->bbox.y1 += ext;
  vertex->bbox.z0 -= ext;  vertex->bbox.z1 += ext;
  cp = vertex->ctlpoints;
  d = 0.0;
  for ( i = 1; i < ncp; i++ ) {
    SubtractPoints3d ( &cp[i], &cp[i-1], &a );
    e = DotProduct3d ( &a, &a );
    if ( e > d )
      d = e;
  }
  vertex->maxder = (double)tree->degree*sqrt ( d );
#undef EPS
} /*FindCBoundingBoxd*/

static void UpdateCBoundingBoxesd ( BezCurveTreeVertexdp vertex )
{
  while ( vertex ) {
    if ( rbez_NarrowBBoxSumd ( &vertex->left->bbox, &vertex->right->bbox,
                               &vertex->bbox ) )
      vertex = vertex->up;
    else
      return;
  }
} /*UpdateCBoundingBoxesd*/

BezCurveTreedp rbez_NewBezCurveTreed ( int object_id, short degree,
                                       double t0, double t1, double ext,
                                       CONST_ point3d *ctlpoints )
{
  BezCurveTreedp tree;

  PKV_MALLOC ( tree, sizeof(BezCurveTreed) );
  if ( tree ) {
    tree->object_id = object_id;
    tree->degree = degree;
    tree->cpsize = (degree+1)*sizeof(point3d);
    tree->ext = ext;
    tree->root = AllocCTreeVertexd ( tree, NULL, t0, t1 );
    if ( tree->root ) {
      memcpy ( tree->root->ctlpoints, ctlpoints, tree->cpsize );
      FindCBoundingBoxd ( tree, tree->root );
    }
    else {
      PKV_FREE ( tree );
      return NULL;
    }
  }
  return tree;
} /*rbez_NewBezCurveTreed*/

static void r_DestroyCTreed ( BezCurveTreeVertexdp vertex )
{
  if ( vertex->right )
    r_DestroyCTreed ( vertex->right );
  if ( vertex->left )
    r_DestroyCTreed ( vertex->left );
  PKV_FREE ( vertex );
} /*r_DestroyCTreed*/

void rbez_DestroyBezCurveTreed ( BezCurveTreedp tree )
{
  r_DestroyCTreed ( tree->root );
  PKV_FREE ( tree );
} /*rbez_DestroyBezCurveTreed*/

static void DivideCVertexd ( BezCurveTreedp tree, BezCurveTreeVertexdp vertex )
{
  double h;

  h = 0.5*(vertex->t0 + vertex->t1);
  vertex->left  = AllocCTreeVertexd ( tree, vertex, vertex->t0, h );
  vertex->right = AllocCTreeVertexd ( tree, vertex, h, vertex->t1 );
  if ( !vertex->left || !vertex->right )
    return;
  memcpy ( vertex->right->ctlpoints, vertex->ctlpoints, tree->cpsize );
  mbs_BisectBC3d ( tree->degree,
                   vertex->right->ctlpoints, vertex->left->ctlpoints );
  FindCBoundingBoxd ( tree, vertex->left );
  FindCBoundingBoxd ( tree, vertex->right );
  UpdateCBoundingBoxesd ( vertex );
  vertex->leaf = false;
} /*DivideCVertexd*/

BezCurveTreeVertexdp rbez_GetBezCurveLeftVertexd ( BezCurveTreedp tree,
                                                   BezCurveTreeVertexdp vertex )
{
  if ( vertex->leaf ) {
    if ( raybez_use_mutex )
      pthread_mutex_lock ( &raybez_mutex );
    if ( !vertex->left )
      DivideCVertexd ( tree, vertex );
    if ( raybez_use_mutex )
      pthread_mutex_unlock ( &raybez_mutex );
  }
  return vertex->left;
} /*rbez_GetBezCurveLeftVertexd*/

BezCurveTreeVertexdp rbez_GetBezCurveRightVertexd ( BezCurveTreedp tree,
                                                    BezCurveTreeVertexdp vertex )
{
  if ( vertex->leaf ) {
    if ( raybez_use_mutex )
      pthread_mutex_lock ( &raybez_mutex );
    if ( !vertex->right )
      DivideCVertexd ( tree, vertex );
    if ( raybez_use_mutex )
      pthread_mutex_unlock ( &raybez_mutex );
  }
  return vertex->right;
} /*rbez_GetBezCurveRightVertexd*/

