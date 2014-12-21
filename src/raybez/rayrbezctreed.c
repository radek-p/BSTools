
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012, 2014                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

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
static RBezCurveTreeVertexdp
    AllocRCTreeVertexd ( RBezCurveTreedp tree,
                         RBezCurveTreeVertexdp up,
                         double t0, double t1 )
{
  RBezCurveTreeVertexdp vertex;

  PKV_MALLOC ( vertex, sizeof(RBezCurveTreeVertexd) + tree->cpsize );
  if ( vertex ) {
    vertex->t0 = t0;
    vertex->t1 = t1;
    vertex->left = vertex->right = NULL;
    vertex->up = up;
    if ( up )
      vertex->level = (short)(up->level + 1);
    else
      vertex->level = 0;
    vertex->ctlpoints = (point4d*)((char*)vertex+sizeof(RBezCurveTreeVertexd));
  }
  return vertex;
} /*AllocRCTreeVertexd*/

static void FindRCBoundingBoxd ( RBezCurveTreedp tree,
                                 RBezCurveTreeVertexdp vertex )
{
#define EPS 1.0e-10
  int ncp;
  point4d *cp, dcp;
  point3d a, b;
  int     i;
  double   ext, d, e;

  ncp = tree->degree + 1;
  cp = vertex->ctlpoints;
  Point4to3d ( cp++, &a );
  if ( ncp & 1 ) {
    vertex->bbox.x0 = vertex->bbox.x1 = a.x;
    vertex->bbox.y0 = vertex->bbox.y1 = a.y;
    vertex->bbox.z0 = vertex->bbox.z1 = a.z;
    i = 2;
  }
  else {
    Point4to3d ( cp++, &b );
    pkv_Sort2d ( &a.x, &b.x );  vertex->bbox.x0 = a.x;  vertex->bbox.x1 = b.x;
    pkv_Sort2d ( &a.y, &b.y );  vertex->bbox.y0 = a.y;  vertex->bbox.y1 = b.y;
    pkv_Sort2d ( &a.z, &b.z );  vertex->bbox.z0 = a.z;  vertex->bbox.z1 = b.z;
    i = 3;
  }
  for ( ; i < ncp; i += 2 ) {
    Point4to3d ( cp++, &a );
    Point4to3d ( cp++, &b );
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
    SubtractPoints4d ( &cp[i], &cp[i-1], &dcp );
    e = DotProduct4d ( &dcp, &dcp );
    if ( e > d )
      d = e;
  }
  vertex->maxder = (double)tree->degree*sqrt ( d );
#undef EPS
} /*FindRCBoundingBoxd*/

static void UpdateRCBoundingBoxesd ( RBezCurveTreeVertexdp vertex )
{
  while ( vertex ) {
    if ( rbez_NarrowBBoxSumd ( &vertex->left->bbox, &vertex->right->bbox,
                               &vertex->bbox ) )
      vertex = vertex->up;
    else
      return;
  }
} /*UpdateRCBoundingBoxesd*/

RBezCurveTreedp rbez_NewRBezCurveTreed ( int object_id, short degree,
                                        double t0, double t1, double ext,
                                        CONST_ point4d *ctlpoints )
{
  RBezCurveTreedp tree;

  PKV_MALLOC ( tree, sizeof(RBezCurveTreed) );
  if ( tree ) {
    tree->object_id = object_id;
    tree->degree = degree;
    tree->cpsize = (degree+1)*sizeof(point4d);
    tree->ext = ext;
    tree->root = AllocRCTreeVertexd ( tree, NULL, t0, t1 );
    if ( tree->root ) {
      memcpy ( tree->root->ctlpoints, ctlpoints, tree->cpsize );
      FindRCBoundingBoxd ( tree, tree->root );
    }
    else {
      PKV_FREE ( tree );
      return NULL;
    }
  }
  return tree;
} /*rbez_NewRBezCurveTreed*/

static void r_DestroyRCTreed ( RBezCurveTreeVertexdp vertex )
{
  if ( vertex->right )
    r_DestroyRCTreed ( vertex->right );
  if ( vertex->left )
    r_DestroyRCTreed ( vertex->left );
  PKV_FREE ( vertex );
} /*r_DestroyRCTreed*/

void rbez_DestroyRBezCurveTreed ( RBezCurveTreedp tree )
{
  r_DestroyRCTreed ( tree->root );
  PKV_FREE ( tree );
} /*rbez_DestroyRBezCurveTreed*/

static void DivideRCVertexd ( RBezCurveTreedp tree, RBezCurveTreeVertexdp vertex )
{
  double h;

  h = 0.5*(vertex->t0 + vertex->t1);
  vertex->left  = AllocRCTreeVertexd ( tree, vertex, vertex->t0, h );
  vertex->right = AllocRCTreeVertexd ( tree, vertex, h, vertex->t1 );
  if ( !vertex->left || !vertex->right )
    return;
  memcpy ( vertex->right->ctlpoints, vertex->ctlpoints, tree->cpsize );
  mbs_BisectBC4d ( tree->degree,
                   vertex->right->ctlpoints, vertex->left->ctlpoints );
  FindRCBoundingBoxd ( tree, vertex->left );
  FindRCBoundingBoxd ( tree, vertex->right );
  UpdateRCBoundingBoxesd ( vertex );
} /*DivideRCVertexd*/

RBezCurveTreeVertexdp rbez_GetRBezCurveLeftVertexd ( RBezCurveTreedp tree,
                                                     RBezCurveTreeVertexdp vertex )
{
  if ( !vertex->left ) {
    if ( raybez_use_mutex )
      pthread_mutex_lock ( &raybez_mutex );
    if ( !vertex->left )
      DivideRCVertexd ( tree, vertex );
    if ( raybez_use_mutex )
      pthread_mutex_unlock ( &raybez_mutex );
  }
  return vertex->left;
} /*rbez_GetRBezCurveLeftVertexd*/

RBezCurveTreeVertexdp rbez_GetRBezCurveRightVertexd ( RBezCurveTreedp tree,
                                                      RBezCurveTreeVertexdp vertex )
{
  if ( !vertex->right ) {
    if ( raybez_use_mutex )
      pthread_mutex_lock ( &raybez_mutex );
    if ( !vertex->right )
      DivideRCVertexd ( tree, vertex );
    if ( raybez_use_mutex )
      pthread_mutex_unlock ( &raybez_mutex );
  }
  return vertex->right;
} /*rbez_GetRBezCurveRightVertexd*/

