
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
static RBezCurveTreeVertexfp
    AllocRCTreeVertexf ( RBezCurveTreefp tree,
                         RBezCurveTreeVertexfp up,
                         float t0, float t1 )
{
  RBezCurveTreeVertexfp vertex;

  PKV_MALLOC ( vertex, sizeof(RBezCurveTreeVertexf) + tree->cpsize );
  if ( vertex ) {
    vertex->t0 = t0;
    vertex->t1 = t1;
    vertex->left = vertex->right = NULL;
    vertex->up = up;
    if ( up )
      vertex->level = (short)(up->level + 1);
    else
      vertex->level = 0;
    vertex->ctlpoints = (point4f*)((char*)vertex+sizeof(RBezCurveTreeVertexf));
  }
  return vertex;
} /*AllocRCTreeVertexf*/

static void FindRCBoundingBoxf ( RBezCurveTreefp tree,
                                 RBezCurveTreeVertexfp vertex )
{
#define EPS 1.0e-10
  int ncp;
  point4f *cp, dcp;
  point3f a, b;
  int     i;
  float   ext, d, e;

  ncp = tree->degree + 1;
  cp = vertex->ctlpoints;
  Point4to3f ( cp++, &a );
  if ( ncp & 1 ) {
    vertex->bbox.x0 = vertex->bbox.x1 = a.x;
    vertex->bbox.y0 = vertex->bbox.y1 = a.y;
    vertex->bbox.z0 = vertex->bbox.z1 = a.z;
    i = 2;
  }
  else {
    Point4to3f ( cp++, &b );
    pkv_Sort2f ( &a.x, &b.x );  vertex->bbox.x0 = a.x;  vertex->bbox.x1 = b.x;
    pkv_Sort2f ( &a.y, &b.y );  vertex->bbox.y0 = a.y;  vertex->bbox.y1 = b.y;
    pkv_Sort2f ( &a.z, &b.z );  vertex->bbox.z0 = a.z;  vertex->bbox.z1 = b.z;
    i = 3;
  }
  for ( ; i < ncp; i += 2 ) {
    Point4to3f ( cp++, &a );
    Point4to3f ( cp++, &b );
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
  ext = tree->ext+EPS;
  vertex->bbox.x0 -= ext;  vertex->bbox.x1 += ext;
  vertex->bbox.y0 -= ext;  vertex->bbox.y1 += ext;
  vertex->bbox.z0 -= ext;  vertex->bbox.z1 += ext;
  cp = vertex->ctlpoints;
  d = 0.0;
  for ( i = 1; i < ncp; i++ ) {
    SubtractPoints4f ( &cp[i], &cp[i-1], &dcp );
    e = DotProduct4f ( &dcp, &dcp );
    if ( e > d )
      d = e;
  }
  vertex->maxder = (float)tree->degree*sqrt ( d );
#undef EPS
} /*FindRCBoundingBoxf*/

static void UpdateRCBoundingBoxesf ( RBezCurveTreeVertexfp vertex )
{
  while ( vertex ) {
    if ( rbez_NarrowBBoxSumf ( &vertex->left->bbox, &vertex->right->bbox,
                               &vertex->bbox ) )
      vertex = vertex->up;
    else
      return;
  }
} /*UpdateRCBoundingBoxesf*/

RBezCurveTreefp rbez_NewRBezCurveTreef ( int object_id, short degree,
                                        float t0, float t1, float ext,
                                        CONST_ point4f *ctlpoints )
{
  RBezCurveTreefp tree;

  PKV_MALLOC ( tree, sizeof(RBezCurveTreef) );
  if ( tree ) {
    tree->object_id = object_id;
    tree->degree = degree;
    tree->cpsize = (degree+1)*sizeof(point4f);
    tree->ext = ext;
    tree->root = AllocRCTreeVertexf ( tree, NULL, t0, t1 );
    if ( tree->root ) {
      memcpy ( tree->root->ctlpoints, ctlpoints, tree->cpsize );
      FindRCBoundingBoxf ( tree, tree->root );
    }
    else {
      PKV_FREE ( tree );
      return NULL;
    }
  }
  return tree;
} /*rbez_NewRBezCurveTreef*/

static void r_DestroyRCTreef ( RBezCurveTreeVertexfp vertex )
{
  if ( vertex->right )
    r_DestroyRCTreef ( vertex->right );
  if ( vertex->left )
    r_DestroyRCTreef ( vertex->left );
  PKV_FREE ( vertex );
} /*r_DestroyRCTreef*/

void rbez_DestroyRBezCurveTreef ( RBezCurveTreefp tree )
{
  r_DestroyRCTreef ( tree->root );
  PKV_FREE ( tree );
} /*rbez_DestroyRBezCurveTreef*/

static void DivideRCVertexf ( RBezCurveTreefp tree, RBezCurveTreeVertexfp vertex )
{
  float h;

  h = 0.5*(vertex->t0 + vertex->t1);
  vertex->left  = AllocRCTreeVertexf ( tree, vertex, vertex->t0, h );
  vertex->right = AllocRCTreeVertexf ( tree, vertex, h, vertex->t1 );
  if ( !vertex->left || !vertex->right )
    return;
  memcpy ( vertex->right->ctlpoints, vertex->ctlpoints, tree->cpsize );
  mbs_BisectBC4f ( tree->degree,
                   vertex->right->ctlpoints, vertex->left->ctlpoints );
  FindRCBoundingBoxf ( tree, vertex->left );
  FindRCBoundingBoxf ( tree, vertex->right );
  UpdateRCBoundingBoxesf ( vertex );
} /*DivideRCVertexf*/

RBezCurveTreeVertexfp rbez_GetRBezCurveLeftVertexf ( RBezCurveTreefp tree,
                                                     RBezCurveTreeVertexfp vertex )
{
  if ( !vertex->left ) {
    if ( raybez_use_mutex )
      pthread_mutex_lock ( &raybez_mutex );
    if ( !vertex->left )
      DivideRCVertexf ( tree, vertex );
    if ( raybez_use_mutex )
      pthread_mutex_unlock ( &raybez_mutex );
  }
  return vertex->left;
} /*rbez_GetRBezCurveLeftVertexf*/

RBezCurveTreeVertexfp rbez_GetRBezCurveRightVertexf ( RBezCurveTreefp tree,
                                                      RBezCurveTreeVertexfp vertex )
{
  if ( !vertex->right ) {
    if ( raybez_use_mutex )
      pthread_mutex_lock ( &raybez_mutex );
    if ( !vertex->right )
      DivideRCVertexf ( tree, vertex );
    if ( raybez_use_mutex )
      pthread_mutex_unlock ( &raybez_mutex );
  }
  return vertex->right;
} /*rbez_GetRBezCurveRightVertexf*/

