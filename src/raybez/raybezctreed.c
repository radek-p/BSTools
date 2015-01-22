
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012, 2015                            */
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

#include "raybezprivate.h"

/* ////////////////////////////////////////////////////////////////////////// */
static BezCurveTreeVertexdp
    AllocCTreeVertexd ( BezCurveTreedp tree,
                        BezCurveTreeVertexdp up,
                        double t0, double t1, boolean spline )
{
  int                  size;
  BezCurveTreeVertexdp vertex;

  if ( spline )
    size = sizeof(BezCurveTreeVertexd);
  else
    size = sizeof(BezCurveTreeVertexd) + tree->cpsize;
  PKV_MALLOC ( vertex, size );
  if ( vertex ) {
    vertex->t0 = t0;
    vertex->t1 = t1;
    vertex->left = vertex->right = NULL;
    vertex->up = up;
    if ( up )
      vertex->level = (short)(up->level + 1);
    else
      vertex->level = 0;
    if ( !spline )
      vertex->ctlpoints = (point3d*)((char*)vertex+sizeof(BezCurveTreeVertexd));
    else
      vertex->ctlpoints = NULL;
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
    i = 1;
  }
  else {
    b = *cp;  cp++;
    pkv_Sort2d ( &a.x, &b.x );  vertex->bbox.x0 = a.x;  vertex->bbox.x1 = b.x;
    pkv_Sort2d ( &a.y, &b.y );  vertex->bbox.y0 = a.y;  vertex->bbox.y1 = b.y;
    pkv_Sort2d ( &a.z, &b.z );  vertex->bbox.z0 = a.z;  vertex->bbox.z1 = b.z;
    i = 2;
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

static void BezCurveInitialDivision ( BezCurveTreedp tree,
                                      BezCurveTreeVertexdp vertex,
                                      int degree, int lastknot, double *knots,
                                      point3d *cp )
{
  int                  ku, kk;
  BezCurveTreeVertexdp left, right;
  double               uu;

  left = right = NULL;
  ku = (lastknot-degree)/(degree+1);
  if ( ku == 1 ) {
    memcpy ( vertex->ctlpoints, cp, tree->cpsize );
    FindCBoundingBoxd ( tree, vertex );
    mbs_BCHornerC3d ( degree, cp, 0.5, &vertex->ccent );
  }
  else {
    kk = ku / 2;
    uu = knots[kk*(degree+1)];
    left = AllocCTreeVertexd ( tree, vertex, vertex->t0, uu, kk > 1 );
    right = AllocCTreeVertexd ( tree, vertex, uu, vertex->t1, ku-kk > 1 );
    if ( !left || !right )
      goto failure;
    vertex->left = left;
    vertex->right = right;
    vertex->tag = 2;
    BezCurveInitialDivision ( tree, left, degree,
                              (kk+1)*(degree+1)-1, knots, cp );
    BezCurveInitialDivision ( tree, right, degree,
                              (ku-kk+1)*(degree+1)-1, &knots[kk*(degree+1)],
                              &cp[kk*(degree+1)]);
    rbez_FindSumBBoxd ( &left->bbox, &right->bbox, &vertex->bbox );
  }
  return;

failure:
  if ( left )  PKV_FREE ( left );
  if ( right ) PKV_FREE ( right );
} /*BezCurveInitialDivision*/

BezCurveTreedp rbez_NewBSCurveTreed ( int object_id,
                                      short degree, int lastknot, double *knots,
                                      point3d *ctlpoints, double ext )
{
  void                 *sp;
  BezCurveTreedp       tree;
  BezCurveTreeVertexdp root;
  int                  _lkn, ku;
  double               *knu;
  point3d              *bcp;

  sp = pkv_GetScratchMemTop ();
  tree = NULL;  root = NULL;
  ku = mbs_NumKnotIntervalsd ( degree, lastknot, knots );
  bcp = pkv_GetScratchMem ( ku*(degree+1)*sizeof(point3d) );
  knu = pkv_GetScratchMemd ( (ku+1)*(degree+1) );
  if ( !bcp || !knu )
    goto failure;
  mbs_BSToBezC3d ( degree, lastknot, knots, ctlpoints, &ku, &_lkn, knu, bcp );
  PKV_MALLOC ( tree, sizeof(BezCurveTreed) );
  if ( !tree )
    goto failure;
  tree->object_id = object_id;
  tree->degree = degree;
  tree->cpsize = (degree+1)*sizeof(point3d);
  tree->ext = ext;
  tree->root = root = AllocCTreeVertexd ( tree, NULL, knu[degree], knu[_lkn-degree],
                                          ku > 1 );
  if ( !root )
    goto failure;
  BezCurveInitialDivision ( tree, root, degree, _lkn, knu, bcp );

  pkv_SetScratchMemTop ( sp );
  return tree;

failure:
  if ( root ) PKV_FREE ( root );
  if ( tree ) PKV_FREE ( tree );
  pkv_SetScratchMemTop ( sp );
  return NULL;
} /*rbez_NewBSCurveTreef*/

BezCurveTreedp rbez_NewBezCurveTreed ( int object_id, short degree,
                                       double t0, double t1,
                                       CONST_ point3d *ctlpoints, double ext )
{
  BezCurveTreedp tree;

  PKV_MALLOC ( tree, sizeof(BezCurveTreed) );
  if ( tree ) {
    tree->object_id = object_id;
    tree->degree = degree;
    tree->cpsize = (degree+1)*sizeof(point3d);
    tree->ext = ext;
    tree->root = AllocCTreeVertexd ( tree, NULL, t0, t1, false );
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
  vertex->left  = AllocCTreeVertexd ( tree, vertex, vertex->t0, h, false );
  vertex->right = AllocCTreeVertexd ( tree, vertex, h, vertex->t1, false );
  if ( !vertex->left || !vertex->right )
    return;
  memcpy ( vertex->right->ctlpoints, vertex->ctlpoints, tree->cpsize );
  mbs_BisectBC3d ( tree->degree,
                   vertex->right->ctlpoints, vertex->left->ctlpoints );
  FindCBoundingBoxd ( tree, vertex->left );
  FindCBoundingBoxd ( tree, vertex->right );
  UpdateCBoundingBoxesd ( vertex );
} /*DivideCVertexd*/

BezCurveTreeVertexdp rbez_GetBezCurveLeftVertexd ( BezCurveTreedp tree,
                                                   BezCurveTreeVertexdp vertex )
{
  if ( vertex->tag < 2 )
    raybez_DivideTreeVertex ( tree, vertex, &vertex->tag,
                              (divide_vertex_proc)&DivideCVertexd );
  return vertex->left;
} /*rbez_GetBezCurveLeftVertexd*/

BezCurveTreeVertexdp rbez_GetBezCurveRightVertexd ( BezCurveTreedp tree,
                                                    BezCurveTreeVertexdp vertex )
{
  if ( vertex->tag < 2 )
    raybez_DivideTreeVertex ( tree, vertex, &vertex->tag,
                              (divide_vertex_proc)&DivideCVertexd );
  return vertex->right;
} /*rbez_GetBezCurveRightVertexd*/

