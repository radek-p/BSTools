
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2015                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "pkvaria.h"
#include "pkvthreads.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "camera.h"
#include "raybez.h"
#include "pkrender.h"

#include "pkrenderprivate.h"

rendertnode* NewLeafNode ( pkRenderer *rend, int k )
{
  rendertnode *node;

  node = malloc ( sizeof(rendertnode) );
  if ( !node )
    exit ( 1 );
  node->left = node->right = NULL;
  node->k = k;
  switch ( rend->obj_tab[k].type ) {
case obj_TRIANGLE:
    memcpy ( &node->bbox, &rend->obj_tab[k].triang.trdata->bbox, sizeof(Box3d) );
    break;
case obj_BSPATCH:
    memcpy ( &node->bbox, &rend->obj_tab[k].bsp.ptree->root->bbox, sizeof(Box3d) );
    break;
case obj_RBSPATCH:
    memcpy ( &node->bbox, &rend->obj_tab[k].rbsp.ptree->root->bbox, sizeof(Box3d) );
    break;
case obj_BEZCURVE:
    memcpy ( &node->bbox, &rend->obj_tab[k].bezc.ctree->root->bbox, sizeof(Box3d) );
    break;
case obj_RBEZCURVE:
    memcpy ( &node->bbox, &rend->obj_tab[k].rbezc.ctree->root->bbox, sizeof(Box3d) );
    break;
default:
    free ( node );
    return NULL;
  }
  return node;
} /*NewLeafNode*/

rendertnode* NewInnerNode ( rendertnode *left, rendertnode *right )
{
  rendertnode *node;

  node = malloc ( sizeof(rendertnode) );
  if ( !node )
    exit ( 1 );
  node->left = left;
  node->right = right;
  node->k = -1;
  rbez_FindSumBBoxd ( &left->bbox, &right->bbox, &node->bbox );
  return node;
} /*NewInnerNode*/

double BoxDiameter ( Box3d *box )
{
  double a, b, c;

  a = box->x1 - box->x0;
  b = box->y1 - box->y0;
  c = box->z1 - box->z0;
  return a*a + b*b + c*c;
} /*BoxDiameter*/

boolean CompNodePrio ( void *node1, void *node2 )
{
  double d1, d2;

  d1 = BoxDiameter ( &((rendertnode*)node1)->bbox );
  d2 = BoxDiameter ( &((rendertnode*)node2)->bbox );
  return d1 < d2;
} /*CompNodePrio*/

boolean BuildObjectTree ( pkRenderer *rend )
{
  void        *sp;
  int         nobjects;
  void        **pqueue;
  int         i, j, k, last;
  rendertnode *tn0, *tn1, *tn2;
  Box3d       box;
  double      d1, d2;

  sp = pkv_GetScratchMemTop ();
  rend->root = NULL;
  nobjects = rend->nobjects;
  if ( !nobjects )
    return false;
  pqueue = malloc ( nobjects*sizeof(void*) );
  if ( !pqueue )
    return false;
  for ( i = 0; i < nobjects; i++ )
    pqueue[i] = NewLeafNode ( rend, i );
  pkv_HeapOrder ( pqueue, nobjects, CompNodePrio );
  for ( i = nobjects-1; i > 0; i-- ) {
    tn0 = pqueue[0];
    last = i;
    pkv_HeapRemove ( pqueue, &last, 0, CompNodePrio );
    tn1 = pqueue[0];
    rbez_FindSumBBoxd ( &tn0->bbox, &tn1->bbox, &box );
    d1 = BoxDiameter ( &box );
    k = 0;
    for ( j = 1; j <= last; j++ ) {
      tn2 = pqueue[j];
      rbez_FindSumBBoxd ( &tn0->bbox, &tn2->bbox, &box );
      d2 = BoxDiameter ( &box );
      if ( d2 < d1 ) {
        k = j;
        d1 = d2;
        tn1 = tn2;
      }
    }
    pkv_HeapRemove ( pqueue, &last, k, CompNodePrio );
    tn2 = NewInnerNode ( tn0, tn1 );
    pkv_HeapInsert ( pqueue, &last, (void*)tn2, CompNodePrio );
  }
  rend->root = pqueue[0];
  pkv_SetScratchMemTop ( sp );
  return true;
} /*BuildObjectTree*/

void DestroyTNodeTree ( rendertnode *node )
{
  if ( node ) {
    DestroyTNodeTree ( node->left );
    DestroyTNodeTree ( node->right );
    free ( node );
  }
} /*DestroyTNodeTree*/

