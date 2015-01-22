
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014, 2015                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "raybez.h"

#include "raybezprivate.h"
#include "raybezprivated.h"

#define BOX_EPS 1.0e-10

/* ////////////////////////////////////////////////////////////////////////// */
static RBezPatchTreeVertexdp
  AllocRBezTreeVertexd ( RBezPatchTreedp tree,
                         RBezPatchTreeVertexdp up,
                         double u0, double u1, double v0, double v1,
                         boolean spline )
{
  int                   size;
  RBezPatchTreeVertexdp vertex;

  
  if ( spline )  /* no control points */
    size = sizeof(RBezPatchTreeVertexd);
  else
    size = sizeof(BezPatchTreeVertexd) + tree->cpsize;
  PKV_MALLOC ( vertex, size );
  if ( vertex ) {
    vertex->u0 = u0;  vertex->u1 = u1;
    vertex->v0 = v0;  vertex->v1 = v1;
    vertex->left = vertex->right = NULL;
    vertex->vertex_colour = RAYBEZ_WHITE;
    vertex->tag = 0;
    if ( (vertex->up = up) != NULL )
      vertex->level = (short)(up->level + 1);
    else
      vertex->level = 0;
    if ( spline )
      vertex->ctlpoints = NULL;
    else
      vertex->ctlpoints = (point4d*)((char*)vertex+sizeof(RBezPatchTreeVertexd));
    vertex->nvcpoints = NULL;
  }
  return vertex;
} /*AllocRBezTreeVertexd*/

static void RBezPatchGetMaxDerd ( int n, int m, int pitch, point4d* cp,
                                  double *du, double *dv )
{
  int      i, j, k;
  vector4d v;
  double    _du, _dv, d;

  _du = _dv = 0.0;
  for ( i = 0; i < n; i++ )
    for ( j = 0, k = i*pitch;  j <= m;  j++, k++ ) {
      SubtractPoints4d ( &cp[k+pitch], &cp[k], &v );
      d = DotProduct4d ( &v, &v );
      _du = max ( _du, d );
    }
  for ( i = 0; i <= n; i++ )
    for ( j = 0, k = i*pitch;  j < m;  j++, k++ ) {
      SubtractPoints4d ( &cp[k+1], &cp[k], &v );
      d = DotProduct4d ( &v, &v );
      _dv = max ( _dv, d );
    }
  *du = (double)n * sqrt ( _du );
  *dv = (double)m * sqrt ( _dv );
} /*RBezPatchGetMaxDerd*/

static void RBezPatchGetDivDird ( unsigned char n, unsigned int lknu,
                                  unsigned char m, unsigned int lknv,
                                  unsigned int pitch, point4d *ctlpoints,
                                  unsigned char *divdir, double *maxdl )
{
  int     ku, kv;
  int     i, j;
  point4d *cp, *cq;
  double   du, dv, maxdlu, maxdlv;

  ku = (lknu-n)/(n+1);
  kv = (lknv-m)/(m+1);
  maxdlu = maxdlv = 0.0;
  for ( i = 0; i < ku; i++ ) {
    cp = &ctlpoints[i*(n+1)*pitch];
    for ( j = 0; j < kv; j++ ) {
      cq = &cp[j*(m+1)];
      RBezPatchGetMaxDerd ( n, m, pitch, cq, &du, &dv );
      if ( du > maxdlu ) maxdlu = du;
      if ( dv > maxdlv ) maxdlv = dv;
    }
  }
  maxdlu *= (double)ku;
  maxdlv *= (double)kv;
  *maxdl = max ( maxdlu, maxdlv );
  if ( ku == 1 ) {
    if ( kv == 1 )
      goto by_length;
      else *divdir = 1;
  }
  else if ( kv == 1 )
    *divdir = 0;
  else {
by_length:
    if ( maxdlu > maxdlv ) { *divdir = 0; }
                      else { *divdir = 1; }
  }
} /*RBezPatchGetDivDird*/

static void RBezPatchInitialDivision ( RBezPatchTreedp tree,
                                       RBezPatchTreeVertexdp vertex,
                                       int n, int lknu, double *knu,
                                       int m, int lknv, double *knv,
                                       int pitch, point4d *cp )
{
  int                   ku, kv, kk;
  unsigned char         divdir;
  double                 maxder, uv;
  RBezPatchTreeVertexdp left, right;

  ku = (lknu-n)/(n+1);
  kv = (lknv-m)/(m+1);
  RBezPatchGetDivDird ( n, lknu, m, lknv, pitch, cp, &divdir, &maxder );
  vertex->maxder = maxder;
  vertex->divdir = divdir;
  switch ( divdir ) {
case 0:  /* divide the "u" variable interval */
    kk = ku / 2;
    uv = knu[kk*(n+1)];
    left = AllocRBezTreeVertexd ( tree, vertex, vertex->u0, uv,
                                  vertex->v0, vertex->v1, kk > 1 || kv > 1 );
    right = AllocRBezTreeVertexd ( tree, vertex, uv, vertex->u1,
                                   vertex->v0, vertex->v1, ku-kk > 1 || kv > 1 );
    if ( !left || !right )
      goto failure;
    if ( kk == 1 && kv == 1 ) {  /* left vertex is a Bezier patch */
      pkv_Selectc ( n+1, (m+1)*sizeof(point4d),
                    pitch*sizeof(point4d), (m+1)*sizeof(point4d),
                    (char*)cp, (char*)left->ctlpoints );
      RBezPatchGetDivDird ( n, 2*n+1, m, 2*m+1, pitch, cp, 
                            &left->divdir, &left->maxder );
      mbs_BCHornerNvP3Rd ( n, m, left->ctlpoints, 0.5, 0.5,
                           &left->pcent, &left->nvcent );
      rbez_FindCPBoundingBox3Rd ( 1, (n+1)*(m+1), 0, left->ctlpoints, BOX_EPS,
                                  &left->bbox );
    }
    else {  /* left vertex consists of more than one Bezier patch */
      RBezPatchInitialDivision ( tree, left, n, (kk+1)*(n+1)-1, knu,
                                 m, lknv, knv, pitch, cp );
    }
    if ( ku-kk == 1 && kv == 1 ) { /* right vertex is a Bezier patch */
      pkv_Selectc ( n+1, (m+1)*sizeof(point4d),
                    pitch*sizeof(point4d), (m+1)*sizeof(point4d),
                    (char*)&cp[kk*(n+1)*pitch], (char*)right->ctlpoints );
      RBezPatchGetDivDird ( n, 2*n+1, m, 2*m+1, pitch, &cp[kk*(n+1)*pitch],
                            &right->divdir, &right->maxder );
      mbs_BCHornerNvP3Rd ( n, m, right->ctlpoints, 0.5, 0.5,
                           &right->pcent, &right->nvcent );
      rbez_FindCPBoundingBox3Rd ( 1, (n+1)*(m+1), 0, right->ctlpoints, BOX_EPS,
                                  &right->bbox );
    }
    else {  /* right vertex consists of more than one Bezier patch */
      RBezPatchInitialDivision ( tree, right, n, (ku-kk+1)*(n+1)-1, &knu[kk*(n+1)],
                                 m, lknv, knv, pitch, &cp[kk*(n+1)*pitch] );
    }
    break;

case 1:  /* divide the "v" variable interval */
    kk = kv/2;
    uv = knv[kk*(m+1)];
    left = AllocRBezTreeVertexd ( tree, vertex, vertex->u0, vertex->u1,
                                  vertex->v0, uv, ku > 1 || kk > 1 );
    right = AllocRBezTreeVertexd ( tree, vertex, vertex->u0, vertex->u1,
                                   uv, vertex->v1, ku > 1 || kv-kk > 1 );
    if ( !left || !right )
      goto failure;
    if ( ku == 1 && kk == 1 ) {  /* left vertex is a Bezier patch */
      pkv_Selectc ( n+1, (m+1)*sizeof(point4d),
                    pitch*sizeof(point4d), (m+1)*sizeof(point4d),
                    (char*)cp, (char*)left->ctlpoints );
      RBezPatchGetDivDird ( n, 2*n+1, m, 2*m+1, pitch, cp,
                            &left->divdir, &left->maxder );
      mbs_BCHornerNvP3Rd ( n, m, left->ctlpoints, 0.5, 0.5,
                           &left->pcent, &left->nvcent );
      rbez_FindCPBoundingBox3Rd ( 1, (n+1)*(m+1), 0, left->ctlpoints, BOX_EPS,
                                  &left->bbox );
    }
    else {
      RBezPatchInitialDivision ( tree, left, n, lknu, knu,
                                 m, (kk+1)*(m+1)-1, knv, pitch, cp );
    }
    if ( ku == 1 && kv-kk == 1 ) { /* right vertex is a Bezier patch */
      pkv_Selectc ( n+1, (m+1)*sizeof(point4d),
                    pitch*sizeof(point4d), (m+1)*sizeof(point4d),
                    (char*)&cp[kk*(m+1)], (char*)right->ctlpoints );
      RBezPatchGetDivDird ( n, 2*n+1, m, 2*m+1, pitch, &cp[kk*(m+1)],
                            &right->divdir, &right->maxder );
      mbs_BCHornerNvP3Rd ( n, m, right->ctlpoints, 0.5, 0.5,
                           &right->pcent, &right->nvcent );
      rbez_FindCPBoundingBox3Rd ( 1, (n+1)*(m+1), 0, right->ctlpoints, BOX_EPS,
                                  &right->bbox );
    }
    else {  /* right vertex consists of more than one Bezier patch */
      RBezPatchInitialDivision ( tree, right, n, lknu, knu,
                                 m, (kv-kk+1)*(m+1)-1, &knv[kk*(m+1)],
                                 pitch, &cp[kk*(m+1)] );
    }
    break;

default:
    goto failure;
  }
  vertex->left = left;
  vertex->right = right;
  vertex->tag = 2;
  rbez_FindSumBBoxd ( &left->bbox, &right->bbox, &vertex->bbox );
  return;

failure:
  if ( left )  PKV_FREE ( left );
  if ( right ) PKV_FREE ( right );
} /*RBezPatchInitialDivision*/

static void RUpdateBoundingBoxesd ( RBezPatchTreeVertexdp vertex )
{
  while ( vertex ) {
    if ( rbez_NarrowBBoxSumd ( &vertex->left->bbox, &vertex->right->bbox,
                               &vertex->bbox ) )
      vertex = vertex->up;
    else
      return;
  }
} /*RUpdateBoundingBoxesd*/

RBezPatchTreedp
    rbez_NewRBSPatchTreed ( int object_id,
                   unsigned char n, unsigned int lknu, CONST_ double *knotsu,
                   unsigned char m, unsigned int lknv, CONST_ double *knotsv,
                   unsigned int pitch, CONST_ point4d *ctlpoints )
{
  void                  *sp;
  RBezPatchTreedp       tree;
  RBezPatchTreeVertexdp root;
  int                   _pitch, _lknu, _lknv, ku, kv;
  double                 *knu, *knv;
  point4d               *bcp;

      /* pitch is expressed in the units sizeof(double) */
      /* _pitch is in units sizeof(point4d) */
  sp = pkv_GetScratchMemTop ();
  tree = NULL;  root = NULL;
  ku = mbs_NumKnotIntervalsd ( n, lknu, knotsu );
  kv = mbs_NumKnotIntervalsd ( m, lknv, knotsv );
  _pitch = (m+1)*kv;
  bcp = pkv_GetScratchMem ( _pitch*ku*(n+1)*sizeof(point4d) );
  knu = pkv_GetScratchMemd ( (ku+1)*(n+1)+(kv+1)*(m+1) );
  if ( !bcp || !knu )
    goto failure;
  knv = &knu[(ku+1)*(n+1)];
  mbs_BSPatchToBezd ( 4, n, lknu, knotsu, m, lknv, knotsv,
                      4*pitch, &ctlpoints[0].x,
                      &ku, &_lknu, knu, &kv, &_lknv, knv, 4*_pitch, &bcp[0].x );
  PKV_MALLOC ( tree, sizeof(RBezPatchTreed) );
  if ( !tree )
    goto failure;
  tree->object_id = object_id;
  tree->n = n;
  tree->m = m;
  tree->cpsize = (n+1)*(m+1)*sizeof(point4d);
  tree->nvsize = 4*n*m*sizeof(vector3d);
  tree->root = root = AllocRBezTreeVertexd ( tree, NULL,
                         knotsu[n], knotsu[lknu-n], knotsv[m], knotsv[lknv-m],
                         lknu > n+n+1 || lknv > m+m+1 );
  if ( !root )
    goto failure;
  if ( ku > 1 || kv > 1 )
    RBezPatchInitialDivision ( tree, root, n, _lknu, knu, m, _lknv, knv,
                               _pitch, bcp );
  else {
    pkv_Selectc ( n+1, (m+1)*sizeof(point4d),
                  pitch*sizeof(point4d), (m+1)*sizeof(point4d),
                  (char*)bcp, (char*)root->ctlpoints );
    RBezPatchGetDivDird ( n, _lknu, m, _lknv, m+1, root->ctlpoints,
                          &root->divdir, &root->maxder );
    mbs_BCHornerNvP3Rd ( n, m, root->ctlpoints, 0.5, 0.5,
                         &root->pcent, &root->nvcent );
    rbez_FindCPBoundingBox3Rd ( 1, (n+1)*(m+1), 0, root->ctlpoints, BOX_EPS,
                                &root->bbox );
  }

  pkv_SetScratchMemTop ( sp );
  return tree;

failure:
  if ( root ) PKV_FREE ( root );
  if ( tree ) PKV_FREE ( tree );
  pkv_SetScratchMemTop ( sp );
  return NULL;
} /*rbez_NewRBSPatchTreed*/

RBezPatchTreedp
  rbez_NewRBezPatchTreed ( int object_id,
                           unsigned char n, unsigned char m,
                           double u0, double u1, double v0, double v1,
                           CONST_ point4d *ctlpoints )
{
  void            *sp;
  RBezPatchTreedp tree;
  double           *knotsu, *knotsv;
  int             i;

  sp = pkv_GetScratchMemTop ();
  knotsu = pkv_GetScratchMemd ( 2*(n+m+2) );
  if ( knotsu ) {
    knotsv = &knotsu[2*n+2];
    for ( i = 0; i <= n; i++ )
      knotsu[i] = u0,  knotsu[i+n+1] = u1;
    for ( i = 0; i <= m; i++ )
      knotsv[i] = v0,  knotsv[i+m+1] = v1;
    tree = rbez_NewRBSPatchTreed ( object_id, n, n+n+1, knotsu, m, m+m+1, knotsv,
                                   m+1, ctlpoints );
  }
  else
    tree = NULL;
  pkv_SetScratchMemTop ( sp );
  return tree;
} /*rbez_NewRBezPatchTreed*/

static void r_DestroyRBezPTreed ( RBezPatchTreeVertexdp vertex )
{
  if ( vertex->right )
    r_DestroyRBezPTreed ( vertex->right );
  if ( vertex->left )
    r_DestroyRBezPTreed ( vertex->left );
  if ( vertex->nvcpoints )
    PKV_FREE ( vertex->nvcpoints );
  PKV_FREE ( vertex );
} /*r_DestroyRBezPTreed*/

void rbez_DestroyRBezPatchTreed ( RBezPatchTreedp tree )
{
  r_DestroyRBezPTreed ( tree->root );
  PKV_FREE ( tree );
} /*rbez_DestroyRBezPatchTreed*/

static void RBezPDivideVertexd ( RBezPatchTreedp tree,
                                 RBezPatchTreeVertexdp vertex )
{
  double h;
  int    n, m, ncp;

  n = tree->n;
  m = tree->m;
  ncp = (n+1)*(m+1);
        /* allocate subtree roots */
  if ( !vertex->divdir ) {
    h = 0.5*(vertex->u0 + vertex->u1);
    vertex->left = AllocRBezTreeVertexd ( tree, vertex,
                          vertex->u0, h, vertex->v0, vertex->v1, false );
    vertex->right = AllocRBezTreeVertexd ( tree, vertex,
                          h, vertex->u1, vertex->v0, vertex->v1, false );
    if ( !vertex->left || !vertex->right )
      goto dealloc;
    memcpy ( vertex->right->ctlpoints, vertex->ctlpoints, tree->cpsize );
    mbs_BisectBP4ud ( n, m, vertex->right->ctlpoints, vertex->left->ctlpoints  );
  }
  else {
    h = 0.5*(vertex->v0 + vertex->v1);
    vertex->left = AllocRBezTreeVertexd ( tree, vertex,
                          vertex->u0, vertex->u1, vertex->v0, h, false );
    vertex->right = AllocRBezTreeVertexd ( tree, vertex,
                          vertex->u0, vertex->u1, h, vertex->v1, false );
    if ( !vertex->left || !vertex->right )
      goto dealloc;
    memcpy ( vertex->right->ctlpoints, vertex->ctlpoints, tree->cpsize );
    mbs_BisectBP4vd ( n, m, vertex->right->ctlpoints, vertex->left->ctlpoints  );
  }
  if ( !mbs_BCHornerNvP3Rd ( n, m, vertex->left->ctlpoints,
                  0.5, 0.5, &vertex->left->pcent, &vertex->left->nvcent ) )
    goto dealloc;
  if ( !mbs_BCHornerNvP3Rd ( n, m, vertex->right->ctlpoints,
                  0.5, 0.5, &vertex->right->pcent, &vertex->right->nvcent ) )
    goto dealloc;
  RBezPatchGetDivDird ( n, n+n+1, m, m+m+1, m+1, vertex->left->ctlpoints,
                        &vertex->left->divdir, &vertex->left->maxder );
  RBezPatchGetDivDird ( n, n+n+1, m, m+m+1, m+1, vertex->right->ctlpoints,
                        &vertex->right->divdir, &vertex->right->maxder );
  rbez_FindCPBoundingBox3Rd ( 1, ncp, 0, vertex->left->ctlpoints, BOX_EPS,
                              &vertex->left->bbox );
  rbez_FindCPBoundingBox3Rd ( 1, ncp, 0, vertex->right->ctlpoints, BOX_EPS,
                              &vertex->right->bbox );
  RUpdateBoundingBoxesd ( vertex );
  return;

dealloc:
  if ( vertex->left ) PKV_FREE ( vertex->left );
  if ( vertex->right ) PKV_FREE ( vertex->right );
} /*RBezPDivideVertexd*/

RBezPatchTreeVertexdp
    rbez_GetRBezLeftVertexd ( RBezPatchTreedp tree,  
                              RBezPatchTreeVertexdp vertex )
{
  if ( vertex->tag < 2 )
    raybez_DivideTreeVertex ( tree, vertex, &vertex->tag,
                              (divide_vertex_proc)&RBezPDivideVertexd );
  return vertex->left;
} /*rbez_GetRBezLeftVertexd*/

RBezPatchTreeVertexdp
    rbez_GetRBezRightVertexd ( RBezPatchTreedp tree, 
                               RBezPatchTreeVertexdp vertex )
{
  if ( vertex->tag < 2 )
    raybez_DivideTreeVertex ( tree, vertex, &vertex->tag,
                              (divide_vertex_proc)&RBezPDivideVertexd );
  return vertex->right;
} /*rbez_GetRBezRightVertexd*/

