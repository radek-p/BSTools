
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

typedef struct {
    pkRenderer *rend;
    int        dens;
    int        *jobs;
    boolean    ok;
    float      fmin, fmax;
  } minmax_struct;

static void _BSPatchMinMaxShapeFunc ( pkRenderer *rend,
                                      int dens, int id, renderobj *obj,
                                      BezPatchTreeVertexd *vertex,
                                      double *fmin, double *fmax )
{
  BezPatchTreed    *tree;
  point3d          *cp;
  RayObjectIntersd roint;
  ray3d            ray;
  int              m, n, j, k;
  double           f;

  if ( vertex->ctlpoints ) {
    tree = obj->bsp.ptree;
    n = tree->n;
    m = tree->m;
    cp = vertex->ctlpoints;

    ray.p = rend->CPos.position;
    roint.object_id = id;
    roint.extra_info = (void*)vertex;
    if ( *fmin > *fmax ) {
      roint.u = vertex->u0;
      roint.v = vertex->v0;
      mbs_BCHornerNvP3d  ( m, n, cp, 0.0, 0.0, &roint.p, &roint.nv );
      SubtractPoints3d ( &roint.p, &ray.p, &ray.v );
      roint.t = sqrt ( DotProduct3d ( &ray.v, &ray.v ) );
      MultVector3d ( 1.0/roint.t, &ray.v, &ray.v );
      *fmin = *fmax = cShapeFunc ( rend, &ray, obj, &roint );
    }
    for ( j = 0; j <= dens; j++ ) {
      roint.u = vertex->u0 + (double)j/(double)dens*(vertex->u1-vertex->u0);
      for ( k = 0; k <= dens; k++ ) {
        roint.v = vertex->v0 + (double)k/(double)dens*(vertex->v1-vertex->v0);
        mbs_BCHornerNvP3d  ( m, n, cp, roint.u, roint.v, &roint.p, &roint.nv );
        SubtractPoints3d ( &roint.p, &ray.p, &ray.v );
        roint.t = sqrt ( DotProduct3d ( &ray.v, &ray.v ) );
        MultVector3d ( 1.0/roint.t, &ray.v, &ray.v );
        f = cShapeFunc ( rend, &ray, obj, &roint );
        if ( f < *fmin )      *fmin = f;
        else if ( f > *fmax ) *fmax = f;
      }
    }
  }
  else {
    if ( vertex->left )
      _BSPatchMinMaxShapeFunc ( rend, dens, id, obj, vertex->left, fmin, fmax );
    if (vertex->right )
      _BSPatchMinMaxShapeFunc ( rend, dens, id, obj, vertex->right, fmin, fmax );
  }
} /*_BSPatchMinMaxShapeFunc*/

static void _RBSPatchMinMaxShapeFunc ( pkRenderer *rend,
                                       int dens, int id, renderobj *obj,
                                       RBezPatchTreeVertexd *vertex,
                                       double *fmin, double *fmax )
{
  RBezPatchTreed   *tree;
  point4d          *cp;
  RayObjectIntersd roint;
  ray3d            ray;
  int              m, n, j, k;
  double           f;

  if ( vertex->ctlpoints ) {
    tree = obj->rbsp.ptree;
    n = tree->n;
    m = tree->m;
    cp = vertex->ctlpoints;

    ray.p = rend->CPos.position;
    roint.object_id = id;
    roint.extra_info = (void*)vertex;
    if ( *fmin > *fmax ) {
      roint.u = vertex->u0;
      roint.v = vertex->v0;
      mbs_BCHornerNvP3Rd  ( m, n, cp, 0.0, 0.0, &roint.p, &roint.nv );
      SubtractPoints3d ( &roint.p, &ray.p, &ray.v );
      roint.t = sqrt ( DotProduct3d ( &ray.v, &ray.v ) );
      MultVector3d ( 1.0/roint.t, &ray.v, &ray.v );
      *fmin = *fmax = cShapeFunc ( rend, &ray, obj, &roint );
    }
    for ( j = 0; j <= dens; j++ ) {
      roint.u = vertex->u0 + (double)j/(double)dens*(vertex->u1-vertex->u0);
      for ( k = 0; k <= dens; k++ ) {
        roint.v = vertex->v0 + (double)k/(double)dens*(vertex->v1-vertex->v0);
        mbs_BCHornerNvP3Rd  ( m, n, cp, roint.u, roint.v, &roint.p, &roint.nv );
        SubtractPoints3d ( &roint.p, &ray.p, &ray.v );
        roint.t = sqrt ( DotProduct3d ( &ray.v, &ray.v ) );
        MultVector3d ( 1.0/roint.t, &ray.v, &ray.v );
        f = cShapeFunc ( rend, &ray, obj, &roint );
        if ( f < *fmin )      *fmin = f;
        else if ( f > *fmax ) *fmax = f;
      }
    }
  }
  else {
    if ( vertex->left )
      _RBSPatchMinMaxShapeFunc ( rend, dens, id, obj, vertex->left, fmin, fmax );
    if (vertex->right )
      _RBSPatchMinMaxShapeFunc ( rend, dens, id, obj, vertex->right, fmin, fmax );
  }
} /*_RBSPatchMinMaxShapeFunc*/

static boolean _FindMinMaxShapeFunc ( void *usrdata, int4 *jobnum )
{
  minmax_struct *mms;
  pkRenderer    *rend;
  int           i, i0, i1;
  double        fmin, fmax;

  mms = (minmax_struct*)usrdata;
  rend = mms->rend;
  i = jobnum->x;
  i0 = mms->jobs[i];  i1 = mms->jobs[i+1];
  fmin = 1.0;  fmax = -1.0;
  for ( i = i0; i < i1; i++ ) {
    switch ( rend->obj_tab[i].type ) {
  case obj_BSPATCH:
      _BSPatchMinMaxShapeFunc ( rend, mms->dens, i, &rend->obj_tab[i],
                                rend->obj_tab[i].bsp.ptree->root, &fmin, &fmax );
      break;
  case obj_RBSPATCH:
      _RBSPatchMinMaxShapeFunc ( rend, mms->dens, i, &rend->obj_tab[i],
                                 rend->obj_tab[i].rbsp.ptree->root, &fmin, &fmax );
      break;
  default:
      break;
    }
  }
        /* store the minimal and maximal values found so far */
  if ( fmin < fmax ) {
    if ( mms->ok ) {
      if ( fmin < mms->fmin ) mms->fmin = fmin;
      if ( fmax > mms->fmax ) mms->fmax = fmax;
    }
    else {
      mms->fmin = fmin;
      mms->fmax = fmax;
      mms->ok = true;
    }
  }
  return true;
} /*_FindMinMaxShapeFunc*/

void FindMinMaxShapeFunc ( pkRenderer *rend, int dens )
{
  void          *sp;
  minmax_struct mms;
  int           i, d, nthr;
  int4          jobsize;
  boolean       success;

  sp = pkv_GetScratchMemTop ();
  pkv_Tic ( &rend->tic );

  nthr = rend->nthr > 1 ? rend->nthr : 1;
  nthr = min ( rend->nobjects, nthr );
  mms.rend = rend;
  mms.dens = dens;
  mms.ok = false;
  mms.fmin = mms.fmax = 0.0;
  mms.jobs = pkv_GetScratchMemi ( nthr+1 );
  if ( !mms.jobs || rend->nobjects <= 0 )
    goto failure;
  mms.jobs[0] = 0;
  mms.jobs[nthr] = rend->nobjects;
  if ( nthr > 1 ) {

printf ( "A: nthr = %d, nobjects = %d\n", nthr, rend->nobjects );

    d = (rend->nobjects+nthr-1)/nthr;
    for ( i = 1; i < nthr; i++ )
      mms.jobs[i] = mms.jobs[i-1]+d;
    mms.jobs[nthr] = rend->nobjects;
    jobsize.x = nthr;
    pkv_SetPThreadsToWork ( 1, &jobsize, nthr, PKRENDER_STACK, PKRENDER_SCRATCHMEM,
                            (void*)&mms, _FindMinMaxShapeFunc, NULL, NULL,
                            &success );
  }
  else {

printf ( "B\n" );

    jobsize.x = 0;
    _FindMinMaxShapeFunc ( (void*)&mms, &jobsize );
  }

failure:
  rend->minsf = mms.fmin + (mms.fmax-mms.fmin)*rend->cfrange[0];
  rend->maxsf = mms.fmin + (mms.fmax-mms.fmin)*rend->cfrange[1];
  if ( rend->minsf >= rend->maxsf ) {
    rend->minsf -= 2.0;
    rend->maxsf += 3.0;
  }
  rend->ticks1 = pkv_Toc ( &rend->tic );
  pkv_SetScratchMemTop ( sp );
} /*FindMinMaxShapeFunction*/

