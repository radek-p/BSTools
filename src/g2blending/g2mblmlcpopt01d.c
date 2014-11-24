
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012, 2014                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "bsmesh.h"
#include "egholed.h"
#include "g2blendingd.h"
#include "g2mblendingd.h"

#include "g2blprivated.h"
#include "g2mblprivated.h"
#include "g2mblmlprivated.h"

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g2mbl_CMPAssignMeshd ( mesh_ml_optdata *d,
                    int fnv, BSMvertex *fmv, int *fmvhei, point3d *fmvcp,
                    int fnhe, BSMhalfedge *fmhe,
                    int fnfac, BSMfacet *fmfac, int *fmfhei, byte *fmkcp,
                    int cnv, int rmnnz, index2 *rmnzi, double *rmnzc )
{
  void *sp;
  int  vi, vb, ei, eb, nvcp, i;
  char *vtag, *ftag;

  sp = pkv_GetScratchMemTop ();
        /* verify the mesh */
  if ( !bsm_CheckMeshIntegrity ( fnv, fmv, fmvhei, fnhe, fmhe,
                                 fnfac, fmfac, fmfhei, NULL, NULL ) )
    goto failure;
        /* all facets must be quadrangular */
  for ( i = 0; i < fnfac; i++ )
    if ( fmfac[i].degree != 4 )
      goto failure;
  PKV_MALLOC ( d->mvtag, fnv );
  if ( !d->mvtag )
    goto failure;
  vtag = d->mvtag;
  ftag = pkv_GetScratchMem ( fnfac );
  if ( !ftag )
    goto failure;
        /* boundary vertices must be of degree less than 4 */
  bsm_TagMesh ( fnv, fmv, fmvhei, fnhe, fmhe, fnfac, fmfac, fmfhei,
                vtag, ftag, &vi, &vb, &ei, &eb );
  for ( i = 0; i < fnv; i++ )
    if ( vtag[i] == 1 ) {
      if ( fmv[i].degree >= 4 )
        goto failure;
      else
        _g2mbl_TagBoundaryCondVert ( fnv, fmv, fmvhei, fnhe, fmhe, i, vtag );
    }
    else if ( fmv[i].degree < 3 || fmv[i].degree > GH_MAX_K )
      goto failure;
        /* tag the vertices fixed by the constraints */
  if ( fmkcp ) {
    for ( i = 0; i < fnv; i++ )
      if ( fmkcp[i] )
        vtag[i] = 1;
  }
        /* count the vertices fixed by constraints */
  for ( i = nvcp = 0; i < fnv; i++ )
    if ( !vtag[i] )
      nvcp ++;
  if ( !nvcp )    /* nothing to optimize */
    goto failure;
  d->nvcp = nvcp;
        /* attach the meshes to the data structure */
  d->nv    = fnv;
  d->nhe   = fnhe;
  d->nfac  = fnfac;
  d->mv    = fmv;
  d->mvhei = fmvhei;
  d->mvcp  = fmvcp;
  d->mhe   = fmhe;
  d->mfac  = fmfac;
  d->mfhei = fmfhei;

  d->cnv    = cnv;
        /* attach the refinement matrix */
  d->rmnnz  = rmnnz;
  d->rmnzi  = rmnzi;
  d->rmnzc  = rmnzc;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_CMPAssignMeshd*/

