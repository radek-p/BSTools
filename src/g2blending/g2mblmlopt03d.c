
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2011, 2012                            */
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
boolean _g2mbl_MLAssignMeshd ( mesh_ml_optdata *d,
                               int nv, BSMvertex *mv, int *mvhei, point3d *mvcp,
                               int nhe, BSMhalfedge *mhe,
                               int nfac, BSMfacet *mfac, int *mfhei,
                               byte *mkcp )
{
  void *sp;
  int  vi, vb, ei, eb, nvcp, i;
  char *mvtag, *ftag;

  sp = pkv_GetScratchMemTop ();
        /* verify the mesh */
  if ( !bsm_CheckMeshIntegrity ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei ) )
    goto failure;
          /* all facets must be quadrangles */
  for ( i = 0; i < nfac; i++ )
    if ( mfac[i].degree != 4 )
      goto failure;
  PKV_MALLOC ( d->mvtag, nv );
  if ( !d->mvtag )
    goto failure;
  mvtag = d->mvtag;
  ftag = pkv_GetScratchMem ( nfac );
  if ( !ftag )
    goto failure;
          /* boundary vertices must be of degree less than 4 */
  bsm_TagMesh ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                mvtag, ftag, &vi, &vb, &ei, &eb );
  for ( i = 0; i < nv; i++ )
    if ( mvtag[i] == 1 ) {
      if ( mv[i].degree >= 4 )
        goto failure;
      else
        _g2mbl_TagBoundaryCondVert ( nv, mv, mvhei, nhe, mhe, i, mvtag );
    }
    else if ( mv[i].degree < 3 || mv[i].degree > GH_MAX_K )
      goto failure;
        /* tag the vertices fixed by the constraints */
  if ( mkcp ) {
    for ( i = 0; i < nv; i++ )
      if ( mkcp[i] )
        mvtag[i] = 1;
  }
        /* count the vertices fixed by constraints */
  for ( i = nvcp = 0;  i < nv;  i++ )
    if ( !mvtag[i] )
      nvcp ++;
  if ( !nvcp )    /* nothing to optimize */
    goto failure;
  d->nvcp = nvcp;
        /* attach the mesh to the data structure */
  d->nv    = nv;
  d->mv    = mv;
  d->mvhei = mvhei;
  d->mvcp  = mvcp;
  d->nhe   = nhe;
  d->mhe   = mhe;
  d->nfac  = nfac;
  d->mfac  = mfac;
  d->mfhei = mfhei;

  pkv_SetScratchMemTop ( sp );  
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_MLAssignMeshd*/

