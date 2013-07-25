
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2012                            */
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

/* ///////////////////////////////////////////////////////////////////////// */
boolean g2mbl_InitBlSurfaceOptLMTd ( int nv, BSMvertex *mv, int *mvhei,
                                     point3d *mvcp, int nhe, BSMhalfedge *mhe,
                                     int nfac, BSMfacet *mfac, int *mfhei,
                                     byte *mkcp,
                                     double C, double dO, double dM,
                                     int nkn1, int nkn2, void **data )
{
  void             *sp;
  mesh_lmt_optdata *d;

  sp = pkv_GetScratchMemTop ();
  *data = d = NULL;
  PKV_MALLOC ( *data, sizeof(mesh_lmt_optdata) );
  d = *data;
  if ( !d )
    goto failure;
  memset ( d, 0, sizeof(mesh_lmt_optdata) );  /* clear all pointer fields */
                                              /* and the eltypes array */
  if ( !_g2mbl_AssignMeshd ( d, nv, mv, mvhei, mvcp, nhe, mhe,
                             nfac, mfac, mfhei, mkcp ) )
    goto failure;
  if ( !_g2mbl_SetupHessianProfiled ( d, false ) )
    goto failure;
  if ( g2mbl_outputnzdistr )
    _g2mbl_OutputNZDistr ( d );
  if ( !_g2mbl_AllocBFArraysd ( d, nkn1, nkn2 ) )
    goto failure;
  if ( !_g2mbl_SetupElemConstd ( d, dM, dO, C ) )
    goto failure;
  d->nbl = 1;
  d->newpoint = true;
  d->accurate = nkn1 == nkn2;
  d->itn = 0;
  d->func = MYINFINITY;
  d->nu[0] = d->nu[1] = -1.0;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  g2mbl_OptLMTDeallocated ( data );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2mbl_InitBlSurfaceOptLMTd*/

