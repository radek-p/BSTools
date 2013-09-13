
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2011, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <sys/times.h>
#include <unistd.h>

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
#include "msgpool.h"

/* ///////////////////////////////////////////////////////////////////////// */
boolean g2mbl_MLSOptInitd ( int nv, BSMvertex *mv, int *mvhei, point3d *mvcp,
                            int nhe, BSMhalfedge *mhe,
                            int nfac, BSMfacet *mfac, int *mfhei,
                            byte *mkcp,
                            int nkn1, int nkn2, short nlevels, void **data )
{
  void            *sp;
  mesh_ml_optdata *d;
  int             i, nbl, nel;
  mlblock_desc    *bd;
  meshdom_elem    *el;

  sp = pkv_GetScratchMemTop ();
  *data = d = NULL;
#ifdef G2MBL_TIME_IT
  pkv_Tic ( NULL );
#endif
  PKV_MALLOC ( *data, sizeof(mesh_ml_optdata) );
  d = *data;
  if ( !d ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_1, ERRMSG_1 );
    goto failure;
  }
  memset ( d, 0, sizeof(mesh_ml_optdata) );  /* clear all pointer fields */
  if ( !_g2mbl_MLAssignMeshd ( d, nv, mv, mvhei, mvcp, nhe, mhe,
                               nfac, mfac, mfhei, mkcp ) ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_14, ERRMSG_14 );
    goto failure;
  }
  d->nvars = d->nvcp;

printf ( "nvcp = %d, nvars = %d\n", d->nvcp, d->nvars );

  if ( !_g2mbl_MLSetupBlocksd ( d, nlevels, 1, SUBBLOCK_STEP ) )
    goto failure;
  if ( !_g2mbl_MLFindElementsd ( d, 1, nkn1 != nkn2 ) ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_15, ERRMSG_15 );
    goto failure;
  }
  if ( !_g2mbl_MLFindBlockElementsd ( d ) )
    goto failure;
  if ( !_g2mbl_MLSSetupBlockHessiansd ( d ) )
    goto failure;
  if ( !_g2mbl_MLOptAllocBFArraysd ( d, nkn1, nkn2, false ) )
    goto failure;
  if ( !_g2mbl_MLSFindCPNormalsd ( d ) ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_16, ERRMSG_16 );
    goto failure;
  }
  nel = d->ndomelems;
  el = d->domelem;
  for ( i = 0; i < nel; i++ )
    el[i].C = 0.0;
  nbl = d->nblocks;
  bd = d->bd;
  for ( i = 0; i < nbl; i++ )
    bd[i].fghflag = 0;
  d->currentblock = nbl-1;
  d->nextlevel = d->nlevels-1;
  d->dirtyblock = 1;
  d->nu[0] = d->nu[1] = -1.0;
#ifdef G2MBL_TIME_IT
  d->time_prep = pkv_Toc ( NULL );
  d->time_h = d->time_cg = 0;
#endif

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  g2mbl_MLOptDeallocated ( data );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2mbl_MLSOptInitd*/

