
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
#include "msgpool.h"

/*#define __DEBUG*/
/* ///////////////////////////////////////////////////////////////////////// */
boolean g2mbl_MLOptInitd ( int nv, BSMvertex *mv, int *mvhei, point3d *mvcp,
                           int nhe, BSMhalfedge *mhe,
                           int nfac, BSMfacet *mfac, int *mfhei,
                           byte *mkcp,
                           double C, double dO, double dM,
                           int nkn1, int nkn2, short nlevels, void **data )
{
  void            *sp;
  mesh_ml_optdata *d;
  int             i, nbl;
  mlblock_desc    *bd;

  sp = pkv_GetScratchMemTop ();
  *data = d = NULL;
#ifdef G2MBL_TIME_IT
  pkv_Tic ( NULL );
#endif
  PKV_MALLOC ( *data, sizeof(mesh_ml_optdata) );
  d = *data;
  if ( !d ) {
printf ( "%s\n", ERRMSG_1 );
    goto failure;
  }
  memset ( d, 0, sizeof(mesh_ml_optdata) );  /* clear all pointer fields */
  if ( !_g2mbl_MLAssignMeshd ( d, nv, mv, mvhei, mvcp, nhe, mhe,
                               nfac, mfac, mfhei, mkcp ) ) {
printf ( "%s\n", ERRMSG_14 );
    goto failure;
  }
  d->nvars = 3*d->nvcp;

printf ( "nvcp = %d, nvars = %d\n", d->nvcp, d->nvars );

  if ( !_g2mbl_MLSetupBlocksd ( d, nlevels, 3, SUBBLOCK_STEP ) )
    goto failure;
  if ( !_g2mbl_MLFindElementsd ( d, 3, nkn1 != nkn2 ) ) {
printf ( "%s\n", ERRMSG_15 );
    goto failure;
  }
  if ( !_g2mbl_MLFindBlockElementsd ( d ) )
    goto failure;
  if ( !_g2mbl_MLSetupBlockHessiansd ( d ) )
    goto failure;
  if ( !_g2mbl_MLOptAllocBFArraysd ( d, nkn1, nkn2, true ) )
    goto failure;
  if ( !_g2mbl_MLSetupElemConstd ( d, dM, dO, C ) )
    goto failure;
  nbl = d->nblocks;
  bd = d->bd;
  for ( i = 0; i < nbl; i++ )
    bd[i].fghflag = 0;
  d->currentblock = d->nblocks-1;
  d->nextlevel = d->nlevels-1;
  d->dirtyblock = 0;
  d->nu[0] = d->nu[1] = -1.0;
#ifdef G2MBL_TIME_IT
  d->time_prep = pkv_Toc ( NULL );
  d->time_h = d->time_cg = 0;
#endif

#ifdef __DEBUG
for ( i = 0; i < nbl; i++ ) {
  printf ( "block %2d, nvcp = %2d, vdomel = %2d",
           i, bd[i].nvcp, bd[i].ndomel );
  if ( bd[i].iHbl )
    printf ( "*" );
  printf ( "\n" );
}
#endif
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  g2mbl_MLOptDeallocated ( data );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2mbl_MLOptInitd*/

