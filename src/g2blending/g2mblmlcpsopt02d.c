
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012, 2013                            */
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

/* ///////////////////////////////////////////////////////////////////////// */
#ifdef __DEBUG
void __g2mbl_WriteOutBlocks ( mesh_ml_optdata *d )
{
  FILE         *f;
  mlblock_desc *bd;
  int          i, nbl;

  f = fopen ( "blocks.txt", "w+" );
  nbl = d->nblocks;
  for ( i = 0; i < nbl; i++ ) {
    bd = &d->bd[i];
    fprintf ( f, "%2d:", i );
    fprintf ( f, " nvcp = %d, nvars = %d, ndomel = %d\n",
              bd->nvcp, bd->nvars, bd->ndomel );
    if ( bd->iHbl ) {
      fprintf ( f, "    nwcp = %d, rmnnz = %d,", bd->nwcp, bd->rmnnz );
      fprintf ( f, " nnz1 = %d, nnz2 = %d\n", bd->nnz1, bd->nnz2 );
      fprintf ( f, "\n" );
    }
    fprintf ( f, "\n" );
  }
  fclose ( f );
} /*__g2mbl_WriteOutBlocks*/
#endif

/* ///////////////////////////////////////////////////////////////////////// */
boolean g2mbl_MLSCMPOptInitd ( /* fine mesh */
                    int fnv, BSMvertex *fmv, int *fmvhei, point3d *fmvcp,
                    int fnhe, BSMhalfedge *fmhe,
                    int fnfac, BSMfacet *fmfac, int *fmhei,
                    byte *fmkcp,
                           /* coarse mesh */
                    int cnv,
                           /* refinement matrix */
                    int rmnnz, index2 *rmnzi, double *rmnzc,
                           /* optimization parameters */
                    int nkn1, int nkn2, short nlevels,
                           /* created data structure */
                    void **data )
{
  void            *sp;
  mesh_ml_optdata *d;
  int             i, nbl, nel;
  mlblock_desc    *bd;
  meshdom_elem    *el;

  sp = pkv_GetScratchMemTop ();
  *data = d = NULL;
  if ( nlevels < 1 )
    goto failure;
#ifdef G2MBL_TIME_IT
  pkv_Tic ( NULL );
#endif
  PKV_MALLOC ( *data, sizeof(mesh_ml_optdata) );
  d = *data;
  if ( !d ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_14, ERRMSG_14 );
    goto failure;
  }
  memset ( d, 0, sizeof(mesh_ml_optdata) );
  if ( !_g2mbl_CMPAssignMeshd ( d, fnv, fmv, fmvhei, fmvcp, fnhe, fmhe,
                                fnfac, fmfac, fmhei, fmkcp,
                                cnv, rmnnz, rmnzi, rmnzc ) ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_15, ERRMSG_15 );
    goto failure;
  }
  d->nvars = d->nvcp;

printf ( "nvcp = %d, nvars = %d\n", d->nvcp, d->nvars );

  if ( !_g2mbl_MLSetupBlocksd ( d, nlevels, 1, SUBBLOCK_STEP_CP ) )
    goto failure;
  if ( !_g2mbl_MLFindElementsd ( d, 1, nkn1 != nkn2 ) ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_15, ERRMSG_15 );
    goto failure;
  }
  if ( !_g2mbl_MLFindBlockElementsd ( d ) )
    goto failure;
  if ( !_g2mbl_CMPSSetupBlockHessiansd ( d ) )
    goto failure;
  if ( !_g2mbl_MLOptAllocBFArraysd ( d, nkn1, nkn2, true ) )
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
  d->currentblock = d->nblocks-1;
  d->nextlevel = d->nlevels-1;
  d->dirtyblock = 0;
  d->nu[0] = d->nu[1] = -1.0;
#ifdef G2MBL_TIME_IT
  d->time_prep = pkv_Toc ( NULL );
  d->time_h = d->time_cg = 0;
#endif
#ifdef __DEBUG
__g2mbl_WriteOutBlocks ( d );
#endif
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  g2mbl_MLOptDeallocated ( data );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2mbl_MLSCMPOptInitd*/

