
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
boolean _g2mbl_MLSetupBlocksd ( mesh_ml_optdata *d,
                                short nlevels, short bsize, int step )
{
  void         *sp;
  short        nblocks;
  mlblock_desc *bd;
  int          nv, nhe, nfac, *mvhei, *mfhei;
  BSMvertex    *mv;
  BSMhalfedge  *mhe;
  BSMfacet     *mfac;
  int          nvcp, *vncpi;
  int          i, k, b, b1, b2;
  char         *mvtag;

  sp = pkv_GetScratchMemTop ();
  if ( nlevels < 1 || nlevels > G2MBL_MAX_LEVELS )
    goto failure;
  d->nlevels = nlevels;
  d->nblocks = nblocks = (0x0001 << nlevels) - 1;

printf ( "nlevels = %d, nblocks = %d\n", nlevels, nblocks );

  PKV_MALLOC ( d->bd, nblocks*sizeof(mlblock_desc) );
  if ( !(bd = d->bd) )
    goto failure;
        /* clear all pointer fields */
  memset ( bd, 0, nblocks*sizeof(mlblock_desc) );
        /* get the mesh */
  nv    = d->nv;
  mv    = d->mv;
  mvhei = d->mvhei;
  nhe   = d->nhe;
  mhe   = d->mhe;
  nfac  = d->nfac;
  mfac  = d->mfac;
  mfhei = d->mfhei;
  mvtag = d->mvtag;
        /* setup the first block, for the entire optimization job */
  bd[0].nvcp = nvcp = d->nvcp;
  bd[0].nvars = bsize*nvcp;
  PKV_MALLOC ( bd[0].vncpi, nvcp*sizeof(int) );
  if ( !(vncpi = bd[0].vncpi) )
    goto failure;
        /* it is important that the vertex indices are stored */
        /* in ascending order */
  for ( i = k = 0;  i < nv;  i++ )
    if ( !mvtag[i] )
      vncpi[k++] = i;
        /* now the recursive block subdivision */
  for ( b = 0; b+b < nblocks-2; b++ ) {
          /* numbers of the block parts */
/*
  printf ( "b%3d: ", b );
*/
    nvcp = bd[b].nvcp;
    vncpi = bd[b].vncpi;
    b1 = b+b+1;
    b2 = b1+1;
    if ( !_g2mbl_MLDivideBlock ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                 nvcp, vncpi,
                 &bd[b1].nvcp, &bd[b1].vncpi, &bd[b1].seed,
                 &bd[b2].nvcp, &bd[b2].vncpi, &bd[b2].seed, step ) )
      goto failure;
        /* make sure that the vertex indices are in ascending order */
    if ( pkv_SortFast ( sizeof(int), ID_UNSIGNED, sizeof(int),
                        0, bd[b1].nvcp, bd[b1].vncpi ) != SORT_OK )
      goto failure;
    if ( pkv_SortFast ( sizeof(int), ID_UNSIGNED, sizeof(int),
                        0, bd[b2].nvcp, bd[b2].vncpi ) != SORT_OK )
      goto failure;
    bd[b1].nvars = bsize*bd[b1].nvcp;
    bd[b2].nvars = bsize*bd[b2].nvcp;
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_MLSetupBlocksd*/

