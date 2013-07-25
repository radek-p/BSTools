
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
boolean _g2mbl_CMPSetupBlockCGPrecondd ( mesh_ml_optdata *d, int bn )
{
  void          *sp;
  int           fnv, cnv;
  int           rmnnz;
  index2        *rmnzi;
  mlblock_desc  *bd;
  int           bnvcp, *bvncpi, bnwcp, *bwncpi;
  int           bnHbl, *cHbl, bnnzH, brmnnz;
  nzHbl_rowdesc *iHbl;
  index2        *bhnzi, *brmnzi, *hri, *rthri;
  int           i, j, k, l, j0, j1;
  unsigned int  *hpermut, *rpermut, *hrpermut;
  int           *hrows, *rrows, *rcols, *hrcols;
  unsigned int  size_hr, nmult_hr, size_rthr, nmult_rthr, hwsize;

  sp = pkv_GetScratchMemTop ();
        /* extract data for the entire mesh */
  fnv    = d->nv;    /* the fine mesh */
  cnv    = d->cnv;    /* number of vertices of the coarse mesh */
  rmnnz  = d->rmnnz;  /* the refinement matrix */
  rmnzi  = d->rmnzi;
        /* extract data for the block */
  bd     = &d->bd[bn];
  bvncpi = bd->vncpi;
  bnvcp  = bd->nvcp;  /* fine Hessian matrix size */
  bnHbl  = bd->nHbl;  /* number of nonzero Hessian blocks */
  iHbl   = bd->iHbl;  /* indices for rows */
  cHbl   = bd->cHbl;  /* indices of columns with nonzero blocks in rows */

  if ( !_g2mbl_CMPFindCoarseMeshBlock ( fnv, cnv, rmnnz, rmnzi, bnvcp, bvncpi,
                                       &bd->nwcp, &bd->wncpi ) )
    goto failure;
  bnwcp = bd->nwcp;
  bwncpi = bd->wncpi;
  if ( !_g2mbl_CMPFindBlockCGPmatrix ( fnv, cnv, rmnnz, rmnzi, bnvcp, bvncpi,
                                      bnwcp, bwncpi, &bd->rmnnz, &bd->rmnzi ) )
    goto failure;
  brmnnz = bd->rmnnz;
  brmnzi = pkv_GetScratchMem ( brmnnz*sizeof(index2) );
  if ( !brmnzi )
    goto failure;
  pkn_SPMindex3to2 ( brmnnz, bd->rmnzi, brmnzi );
        /* find the nonzero Hessian blocks distribution in the sparse representation */
  bnnzH = 2*bnHbl - bnvcp;
  bhnzi = pkv_GetScratchMem ( bnnzH*sizeof(index2) );
  if ( !bhnzi )
    goto failure;
  for ( i = l = 0; i < bnvcp; i++ ) {
    j0 = iHbl[i].firsthbl;
    j1 = j0 + iHbl[i].nhbl - 1;
    for ( j = j0; j < j1; j++ ) {
      k = cHbl[j];
      bhnzi[l  ].i = i;
      bhnzi[l++].j = k;
      bhnzi[l  ].i = k;
      bhnzi[l++].j = i;
    }
    bhnzi[l  ].i = i;
    bhnzi[l++].j = i;
  }
        /* find the distribution of nonzero coefficients of the matrix */
          /* B = H*R, where H is the block Hessian and R is the */
          /* block refinement matrix */
  hpermut = (unsigned int*)pkv_GetScratchMemi ( bnnzH+brmnnz+bnvcp+bnwcp+2 );
  if ( !hpermut )
    goto failure;
  rpermut = &hpermut[bnnzH];
  hrows = (int*)&rpermut[brmnnz];
  rcols = &hrows[bnvcp+1];
  if ( !pkn_SPMCountMMnnzC ( bnvcp, bnvcp, bnwcp,
                             bnnzH, bhnzi, hpermut, hrows, false,
                             brmnnz, brmnzi, rpermut, rcols, false,
                             &size_hr, &nmult_hr ) )
    goto failure;
  bd->nnz1 = size_hr;
  hri = pkv_GetScratchMem ( size_hr*sizeof(index2) );
  if ( !hri )
    goto failure;
  if ( !pkn_SPMmultMMCempty ( bnvcp, bnvcp, bnwcp,
                              bnnzH, bhnzi, hpermut, hrows, true,
                              brmnnz, brmnzi, rpermut, rcols, true,
                              hri ) )
    goto failure;
          /* find the distribution of R^T*B */
  hrpermut = (unsigned int*)pkv_GetScratchMemi ( size_hr );
  if ( !hrpermut )
    goto failure;
  rrows = hrows;
  hrcols = rcols;
  if ( !pkn_SPMCountMTMnnzC ( bnvcp, bnwcp, bnwcp,
                              brmnnz, brmnzi, rpermut, rrows, false,
                              size_hr, hri, hrpermut, hrcols, false /*true?*/,
                              &size_rthr, &nmult_rthr ) )
    goto failure;
  bd->nnz2 = size_rthr;
  rthri = pkv_GetScratchMem ( size_rthr*sizeof(index2) );
  if ( !rthri )
    goto failure;
  if ( !pkn_SPMmultMTMCempty ( bnvcp, bnwcp, bnwcp,
                               brmnnz, brmnzi, rpermut, rrows, true,
                               size_hr, hri, hrpermut, hrcols, true,
                               rthri ) )
    goto failure;

        /* reorder the coarse mesh vertices and prepare the matrix profile */
  PKV_MALLOC ( bd->hprof, 3*bnwcp*sizeof(int) );
  if ( !bd->hprof )
    goto failure;
  if ( !_g2mbl_CMPOrderCoarsePoints ( bnwcp, size_rthr, rthri, brmnnz, bd->rmnzi,
                                      3, &bd->hsize, bd->hprof ) )
    goto failure;

  hwsize = 3*bnwcp;
  PKV_MALLOC ( bd->hrows, 2*(hwsize*sizeof(double*)+bd->hsize*sizeof(double)) );
  if ( !bd->hrows ) {
printf ( "%s\n", ERRMSG_1 );
    goto failure;
  }
  bd->lhrows = &bd->hrows[hwsize];
  pkn_NRBFindRowsd ( hwsize, bd->hprof, (double*)&bd->lhrows[hwsize], bd->hrows );
  pkn_NRBFindRowsd ( hwsize, bd->hprof, &bd->hrows[0][bd->hsize], bd->lhrows );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_CMPSetupBlockCGPrecondd*/

boolean _g2mbl_CMPSetupBlockHessiansd ( mesh_ml_optdata *d )
{
  int          i, nblocks;
  mlblock_desc *bd;

  nblocks = d->nblocks;
  for ( i = 0; 2*i+1 < nblocks; i++ ) {
    bd = &d->bd[i];
    if ( bd->nvcp <= MAX_NVCP ) {
      if ( !_g2mbl_MLSetupBlockCholHessiand ( d, i ) ) {
printf ( "%s\n", ERRMSG_12 );
        return false;
      }
    }
    else {
      if ( !_g2mbl_MLSetupBlockCGHessiand ( d, i ) ) {
printf ( "%s\n", ERRMSG_13 );
        return false;
      }
      if ( d->bd[2*i+1].nvcp > MAX_NVCP ) {
        if ( !_g2mbl_CMPSetupBlockCGPrecondd ( d, i ) )
          continue; /*return false;*/ /* this error is ignored */
      }
    }
  }
  for ( ; i < nblocks; i++ )
    if ( !_g2mbl_MLSetupBlockCholHessiand ( d, i ) ) {
printf ( "%s\n", ERRMSG_12 );
      return false;
    }
  return true;
} /*_g2mbl_CPSetupBlockHessiansd*/

