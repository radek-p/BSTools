
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

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g2mbl_CMPSetupCGPrecondd ( mesh_lmt_optdata *d )
{
  void          *sp;
  int           rmnnz, bnvcp, *bvncpi, bnwcp, *bwncpi;
  int           bnHbl, *cHbl, bnnzH, brmnnz;
  nzHbl_rowdesc *iHbl;
  index2        *rmnzi, *bhnzi, *brmnzi, *hri, *rthri;
  int           i, j, k, l, j0, j1;
  unsigned int  *hpermut, *rpermut, *hrpermut;
  int           *hrows, *rrows, *rcols, *hrcols;
  unsigned int  size_hr, nmult_hr, size_rthr, nmult_rthr, hwsize;

  sp = pkv_GetScratchMemTop ();
  rmnnz  = d->rmnnz;
  rmnzi  = d->rmnzi;
  bvncpi = d->vncpi;
  bnvcp  = d->nvcp;
  bnHbl  = d->Hblsize;
  iHbl   = d->iHbl;
  cHbl   = d->cHbl;

  if ( !_g2mbl_CMPFindCoarseMeshBlock ( d->nv, d->cnv, rmnnz, rmnzi,
                                bnvcp, bvncpi, &d->nwcp, &d->wncpi ) )
    goto failure;
  bnwcp = d->nwcp;
  bwncpi = d->wncpi;
  if ( !_g2mbl_CMPFindBlockCGPmatrix ( d->nv, d->cnv, rmnnz, rmnzi,
              d->nvcp, d->vncpi, bnwcp, bwncpi, &d->pmnnz, &d->pmnzi ) )
    goto failure;
  brmnnz = d->pmnnz;
  brmnzi = pkv_GetScratchMem ( brmnnz*sizeof(index2) );
  if ( !brmnzi )
    goto failure;
  pkn_SPMindex3to2 ( brmnnz, d->pmnzi, brmnzi );
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
  d->nnz1 = size_hr;
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
  d->nnz2 = size_rthr;
  rthri = pkv_GetScratchMem ( size_rthr*sizeof(index2) );
  if ( !rthri )
    goto failure;
  if ( !pkn_SPMmultMTMCempty ( bnvcp, bnwcp, bnwcp,
                               brmnnz, brmnzi, rpermut, rrows, true,
                               size_hr, hri, hrpermut, hrcols, true,
                               rthri ) )
    goto failure;

        /* reorder the coarse mesh vertices and prepare the matrix profile */
  PKV_MALLOC ( d->phprof, 3*bnwcp*sizeof(int) );
  if ( !d->phprof )
    goto failure;
  if ( !_g2mbl_CMPOrderCoarsePoints ( d->nwcp, size_rthr, rthri,
                              d->pmnnz, d->pmnzi, 3, &d->phsize, d->phprof ) )
    goto failure;
  hwsize = 3*d->nwcp;
  PKV_MALLOC ( d->phrows, hwsize*sizeof(double*) + d->phsize*sizeof(double) );
  if ( !d->phrows ) {
printf ( "%s\n", ERRMSG_9 );
    goto failure;
  }
  pkn_NRBFindRowsd ( hwsize, d->phprof, (double*)&d->phrows[hwsize], d->phrows );

#ifdef _DEBUG
printf ( "coarse mesh block: ncp = %d, nvars = %d, hsize = %d\n",
         d->nwcp, 3*d->nwcp, d->phsize );
printf ( "nnz1 = %d, nnz2 = %d\n", d->nnz1, d->nnz2 );
#endif

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_CMPSetupCGPrecondd*/

boolean g2mbl_InitBlCMPSurfaceOptd ( /* fine mesh */ 
                    int fnv, BSMvertex *fmv, int *fmvhei, point3d *fmvcp,
                    int fnhe, BSMhalfedge *fmhe,
                    int fnfac, BSMfacet *fmfac, int *fmfhei,
                    byte *fmkcp,
                           /* number of vertices of the coarse mesh */
                    int cnv,
                           /* refinement matrix */
                    int rmnnz, index2 *rmnzi, double *rmnzc,
                           /* optimization parameters */
                    double C, double dO, double dM,
                    int nkn1, int nkn2, int nbl,
                           /* created data structure */
                    void **data )
{
  void             *sp;
  mesh_lmt_optdata *d;
  int              blseed[G2MBL_MAX_BLOCKS];

  sp = pkv_GetScratchMemTop ();
  *data = d = NULL;
  PKV_MALLOC ( *data, sizeof(mesh_lmt_optdata) );
  d = *data;
  if ( !d )
    goto failure;
  memset ( d, 0, sizeof(mesh_lmt_optdata) );  /* clear all pointer fields */
                                              /* and the eltypes array */
  if ( !_g2mbl_AssignMeshd ( d, fnv, fmv, fmvhei, fmvcp, fnhe, fmhe,
                             fnfac, fmfac, fmfhei, fmkcp ) )
    goto failure;
        /* record the coarse mesh preconditioner */
  d->cnv   = cnv;
  d->rmnnz = rmnnz;
  d->rmnzi = rmnzi;
  d->rmnzc = rmnzc;
  if ( !_g2mbl_SetupHessianProfiled ( d, true ) )
    goto failure;
  if ( g2mbl_outputnzdistr )
    _g2mbl_OutputNZDistr ( d );
  if ( !_g2mbl_SetupHbl3x3d ( d ) )
    goto failure;
  if ( !_g2mbl_AllocBFArraysd ( d, nkn1, nkn2 ) )
    goto failure;
  if ( !_g2mbl_SetupElemConstd ( d, dM, dO, C ) )
    goto failure;
  PKV_MALLOC ( d->bltag, fnv*sizeof(int) );
  if ( !d->bltag )
    goto failure;
  if ( !_g2mbl_SetupAltBlocks ( fnv, fmv, fmvhei, fnhe, fmhe, fnfac, fmfac, fmfhei,
                                d->nncpi, nbl, d->bltag, blseed ) )
    goto failure;
  if ( !_g2mbl_SetupAltBlockDescription ( d, nbl ) )
    goto failure;
  if ( !_g2mbl_SetupBlockRowsd ( d ) )
    goto failure;
  if ( !_g2mbl_CMPSetupCGPrecondd ( d ) )
    goto failure;
                                            
  d->nbl = nbl;
  d->newpoint = true;
  if ( (d->accurate = nkn1 == nkn2) )
    d->ibl = 0;
  else
    d->ibl = nbl;
  d->itn = 0;
  d->next_entire = nbl;
  d->func = MYINFINITY;
  d->nu[0] = d->nu[1] = -1.0;
  d->nextblock = 0;
  d->last_ntn = false;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  g2mbl_OptLMTDeallocated ( data );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2mbl_InitBlCMPSurfaceOptd*/

