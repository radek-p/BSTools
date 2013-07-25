
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2011                                  */
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
boolean _g2mbl_MLOptAllocBFArraysd ( mesh_ml_optdata *d, int nkn1, int nkn2,
                                     boolean reparam )
{
  int i, j;
  int size;
  double *aqcoeff, *aqknots, *bqcoeff, *bqknots;
  double *abf, *adbf, *addbf, *adddbf, *bbf, *bdbf, *bddbf, *bdddbf;

        /* conpute the size and allocate arrays for basis functions */
  if ( nkn2 != nkn1 ) {
    size = nkn1+nkn2;
    for ( i = 3, j = 0;  i <= GH_MAX_K;  i++, j++ )
      if ( d->eltypes[j] )
        size += g2mbl_NiSize ( nkn1, i ) + g2mbl_NijSize ( nkn1, i ) +
                g2mbl_MijSize ( nkn1, i ) + g2mbl_NiSize ( nkn2, i ) +
                nkn1*nkn1 + nkn2*nkn2;
  }
  else {
    size = nkn1;
    for ( i = 3, j = 0;  i <= GH_MAX_K;  i++, j++ )
      if ( d->eltypes[j] )
        size += g2mbl_NiSize ( nkn1, i ) + g2mbl_NijSize ( nkn1, i ) +
                g2mbl_MijSize ( nkn1, i ) + nkn1*nkn1;
  }
#ifdef _DEBUG
printf ( "Nsize = %d\n", size );
#endif
  PKV_MALLOC ( d->aqcoeff, size*sizeof(double) );
  if ( !d->aqcoeff )
    goto failure;
  if ( nkn2 != nkn1 )
    d->bqcoeff = &d->aqcoeff[nkn1];
  else
    d->bqcoeff = d->aqcoeff;
  bqcoeff = &d->bqcoeff[nkn2];
  for ( i = 3, j = 0;  i <= GH_MAX_K;  i++, j++ )
    if ( d->eltypes[j] ) {
      d->aJac[j] = bqcoeff;
      d->aNitabs[j] = d->aJac[j] + nkn1*nkn1;
      d->aNijtabs[j] = d->aNitabs[j] + g2mbl_NiSize ( nkn1, i );
      d->aMijtabs[j] = d->aNijtabs[j] + g2mbl_NijSize ( nkn1, i );
      if ( nkn2 != nkn1 ) {
        d->bJac[j] = d->aMijtabs[j] + g2mbl_MijSize ( nkn1, i );
        d->bNitabs[j] = d->bJac[j] + nkn2*nkn2;
        bqcoeff = d->bNitabs[j] + g2mbl_NiSize ( nkn2, i );
      }
      else {
        d->bJac[j] = d->aJac[j];
        d->bNitabs[j] = d->aNitabs[j];
        bqcoeff = d->aMijtabs[j] + g2mbl_MijSize ( nkn1, i );
      }
    }
        /* evaluate the basis functions for regular squares */
  if ( !_g2bl_TabBasisFuncd ( nkn1, &aqknots, &aqcoeff,
                              &abf, &adbf, &addbf, &adddbf ) )
    goto failure;
  memcpy ( d->aqcoeff, aqcoeff, nkn1*sizeof(double) );
  g2bl_TabNid ( nkn1, abf, adbf, addbf, adddbf, d->aNitabs[1] );
  g2bl_TabNijd ( nkn1, abf, adbf, addbf, d->aNijtabs[1] );
  g2bl_TabMijd ( nkn1, abf, adbf, addbf, adddbf, d->aMijtabs[1] );
  for ( i = 0; i < nkn1*nkn1; i++ )
    d->aJac[1][i] = 1.0;
  if ( nkn2 != nkn1 ) {
    if ( !_g2bl_TabBasisFuncd ( nkn2, &bqknots, &bqcoeff,
                                &bbf, &bdbf, &bddbf, &bdddbf ) )
      goto failure;
    memcpy ( d->bqcoeff, bqcoeff, nkn2*sizeof(double) );
    g2bl_TabNid ( nkn2, bbf, bdbf, bddbf, bdddbf, d->bNitabs[1] );
    for ( i = 0; i < nkn2*nkn2; i++ )
      d->bJac[1][i] = 1.0;
  }
  d->nkn1 = nkn1;
  d->nkn2 = nkn2;
        /* evaluate the basis functions for special elements */
  if ( d->eltypes[0] ) {
    if ( !g2mbl_SetupHolePatchMatrixd ( 3 ) )
      goto failure;
    if ( !g2mbl_TabNid ( 3, nkn1, aqknots, d->aNitabs[0], d->aJac[0],
                         reparam ) )
      goto failure;
    g2mbl_TabNijd ( 3*6+1, nkn1, d->aNitabs[0], d->aNijtabs[0] );
    g2mbl_TabMijd ( 3*6+1, nkn1, d->aNitabs[0], d->aMijtabs[0] );
    if ( nkn2 != nkn1 ) {
      if ( !g2mbl_TabNid ( 3, nkn2, bqknots, d->bNitabs[0], d->bJac[0],
                           reparam ) )
        goto failure;
    }
  }
  for ( i = 5; i <= GH_MAX_K; i++ )
    if ( d->eltypes[i-3] ) {
      if ( !g2mbl_SetupHolePatchMatrixd ( i ) )
        goto failure;
      if ( !g2mbl_TabNid ( i, nkn1, aqknots, d->aNitabs[i-3], d->aJac[i-3],
                           reparam ) )
        goto failure;
      g2mbl_TabNijd ( i*6+1, nkn1, d->aNitabs[i-3], d->aNijtabs[i-3] );
      g2mbl_TabMijd ( i*6+1, nkn1, d->aNitabs[i-3], d->aMijtabs[i-3] );
      if ( nkn2 != nkn1 ) {
        if ( !g2mbl_TabNid ( i, nkn2, bqknots, d->bNitabs[i-3], d->bJac[i-3],
                             reparam ) )
          goto failure;
      }
    }
  return true;

failure:
  return false;
} /*_g2mbl_MLOptAllocBFArraysd*/

