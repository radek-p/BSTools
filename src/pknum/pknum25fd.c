
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"

/* /////////////////////////////////////////////////////////////////////////// */
boolean pkn_SPMmultMTMCempty ( int nra, int nca, int ncb,
                               unsigned int nnza, index2 *ai,
                               unsigned int *apermut, int *arows, boolean ra,
                               unsigned int nnzb, index2 *bi,
                               unsigned int *bpermut, int *bcols, boolean cb,
                               index2 *abi )
{
  void         *sp, *sp1;
  unsigned int nnz;
  int          *wsp, wsps, wspn;
  unsigned int *auxpermut;
  index2       *aikbkj;
  int          i, j, k, k0, k1, l, l0, l1, cnt, m, n;

  sp = pkv_GetScratchMemTop ();
  wsp = NULL;
  aikbkj = NULL;
  if ( !apermut ) {
    apermut = (unsigned int*)pkv_GetScratchMemi ( nnza+nra+1 );
    if ( !apermut )
      goto failure;
    ra = false;
    arows = (int*)&apermut[nnza];
  }
  if ( !bpermut ) {
    bpermut = (unsigned int*)pkv_GetScratchMemi ( nnzb+ncb+1 );
    if ( !bpermut )
      goto failure;
    cb = false;
    bcols = (int*)&bpermut[nnzb];
  }
  if ( !ra ) {
    if ( !pkn_SPMFindRows ( nra, nca, nnza, ai, apermut, false, arows ) )
      goto failure;
  }
  if ( !cb ) {
    if ( !pkn_SPMFindCols ( nca, ncb, nnzb, bi, bpermut, false, bcols ) )
      goto failure;
  }

  sp1 = pkv_GetScratchMemTop ();
  wsps = 0;
  nnz = 0;
  for ( i = 0; i < ncb; i++ ) {
    k0 = bcols[i];
    k1 = bcols[i+1];
        /* find the workspace size and allocate the workspace */
    wspn = 0;
    for ( k = k0; k < k1; k++ ) {
      j = bi[bpermut[k]].i;
      wspn += arows[j+1]-arows[j];
    }
    if ( wspn > wsps ) {
      pkv_SetScratchMemTop ( sp1 );
      wsp = pkv_GetScratchMemi ( 4*wspn );
      if ( !wsp )
        goto failure;
      aikbkj = (index2*)&wsp[wspn];
      wsps = wspn;
    }
        /* find the positions (row numbers) of nonzero product coefficients */
    cnt = 0;
    for ( k = k0; k < k1; k++ ) {
      m = bpermut[k];
      j = bi[m].i;
      l0 = arows[j];
      l1 = arows[j+1];
      for ( l = l0; l < l1; l++ ) {
        aikbkj[cnt].i = n = apermut[l];
        aikbkj[cnt].j = m;
        wsp[cnt++] = ai[n].j;
      }
    }
    if ( cnt ) {
          /* set up the identity permutation */
      auxpermut = (unsigned int*)&aikbkj[wspn];
      for ( l = 0; l < cnt; l++ )
        auxpermut[l] = l;
          /* sort by columns */
      if ( pkv_SortKernel ( sizeof(int), ID_UNSIGNED, sizeof(int),
                            0, cnt, wsp, auxpermut ) != SORT_OK )
        goto failure;
      k = auxpermut[0];
      l0 = wsp[k];
      abi[nnz  ].j = i;
      abi[nnz++].i = l0;
      for ( l = 1; l < cnt; l++ ) {
        k = auxpermut[l];
        if (wsp[k] != l0 ) {
          l0 = wsp[k];
          abi[nnz  ].j = i;
          abi[nnz++].i = l0;
        }
      }
    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_SPMmultMTMCempty*/

