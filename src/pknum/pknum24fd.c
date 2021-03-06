
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
boolean pkn_SPMmultMMTCempty ( int nra, int nca, int nrb,
                               unsigned int nnza, index2 *ai,
                               unsigned int *apermut, int *acols, boolean ca,
                               unsigned int nnzb, index2 *bi,
                               unsigned int *bpermut, int *brows, boolean rb,
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
    apermut = (unsigned int*)pkv_GetScratchMemi ( nnza+nca+1 );
    if ( !apermut )
      goto failure;
    ca = false;
    acols = (int*)&apermut[nnza];
  }
  if ( !bpermut ) {
    bpermut = (unsigned int*)pkv_GetScratchMemi ( nnzb+nrb+1 );
    if ( !bpermut )
      goto failure;
    rb = false;
    brows = (int*)&bpermut[nnzb];
  }
  if ( !ca ) {
    if ( !pkn_SPMFindCols ( nra, nca, nnza, ai, apermut, false, acols ) )
      goto failure;
  }
  if ( !rb ) {
    if ( !pkn_SPMFindRows ( nrb, nca, nnzb, bi, bpermut, false, brows ) )
      goto failure;
  }

  sp1 = pkv_GetScratchMemTop ();
  wsps = 0;
  nnz = 0;
  for ( i = 0; i < nrb; i++ ) {
    k0 = brows[i];
    k1 = brows[i+1];
        /* find the workspace size and allocate the workspace */
    wspn = 0;
    for ( k = k0; k < k1; k++ ) {
      j = bi[bpermut[k]].j;
      wspn += acols[j+1]-acols[j];
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
      j = bi[m].j;
      l0 = acols[j];
      l1 = acols[j+1];
      for ( l = l0; l < l1; l++ ) {
        aikbkj[cnt].i = n = apermut[l];
        aikbkj[cnt].j = m;
        wsp[cnt++] = ai[n].i;
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
} /*pkn_SPMmultMMTCempty*/

