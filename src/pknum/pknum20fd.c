
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
boolean pkn_SPsubMCountMMTnnzC ( int nra, int nca, int nrb,
                                 unsigned int nnza, index3 *ai,
                                 unsigned int *apermut, int *acols, boolean ca,
                                 unsigned int nnzb, index3 *bi,
                                 unsigned int *bpermut, int *brows, boolean rb,
                                 unsigned int *nnzab, unsigned int *nmultab )
{
  void         *sp;
  unsigned int nnz, nmult;
  int          i, j, k, k0, k1, l, l0, l1, cnt;
  int          *wsp, wsps, wspn;

  sp = pkv_GetScratchMemTop ();
  if ( !ca ) {
    if ( !pkn_SPsubMFindCols ( nra, nca, nnza, ai, apermut, ca, acols ) )
      goto failure;
  }
  if ( !rb ) {
    if ( !pkn_SPsubMFindRows ( nca, nrb, nnzb, bi, bpermut, rb, brows ) )
      goto failure;
  }

  wsp = NULL;
  wsps = 0;
  nnz = nmult = 0;
  for ( i = 0; i < nrb; i++ ) {
    k0 = brows[i];
    k1 = brows[i+1];
        /* find the workspace size and allocate the workspace */
    wspn = 0;
    for ( k = k0; k < k1; k++ ) {
      j = bi[bpermut[k]].j;
      wspn += acols[j+1]-acols[j];
    }
    nmult += wspn;
    if ( wspn > wsps ) {
      pkv_SetScratchMemTop ( sp );
      wsp = pkv_GetScratchMemi ( wspn );
      if ( !wsp )
        goto failure;
      wsps = wspn;
    }
        /* find the positions (row numbers) of nonzero product coefficients */
        /* in the i-th column */
    cnt = 0;
    for ( k = k0; k < k1; k++ ) {
      j = bi[bpermut[k]].j;
      l0 = acols[j];
      l1 = acols[j+1];
      for ( l = l0; l < l1; l++ )
        wsp[cnt++] = ai[apermut[l]].i;
    }
        /* count the different positions */
    if ( cnt ) {
      if ( pkv_SortFast ( sizeof(int), ID_UNSIGNED, sizeof(int),
                          0, cnt, wsp ) != SORT_OK )
        goto failure;
      l0 = wsp[0];
      nnz ++;
      for ( l = 1; l < cnt; l++ )
        if ( wsp[l] != l0 ) {
          l0 = wsp[l];
          nnz ++;
        }
    }
  }

  *nnzab = nnz;
  if (nmultab)
    *nmultab = nmult;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_SPsubMCountMMTnnzC*/

boolean pkn_SPsubMFindMMTnnzC ( int nra, int nca, int nrb,
                                unsigned int nnza, index3 *ai,
                                unsigned int *apermut, int *acols,
                                unsigned int nnzb, index3 *bi,
                                unsigned int *bpermut, int *brows,
                                index2 *abi, int *abpos, index2 *aikbkj )
{
  void         *sp;
  unsigned int nnz, nmult;
  int          i, j, k, k0, k1, l, l0, l1, m, n, cnt, nn;
  int          *wsp, wsps, wspn;
  unsigned int *auxpermut;

  sp = pkv_GetScratchMemTop ();
  wsp = NULL;
  wsps = 0;
  nnz = nmult = nn = 0;
  abpos[0] = 0;
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
      pkv_SetScratchMemTop ( sp );
      wsp = pkv_GetScratchMemi ( wspn );
      if ( !wsp )
        goto failure;
      wsps = wspn;
    }
        /* find the positions (row numbers) of nonzero product coefficients */
        /* in the i-th column */
    cnt = 0;
    for ( k = k0; k < k1; k++ ) {
      m = bpermut[k];
      j = bi[m].j;
      l0 = acols[j];
      l1 = acols[j+1];
      for ( l = l0; l < l1; l++ ) {
        n = apermut[l];
        wsp[cnt++] = ai[n].i;
        aikbkj[nmult].j = bi[m].k;
        aikbkj[nmult++].i = ai[n].k;
      }
    }
        /* count the different positions */
    if ( cnt ) {
          /* set up the identity permutation */
      auxpermut = (unsigned int*)pkv_GetScratchMemi ( cnt );
      if ( !auxpermut )
        goto failure;
      for ( l = 0; l < cnt; l++ )
        auxpermut[l] = l;
          /* sort by columns */
      if ( pkv_SortKernel ( sizeof(int), ID_UNSIGNED, sizeof(int),
                            0, cnt, wsp, auxpermut ) != SORT_OK )
        goto failure;
      k1 = nn;
      l0 = wsp[auxpermut[0]];
      abpos[nnz] = nn ++;
      abi[nnz].j = i;
      abi[nnz++].i = l0;
      for ( l = 1; l < cnt; l++ ) {
        k0 = auxpermut[l];
        if ( wsp[k0] != l0 ) {
          l0 = wsp[k0];
          abpos[nnz] = nn;
          abi[nnz].j = i;
          abi[nnz++].i = l0;
        }
        nn ++;
      }
      pkv_SortPermute ( sizeof(index2), cnt, &aikbkj[k1], auxpermut );
    }
  }
  abpos[nnz] = nn;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_SPsubMFindMMTnnzC*/

