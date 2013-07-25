
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pkvaria.h"
#include "pknum.h"
#include "bsmesh.h"
#include "bsmprivate.h"

/* ////////////////////////////////////////////////////////////////////////// */
boolean bsm_FindVertexDistances2 ( int nv, BSMvertex *mv, int *mvhei,
                                   int nhe, BSMhalfedge *mhe,
                                   int nfac, BSMfacet *mfac, int *mfhei,
                                   int v, int *dist )
{
  void      *sp;
  pkv_queue *q;
  int       i, j, k, f, fhe, deg, ffhe, fdeg;
  int       d;
  char      *vtag;

  sp = pkv_GetScratchMemTop ();
/* integer numbers of vertices will be stored in the queue */
  q = pkv_InitQueue ( nv, sizeof(int) );
  if ( !q )
    return false;
  vtag = pkv_GetScratchMem ( nv );
  if ( !vtag )
    goto failure;
  memset ( vtag, 0, nv );
  for ( i = 0; i < nv; i++ )
    dist[i] = nv;  /* initially infinite */
  vtag[v] = 1;
  dist[v] = 0;
  pkv_QueueInsert ( q, &v );
  do {
    pkv_QueueRemoveFirst ( q, &i );
    d = dist[i]+1;
    deg = mv[i].degree;
    fhe = mv[i].firsthalfedge;
    for ( j = 0; j < deg; j++ ) {
      f = mhe[mvhei[fhe+j]].facetnum;
      fdeg = mfac[f].degree;
      ffhe = mfac[f].firsthalfedge;
      for ( k = 0; k < fdeg; k++ ) {
        v = mhe[mfhei[ffhe+k]].v0;
        if ( !vtag[v] ) {
          vtag[v] = 1;
          dist[v] = d;
          pkv_QueueInsert ( q, &v );
        }
      }
    }
  } while ( !pkv_QueueEmpty ( q ) );
  PKV_FREE ( q );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  PKV_FREE ( q );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*bsm_FindVertexDistances2*/

