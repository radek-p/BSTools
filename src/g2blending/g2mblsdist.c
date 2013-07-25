
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2012                            */
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

/* ////////////////////////////////////////////////////////////////////////// */
typedef struct {
    int vertnum, dist;
  } queue_elem;

/* ////////////////////////////////////////////////////////////////////////// */
boolean g2mbl_FindDistances ( int nv, BSMvertex *mv, int *mvhei,
                              int nhe, BSMhalfedge *mhe,
                              int nfac, BSMfacet *mfac, int *mfhei,
                              int *fdist )
{
  void       *sp;
  pkv_queue  *q;
  queue_elem qe;
  int        i, j, vfhe, vd, ffhe, fd;
  int        vn, dist, fn, vn1;
  char       *vtag, *ftag;

  sp = pkv_GetScratchMemTop ();
  q = pkv_InitQueue ( nv, sizeof(queue_elem) );
  vtag = pkv_GetScratchMem ( nv );
  ftag = pkv_GetScratchMem ( nfac );
  if ( !q || !vtag || !ftag )
    goto failure;

  memset ( ftag, 0, nfac );
  for ( i = 0; i < nfac; i++ )
    fdist[i] = INT_MAX;    /* stands for infinity */
  for ( i = 0; i < nv; i++ ) {
    vd = mv[i].degree;
    vfhe = mv[i].firsthalfedge;
    if ( vd != 4 && mhe[mvhei[vfhe+vd-1]].otherhalf != -1 ) {
      vtag[i] = 1;
      qe.vertnum = i;
      qe.dist = 0;
      pkv_QueueInsert ( q, &qe );
    }
    else
      vtag[i] = 0;
  }
  while ( !pkv_QueueEmpty ( q ) ) {
    pkv_QueueRemoveFirst ( q, &qe );
    vn = qe.vertnum;
    dist = qe.dist;
    vd = mv[vn].degree;
    vfhe = mv[vn].firsthalfedge;
    for ( i = 0; i < vd; i++ ) {
      fn = mhe[mvhei[vfhe+i]].facetnum;
      ftag[fn] = 1;
      fdist[fn] = min ( fdist[fn], dist );
      fd = mfac[fn].degree;
      ffhe = mfac[fn].firsthalfedge;
      for ( j = 0; j < fd; j++ ) {
        vn1 = mhe[mfhei[ffhe+j]].v0;
        if ( !vtag[vn1] ) {
          vtag[vn1] = 1;
          qe.vertnum = vn1;
          qe.dist = dist+1;
          pkv_QueueInsert ( q, &qe );
        }
      }
    }
  }

  pkv_SetScratchMemTop ( sp );
  if ( q ) PKV_FREE ( q );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2mbl_FindDistances*/

