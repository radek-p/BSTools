
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2014                            */
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

boolean bsm_CheckMeshIntegrity ( int nv, const BSMvertex *mv, const int *mvhei,
                                 int nhe, const BSMhalfedge *mhe,
                                 int nfac, const BSMfacet *mfac, const int *mfhei,
                                 int *errcode, int *errelem )
{
#define FAILURE(c,n) \
  { error_code = c;  error_el = n;  goto failure; }
#define PREVVERT_HEDGE(vn,en) \
  ((int)(en) > 0 ? \
    mvhei[mv[vn].firsthalfedge+(int)(en)-1] : \
    (mhe[mvhei[mv[vn].firsthalfedge+mv[vn].degree-1]].otherhalf >= 0 ? \
      mvhei[mv[vn].firsthalfedge+mv[vn].degree-1] : \
      -1))
#define NEXTFAC_HEDGE(fn,en) \
  ((int)(en) < mfac[fn].degree-1 ? \
    mfhei[mfac[fn].firsthalfedge+(int)(en)+1] : \
    mfhei[mfac[fn].firsthalfedge])

  void *sp;
  int  i, j, k, l, d;
  int  error_code, error_el;
  char *iset, *wlv, *wlf;
  int  he0, he1, he2;

  sp = pkv_GetScratchMemTop ();
  if ( nv < 3 || nhe < 3 || nfac < 1 )
    FAILURE ( BSM_ERROR_NOT_ENOUGH_ELEMENTS, -1 );

        /* the arrays mvhei and mfhei have to hold permutations of */
        /* integers from 0 to nhe-1; verify it */
  iset = pkv_GetScratchMem ( nhe );
  if ( !iset )
    FAILURE ( BSM_ERROR_NO_SCRATCH_MEMORY, -1 );
  memset ( iset, 0, nhe );
  for ( i = 0; i < nhe; i++ ) {
    if ( mvhei[i] < 0 || mvhei[i] >= nhe )
      FAILURE ( BSM_ERROR_INVALID_MVHEI, i )
    else
      iset[mvhei[i]] |= 0x01;
    if ( mfhei[i] < 0 || mfhei[i] >= nhe )
      FAILURE ( BSM_ERROR_INVALID_MFHEI, i )
    else
      iset[mfhei[i]] |= 0x02;
  }
  for ( i = 0; i < nhe; i++ )
    if ( iset[i] != 0x03 )
      FAILURE ( BSM_ERROR_INVALID_PERMUTATION, i )
  pkv_SetScratchMemTop ( sp );

        /* verify the edges */
  for ( i = 0; i < nhe; i++ ) {
    if ( mhe[i].v0 < 0 || mhe[i].v0 >= nv ||
         mhe[i].v1 < 0 || mhe[i].v1 >= nv ||
         mhe[i].v0 == mhe[i].v1 )
      FAILURE ( BSM_ERROR_INVALID_HE_V0V1, i )
    if ( mhe[i].facetnum < 0 || mhe[i].facetnum >= nfac )
      FAILURE ( BSM_ERROR_INVALID_HE_FACETNUM, i )
    j = mhe[i].otherhalf;
    if ( j < -1 || j >= nhe || j == i )
      FAILURE ( BSM_ERROR_INVALID_HE_OTHERHALF, i )
    if ( j >= 0 ) {  /* an internal edge */
      if ( mhe[j].otherhalf != i )
        FAILURE ( BSM_ERROR_INVALID_HE_PAIR, i )
      if ( mhe[j].v0 != mhe[i].v1 || mhe[j].v1 != mhe[i].v0 )
        FAILURE ( BSM_ERROR_INVALID_HE_PAIR_V0V1, i )
    }
  }

        /* verify the verticess */
  for ( i = 0; i < nv; i++ ) {
    d = mv[i].degree;
    j = mv[i].firsthalfedge;
    if ( d < 1 || j < 0 || j+d > nhe )
      FAILURE ( BSM_ERROR_INVALID_VERTEX_DEGREE, i )
    for ( k = j; k < j+d-1; k++ )
      if ( mhe[mvhei[k]].v0 != i || mhe[mvhei[k]].otherhalf < 0 )
        FAILURE ( BSM_ERROR_INVALID_VERTEX_INNERHE, i )
          /* only the last edge leaving a point may be represented by */
          /* only one halfedge - at the mesh boundary */
    if ( mhe[mvhei[j+d-1]].v0 != i || mhe[mvhei[j+d-1]].otherhalf < -1 )
      FAILURE ( BSM_ERROR_INVALID_VERTEX_LASTHE, i )
          /* an inner vertex degree must be at least 3 */
    if ( mhe[mvhei[j+d-1]].otherhalf >= 0 && d < 3 )
      FAILURE ( BSM_ERROR_INVALID_INNER_VERTEX_DEG, i );
  }

        /* verify the facets */
  for ( i = 0; i < nfac; i++ ) {
    d = mfac[i].degree;
    j = mfac[i].firsthalfedge;
    if ( d < 3 || j < 0 || j+d > nhe )
      FAILURE ( BSM_ERROR_INVALID_FACET_DEGREE, i )
    for ( k = j; k < j+d-1; k++ ) {
      if ( mhe[mfhei[k]].v1 != mhe[mfhei[k+1]].v0 )
        FAILURE ( BSM_ERROR_INVALID_FACET_VERT, i )
      if ( mhe[mfhei[k]].facetnum != i )
        FAILURE ( BSM_ERROR_INVALID_FACET_HEDGE, i )
    }
    if ( mhe[mfhei[j+d-1]].v1 != mhe[mfhei[j]].v0 )
      FAILURE ( BSM_ERROR_INVALID_FACET_VERT, i )
    if ( mhe[mfhei[j+d-1]].facetnum != i )
      FAILURE ( BSM_ERROR_INVALID_FACET_HEDGE, i )
  }

        /* verify the orientation */
  wlv = pkv_GetScratchMem ( nhe );
  wlf = pkv_GetScratchMem ( nhe );
  if ( !wlv || !wlf )
    FAILURE ( BSM_ERROR_NO_SCRATCH_MEMORY, -1 )
  for ( i = 0; i < nv; i++ ) {
    d = mv[i].degree;
    j = mv[i].firsthalfedge;
    for ( k = 0; k < d; k++ )
      wlv[mvhei[j+k]] = k;
  }
  for ( i = 0; i < nfac; i++ ) {
    d = mfac[i].degree;
    j = mfac[i].firsthalfedge;
    for ( k = 0; k < d; k++ )
      wlf[mfhei[j+k]] = k;
  }
          /* orientation for facets */
  for ( i = 0; i < nfac; i++ ) {
    d = mfac[i].degree;
    j = mfac[i].firsthalfedge;
    he0 = mfhei[j+d-1];
    for ( k = 0; k < d; k++ ) {
      he1 = mfhei[j+k];
      l = wlv[he1];
      he2 = PREVVERT_HEDGE ( mhe[he1].v0, l );
      if ( he2 != mhe[he0].otherhalf )
        FAILURE ( BSM_ERROR_INVALID_FACET_ORIENTATON, i )
      he0 = he1;
    }
  }
          /* orientation for vertices */
  for ( i = 0; i < nv; i++ ) {
    d = mv[i].degree;
    j = mv[i].firsthalfedge;
    he0 = PREVVERT_HEDGE ( i, 0 );
    if ( he0 >= 0 ) {
      for ( k = 0; k < d; k++ ) {
        he1 = mhe[he0].otherhalf;
        he2 = NEXTFAC_HEDGE ( mhe[he1].facetnum, wlf[he1] );
        if ( wlv[he2] != k )
          FAILURE ( BSM_ERROR_INVALID_VERTEX_ORIENTATION, i )
        he0 = he2;
      }
    }
    else {
      he0 = mvhei[j];
      for ( k = 0; k < d-1; k++ ) {
        he1 = mhe[he0].otherhalf;
        he2 = NEXTFAC_HEDGE ( mhe[he1].facetnum, wlf[he1] );
        he0 = mvhei[j+k+1];
        if ( he2 != he0 )
          FAILURE ( BSM_ERROR_INVALID_VERTEX_ORIENTATION, i )
      }
    }
  }
  pkv_SetScratchMemTop ( sp );

        /* vertices of each facet must not appear more than once */
  wlv = pkv_GetScratchMem ( nv );
  if ( !wlv )
    FAILURE ( BSM_ERROR_NO_SCRATCH_MEMORY, -1 );
  memset ( wlv, 0, nv );
  for ( i = 0; i < nfac; i++ ) {
    d = mfac[i].degree;
    j = mfac[i].firsthalfedge;
    for ( k = 0; k < d; k++ ) {
      l = mhe[mfhei[j+k]].v0;
      if ( wlv[l] )
        FAILURE ( BSM_ERROR_DOUBLE_FACET_VERTEX, i );
      wlv[l] = 1;
    }
    for ( k = 0; k < d; k++ )
      wlv[mhe[mfhei[j+k]].v0] = 0;
  }
  pkv_SetScratchMemTop ( sp );

        /* facets surrounding each facet must not appear more than once */
  wlf = pkv_GetScratchMem ( nfac );
  if ( !wlf )
    FAILURE ( BSM_ERROR_NO_SCRATCH_MEMORY, -1 );
  memset ( wlf, 0, nfac );
  for ( i = 0; i < nv; i++ ) {
    d = mv[i].degree;
    j = mv[i].firsthalfedge;
    for ( k = 0; k < d; k++ ) {
      l = mhe[mvhei[j+k]].facetnum;
      if ( wlf[l] )
        FAILURE ( BSM_ERROR_DOUBLE_VERTEX_FACET, i );
      wlf[l] = 1;
    }
    for ( k = 0; k < d; k++ )
      wlf[mhe[mvhei[j+k]].facetnum] = 0;
  }
  pkv_SetScratchMemTop ( sp );

  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  if ( errcode && errelem ) {
    *errcode = error_code;
    *errelem = error_el;
  }
  else
    printf ( "error: %d, %d\n", error_code, error_el );
  return false;
#undef FAILURE
#undef PREVVERT_HEDGE
#undef NEXTFAC_HEDGE
} /*bsm_CheckMeshIntegrity*/

