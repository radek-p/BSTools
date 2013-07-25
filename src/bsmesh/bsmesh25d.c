
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
boolean bsm_GlueHalfedgeLoopsd ( int spdimen,
                                 int inv, BSMvertex *imv, int *imvhei, double *ivc,
                                 int inhe, BSMhalfedge *imhe,
                                 int infac, BSMfacet *imfac, int *imfhei,
                                 int he1, int he2,
                                 int *onv, BSMvertex *omv, int *omvhei, double *ovc,
                                 int *onhe, BSMhalfedge *omhe,
                                 int *onfac, BSMfacet *omfac, int *omfhei )
{
  void *sp;
  int  llgt, i, j, k, e, v, vfhe, vd, w, wfhe, wd, _onv;
  int  *lhe1, *lhe2, *vertn;
  char *vtag;

  sp = pkv_GetScratchMemTop ();
  vtag = pkv_GetScratchMem ( inv );
  if ( !vtag )
    goto failure;
        /* both loops must have the same length - find and check it */
  llgt = bsm_HalfedgeLoopLength ( inv, imv, imvhei, inhe, imhe, he1 );
  if ( llgt < 3 )
    goto failure;
  if ( bsm_HalfedgeLoopLength ( inv, imv, imvhei, inhe, imhe, he2 ) != llgt )
    goto failure;
        /* allocate the arrays for the loop halfedge indices */
  lhe1 = pkv_GetScratchMemi ( 2*llgt );
  if ( !lhe1 )
    goto failure;
  lhe2 = &lhe1[llgt];
        /* find the loop halfedges */
          /* first loop */
  i = 0;
  e = he1;
  do {
    lhe1[i++] = e;
    v = imhe[e].v1;
    vd = imv[v].degree;
    vfhe = imv[v].firsthalfedge;
    e = imvhei[vfhe+vd-1];
  } while ( i < llgt );
          /* second loop */
  lhe2[0] = he2;
  v = imhe[he2].v1;
  vd = imv[v].degree;
  vfhe = imv[v].firsthalfedge;
  e = imvhei[vfhe+vd-1];
  i = llgt-1;
  do {
    lhe2[i--] = e;
    v = imhe[e].v1;
    vd = imv[v].degree;
    vfhe = imv[v].firsthalfedge;
    e = imvhei[vfhe+vd-1];
  } while ( i > 0 );
        /* find out, whether the loops are disjoint */
  for ( i = 0; i < llgt; i++ ) {
    e = lhe1[i];
    v = imhe[e].v0;
    vd = imv[v].degree;
    vfhe = imv[v].firsthalfedge;
    vtag[v] = 1;
    for ( j = 0; j < vd; j++ ) {
      e = imvhei[vfhe+j];
      v = imhe[e].v1;
      vtag[v] = 1;
    }
  }
  for ( i = 0; i < llgt; i++ ) {
    e = lhe2[i];
    v = imhe[e].v0;
    vtag[v] = 0;
  }
  for ( i = 0; i < llgt; i++ ) {
    e = lhe1[i];
    v = imhe[e].v0;
    vd = imv[v].degree;
    vfhe = imv[v].firsthalfedge;
    if ( !vtag[v] )
      goto failure;
    for ( j = 0; j < vd; j++ ) {
      e = imvhei[vfhe+j];
      v = imhe[e].v1;
      if ( !vtag[v] )
        goto failure;
    }
  }
        /* o.k., find the numbering of the vertices to remain */
  vertn = pkv_GetScratchMemi ( inv );
  if ( !vertn )
    goto failure;
  memset ( vertn, 0, inv*sizeof(int) );
  for ( i = 0; i < llgt; i++ ) {
    e = lhe2[i];
    v = imhe[e].v0;
    vertn[v] = -1;
  }
  for ( i = j = 0; i < inv; i++ )
    if ( vertn[i] >= 0 ) {
      vertn[i] = j;
      omv[j].degree = imv[i].degree;
      memcpy ( &ovc[j*spdimen], &ivc[i*spdimen], spdimen*sizeof(double) );
      j ++;
    }
  *onv = _onv = j;
  memcpy ( omhe, imhe, inhe*sizeof(BSMhalfedge) );
  *onhe = inhe;
  memcpy ( omfac, imfac, infac*sizeof(BSMfacet) );
  memcpy ( omfhei, imfhei, inhe*sizeof(int) );
  *onfac = infac;
        /* glue the halfedges */
  for ( i = 0; i < llgt; i++ ) {
    j = lhe1[i];
    k = lhe2[i];
    omhe[j].otherhalf = k;
    omhe[k].otherhalf = j;
    omhe[k].v0 = omhe[j].v1;
    omv[vertn[omhe[j].v0]].degree += imv[omhe[k].v1].degree;
  }
        /* glue the vertices */
  for ( i = k = 0; i < inv; i++ )
    if ( vertn[i] >= 0 ) {
      j = vertn[i];
      omv[j].firsthalfedge = k;
      memcpy ( &omvhei[k], &imvhei[imv[i].firsthalfedge],
               imv[i].degree*sizeof(int) );
      k += omv[j].degree;
    }
  for ( i = 0; i < llgt; i++ ) {
    v = vertn[imhe[lhe1[i]].v0];
    vd = omv[v].degree;
    vfhe = omv[v].firsthalfedge;
    w = imhe[lhe2[i]].v1;
    wd = imv[w].degree;
    wfhe = imv[w].firsthalfedge;
    memcpy ( &omvhei[vfhe+vd-wd], &imvhei[wfhe], wd*sizeof(int) );
    for ( j = 0; j < vd; j++ ) {
      e = omvhei[vfhe+j];
      k = omhe[e].otherhalf;
      omhe[e].v0 = omhe[k].v1 = v;
    }
    for ( j = 0; j < spdimen; j++ )
      ovc[v*spdimen+j] = 0.5*(ovc[v*spdimen+j]+ivc[w*spdimen+j] );
  }
  for ( i = 0; i < _onv; i++ ) {
    vd = omv[i].degree;
    vfhe = omv[i].firsthalfedge;
    for ( j = 0; j < vd; j++ ) {
      e = omvhei[vfhe+j];
      omhe[e].v0 = i;
      v = vertn[imhe[e].v1];
      if ( v >= 0 )
        omhe[e].v1 = v;
      e = omhe[e].otherhalf;
      if ( e >= 0 )
        omhe[e].v1 = i;

    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*bsm_GlueHalfedgeLoopsd*/

