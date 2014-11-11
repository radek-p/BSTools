
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */
/* This file was written by Anna Sierhej                                     */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pkvaria.h"
#include "pknum.h"
#include "bsmesh.h"
#include "bsmprivate.h"

/* ////////////////////////////////////////////////////////////////////////// */
void bsm_TriangulateFacetsNum ( int inv, BSMvertex *imv, int *imvhei,
                                int inhe, BSMhalfedge *imhe,
                                int infac, BSMfacet *imfac, int *imfhei,
                                int *onv, int *onhe, int *onfac )
{
  int i, d;

  *onv = inv;
  *onhe = inhe;
  *onfac = infac;
  for ( i = 0; i < infac; i++ ) {
    if ( imfac[i].degree > 3 ) {
      d = imfac[i].degree;
      *onhe = *onhe+(d-3)*2;
      *onfac = *onfac+(d-3);
    }
  }
} /*bsm_TriangulateFacetsNum*/

boolean bsm_TriangulateFacetsd ( int spdimen, int inv,
                                 const BSMvertex *imv, const int *imvhei,
                                 double *iptc, int inhe, const BSMhalfedge *imhe,
                                 int infac, const BSMfacet *imfac, const int *imfhei,
                                 int *onv, BSMvertex *omv, int *omvhei,
                                 double *optc, int *onhe, BSMhalfedge *omhe,
                                 int *onfac, BSMfacet *omfac, int *omfhei )
{
  void        *sp;
  BSMvertex   *tmv,  *t2mv;
  BSMhalfedge *tmhe, *t2mhe;
  BSMfacet    *tmfac, *t2mfac;
  int         tnv, tnhe, tnfac, t2nv, t2nhe, t2nfac;
  int         *tmvhei, *tmfhei, *t2mvhei, *t2mfhei;
  int         i, j, k, d, fhe, v0, v1;
  double      *tptc, *t2ptc;
  boolean     isT, isK, isJ;

  sp = pkv_GetScratchMemTop ();

      /*zarezerwowanie miejsca na tablice na t*/
  tmfac = pkv_GetScratchMem ( *onfac*sizeof(BSMfacet) );
  tmfhei = pkv_GetScratchMemi ( *onhe );
  tmhe = pkv_GetScratchMem ( *onhe*sizeof(BSMhalfedge) );
  tmv = pkv_GetScratchMem ( *onv*sizeof(BSMvertex) );
  tmvhei = pkv_GetScratchMemi ( *onhe );
  tptc = pkv_GetScratchMemd ( *onv*spdimen );
  t2mfac = pkv_GetScratchMem ( *onfac*sizeof(BSMfacet) );
  t2mfhei = pkv_GetScratchMemi ( *onhe );
  t2mhe = pkv_GetScratchMem ( *onhe*sizeof(BSMhalfedge) );
  t2mv = pkv_GetScratchMem ( *onv*sizeof(BSMvertex) );
  t2mvhei = pkv_GetScratchMemi ( *onhe );
  t2ptc = pkv_GetScratchMemd ( *onv*spdimen );
  if ( !tmfac || !tmfhei || !tmhe || !tmv || !tmvhei || !tptc ||
       !t2mfac || !t2mfhei || !t2mhe || !t2mv || !t2mvhei || !t2ptc )
    goto failure;

  memcpy ( tmfac, imfac, infac*sizeof(BSMfacet) );
  memcpy ( tmfhei, imfhei, inhe*sizeof(int) );
  memcpy ( tmhe, imhe, inhe*sizeof(BSMhalfedge) );
  memcpy ( tmv, imv, inv*sizeof(BSMvertex) );
  memcpy ( tmvhei, imvhei, inhe*sizeof(int) );
  memcpy ( tptc, iptc, inv*spdimen*sizeof(double) );

  tnv = inv;
  tnhe = inhe;
  tnfac = infac;
  isT = true;

  for ( i = 0; i < infac; i++ ) {
    if ( imfac[i].degree > 3 ) {
      isK = true;
      isJ = false;
      d = imfac[i].degree;
      fhe = imfac[i].firsthalfedge;
      for ( j = 0, k = 1;  j+k < d-2; ) {
        if ( isT ) {
          v0 = imhe[imfhei[j+fhe]].v0;
          v1 = imhe[imfhei[fhe+d-k-1]].v0;
          if ( !bsm_DivideFacetd ( spdimen, tnv, tmv, tmvhei, tptc, tnhe, tmhe,
                             tnfac, tmfac, tmfhei, v0, v1,
                             &t2nv, t2mv, t2mvhei, t2ptc, &t2nhe, t2mhe,
                             &t2nfac, t2mfac, t2mfhei ) )
            goto failure;
        }
        else{
          v0 = imhe[imfhei[j+fhe]].v0;
          v1 = imhe[imfhei[fhe+d-k-1]].v0;
          if ( !bsm_DivideFacetd ( spdimen, t2nv, t2mv, t2mvhei, t2ptc, t2nhe, t2mhe,
                                   t2nfac, t2mfac, t2mfhei, v0, v1,
                                   &tnv, tmv, tmvhei, tptc, &tnhe, tmhe,
                                   &tnfac, tmfac, tmfhei ) )
            goto failure;
        }
        isT = !isT;
        isK = !isK;
        isJ = !isJ;
        if ( isJ )
          j ++;
        if ( isK )
          k ++;
      }
    }
  }

  if ( isT ) {
        /*najnowsza siatka znajduje sie w t*/
    *onv = tnv;
    *onhe = tnhe;
    *onfac = tnfac;
    memcpy ( omfac, tmfac, *onfac*sizeof(BSMfacet) );
    memcpy ( omfhei, tmfhei, *onhe*sizeof(int) );
    memcpy ( omhe, tmhe, *onhe*sizeof(BSMhalfedge) );
    memcpy ( omv, tmv, *onv*sizeof(BSMvertex) );
    memcpy ( omvhei, tmvhei, *onhe*sizeof(int) );
    memcpy ( optc, tptc, *onv*spdimen*sizeof(double) );
  }
  else{
        /*najnowsza siatka znajduje sie w t2*/
    *onv = t2nv;
    *onhe = t2nhe;
    *onfac = t2nfac;
    memcpy ( omfac, t2mfac, *onfac*sizeof(BSMfacet) );
    memcpy ( omfhei, t2mfhei, *onhe*sizeof(int) );
    memcpy ( omhe, t2mhe, *onhe*sizeof(BSMhalfedge) );
    memcpy ( omv, t2mv, *onv*sizeof(BSMvertex) );
    memcpy ( omvhei, t2mvhei, *onhe*sizeof(int) );
    memcpy ( optc, t2ptc, *onv*spdimen*sizeof(double) );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*bsm_TriangulateFacetsd*/

