
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
void bsm_EdgeLoopDoublingNum ( int inv, int inhe, int infac,
                               int EdgeLoopLength,
                               int *onv, int *onhe, int *onfac )
{
  *onv = inv + EdgeLoopLength;
  *onhe = inhe + 4*EdgeLoopLength;
  *onfac = infac + EdgeLoopLength;
} /*bsm_EdgeLoopDoublingNum*/

boolean bsm_EdgeLoopDoublingd ( int spdimen,
                                int inv, BSMvertex *imv, int *imvhei, double *iptc,
                                int inhe, BSMhalfedge *imhe,
                                int infac, BSMfacet *imfac, int *imfhei,
                                int EdgeLoopLength, int *EdgeLoop,
                                int *onv, BSMvertex *omv, int *omvhei, double *optc,
                                int *onhe, BSMhalfedge *omhe,
                                int *onfac, BSMfacet *omfac, int *omfhei )
{
  void    *sp;
  int     i, j, k, l, m, n, vd, vf, vf2, fhe,
          e0, e1, e2, e3, e4, e5, e6, e7, e8, e9, v0, v1, v2, v3;
  int     pfn, pffhe, phe;
  boolean ce8, ce9, incycle;
  boolean *fiel;
  pkv_queue *q;

  sp = pkv_GetScratchMemTop ();
  q = NULL;
  fiel = pkv_GetScratchMem ( infac*sizeof(boolean) );
  if ( !fiel )
    goto failure;
  memset ( fiel, false, infac*sizeof(boolean) );

  q = pkv_InitQueue ( infac, sizeof(int) );
  if ( !q )
    goto failure;
  vf = imhe[EdgeLoop[0]].facetnum;
  fiel[vf] = true;
  pkv_QueueInsert ( q, &vf );
  do {
    pkv_QueueRemoveFirst(q, &vf2);
    for ( i = imfac[vf2].firsthalfedge;
          i < imfac[vf2].degree+imfac[vf2].firsthalfedge;
          i++ ) {
      fhe = imfhei[i];
      incycle = false;
      for ( j = 0; j < EdgeLoopLength; j++ ) {
        if ( EdgeLoop[j] == fhe )
          incycle = true;
      }
      if ( !incycle && imhe[fhe].otherhalf != -1 ) {
        vf = imhe[imhe[fhe].otherhalf].facetnum;
        if ( !fiel[vf] ) {
          fiel[vf] = true;
          pkv_QueueInsert ( q, &vf );
        }
      }
    }
  } while ( !pkv_QueueEmpty ( q ) );

  *onv = inv + EdgeLoopLength;
  *onhe = inhe + 4*EdgeLoopLength;
  *onfac = infac + EdgeLoopLength;

        /* copy all that you can to the output arrays */
  memcpy ( omfac, imfac, infac*sizeof(BSMfacet) );
  memcpy ( omfhei, imfhei, inhe*sizeof(int) );
  memcpy ( omhe, imhe, inhe*sizeof(BSMhalfedge) );
  memcpy ( omv, imv, inv*sizeof(BSMvertex) );
  memcpy ( optc, iptc, inv*spdimen*sizeof(double) );

  for ( i = 0; i < EdgeLoopLength; i++ ) {
    v1 = imhe[EdgeLoop[i]].v1;
        /* the vertex degree is increased only if there is no halfedge not */
        /* being the loop element, between the loop halfedges. Otherwise the */
        /* vertex degree may even decrease */
    vf = imhe[EdgeLoop[i]].facetnum;
    vf2 = imhe[EdgeLoop[(i+1)%EdgeLoopLength]].facetnum;
    if ( vf == vf2 )
      omv[v1].degree ++;
    else {
      omv[v1].degree++;
      while ( vf != vf2 ) {
        pffhe = imfac[vf].firsthalfedge;
        while ( imhe[imfhei[pffhe]].v0 != v1 )
          pffhe ++;
        vf = imhe[imhe[imfhei[pffhe]].otherhalf].facetnum;
        omv[v1].degree --;
      }
    }
    if ( omv[v1].degree <= 0 )
      goto failure;
  }

  for ( i = m = 0;  i < inv;  i++ ) {
    memcpy ( &omvhei[m], &imvhei[imv[i].firsthalfedge],
             imv[i].degree*sizeof(int) );
    omv[i].firsthalfedge = m;
    m += omv[i].degree;
  }
  k = m - EdgeLoopLength;

        /* generate the new data */
  for ( i = 0, k += EdgeLoopLength, l = inhe,
        v3 = imhe[EdgeLoop[0]].v0, v2 = *onv-1,
        e6 = *onhe-2, e7 = *onhe-1;
        i < EdgeLoopLength;
        i++, l += 4,
        v3 = v0, v2 = v1,
        e6 = e4, e7 = e5 ) {
    e0 = EdgeLoop[i];
    e1 = imhe[e0].otherhalf;
    e2 = inhe + 4*i;  e3 = e2 + 1;  e4 = e3 + 1;  e5 = e4 + 1;
    v0 = imhe[e0].v1;
    v1 = inv+i;
    omhe[e2].otherhalf = e0;
    omhe[e0].otherhalf = e2;
    omhe[e4].otherhalf = e5;
    omhe[e5].otherhalf = e4;
    omhe[e3].otherhalf = e1;
    if ( e1 >= 0)
      omhe[e1].otherhalf = e3;
    omhe[e0].v0 = omhe[e2].v1 = v2;
    omhe[e0].v1 = omhe[e4].v1 = omhe[e2].v0 = omhe[e5].v0 = v1;
    omhe[e3].v0 = v3;
    omhe[e4].v0 = omhe[e5].v1 = omhe[e3].v1 = v0;
    omhe[e2].facetnum = omhe[e3].facetnum =
    omhe[e4].facetnum = omhe[e7].facetnum = infac+i;
      /* here degree = 3 + number of unmarked halfedges between e0 */
      /* and the next element of the loop */
    omv[v1].degree = 3;
    omv[v1].firsthalfedge = k;
    omvhei[k] = e5;
    k ++;
    omvhei[k] = e2;
    k ++;
    pfn = imhe[e0].facetnum;
    while ( pfn != imhe[EdgeLoop[(i+1) % EdgeLoopLength]].facetnum ) {
      pffhe = imfac[pfn].firsthalfedge;
      phe = imfac[pfn].firsthalfedge;
      while ( imhe[imfhei[phe]].v0 != v0 )
        phe ++;
      omvhei[k] = imfhei[phe];
      omv[v1].degree ++;
      k ++;
      pfn = imhe[imhe[imfhei[phe]].otherhalf].facetnum;
    }
    omvhei[k] = EdgeLoop[(i+1)%EdgeLoopLength];
    k ++;

    memcpy ( &optc[v1*spdimen], &iptc[v0*spdimen], spdimen*sizeof(double) );
    omfac[infac+i].degree = 4;
    omfac[infac+i].firsthalfedge = l;
    omfhei[l] = e2;
    omfhei[l+1] = e7;
    omfhei[l+2] = e3;
    omfhei[l+3] = e4;
    vd = imv[v3].degree;
    vf = omv[v3].firsthalfedge;
    m = imv[v3].firsthalfedge;
    for ( j = n = 0;  j < vd;  j++) {
      if ( imvhei[m+j] == e0 ) {
        omvhei[vf+n] = e6;
        n ++;
        omvhei[vf+n] = e3;
        n ++;
      }
      else {
        if ( fiel[imhe[imvhei[m+j]].facetnum] ) {
            /* is there another edge from the loop ? */
        }
        else {
          omvhei[vf+n] = imvhei[m+j];
          n ++;
        }
      }
    }
  }

    /* exchanging e8 and e9 in appropriate facets */
  for ( i = 0; i < EdgeLoopLength; i++ ) {
      /* searching for the location of our halfedge */
    vf = imhe[EdgeLoop[i]].facetnum;
    vd = imfac[vf].degree;
    fhe = imfac[vf].firsthalfedge;
    for ( j = fhe; j < fhe+vd; j++ ) {
      if ( imfhei[j] == EdgeLoop[i] ) {
        k = (j-fhe-1) % vd;
        if ( k < 0 )
          k += vd;
        e8 = imfhei[fhe+k];
        e9 = imfhei[fhe + (j-fhe+1) % vd];
          /* if e8 is not an element of the loop, it must be renumbered */
          /* similarly e9 */
        ce8 = false;
        ce9 = false;
        for ( l = 0; l < EdgeLoopLength; l++ ) {
          if ( EdgeLoop[l] == e8 )
            ce8 = true;
          if ( EdgeLoop[l] == e9 )
            ce9 = true;
        }
        if ( !ce8 ) {
          omhe[e8].v1 = omhe[EdgeLoop[i]].v0;
          omhe[omhe[e8].otherhalf].v0 = omhe[e8].v1;
        }
        if ( !ce9 ) {
          omhe[e9].v0 = omhe[EdgeLoop[i]].v1;
          omhe[omhe[e9].otherhalf].v1 = omhe[e9].v0;
        }
      }
    }
  }
  if ( q ) PKV_FREE ( q );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( q ) PKV_FREE ( q );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*bsm_EdgeLoopDoublingd*/

