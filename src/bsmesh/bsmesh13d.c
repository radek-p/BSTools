
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010                                  */
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
boolean bsm_FacetEdgeDoublingd ( int spdimen,
                                 int inv, BSMvertex *imv, int *imvhei, double *iptc,
                                 int inhe, BSMhalfedge *imhe,
                                 int infac, BSMfacet *imfac, int *imfhei,
                                 int fn,
                                 int *onv, BSMvertex *omv, int *omvhei, double *optc,
                                 int *onhe, BSMhalfedge *omhe,
                                 int *onfac, BSMfacet *omfac, int *omfhei )
{
  int d, fhe, i, j, k, l, m, n, vd, vf,
      e0, e1, e2, e3, e4, e5, e6, e7, v0, v1, v2, v3;

  if ( fn < 0 || fn >= infac ) {
    *onv = *onhe = *onfac = 0;
    return false;
  }
  d = imfac[fn].degree;
  *onv = inv + d;
  *onhe = inhe + 4*d;
  *onfac = infac + d;

      /* copy all that you can to the output arrays */
  memcpy ( omfac, imfac, infac*sizeof(BSMfacet) );
  memcpy ( omfhei, imfhei, inhe*sizeof(int) );
  memcpy ( omhe, imhe, inhe*sizeof(BSMhalfedge) );
  memcpy ( omv, imv, inv*sizeof(BSMvertex) );
  memcpy ( optc, iptc, inv*spdimen*sizeof(double) );
  fhe = imfac[fn].firsthalfedge;
  for ( i = 0; i < d; i++ ) {
    v0 = imhe[imfhei[fhe+i]].v0;
    omv[v0].degree ++;
    if ( omv[v0].degree <= 0 )
      return false;
  }
  for ( i = k = 0;  i < inv;  i++ ) {
    memcpy ( &omvhei[k], &imvhei[imv[i].firsthalfedge],
             imv[i].degree*sizeof(int) );
    omv[i].firsthalfedge = k;
    k += omv[i].degree;
  }
      /* generate the new data */
  for ( i = 0, k = inhe+d, l = inhe,
        v3 = imhe[imfhei[fhe]].v0, v2 = *onv-1,
        e6 = *onhe-2, e7 = *onhe-1;
        i < d;
        i++, k += 3, l += 4,
        v3 = v0, v2 = v1,
        e6 = e4, e7 = e5 ) {
    e0 = imfhei[fhe+i];
    e1 = imhe[e0].otherhalf;
    e2 = inhe+4*i;  e3 = e2+1;  e4 = e3+1;  e5 = e4+1;
    v0 = imhe[e0].v1;
    v1 = inv+i;
    omhe[e2].otherhalf = e0;
    omhe[e0].otherhalf = e2;
    omhe[e4].otherhalf = e5;
    omhe[e5].otherhalf = e4;
    omhe[e3].otherhalf = e1;
    if ( e1 >= 0) omhe[e1].otherhalf = e3;
    omhe[e0].v0 = omhe[e2].v1 = v2;
    omhe[e0].v1 = omhe[e4].v1 = omhe[e2].v0 = omhe[e5].v0 = v1;
    omhe[e3].v0 = v3;
    omhe[e4].v0 = omhe[e5].v1 = omhe[e3].v1 = v0;
    omhe[e2].facetnum = omhe[e3].facetnum =
    omhe[e4].facetnum = omhe[e7].facetnum = infac+i;
    omv[v1].degree = 3;
    omv[v1].firsthalfedge = k;
    omvhei[k] = e5;
    omvhei[k+1] = e2;
    omvhei[k+2] = imfhei[fhe + (i+1)%d];
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
    for ( j = n = 0;  j < vd;  j++, n++ ) {
      if ( imvhei[m+j] == e0 ) {
        omvhei[vf+n] = e6;
        n ++;
        omvhei[vf+n] = e3;
      }
      else
        omvhei[vf+n] = imvhei[m+j];
    }
  }

  return true;
} /*bsm_FacetEdgeDoublingd*/

