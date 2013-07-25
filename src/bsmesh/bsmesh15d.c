
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
boolean bsm_RemoveVertexd ( int spdimen,
                            int inv, BSMvertex *imv, int *imvhei, double *iptc,
                            int inhe, BSMhalfedge *imhe,
                            int infac, BSMfacet *imfac, int *imfhei,
                            int nvr,
                            int *onv, BSMvertex *omv, int *omvhei, double *optc,
                            int *onhe, BSMhalfedge *omhe,
                            int *onfac, BSMfacet *omfac, int *omfhei )
{
  void *sp;
  int  d, fhe, i, j, k, l, m, n, o, p, q, fd, v0, v1;
  int  *newhenum;
  char *vtag, *ftag;

  sp = pkv_GetScratchMemTop ();
  if ( nvr < 0 || nvr >= inv )  /* invalid vertex number */
    goto failure;

  newhenum = pkv_GetScratchMemi ( inhe );
  vtag = pkv_GetScratchMem ( inv );
  ftag = pkv_GetScratchMem ( infac );
  if ( !newhenum || !vtag || !ftag )
    goto failure;

  d = imv[nvr].degree;
  fhe = imv[nvr].firsthalfedge;
  *onv = inv-1;
  *onfac = infac-d+1;
  memset ( newhenum, 0, inhe*sizeof(int) );
  memset ( vtag, 0, inv );
  vtag[nvr] = 2;
        /* tag the facets to be removed */
  memset ( ftag, 0, infac );
  for ( i = 0; i < d; i++ ) {
    j = imhe[imvhei[fhe+i]].facetnum;
    ftag[j] = 1;
  }

  if ( imhe[imvhei[fhe+d-1]].otherhalf < 0 ) { /* a boundary vertex is to be removed */
    if ( imv[nvr].degree == 1 && imfac[imhe[fhe].facetnum].degree == 3 )
      goto failure;
    *onhe = inhe-2*d+1;
          /* tag the vertices, whose degree will change */
    for ( i = 0; i < d-1; i++ )
      vtag[imhe[imvhei[fhe+i]].v1] = 1;
          /* set numbers of remaining halfedges */
    for ( i = k = 0; i < inhe; i++ )
      if ( imhe[i].v0 == nvr || imhe[i].v1 == nvr )
        newhenum[i] = -1;
      else
        newhenum[i] = k ++;
          /* copy the remaining halfedges */
    for ( i = 0; i < inhe; i++ )
      if ( newhenum[i] >= 0 ) {
        k = newhenum[i];
        omhe[k].v0 = imhe[i].v0 < nvr ? imhe[i].v0 : imhe[i].v0-1;
        omhe[k].v1 = imhe[i].v1 < nvr ? imhe[i].v1 : imhe[i].v1-1;
        if ( imhe[i].otherhalf >= 0 )
          omhe[k].otherhalf = newhenum[imhe[i].otherhalf];
        else
          omhe[k].otherhalf = -1;
      }
          /* generate the new halfedge */
    k = *onhe-1;
    i = imhe[imvhei[fhe]].facetnum;
    j = imfac[i].firsthalfedge;
    fd = imfac[i].degree;
    for ( l = 0; l < fd; l++ )
      if ( imhe[imfhei[j+l]].v1 == nvr )
        break;
    if ( l >= fd )
      goto failure;
    v0 = imhe[imfhei[j+l]].v0;
    v1 = imhe[imvhei[fhe+d-1]].v1;
    omhe[k].v0 = v0 = v0 < nvr ? v0 : v0-1;
    omhe[k].v1 = v1 < nvr ? v1 : v1-1;
    omhe[k].facetnum = *onfac-1;
    omhe[k].otherhalf = -1;
          /* copy the remaining facets */
    for ( i = j = k = 0; i < infac; i++ )
      if ( !ftag[i] ) {
        omfac[j].degree = fd = imfac[i].degree;
        omfac[j].firsthalfedge = k;
        for ( l = 0; l < fd; l++ ) {
          omfhei[k+l] = newhenum[imfhei[imfac[i].firsthalfedge+l]];
          omhe[omfhei[k+l]].facetnum = j;
        }
        k += fd;
        j++;
      }
          /* setup the new facet */
    for ( i = fd = 0; i < d; i++ )
      fd += imfac[imhe[imvhei[fhe+i]].facetnum].degree;
    fd -= 2*d-1;
    if ( fd > 127 )  /* the facet would have too many edges */
      goto failure;
        /* at this point j ought to be equal to *obfac-1 */
        /* and k ought to be the index of the first halfedge number */
        /* of the new facet */
    omfac[j].degree = fd;
    omfac[j].firsthalfedge = k;
    for ( i = d-1, l = 0; i >= 0; i-- ) {
      m = imhe[imvhei[fhe+i]].facetnum;
      n = imfac[m].degree;
      o = imfac[m].firsthalfedge;
      for ( p = 0; p < n; p++ )
        if ( imhe[imfhei[o+p]].v0 == nvr )
          break;
        if ( p >= n )
          goto failure;
        for ( q = 0; q < n-2; q++ ) {
          p = (p+1) % n;
          omfhei[k] = newhenum[imfhei[o+p]];
          omhe[omfhei[k]].facetnum = j;
          k ++;
        }
    }
    omfhei[k] = k;
          /* copy the vertices */
    for ( i = j = k = 0; i < inv; i++ )
      if ( i != nvr ) {
        memcpy ( &optc[j*spdimen], &iptc[i*spdimen], spdimen*sizeof(double) );
        omv[j].firsthalfedge = k;
        omv[j].degree = d = vtag[i] ? imv[i].degree-1 : imv[i].degree;
        if ( imhe[imvhei[imv[i].firsthalfedge+imv[i].degree-1]].otherhalf == -1 ) {
          if ( d < 1 )
            goto failure;
        }
        else {
          if ( d <= 2 )
            goto failure;
        }
        for ( n = o = 0; n < imv[i].degree; n++ )
          if ( newhenum[imvhei[imv[i].firsthalfedge+n]] >= 0 ) {
            omvhei[k+o] = newhenum[imvhei[imv[i].firsthalfedge+n]];
            o ++;
          }
        j ++;
        k += d;
      }
    omvhei[omv[v0].firsthalfedge+omv[v0].degree-1] = *onhe-1;
  }
  else { /* an inner vertex is to be removed */
    *onhe = inhe-2*d;
          /* tag the vertices, whose degree will change */
    for ( i = 0; i < d; i++ ) {
      newhenum[imvhei[fhe+i]] = -1;
      newhenum[imhe[imvhei[fhe+i]].otherhalf] = -1;
      vtag[imhe[imvhei[fhe+i]].v1] = 1;
    }
          /* set numbers of remaining halfedges */
    for ( i = k = 0; i < inhe; i++ )
      if ( newhenum[i] >= 0 )
        newhenum[i] = k ++;
          /* copy the remaining halfedges */
    for ( i = 0; i < inhe; i++ )
      if ( newhenum[i] >= 0 ) {
        k = newhenum[i];
        omhe[k].v0 = imhe[i].v0 < nvr ? imhe[i].v0 : imhe[i].v0-1;
        omhe[k].v1 = imhe[i].v1 < nvr ? imhe[i].v1 : imhe[i].v1-1;
        if ( imhe[i].otherhalf >= 0 )
          omhe[k].otherhalf = newhenum[imhe[i].otherhalf];
        else
          omhe[k].otherhalf = -1;
      }
          /* copy the remaining facets */
    for ( i = j = k = 0; i < infac; i++ )
      if ( !ftag[i] ) {
        omfac[j].degree = fd = imfac[i].degree;
        omfac[j].firsthalfedge = k;
        for ( l = 0; l < fd; l++ ) {
          omfhei[k+l] = newhenum[imfhei[imfac[i].firsthalfedge+l]];
          omhe[omfhei[k+l]].facetnum = j;
        }
        k += fd;
        j++;
      }
        /* setup the new facet */
    for ( i = fd = 0; i < d; i++ )
      fd += imfac[imhe[imvhei[fhe+i]].facetnum].degree;
    fd -= 2*d;
    if ( fd > 127 )  /* the facet would have too many edges */
      goto failure;
        /* at this point j ought to be equal to *onfac-1 */
        /* and k ought to be the index of the first halfedge number */
        /* of the new facet */
    omfac[j].degree = fd;
    omfac[j].firsthalfedge = k;
    for ( i = d-1, l = 0; i >= 0; i-- ) {
      m = imhe[imvhei[fhe+i]].facetnum;
      n = imfac[m].degree;
      o = imfac[m].firsthalfedge;
      for ( p = 0; p < n; p++ )
        if ( imhe[imfhei[o+p]].v0 == nvr )
          break;
      if ( p >= n )
        goto failure;
      for ( q = 0; q < n-2; q++ ) {
        p = (p+1) % n;
        omfhei[k] = newhenum[imfhei[o+p]];
        omhe[omfhei[k]].facetnum = j;
        k ++;
      }
    }
        /* copy the vertices */
    for ( i = j = k = 0; i < inv; i++ )
      if ( i != nvr ) {
        memcpy ( &optc[j*spdimen], &iptc[i*spdimen], spdimen*sizeof(double) );
        omv[j].firsthalfedge = k;
        omv[j].degree = d = vtag[i] ? imv[i].degree-1 : imv[i].degree;
        if ( imhe[imvhei[imv[i].firsthalfedge+imv[i].degree-1]].otherhalf == -1 ) {
          if ( d < 1 )
            goto failure;
        }
        else {
          if ( d <= 2 )
            goto failure;
        }
        for ( n = o = 0; n < imv[i].degree; n++ )
          if ( newhenum[imvhei[imv[i].firsthalfedge+n]] >= 0 ) {
            omvhei[k+o] = newhenum[imvhei[imv[i].firsthalfedge+n]];
            o ++;
          }
        j++;
        k += d;
      }
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  *onv = *onhe = *onfac = -1;
  pkv_SetScratchMemTop ( sp );
  return false;
} /*bsm_RemoveVertexd*/

