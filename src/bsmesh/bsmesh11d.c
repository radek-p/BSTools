
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
boolean bsm_RemoveFacetd ( int spdimen,
                           int inv, BSMvertex *imv, int *imvhei, double *iptc,
                           int inhe, BSMhalfedge *imhe,
                           int infac, BSMfacet *imfac, int *imfhei,
                           int nfr,
                           int *onv, BSMvertex *omv, int *omvhei, double *optc,
                           int *onhe, BSMhalfedge *omhe,
                           int *onfac, BSMfacet *omfac, int *omfhei )
{
  void *sp;
  BSMvertex   *wmv;
  BSMhalfedge *wmhe;
  BSMfacet    *wmfac;
  int         *wmvhei, *wmfhei;
  int         *newvi, *newhei, *newfi;
  double      *wptc;
  int         i, j, k, fd, ffhe, v0, vd, vfhe, hn, iihei, newvn;
  char        *vtag, *ftag;

  sp = pkv_GetScratchMemTop ();
        /* make a working copy of the mesh */
          /* in case there are vertices to be split, allocate arrays long enough */
  fd = imfac[nfr].degree;
  ffhe = imfac[nfr].firsthalfedge;
  for ( i = k = 0;  i < fd;  i++ )
    k += imv[imhe[imfhei[ffhe+i]].v0].degree - 2;
  k = max ( k, 0 );
  wmv = pkv_GetScratchMem ( (inv+fd)*sizeof(BSMvertex) );
  wmhe = pkv_GetScratchMem ( (inhe+k)*sizeof(BSMhalfedge) );
  wmfac = pkv_GetScratchMem ( infac*sizeof(BSMfacet) );
  wmvhei = pkv_GetScratchMemi ( inhe+k );
  wmfhei = pkv_GetScratchMemi ( inhe+k );
  wptc = pkv_GetScratchMemd ( fd*spdimen );
  newvi = pkv_GetScratchMemi ( inv+fd );
  newhei = pkv_GetScratchMemi ( inhe+k );
  newfi = pkv_GetScratchMemi ( infac );
  vtag = pkv_GetScratchMem ( inv );
  ftag = pkv_GetScratchMem ( infac );
  if ( !wmv || !wmhe || !wmfac || !wmvhei || !wmfhei || !wptc ||
       !newvi || !newhei || !newfi || !inv || !infac )
    goto failure;
  memcpy ( wmv, imv, inv*sizeof(BSMvertex) );
  memcpy ( wmhe, imhe, inhe*sizeof(BSMhalfedge) );
  memcpy ( wmfac, imfac, infac*sizeof(BSMfacet) );
  memcpy ( wmvhei, imvhei, inhe*sizeof(int) );
  memcpy ( wmfhei, imfhei, inhe*sizeof(int) );
        /* mark the vertices, halfedges and the facet to be removed */
  for ( i = 0; i < inv; i++ ) {
    vd = wmv[i].degree;
    vfhe = wmv[i].firsthalfedge;
    if ( wmhe[wmvhei[vfhe+vd-1]].otherhalf >= 0 )
      vtag[i] = 0;  /* inner vertex */
    else
      vtag[i] = 1;  /* boundary vertex */
  }
  iihei = inhe;
  memset ( newhei, 0, inhe*sizeof(int) );
  newvn = 0;
  for ( i = 0; i < infac; i++ )
    ftag[i] = 0;
  ftag[nfr] = 1;
  fd = wmfac[nfr].degree;
  ffhe = wmfac[nfr].firsthalfedge;
  for ( i = 0; i < fd; i++ ) {
    hn = wmfhei[ffhe+i];
    newhei[hn] = -1;  /* halfedge to be removed */
    v0 = wmhe[hn].v0;
    vd = wmv[v0].degree;
    vfhe = wmv[v0].firsthalfedge;
    if ( vtag[v0] > 1 )
      goto failure;
    if ( wmv[v0].degree == 1 )
      vtag[v0] = 2;  /* vertex to be removed */
    else {
      vtag[v0] = 3;  /* vertex to remain, with less edges */
      if ( wmhe[wmvhei[vfhe+vd-1]].otherhalf < 0 ) {
        if ( wmvhei[vfhe] == hn )
          memmove ( &wmvhei[vfhe], &wmvhei[vfhe+1], (vd-1)*sizeof(int) );
        else if ( wmvhei[vfhe+vd-1] != hn ) { /* vertex to be split */
          j = vd-2;
          while ( j > 0 && wmvhei[vfhe+j] != hn )
            j --;
          if ( j <= 0 )
            goto failure;
          wmv[inv].degree = vd = vd-j-1;
          vtag[inv] = 4;             /* a new vertex copy */
          wmv[inv].firsthalfedge = iihei;
          memcpy ( &wmvhei[inhe], &wmvhei[vfhe+j+1], vd*sizeof(int) );
          wmhe[wmhe[hn].otherhalf].v1 = inv;
          for ( k = 0; k < vd-1; k++ ) {
            wmhe[wmvhei[inhe+k]].v0 = inv;
            wmhe[wmhe[inhe+k].otherhalf].v1 = inv;
          }
          wmhe[wmvhei[inhe+vd-1]].v0 = inv;
          memcpy ( &wptc[newvn*spdimen], &iptc[v0*spdimen], spdimen*sizeof(double) );
          newvn ++;
          inv ++;
          iihei += vd;
          wmv[v0].degree = j+1;
        }
      }
      else {
        if ( !_bsm_RotateHalfedgei ( vd, &wmvhei[vfhe], hn ) )
          goto failure;
      }
      wmv[v0].degree --;
    }
  }
        /* renumber the vertices */
  for ( i = k = 0; i < inv; i++ )
    if ( vtag[i] != 2 )
      newvi[i] = k++;
    else
      newvi[i] = -1;
  *onv = k;
        /* renumber the halfedges */
  for ( i = k = 0; i < inhe; i++ )
    if ( newhei[i] >= 0 )
      newhei[i] = k++;
  *onhe = k;
        /* renumber the facets */
  for ( i = k = 0; i < infac; i++ )
    if ( i != nfr )
      newfi[i] = k++;
    else
      newfi[i] = -1;
  *onfac = k;
        /* store the vertices in the final arrays */
  vfhe = newvn = 0;
  for ( i = 0; i < inv; i++ ) {
    k = newvi[i];
    if ( k >= 0 ) {
      omv[k].degree = vd = wmv[i].degree;
      omv[k].firsthalfedge = vfhe;
      for ( j = 0; j < vd; j++ )
        omvhei[vfhe+j] = newhei[wmvhei[imv[i].firsthalfedge+j]];
      if ( vtag[i] == 4 ) {
        memcpy ( &optc[k*spdimen], &wptc[newvn*spdimen], spdimen*sizeof(double) );
        newvn ++;
      }
      else
        memcpy ( &optc[k*spdimen], &iptc[i*spdimen], spdimen*sizeof(double) );
      vfhe += vd;
    }
  }

  _bsm_OutputWorkMeshd ( spdimen, inv, wmv, wmvhei, NULL, inhe, wmhe,
                         infac, wmfac, wmfhei, newvi, newhei, newfi,
                         onv, omv, omvhei, optc, onhe, omhe, onfac, omfac, omfhei );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*bsm_RemoveFacetd*/

