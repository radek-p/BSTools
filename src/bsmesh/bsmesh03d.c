
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

/* This procedure uses the same doubling algorithm that the bsm_DoublingMatd */  
/* procedure, which is intended to obtain the same correspondence between    */
/* the vertices of the input and output mesh - the only difference is        */
/* that itcopies the input vertices instead of assembling the doubling       */
/* Should any modifications be necessary, the same have to be made in the    */
/* code of bsm_DoublingMatd.                                                 */
boolean bsm_Doublingd ( int spdimen,
                        int inv, BSMvertex *imv, int *imvhei, double *ivc,
                        int inhe, BSMhalfedge *imhe,
                        int infac, BSMfacet *imfac, int *imfhei,
                        int *onv, BSMvertex *omv, int *omvhei, double *ovc,
                        int *onhe, BSMhalfedge *omhe,
                        int *onfac, BSMfacet *omfac, int *omfhei )
{
#define PREVIFAC_HEDGE(fn,en) \
  ((int)(en) > 0 ? \
    imfhei[imfac[fn].firsthalfedge+(int)(en)-1] : \
    imfhei[imfac[fn].firsthalfedge+imfac[fn].degree-1])

  void *sp;
  int  vi, vb, ei, eb, _onfac, _onhe, _onv, fvf, fvhe;
  int  *vcn, *ecn, *efn;
  char *wlv, *wlf, *vtag, *ftag;
  int  i, j, k, l, m, f, d, v0, v1, p, ecni, ecnl;

  sp = pkv_GetScratchMemTop ();
  vtag = pkv_GetScratchMem ( inv );
  ftag = pkv_GetScratchMem ( infac );
  if ( !vtag || !ftag )
    goto failure;

  bsm_TagMesh ( inv, imv, imvhei, inhe, imhe, infac, imfac, imfhei,
                vtag, ftag, &vi, &vb, &ei, &eb );
        /* the number of output facets */
  *onfac = _onfac = infac + ei + eb + inv;
        /* the number of output halfedges */
  *onhe = _onhe = 8*ei + 6*eb + 2*vb;
  fvhe = _onhe - 2*vb;
        /* the number of output verticess */
  _onv = 0;
  for ( i = 0; i < inv; i++ ) {
    j = imv[i].firsthalfedge;
    d = imv[i].degree;
    if ( vtag[i] )
      _onv += d + 2;
    else
      _onv += d;
  }
  *onv = _onv;

  ecn = pkv_GetScratchMemi ( inhe+1 );
  vcn = pkv_GetScratchMemi ( inv+1 );
  efn = pkv_GetScratchMemi ( inhe+1 );
  wlv = pkv_GetScratchMem ( inhe );
  wlf = pkv_GetScratchMem ( inhe );
  if ( !ecn || !vcn || !efn || !wlv || !wlf )
    goto failure;

        /* determine the number of output halfedges corresponding to */
        /* each input halfedge and compute the prefix sums */
  ecn[0] = 0;
  for ( i = 0; i < inhe; i++ )
    if ( imhe[i].otherhalf == -1 )
      ecn[i+1] = ecn[i] + 6;
    else
      ecn[i+1] = ecn[i] + 4;
        /* determine the number of output vertices corresponding to */
        /* each input vertex and compute the prefix sums; */
        /* then copy the point coordinates */
  vcn[0] = 0;
  for ( i = 0; i < inv; i++ ) {
    if ( vtag[i] )
      d = imv[i].degree + 2;
    else
      d = imv[i].degree;
    vcn[i+1] = vcn[i] + d;
    for ( j = 0; j < d; j++ ) {
      memcpy ( &ovc[(vcn[i]+j)*spdimen], &ivc[i*spdimen],
               spdimen*sizeof(double) );
      omv[vcn[i]+j].degree = 4;
    }
          /* all inner vertices are of degree 4 and boundary vertices are of */
          /* degree 2 */
    if ( vtag[i] )
      omv[vcn[i]].degree = omv[vcn[i]+d-1].degree = 2;
  }
  omv[0].firsthalfedge = 0;
  for ( i = 1; i < _onv; i++ )
    omv[i].firsthalfedge = omv[i-1].firsthalfedge + omv[i-1].degree;
        /* for each input halfedge find its number in the list */
        /* of halfedges of the vertex */
  for ( i = 0; i < inv; i++ ) {
    d = imv[i].degree;
    j = imv[i].firsthalfedge;
    for ( k = 0; k < d; k++ )
      wlv[imvhei[j+k]] = k;
  }
        /* and in the list of halfedges of the facet */
  for ( i = 0; i < infac; i++ ) {
    d = imfac[i].degree;
    j = imfac[i].firsthalfedge;
    for ( k = 0; k < d; k++ )
      wlf[imfhei[j+k]] = k;
  }
        /* for each input halfedge find the number of the corresponding */
        /* output facet */
  for ( i = 0, k = infac;  i < inhe; i++ ) {
    j = imhe[i].otherhalf;
    if ( j < 0 || j > i )
      efn[i] = k++;
    else
      efn[i] = efn[j];
  }
        /* the first infac facets correspond to the facets of the input mesh */
  memcpy ( omfac, imfac, infac*sizeof(BSMfacet) );
        /* the next ei + eb facets correspond to the edges of the input mesh */
  fvf = infac+ei+eb;
  for ( i = infac; i < fvf; i++ ) {
    omfac[i].degree = 4;
    omfac[i].firsthalfedge = omfac[i-1].firsthalfedge + omfac[i-1].degree;
  }
        /* the last inv facets correspond to the vertices of the input mesh */
  for ( i = fvf; i < _onfac; i++ ) {
    if ( vtag[i-fvf] )
      omfac[i].degree = imv[i-fvf].degree + 2;
    else
      omfac[i].degree = imv[i-fvf].degree;
    omfac[i].firsthalfedge = omfac[i-1].firsthalfedge + omfac[i-1].degree;
  }
        /* bind new halfedges in pairs */
  for ( i = k = 0; i < inhe; i++ ) {
    j = imhe[i].otherhalf;
    ecni = ecn[i];
    omhe[ecni].otherhalf = ecni+1;
    omhe[ecni+1].otherhalf = ecni;
    omhe[ecni+2].otherhalf = ecni+3;
    omhe[ecni+3].otherhalf = ecni+2;
    omhe[ecni].facetnum = imhe[i].facetnum;
    omhe[ecni+3].facetnum = fvf + imhe[i].v0;
    if ( j < 0 ) {
      omhe[ecni+4].otherhalf = omhe[ecni+5].otherhalf = -1;
      omhe[ecni+1].facetnum = omhe[ecni+2].facetnum =
      omhe[ecni+4].facetnum = omhe[ecni+5].facetnum = infac+k;
      k ++;
    }
    else if ( i < j ) {
      omhe[ecni+1].facetnum = omhe[ecni+2].facetnum = infac+k;
      k ++;
    }
    else
      omhe[ecni+1].facetnum = omhe[ecni+2].facetnum =
        omhe[omhe[ecn[j]].otherhalf].facetnum;
  }

        /* find halfedge sequences for output facets */
          /* corresponding to the input facets */
  for ( i = k = 0;  i < infac;  i++ ) {
    d = imfac[i].degree;
    for ( j = 0;  j < d;  j++, k++ )
      omfhei[k] = ecn[imfhei[imfac[i].firsthalfedge+j]];
  }
          /* corresponding to the input edges */
  for ( i = 0; i < inhe; i++ ) {
    k = omfac[efn[i]].firsthalfedge;
    ecni = ecn[i];
    if ( imhe[i].otherhalf < 0 ) {
      omfhei[k]   = ecni+1;
      omfhei[k+1] = ecni+2;
      omfhei[k+2] = ecni+4;
      omfhei[k+3] = ecni+5;
    }
    else {
      if ( imhe[i].otherhalf > i ) {
        omfhei[k] = ecni+1;
        omfhei[k+1] = ecni+2;
      }
      else {
        omfhei[k+2] = ecni+1;
        omfhei[k+3] = ecni+2;
      }
    }
  }
          /* corresponding to the input vertices */
  for ( i = m = 0; i < inv; i++ ) {
    d = imv[i].degree;
    j = imv[i].firsthalfedge;
    if ( vtag[i] ) {
      v0 = vcn[i];
      l = imvhei[imv[i].firsthalfedge];
      f = imhe[l].facetnum;
      l = PREVIFAC_HEDGE ( f, wlf[l] );
      ecnl = ecn[l];
      omhe[ecnl+4].v1 = v0;
      omhe[ecnl+5].otherhalf = p = fvhe + 2*m + 1;
      omhe[fvhe + 2*m + 1].otherhalf = ecnl+5;
      omhe[ecnl+5].v0 = omhe[fvhe + 2*m + 1].v1 = v0;
      omhe[ecnl+5].v1 = omhe[fvhe + 2*m + 1].v0 = v0+1;
      omvhei[omv[v0].firsthalfedge] = ecnl+5;
      omhe[fvhe+2*m].otherhalf = -1;
      omhe[fvhe+2*m].v0 = v0;
      omhe[fvhe+2*m].v1 = v0+d+1;
      omvhei[omv[v0].firsthalfedge+1] = fvhe + 2*m;
      omhe[fvhe+2*m].facetnum = omhe[fvhe+2*m+1].facetnum = fvf + i;
      omfhei[omfac[fvf+i].firsthalfedge] = fvhe + 2*m;
      for ( k = 0; k < d; k++ ) {
        v0 = vcn[i] + k + 1;
        l = imvhei[imv[i].firsthalfedge + k];
        f = imhe[l].facetnum;
        ecnl = ecn[l];
        omvhei[omv[v0].firsthalfedge] = p;
        omfhei[omfac[fvf+i].firsthalfedge+d+1-k] = p;
        omhe[ecnl].v0 = omhe[ecnl+1].v1 = v0;
        omhe[omhe[ecnl].otherhalf].v1 = v0;
        omhe[ecnl+2].v0 = omhe[ecnl+3].v1 = v0;
        omhe[ecnl+2].v1 = omhe[ecnl+3].v0 = v0+1;
        omvhei[omv[v0].firsthalfedge+2] = ecnl;
        omvhei[omv[v0].firsthalfedge+3] = ecnl+2;
        p = ecnl+3;
        l = PREVIFAC_HEDGE ( f, wlf[l] );
        ecnl = ecn[l];
        omhe[ecnl].v1 = v0;
        omhe[omhe[ecnl].otherhalf].v0 = v0;
        omvhei[omv[v0].firsthalfedge+1] = ecnl+1;
      }
      omfhei[omfac[fvf+i].firsthalfedge+1] = p;
      l = imvhei[imv[i].firsthalfedge+d-1];
      ecnl = ecn[l];
      omhe[ecnl+4].v0 = vcn[i] + d + 1;
      omvhei[omv[vcn[i]+d+1].firsthalfedge] = ecnl+3;
      omvhei[omv[vcn[i]+d+1].firsthalfedge+1] = ecnl+4;
      m ++;
    }
    else {
      for ( k = 0; k < d; k++ ) {
        v0 = vcn[i] + k;
        v1 = (k < d-1) ? v0+1 : vcn[i];
        l = imvhei[imv[i].firsthalfedge + k];
        f = imhe[l].facetnum;
        ecnl = ecn[l];
        omhe[ecnl].v0 = v0;
        omhe[ecnl+2].v0 = omhe[ecnl+3].v1 = v0;
        omhe[ecnl+2].v1 = omhe[ecnl+3].v0 = v1;
        omhe[omhe[ecnl].otherhalf].v1 = v0;
        omvhei[omv[v0].firsthalfedge] = ecnl;
        omvhei[omv[v0].firsthalfedge+1] = ecnl+2;
        omvhei[omv[v1].firsthalfedge+2] = ecnl+3;
        omfhei[omfac[fvf+i].firsthalfedge+d-1-k] = ecnl+3;
        l = PREVIFAC_HEDGE ( f, wlf[l] );
        ecnl = ecn[l];
        omhe[ecnl].v1 = v0;
        omhe[omhe[ecnl].otherhalf].v0 = v0;
        omvhei[omv[v0].firsthalfedge+3] = omhe[ecnl].otherhalf;
      }
    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef PREVIFAC_HEDGE
} /*bsm_Doublingd*/

