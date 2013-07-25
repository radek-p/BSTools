
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

/*#define DEBUG*/
#ifdef DEBUG
#define ERROR(n) { printf ( "Averaging error %d\n", n ); goto failure; }
#endif

/* This procedure uses the same averaging algorithm that the              */
/* bsm_AveragingMatd procedure, which is intended to obtain the same      */
/* correspondence between the vertices of the input and output mesh - the */
/* only difference is that it averages the vertices instead of producing  */
/* the averaging matrix. Should any modifications be necessary, te same   */
/* have to be made in the code of bsm_AveragingMatd.                      */

boolean bsm_Averagingd ( int spdimen,
                         int inv, BSMvertex *imv, int *imvhei, double *ivc,
                         int inhe, BSMhalfedge *imhe,
                         int infac, BSMfacet *imfac, int *imfhei,
                         int *onv, BSMvertex *omv, int *omvhei, double *ovc,
                         int *onhe, BSMhalfedge *omhe,
                         int *onfac, BSMfacet *omfac, int *omfhei )
{
  void   *sp;
  int    vi, vb, ei, eb, _onv, _onhe, _onfac, fvn;
  int    *nhei, *nfi, *nvi, *fvnum;
  int    i, j, k, l, m, n, d, e, r, s, t, v0, v1, fhe;
  double *av;
  char   *vtag, *ftag;

  sp = pkv_GetScratchMemTop ();

  nvi = pkv_GetScratchMemi ( infac );
  nhei = pkv_GetScratchMemi ( inhe );
  nfi = pkv_GetScratchMemi ( inv );
  fvnum = pkv_GetScratchMemi ( infac );
  vtag = pkv_GetScratchMem ( inv );
  ftag = pkv_GetScratchMem ( infac );
  if ( !nvi || !nhei || !nfi || !fvnum || !vtag || !ftag )
    goto failure;

  bsm_TagMesh ( inv, imv, imvhei, inhe, imhe, infac, imfac, imfhei,
                vtag, ftag, &vi, &vb, &ei, &eb );
  if ( vb == 0 ) { /* there are no boundary vertices and edges */
    *onv = _onv = infac;
    *onhe = _onhe = inhe;
    *onfac = _onfac = inv;
    for ( i = 0; i < infac; i++ ) {
      nvi[i] = i;
      fvnum[i] = 1;
    }
    for ( i = 0; i < inhe; i++ )
      nhei[i] = i;
    for ( i = 0; i < inv; i++ )
      nfi[i] = i;
  }
  else {
        /* number of output vertices */
    for ( i = k = _onv = 0;  i < infac;  i++ ) {
      fvnum[i] = fvn = _bsm_AveragingFVN ( inv, imv, imvhei, inhe, imhe,
                                           infac, imfac, imfhei, i, vtag );
      if ( fvn ) {
        _onv += fvn;
        nvi[i] = k;
        k += fvn;
      }
      else
        nvi[i] = -1;
    }
    *onv = _onv;
        /* number of output facets */
    for ( i = k = 0;  i < inv;  i++ )
      if ( !vtag[i] )
        nfi[i] = k++;
      else
        nfi[i] = -1;
    *onfac = _onfac = vi;
        /* number of output halfedges */
    for ( i = _onhe = 0;  i < inhe;  i++ )
      if ( imhe[i].otherhalf >= 0 && !vtag[imhe[i].v1] )
        nhei[i] = _onhe ++;
      else
        nhei[i] = -1;
    *onhe = _onhe;
  }
#ifdef DEBUG
memset ( omv, 0, _onv*sizeof(BSMvertex) );
memset ( omhe, 0, _onhe*sizeof(BSMhalfedge) );
memset ( omfac, 0, _onfac*sizeof(BSMfacet) );
memset ( omvhei, 0, _onhe*sizeof(int) );
memset ( omfhei, 0, _onhe*sizeof(int) );
printf ( "*\n" );
#endif

        /* the number of edges of each output facet is the degree */
        /* of the corresponding inner input vertex */
  for ( i = k = j = 0;  i < inv;  i++ ) {
    k = nfi[i];
    if ( k >= 0 ) {
#ifdef DEBUG
if ( k < 0 || k >= _onfac )
  ERROR ( 1 )
#endif
      d = imv[i].degree;
      omfac[k].degree = d;
      omfac[k].firsthalfedge = j;
      j += d;
    }
  }
        /* generate the vertices for each facet */
  for ( i = j = 0; i < infac; i++ )
    if ( fvnum[i] > 0 ) {
      n = nvi[i];
      d = imfac[i].degree;
      fhe = imfac[i].firsthalfedge;
      r = fvnum[i];
      for ( k = 0; k < d; k++ ) {
        v1 = imhe[imfhei[fhe+k]].v1;
        if ( vtag[v1] )
          break;
      }
        /* generate the vertices for the facet */
      for ( s = 0;  s < r;  s++, n++ ) {
          /* find the first inner halfedge */
        do {
          k = k >= d-1 ? 0 : k+1;
          v1 = imhe[imfhei[fhe+k]].v1;
        } while ( vtag[v1] );
          /* determine the vertex degree */
        for ( m = 0, t = (k+1) % d;  m < d;  m++, t = (t+1) % d ) {
          v0 = imhe[imfhei[fhe+t]].v0;
          if ( vtag[v0] )
            break;
        }
          /* create the vertex description */
#ifdef DEBUG
if ( n < 0 || n > _onv )
  ERROR ( 2 )
#endif
        omv[n].degree = m;
        omv[n].firsthalfedge = j;
        for ( l = m-1, t = k;  l >= 0;  l--, t = (t+1) % d ) {
          v1 = imhe[imfhei[fhe+t]].v1;
#ifdef DEBUG
if ( j+l < 0 || j+l > _onhe )
  ERROR ( 3 )
#endif
          omvhei[j+l] = e = nhei[imfhei[fhe+t]];
#ifdef DEBUG
if ( e < 0 || e > _onhe )
  ERROR ( 4 )
#endif
          omhe[e].v0 = n;
          omhe[e].facetnum = nfi[v1];
        }
        j += m;
      }
    }
        /* bind the output halfedges in pairs */
  for ( i = 0;  i < inhe;  i++ ) {
    k = nhei[i];
    if ( k >= 0 )
      omhe[k].otherhalf = nhei[imhe[i].otherhalf];
  }
        /* for output facets find their halfedge sequences */
  for ( i = 0;  i < inv;  i++ ) {
    k = nfi[i];
    if ( k >= 0 ) {
      d = imv[i].degree;
      j = imv[i].firsthalfedge;
      l = omfac[k].firsthalfedge;
      for ( m = 0; m < d; m++ )
        omfhei[l+m] = nhei[imhe[imvhei[j+d-1-m]].otherhalf];
          /* find the end vertex of each halfedge */
      for ( m = d-1, v1 = omhe[omfhei[l]].v0;
            m >= 0;
            v1 = omhe[omfhei[l+m]].v0, m-- )
        omhe[omfhei[l+m]].v1 = v1;
    }
  }
        /* the numerical computation - actual averaging of vertices */
  for ( i = 0; i < infac; i++ )
    if ( fvnum[i] > 0 ) {
      n = nvi[i];
      av = &ovc[n*spdimen];
      d = imfac[i].degree;
      j = imfac[i].firsthalfedge;
      pkn_AddMatrixd ( 1, spdimen, 0, &ivc[imhe[imfhei[j]].v0*spdimen],
                       0, &ivc[imhe[imfhei[j+1]].v0*spdimen], 0, av );
      for ( k = 2; k < d; k++ )
        pkn_AddMatrixd ( 1, spdimen, 0, av,
                         0, &ivc[imhe[imfhei[j+k]].v0*spdimen], 0, av );
      pkn_MultMatrixNumd ( 1, spdimen, 0, av, 1.0/(double)d,
                           0, av );
      for ( s = 1; s < fvnum[i]; s++ )
        memcpy ( &ovc[(n+s)*spdimen], av, spdimen*sizeof(double) );
    }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*bsm_Averagingd*/

