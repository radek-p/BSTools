
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
boolean bsm_RemoveFacetNum ( int inv, BSMvertex *imv, int *imvhei,
                             int inhe, BSMhalfedge *imhe,
                             int infac, BSMfacet *imfac, int *imfhei,
                             int nfr,
                             int *onv, int *onhe, int *onfac )
{
  void *sp;
  int  ei, eb, vi, vb;
  int  d, fhe, i, k, en, v0, fve, v0d;
  char *vtag, *ftag;

  sp = pkv_GetScratchMemTop ();
        /* compute the numbers of mesh vertices, halfedges and facets */
        /* after facet removal */
  if ( nfr < 0 || nfr >= infac ) {  /* invalid facet number */
    *onv = *onhe = *onfac = -1;
    return false;
  }

  vtag = pkv_GetScratchMem ( inv );
  ftag = pkv_GetScratchMem ( infac );
  if ( vtag && ftag ) {
    bsm_TagMesh ( inv, imv, imvhei, inhe, imhe, infac, imfac, imfhei,
                  vtag, ftag, &vi, &vb, &ei, &eb );
        /* the number of facets is decreased by one */
    *onfac = infac - 1;
        /* the number of halfedges is decreased by the number of facet edges */
    d = imfac[nfr].degree;
    *onhe = inhe - d;
        /* the number of vertices is decreased by the number of facet */
        /* vertices of degree 1 and increased by the number of vertices */
        /* to be split */
    fhe = imfac[nfr].firsthalfedge;
    for ( i = k = 0;  i < d;  i++ ) {
      en = imfhei[fhe+i];
      v0 = imhe[en].v0;
      if ( imv[v0].degree == 1 )
        k ++;
      else if ( vtag[v0] ) {
        v0d = imv[v0].degree;
        fve = imv[v0].firsthalfedge;
        if ( en != imvhei[fve] && en != imvhei[fve+v0d-1] )
          k --;  /* the vertex has to be split */
      }
    }
    *onv = inv - k;
    pkv_SetScratchMemTop ( sp );
    return true;
  }
  else {
    pkv_SetScratchMemTop ( sp );
    return false;
  }
} /*bsm_RemoveFacetNum*/

