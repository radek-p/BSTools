
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2015                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */
/* this file was written by Anna Sierhej                                     */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "bsmesh.h"

#include "bsmprivate.h"

boolean _bsm_DivideFacetd ( int spdimen, int *inv,
                            BSMvertex *imv, int *imvhei,
                            double *iptc, int *inhe, BSMhalfedge *imhe,
                            int *infac, BSMfacet *imfac, int *imfhei,
                            int nV0, int nV1 )
{
  void    *sp;
  int     i, j, d, fhe0, d0, fhe1, d1, nf, lastFacetStart;
  boolean rightFace;
  int     face2d;
  int     *tmfhei, *tmvhei;

  sp = pkv_GetScratchMemTop ();
  if ( nV0 < 0 || nV0 >= *inv || nV1 < 0 || nV1 >= *inv || nV0 == nV1 )
    goto failure;

  tmfhei = pkv_GetScratchMemi ( *inhe );
  tmvhei = pkv_GetScratchMemi ( *inhe );
  if ( !tmfhei || !tmvhei )
    goto failure;
  memcpy ( tmfhei, imfhei, *inhe*sizeof(int) );
  memcpy ( tmvhei, imvhei, *inhe*sizeof(int) );

  if ( nV0 > nV1 ) {
    i = nV0;
    nV0 = nV1;
    nV1 = i;
  }

  face2d = 0;
  nf = -1;

  /* zorientowanie sie w ktorej scianie ma lezec polkrawedz */
  fhe0 = imv[nV0].firsthalfedge;
  d0 = imv[nV0].degree;
  fhe1 = imv[nV1].firsthalfedge;
  d1 = imv[nV1].degree;
  for ( i = fhe0; i < fhe0+d0; i++ ) {
    for ( j = fhe1; j < fhe1+d1; j++ ) {
    /* jesli wierzcholki leza na tej samej halfedge - failure */
      if ( imhe[imvhei[j]].v0 == imhe[imvhei[i]].v1 ||
           imhe[imvhei[i]].v0 == imhe[imvhei[j]].v1)
        goto failure;

      if ( imhe[imvhei[j]].facetnum == imhe[imvhei[i]].facetnum ) {
        nf = imhe[imvhei[j]].facetnum;
        goto facet_found;
      }
    }
  }
  goto failure;

facet_found:
  /* wiemy juz ktora sciane bedziemy dzielic. teraz trzeba zdecydowac
     jak bedziemy ja dzielic. przyjmuje zalozenie, ze z V0 bedzie
     wychodzila polkrawedz nalezaca do starej sciany, a z V1 polkrawedz
     nalezaca do nowej sciany. przyjmuje tez, ze n-1 - sza polkrawedz w tablicy
     omhe bedzie z V0 do V1, a ostatnia - z V1 do V0 */

  /*przepisanie punktow*/
  *inhe += 2;
  *infac += 1;

  /*omf*/
  d = imfac[nf].degree;
  imfhei[imfac[nf].firsthalfedge] = *inhe-2;
  i = 0;
  while ( imhe[tmfhei[imfac[nf].firsthalfedge + i]].v0 != nV1 )
    i++;

  imfac[nf].degree = 1;
  j = 1;
  do {
    imfhei[imfac[nf].firsthalfedge+j] = tmfhei[imfac[nf].firsthalfedge + i];
    i++;
    j++;
    i = i % d;
    imfac[nf].degree ++;
  } while ( imhe[tmfhei[imfac[nf].firsthalfedge + i]].v0 != nV0 );

  lastFacetStart = i;
  face2d = d-imfac[nf].degree+1;

  for ( i = imfac[nf].degree+imfac[nf].firsthalfedge; i < *inhe-2-face2d+1; i++ )
    imfhei[i] = tmfhei[i+face2d-1];
  for ( i = nf+1; i < *infac-1; i++ )
    imfac[i].firsthalfedge = imfac[i].firsthalfedge-face2d+1;

  i = *infac-1;
  imfac[i].firsthalfedge = *inhe-2-face2d+1;
  imfac[i].degree = face2d+1;
  j = 0;
  i = lastFacetStart;

  do {
    imfhei[imfac[*infac-1].firsthalfedge+j] = tmfhei[imfac[nf].firsthalfedge+i];
    i++;
    j++;
    i = i%d;
  } while ( imhe[tmfhei[imfac[nf].firsthalfedge+i]].v0 != nV1 );

  imfhei[*inhe-1] = *inhe-1;

  /*omv*/
  rightFace = false;

  imv[nV0].degree ++;
  i = imv[nV0].firsthalfedge;
  while ( !rightFace ) {
    if ( imhe[tmvhei[i]].facetnum == nf ) {
      rightFace = true;
      imvhei[i] = *inhe-2;
      i ++;
      imvhei[i] = tmvhei[i-1];
    }
    i ++;
  }

  while ( i != imv[nV0].firsthalfedge+imv[nV0].degree) {
    imvhei[i] = tmvhei[i-1];
    i ++;
  }

  for ( i = nV0+1; i < nV1; i++ ) {
    imv[i].firsthalfedge ++;
    for ( j = imv[i].firsthalfedge; j < imv[i].firsthalfedge+imv[i].degree; j++ )
      imvhei[j] = tmvhei[j-1];
  }

  /*omvhei*/
  rightFace = false;

  imv[nV1].degree ++;
  imv[nV1].firsthalfedge ++;
  i = imv[nV1].firsthalfedge;

  while ( !rightFace ) {
    if ( imhe[tmvhei[i-1]].facetnum == nf ) {
      rightFace = true;
      imvhei[i] = *inhe-1;
      i++;
      imvhei[i] = tmvhei[i-2];
    }
    else
      imvhei[i] = tmvhei[i-1];
    i ++;
  }

  while ( i != imv[nV1].firsthalfedge+imv[nV1].degree ) {
    imvhei[i] = tmvhei[i-2];
    i ++;
  }

  for ( i = nV1+1; i < *inv; i++ ) {
    imv[i].firsthalfedge += 2;
    for ( j = imv[i].firsthalfedge; j < imv[i].firsthalfedge+imv[i].degree; j++ )
      imvhei[j] = tmvhei[j-2];
  }

/*omhe*/
  imhe[*inhe-2].facetnum = nf;
  imhe[*inhe-2].otherhalf = *inhe-1;
  imhe[*inhe-2].v0 = nV0;
  imhe[*inhe-2].v1 = nV1;

  imhe[*inhe-1].facetnum = *infac-1;
  imhe[*inhe-1].otherhalf = *inhe-2;
  imhe[*inhe-1].v0 = nV1;
  imhe[*inhe-1].v1 = nV0;

  for ( i = imfac[*infac-1].firsthalfedge;
        i < imfac[*infac-1].degree+imfac[*infac-1].firsthalfedge;
        i++ )
    imhe[imfhei[i]].facetnum = *infac-1;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_bsm_DivideFacetd*/

