
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */
/* this file was written by Anna Sierhej                                     */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "pkvaria.h"
#include "pknum.h"
#include "bsmesh.h"
#include "bsmprivate.h"

boolean bsm_DivideFacetd ( int spdimen, int inv,
                           const BSMvertex *imv, const int *imvhei,
                           double *iptc, int inhe, const BSMhalfedge *imhe,
                           int infac, const BSMfacet *imfac, const int *imfhei,
                           int nV0, int nV1,
                           int *onv, BSMvertex *omv, int *omvhei,
                           double *optc, int* onhe, BSMhalfedge *omhe,
                           int *onfac, BSMfacet *omfac, int *omfhei )
{
  int       i, j, d, fhe0, d0, fhe1, d1, nf, lastFacetStart;
  boolean   rightFace;
  int       face2d;

  if ( nV0 < 0 || nV0 >= inv || nV1 < 0 || nV1 >= inv || nV0 == nV1 )
    goto failure;

  if(nV0>nV1){
    i=nV0;
    nV0=nV1;
    nV1=i;
  }

  *onv = *onfac = *onhe = 0;
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
  *onhe = inhe+2;
  *onfac = infac+1;
  *onv = inv;

    /*optc*/
  memcpy ( optc, iptc, inv*spdimen*sizeof(double) );
    /*omf*/

  for ( i = 0; i < nf; i++ ) {
    omfac[i] = imfac[i];
    for ( j = omfac[i].firsthalfedge; j < omfac[i].firsthalfedge+omfac[i].degree; j++ )
      omfhei[j] = imfhei[j];
  }

 omfac[nf]=imfac[nf];

    d=omfac[nf].degree;
    omfhei[omfac[nf].firsthalfedge]=inhe;

    i=0;
    while(imhe[imfhei[imfac[nf].firsthalfedge + i]].v0!=nV1){
        i++;
    }

    omfac[nf].degree=1;
    j=1;
    do{
        omfhei[omfac[nf].firsthalfedge+j]=imfhei[imfac[nf].firsthalfedge + i];
        i++;
        j++;
        i=i%d;
        omfac[nf].degree++;
    }while(imhe[imfhei[imfac[nf].firsthalfedge + i]].v0!=nV0);

    lastFacetStart=i;
    face2d=d-omfac[nf].degree+1;

    for(i=omfac[nf].degree+omfac[nf].firsthalfedge; i<inhe-face2d+1;i++){
        omfhei[i]=imfhei[i+face2d-1];
    }

    for(i=nf+1; i<infac; i++){
        omfac[i]=imfac[i];
        omfac[i].firsthalfedge=omfac[i].firsthalfedge-face2d+1;
    }

    i=infac;
    omfac[i].firsthalfedge=inhe-face2d+1;
    omfac[i].degree=face2d+1;
    j=0;
    i=lastFacetStart;

    do{
        omfhei[omfac[infac].firsthalfedge+j]=imfhei[imfac[nf].firsthalfedge + i];
        i++;
        j++;
        i=i%d;
    }while(imhe[imfhei[imfac[nf].firsthalfedge + i]].v0!=nV1);

    omfhei[inhe+1]=inhe+1;

    /*omhe*/
  for ( i = 0; i < inhe; i++ )
    omhe[i] = imhe[i];

  omhe[inhe].facetnum = nf;
  omhe[inhe].otherhalf = inhe+1;
  omhe[inhe].v0 = nV0;
  omhe[inhe].v1 = nV1;

  omhe[inhe+1].facetnum = infac;
  omhe[inhe+1].otherhalf = inhe;
  omhe[inhe+1].v0 = nV1;
  omhe[inhe+1].v1 = nV0;

  for ( i = omfac[infac].firsthalfedge;
        i < omfac[infac].degree+omfac[infac].firsthalfedge;
        i++ )
    omhe[omfhei[i]].facetnum = infac;

    /*omv*/
  for ( i = 0; i < nV0; i++ ) {
    omv[i] = imv[i];
    for ( j = omv[i].firsthalfedge; j < omv[i].firsthalfedge+omv[i].degree; j++ )
      omvhei[j]=imvhei[j];
  }

  rightFace = false;
  omv[nV0] = imv[nV0];
  omv[nV0].degree++;
  i = omv[nV0].firsthalfedge;

  while ( !rightFace ) {
    if ( imhe[imvhei[i]].facetnum == nf ) {
      rightFace = true;
      omvhei[i] = inhe;
      i++;
      omvhei[i] = imvhei[i-1];
      i++;
    }
    else {
      omvhei[i] = imvhei[i];
      i++;
    }
  }

  while ( i != omv[nV0].firsthalfedge+omv[nV0].degree ) {
    omvhei[i] = imvhei[i-1];
    i++;
  }

  for ( i = nV0+1; i <= nV1; i++ ) {
    omv[i] = imv[i];
    omv[i].firsthalfedge ++;
    for ( j = omv[i].firsthalfedge; j < omv[i].firsthalfedge+omv[i].degree; j++ )
      omvhei[j] = imvhei[j-1];
  }

    /*omvhei*/
  rightFace = false;

  omv[nV1] = imv[nV1];
  omv[nV1].degree ++;
  omv[nV1].firsthalfedge ++;
  i = omv[nV1].firsthalfedge;

  while ( !rightFace ) {
    if ( imhe[imvhei[i-1]].facetnum == nf ) {
      rightFace = true;
      omvhei[i] = inhe+1;
      i++;
      omvhei[i] = imvhei[i-2];
      i++;
    }
    else {
      omvhei[i] = imvhei[i-1];
      i++;
    }
  }

  while ( i != omv[nV1].firsthalfedge+omv[nV1].degree ) {
    omvhei[i] = imvhei[i-2];
    i++;
  }

  for ( i = nV1+1; i < inv; i++ ) {
    omv[i] = imv[i];
    omv[i].firsthalfedge += 2;
    for ( j = omv[i].firsthalfedge; j < omv[i].firsthalfedge+omv[i].degree; j++ )
      omvhei[j] = imvhei[j-2];
  }
  return true;

failure:
  return false;
} /*bsm_DivideFacetd*/

