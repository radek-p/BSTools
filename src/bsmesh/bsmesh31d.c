
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2013                                  */
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

/*#define LZZ*/

#ifdef LZZ
static int liczbaZmianZnaku ( int *tablica, int dlugoscTablicy )
{
  int pom, i;

  if ( tablica[0] == tablica[dlugoscTablicy-1] )
    pom = -1;
  else
    pom=0;
  for ( i = 0; i < dlugoscTablicy-1; i++ ) {
    if ( tablica[i] != tablica[i+1] )
      pom++;
  }
  if ( pom == -1 )
    pom = 0;
  return pom;
} /*liczbaZmianZnaku*/
#endif

boolean bsm_ExtractSubmeshVNum ( int inv,
                                 const BSMvertex *imv, const int *imvhei,
                                 int inhe, const BSMhalfedge *imhe,
                                 int infac, const BSMfacet *imfac, const int *imfhei,
                                 boolean *vtag,
                                 int *onv, int *onhe, int *onfac )
{
  void *sp;
  int  d, f, i, j, k, l;
  int *vertexCopiesCounter;
  int *isHalfedgeInResultingNet;
  int *isFacetInResultingNet;
  int *isVertexInResultingNet;
  int *tablicaPlusowIMinusow;
  int wielkoscomv, wielkoscomhe, wielkoscomfac;

  sp = pkv_GetScratchMemTop ();
  isHalfedgeInResultingNet = pkv_GetScratchMemi ( inhe );
  isFacetInResultingNet = pkv_GetScratchMemi ( infac );
  isVertexInResultingNet = pkv_GetScratchMemi ( inv );
  vertexCopiesCounter = pkv_GetScratchMemi ( inv );
  if ( !isHalfedgeInResultingNet || !isFacetInResultingNet ||
       !isVertexInResultingNet || !vertexCopiesCounter )
    goto failure;
  for ( i = 0; i < inhe; i++ )
    isHalfedgeInResultingNet[i] = -1;
  for ( i = 0; i < infac; i++ )
    isFacetInResultingNet[i] = -1;
  for ( i = 0; i < inv; i++ )
    isVertexInResultingNet[i] = -1;

  wielkoscomv = wielkoscomfac = wielkoscomhe = 0;
  for ( i = k = l = 0; i < infac; i++ ) {
    f = imfac[i].firsthalfedge;
    d = imfac[i].degree;
    for ( j = 0; j < d; j++ ) {
      if ( !vtag[imhe[imfhei[f+j]].v0] )
        goto skip;
    }
    isFacetInResultingNet[i] = k++;
    for ( j = 0; j < d; j++ ) {
      isHalfedgeInResultingNet[imfhei[f+j]] = l++;
      isVertexInResultingNet[imhe[imfhei[f+j]].v0] = 0;
    }
    wielkoscomfac ++;
    wielkoscomhe += d;
skip:
    continue;
  }
  for ( i = j = 0; i < inv; i++ ) {
    if ( isVertexInResultingNet[i] == 0 )
      isVertexInResultingNet[i] = j++;
  }
#ifndef LZZ
  for ( i = 0; i < inv; i++ ) {
    vertexCopiesCounter[i] = 0;
    if ( isVertexInResultingNet[i] >= 0 ) {
      f = imv[i].firsthalfedge;
      d = imv[i].degree;
/* tu jest blad; dla wierzcholka wewnetrznego, ktory pozostaje wewnetrzny, */
/* vertexCopiesCounter dostaje wartosc 0 zamiast 1. Trzeba najpierw sprawdzic, */
/* czy wierzcholek jest wewnetrzny, a nastepnie osobno policzyc zmiany znaku */
/* dla wewnetrznego i brzegowego */
      for ( j = 0; j < d; j++ ) {
        if ( isHalfedgeInResultingNet[imvhei[f+j]] >= 0 ) {
          if ( imhe[imvhei[f+j]].otherhalf==-1 ||
               isHalfedgeInResultingNet[imhe[imvhei[f+j]].otherhalf] < 0 )
            vertexCopiesCounter[i] ++;
        }
      }
      if ( vertexCopiesCounter[i] > 1 )
        goto failure;
      wielkoscomv ++;
    }
  }
#else
        /* tutaj kod z liczeniem plusow i minusow */
  for ( i = 0; i < inv; i++ ) {
            /* wybieram podciag polkrawedzi takich ktore wychodza z danego wierzcholka */
    f = imv[i].firsthalfedge;
    d = imv[i].degree;
    tablicaPlusowIMinusow = pkv_GetScratchMemi( d );
    for ( j = 0; j < d; j++ ) {
      if ( isHalfedgeInResultingNet[imvhei[f+j]] >= 0 )
        tablicaPlusowIMinusow[j] = 1;
      else
        tablicaPlusowIMinusow[j] = -1;
    }
    if ( liczbaZmianZnaku ( tablicaPlusowIMinusow, d ) <= 1 )
      goto failure;
    wielkoscomv ++;
  }
#endif

  *onv = wielkoscomv;
  *onhe = wielkoscomhe;
  *onfac = wielkoscomfac;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*bsm_ExtractSubmeshVNum*/

boolean bsm_ExtractSubmeshVd ( int spdimen, int inv,
                               const BSMvertex *imv, const int *imvhei,
                               double *iptc, int inhe, const BSMhalfedge *imhe,
                               int infac, const BSMfacet *imfac, const int *imfhei,
                               boolean *vtag,
                               int *onv,  BSMvertex *omv,  int *omvhei,
                               double *optc, int *onhe,  BSMhalfedge *omhe,
                               int *onfac,  BSMfacet *omfac,  int *omfhei )
{
  void *sp;
  int  d, f, i, j, k, l, m, ifhe, ofhe;
  int *vertexCopiesCounter;
  int *isHalfedgeInResultingNet;
  int *isFacetInResultingNet;
  int *isVertexInResultingNet;
  int *tablicaPlusowIMinusow;
  int wielkoscomv, wielkoscomhe, wielkoscomfac;

  sp = pkv_GetScratchMemTop ();
  isHalfedgeInResultingNet = pkv_GetScratchMemi ( inhe );
  isFacetInResultingNet = pkv_GetScratchMemi ( infac );
  isVertexInResultingNet = pkv_GetScratchMemi ( inv );
  vertexCopiesCounter = pkv_GetScratchMemi ( inv );
  if ( !isHalfedgeInResultingNet || !isFacetInResultingNet ||
       !isVertexInResultingNet || !vertexCopiesCounter )
    goto failure;
  for ( i = 0; i < inhe; i++ )
    isHalfedgeInResultingNet[i] = -1;
  for ( i = 0; i < infac; i++ )
    isFacetInResultingNet[i] = -1;
  for ( i = 0; i < inv; i++ )
    isVertexInResultingNet[i] = -1;

  wielkoscomv = wielkoscomfac = wielkoscomhe = 0;
  for ( i = k = l = 0; i < infac; i++ ) {
    f = imfac[i].firsthalfedge;
    d = imfac[i].degree;
    for ( j = 0; j < d; j++ ) {
      if ( !vtag[imhe[imfhei[f+j]].v0] )
        goto skip;
    }
    isFacetInResultingNet[i] = k++;
    for ( j = 0; j < d; j++ ) {
      isHalfedgeInResultingNet[imfhei[f+j]] = l++;
      isVertexInResultingNet[imhe[imfhei[f+j]].v0] = 0;
    }
    wielkoscomfac ++;
    wielkoscomhe += d;
skip:
    continue;
  }
  for ( i = j = 0; i < inv; i++ ) {
    if ( isVertexInResultingNet[i] == 0 )
      isVertexInResultingNet[i] = j++;
  }
#ifndef LZZ
  for ( i = 0; i < inv; i++ ) {
    vertexCopiesCounter[i] = 0;
    if ( isVertexInResultingNet[i] >= 0 ) {
      f = imv[i].firsthalfedge;
      d = imv[i].degree;
/* tu jest blad; dla wierzcholka wewnetrznego, ktory pozostaje wewnetrzny, */
/* vertexCopiesCounter dostaje wartosc 0 zamiast 1. Trzeba najpierw sprawdzic, */
/* czy wierzcholek jest wewnetrzny, a nastepnie osobno policzyc zmiany znaku */
/* dla wewnetrznego i brzegowego */
      for ( j = 0; j < d; j++ ) {
        if ( isHalfedgeInResultingNet[imvhei[f+j]] >= 0 ) {
          if ( imhe[imvhei[f+j]].otherhalf==-1 ||
               isHalfedgeInResultingNet[imhe[imvhei[f+j]].otherhalf] < 0 )
            vertexCopiesCounter[i] ++;
        }
      }
      if ( vertexCopiesCounter[i] > 1 )
        goto failure;
      wielkoscomv ++;
    }
  }
#else
        /* tutaj kod z liczeniem plusow i minusow */
  for ( i = 0; i < inv; i++ ) {
            /* wybieram podciag polkrawedzi takich ktore wychodza z danego wierzcholka */
    f = imv[i].firsthalfedge;
    d = imv[i].degree;
    tablicaPlusowIMinusow = pkv_GetScratchMemi( d );
    for ( j = 0; j < d; j++ ) {
      if ( isHalfedgeInResultingNet[imvhei[f+j]] >= 0 )
        tablicaPlusowIMinusow[j] = 1;
      else
        tablicaPlusowIMinusow[j] = -1;
    }
    if ( liczbaZmianZnaku ( tablicaPlusowIMinusow, d ) <= 1 )
      goto failure;
    wielkoscomv ++;
  }
#endif

/* wpisuje wierzcholki, w odpowiednim momencie je podwajajac */

  /* przepisanie elementow nowej siatki */
    /* omv, optc, omvhei */
  for ( i = j = 0; i < inv; i++ ) {
    j = isVertexInResultingNet[i];
    if ( j >= 0 ) {
      memcpy ( optc+j*spdimen, iptc+i*spdimen, spdimen*sizeof(double) );
      ifhe = imv[i].firsthalfedge;
      omv[j].degree = 0;
      for ( k = 0; k < imv[i].degree; k++ ) {
        if ( isHalfedgeInResultingNet[imvhei[ifhe+k]] >= 0 )
          omv[j].degree ++;
      }
      omv[j].firsthalfedge = ofhe =
                  j ? omv[j-1].degree + omv[j-1].firsthalfedge : 0;
/* tu tez jest delikatny moment. Jesli stopien wierzcholka sie zmienil, to */
/* znaczy, ze jakies jgo polkrawedzie wypadly. Wtedy ostatnia wychodzaca */
/* polkrawedz ma byc bez pary, wiec przepisywac na ogol trzeba "od srodka", */
/* nie od poczatku. Proponuje rozwiazanie, ktore w wersji dopuszczajacej */
/* powielanie wierzcholkow trzeba bedzie zmodyfikowac. */
      d = imv[i].degree;
      if ( omv[j].degree < d ) {
        m = d;
        while ( isHalfedgeInResultingNet[imvhei[ifhe+m-1]] >= 0 )
          m--;
      }
      else
        m = 0;
      for ( k = l = 0; k < d; k++ ) {
        if ( m+k == d )
          m -= d;
        if ( isHalfedgeInResultingNet[imvhei[ifhe+m+k]] >= 0 ) {
          omvhei[ofhe+l] = isHalfedgeInResultingNet[imvhei[ifhe+m+k]];
          l ++;
        }
      }
    }
  }

    /* omf, omfhei */
  for ( i = j = 0; i < infac; i++ ) {
    j = isFacetInResultingNet[i];
    if ( j >= 0 ) {
      d = omfac[j].degree = imfac[i].degree;
      ifhe = imfac[i].firsthalfedge;
      omfac[j].firsthalfedge = ofhe =
                  j ? omfac[j-1].degree + omfac[j-1].firsthalfedge : 0;
      for ( k = l = 0;  k < imfac[i].degree;  k++, l++ )
        omfhei[ofhe+l] = isHalfedgeInResultingNet[imfhei[ifhe+k]];
    }
  }

    /* omhe */
  for ( i = 0; i < inhe; i++ ) {
    j = isHalfedgeInResultingNet[i];
    if ( j >= 0 ) {
      omhe[j].facetnum = isFacetInResultingNet[imhe[i].facetnum];
      k = imhe[i].otherhalf;
      omhe[j].otherhalf = k >= 0 ? isHalfedgeInResultingNet[k] : -1;
      omhe[j].v0 = isVertexInResultingNet[imhe[i].v0];
      omhe[j].v1 = isVertexInResultingNet[imhe[i].v1];
    }
  }

  *onv = wielkoscomv;
  *onhe = wielkoscomhe;
  *onfac = wielkoscomfac;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*bsm_ExtractSubmeshVd*/

