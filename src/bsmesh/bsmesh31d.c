
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2013, 2014                            */
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

/*#define DEBUG*/

boolean bsm_ExtractSubmeshVNum ( int inv, const BSMvertex *imv, const int *imvhei,
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
  int wielkoscomv, wielkoscomhe, wielkoscomfac;

  sp = pkv_GetScratchMemTop ();
  isHalfedgeInResultingNet = pkv_GetScratchMemi ( inhe );
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
      if ( vertexCopiesCounter[i] == 0 )
        vertexCopiesCounter[i]=1;
      wielkoscomv+=vertexCopiesCounter[i];
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
} /*bsm_ExtractSubmeshVNum*/

boolean bsm_ExtractSubmeshVd ( int spdimen,
                               int inv, const BSMvertex *imv, const int *imvhei,
                               double *iptc, int inhe, const BSMhalfedge *imhe,
                               int infac, const BSMfacet *imfac, const int *imfhei,
                               boolean *vtag,
                               int *onv, BSMvertex *omv, int *omvhei,
                               double *optc, int *onhe, BSMhalfedge *omhe,
                               int *onfac, BSMfacet *omfac, int *omfhei )
{
  void        *sp;
  int         d, f, i, j, k, l, omvLastFreeSlot, ifhe, ofhe, face,
              currentVertex, heIter, omvheiIter, boundaryHalfedge, omvheiEntry;
  int         *vertexCopiesCounter;
  int         *isHalfedgeInResultingNet;
  int         *isFacetInResultingNet;
  int         *isVertexInResultingNet;
  int         *V0, *V1;
  boolean     isSecondOrNext;
  BSMhalfedge heIn, heOut;
  int         zDrugiejStrony, a, pom;

    /*przypisanie pamieci na tablice*/
  sp = pkv_GetScratchMemTop ();
  V1 = pkv_GetScratchMemi ( inhe );
  V0 = pkv_GetScratchMemi ( inhe );
  isHalfedgeInResultingNet = pkv_GetScratchMemi ( inhe );
  isFacetInResultingNet = pkv_GetScratchMemi ( infac );
  isVertexInResultingNet = pkv_GetScratchMemi ( inv );
  vertexCopiesCounter = pkv_GetScratchMemi ( inv );
  if ( !isHalfedgeInResultingNet || !isFacetInResultingNet ||
       !isVertexInResultingNet || !vertexCopiesCounter )
    goto failure;

  if ( !V0 || !V1 )
    goto failure;

    /*zainicjowanie tablic i zmiennych*/
  for ( i = 0; i < inhe; i++ ) {
    V0[i] = -10;
    V1[i] = -10;
  }
  for ( i = 0; i < *onv; i++ )
    omv[i].degree=0;
  for ( i = 0; i < inhe; i++ )
    isHalfedgeInResultingNet[i] = -1;
  for ( i = 0; i < infac; i++ )
    isFacetInResultingNet[i] = -1;
  for ( i = 0; i < inv; i++ )
    isVertexInResultingNet[i] = -1;

  *onv = *onfac = *onhe = 0;

    /*oznaczenie scian wierzcholkow i polkrawedzi wchodzacych do siatki wynikowej*/
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
    (*onfac) ++;
    (*onhe) += d;
skip:
    continue;
  }

    /*przenumerowanie wierzcholkow wchodzacych do siatki wynikowej*/
  for ( i = j = 0; i < inv; i++ ) {
    if ( isVertexInResultingNet[i] == 0 )
      isVertexInResultingNet[i] = j++;
  }
    /*liczenie kopii wierzcholkow*/
  for ( i = 0; i < inv; i++ ) {
    vertexCopiesCounter[i] = -1;
    if ( isVertexInResultingNet[i] >= 0 ) {
      vertexCopiesCounter[i]++;
      f = imv[i].firsthalfedge;
      d = imv[i].degree;
      for ( j = 0; j < d; j++ ) {
        if ( isHalfedgeInResultingNet[imvhei[f+j]] >= 0 ) {
          if ( imhe[imvhei[f+j]].otherhalf==-1 ||
               isHalfedgeInResultingNet[imhe[imvhei[f+j]].otherhalf] < 0 )
            vertexCopiesCounter[i] ++;
        }
      }
      if ( vertexCopiesCounter[i] == 0 )
        (*onv) ++;
      (*onv) += vertexCopiesCounter[i];
    }
  }

#ifdef DEBUG
printf ( "VertexCopiesCounter\n" );
for ( i = 0; i < inv; i++ ) {
  printf ( "%d", vertexCopiesCounter[i] );
  printf ( " " );
}
#endif

    /*wierzcholki podwojone bede numerowac 'od tylu'*/
  omvLastFreeSlot=*onv;
  omvheiIter=0;
  for ( i = 0; i < inv; i++ ) {
    if ( vertexCopiesCounter[i] >= 1 ) {
      isSecondOrNext = false;
      f = imv[i].firsthalfedge;
      d = imv[i].degree;
      for ( j = 0; j < d; j++ ) {
        if ( isHalfedgeInResultingNet[imvhei[f+j]] >= 0 &&
             ( imhe[imvhei[f+j]].otherhalf == -1 ||
               isHalfedgeInResultingNet[imhe[imvhei[f+j]].otherhalf] < 0 ) ) {
          if ( !isSecondOrNext ) {
            isSecondOrNext = true;
            currentVertex = isVertexInResultingNet[i];
          }
          else {
            omvLastFreeSlot --;
            currentVertex = omvLastFreeSlot;
          }
    /*iteruje po polkrawedziach taka jakby harmonijka, po drodze zmieniajac odpowiednie wartosci*/
          omv[currentVertex].degree ++;
          boundaryHalfedge = omvheiIter;
          omvheiEntry = isHalfedgeInResultingNet[imvhei[f+j]];
          heOut = imhe[imvhei[f+j]];
          V0[imvhei[f+j]] = currentVertex;
          face = heOut.facetnum;
          heIn.otherhalf = -1;
          for ( heIter = 0; heIter < imfac[face].degree; heIter++ ) {
            if ( imhe[imfhei[imfac[face].firsthalfedge+heIter]].v1 == i ) {
              heIn = imhe[imfhei[imfac[face].firsthalfedge+heIter]];
              V1[imfhei[imfac[face].firsthalfedge+heIter]] = currentVertex;
            }
          }
          if ( heIn.otherhalf >= 0 ) {
            V0[heIn.otherhalf] = currentVertex;
            heOut = imhe[heIn.otherhalf];
            face = heOut.facetnum;
          }
          while ( heIn.otherhalf >= 0 && isFacetInResultingNet[face] >= 0 ) {
            omv[currentVertex].degree ++;
            omvhei[omvheiIter] = isHalfedgeInResultingNet[heIn.otherhalf];
            omvheiIter ++;
            for ( heIter = 0; heIter < imfac[face].degree; heIter++ ) {
              if ( imhe[imfhei[imfac[face].firsthalfedge+heIter]].v1 == i ) {
                heIn = imhe[imfhei[imfac[face].firsthalfedge+heIter]];
                V1[imfhei[imfac[face].firsthalfedge+heIter]] = currentVertex;
              }
            }
            if ( heIn.otherhalf < 0 )
              break;
            heOut = imhe[heIn.otherhalf];
            V0[heIn.otherhalf] = currentVertex;
            face = heOut.facetnum;
          }
          omv[currentVertex].firsthalfedge = boundaryHalfedge;
          omvhei[omvheiIter] = omvheiEntry;
          zDrugiejStrony = 2;
          for ( a = omv[currentVertex].firsthalfedge;
                a < omv[currentVertex].firsthalfedge+omv[currentVertex].degree/2;
                a++ ) {
            pom = omvhei[a];
            omvhei[a] = omvhei[omv[currentVertex].firsthalfedge+omv[currentVertex].degree-zDrugiejStrony];
            omvhei[omv[currentVertex].firsthalfedge+omv[currentVertex].degree-zDrugiejStrony] = pom;
            zDrugiejStrony++;
          }
          omvheiIter++;
        }
      }
    }
    else if ( vertexCopiesCounter[i] == 0 ) {
      f = imv[i].firsthalfedge;
      d = imv[i].degree;
      currentVertex = isVertexInResultingNet[i];
      isSecondOrNext = false;
      for ( j = 0; j < d; j++ ) {
        if ( isHalfedgeInResultingNet[imvhei[f+j]] >= 0 ) {
          if ( !isSecondOrNext ) {
            omv[currentVertex].firsthalfedge = omvheiIter;
            isSecondOrNext = true;
          }
          omv[currentVertex].degree ++;
          heOut = imhe[imvhei[f+j]];
          V0[imvhei[f+j]] = currentVertex;
          omvhei[omvheiIter] = isHalfedgeInResultingNet[imvhei[f+j]];
          omvheiIter ++;
          face = heOut.facetnum;
          for ( heIter = 0; heIter < imfac[face].degree; heIter++ ) {
            if ( imhe[imfhei[imfac[face].firsthalfedge+heIter]].v1 == i ) {
              heIn = imhe[imfhei[imfac[face].firsthalfedge+heIter]];
              V1[imfhei[imfac[face].firsthalfedge+heIter]] = currentVertex;
            }
          }
        }
      }
    }
  }

  for ( i = 0; i < inv; i++ ) {
    if ( vertexCopiesCounter[i] == 0 )
      vertexCopiesCounter[i]++;
  }

    /* wpisuje wierzcholki, w odpowiednim momencie je podwajajac */

    /* przepisanie elementow nowej siatki */
    /*optc*/
  omvLastFreeSlot = *onv;
  for ( i = 0; i < inv; i++ ) {
    j = isVertexInResultingNet[i];
    isSecondOrNext = false;
    for ( k = 0; k < vertexCopiesCounter[i]; k++ ) {
      if ( !isSecondOrNext ) {
        memcpy ( optc+j*spdimen, iptc+i*spdimen, spdimen*sizeof(double) );
        isSecondOrNext = true;
      }
      else {
        omvLastFreeSlot--;
        memcpy ( optc+omvLastFreeSlot*spdimen, iptc+i*spdimen, spdimen*sizeof(double) );
      }
    }
  }

    /*omhe*/
  for ( i = 0; i < inhe; i++ ) {
    j = isHalfedgeInResultingNet[i];
    if ( j >= 0 ) {
      omhe[j].facetnum = isFacetInResultingNet[imhe[i].facetnum];
      k = imhe[i].otherhalf;
      omhe[j].otherhalf = k >= 0 ? isHalfedgeInResultingNet[k] : -1;
      omhe[j].v0 = V0[i];
      omhe[j].v1 = V1[i];
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

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*bsm_ExtractSubmeshVd*/

