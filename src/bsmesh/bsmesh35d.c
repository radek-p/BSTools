
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

boolean bsm_SplitBoundaryEdged ( int spdimen, int inv,
                                 const BSMvertex *imv, const int *imvhei,
                                 double *iptc, int inhe, const BSMhalfedge *imhe,
                                 int infac, const BSMfacet *imfac, const int *imfhei,
                                 int splithe,
                                 int *onv, BSMvertex *omv, int *omvhei,
                                 double *optc, int *onhe, BSMhalfedge *omhe,
                                 int *onfac, BSMfacet *omfac, int *omfhei )
{
  void   *sp;
  int    i;
  int    v0, v1, fac;

  sp = pkv_GetScratchMemTop ();

    /* krok pierwszy - sprawdzamy czy jest po polkrawedz brzegowa jesli nie - error */
  if ( imhe[splithe].otherhalf != -1 )
    goto failure;
  v0 = imhe[splithe].v0;
  v1 = imhe[splithe].v1;
  fac = imhe[splithe].facetnum;

    /* powstana:
      -nowy wierzcholek o numerze inv;
      -nowa polkrawedz o numerze inhe; */

    /* przeliczenie dlugosci tablic */
  *onhe = inhe+1;
  *onv = inv+1;
  *onfac = infac;

    /* przepisujemy polkrawedzie, az do konca */
  for ( i = 0; i < inhe; i++ )
    omhe[i] = imhe[i];

    /* zmieniamy zaznaczona polkrawedz */
  omhe[splithe].v1 = inv;

    /* dodajemy nowa polkrawedz */
  omhe[inhe].otherhalf = -1;
  omhe[inhe].facetnum = fac;
  omhe[inhe].v0 = inv;
  omhe[inhe].v1 = v1;

    /* przepisujemy wierzcholki, az do konca */
  for ( i = 0; i < inv; i++ )
    omv[i] = imv[i];

    /* dodajemy nowy wierzcholek */
  omv[inv].firsthalfedge = inhe;
  omv[inv].degree = 1;

    /* przepisujemy wspolrzedne */
  memcpy ( optc, iptc, inv*spdimen*sizeof(double) );

    /*dodajemy nowa wspolrzedna jako prosta srednia*/
  for ( i = 0; i < spdimen; i++ )
    optc[inv*spdimen+i] = (iptc[spdimen*v0+i] +
                           iptc[spdimen*v1+i])/2.0;

    /* omfac */
  for ( i = 0; i < infac; i++ )
    omfac[i] = imfac[i];

  omfac[fac].degree++;
  for ( i = fac+1; i < infac; i++ )
    omfac[i].firsthalfedge = imfac[i].firsthalfedge+1;

    /* omvhei */
  for ( i = 0; i < inhe; i++ )
    omvhei[i] = imvhei[i];

  omvhei[inhe] = inhe;

    /* omfhei */

    /* czesc przed nasza sciana */
  for ( i = 0; i < imfac[fac].firsthalfedge; i++ )
    omfhei[i] = imfhei[i];
    /* nasza sciana */
  do {
    omfhei[i] = imfhei[i];
  } while ( imfhei[i++] != splithe );

  omfhei[i++] = inhe;
    /* czesc po naszej scianie */
  for( ; i < inhe+1; i++ )
    omfhei[i] = imfhei[i-1];

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*bsm_SplitBoundaryEdged*/

