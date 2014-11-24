
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2014                            */
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
void bsm_ContractEdgeNum ( int inv, BSMvertex *imv, int *imvhei,
                           int inhe, BSMhalfedge *imhe,
                           int infac, BSMfacet *imfac, int *imfhei,
                           int nche,
                           int *onv, int *onhe, int *onfac )
{
  int v0, v1;

  if ( nche < 0 || nche >= inhe ) {
    *onv = *onhe = *onfac = -1;
    return;
  }

  *onv = inv-1;
  if ( imhe[nche].otherhalf == -1 ) {  /* a boundary edge to be contracted */
    if ( imfac[imhe[nche].facetnum].degree > 3 ) {
      *onhe = inhe-1;
      *onfac = infac;
    }
    else {
      *onhe = inhe-3;
      *onfac = infac-1;
    }
  }
  else {  /* an inner edge to be contracted */
    *onhe = inhe-2;
    *onfac = infac;
    if ( imfac[imhe[nche].facetnum].degree == 3 ) {
      *onhe -= 2;
      *onfac -= 1;
    }
    nche = imhe[nche].otherhalf;
    if ( imfac[imhe[nche].facetnum].degree == 3 ) {
      *onhe -= 2;
      *onfac -= 1;
    }
        /* if both vertices of the halfedge are boundary, */
        /* the contracted vertex should be split - in the future */
    v0 = imhe[nche].v0;
    v1 = imhe[nche].v1;
    if ( imhe[imvhei[imv[v0].firsthalfedge+imv[v0].degree-1]].otherhalf < 0 &&
         imhe[imvhei[imv[v1].firsthalfedge+imv[v1].degree-1]].otherhalf < 0 )
      *onv = *onhe = *onfac = -1;
  }
} /*bsm_ContractEdgeNum*/

