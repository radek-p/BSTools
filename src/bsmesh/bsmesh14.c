
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010                                  */
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
void bsm_RemoveVertexNum ( int inv, BSMvertex *imv, int *imvhei,
                           int inhe, BSMhalfedge *imhe,
                           int infac, BSMfacet *imfac, int *imfhei,
                           int nvr,
                           int *onv, int *onhe, int *onfac )
{
  int d, fhe;

  if ( nvr < 0 || nvr >= inv ) { /* invalid vertex number */
    *onv = *onhe = *onfac = -1;
    return;
  }
  d = imv[nvr].degree;
  fhe = imv[nvr].firsthalfedge;
  *onv = inv-1;
  *onfac = infac-d+1;
  if ( imhe[imvhei[fhe+d-1]].otherhalf < 0 ) { /* a boundary vertex to be removed */
    if ( imv[nvr].degree == 1 && imfac[imhe[fhe].facetnum].degree == 3 )
      *onv = *onhe = *onfac = -1;
    else
      *onhe = inhe-2*d+1;
  }
  else /* an inner vertex to be removed */
    *onhe = inhe-2*d;
} /*bsm_RemoveVertexNum*/

