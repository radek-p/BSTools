
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
void bsm_FacetEdgeDoublingNum ( int inv, BSMvertex *imv, int *imvhei,
                                int inhe, BSMhalfedge *imhe,
                                int infac, BSMfacet *imfac, int *imfhei,
                                int fn,
                                int *onv, int *onhe, int *onfac )
{
  int d;

  if ( fn < 0 || fn >= infac ) {
    *onv = *onhe = *onfac = 0;
    return;
  }
  d = imfac[fn].degree;
  *onv = inv + d;
  *onhe = inhe + 4*d;
  *onfac = infac + d;
} /*bsm_FacetEdgeDoublingNum*/

