
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

void bsm_MergeMeshesd ( int spdimen,
                        int nv1, BSMvertex *mv1, int *mvhei1, double *vpc1,
                        int nhe1, BSMhalfedge *mhe1,
                        int nfac1, BSMfacet *mfac1, int *mfhei1,
                        int nv2, BSMvertex *mv2, int *mvhei2, double *vpc2,
                        int nhe2, BSMhalfedge *mhe2,
                        int nfac2, BSMfacet *mfac2, int *mfhei2,
                        int *onv, BSMvertex *omv, int *omvhei, double *ovpc,
                        int *onhe, BSMhalfedge *omhe,
                        int *onfac, BSMfacet *omfac, int *omfhei )
{
  int i, j, k, d, fh, fho;

        /* the output mesh may be stored in the same arrays as the first */
        /* of the meshes merged - copy the data only if it is not the case */
  if ( mv1 != omv ) {
    memcpy ( omv, mv1, nv1*sizeof(BSMvertex) );
    memcpy ( omvhei, mvhei1, nhe1*sizeof(int) );
    memcpy ( ovpc, vpc1, nv1*spdimen*sizeof(double) );
    memcpy ( omhe, mhe1, nhe1*sizeof(BSMhalfedge) );
    memcpy ( omfac, mfac1, nfac1*sizeof(BSMfacet) );
    memcpy ( omfhei, mfhei1, nhe1*sizeof(int) );
  }
  memcpy ( &ovpc[nv1*spdimen], vpc2, nv2*spdimen*sizeof(double) );
  for ( i = 0, j = nv1;  i < nv2;  i++, j++ ) {
    omv[j].degree = d = mv2[i].degree;
    fh = mv2[i].firsthalfedge;
    omv[j].firsthalfedge = fho = fh + nhe1;
    for ( k = 0; k < d; k++ )
      omvhei[fho+k] = mvhei2[fh+k] + nhe1;
  }
  for ( i = 0, j = nfac1;  i < nfac2;  i++, j++ ) {
    omfac[j].degree = d = mfac2[i].degree;
    fh = mfac2[i].firsthalfedge;
    omfac[j].firsthalfedge = fho = fh + nhe1;
    for ( k = 0; k < d; k++ )
      omfhei[fho+k] = mfhei2[fh+k] + nhe1;
  }
  for ( i = 0, j = nhe1;  i < nhe2;  i++, j++ ) {
    omhe[j].v0 = mhe2[i].v0 + nv1;
    omhe[j].v1 = mhe2[i].v1 + nv1;
    omhe[j].facetnum = mhe2[i].facetnum + nfac1;
    if ( mhe2[i].otherhalf >= 0 )
      omhe[j].otherhalf = mhe2[i].otherhalf + nhe1;
    else
      omhe[j].otherhalf = -1;
  }
  *onv = nv1 + nv2;
  *onhe = nhe1 + nhe2;
  *onfac = nfac1 + nfac2;
} /*bsm_MergeMeshesd*/

