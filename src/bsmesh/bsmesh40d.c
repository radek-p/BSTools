
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
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "bsmesh.h"

#include "bsmprivate.h"

/* ////////////////////////////////////////////////////////////////////////// */
void bsm_TriangulateFacetNum ( int inv, BSMvertex *imv, int *imvhei,
                               int inhe, BSMhalfedge *imhe,
                               int facetNum,
                               int infac, BSMfacet *imfac, int *imfhei,
                               int *onv, int *onhe, int *onfac )
{
  int d;

  *onv = inv;
  *onhe = inhe;
  *onfac = infac;
  if ( imfac[facetNum].degree > 3 ) {
    d = imfac[facetNum].degree;
    *onhe = *onhe+(d-3)*2;
    *onfac = *onfac+(d-3);
  }
} /*bsm_TriangulateFacetNum*/

static double point3dist ( const point3d *p1, const point3d *p2 )
{
  vector3d v;

  SubtractPoints3d ( p1, p2, &v );
  return sqrt ( DotProduct3d ( &v, &v ) );
} /*point3dist*/

void bsm_TriangulateFacetd ( int spdimen, int *inv,
                             BSMvertex *imv, int *imvhei,
                             double *iptc, int *inhe,BSMhalfedge *imhe,
                             int *infac, BSMfacet *imfac, int *imfhei,
                             int facetNum )
{
  int     j, d, fhe, v0, v1, v2, facetNum2, minj;
  point3d a, b, c;
  double  miaraj, miaramin;

  if ( imfac[facetNum].degree > 3 ) {
    d = imfac[facetNum].degree;
    fhe = imfac[facetNum].firsthalfedge;
    miaramin = 1;
    minj = 0;
            /*szukamy takien polkrawedzi, za utworzony trojkat bedzie najblizej rownobocznego*/
    for ( j = 0; j < d; j++ ) {
      v0 = imhe[imfhei[fhe+j]].v0;
      v1 = imhe[imfhei[fhe+j]].v1;
      v2 = imhe[imfhei[fhe+(j+1)%d]].v1;
      SetPoint3d ( &a, iptc[3*v0], iptc[3*v0+1], iptc[3*v0+2] );
      SetPoint3d ( &b, iptc[3*v1], iptc[3*v1+1], iptc[3*v1+2] );
      SetPoint3d ( &c, iptc[3*v2], iptc[3*v2+1], iptc[3*v2+2] );
      miaraj = fabs ( point3dist ( &a, &c )/
                      (point3dist ( &a, &b ) + point3dist ( &b, &c )) - 0.5 );
      if ( miaraj < miaramin ) {
        miaramin = miaraj;
        minj = j;
      }
    }
    v0 = imhe[imfhei[minj+fhe]].v0;
    v1 = imhe[imfhei[fhe+(minj+1)%d]].v1;
    _bsm_DivideFacetd ( spdimen, inv, imv, imvhei, iptc, inhe, imhe,
                        infac, imfac, imfhei, v0, v1 );
    facetNum2 = *infac-1;
    bsm_TriangulateFacetd ( spdimen, inv, imv, imvhei, iptc, inhe, imhe,
                            infac, imfac, imfhei, facetNum2 );
    bsm_TriangulateFacetd ( spdimen, inv, imv, imvhei, iptc, inhe, imhe,
                            infac, imfac, imfhei, facetNum );
  }
} /*bsm_TriangulateFacetd*/

