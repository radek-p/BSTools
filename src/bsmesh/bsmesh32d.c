
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
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

boolean bsm_SealMeshHoleNum ( int inv, const BSMvertex *imv, const int *imvhei,
                              int inhe, const BSMhalfedge *imhe,
                              int infac, const BSMfacet *imfac, const int *imfhei,
                              int nbhe,
                              int *onv, int *onhe, int *onfac )
{
  int cnt, v, vfhe, vd, e;

  if ( nbhe < 0 || nbhe >= inhe )
    goto failure;
  if ( imhe[nbhe].otherhalf >= 0 )
    goto failure;

  cnt = 0;
  e = nbhe;
  do {
    v = imhe[e].v1;
    vd = imv[v].degree;
    if ( vd < 2 )
      goto failure;
    vfhe = imv[v].firsthalfedge;
    e = imvhei[vfhe+vd-1];
    cnt ++;  /* count the halfedge */
  } while ( e != nbhe );

  *onv = inv;
  *onhe = inhe + cnt;
  *onfac = infac + 1;
  return true;

failure:
  *onv = *onhe = *onfac = -1;
  return false;
} /*bsm_SealMeshHoleNum*/

boolean bsm_SealMeshHoled ( int spdimen,
                            int inv, const BSMvertex *imv, const int *imvhei,
                            const double *iptc,
                            int inhe, const BSMhalfedge *imhe,
                            int infac, const BSMfacet *imfac, const int *imfhei,
                            int nbhe,
                            int *onv, BSMvertex *omv, int *omvhei, double *optc,
                            int *onhe, BSMhalfedge *omhe,
                            int *onfac, BSMfacet *omfac, int *omfhei )
{
  int cnt, v, ivfhe, ovfhe, vd, e, newe, ffhe;

  if ( nbhe < 0 || nbhe >= inhe )
    goto failure;
  if ( imhe[nbhe].otherhalf >= 0 )
    goto failure;

        /* copy the entire mesh, except for vertex halfedge numbers */
  memcpy ( omv, imv, inv*sizeof(BSMvertex) );
  memcpy ( optc, iptc, spdimen*inv*sizeof(double) );
  memcpy ( omhe, imhe, inhe*sizeof(BSMhalfedge) );
  memcpy ( omfac, imfac, infac*sizeof(BSMfacet) );
  memcpy ( omfhei, imfhei, inhe*sizeof(int) );
        /* add the new facet and seal the hole */
  omfac[infac].firsthalfedge = ffhe = inhe;
  cnt = 0;
  e = nbhe;
  do {
        /* add the other half of the halfedge */
    newe = ffhe+cnt;
    omhe[newe].v1 = imhe[e].v0;
    omhe[newe].v0 = v = imhe[e].v1;
    omhe[newe].facetnum = infac;
    omhe[newe].otherhalf = e;
    omhe[e].otherhalf = newe;
    vd = imv[v].degree;
    if ( vd < 2 )
      goto failure;
    omv[v].degree = vd + 1;
    ivfhe = imv[v].firsthalfedge;
    omvhei[v] = newe;  /* store it here for a while */
    e = imvhei[ivfhe+vd-1];
    cnt ++;  /* count the halfedge */
  } while ( e != nbhe );
  omfac[infac].degree = cnt;
  for ( e = 0; e < cnt; e++ )
    omfhei[ffhe+e] = ffhe+cnt-1-e;
        /* update the numbers of halfedges of the vertices */
  omv[0].firsthalfedge = 0;
  for ( v = 1; v < inv; v++ )
    omv[v].firsthalfedge = omv[v-1].firsthalfedge + omv[v-1].degree;
  for ( v = inv-1; v >= 0; v-- ) {
    vd = imv[v].degree;
    ivfhe = imv[v].firsthalfedge;
    ovfhe = omv[v].firsthalfedge;
    if ( vd < omv[v].degree )
      omvhei[ovfhe++] = omvhei[v];
    memcpy ( &omvhei[ovfhe], &imvhei[ivfhe], vd*sizeof(int) );
  }

  *onv = inv;
  *onhe = inhe + cnt;
  *onfac = infac + 1;
  return true;

failure:
  *onv = *onhe = *onfac = -1;
  return false;
} /*bsm_SealMeshHoled*/

