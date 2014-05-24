
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

boolean bsm_GlueTwoHalfedgesNum ( int inv, const BSMvertex *imv, const int *imvhei,
                                  int inhe, const BSMhalfedge *imhe,
                                  int infac, const BSMfacet *imfac, const int *imfhei,
                                  int he1, int he2,
                                  int *onv, int *onhe, int *onfac )
{
  int f1, f2, e, e1, e2, fd, ffhe;
  int v10, v11, v20, v21, vd, vfhe;

  if ( he1 < 0 || he1 >= inhe || he2 < 0 || he2 > inhe )
    goto failure;
        /* the halfedges to be glued must be boundary */
  if ( imhe[he1].otherhalf >= 0 || imhe[he2].otherhalf >= 0 )
    goto failure;
        /* the halfedges to be glued must belong to different facets */
  f1 = imhe[he1].facetnum;
  f2 = imhe[he2].facetnum;
  if ( f1 == f2 )
    goto failure;
        /* two different facets may have at most one common edge */
  fd = imfac[f1].degree;
  ffhe = imfac[f1].firsthalfedge;
  for ( e = 0; e < fd; e++ ) {
    e1 = imfhei[ffhe+e];
    if ( imhe[e1].otherhalf >= 0 ) {
      if ( imhe[imhe[e1].otherhalf].facetnum == f2 )
        goto failure;
    }
  }
  fd = imfac[f2].degree;
  ffhe = imfac[f2].firsthalfedge;
  for ( e = 0; e < fd; e++ ) {
    e1 = imfhei[ffhe+e];
    if ( imhe[e1].otherhalf >= 0 ) {
      if ( imhe[imhe[e1].otherhalf].facetnum == f1 )
        goto failure;
    }
  }
        /* glueing the halfedges must not produce a loop */
  v10 = imhe[he1].v0;
  v11 = imhe[he1].v1;
  v20 = imhe[he2].v0;
  v21 = imhe[he2].v1;
  vd = imv[v11].degree;
  vfhe = imv[v11].firsthalfedge;
  e1 = imvhei[vfhe+vd-1];
  if ( imhe[e1].v1 == v20 )
    goto failure;
  vd = imv[v21].degree;
  vfhe = imv[v21].firsthalfedge;
  e2 = imvhei[vfhe+vd-1];
  if ( imhe[e2].v1 == v10 )
    goto failure;
        /* now determine the case, i.e. how many vertices will be removed */
  if ( v21 != v10 )
    inv --;
  if ( v11 != v20 )
    inv --;
  *onv = inv;  *onhe = inhe;  *onfac = infac;
  return true;

failure:
  *onv = *onhe = *onfac = -1;
  return false;
} /*bsm_GlueTwoHalfedgesNum*/

boolean bsm_GlueTwoHalfedgesd ( int spdimen,
                                int inv, const BSMvertex *imv, const int *imvhei,
                                const double *ivc,
                                int inhe, const BSMhalfedge *imhe,
                                int infac, const BSMfacet *imfac, const int *imfhei,
                                int he1, int he2,
                                int *onv, BSMvertex *omv, int *omvhei, double *ovc,
                                int *onhe, BSMhalfedge *omhe,
                                int *onfac, BSMfacet *omfac, int *omfhei )
{
  int f1, f2, e, e1, e2, fd, ffhe;
  int v10, v11, v20, v21, vd1, vfhe1, vd2, vfhe2, vd, vfhe;
  int iv, ov, j;

  if ( he1 < 0 || he1 >= inhe || he2 < 0 || he2 > inhe )
    goto failure;
        /* the halfedges to be glued must be boundary */
  if ( imhe[he1].otherhalf >= 0 || imhe[he2].otherhalf >= 0 )
    goto failure;
        /* the halfedges to be glued must belong to different facets */
  f1 = imhe[he1].facetnum;
  f2 = imhe[he2].facetnum;
  if ( f1 == f2 )
    goto failure;
        /* two different facets may have at most one common edge */
  fd = imfac[f1].degree;
  ffhe = imfac[f1].firsthalfedge;
  for ( e = 0; e < fd; e++ ) {
    e1 = imfhei[ffhe+e];
    if ( imhe[e1].otherhalf >= 0 ) {
      if ( imhe[imhe[e1].otherhalf].facetnum == f2 )
        goto failure;
    }
  }
  fd = imfac[f2].degree;
  ffhe = imfac[f2].firsthalfedge;
  for ( e = 0; e < fd; e++ ) {
    e1 = imfhei[ffhe+e];
    if ( imhe[e1].otherhalf >= 0 ) {
      if ( imhe[imhe[e1].otherhalf].facetnum == f1 )
        goto failure;
    }
  }
        /* glueing the halfedges must not produce a loop */
  v10 = imhe[he1].v0;
  v11 = imhe[he1].v1;
  v20 = imhe[he2].v0;
  v21 = imhe[he2].v1;
  vd1 = imv[v11].degree;
  vfhe1 = imv[v11].firsthalfedge;
  e1 = imvhei[vfhe1+vd1-1];
  if ( imhe[e1].v1 == v20 )
    goto failure;
  vd2 = imv[v21].degree;
  vfhe2 = imv[v21].firsthalfedge;
  e2 = imvhei[vfhe2+vd2-1];
  if ( imhe[e2].v1 == v10 )
    goto failure;
        /* copy the halfedges and facets */
  memcpy ( omhe, imhe, inhe*sizeof(BSMhalfedge) );
  memcpy ( omfac, imfac, infac*sizeof(BSMfacet) );
  memcpy ( omfhei, imfhei, inhe*sizeof(int) );
        /* actually glue the halfedges */
  omhe[he1].otherhalf = he2;
  omhe[he2].otherhalf = he1;
        /* now determine the case, i.e. how many vertices will be removed */
  if ( e1 == he2 ) {
    if ( e2 == he1 ) {  /* remove a boundary cycle of length 2, */
                        /* all vertices remain */
      memcpy ( omv, imv, inv*sizeof(BSMvertex) );
      memcpy ( omvhei, imvhei, inhe*sizeof(int) );
      memcpy ( ovc, ivc, inv*spdimen*sizeof(double) );
      *onv = inv;
    }
    else {
        /* copy the vertices, skipping v21 */
      vfhe = 0;
      for ( iv = ov = 0;  iv < inv;  iv++ )
        if ( iv != v21 ) {
          vd = imv[iv].degree;
          memcpy ( &omvhei[vfhe], &imvhei[imv[iv].firsthalfedge], vd*sizeof(int) );
          if ( iv == v10 ) {    /* merge v21 to v10 */
            memcpy ( &omvhei[vfhe+vd], &imvhei[imv[v21].firsthalfedge], vd2*sizeof(int) );
            vd += vd2;
            for ( j = 0; j < spdimen; j++ )
              ovc[spdimen*ov+j] = 0.5*(ivc[spdimen*iv+j]+ivc[spdimen*v21+j]);
          }
          else  /* just copy the vertex position */
            memcpy ( &ovc[ov*spdimen], &ivc[iv*spdimen], spdimen*sizeof(double) );
          omv[ov].firsthalfedge = vfhe;
          omv[ov].degree = vd;
          vfhe += vd;
          ov ++;
        }
        /* renumber the vertices in the output halfedges */
      for ( e = 0; e < inhe; e++ ) {
        if ( omhe[e].v0 == v21 ) omhe[e].v0 = v10;
        if ( omhe[e].v0 >= v21 ) omhe[e].v0 --;
        if ( omhe[e].v1 == v21 ) omhe[e].v1 = v10;
        if ( omhe[e].v1 >= v21 ) omhe[e].v1 --;
      }
      *onv = inv-1;
    }
  }
  else if ( e2 == he1 ) {
        /* copy the vertices, skipping v11 */
    vfhe = 0;
    for ( iv = ov = 0;  iv < inv;  iv++ )
      if ( iv != v11 ) {
        vd = imv[iv].degree;
        memcpy ( &omvhei[vfhe], &imvhei[imv[iv].firsthalfedge], vd*sizeof(int) );
        if ( iv == v20 ) {    /* merge v11 to v20 */
          memcpy ( &omvhei[vfhe+vd], &imvhei[imv[v11].firsthalfedge], vd1*sizeof(int) );
          vd += vd1;
          for ( j = 0; j < spdimen; j++ )
            ovc[spdimen*ov+j] = 0.5*(ivc[spdimen*iv+j]+ivc[spdimen*v11+j]);
        }
        else  /* just copy the vertex position */
          memcpy ( &ovc[ov*spdimen], &ivc[iv*spdimen], spdimen*sizeof(double) );
        omv[ov].firsthalfedge = vfhe;
        omv[ov].degree = vd;
        vfhe += vd;
        ov ++;
      }
        /* renumber the vertices in the output halfedges */
    for ( e = 0; e < inhe; e++ ) {
      if ( omhe[e].v0 == v11 ) omhe[e].v0 = v20;
      if ( omhe[e].v0 >= v11 ) omhe[e].v0 --;
      if ( omhe[e].v1 == v11 ) omhe[e].v1 = v20;
      if ( omhe[e].v1 >= v11 ) omhe[e].v1 --;
    }
    *onv = inv-1;
  }
  else {
        /* copy the vertices, skipping v11 and v21 */
    vfhe = 0;
    for ( iv = ov = 0;  iv < inv;  iv++ )
      if ( iv != v11 && iv != v21 ) {
        vd = imv[iv].degree;
        memcpy ( &omvhei[vfhe], &imvhei[imv[iv].firsthalfedge], vd*sizeof(int) );
        if ( iv == v10 ) {      /* merge v21 with v10 */
          memcpy ( &omvhei[vfhe+vd], &imvhei[imv[v21].firsthalfedge], vd2*sizeof(int) );
          vd += vd2;
          for ( j = 0; j < spdimen; j++ )
            ovc[spdimen*ov+j] = 0.5*(ivc[spdimen*iv+j]+ivc[spdimen*v21+j]);
        }
        else if ( iv == v20 ) { /* merge v11 with v20 */
          memcpy ( &omvhei[vfhe+vd], &imvhei[imv[v11].firsthalfedge], vd1*sizeof(int) );
          vd += vd1;
          for ( j = 0; j < spdimen; j++ )
            ovc[spdimen*ov+j] = 0.5*(ivc[spdimen*iv+j]+ivc[spdimen*v11+j]);
        }
        else  /* just copy the vertex position */
          memcpy ( &ovc[ov*spdimen], &ivc[iv*spdimen], spdimen*sizeof(double) );
        omv[ov].firsthalfedge = vfhe;
        omv[ov].degree = vd;
        vfhe += vd;
        ov ++;
      }
        /* renumber the vertices in the output halfedges */
    for ( e = 0; e < inhe; e++ ) {
      if ( omhe[e].v0 == v11 ) omhe[e].v0 = v20;
      else if ( omhe[e].v0 == v21 ) omhe[e].v0 = v10;
      if ( omhe[e].v0 >= v11 ) omhe[e].v0 --;
      if ( omhe[e].v0 >= v21 ) omhe[e].v0 --;
      if ( omhe[e].v1 == v11 ) omhe[e].v1 = v20;
      else if ( omhe[e].v1 == v21 ) omhe[e].v1 = v10;
      if ( omhe[e].v1 >= v11 ) omhe[e].v1 --;
      if ( omhe[e].v1 >= v21 ) omhe[e].v1 --;
    }
    *onv = inv-2;
  }
  *onhe = inhe;  *onfac = infac;
  return true;

failure:
  *onv = *onhe = *onfac = -1;
  return false;
} /*bsm_GlueTwoHalfedgesd*/

