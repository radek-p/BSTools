
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

/* ////////////////////////////////////////////////////////////////////////// */
void _bsm_OutputWorkMeshd ( int spdimen,
                            int wnv, BSMvertex *wmv, int *wmvhei, double *wptc,
                            int wnhe, BSMhalfedge *wmhe,
                            int wnfac, BSMfacet *wmfac, int *wmfhei,
                            int *newvi, int *newhei, int *newfi,
                            int *onv, BSMvertex *omv, int *omvhei, double *optc,
                            int *onhe, BSMhalfedge *omhe,
                            int *onfac, BSMfacet *omfac, int *omfhei )
{
  int vd, vfhe, fd, ffhe, i, j, k;

        /* store the vertices in the final arrays */
  vfhe = 0;
  for ( i = 0; i < wnv; i++ ) {
    k = newvi[i];
    if ( k >= 0 ) {
      omv[k].degree = vd = wmv[i].degree;
      omv[k].firsthalfedge = vfhe;
      for ( j = 0; j < vd; j++ )
        omvhei[vfhe+j] = newhei[wmvhei[wmv[i].firsthalfedge+j]];
      if ( wptc )
        memcpy ( &optc[k*spdimen], &wptc[i*spdimen], spdimen*sizeof(double) );
      vfhe += vd;
    }
  }
        /* store the halfedges in the final arrays */
  for ( i = 0; i < wnhe; i++ ) {
    k = newhei[i];
    if ( k >= 0 ) {
      omhe[k].v0 = newvi[wmhe[i].v0];
      omhe[k].v1 = newvi[wmhe[i].v1];
      omhe[k].facetnum = newfi[wmhe[i].facetnum];
      if ( wmhe[i].otherhalf >= 0 )
        omhe[k].otherhalf = newhei[wmhe[i].otherhalf];
      else
        omhe[k].otherhalf = -1;
    }
  }
        /* store the facets in the final arrays */
  ffhe = 0;
  for ( i = 0; i < wnfac; i++ ) {
    k = newfi[i];
    if ( k >= 0 ) {
      omfac[k].degree = fd = wmfac[i].degree;
      omfac[k].firsthalfedge = ffhe;
      for ( j = 0; j < fd; j++ )
        omfhei[ffhe+j] = newhei[wmfhei[wmfac[i].firsthalfedge+j]];
      ffhe += fd;
    }
  }
} /*_bsm_OutputWorkMeshd*/

