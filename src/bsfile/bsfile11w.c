
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
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camerad.h"
#include "multibs.h"
#include "bsfile.h"

boolean bsf_WriteBSMeshd ( int spdimen, int cpdimen, boolean rational, int degree,
                           int nv, const BSMvertex *mv, const int *mvhei,
                           const double *vc,
                           int nhe, const BSMhalfedge *mhe,
                           int nfac, const BSMfacet *mfac, const int *mfhei,
                           const char *name,
                           bsf_WriteAttr_fptr WriteAttr, void *userData )
{
  int i, j, d, l;

  fprintf ( bsf_output, "%s\n", "%B-spline mesh" );
  fprintf ( bsf_output, "%s {\n", bsf_keyword[BSF_SYMB_BSMESH-BSF_FIRST_KEYWORD] );
  if ( name && *name )
    fprintf ( bsf_output, "  %s \"%s\"\n",
              bsf_keyword[BSF_SYMB_NAME-BSF_FIRST_KEYWORD], name );
  bsf_WriteSpaceDim ( spdimen, rational );
  if ( degree > 0 )
    bsf_WriteCurveDegree ( degree );
        /* write out the mesh vertices */
  fprintf ( bsf_output, "  %s %d {\n",
            bsf_keyword[BSF_SYMB_VERTICES-BSF_FIRST_KEYWORD],
            nv );
  for ( i = 0; i < nv; i++ ) {
    fprintf ( bsf_output, "    {" );
    d = mv[i].degree;
    l = mv[i].firsthalfedge;
    for ( j = 0; j < d; j++ ) {
      fprintf ( bsf_output, "%d", mvhei[l+j] );
      if ( j < d-1 )
        fprintf ( bsf_output, ", " );
    }
    fprintf ( bsf_output, " " );
    bsf_WritePointd ( cpdimen, &vc[i*spdimen] );
    fprintf ( bsf_output, "}" );
    if ( i < nv-1 )
      fprintf ( bsf_output, ",\n" );
  }
  fprintf ( bsf_output, "}\n" );
        /* write out the mesh halfedges */
  fprintf ( bsf_output, "  %s %d {\n",
            bsf_keyword[BSF_SYMB_HALFEDGES-BSF_FIRST_KEYWORD],
            nhe );
  for ( i = 0; i < nhe; i++ ) {
    fprintf ( bsf_output, "    {%d, %d, %d, %d}",
              mhe[i].v0, mhe[i].v1, mhe[i].facetnum, mhe[i].otherhalf );
    if ( i < nhe-1 )
      fprintf ( bsf_output, ",\n" );
  }
  fprintf ( bsf_output, "}\n");
        /* write out the mesh facets */
  fprintf ( bsf_output, "  %s %d {\n",
            bsf_keyword[BSF_SYMB_FACETS-BSF_FIRST_KEYWORD],
            nfac );
  for ( i = 0; i < nfac; i++ ) {
    fprintf ( bsf_output, "    {" );
    d = mfac[i].degree;
    l = mfac[i].firsthalfedge;
    for ( j = 0; j < d; j++ ) {
      fprintf ( bsf_output, "%d", mfhei[l+j] );
      if ( j < d-1 )
        fprintf ( bsf_output, ", " );
    }
    fprintf ( bsf_output, "}" );
    if ( i < nfac-1 )
      fprintf ( bsf_output, ",\n" );
  }
  fprintf ( bsf_output, "}\n" );
  if ( WriteAttr ) {
    bsf_current_indentation = 2;
    WriteAttr ( userData );
    bsf_current_indentation = 0;
  }
  fprintf ( bsf_output, "}\n\n" );
  return true;
} /*bsf_WriteBSMeshd*/

