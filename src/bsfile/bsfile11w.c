
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2015                            */
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

#include "bsfprivate.h"

boolean bsf_WriteBSMeshd ( int spdimen, int cpdimen, boolean rational, int degree,
                           int nv, const BSMvertex *mv, const int *mvhei,
                           const double *vc,
                           int nhe, const BSMhalfedge *mhe,
                           int nfac, const BSMfacet *mfac, const int *mfhei,
                           const char *name, int ident,
                           bsf_WriteAttr_fptr WriteAttr, void *userData )
{
  int sci;
  int i, j, d, l;

  sci = bsf_current_indentation;
  bsf_current_length += fprintf ( bsf_output, "%s\n", "%B-spline mesh" );
  bsf_current_length += fprintf ( bsf_output, "%s {",
          bsf_keyword[BSF_SYMB_BSMESH-BSF_FIRST_KEYWORD] );
  BSFeol
  bsf_current_indentation += 2;
  if ( name && *name ) {
    BSFwci
    bsf_current_length += fprintf ( bsf_output, "%s \"%s\"",
              bsf_keyword[BSF_SYMB_NAME-BSF_FIRST_KEYWORD], name );
    BSFeol
  }
  BSFwci
  bsf_WriteSpaceDim ( spdimen, rational );
  if ( degree >= -1 )
    bsf_WriteCurveDegree ( degree );
        /* write out the mesh vertices */
  BSFwci
  bsf_current_length += fprintf ( bsf_output, "%s %d {",
            bsf_keyword[BSF_SYMB_VERTICES-BSF_FIRST_KEYWORD],
            nv );
  BSFeol
  bsf_current_indentation += 2;
  BSFwci
  for ( i = 0; i < nv; i++ ) {
    bsf_current_length += fprintf ( bsf_output, "{" );
    d = mv[i].degree;
    l = mv[i].firsthalfedge;
    for ( j = 0; j < d; j++ ) {
      bsf_current_length += fprintf ( bsf_output, "%d", mvhei[l+j] );
      if ( j < d-1 )
        bsf_current_length += fprintf ( bsf_output, ", " );
    }
    bsf_current_length += fprintf ( bsf_output, " " );
    bsf_WritePointd ( cpdimen, &vc[i*spdimen] );
    bsf_current_length += fprintf ( bsf_output, "}" );
    if ( i < nv-1 ) {
      bsf_current_length += fprintf ( bsf_output, "," );
      BSFdol
      BSFwci
    }
  }
  bsf_current_length += fprintf ( bsf_output, "}" );
  BSFeol
  bsf_current_indentation -= 2;
        /* write out the mesh halfedges */
  BSFwci
  bsf_current_length += fprintf ( bsf_output, "%s %d {",
            bsf_keyword[BSF_SYMB_HALFEDGES-BSF_FIRST_KEYWORD],
            nhe );
  BSFeol
  bsf_current_indentation += 2;
  BSFwci
  for ( i = 0; i < nhe; i++ ) {
    bsf_current_length += fprintf ( bsf_output, "{%d, %d, %d, %d}",
              mhe[i].v0, mhe[i].v1, mhe[i].facetnum, mhe[i].otherhalf );
    if ( i < nhe-1 ) {
      bsf_current_length += fprintf ( bsf_output, "," );
      BSFdol
      BSFwci
    }
  }
  bsf_current_length += fprintf ( bsf_output, "}");
  BSFeol
  bsf_current_indentation -= 2;
  BSFwci
        /* write out the mesh facets */
  bsf_current_length += fprintf ( bsf_output, "%s %d {",
            bsf_keyword[BSF_SYMB_FACETS-BSF_FIRST_KEYWORD],
            nfac );
  BSFeol
  bsf_current_indentation += 2;
  BSFwci
  for ( i = 0; i < nfac; i++ ) {
    bsf_current_length += fprintf ( bsf_output, "{" );
    d = mfac[i].degree;
    l = mfac[i].firsthalfedge;
    for ( j = 0; j < d; j++ ) {
      bsf_current_length += fprintf ( bsf_output, "%d", mfhei[l+j] );
      if ( j < d-1 )
        bsf_current_length += fprintf ( bsf_output, ", " );
    }
    bsf_current_length += fprintf ( bsf_output, "}" );
    if ( i < nfac-1 ) {
      bsf_current_length += fprintf ( bsf_output, "," );
      BSFdol
      BSFwci
    }
  }
  bsf_current_length += fprintf ( bsf_output, "}" );
  BSFeol
  bsf_current_indentation -= 2;
  BSFwci
  if ( WriteAttr )
    WriteAttr ( userData );
  bsf_current_indentation = sci;
  BSFwci
  bsf_current_length += fprintf ( bsf_output, "}" );
  BSFeol
  return true;
} /*bsf_WriteBSMeshd*/

