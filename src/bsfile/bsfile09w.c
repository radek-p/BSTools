
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009, 2014                            */
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

boolean bsf_WriteBSplineHoled ( int spdimen, int cpdimen, boolean rational,
                                int hole_k, const double *knots,
                                const point2d *domain_cp, const double *hole_cp,
                                const char *name, int ident,
                                bsf_WriteAttr_fptr WriteAttr, void *userData )
{
  int i;
  int sci;

  sci = bsf_current_indentation;
  BSFwci
  bsf_current_length += fprintf ( bsf_output,
          "%s", "%Hole in a piecewise bicubic B-spline surface" );
  BSFeol
  BSFwci
  bsf_current_length += fprintf ( bsf_output,
          "%s {", bsf_keyword[BSF_SYMB_BSHOLE-BSF_FIRST_KEYWORD] );
  BSFeol
  bsf_current_indentation += 2;
  if ( name && *name ) {
    BSFwci
    bsf_current_length += fprintf ( bsf_output, "%s \"%s\"",
              bsf_keyword[BSF_SYMB_NAME-BSF_FIRST_KEYWORD], name );
    BSFeol
  }
  bsf_WriteSpaceDim ( spdimen, rational );
  BSFwci
  bsf_current_length += fprintf ( bsf_output, "%s %d",
            bsf_keyword[BSF_SYMB_SIDES-BSF_FIRST_KEYWORD], hole_k );
  BSFeol
  for ( i = 0; i < hole_k; i++ ) {
    BSFwci
    bsf_current_length += fprintf ( bsf_output, "%s",
              bsf_keyword[BSF_SYMB_KNOTS-BSF_FIRST_KEYWORD] );
    bsf_current_indentation += 2;
    bsf_WriteKnotSequenced ( 10, &knots[11*i], false );
    bsf_current_indentation -= 2;
    BSFeol
  }
  BSFwci
  bsf_current_length += fprintf ( bsf_output, "%s %s ",
            bsf_keyword[BSF_SYMB_DOMAIN-BSF_FIRST_KEYWORD],
            bsf_keyword[BSF_SYMB_CPOINTS-BSF_FIRST_KEYWORD] );
  bsf_current_indentation += 2;
  bsf_WritePointsd ( 2, 1, 12*hole_k+1, 0, &domain_cp[0].x );
  bsf_current_indentation -= 2;
  BSFwci
  bsf_current_length += fprintf ( bsf_output, "%s ",
            bsf_keyword[BSF_SYMB_CPOINTS-BSF_FIRST_KEYWORD] );
  BSFeol
  bsf_current_indentation += 2;
  bsf_WritePointsd ( cpdimen, 1, 12*hole_k+1, 0, hole_cp );
  bsf_current_indentation -= 2;
  if ( WriteAttr )
    WriteAttr ( userData );
  bsf_current_indentation = sci;
  BSFwci
  bsf_current_length += fprintf ( bsf_output, "}" );
  BSFeol
  return true;
} /*bsf_WriteBSplineHoled*/

