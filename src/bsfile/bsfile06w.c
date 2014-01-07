
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

boolean bsf_WriteBSplineCurved ( int spdimen, int cpdimen, boolean rational,
                                 int deg, int lastknot, const double *knots,
                                 boolean closed, const double *cpoints,
                                 const char *name,
                                 bsf_WriteAttr_fptr WriteAttr, void *userData )
{
  fprintf ( bsf_output, "%s\n", "%B-spline curve" );
  fprintf ( bsf_output, "%s {\n", bsf_keyword[BSF_SYMB_BSCURVE-BSF_FIRST_KEYWORD] );
  if ( name && *name )
    fprintf ( bsf_output, "  %s \"%s\"\n",
              bsf_keyword[BSF_SYMB_NAME-BSF_FIRST_KEYWORD], name );
  bsf_WriteSpaceDim ( spdimen, rational );
  bsf_WriteCurveDegree ( deg );
  fprintf ( bsf_output, "  %s", bsf_keyword[BSF_SYMB_KNOTS-BSF_FIRST_KEYWORD] );
  bsf_WriteKnotSequenced ( lastknot, knots, closed );
  fprintf ( bsf_output, "  %s\n", bsf_keyword[BSF_SYMB_CPOINTS-BSF_FIRST_KEYWORD] );
  bsf_WritePointsd ( cpdimen, 1, lastknot-deg, 0, cpoints );
  if ( WriteAttr ) {
    bsf_current_indentation = 2;
    WriteAttr ( userData );
    bsf_current_indentation = 0;
  }
  fprintf ( bsf_output, "}\n\n" );
  return true;
} /*bsf_WriteBSplineCurved*/

