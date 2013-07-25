
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009, 2010                            */
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

boolean bsf_WriteBezierCurved ( int spdimen, int cpdimen, boolean rational,
                                int deg, const double *cpoints,
                                const byte *mk, const char *name )
{
  fprintf ( bsf_output, "%s\n", "%Bezier curve" );
  fprintf ( bsf_output, "%s {\n", bsf_keyword[BSF_SYMB_BCURVE-BSF_FIRST_KEYWORD] );
  if ( name && *name )
    fprintf ( bsf_output, "  %s \"%s\"\n",
              bsf_keyword[BSF_SYMB_NAME-BSF_FIRST_KEYWORD], name );
  bsf_WriteSpaceDim ( spdimen, rational );
  bsf_WriteCurveDegree ( deg );
  fprintf ( bsf_output, "  %s\n", bsf_keyword[BSF_SYMB_CPOINTS-BSF_FIRST_KEYWORD] );
  bsf_WritePointsd ( cpdimen, 1, deg+1, 0, cpoints );
  if ( mk )
    bsf_WritePointsMK ( deg+1, mk );
  fprintf ( bsf_output, "}\n\n" );
  return true;
} /*bsf_WriteBezierCurved*/

