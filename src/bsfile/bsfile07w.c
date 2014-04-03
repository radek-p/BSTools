
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

boolean bsf_WriteBezierPatchd ( int spdimen, int cpdimen, boolean rational,
                                int udeg, int vdeg, int pitch, const double *cpoints,
                                const char *name, int ident,
                                bsf_WriteAttr_fptr WriteAttr, void *userData )
{
  fprintf ( bsf_output, "%s\n", "%Bezier patch" );
  fprintf ( bsf_output, "%s {\n", bsf_keyword[BSF_SYMB_BPATCH-BSF_FIRST_KEYWORD] );
  if ( name && *name )
    fprintf ( bsf_output, "  %s \"%s\"\n",
              bsf_keyword[BSF_SYMB_NAME-BSF_FIRST_KEYWORD], name );
  bsf_WriteIdent ( ident );
  bsf_WriteSpaceDim ( spdimen, rational );
  bsf_WritePatchDegree ( udeg, vdeg );
  fprintf ( bsf_output, "  %s\n", bsf_keyword[BSF_SYMB_CPOINTS-BSF_FIRST_KEYWORD] );
  bsf_WritePointsd ( cpdimen, udeg+1, vdeg+1, pitch, cpoints );
  if ( WriteAttr ) {
    bsf_current_indentation = 2;
    WriteAttr ( userData );
    bsf_current_indentation = 0;
  }
  fprintf ( bsf_output, "}\n\n" );
  return true;
} /*bsf_WriteBezierPatchd*/

