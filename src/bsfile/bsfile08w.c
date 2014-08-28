
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

boolean bsf_WriteBSplinePatchd ( int spdimen, int cpdimen, boolean rational,
                                 int udeg, int lastknotu, const double *knotsu,
                                 int vdeg, int lastknotv, const double *knotsv,
                                 boolean closed_u, boolean closed_v,
                                 int pitch, const double *cpoints,
                                 const char *name, int ident,
                                 bsf_WriteAttr_fptr WriteAttr, void *userData )
{
  int sci;

  sci = bsf_current_indentation;
  BSFwci
  bsf_current_length += fprintf ( bsf_output, "%s", "%B-spline patch" );
  BSFeol
  BSFwci
  bsf_current_length += fprintf ( bsf_output,
          "%s {", bsf_keyword[BSF_SYMB_BSPATCH-BSF_FIRST_KEYWORD] );
  BSFeol
  bsf_current_indentation += 2;
  if ( name && *name ) {
    BSFwci
    bsf_current_length += fprintf ( bsf_output, "%s \"%s\"",
              bsf_keyword[BSF_SYMB_NAME-BSF_FIRST_KEYWORD], name );
    BSFeol
  }
  bsf_WriteIdent ( ident );
  bsf_WriteSpaceDim ( spdimen, rational );
  bsf_WritePatchDegree ( udeg, vdeg );
  BSFwci
  bsf_current_length += fprintf ( bsf_output,
          "%s", bsf_keyword[BSF_SYMB_KNOTS_U-BSF_FIRST_KEYWORD] );
  bsf_WriteKnotSequenced ( lastknotu, knotsu, closed_u );
  BSFwci
  bsf_current_length += fprintf ( bsf_output,
          "%s", bsf_keyword[BSF_SYMB_KNOTS_V-BSF_FIRST_KEYWORD] );
  bsf_WriteKnotSequenced ( lastknotv, knotsv, closed_v );
  BSFwci
  bsf_current_length += fprintf ( bsf_output,
          "%s ", bsf_keyword[BSF_SYMB_CPOINTS-BSF_FIRST_KEYWORD] );
  bsf_current_indentation += 2;
  bsf_WritePointsd ( cpdimen, lastknotu-udeg, lastknotv-vdeg, pitch, cpoints );
  bsf_current_indentation -= 2;
  if ( WriteAttr )
    WriteAttr ( userData );
  bsf_current_length += fprintf ( bsf_output, "}" );
  BSFeol
  bsf_current_indentation = sci;
  return true;
} /*bsf_WriteBSplinePatchd*/

