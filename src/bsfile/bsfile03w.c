
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

static boolean IsUniform ( int lastknot, const double *knots )
{
  int i;

  for ( i = 0; i <= lastknot; i++ )
    if ( knots[i] != (double)i )
      return false;
  return true;
} /*IsUniform*/

void bsf_WriteKnotSequenced ( int lastknot, const double *knots, boolean closed )
{
  int i;
  int sci;

  if ( lastknot < 1 )
    return;

  sci = bsf_current_indentation;
  bsf_current_indentation += 2;
  if ( closed ) {
    bsf_current_length +=
          fprintf ( bsf_output, " %s",
                    bsf_keyword[BSF_SYMB_CLOSED-BSF_FIRST_KEYWORD] );
  }

  if ( IsUniform ( lastknot, knots ) ) {
    bsf_current_length +=
          fprintf ( bsf_output, " %s %d",
                    bsf_keyword[BSF_SYMB_UNIFORM-BSF_FIRST_KEYWORD], lastknot );
  }
  else {
    bsf_current_length += fprintf ( bsf_output, " {" );
    for ( i = 0; i < lastknot; i++ ) {
      bsf_WriteDoubleNumber ( knots[i] );
      bsf_current_length += fprintf ( bsf_output, "," );
      BSFdol
      BSFwci
    }
    bsf_WriteDoubleNumber ( knots[lastknot] );
    bsf_current_length += fprintf ( bsf_output, "}" );
  }
  BSFeol
  bsf_current_indentation = sci;
} /*WriteKnotSequence*/

