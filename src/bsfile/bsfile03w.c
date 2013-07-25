
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009                                  */
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

  if ( lastknot < 1 )
    return;

  if ( closed ) {
    fprintf ( bsf_output, " %s",
              bsf_keyword[BSF_SYMB_CLOSED-BSF_FIRST_KEYWORD] );
  }

  if ( IsUniform ( lastknot, knots ) ) {
    fprintf ( bsf_output, " %s %d\n",
              bsf_keyword[BSF_SYMB_UNIFORM-BSF_FIRST_KEYWORD], lastknot );
  }
  else {
    fprintf ( bsf_output, "\n   {" );
    bsf_WriteDoubleNumber ( knots[0] );
    fprintf ( bsf_output, ",\n" );
    for ( i = 1; i < lastknot; i++ ) {
      fprintf ( bsf_output, "    " );
      bsf_WriteDoubleNumber ( knots[i] );
      fprintf ( bsf_output, ",\n" );
    }
    fprintf ( bsf_output, "    "  );
    bsf_WriteDoubleNumber ( knots[lastknot] );
    fprintf ( bsf_output, "}\n" );
  }
} /*WriteKnotSequence*/

