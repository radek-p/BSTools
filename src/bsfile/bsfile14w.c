
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
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

void bsf_WriteDependencies ( int depname, int ndep, const int *dep )
{
  int i, k;

/* ndep must be at least 1 */
  bsf_WriteCurrentIndentation ();
  fprintf ( bsf_output, "%s { ", bsf_keyword[depname-BSF_FIRST_KEYWORD] );
  for ( i = k = 0;  i < ndep-1;  i++, k++ ) {
    fprintf ( bsf_output, "%d,", dep[i] );
    if ( k > 10 ) {
      fprintf ( bsf_output, "\n" );
      bsf_WriteCurrentIndentation ();
      fprintf ( bsf_output, "    " );
    }
    else
      fprintf ( bsf_output, " " );
  }
  fprintf ( bsf_output, "%d }\n", dep[ndep-1] );
} /*bsf_WriteDependencies*/

