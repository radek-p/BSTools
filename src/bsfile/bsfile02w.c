
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

void bsf_WritePointd ( int spdimen, const double *point )
{
  double defaultp[4] = {0.0,0.0,0.0,1.0};
  int i, k;

        /* do not write out default coordinates - y == 0, z == 0, w == 1 */
  k = spdimen;
  if ( k <= 4 ) {
    while ( k > 1 && point[k-1] == defaultp[k-1] )
      k--;
  }
  bsf_current_length += fprintf ( bsf_output, "{" );
  for ( i = 0; i < k-1; i++ ) {
    bsf_WriteDoubleNumber ( point[i] );
    bsf_current_length += fprintf ( bsf_output, "," );
  }
  bsf_WriteDoubleNumber ( point[k-1] );
  bsf_current_length += fprintf ( bsf_output, "}" );
} /*bsf_WritePointd*/

void bsf_WritePointsd ( int spdimen, int cols, int rows, int pitch,
                        const double *points )
{
  int i, j, k, l;
  int sci;

  sci = bsf_current_indentation;
  BSFwci
  bsf_current_length += fprintf ( bsf_output, "{" );
  BSFeol
/*  bsf_current_indentation += 2;*/
  BSFwci
  for ( i = k = 0;  i < cols;  i++, k += pitch ) {
    for ( j = 0, l = k;  j < rows;  j++, l += spdimen ) {
      bsf_WritePointd ( spdimen, &points[l] );
      if ( j < rows-1 ) {
        bsf_current_length += fprintf ( bsf_output, "," );
        BSFdol
        BSFwci
      }
    }
    if ( i < cols-1 ) {
      bsf_current_length += fprintf ( bsf_output, "," );
      BSFeol
      bsf_current_indentation -= 2;
      BSFwci
      bsf_current_length += fprintf ( bsf_output, "%c %d:", '%', i+1 );
      bsf_current_indentation += 2;
      BSFeol
      BSFwci
    }
    else {
      bsf_current_length += fprintf ( bsf_output, "}" );
      BSFeol
    }
  }
  bsf_current_indentation = sci;
} /*bsf_WritePointsd*/

void bsf_WritePointsMK ( int npoints, const unsigned int *mk )
{
  int i;
  int sci;

  if ( mk && npoints > 0 ) {
    sci = bsf_current_indentation;
    BSFwci
    bsf_current_indentation += 2;
    bsf_current_length +=
          fprintf ( bsf_output, "%s {",
                    bsf_keyword[BSF_SYMB_CPOINTSMK-BSF_FIRST_KEYWORD] );
    for ( i = 0; i < npoints-1; i++ ) {
      bsf_current_length += fprintf ( bsf_output, "%d,", mk[i] );
      BSFdol
      BSFwci
    }
    bsf_current_length += fprintf ( bsf_output, "%d}", mk[npoints-1] );
    BSFeol
    bsf_current_indentation = sci;
  }
} /*bsf_WritePointsMK*/

