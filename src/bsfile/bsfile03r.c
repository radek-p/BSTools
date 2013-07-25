
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

boolean bsf_ReadKnotSequenced ( int maxlastknot, int *lastknot, double *knots,
                                boolean *closed )
{
  int i;

  if ( closed )
    *closed = false;
  if ( bsf_nextsymbol == BSF_SYMB_CLOSED ) {
    if ( closed )
      *closed = true;
    bsf_GetNextSymbol ();
  }
  switch ( bsf_nextsymbol ) {
case BSF_SYMB_UNIFORM:
    bsf_GetNextSymbol ();
    if ( bsf_nextsymbol != BSF_SYMB_INTEGER || bsf_nextint < 1 )
      return false;
    if ( bsf_nextint > maxlastknot ) 
      return false;
    *lastknot = bsf_nextint;
    for ( i = 0; i <= bsf_nextint; i++ )
      knots[i] = (double)i;
    bsf_GetNextSymbol ();
    return true;

case BSF_SYMB_LBRACE:
    bsf_GetNextSymbol ();
    i = 0;
    for (;;) {
      if ( !bsf_ReadDoubleNumber ( &knots[i] ) )
        return false;
      i ++;
      if ( bsf_nextsymbol == BSF_SYMB_COMMA )
        bsf_GetNextSymbol ();
      else
        break;
    }
    if ( bsf_nextsymbol != BSF_SYMB_RBRACE )
      return false;
    bsf_GetNextSymbol ();
    i--;
    if ( i < 1 )
      return false;
    *lastknot = i;
        /* the knot sequence must be nondecreasing */
    while ( i > 0 )
      if ( knots[i] >= knots[i-1] )
        i--;
      else
        return false;
    return true;

default:
    return false;
  }
} /*bsf_ReadKnotSequenced*/

