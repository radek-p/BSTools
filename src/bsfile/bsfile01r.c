
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

boolean bsf_ReadIntNumber ( int *number )
{
  switch ( bsf_nextsymbol ) {
case BSF_SYMB_MINUS:
    bsf_GetNextSymbol ();
    if ( bsf_nextsymbol != BSF_SYMB_INTEGER )
      return false;
    else {
      *number = -bsf_nextint;
      bsf_GetNextSymbol ();
      return true;
    }

case BSF_SYMB_PLUS:
    bsf_GetNextSymbol ();
    if ( bsf_nextsymbol != BSF_SYMB_INTEGER )
      return false;
    else {
      *number = bsf_nextint;
      bsf_GetNextSymbol ();
      return true;
    }

case BSF_SYMB_INTEGER:
    *number = bsf_nextint;
    bsf_GetNextSymbol ();
    return true;

default:
    return false;
  }
} /*bsf_ReadIntNumber*/

boolean bsf_ReadDoubleNumber ( double *number )
{
  switch ( bsf_nextsymbol ) {
case BSF_SYMB_MINUS:
    bsf_GetNextSymbol ();
    if ( bsf_nextsymbol == BSF_SYMB_INTEGER || bsf_nextsymbol == BSF_SYMB_FLOAT ) {
      *number = -bsf_nextfloat;
      bsf_GetNextSymbol ();
      return true;
    }
    else
      return false;

case BSF_SYMB_PLUS:
    bsf_GetNextSymbol ();
    if ( bsf_nextsymbol != BSF_SYMB_INTEGER && bsf_nextsymbol != BSF_SYMB_FLOAT )
      return false;
case BSF_SYMB_INTEGER:
case BSF_SYMB_FLOAT:
    *number = bsf_nextfloat;
    bsf_GetNextSymbol ();
    return true;

default:
    return false;
  }
} /*bsf_ReadDoubleNumber*/

boolean bsf_ReadIdent ( int *ident )
{
  if ( bsf_nextsymbol != BSF_SYMB_IDENT )
    return false;
  bsf_GetNextSymbol ();
  if ( bsf_nextsymbol != BSF_SYMB_INTEGER )
    return false;
  if ( ident )
    *ident = bsf_nextint;
  bsf_GetNextSymbol ();
  return true;
} /*bsf_ReadIdent*/

