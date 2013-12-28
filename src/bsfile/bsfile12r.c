
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2013                                  */
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
#include "bsmesh.h"
#include "bsfile.h"

#define CLAMP(a) \
  if ( a < 0.0 ) a = 0.0; else if ( a > 1.0 ) a = 1.0;

boolean bsf_ReadColour ( point3d *colour )
{
  int dim;

  if ( bsf_nextsymbol != BSF_SYMB_COLOR &&
       bsf_nextsymbol != BSF_SYMB_COLOUR )
    goto failure;
  bsf_GetNextSymbol ();
  if ( !bsf_ReadPointd ( 3, &colour->x, &dim ) )
    goto failure;
  CLAMP ( colour->x )
  CLAMP ( colour->y )
  CLAMP ( colour->z )
  return true;

failure:
  return false;
} /*bsf_ReadColour*/

