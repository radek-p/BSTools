
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

#include "bsfprivate.h"

boolean bsf_WriteColour ( point3d *colour )
{
  BSFwci
        /* the keyword we write is in the British spelling */
  bsf_current_length += fprintf ( bsf_output,
          "%s { ", bsf_keyword[BSF_SYMB_COLOUR-BSF_FIRST_KEYWORD] );
  bsf_WriteAltDoubleNumber ( colour->x, 3 );
  bsf_current_length += fprintf ( bsf_output, ", " );
  bsf_WriteAltDoubleNumber ( colour->y, 3 );
  bsf_current_length += fprintf ( bsf_output, ", " );
  bsf_WriteAltDoubleNumber ( colour->z, 3 );
  bsf_current_length += fprintf ( bsf_output, " }" );
  BSFeol
  return true;
} /*bsf_WriteColour*/

