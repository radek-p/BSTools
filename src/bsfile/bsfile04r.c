
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

boolean bsf_ReadSpaceDim ( int maxdim, int *spdimen )
{
  if ( bsf_nextsymbol != BSF_SYMB_INTEGER || bsf_nextint > maxdim )
    return false;   
  *spdimen = bsf_nextint;
  bsf_GetNextSymbol ();
  return true;
} /*bsf_ReadSpaceDim*/

boolean bsf_ReadCurveDegree ( int maxdeg, int *degree )
{
  if ( bsf_nextsymbol != BSF_SYMB_INTEGER || bsf_nextint > maxdeg )
    return false;   
  *degree = bsf_nextint;
  bsf_GetNextSymbol ();
  return true;
} /*bsf_ReadCurveDegree*/

boolean bsf_ReadPatchDegree ( int maxdeg, int *udeg, int *vdeg )
{
  if ( bsf_nextsymbol != BSF_SYMB_INTEGER || bsf_nextint > maxdeg )
    return false;  
  *udeg = bsf_nextint;
  bsf_GetNextSymbol ();
  if ( bsf_nextsymbol != BSF_SYMB_INTEGER || bsf_nextint > maxdeg )
    return false;   
  *vdeg = bsf_nextint;
  bsf_GetNextSymbol ();
  return true;
} /*bsf_ReadPatchDegree*/

