
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */ 
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009, 1020                            */
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

boolean bsf_ReadPointd ( int maxcpdimen, double *point, int *cpdimen )
{
  int i;

  if ( bsf_nextsymbol != BSF_SYMB_LBRACE )
    return false;
  memset ( point, 0, maxcpdimen*sizeof(double) );
  if ( maxcpdimen == 4)
    point[3] = 1.0;
  bsf_GetNextSymbol ();
  for ( i = 0; i < maxcpdimen; i++ ) {
    if ( !bsf_ReadDoubleNumber ( &point[i] ) )
      return false;
    if ( bsf_nextsymbol != BSF_SYMB_COMMA )
      break;
    bsf_GetNextSymbol ();
  }
  if ( bsf_nextsymbol != BSF_SYMB_RBRACE )
    return false;
  *cpdimen = min ( 4, i+1 );
  bsf_GetNextSymbol ();
  return true;
} /*bsf_ReadPointd*/

int bsf_ReadPointsd ( int maxcpdimen, int maxnpoints,
                      double *points, int *cpdimen )
{
  int    ncp, spd0, spd1;
  double *p;

  ncp = spd0 = 0;
  if ( bsf_nextsymbol != BSF_SYMB_LBRACE )
    return 0;
  bsf_GetNextSymbol ();
  for ( p = points; ; p += maxcpdimen ) {
    if ( ncp >= maxnpoints )
      return 0;
    if ( bsf_ReadPointd ( maxcpdimen, p, &spd1 ) ) {
      ncp ++;
      spd0 = max ( spd0, spd1 );
    }
    else
      return 0;
    if ( bsf_nextsymbol != BSF_SYMB_COMMA )
      break;
    bsf_GetNextSymbol ();
  }
  if ( bsf_nextsymbol != BSF_SYMB_RBRACE )
    return 0;
  bsf_GetNextSymbol ();
  *cpdimen = spd0;
  return ncp;
} /*bsf_ReadPointsd*/

int bsf_ReadPointsMK ( int maxnpoints, byte *mk )
{
  int np;
  int num;

  if ( bsf_nextsymbol != BSF_SYMB_LBRACE )
    return 0;
  bsf_GetNextSymbol ();
  for ( np = 0; np < maxnpoints; ) {
    if ( !bsf_ReadIntNumber ( &num ) )
      return 0;
    if ( mk )  /* NULL is allowed here */
      mk[np++] = (byte)num;
    else
      np ++;
    if ( bsf_nextsymbol != BSF_SYMB_COMMA )
      break;
    bsf_GetNextSymbol ();
  }
  if ( bsf_nextsymbol != BSF_SYMB_RBRACE )
    return 0;
  bsf_GetNextSymbol ();
  return np;
} /*bsf_ReadPointsMK*/

