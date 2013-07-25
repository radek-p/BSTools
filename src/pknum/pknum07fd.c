
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"

/* ///////////////////////////////////////////////////////////////////////// */
int pkn_TMBSize ( int n )
{
  n = pkn_SymMatIndex ( n, 0 );
  return (n + 0x07) >> 0x03;
} /*pkn_TMBSize*/

boolean pkn_TMBElem ( byte *bittm, int i, int j )
{
  byte mask;

  i = pkn_SymMatIndex ( i, j );
  mask = 0x01 << (i & 0x07);
  i >>= 3;
  return (bittm[i] & mask) != 0;
} /*pkn_TMBElem*/

void pkn_TMBElemSet ( byte *bittm, int i, int j )
{
  byte mask;

  i = pkn_SymMatIndex ( i, j );
  mask = 0x01 << (i & 0x07);
  i >>= 3;
  bittm[i] |= mask;   /* set the bit */
} /*pkn_TMBElemSet*/

void pkn_TMBElemClear ( byte *bittm, int i, int j )
{
  byte mask;

  i = pkn_SymMatIndex ( i, j );
  mask = 0x01 << (i & 0x07);
  i >>= 3;
  bittm[i] &= ~mask;
} /*pkn_TMBElemClear*/

boolean pkn_TMBTestAndSet ( byte *bittm, int i, int j )
{
  byte mask;

  i = pkn_SymMatIndex ( i, j );
  mask = 0x01 << (i & 0x07);
  i >>= 3;
  if ( bittm[i] & mask )
    return true;        /* the bit was already set */
  else {
    bittm[i] |= mask;
    return false;       /* the bit has just been set */
  }
} /*pkn_TMBTestAndSet*/

boolean pkn_TMBTestAndClear ( byte *bittm, int i, int j )
{
  byte mask;

  i = pkn_SymMatIndex ( i, j );
  mask = 0x01 << (i & 0x07);
  i >>= 3;
  if ( bittm[i] & mask ) {
    bittm[i] &= ~mask;
    return true;        /* the bit has just been cleared */
  }
  else
    return false;       /* the bit was already zero */
} /*pkn_TMBTestAndClear*/

