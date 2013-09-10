
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2013                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

#undef CONST_
#define CONST_

#include "eg1hprivate.h"

/* ///////////////////////////////////////////////////////////////////////// */
unsigned short _g1h_ExtendSupport ( int hole_k, unsigned short supp )
{
  unsigned short s, m;

  m = (unsigned short)(0x0001 << (hole_k-1));
  s = (unsigned short)(supp | ( supp << 1 ) | (supp >> 1));
  if ( supp & 0x0001 ) s = (unsigned short)(s | m);
  if ( supp & m ) s |= 0x0001;
  return s;
} /*_g1h_ExtendSupport*/

