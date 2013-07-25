
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2008                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "pkvaria.h"
#include "pknum.h"

#include "msgpool.h"

int pkn_NRBArraySize ( int n, const int *prof )
{
  int i, size;

  for ( i = size = 0;  i < n; i++ ) {
    if ( prof[i] < 0 || prof[i] > i )
      return 0;
    size += i-prof[i]+1;
  }
  return size;
} /*pkn_NRBArraySize*/

