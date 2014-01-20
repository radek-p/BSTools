
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>

#include "pkvaria.h"

boolean pkv_IncMultiCounter ( int dim, const int *limits, int *cnt )
{
  int i;

  for ( i = 0; i < dim; i++ ) {
    cnt[i] ++;
    if ( cnt[i] >= limits[i] )
      cnt[i] = 0;
    else
      return true;
  }
  return false;
} /*pkv_IncMultiCounter*/

