
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>

#include "pkvaria.h"

/* /////////////////////////////////////////// */
/* computing integer powers of a real numbers  */

double pkv_rpower ( double x, int e )
{
  double z;

  if ( e > 0 ) {
    z = 1.0;
    for ( ; ; ) {
      if ( e & 0x1 )
        z *= x;
      e = e >> 1;
      if ( !e )
        return z;
      x *= x;
    }
  }
  else if ( e < 0 )
    return 1.0 / pkv_rpower ( x, -e );
  else
    return 1.0;
} /*pkv_rpower*/

