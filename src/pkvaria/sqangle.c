
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <memory.h>
#include <math.h>

#include "pkvaria.h"

double pkv_SqAngle ( double x, double y )
{
  if ( x > 0.0 ) {
    if ( y >= 0.0 )     return y/(x+y);
    else                return 3.0 + x/(x - y);
  }
  else if ( x < 0.0 ) {
    if ( y >= 0.0 )     return 1.0 + x/(x - y);
    else                return 2.0 + y/(x + y);
  }
  else {
    if ( y > 0.0 )      return 1.0;
    else if ( y < 0.0 ) return 3.0;
    else                return -1.0;  /* zero vector - angle undefined */
  }
} /*pkv_SqAngle*/

