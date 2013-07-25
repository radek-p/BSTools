
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>

#include "pkvaria.h"


/* /////////////////////////////////////////// */
/* initialization of a matrix to zero */

void pkv_ZeroMatc ( int nrows, int rowlen, int pitch, char *data )
{
  int i;

  for ( i = 0; i < nrows; i++, data += pitch )
    memset ( data, 0, rowlen );
} /*pkv_Selectdf*/

