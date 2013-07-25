
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
void pkv_ReverseMatc ( int nrows, int rowlen, int pitch, char *data )
{
  int  i, h;
  char *p;

  h = nrows/2;
  for ( i = 0, p = data+(nrows-1)*pitch;
        i < h;
        i++, data += pitch, p -= pitch )
    pkv_Exchange ( data, p, rowlen );
} /*pkv_ReverseMatc*/

