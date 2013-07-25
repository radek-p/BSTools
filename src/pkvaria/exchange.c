
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

void pkv_Exchange ( void *x, void *y, int size )
{
  char *a, *b, *buf;
  int  bufsize, i;

  if ( size > 0 ) {
    bufsize = min ( size, 1024 );
    if ( (buf = pkv_GetScratchMem ( bufsize )) ) {
      a = (char*)x;
      b = (char*)y;
      while ( (i = min ( bufsize, size )) ) {
        memcpy ( buf, a, i );
        memcpy ( a, b,   i );
        memcpy ( b, buf, i );
        a += i;
        b += i;
        size -= i;
      }
      pkv_FreeScratchMem ( bufsize );
    }
  }
} /*pkv_Exchange*/

