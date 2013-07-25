
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>

#include "pkvaria.h"
#include "pkgeom.h"


/* /////////////////////////////////////////// */
/* procedures rearranging data in rectangular arrays */

void pkv_Rearrangec ( int nrows, int rowlen, int inpitch, int outpitch,
                      char *data )
{
  int i, j, k;

  if ( nrows <= 0 || rowlen <= 0 )
    return;
    
  if ( inpitch < outpitch ) {
    for ( i = nrows-1,  j = i*inpitch,  k = i*outpitch;
          i > 0;
          i--,  j -= inpitch, k -= outpitch )
      memmove ( &data[k], &data[j], rowlen );
  }
  else if ( inpitch > outpitch ) {
    for ( i = 1,  j = inpitch,  k = outpitch;
          i < nrows;
          i++,  j += inpitch,  k += outpitch )
      memmove ( &data[k], &data[j], rowlen );
  }
} /*pkv_Rearrangec*/

void pkv_Selectc ( int nrows, int rowlen, int inpitch, int outpitch,
                   const char *indata, char *outdata )
{
  int i, j, k;

  for ( i = j = k = 0; i < nrows; i++, j+= inpitch, k += outpitch )
    memcpy ( &outdata[k], &indata[j], rowlen );
} /*pkv_Selectc*/

void pkv_Movec ( int nrows, int rowlen, int pitch, int shift, char *data )
{
  int i;
  char *d;

  if ( shift > 0 ) {
    for ( i = nrows-1, d = &data[(nrows-1)*pitch]; i >= 0; i--, d -= pitch )
      memmove ( &d[shift], d, rowlen );
  }
  else if ( shift < 0 ) {
    for ( i = 0, d = data; i < nrows; i++, d += pitch )
      memmove ( &d[shift], d, rowlen );
  }
} /*pkv_Movec*/


