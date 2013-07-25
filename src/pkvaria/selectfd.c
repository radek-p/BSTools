
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
/* procedures rearranging floating-point data in rectangular arrays; */
/* conversion between float and double */

void pkv_Selectfd ( int nrows, int rowlen, int inpitch, int outpitch,
                    const float *indata, double *outdata )
{
  int    i, j, k, l;

  for ( i = k = l = 0;
        i < nrows;
        i++, k += inpitch, l += outpitch ) {
    for ( j = 0; j < rowlen; j++ )
      outdata[l+j] = (double)indata[k+j];
  }
} /*pkv_Selectfd*/

void pkv_Selectdf ( int nrows, int rowlen, int inpitch, int outpitch,
                    const double *indata, float *outdata )
{
  int    i, j, k, l;

  for ( i = k = l = 0;
        i < nrows;
        i++, k += inpitch, l += outpitch ) {
    for ( j = 0; j < rowlen; j++ )
      outdata[l+j] = (float)indata[k+j];
  }
} /*pkv_Selectdf*/

