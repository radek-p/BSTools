
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
#include "pknum.h"


/* /////////////////////////////////////////// */

void pkn_MultMatrixNumd ( int nrows, int rowlen,
                          int inpitch, const double *indata,
                          double a,
                          int outpitch, double *outdata )
{
  int i, j;

  for ( i = 0;  i < nrows;  i++, indata += inpitch, outdata += outpitch )
    for ( j = 0; j < rowlen; j++ )
      outdata[j] = a*indata[j];
} /*pkn_MultMatrixNumd*/

