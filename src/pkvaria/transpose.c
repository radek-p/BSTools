
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

void pkv_TransposeMatrixc ( int nrows, int ncols, int elemsize,
                            int inpitch, const char *indata,
                            int outpitch, char *outdata )
{
  int i;

  for ( i = 0; i < nrows; i++ )
    pkv_Selectc ( ncols, elemsize, elemsize, outpitch,
                  indata+i*inpitch, outdata+i*elemsize );
} /*pkv_TransposeMatrix*/

