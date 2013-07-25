
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"

/* /////////////////////////////////////////////////////////////////////////// */
void pkn_SPMFastMultMMd ( double *ac, double *bc,
                          int nnzab, int *abpos, index2 *aikbkj,
                          double *abc )
{
  int    i, j, j1;
  double s;

  j = abpos[0];
  for ( i = 0; i < nnzab; i++ ) {
    j1 = abpos[i+1];
    s = ac[aikbkj[j].i]*bc[aikbkj[j].j];
    for ( j++; j < j1; j++ )
      s += ac[aikbkj[j].i]*bc[aikbkj[j].j];
    abc[i] = s;
  }
} /*pkn_SPMFastMultMMd*/

