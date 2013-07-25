
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"


void pkn_BandmFindQRMSizes ( int ncols, const bandm_profile *aprof,
                             int *qsize, int *rsize )
{
  int i, j, k, fr, lr, qs, rs;

  qs = ncols;
  rs = 0;
  for ( i = 0; i < ncols; i++ ) {
    k = fr = aprof[i].firstnz;
    lr = fr+aprof[i+1].ind-aprof[i].ind;
    qs += lr-i;
    for ( j = 0; j < i; j++ ) {
      if ( k < aprof[j].firstnz+aprof[j+1].ind-aprof[j].ind ) {
        k = min ( k, j );
        break;
      }
    }
    if ( k > i )
      rs += 1;
    else
      rs += i+1-k;
  }
  *qsize = qs;
  *rsize = rs;
} /*pkn_BandmFindQRMSizes*/

