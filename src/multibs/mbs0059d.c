
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2015                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

boolean mbs_TabCubicHFuncDer2d ( double a, double b, int nkn, const double *kn,
                                 double *hfunc, double *dhfunc, double *ddhfunc )
{
  double HFunc[16] = {1.0, 1.0, 0.0, 0.0,                  /* h00 */
                      0.0, 0.0, 1.0, 1.0,                  /* h10 */
                      0.0, (double)(1.0/3.0), 0.0, 0.0,    /* h01 */
                      0.0, 0.0, (double)(-1.0/3.0), 0.0};  /* h11 */
  int    i, j;
  double h, h_1, h_2, wsp[3];

  if ( a == 0.0 && b == 1.0 ) {
     for ( i = j = 0;  i < nkn;  i++, j += 4 )
       if ( !_mbs_multiBCHornerDer2d ( 3, 4, 1, 4, HFunc, kn[i],
                       &hfunc[j], &dhfunc[j], &ddhfunc[j], wsp ) )
         return false;
  }
  else {
    h = b - a;
    h_1 = 1.0/h;
    h_2 = h_1*h_1;
    for ( i = j = 0;  i < nkn;  i++, j += 4 ) {
      if ( !_mbs_multiBCHornerDer2d ( 3, 4, 1, 4, HFunc, (kn[i]-a)/h,
                      &hfunc[j], &dhfunc[j], &ddhfunc[j], wsp ) )
        return false;
      hfunc[j+2]   *= h;      hfunc[j+3]   *= h;
      dhfunc[j]    *= h_1;    dhfunc[j+1]  *= h_1;  
      ddhfunc[j]   *= h_2;    ddhfunc[j+1] *= h_2;  
      ddhfunc[j+2] *= h_1;    ddhfunc[j+3] *= h_1;  
    }
  }
  return true;
} /*mbs_TabCubicHFuncDer2d*/

