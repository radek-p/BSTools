
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

boolean mbs_TabQuinticHFuncDer3d ( double a, double b, int nkn, const double *kn,
                                   double *hfunc, double *dhfunc,
                                   double *ddhfunc, double *dddhfunc )
{
  double HFunc[36] = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0,   /* h00 */
                      0.0, 0.0, 0.0, 1.0, 1.0, 1.0,   /* h10 */
                      0.0, 0.2, 0.4, 0.0, 0.0, 0.0,   /* h01 */
                      0.0, 0.0, 0.0,-0.4,-0.2, 0.0,   /* h11 */
                      0.0, 0.0,0.05, 0.0, 0.0, 0.0,   /* h02 */
                      0.0, 0.0, 0.0,0.05, 0.0, 0.0};  /* h12 */
  int    i, j;
  double h, h2, h_1, h_2, h_3;

  if ( a == 0.0 && b == 1.0 ) {
    for ( i = j = 0;  i < nkn;  i++, j += 6 )
      if ( !mbs_multiBCHornerDer3d ( 5, 6, 1, 6, HFunc, kn[i],
              &hfunc[j], &dhfunc[j], &ddhfunc[j], &dddhfunc[j] ) )
        return false;
  }
  else {
    h = b - a;    
    h2 = h*h;     
    h_1 = 1.0/h;  
    h_2 = 1.0/h2;
    h_3 = h_1*h_2;
    for ( i = j = 0;  i < nkn;  i++, j += 6 ) {
      if ( !mbs_multiBCHornerDer3d ( 5, 6, 1, 6, HFunc, kn[i],
               &hfunc[j], &dhfunc[j], &ddhfunc[j], &dddhfunc[j] ) )
        return false;
      hfunc[j+2] *= h;       hfunc[j+3] *= h;
      hfunc[j+4] *= h2;      hfunc[j+5] *= h2;
      dhfunc[j] *= h_1;      dhfunc[j+1] *= h_1;
      dhfunc[j+4] *= h;      dhfunc[j+5] *= h;  
      ddhfunc[j] *= h_2;     ddhfunc[j+1] *= h_2;
      ddhfunc[j+2] *= h_1;   ddhfunc[j+3] *= h_1;
      dddhfunc[j] *= h_3;    dddhfunc[j+1] *= h_3;
      dddhfunc[j+2] *= h_2;  dddhfunc[j+3] *= h_2;
      dddhfunc[j+4] *= h_1;  dddhfunc[j+5] *= h_1;
    }
  }
  return true;
} /*mbs_TabQuinticHFuncDer3d*/

