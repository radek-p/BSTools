
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2014                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>

#include "pkvaria.h"
#include "pknum.h"


/* ///////////////////////////////////////////// */
/* solving nonlinear equations with one variable */
/* Illinois algorithm                            */

float pkn_Illinoisf ( float (*f)(void*,float), void *usrptr,
                      float a, float b, float eps, boolean *error )
{
  float c, cb, fa, fb, fc, fab, fbc, beta, d;
  char  sa, sb, sc, sab, sbc;

  *error = false;
  fa = f ( usrptr, a );  if ( !fa ) return a;
  fb = f ( usrptr, b );  if ( !fb ) return b;
  sa = pkv_signf ( fa );
  sb = pkv_signf ( fb );
  if ( sa*sb < 0 ) {
    fab = (fb-fa)/(b-a);
    sab = pkv_signf ( fab );
    beta = 1.0;
    d = fab;
    while ( fabs(b-a) > eps ) {
      c = b - fb/d;  cb = c-b;  if ( !cb ) return c;
      fc = f ( usrptr, c );  if ( !fc ) return c;
      sc = pkv_signf ( fc );
      fbc = (fc-fb)/cb;
      sbc = pkv_signf ( fbc );
      if ( sc*sb < 0 ) {
        a = b;  fa = fb;  sa = sb;
        beta = 1.0;  d = fbc;
      }
      else {
        if ( sab*sbc > 0 )
          beta = fbc/fab;
          else beta *= 0.5;
        d = (fc-beta*fa)/(c-a);
      }
      b = c;  fb = fc;  sb = sc;
      fab = fbc;  sab = sbc;
    }
    return b;
  }
  else {
    *error = true;
    return 0.0;
  }
} /*pkn_Illinoisf*/

