
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkgeom.h"
#include "multibs.h"
#include "convh.h"
#include "camerad.h"

#include "xgedit.h"
#include "spl3d.h"

void SelectDomPoints ( const Box2s *sbox, boolean mk )
{
  short a, b, xi, eta;
  int   i, j, k, pitch;

  pitch = lastknot_v-degree_v;
  i = j = 0;
  if ( sbox->x0 == sbox->x1 && sbox->y0 == sbox->y1 ) {
    a = xge_MINDIST;
    for ( k = 0; k < lastknot_u-degree_u; k++ ) {
      xi = xge_T2KnotWindMapKnotU ( &kwind,
                  GrevilleAbscissa ( degree_u, knots_u, k ) );
      b = (short)(abs ( xi - sbox->x0 ));
      if ( b < a ) { a = b;  i = k; }
    }
    if ( a < xge_MINDIST ) {
      a = xge_MINDIST;
      for ( k = 0; k < lastknot_v-degree_v; k++ ) {
        eta = xge_T2KnotWindMapKnotV ( &kwind,
                     GrevilleAbscissa ( degree_v, knots_v, k ) );
        b = (short)(abs ( eta - sbox->y0 ));
        if ( b < a ) { a = b;  j = k; }
      }
      if ( a < xge_MINDIST )
        MkPoint ( i, j, mk );
    }
  }
  else {
    for ( i = 0; i < lastknot_u-degree_u; i++ ) {
      xi = xge_T2KnotWindMapKnotU ( &kwind,
                  GrevilleAbscissa ( degree_u, knots_u, i ) );
      if ( xi >= sbox->x0 && xi <= sbox->x1 ) {
        for ( j = 0; j < pitch; j++ ) {
          eta = xge_T2KnotWindMapKnotV ( &kwind,
                       GrevilleAbscissa ( degree_v, knots_v, j ) );
          if ( eta >= sbox->y0 && eta <= sbox->y1 )
            MkPoint ( i, j, mk );
        }
      }
    }
  }
} /*SelectDomPoints*/

