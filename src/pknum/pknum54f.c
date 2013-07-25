
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2013                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"

void pkn_FindGivensRotationf ( float a, float b, float *c, float *s )
{
  float d, r, xi;

  if ( fabs(a) >= fabs(b) ) {
    if ( a ) {
      d = b/a;  r = sqrt ( 1.0+d*d );
      *s = xi = -d/r;  *c = sqrt ( 1.0-xi*xi );
    }
    else { *s = 0.0;  *c = 1.0; }
  }
  else {
    if ( b ) {
      d = a/b;  r = sqrt ( 1.0+d*d );
      *c = xi = d/r;  *s = -sqrt ( 1.0-xi*xi );
    }
    else { *s = 1.0;  *c = 0.0; }
  }
} /*pkn_FindGivensRotationf*/

void pkn_FindGivensRotXif ( float a, float b, float *c, float *s, float *xi )
{
  float d, r, _xi;

  if ( fabs(a) >= fabs(b) ) {
    if ( a ) {
      d = b/a;  r = sqrt ( 1.0+d*d );
      *s = _xi = -d/r;  *c = sqrt ( 1.0-_xi*_xi );
    }
    else { *s = _xi = 0.0;  *c = 1.0; }
  }
  else {
    if ( b ) {
      d = a/b;  r = sqrt ( 1.0+d*d );
      *c = _xi = d/r;  *s = -sqrt ( 1.0-_xi*_xi );
      _xi = r/d;
    }
    else { *s = _xi = 1.0;  *c = 0.0; }
  }
  *xi = _xi;
} /*pkn_FindGivensRotXif*/

void pkn_FindXiGivensRotf ( float xi, float *c, float *s )
{
  if ( fabs(xi) < 1.0 ) { *c = sqrt(1.0-xi*xi);  *s = -xi; }
  else if ( xi == 1.0 ) { *c = 1.0;  *s= 0.0; }
  else                  { *c = xi = 1.0/xi;  *s = -sqrt(1.0-xi*xi); }
} /*pkn_FindXiGivensRotf*/

void pkn_ApplyGivensRotationf ( float c, float s, float *a, float *b )
{
  float d;

  d  = c*(*a)-s*(*b);
  *b = s*(*a)+c*(*b);
  *a = d;
} /*pkn_ApplyGivensRotationf*/

void pkn_ApplySymGivensRotationf ( float c, float s,
                                   float *d, float *e, float *f )
{
  float u, v, w;

  u = c*s;
  v = c*c*(*d) - (u+u)*(*e) + s*s*(*f);
  w = u*(*d-*f) + (c*c-s*s)*(*e);
  *f = (u+u)*(*e) + s*s*(*d) + c*c*(*f);
  *d = v;
  *e = w;
} /*pkn_ApplySymGivensRotationf*/

