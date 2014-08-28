
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009, 2014                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camerad.h"
#include "multibs.h"
#include "bsfile.h"

#include "bsfprivate.h"

static void DeleteChar ( char *s )
{
  while ( *s ) {
    *s = *(s+1);
    s ++;
  }
} /*DeleteChar*/

static char *FindChar ( char c, char *s )
{
  while ( *s && *s != c )
    s++;
  return s;
} /*FindChar*/

void bsf_WriteDoubleNumber ( double x )
{
#define EPS 1.0e-15
#define UX  1.0e+12
#define LX  1.0e-3
  double ax;
  char   s[50], *d, *e;
  int    l;

  ax = fabs(x);
  if ( ax >= UX ) {      /* large number with exponent */
    sprintf ( s, "%15.12e", x );
    e = FindChar ( 'e', s ) + 1;
    if ( *e == '+' ) e++;
    goto trunc_fract;
  }
  else if ( ax < EPS )
    memcpy ( s, "0.0\000", 4 );
  else if ( ax < LX ) {  /* small number with exponent */
    sprintf ( s, "%15.12e", x );
    e = FindChar ( 'e', s ) + 1;
    e = FindChar ( '-', e ) + 1;
                           /* delete exponent leading zeros */
trunc_fract:
    while ( *e == '0' && *(e+1) )
      DeleteChar ( e );
                           /* delete trailing zeros of the fractional part */
    e = FindChar ( 'e', s )-1;
    while ( *e == '0' )
      DeleteChar ( e-- );
  }
  else {                 /* no exponent */
    sprintf ( s, "%14.12f", x );
    d = FindChar ( '.', s );
    if ( d-s < 14 )       /* truncate the fractional part */
      s[15] = 0;
                          /* delete trailing zeros */
    l = strlen ( s ) - 1;
    while ( s[l] == '0' && s[l-1] != '.' )
      s[l--] = 0;
  }
  bsf_current_length += fprintf ( bsf_output, "%s", s );
#undef LX
#undef UX
#undef EPS
} /*bsf_WriteDoubleNumber*/

void bsf_WriteAltDoubleNumber ( double x, int nzf )
{
  char s[50];
  int  l;

  if ( nzf > 47 ) nzf = 47; else if ( nzf < 1 ) nzf = 1;
  sprintf ( s, "%*.*f", nzf+2, nzf, x );
                          /* delete trailing zeros */
  l = strlen ( s ) - 1;
  while ( s[l] == '0' && s[l-1] != '.' )
    s[l--] = 0;
  bsf_current_length += fprintf ( bsf_output, "%s", s );
} /*bsf_WriteAltDoubleNumber*/

