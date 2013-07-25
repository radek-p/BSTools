
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
#include <stdio.h>

#include "pkvaria.h"


void WriteArrayf ( const char *name, int lgt, const float *tab )
{
  int i;

  printf ( "%s, n = %d\n", name, lgt );
  for ( i = 0; i < lgt; i++ )
    printf ( "%6.3f ", tab[i] );
  printf ( "\n" );
} /*WriteArrayf*/

void WriteArrayd ( const char *name, int lgt, const double *tab )
{
  int i;

  printf ( "%s", name );
  for ( i = 0; i < lgt; i++ )
    printf ( "%6.3f ", tab[i] );
  printf ( "\n" );
} /*WriteArrayd*/

