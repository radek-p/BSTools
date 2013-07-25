
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

/* /////////////////////////////////////////// */
void pkv_Sort2f ( float *a, float *b )
{
  float c;

  if ( *a > *b ) { c = *a;  *a = *b;  *b = c; }
} /*pkv_Sort2f*/
