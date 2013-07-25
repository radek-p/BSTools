
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"


/* /////////////////////////////////////////// */
void mbs_TransformAffKnotsf ( int degree, int lastknot, const float *inknots,
                              float a, float b, float *outknots )
{
  int   i;
  float un, lu, lv;

  un = inknots[degree];
  lu = inknots[lastknot-degree]-un;
  lv = b-a;
  for ( i = 0; i <= lastknot; i++ )
    outknots[i] = a + (inknots[i]-un)/lu*lv;
} /*mbs_TransformAffKnotsf*/

