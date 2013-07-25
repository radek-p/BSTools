
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
/* computing degrees of normal vector Bezier patches      */

void mbs_BezP3NormalDeg ( int degreeu, int degreev, int *ndegu, int *ndegv )
{
  *ndegu = 2*degreeu-1;
  *ndegv = 2*degreev-1;
} /*mbs_BezP3NormalDeg*/

void mbs_BezP3RNormalDeg ( int degreeu, int degreev, int *ndegu, int *ndegv )
{
  *ndegu = 3*degreeu-2;
  *ndegv = 3*degreev-2;
} /*mbs_BezP3RNormalDeg*/

