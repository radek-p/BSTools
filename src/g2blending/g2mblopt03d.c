
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "bsmesh.h"
#include "egholed.h"
#include "g2blendingd.h"
#include "g2mblendingd.h"

#include "g2blprivated.h"
#include "g2mblprivated.h"

/* ///////////////////////////////////////////////////////////////////////// */
void _g2mbl_AddCPIncrement ( int nvcp, int *vncpi, point3d *mvcp,
                             vector3d *incr, point3d *omvcp )
{
  int i, j;

  for ( i = 0; i < nvcp; i++ ) {
    j = vncpi[i];
    SubtractPoints3d ( &mvcp[j], &incr[i], &omvcp[j] );
  }
} /*_g2mbl_AddCPIncrement*/

