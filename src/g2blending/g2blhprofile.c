
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

#include "g2blendingd.h"

/* ////////////////////////////////////////////////////////////////////////// */
int _g2bl_SetupHessian1Profile ( int lastknotu, int lastknotv, int *prof )
{
  int bldim, blnum;
  int i, j, i3;

  bldim = lastknotv-9;
  blnum = lastknotu-9;
  for ( i = i3 = 0; i < blnum; i++ )
    for ( j = 0; j < bldim; j++, i3 += 3 ) {
      prof[i3] = max ( 0, j-3 );
      if ( i > 3 )
        prof[i3] += (i-3)*bldim;
      prof[i3+2] = prof[i3+1] = prof[i3] *= 3;
    }
  return pkn_NRBArraySize ( i3, prof );
} /*_g2bl_SetupHessian1Profile*/

