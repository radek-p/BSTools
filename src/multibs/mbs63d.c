
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2015                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"

void mbs_DeallocatePolycurved ( int nelem, mbs_polycurved *poly )
{
  int i;

  if ( poly ) {
    for ( i = 0; i < nelem; i++ ) {
      if ( poly[i].knots )  PKV_FREE ( poly[i].knots );
      if ( poly[i].points ) PKV_FREE ( poly[i].points );
    }
    PKV_FREE ( poly );
  }
} /*mbs_DeallocatePolycurved*/

