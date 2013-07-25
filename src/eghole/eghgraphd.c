
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

#include "eg1holed.h"
#include "eg1hprivated.h"
#include "eg1herror.h"

/* ////////////////////////////////////////////////////////////////////////// */
void gh_DrawDomSurrndPatchesd ( GHoleDomaind *domain,
               void (*drawpatch) ( int n, int m, const point2d *cp ) )
{
  void    *sp;
  int     hole_k, i;
  point2d *cp, *cp0;

  sp = pkv_GetScratchMemTop ();
  cp0 = (point2d*)pkv_GetScratchMem ( 16*sizeof(point2d) );
  if ( !cp0 )
    goto wayout;

  hole_k = domain->hole_k;
  for ( i = 0; i < hole_k; i++ ) {
    if ( _gh_FindDomSurrndPatchd ( domain, i, 0, cp0 ) )
      drawpatch ( 3, 3, cp0 );
    cp = _gh_GetDomSurrndPatchd ( domain, i, 1 );
    if ( cp )
      drawpatch ( 3, 3, cp );
    cp = _gh_GetDomSurrndPatchd ( domain, i, 2 );
    if ( cp )
      drawpatch ( 3, 3, cp );
  }
wayout:
  pkv_SetScratchMemTop ( sp );
} /*gh_DrawDomSurrndPatchesd*/

