
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

/* diagnostic functions */

boolean _g1h_VerifyJunctionFunctionsd ( GHoleDomaind *domain )
{
  G1HolePrivateRecd *privateG1;
  double *b01, *c01, *f01, *g01, *b11, *c11, *f11, *g11;
  int    hole_k, i, j;

  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  G1GetPolyAddr ( privateG1->jfunc, b01, c01, f01, g01, b11, c11, f11, g11 );
  for ( i = 0; i < hole_k; i++ ) {
    for ( j = 1; j <= G1_CG01DEG; j++ ) {
      if ( c01[i*(G1_CG01DEG+1)+j]*c01[i*(G1_CG01DEG+1)] <= 0.0 )
        return false;
      if ( g01[i*(G1_CG01DEG+1)+j]*g01[i*(G1_CG01DEG+1)] <= 0.0 )
        return false;
      if ( c11[i*(G1_CG11DEG+1)+j]*c11[i*(G1_CG11DEG+1)] <= 0.0 )
        return false;
      if ( g11[i*(G1_CG11DEG+1)+j]*g11[i*(G1_CG11DEG+1)] <= 0.0 )
        return false;
    }
  }
  return true;
} /*_g1h_VerifyJunctionFunctionsd*/

boolean _g1h_VerifyDomPatchesd ( GHoleDomaind *domain )
{
  /* ****************** */
  return true;
} /*_g1h_VerifyDomPatchesd*/

