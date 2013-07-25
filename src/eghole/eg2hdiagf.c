
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
#include "eg2holef.h"

#include "eg2hprivatef.h"
#include "eg2herror.h"

/* diagnostic functions */

boolean _g2h_VerifyJunctionFunctionsf ( GHoleDomainf *domain )
{
  G2HolePrivateRecf *privateG2;
  float *b01, *c01, *f01, *g01, *b11, *c11, *f11, *g11;
  int   hole_k, i, j;

  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  G2GetPolyAddr ( privateG2->jfunc, b01, c01, f01, g01, b11, c11, f11, g11 );
  for ( i = 0; i < hole_k; i++ ) {
    for ( j = 1; j <= G2_CG01DEG; j++ ) {
      if ( c01[i*(G2_CG01DEG+1)+j]*c01[i*(G2_CG01DEG+1)] <= 0.0 )
        return false;
      if ( g01[i*(G2_CG01DEG+1)+j]*g01[i*(G2_CG01DEG+1)] <= 0.0 )
        return false;
      if ( c11[i*(G2_CG11DEG+1)+j]*c11[i*(G2_CG11DEG+1)] <= 0.0 )
        return false;
      if ( g11[i*(G2_CG11DEG+1)+j]*g11[i*(G2_CG11DEG+1)] <= 0.0 )
        return false;
    }
  }
  return true;
} /*_g2h_VerifyJunctionFunctionsf*/

boolean _g2h_VerifyDomPatchesf ( GHoleDomainf *domain )
{
  /* ****************** */
  return true;
} /*_g2h_VerifyDomPatchesf*/

