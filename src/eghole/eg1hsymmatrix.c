
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010                                  */
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

#undef CONST_
#define CONST_

#include "eg1holef.h"
#include "eg1holed.h"
#include "eg1hprivatef.h"
#include "eg1hprivated.h"
#include "eg1herror.h"

/* ////////////////////////////////////////////////////////////////////////// */
int g1h_SymPatchMatrixSize ( int hole_k )
{
  return (G1H_FINALDEG+1)*(G1H_FINALDEG+1)*(6*hole_k+1);
} /*g1h_SymPatchMatrixSize*/

