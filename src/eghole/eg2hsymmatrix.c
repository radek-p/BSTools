
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

#include "eg2holef.h"
#include "eg2holed.h"
#include "eg2hprivatef.h"
#include "eg2hprivated.h"
#include "eg2herror.h"

/* ////////////////////////////////////////////////////////////////////////// */
int g2h_SymPatchMatrixSize ( int hole_k )
{
  return (G2H_FINALDEG+1)*(G2H_FINALDEG+1)*(6*hole_k+1);
} /*g2h_SymPatchMatrixSize*/

