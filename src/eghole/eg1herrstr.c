
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

#include "eg1holef.h"
#include "eg1herror.h"

char *_g1h_GetErrorString ( int errorcode )
{
  switch ( errorcode ) {
case 0:
    return "Diagnostic of this error still ought to be improved";

case G1H_ERROR_NO_SCRATCH_MEMORY:
    return "Not enough scratch memory";

case G1H_ERROR_CANNOT_MALLOC:
    return "Cannot malloc";

case G1H_ERROR_INVALID_OPTION:
    return "Invalid option";

case G1H_ERROR_INVALID_PARTITION:
    return "Invalid partition";

case G1H_ERROR_INVALID_JUNC_FUNC:
    return "Invalid junction functions";

case G1H_ERROR_NONPOSITIVE_MATRIX:
    return "System matrix not positive-definite";

case G1H_ERROR_NONPOSITIVE_EXT_MATRIX:
    return "Extended system matrix not positive-definite";

case G1H_ERROR_UNDEFINED_CONSTR:
    return "Undefined constraints";

case G1H_ERROR_INCONSISTENT_CONSTR:
    return "Inconsistent constraints";

case G1H_ERROR_NL_CANNOT_PROJECT:
    return "NL: Cannot project this surface";

case G1H_ERROR_NL_JACOBIAN:
    return "NL: Jacobian has various signs";

case G1H_ERROR_NL_MINIMIZATION:
    return "NL: Cannot find the minimum";

case G1H_ERROR_SPLINE_BASIS_NOT_READY:
    return "Spline basis not ready";

default:
    return "If there was an error, it is not properly signalled yet";
  }
} /*_g1h_GetErrorString*/

