
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
#include "eg2holed.h"

#include "eg2hprivated.h"
#include "eg2herror.h"

int g2h_GetErrorCoded ( GHoleDomaind *domain, char **ErrorString )
{
  int code;

  code = domain->error_code;
  if ( ErrorString )
    *ErrorString = _g2h_GetErrorString ( code );
  return code;
} /*g2h_GetErrorCoded*/

