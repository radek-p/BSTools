
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "psout.h"

#include "psprivate.h"

boolean _psl_trmk, _psl_btrmk, _psl_htrmk, _psl_bhtrmk;

void _psl_InitPSLib ( void )
{
  _psl_trmk = _psl_btrmk = _psl_htrmk = _psl_bhtrmk = false;
} /*_psl_InitPSLib*/

