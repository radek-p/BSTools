
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camerad.h"
#include "multibs.h"
#include "bsfile.h"

#include "bsfprivate.h"

void bsf_WriteDependencies ( int depname, int ndep, const int *dep )
{
  int i, k;
  int depid;

  switch ( depname ) {
case BSF_DEP_SPHERICAL:
    depid = BSF_SYMB_SPHERICAL_PRODUCT;
    break;
default:
    return;
  }
/* ndep must be at least 1 */
  BSFwci
  bsf_current_length += fprintf ( bsf_output,
          "%s { ", bsf_keyword[depid-BSF_FIRST_KEYWORD] );
  for ( i = k = 0;  i < ndep-1;  i++, k++ ) {
    bsf_current_length += fprintf ( bsf_output, "%d,", dep[i] );
    BSFdol
    BSFwci
  }
  bsf_current_length += fprintf ( bsf_output, "%d }", dep[ndep-1] );
  BSFeol
} /*bsf_WriteDependencies*/

