
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009                                  */
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

/* ///////////////////////////////////////////////////////////////////////// */
FILE *bsf_output = NULL;

/* ///////////////////////////////////////////////////////////////////////// */
boolean bsf_OpenOutputFile ( char *filename, boolean append )
{
  if ( append )
    bsf_output = fopen ( filename, "a+" );
  else
    bsf_output = fopen ( filename, "w+" );
  if ( !bsf_output )
    return false;
  return true;
} /*bsf_OpenOutputFile*/

void bsf_CloseOutputFile ( void )
{
  if ( bsf_output ) {
    fclose ( bsf_output );
    bsf_output = NULL;
  }
} /*bsf_CloseOutputFile*/

