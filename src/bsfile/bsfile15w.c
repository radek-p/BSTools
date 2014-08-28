
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

/* ///////////////////////////////////////////////////////////////////////// */
boolean bsf_newline;
int     bsf_current_length;
int     bsf_current_indentation;

/* ///////////////////////////////////////////////////////////////////////// */
void bsf_EndOutputLine ( void )
{
  if ( bsf_current_length > 0 )
    fprintf ( bsf_output, "\n" );
  bsf_current_length = 0;
  bsf_newline = true;
} /*bsf_EndOutputLine*/

void bsf_DivideOutputLine ( void )
{
  if ( bsf_current_length >= BSF_OUTPUT_LINE_LENGTH+bsf_current_indentation )
    bsf_EndOutputLine ();
  else {
    bsf_current_length += fprintf ( bsf_output, " " );
    bsf_newline = false;
  }
} /*bsf_DivideOutputLine*/

void bsf_WriteCurrentIndentation ( void )
{
  if ( bsf_newline ) {
    if ( bsf_current_indentation > 0 ) {
      bsf_current_length +=
        fprintf ( bsf_output, "%*s", bsf_current_indentation, "" );
      bsf_newline = false;
    }
  }
} /*bsf_WriteCurrentIndentation*/

