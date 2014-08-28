
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009, 2014                            */
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

void bsf_WriteSpaceDim ( int spdimen, boolean rational )
{
  BSFwci
  bsf_current_length +=
          fprintf ( bsf_output, "%s %d",
                    bsf_keyword[BSF_SYMB_DIM-BSF_FIRST_KEYWORD], spdimen );
  if ( rational )
    bsf_current_length +=
          fprintf ( bsf_output, " %s",
                    bsf_keyword[BSF_SYMB_RATIONAL-BSF_FIRST_KEYWORD] );
  BSFeol
} /*bsf_WriteSpaceDim*/

void bsf_WriteCurveDegree ( int degree )
{
  BSFwci
  bsf_current_length +=
          fprintf ( bsf_output, "%s %d",
                    bsf_keyword[BSF_SYMB_DEGREE-BSF_FIRST_KEYWORD], degree );
  BSFeol
} /*bsf_WriteCurveDegree*/

void bsf_WritePatchDegree ( int udeg, int vdeg )
{
  BSFwci
  bsf_current_length +=
          fprintf ( bsf_output, "%s %d %d",
                    bsf_keyword[BSF_SYMB_DEGREE-BSF_FIRST_KEYWORD], udeg, vdeg );
  BSFeol
} /*bsf_WritePatchDegree*/

void bsf_WriteIdent ( int ident )
{
  if ( ident >= 0 ) {
    BSFwci
    bsf_current_length +=
          fprintf ( bsf_output, "%s %d",
                    bsf_keyword[BSF_SYMB_IDENT-BSF_FIRST_KEYWORD], ident );
    BSFeol
  }
} /*bsf_WriteIdent*/

