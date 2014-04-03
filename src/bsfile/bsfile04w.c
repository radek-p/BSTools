
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

void bsf_WriteSpaceDim ( int spdimen, boolean rational )
{
  fprintf ( bsf_output, "  %s %d\n",
            bsf_keyword[BSF_SYMB_DIM-BSF_FIRST_KEYWORD], spdimen );
  if ( rational )
    fprintf ( bsf_output, "  %s\n",
              bsf_keyword[BSF_SYMB_RATIONAL-BSF_FIRST_KEYWORD] );
} /*bsf_WriteSpaceDim*/

void bsf_WriteCurveDegree ( int degree )
{
  fprintf ( bsf_output, "  %s %d\n",
            bsf_keyword[BSF_SYMB_DEGREE-BSF_FIRST_KEYWORD], degree );
} /*bsf_WriteCurveDegree*/

void bsf_WritePatchDegree ( int udeg, int vdeg )
{
  fprintf ( bsf_output, "  %s %d %d\n",
            bsf_keyword[BSF_SYMB_DEGREE-BSF_FIRST_KEYWORD], udeg, vdeg );
} /*bsf_WritePatchDegree*/

void bsf_WriteIdent ( int ident )
{
  if ( ident >= 0 )
    fprintf ( bsf_output, "  %s %d\n",
              bsf_keyword[BSF_SYMB_IDENT-BSF_FIRST_KEYWORD], ident );
} /*bsf_WriteIdent*/

