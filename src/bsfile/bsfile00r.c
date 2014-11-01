
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

#define INBUFSIZE    1024
#define NAMEBUFSIZE  (BSF_MAX_NAME_LENGTH+1)

/* ///////////////////////////////////////////////////////////////////////// */
typedef struct {
    FILE *bsf_input;
  } bsf_input_struct;

static bsf_input_struct bsf_input_data;
static pkv_scanner      bsf_scanner;

char                    *bsf_name = NULL;
int                     bsf_nextsymbol;
int                     bsf_nextint;
double                  bsf_nextfloat;

const char *bsf_keyword[BSF_NKEYWORDS] = /* this table must be sorted */
  { "BCurve",           /* alphabetically, uppercase before lowercase */
    "BPatch",           /* and the contents must match numbers in     */
    "BSCurve",          /* bsfprivate.h */
    "BSHole",
    "BSMesh",
    "BSPatch",
    "Euler_angles",
    "camera",
    "closed",
    "color",
    "colour",
    "cpoints",
    "cpointsmk",
    "degree",
    "depth",
    "dim",
    "domain",
    "facetmk",
    "facets",
    "frame",
    "halfedgemk",
    "halfedges",
    "ident",
    "knots",
    "knots_u",
    "knots_v",
    "name",
    "parallel",
    "perspective",
    "points",
    "polyline",
    "position",
    "rational",
    "sides",
    "spherical_product",
    "trimmed",
    "uniform",
    "vertices" };

/* ///////////////////////////////////////////////////////////////////////// */
void bsf_GetNextSymbol ( void )
{
  bsf_nextsymbol = pkv_GetNextSymbol ( &bsf_scanner );
  switch ( bsf_nextsymbol ) {
case BSF_SYMB_INTEGER:
    bsf_nextint = bsf_scanner.nextinteger;
    break;
case BSF_SYMB_FLOAT:
    bsf_nextint = bsf_scanner.nextinteger;
    bsf_nextfloat = bsf_scanner.nextfloat;
    break;
case PKV_SYMB_IDENT:
    bsf_nextsymbol = pkv_BinSearchNameTable ( BSF_NKEYWORDS, BSF_FIRST_KEYWORD,
                                              bsf_keyword, bsf_scanner.nextname );
    break;
default:
    break;
  }
} /*bsf_GetNextSymbol*/

void bsf_CloseInputFile ( void )
{
  if ( bsf_input_data.bsf_input ) {
    fclose ( bsf_input_data.bsf_input );
    bsf_input_data.bsf_input = NULL;
    bsf_name = NULL;
  }
  pkv_ShutDownScanner ( &bsf_scanner );
} /*bsf_CloseInputFile*/

static int _bsf_ReadInputFile ( void *userdata, int buflength, char *buffer )
{
  bsf_input_struct *instr;
  int              inbufcount;

  instr = (bsf_input_struct*)userdata;
  inbufcount = fread ( buffer, 1, buflength, instr->bsf_input );
  return inbufcount;
} /*_bsf_ReadInputFile*/

boolean bsf_OpenInputFile ( const char *filename )
{
  bsf_input_data.bsf_input = fopen ( filename, "r+" );
  if ( !bsf_input_data.bsf_input )
    return false;
  if ( pkv_InitScanner ( &bsf_scanner, 0, NULL, BSF_MAX_NAME_LENGTH, NULL,
                         '%', '\n', _bsf_ReadInputFile, &bsf_input_data ) ) {
    bsf_GetNextSymbol ();
    bsf_name = bsf_scanner.nextname;
    return true;
  }
  else {
    bsf_CloseInputFile ();
    return false;
  }
} /*bsf_OpenInputFile*/

void bsf_PrintErrorLocation ( void )
{
  fprintf ( stderr, "Error detected at line %d, column %d\n",
            bsf_scanner.linenum+1, bsf_scanner.colnum+1 );
} /*bsf_PrintErrorLocation*/

