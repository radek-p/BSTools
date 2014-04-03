
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

#define INBUFSIZE   1024
#define NAMEBUFSIZE   (BSF_MAX_NAME_LENGTH+1)

/* ///////////////////////////////////////////////////////////////////////// */
FILE   *bsf_input = NULL;
static char   *inbuffer = NULL;
char          *bsf_namebuffer = NULL;
static int    inbufpos, inbufcont, namebufcont;
static int    nextchar;
int           bsf_nextsymbol;
int           bsf_nextint;
double        bsf_nextfloat;
static int    bsf_linenum, bsf_colnum; /* for locating errors */

const char *bsf_keyword[BSF_NKEYWORDS] = /* this table must be sorted alphabetically */
  { "BCurve",
    "BPatch",
    "BSCurve",
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
    "facets",
    "frame",
    "halfedges",
    "ident",
    "knots",
    "knots_u",
    "knots_v",
    "name",
    "parallel",
    "perspective",
    "position",
    "rational",
    "sides",
    "uniform",
    "vertices" };

/* ///////////////////////////////////////////////////////////////////////// */
static void GetNextChar ( void )
{
  if ( inbufpos >= inbufcont ) {
    inbufcont = fread ( inbuffer, 1, INBUFSIZE, bsf_input );
    if ( !inbufcont ) {
      nextchar = EOF;
      return;
    }
    inbufpos = 0;
  }
  nextchar = (int)inbuffer[inbufpos++];
  if ( nextchar == '\n' )
    bsf_linenum ++, bsf_colnum = 0;
  else
    bsf_colnum ++;
} /*GetNextChar*/

void bsf_GetNextSymbol ( void )
{
  int i, j, k, l;

#define ADDTONBUF \
{ if ( namebufcont < NAMEBUFSIZE-1 ) \
    bsf_namebuffer[namebufcont++] = (char)nextchar; \
  else { \
    bsf_nextsymbol = BSF_SYMB_ERROR; \
    return; \
  }}

repeat_scan:
  if ( isdigit ( nextchar ) ) {
        /* integral part */
    bsf_nextsymbol = BSF_SYMB_INTEGER;
    namebufcont = 0;
    do {
      ADDTONBUF
      GetNextChar ();
    } while ( isdigit ( nextchar ) );
    if ( nextchar == '.' ) {
        /* fractional part */
      ADDTONBUF
      GetNextChar ();
      while ( isdigit ( nextchar ) ) {
        ADDTONBUF
        GetNextChar ();
      }
      bsf_nextsymbol = BSF_SYMB_FLOAT;
    }
    if ( nextchar == 'e' || nextchar == 'E' ) {
        /* exponent */
      ADDTONBUF
      GetNextChar ();
      if ( nextchar == '-' || nextchar == '+' ) {
        ADDTONBUF
        GetNextChar ();
      }
      if ( !isdigit ( nextchar ) ) {
        bsf_nextsymbol = BSF_SYMB_ERROR;
        return;
      }
      do {
        ADDTONBUF
        GetNextChar ();
      } while ( isdigit ( nextchar ) );
      bsf_nextsymbol = BSF_SYMB_FLOAT;
    }
        /* convert the character string to a number */
    bsf_namebuffer[namebufcont] = 0;
    if ( bsf_nextsymbol == BSF_SYMB_INTEGER ) {
      sscanf ( bsf_namebuffer, "%d", &bsf_nextint );
      bsf_nextfloat = (double)bsf_nextint;
    }
    else
      sscanf ( bsf_namebuffer, "%lf", &bsf_nextfloat );
  }
  else if ( isalpha ( nextchar ) ) {
        /* a keyword */
    namebufcont = 0;
    do {
      ADDTONBUF
      GetNextChar ();
    } while ( isalpha ( nextchar ) || nextchar == '_' );
    bsf_namebuffer[namebufcont] = 0;
        /* binary search in the keyword table */
    i = k = 0;  j = BSF_NKEYWORDS;
    do {
      k = (i+j)/2;
      l = strcmp ( bsf_namebuffer, bsf_keyword[k] );
      if ( l == 0 ) {
        bsf_nextsymbol = BSF_FIRST_KEYWORD + k;
        return;
      }
      else if ( l < 0 )
        j = k;
      else
        i = k+1;
    } while ( j-i > 0 );
    bsf_nextsymbol = BSF_SYMB_ERROR;  /* keyword not found */
  }
  else {
    switch ( nextchar ) {
case ' ':    /* blank space */
case '\t':   /* tab */
case '\n':   /* end of line */
      GetNextChar ();
      goto repeat_scan;

case '%':    /* comment; skip the entire line */
      do {
        GetNextChar ();
      } while ( nextchar != '\n' && nextchar != EOF );
      if ( nextchar != EOF )
        GetNextChar ();
      goto repeat_scan;

case ',':
      bsf_nextsymbol = BSF_SYMB_COMMA;
      GetNextChar ();
      break;

case '+':
      bsf_nextsymbol = BSF_SYMB_PLUS;
      GetNextChar ();
      break;

case '-':
      bsf_nextsymbol = BSF_SYMB_MINUS;
      GetNextChar ();
      break;

case '{':
      bsf_nextsymbol = BSF_SYMB_LBRACE;
      GetNextChar ();
      break;

case '}':
      bsf_nextsymbol = BSF_SYMB_RBRACE;
      GetNextChar ();
      break;

case '"':
      GetNextChar ();
      namebufcont = 0;
      while ( nextchar != '"' && nextchar != EOF ) {
        ADDTONBUF
        GetNextChar ();
      }
      if ( nextchar == '"' ) {
        bsf_namebuffer[namebufcont] = 0;
        GetNextChar ();
        bsf_nextsymbol = BSF_SYMB_STRING;
      }
      else
        bsf_nextsymbol = BSF_SYMB_ERROR;
      break;

case EOF:
      bsf_nextsymbol = BSF_SYMB_EOF;
      break;

default:    /* skip anything else */
      GetNextChar ();
      goto repeat_scan;
    }
  }
} /*bsf_GetNextSymbol*/

static boolean InitInputBuffer ( void )
{
  inbuffer = malloc ( INBUFSIZE );
  bsf_namebuffer = malloc ( NAMEBUFSIZE );
  if ( inbuffer && bsf_namebuffer ) {
    inbufpos = inbufcont = 0;
    bsf_linenum = 1;  bsf_colnum = 0;
    GetNextChar ();
    bsf_GetNextSymbol ();
    return true;
  }
  else {
    if ( inbuffer ) free ( inbuffer );
    if ( bsf_namebuffer ) free ( bsf_namebuffer );
    inbuffer = bsf_namebuffer = NULL;
    return false;
  }
} /*InitInputBuffer*/

void bsf_CloseInputFile ( void )
{
  if ( bsf_input ) {
    fclose ( bsf_input );
    bsf_input = NULL;
  }
  if ( inbuffer ) {
    free ( inbuffer );
    inbuffer = NULL;
  }
  if ( bsf_namebuffer ) {
    free ( bsf_namebuffer );
    bsf_namebuffer = NULL;
  }
} /*bsf_CloseInputFile*/

boolean bsf_OpenInputFile ( const char *filename )
{
  bsf_input = fopen ( filename, "r+" );
  if ( !bsf_input )
    return false;
  if ( InitInputBuffer () )
    return true;
  else {
    bsf_CloseInputFile ();
    return false;
  }
} /*bsf_OpenInputFile*/

void bsf_PrintErrorLocation ( void )
{
  fprintf ( stderr, "Error detected at line %d, column %d\n",
            bsf_linenum, bsf_colnum );
} /*bsf_PrintErrorLocation*/

