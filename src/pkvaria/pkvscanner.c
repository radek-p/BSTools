
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
#include <string.h>
#include <ctype.h>

#include "pkvaria.h"
#include "pkvscanner.h"

boolean pkv_InitScanner ( pkv_scanner *sc,
                          int maxnamelength,
                          char bcomment, char ecomment,
                          int (*InputData)(void *usrdata,int buflength,char *buffer),
                          void *userdata )
{
  memset ( sc, 0, sizeof(pkv_scanner) );
  PKV_MALLOC ( sc->inbuffer, PKV_SCANNER_BUFLENGTH );
  PKV_MALLOC ( sc->nextname, maxnamelength+1 )
  if ( !sc->inbuffer || !sc->nextname )
    goto failure;
  sc->inbuflength = PKV_SCANNER_BUFLENGTH;
  sc->namebuflength = maxnamelength+1;
  sc->inbufpos = sc->inbufcount = sc->namebufcount = 0;
  sc->linenum = sc->colnum = 0;
  sc->nextsymbol = PKV_SYMB_NONE;
  sc->bcomment = bcomment;
  sc->ecomment = ecomment;
  sc->InputData = InputData;
  sc->userdata = userdata;
  pkv_GetNextChar ( sc );
  return true;

failure:
  pkv_ShutDownScanner ( sc );
  return false;
} /*pkv_InitScanner*/

void pkv_GetNextChar ( pkv_scanner *sc )
{
  if ( sc->inbufpos >= sc->inbufcount ) {
    sc->inbufcount = sc->InputData ( sc->userdata, sc->inbuflength, sc->inbuffer );
    if ( sc->inbufcount <= 0 ) {
      sc->nextchar = EOF;
      return;
    }
    sc->inbufpos = 0;
  }
  sc->nextchar = (int)sc->inbuffer[sc->inbufpos++];
  if ( sc->nextchar == '\n' )
    sc->linenum ++, sc->colnum = 0;
  else
    sc->colnum ++;
} /*pkv_GetNextChar*/

int pkv_GetNextSymbol ( pkv_scanner *sc )
{
#define ADDTONBUF \
{ if ( sc->namebufcount < sc->namebuflength-1 ) \
    sc->nextname[sc->namebufcount++] = (char)sc->nextchar; \
  else \
    return (sc->nextsymbol = PKV_SYMB_ERROR); \
  pkv_GetNextChar ( sc ); \
}

repeat_scan:
  if ( sc->nextchar == sc->bcomment ) {
    do {
      pkv_GetNextChar ( sc );
    } while ( sc->nextchar != sc->ecomment && sc->nextchar != EOF );
    if ( sc->nextchar != EOF ) {
      pkv_GetNextChar ( sc );
      goto repeat_scan;
    }
  }
  else if ( isdigit ( sc->nextchar ) ) {
        /* integral part */
    sc->nextsymbol = PKV_SYMB_INTEGER;
    sc->namebufcount = 0;
    do {
      ADDTONBUF
    } while ( isdigit ( sc->nextchar ) );
    if ( sc->nextchar == '.' ) {
        /* fractional part */
      ADDTONBUF
      while ( isdigit ( sc->nextchar ) ) {
        ADDTONBUF
      }
      sc->nextsymbol = PKV_SYMB_FLOAT;
    }
    if ( sc->nextchar == 'e' || sc->nextchar == 'E' ) {
        /* exponent */
      ADDTONBUF
      if ( sc->nextchar == '-' || sc->nextchar == '+' ) {
        ADDTONBUF
      }
      if ( !isdigit ( sc->nextchar ) )
        return (sc->nextsymbol = PKV_SYMB_ERROR);
      do {
        ADDTONBUF
      } while ( isdigit ( sc->nextchar ) );
      sc->nextsymbol = PKV_SYMB_FLOAT;
    }
        /* convert the character string to a number */
    sc->nextname[sc->namebufcount] = 0;
    if ( sc->nextsymbol == PKV_SYMB_INTEGER ) {
      sscanf ( sc->nextname, "%d", &sc->nextinteger );
      sc->nextfloat = (double)sc->nextinteger;
    }
    else
      sscanf ( sc->nextname, "%lf", &sc->nextfloat );
  }
  else if ( isalpha ( sc->nextchar ) ) {
        /* a keyword */
    sc->namebufcount = 0;
    do {
      ADDTONBUF
    } while ( isalpha ( sc->nextchar ) || isdigit (sc->nextchar) ||
              sc->nextchar == '_' );
    sc->nextname[sc->namebufcount] = 0;
    sc->nextsymbol = PKV_SYMB_IDENT;
  }
  else {
    switch ( sc->nextchar ) {
  case ' ':
  case '\n':
  case '\t':
      pkv_GetNextChar ( sc );
      goto repeat_scan;
  case '"':
      pkv_GetNextChar ( sc );
      sc->namebufcount = 0;
      do {
        ADDTONBUF
      } while ( sc->nextchar != '"' );
      pkv_GetNextChar ( sc );
      sc->nextsymbol = PKV_SYMB_STRING;
      break;
  case EOF:
      pkv_ShutDownScanner ( sc );
      sc->nextsymbol = PKV_SYMB_EOF;
      break;
  case '+':
      sc->nextsymbol = PKV_SYMB_PLUS;      goto get_next_char;
  case '-':
      sc->nextsymbol = PKV_SYMB_MINUS;     goto get_next_char;
  case '*':
      sc->nextsymbol = PKV_SYMB_STAR;      goto get_next_char;
  case '/':
      sc->nextsymbol = PKV_SYMB_SLASH;     goto get_next_char;
  case '(':
      sc->nextsymbol = PKV_SYMB_LPAREN;    goto get_next_char;
  case ')':
      sc->nextsymbol = PKV_SYMB_RPAREN;    goto get_next_char;
  case '[':
      sc->nextsymbol = PKV_SYMB_LBRACKET;  goto get_next_char;
  case ']':
      sc->nextsymbol = PKV_SYMB_RBRACKET;  goto get_next_char;
  case '{':
      sc->nextsymbol = PKV_SYMB_LBRACE;    goto get_next_char;
  case '}':
      sc->nextsymbol = PKV_SYMB_RBRACE;    goto get_next_char;
  case '<':
      sc->nextsymbol = PKV_SYMB_LANGLE;    goto get_next_char;
  case '>':
      sc->nextsymbol = PKV_SYMB_RANGLE;    goto get_next_char;
  case '=':
      sc->nextsymbol = PKV_SYMB_EQUAL;     goto get_next_char;
  case '\'':
      sc->nextsymbol = PKV_SYMB_CHAR;      goto get_next_char;
  case ',':
      sc->nextsymbol = PKV_SYMB_COMMA;     goto get_next_char;
  case '.':
      sc->nextsymbol = PKV_SYMB_DOT;       goto get_next_char;
  case ';':
      sc->nextsymbol = PKV_SYMB_SEMICOLON; goto get_next_char;
  case '%':
      sc->nextsymbol = PKV_SYMB_PERCENT;   goto get_next_char;
  case '#':
      sc->nextsymbol = PKV_SYMB_HASH;      goto get_next_char;
  case '!':
      sc->nextsymbol = PKV_SYMB_EXCLAM;    goto get_next_char;
  case '~':
      sc->nextsymbol = PKV_SYMB_TILDE;     goto get_next_char;
  case '@':
      sc->nextsymbol = PKV_SYMB_AT;        goto get_next_char;
  case '$':
      sc->nextsymbol = PKV_SYMB_DOLLAR;    goto get_next_char;
  case '&':
      sc->nextsymbol = PKV_SYMB_ET;        goto get_next_char;
  case '\\':
      sc->nextsymbol = PKV_SYMB_BACKSLASH; goto get_next_char;
  case '?':
      sc->nextsymbol = PKV_SYMB_QUESTION;  goto get_next_char;
  case '|':
      sc->nextsymbol = PKV_SYMB_VERTBAR;   goto get_next_char;
  default:
      sc->nextsymbol = PKV_SYMB_OTHER;
get_next_char:
      pkv_GetNextChar ( sc );
      break;
    }
  }

  return sc->nextsymbol;
#undef ADDTONBUF
} /*pkv_GetNextSymbol*/

void pkv_ShutDownScanner ( pkv_scanner *sc )
{
  if ( sc->inbuffer ) PKV_FREE ( sc->inbuffer );
  if ( sc->nextname ) PKV_FREE ( sc->nextname );
  memset ( sc, 0, sizeof(pkv_scanner) );
  sc->nextsymbol = PKV_SYMB_ERROR;
} /*pkv_ShutDownScanner*/

int pkv_BinSearchNameTable ( int nnames, int firstname, const char **names,
                             const char *name )
{
  int i, j, k, l;

  i = k = 0;  j = nnames;
  do {
    k = (i+j)/2;
    l = strcmp ( name, names[k] );
    if ( l == 0 )
      return firstname + k;
    else if ( l < 0 )
      j = k;
    else
      i = k+1;
  } while ( j-i > 0 );
  return PKV_SYMB_ERROR;
} /*pkv_BinSearchNameTable*/

