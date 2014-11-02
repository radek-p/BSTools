
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>

#include "pkvaria.h"
#include "pkvscanner.h"

#include "calc.h"

  /* built-in function identifiers, must be consecutive and match */
  /* with the contents of the name table */
#define f_PI      0
#define f_ABS     1
#define f_ARCCOS  2
#define f_ARCSIN  3
#define f_ARCTAN  4
#define f_COS     5 
#define f_EXP     6 
#define f_LN      7
#define f_SIN     8
#define f_SQRT    9
#define f_TAN    10
#define NUMNAMES (f_TAN+1)

#define MAXPARAMS 10  /* maximal number of function parameters */

static const char *names[NUMNAMES] = /* this array must be sorted alphabetically */
  { "PI",
    "abs",
    "arccos",
    "arcsin",
    "arctan",
    "cos",
    "exp",
    "ln",
    "sin",
    "sqrt",
    "tan" };

static void GetNextSymbol ( pkv_scanner *sc )
{
  pkv_GetNextSymbol ( sc );
  if ( sc->nextsymbol == PKV_SYMB_IDENT )
    sc->nextinteger = pkv_BinSearchNameTable ( NUMNAMES, 0, names, sc->nextname );
} /*GetNextSymbol*/

static double Expression ( pkv_scanner *sc );

static double EvaluateFunction ( pkv_scanner *sc )
{
  int    fname, nparams;
  double param[MAXPARAMS];

  fname = sc->nextinteger;
  GetNextSymbol ( sc );
        /* compute the parameters */
  nparams = 0;
  memset ( param, 0, MAXPARAMS*sizeof(double) );
  if ( sc->nextsymbol == PKV_SYMB_LPAREN ) {
    for (;;) {
      GetNextSymbol ( sc );  /* skip the opening parenthesis or comma */
      param[nparams++] = Expression ( sc );
      if ( sc->nextsymbol == PKV_SYMB_RPAREN )
        break;
      else if ( sc->nextsymbol != PKV_SYMB_COMMA || nparams >= MAXPARAMS ) {
        sc->nextsymbol = PKV_SYMB_ERROR;
        return 0.0;
      }
    }
    GetNextSymbol ( sc );  /* skip the closing parenthesis */
  }
        /* compute the proper function value */
  switch ( fname ) {
case f_PI:      return PI;
case f_ABS:     return fabs ( param[0] );
case f_ARCCOS:  return acos ( param[0] );
case f_ARCSIN:  return asin ( param[0] );
case f_ARCTAN:  return atan ( param[0] );
case f_COS:     return cos ( param[0] );
case f_EXP:     return exp ( param[0] );
case f_LN:      return log ( param[0] );
case f_SIN:     return sin ( param[0] );
case f_SQRT:    return sqrt ( param[0] );
case f_TAN:     return tan ( param[0] );
default:
    sc->nextsymbol = PKV_SYMB_ERROR;
    return 0.0;
  }
} /*EvaluateFunction*/

static double Factor ( pkv_scanner *sc )
{
  double c;

  switch ( sc->nextsymbol ) {
case PKV_SYMB_INTEGER:
case PKV_SYMB_FLOAT:
    c = sc->nextfloat;
    GetNextSymbol ( sc );
    return c;
case PKV_SYMB_IDENT:
    c = EvaluateFunction ( sc );
    return c;
case PKV_SYMB_LPAREN:
    GetNextSymbol ( sc );
    c = Expression ( sc );
    if ( sc->nextsymbol == PKV_SYMB_RPAREN ) {
      GetNextSymbol ( sc );
      return c;
    }
    else
      return 0.0;
default:
    return 0.0;
  }
} /*Factor*/

static double Term ( pkv_scanner *sc )
{
  double s;

  s = Factor ( sc );
  for (;;) {
    switch ( sc->nextsymbol ) {
  case PKV_SYMB_STAR:
      GetNextSymbol ( sc );
      s *= Factor ( sc );
      break;
  case PKV_SYMB_SLASH:
      GetNextSymbol ( sc );
      s /= Factor ( sc );
      break;
  default:
      return s;
    }
  }
} /*Term*/

static double Expression ( pkv_scanner *sc )
{
  double w;
  int zn;

  zn = 1;
  if ( sc->nextsymbol == PKV_SYMB_PLUS )
    GetNextSymbol ( sc );
  else if ( sc->nextsymbol == PKV_SYMB_MINUS ) {
    zn = -1;
    GetNextSymbol ( sc );
  }
  w = zn*Term ( sc );
  for (;;) {
    switch ( sc->nextsymbol ) {
  case PKV_SYMB_PLUS:
      GetNextSymbol ( sc );
      w += Term ( sc );
      break;
  case PKV_SYMB_MINUS:
      GetNextSymbol ( sc );
      w -= Term ( sc );
      break;
  default:
      return w;
    }
  }
} /*Expression*/

typedef struct {
    boolean all;
  } expr_data;

static int InputExprText ( void *userdata, int buflength, char *buffer )
{
  expr_data *mydata;

  mydata = (expr_data*)userdata;
  if ( mydata->all )
    return 0;
  else {  /* the text is already in the buffer */
    mydata->all = true;
    return strlen ( buffer );
  }
} /*InputExprText*/

boolean EvaluateExpression ( char *text, double *value )
{
  pkv_scanner sc;
  expr_data   mydata;
  char        name[65];

  mydata.all = false;
  if ( !pkv_InitScanner ( &sc, strlen(text), text, 64, name, 256, 256,
                          InputExprText, (void*)&mydata ) )
    return false;
  GetNextSymbol ( &sc );
  *value = Expression ( &sc );
  GetNextSymbol ( &sc );
  return sc.nextsymbol == PKV_SYMB_EOF;
} /*EvaluateExpression*/

