
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

#include "calc.h"

  /* lexical symbols */
#define sEND      0
#define sNUMBER   1
#define sPLUS     2
#define sMINUS    3
#define sMULT     4
#define sDIV      5
#define sLPAREN   6
#define sRPAREN   7
#define sNAME     8
#define sCOMMA    9
#define sERROR   10

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

static const char *name[NUMNAMES] = /* this array must be sorted alphabetically */
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

static char *buf;           /* expression string */
static int  pos;            /* current string position */

static int    nextsymb;     /* next symbol */
static double nextnum;      /* next number */
static int    nextint;      /* next integer or function number */
static char   nextname[21]; /* next name string */

static void GetNextSymbol ( void )
{
  double fct;
  int    e, es;
  int    i, j, k, l;

  for (;;)
    if ( isdigit ( buf[pos] ) ) {
            /* integral part */
      nextnum = (double)(buf[pos++]-'0');
      while ( isdigit(buf[pos]) )
        nextnum = 10.0*nextnum + (double)(buf[pos++]-'0');
      if ( buf[pos] == '.' ) { /* fractional part */
        pos ++;
        fct = 0.1;
        while ( isdigit(buf[pos]) ) {
          nextnum += fct*(double)(buf[pos++]-'0');
          fct *= 0.1;
        }
      }
      if ( buf[pos] == 'e' || buf[pos] == 'E' ) { /* exponent */
        pos ++;
        e = 0;
        es = 1;
        if ( buf[pos] == '+' )
          pos ++;
        else if ( buf[pos] == '-' ) {
          pos ++;
          es = -1;
        }
        if ( !isdigit(buf[pos]) ) {
          nextsymb = sERROR;
          return;
        }
        while ( isdigit(buf[pos]) )
          e = 10*e + (buf[pos++]-'0');
        nextnum *= pow ( 10.0, es*e );
      }
      nextint = (int)nextnum;
      nextsymb = sNUMBER;
      return;
    }
    else if ( isalpha ( buf[pos] ) ) {
        /* collect the name characters */
      e = 0;
      do {
        nextname[e++] = buf[pos++];
      } while ( isalpha ( buf[pos] ) && e < 20 );
      nextname[e] = 0;
      nextsymb = sNAME;
        /* search the nametable, binary search */
      i = k = 0;
      j = NUMNAMES;
      do {
        k = (i+j)/2;
        l = strcmp ( nextname, name[k] );
        if ( l == 0 ) {  /* name found */
          nextint = k;
          return;
        }
        else if ( l < 0 )
          j = k;
        else
          i = k+1;
      } while ( j-i > 0 );
      nextint = -1;  /* name not found */
      return;
    }
    else {
      switch ( buf[pos] ) {
    case 0:
        nextsymb = sEND;
        return;
    case '+':
        nextsymb = sPLUS;
        pos ++;
        return;
    case '-':
        nextsymb = sMINUS;
        pos ++;
        return;
    case '*':
        nextsymb = sMULT;
        pos ++;
        return;
    case '/':
        nextsymb = sDIV;
        pos ++;
        return;
    case '(':
        nextsymb = sLPAREN;
        pos ++;
        return;
    case ')':
        nextsymb = sRPAREN;
        pos ++;
        return;
    case ',':
        nextsymb = sCOMMA;
        pos++;
        return;
    case ' ':    /* spaces are skipped - this is the loop purpose */
        pos ++;
        break;
    default:
        nextsymb = sERROR;
        return;
      }
    }
} /*GetNextSymbol*/

static double Expression ( void );

static double EvaluateFunction ( void )
{
  int    fname, nparams;
  double param[MAXPARAMS];

  fname = nextint;
  GetNextSymbol ();
        /* compute the parameters */
  nparams = 0;
  memset ( param, 0, MAXPARAMS*sizeof(double) );
  if ( nextsymb == sLPAREN ) {
    for (;;) {
      GetNextSymbol ();  /* skip the opening parenthesis or comma */
      param[nparams++] = Expression ();
      if ( nextsymb == sRPAREN )
        break;
      else if ( nextsymb != sCOMMA || nparams >= MAXPARAMS ) {
        nextsymb = sERROR;
        return 0.0;
      }
    }
    GetNextSymbol ();  /* skip the closing parenthesis */
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
    nextsymb = sERROR;
    return 0.0;
  }
} /*EvaluateFunction*/

static double Factor ( void )
{
  double c;

  switch ( nextsymb ) {
case sNUMBER:
    c = nextnum;
    GetNextSymbol ();
    return c;
case sNAME:
    c = EvaluateFunction ();
    return c;
case sLPAREN:
    GetNextSymbol ();
    c = Expression ();
    if ( nextsymb == sRPAREN ) {
      GetNextSymbol ();
      return c;
    }
    else
      return 0.0;
default:
    return 0.0;
  }
} /*Factor*/

static double Term ( void )
{
  double s;

  s = Factor ();
  for (;;) {
    switch ( nextsymb ) {
  case sMULT:
      GetNextSymbol ();
      s *= Factor ();
      break;
  case sDIV:
      GetNextSymbol ();
      s /= Factor ();
      break;
  default:
      return s;
    }
  }
} /*Term*/

static double Expression ( void )
{
  double w;
  int zn;

  zn = 1;
  if ( nextsymb == sPLUS )
    GetNextSymbol ();
  else if ( nextsymb == sMINUS ) {
    zn = -1;
    GetNextSymbol ();
  }
  w = zn*Term ();
  for (;;) {
    switch ( nextsymb ) {
  case sPLUS:
      GetNextSymbol ();
      w += Term ();
      break;
  case sMINUS:
      GetNextSymbol ();
      w -= Term ();
      break;
  default:
      return w;
    }
  }
} /*Expression*/

boolean EvaluateExpression ( char *text, double *value )
{
  buf = text;
  pos = 0;
  GetNextSymbol ();
  *value = Expression ();
  return nextsymb == sEND;
} /*EvaluateExpression*/

