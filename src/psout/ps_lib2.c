
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "psout.h"

#include "psprivate.h"

void psl_TrMark ( float x, float y )
{
  char sx[8], sy[8];
  char s[21];

  if ( !_psl_trmk ) {
    ps_Write_Command("/trm {");
    ps_Write_Command("  /y exch def");
    ps_Write_Command("  /x exch def");
    ps_Write_Command("  1 setgray");
    ps_Write_Command(
      "  newpath x y moveto 15 -35 rlineto -30 0 rlineto closepath fill");
    ps_Write_Command("  2 setlinewidth");
    ps_Write_Command("  0 setgray");
    ps_Write_Command(
      "  newpath x y moveto 15 -35 rlineto -30 0 rlineto closepath stroke");
    ps_Write_Command("} def");
    _psl_trmk = true;
  }
  sprintf(sx, "%7.2f", x);
  sprintf(sy, "%7.2f", y);
  sprintf(s, "%s %s trm", sx, sy);
  ps_Write_Command(s);
}  /*psl_TrMark*/

void psl_BlackTrMark ( float x, float y )
{
  char sx[8], sy[8];
  char s[23];

  if ( !_psl_btrmk ) {
    ps_Write_Command("/btrm {");
    ps_Write_Command("  /y exch def");
    ps_Write_Command("  /x exch def");
    ps_Write_Command(
      "  newpath x y moveto 15 -35 rlineto -30 0 rlineto closepath fill");
    ps_Write_Command("} def");
    _psl_btrmk = true;
  }
  sprintf(sx, "%7.2f", x);
  sprintf(sy, "%7.2f", y);
  sprintf(s, "%s %s btrm", sx, sy);
  ps_Write_Command(s);
}  /*psl_BlackTrMark*/

void psl_HighTrMark ( float x, float y )
{
  char sx[8], sy[8];
  char s[21];

  if ( !_psl_htrmk ) {
    ps_Write_Command("/htrm {");
    ps_Write_Command("  /y exch def");
    ps_Write_Command("  /x exch def");
    ps_Write_Command("  1 setgray");
    ps_Write_Command(
      "  newpath x y moveto 15 -70 rlineto -30 0 rlineto closepath fill");
    ps_Write_Command("  2 setlinewidth");
    ps_Write_Command("  0 setgray");
    ps_Write_Command(
      "  newpath x y moveto 15 -70 rlineto -30 0 rlineto closepath stroke");
    ps_Write_Command("} def");
    _psl_htrmk = true;
  }
  sprintf(sx, "%7.2f", x);
  sprintf(sy, "%7.2f", y);
  sprintf(s, "%s %s htrm", sx, sy);
  ps_Write_Command(s);
}  /*psl_HighTrMark*/

void psl_BlackHighTrMark ( float x, float y )
{
  char sx[8], sy[8];
  char s[23];

  if ( !_psl_bhtrmk ) {
    ps_Write_Command("/bhtrm {");
    ps_Write_Command("  /y exch def");
    ps_Write_Command("  /x exch def");
    ps_Write_Command(
      "  newpath x y moveto 15 -70 rlineto -30 0 rlineto closepath fill");
    ps_Write_Command("} def");
    _psl_bhtrmk = true;
  }
  sprintf(sx, "%7.2f", x);
  sprintf(sy, "%7.2f", y);
  sprintf(s, "%s %s bhtrm", sx, sy);
  ps_Write_Command(s);
}  /*psl_BlackHighTrMark*/

void psl_LTrMark ( float t )
{
  float x, y;

  psl_GetPointf ( t, &x, &y );
  psl_TrMark ( x, y );
}  /*psl_LTrMark*/

void psl_BlackLTrMark ( float t )
{
  float x, y;

  psl_GetPointf ( t, &x, &y );
  psl_BlackTrMark ( x, y );
}  /*psl_BlackLTrMark*/

void psl_HighLTrMark ( float t )
{
  float x, y;

  psl_GetPointf ( t, &x, &y );
  psl_HighTrMark ( x, y );
}  /*psl_HighLTrMark*/

void psl_BlackHighLTrMark ( float t )
{
  float x, y;

  psl_GetPointf ( t, &x, &y );
  psl_BlackHighTrMark ( x, y );
}  /*psl_BlackHighLTrMark*/

