
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2010                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <memory.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "psout.h"

#include "psprivate.h"

void ps_Set_Gray ( float gray )
{
  gray = (float)max ( 0.0, gray );  gray = (float)min ( 1.0, gray );
  if ( gray != _ps_cgray ) {
    fprintf ( _ps_f, "%5.3f setgray\n", gray );
    _ps_cgray = gray;
    _ps_cred = _ps_cgreen = _ps_cblue = -1.0;
  }
} /*ps_Set_Gray*/

void ps_Set_RGB ( float red, float green, float blue )
{
  red = (float)max ( 0.0, red );      red = (float)min ( 1.0, red );
  green = (float)max ( 0.0, green );  green = (float)min ( 1.0, green );
  blue = (float)max ( 0.0, blue );    blue = (float)min ( 1.0, blue );
  if ( red != _ps_cred || green != _ps_cgreen || blue != _ps_cblue ) {
    fprintf ( _ps_f, "%5.3f %5.3f %5.3f setrgbcolor\n", red, green, blue );
    _ps_cred = red;  _ps_cgreen = green;  _ps_cblue = blue;
    _ps_cgray = -1.0;
  }
} /*ps_Set_RGB*/

void ps_Set_Line_Width ( float w )
{
  if (w != _ps_cwidth) {
    fprintf(_ps_f, "%5.*f setlinewidth\n", ps_dec_digits, w);
    _ps_cwidth = w;
  }
} /*ps_Set_Line_Width*/

void ps_Draw_Line ( float x1, float y1, float x2, float y2 )
{
  _ps_OutProc(_drawline);
  fprintf(_ps_f, "%5.*f %5.*f %5.*f %5.*f _dl\n",
	  ps_dec_digits, x1, ps_dec_digits, y1, ps_dec_digits, x2,
	  ps_dec_digits, y2);
  _ps_ExtendRect(x1, y1);
  _ps_ExtendRect(x2, y2);
} /*ps_Draw_Line*/

void ps_Set_Clip_Rect ( float w, float h, float x, float y )
{
  _ps_OutProc(_moveto);
  _ps_OutProc(_lineto);
  fprintf(_ps_f,
    "%5.*f %5.*f %5.*f %5.*f %5.*f %5.*f %5.*f %5.*f newpath _mt _lt _lt _lt closepath clip\n",
    ps_dec_digits, x, ps_dec_digits, y, ps_dec_digits, x + w, ps_dec_digits,
    y, ps_dec_digits, x + w, ps_dec_digits, y + h, ps_dec_digits, x,
    ps_dec_digits, y + h);
} /*ps_Set_Clip_Rect*/

void ps_Draw_Rect ( float w, float h, float x, float y )
{
  _ps_OutProc(_drawrect);
  fprintf(_ps_f, "%5.*f %5.*f %5.*f %5.*f _dr\n",
	  ps_dec_digits, w, ps_dec_digits, h, ps_dec_digits, x, ps_dec_digits,
	  y);
  _ps_ExtendRect(x, y);
  _ps_ExtendRect(x + w, y + h);
} /*ps_Draw_Rect*/

void ps_Fill_Rect ( float w, float h, float x, float y )
{
  _ps_OutProc(_fillrect);
  fprintf(_ps_f, "%5.*f %5.*f %5.*f %5.*f _fr\n",
	  ps_dec_digits, w, ps_dec_digits, h, ps_dec_digits, x, ps_dec_digits,
	  y);
  _ps_ExtendRect(x, y);
  _ps_ExtendRect(x + w, y + h);
} /*ps_Draw_Rect*/

void ps_Hatch_Rect ( float w, float h, float x, float y, float ang, float d )
{
  _ps_OutProc ( _hatchrect );
  fprintf(_ps_f, "%5.*f %5.*f %5.*f %5.*f %5.*f %5.*f _hatchrect\n",
          ps_dec_digits, w, ps_dec_digits, h,
          ps_dec_digits, x, ps_dec_digits, y,
          ps_dec_digits, ang*(180.0/PI), ps_dec_digits, d);
  _ps_ExtendRect ( x, y );
  _ps_ExtendRect ( x+w, y+h );
} /*ps_Hatch_Rect*/

