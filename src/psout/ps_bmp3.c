
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
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

void ps_Init_BitmapRGB ( int w, int h, int x, int y )
{
  _ps_bmp_w = (unsigned short)w;
  _ps_bmp_p = 8;
  fprintf ( _ps_f, "gsave\n");
  fprintf ( _ps_f, "/picstr 3 string def\n");
  fprintf ( _ps_f, "/displayimagergb {\n");
  fprintf ( _ps_f, "  %d %d 8 [ %d 0 0 %d 0 %d ]\n", w, h, w, -h, h);
  fprintf ( _ps_f, "  { currentfile picstr readhexstring pop }\n");
  fprintf ( _ps_f, "  false 3 colorimage\n");
  fprintf ( _ps_f, "} def\n");
  fprintf ( _ps_f, "%d %d scale\n", w, h);
  fprintf ( _ps_f, "%d %d div %d %d div translate\n", x, w, y, h);
  fprintf ( _ps_f, "displayimagergb\n");
  _ps_bmp_y = (short)h;
  _ps_ExtendRect ( (float)x, (float)y );
  _ps_ExtendRect ( (float)(x + w), (float)(y + h) );
} /*ps_Init_BitmapRGB*/

void ps_Out_LineRGB ( byte *data )
{
  int x;
  char s[2];

  if ( _ps_bmp_y <= 0 )
    return;
  for ( x = 0; x < 3*_ps_bmp_w; x++ ) {
    pkv_HexByte ( data[x], s );
    fputc ( s[0], _ps_f );  fputc ( s[1], _ps_f );
  }
  fputc ( '\n', _ps_f );
  _ps_bmp_y--;
  if ( _ps_bmp_y != 0 )
    return;
  _ps_bmp_y = -1;
  fprintf ( _ps_f, "grestore\n" );
} /*ps_Out_LineRGB*/

