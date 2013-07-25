
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

void ps_Init_BitmapRGBP ( int w, int h, int x, int y )
{
  _ps_bmp_w = (unsigned short)w;
  _ps_bmp_p = 8;
  fprintf ( _ps_f, "gsave\n");
  fprintf ( _ps_f, "/countstr 1 string def\n");
  fprintf ( _ps_f, "/picstr 3 string def\n");
  fprintf ( _ps_f, "/count 0 def\n");
  fprintf ( _ps_f, "/pack false def\n");
  fprintf ( _ps_f, "/uncomprgb {\n");
  fprintf ( _ps_f, "  count 0 eq\n");
  fprintf ( _ps_f, "  { currentfile countstr readhexstring pop pop\n");
  fprintf ( _ps_f, "    /count countstr 0 get def\n");
  fprintf ( _ps_f, "    count 127 gt\n");
  fprintf ( _ps_f, "    { /pack true def /count count 128 sub def }\n");
  fprintf ( _ps_f, "    { /pack false def }\n");
  fprintf ( _ps_f, "    ifelse\n");
  fprintf ( _ps_f, "    currentfile picstr readhexstring pop\n");
  fprintf ( _ps_f, "  }\n");
  fprintf ( _ps_f, "  { /count count 1 sub def\n");
  fprintf ( _ps_f, "    pack { picstr } { currentfile picstr readhexstring pop }\n");
  fprintf ( _ps_f, "    ifelse\n");
  fprintf ( _ps_f, "  }\n");
  fprintf ( _ps_f, "  ifelse\n");
  fprintf ( _ps_f, "} bind def\n");
  fprintf ( _ps_f, "/displayimagergb {\n");
  fprintf ( _ps_f, "  %d %d 8 [ %d 0 0 %d 0 %d ]\n", w, h, w, -h, h);
  fprintf ( _ps_f, "  { uncomprgb }\n");
  fprintf ( _ps_f, "  false 3 colorimage\n");
  fprintf ( _ps_f, "} def\n");
  fprintf ( _ps_f, "%d %d scale\n", w, h);
  fprintf ( _ps_f, "%d %d div %d %d div translate\n", x, w, y, h);
  fprintf ( _ps_f, "displayimagergb\n");
  _ps_bmp_y = (short)h;
  _ps_ExtendRect ( (float)x, (float)y );
  _ps_ExtendRect ( (float)(x + w), (float)(y + h) );
} /*ps_Init_BitmapRGBP*/

void ps_Out_LineRGBP ( byte *data )
{
  int x, i, count;
  int p = 0;
  byte r, g, b;
  char s[2];

  if ( _ps_bmp_y <= 0 )
    return;
  while (p < 3*_ps_bmp_w) {
    r = data[p];  g = data[p+1];  b = data[p+2];
    x = p + 3;
    if (x == 3*_ps_bmp_w) {
      fprintf ( _ps_f, "00");
      pkv_HexByte ( r, s );  fputc ( s[0], _ps_f );  fputc ( s[1], _ps_f );
      pkv_HexByte ( g, s );  fputc ( s[0], _ps_f );  fputc ( s[1], _ps_f );
      pkv_HexByte ( b, s );  fputc ( s[0], _ps_f );  fputc ( s[1], _ps_f );
    }
    else if (data[x] != r || data[x+1] != g || data[x+2] != b) {  /* unpacked run */
      count = 1;
      x += 3;
      while (x < 3*_ps_bmp_w &&
             (data[x] != data[x-3] || data[x+1] != data[x-2] || data[x+2] != data[x-1] ) &&
             count < 127) {
	x += 3;
	count++;
      }
      if (x < 3*_ps_bmp_w &&
          data[x] == data[x-3] && data[x+1] == data[x-2] && data[x+2] == data[x-1] ) {
	x -= 3;
	count--;
      }
      pkv_HexByte ( (byte)count, s) ; fputc ( s[0], _ps_f );  fputc ( s[1], _ps_f );
      for (i = p; i <= x - 1; i++) {
	pkv_HexByte ( data[i], s );  fputc ( s[0], _ps_f );  fputc ( s[1], _ps_f );
      }
    } else {
      count = 1;
      x += 3;
      while (x < 3*_ps_bmp_w &&
             data[x] == r && data[x+1] == g && data[x+2] == b &&
             count < 127) {
	count++;
	x += 3;
      }
      pkv_HexByte((byte)(count + 0x80), s);  fputc ( s[0], _ps_f );  fputc ( s[1], _ps_f );
      pkv_HexByte ( r, s );  fputc ( s[0], _ps_f );  fputc ( s[1], _ps_f );
      pkv_HexByte ( g, s );  fputc ( s[0], _ps_f );  fputc ( s[1], _ps_f );
      pkv_HexByte ( b, s );  fputc ( s[0], _ps_f );  fputc ( s[1], _ps_f );
      /*  packed run */
    }
    p = x;
  }
  fputc ( '\n', _ps_f );
  _ps_bmp_y--;
  if (_ps_bmp_y == 0) {
    _ps_bmp_y = -1;
    fprintf ( _ps_f, "grestore\n");
  }
} /*ps_Out_LineRGBP*/

