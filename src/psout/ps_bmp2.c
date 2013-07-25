
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

void ps_Init_BitmapP ( int w, int h, int x, int y )
{
  _ps_bmp_w = (unsigned short)w;
  _ps_bmp_p = 8;
  _ps_bmp_y = (short)h;
  fprintf ( _ps_f, "gsave\n");
  fprintf ( _ps_f, "/picstr 1 string def\n");
  fprintf ( _ps_f, "/count 0 def\n");
  fprintf ( _ps_f, "/pack false def\n");
  fprintf ( _ps_f, "/uncomp {\n");
  fprintf ( _ps_f, "  count 0 eq\n");
  fprintf ( _ps_f, "  {\n");
  fprintf ( _ps_f,
    "    currentfile picstr readhexstring pop pop /count picstr 0 get def\n");
  fprintf ( _ps_f, "    count 127 gt\n");
  fprintf ( _ps_f, "    { /pack true def /count count 128 sub def }\n");
  fprintf ( _ps_f, "    { /pack false def }\n");
  fprintf ( _ps_f, "    ifelse\n");
  fprintf ( _ps_f, "    currentfile picstr readhexstring pop\n");
  fprintf ( _ps_f, "  }\n");
  fprintf ( _ps_f, "  {\n");
  fprintf ( _ps_f, "    /count count 1 sub def\n");
  fprintf ( _ps_f, "    pack\n");
  fprintf ( _ps_f, "    { picstr }\n");
  fprintf ( _ps_f, "    { currentfile picstr readhexstring pop }\n");
  fprintf ( _ps_f, "    ifelse\n");
  fprintf ( _ps_f, "  }\n");
  fprintf ( _ps_f, "  ifelse\n");
  fprintf ( _ps_f, "} bind def\n");
  fprintf ( _ps_f, "/paint {\n");
  fprintf ( _ps_f, "  %d %d %d [%d 0 0 %d 0 %d] {uncomp} image\n",
	    w, h, 8, w, -h, h);
  fprintf ( _ps_f, "} bind def\n");
  fprintf ( _ps_f, "%d %d scale\n", w, h);
  fprintf ( _ps_f, "%d %d div %d %d div translate\n", x, w, y, h);
  fprintf ( _ps_f, "paint\n");
} /*ps_Init_BitmapP*/

void ps_Out_LineP ( byte *data )
{
  short x;
  short p = 0;
  short i;
  byte b;
  char s[2];
  short count;

  if (_ps_bmp_y <= 0)
    return;
  while (p < _ps_bmp_w) {
    b = data[p];
    x = (short)(p + 1);
    if (x == _ps_bmp_w) {
      fprintf(_ps_f, "00");
      pkv_HexByte(b, s);
      fputc ( s[0], _ps_f );  fputc ( s[1], _ps_f );
    } else if (data[x] != b) {  /* unpacked run */
      count = 1;
      x++;
      while (x < _ps_bmp_w && data[x] != data[x-1] && count < 127) {
	x++;
	count++;
      }
      if (x < _ps_bmp_w && data[x] == data[x-1]) {
	x--;
	count--;
      }
      pkv_HexByte((byte)count, s);
      fputc ( s[0], _ps_f );  fputc ( s[1], _ps_f );
      for (i = p; i <= x - 1; i++) {
	pkv_HexByte(data[i], s);
        fputc ( s[0], _ps_f );  fputc ( s[1], _ps_f );
      }
    } else {
      count = 1;
      x++;
      while (x < _ps_bmp_w && data[x] == b && count < 127) {
	count++;
	x++;
      }
      pkv_HexByte((byte)(count + 0x80), s);
      fputc ( s[0], _ps_f );  fputc ( s[1], _ps_f );
      pkv_HexByte(b, s);
      fputc ( s[0], _ps_f );  fputc ( s[1], _ps_f );
      /*  packed run */
    }
    p = x;
  }
  fputc ( '\n', _ps_f );
  _ps_bmp_y--;
  if (_ps_bmp_y == 0) {
    _ps_bmp_y = -1;
    fprintf(_ps_f, "grestore\n");
  }
} /*ps_Out_LineP*/

