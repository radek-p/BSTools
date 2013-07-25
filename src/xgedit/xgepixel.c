
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2011                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>

#include "pkvaria.h"
#include "xgedit.h"

#include "xgeprivate.h"

xgecolour_int xge_PixelColourf ( float r, float g, float b )
{
  unsigned long pix;

  r = (float)(r < 0.0 ? 0.0 : ( r > 1.0 ? 1.0 : r ));
  g = (float)(g < 0.0 ? 0.0 : ( g > 1.0 ? 1.0 : g ));
  b = (float)(b < 0.0 ? 0.0 : ( b > 1.0 ? 1.0 : b ));
  pix = (((int)(r*(float)xge_rgbmap.r_bits)) << xge_rgbmap.r_shift) +
        (((int)(g*(float)xge_rgbmap.g_bits)) << xge_rgbmap.g_shift) +
        (((int)(b*(float)xge_rgbmap.b_bits)) << xge_rgbmap.b_shift);
  return pix;
} /*xge_PixelColourf*/

xgecolour_int xge_PixelColour ( byte r, byte g, byte b )
{
  unsigned long pix;

  if (xge_rgbmap.nr_bits < 8) r = (byte)(r >> (8-xge_rgbmap.nr_bits));
  if (xge_rgbmap.ng_bits < 8) g = (byte)(g >> (8-xge_rgbmap.ng_bits));
  if (xge_rgbmap.nb_bits < 8) b = (byte)(b >> (8-xge_rgbmap.nb_bits));
  pix = ((int)r << xge_rgbmap.r_shift) +
        ((int)g << xge_rgbmap.g_shift) +
        ((int)b << xge_rgbmap.b_shift); 
  return pix;
} /*xge_PixelColour*/

void xge_GetPixelColour ( xgecolour_int pixel, byte *r, byte *g, byte *b )
{
  unsigned short rr, gg, bb;

  rr = (unsigned short)((pixel & xgevisual->red_mask)   >> xge_rgbmap.r_shift);
  gg = (unsigned short)((pixel & xgevisual->green_mask) >> xge_rgbmap.g_shift);
  bb = (unsigned short)((pixel & xgevisual->blue_mask)  >> xge_rgbmap.b_shift);
  if ( xge_rgbmap.nr_bits == 8 )
    *r = (byte)rr;
  else
    *r = (byte)((rr << (8-xge_rgbmap.nr_bits)) | (rr << (8-2*xge_rgbmap.nr_bits)));
  if ( xge_rgbmap.ng_bits == 8 )
    *g = (byte)gg;
  else
    *g = (byte)((gg << (8-xge_rgbmap.ng_bits)) | (gg << (8-2*xge_rgbmap.ng_bits)));
  if ( xge_rgbmap.nb_bits == 8 )
    *b = (byte)bb;
  else
    *b = (byte)((bb << (8-xge_rgbmap.nb_bits)) | (bb << (8-2*xge_rgbmap.nb_bits)));
} /*xge_GetPixelColour*/

void xge_GetPixelColourf ( xgecolour_int pixel, float *r, float *g, float *b )
{
  *r = (float)((pixel & xgevisual->red_mask)   >> xge_rgbmap.r_shift)/xge_rgbmap.ar;
  *g = (float)((pixel & xgevisual->green_mask) >> xge_rgbmap.g_shift)/xge_rgbmap.ag;
  *b = (float)((pixel & xgevisual->blue_mask)  >> xge_rgbmap.b_shift)/xge_rgbmap.ab;
} /*xge_GetPixelColourf*/

