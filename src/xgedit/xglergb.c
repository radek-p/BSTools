
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2008, 2011                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>

#include <GL/gl.h>
#include <GL/glx.h>

#include "pkvaria.h"
#include "xgledit.h"

#include "xgeprivate.h"

#include "xge_rgb.c"

static void parse_colourmask ( int mask, unsigned short *bits,
                               char *nbits, unsigned char *shift,
                               float *a )
{
  unsigned char sh;
  
  if ( !mask )
    *bits = *nbits = *shift = 0;
  else {
    for ( sh = 0;  !(mask & 0x01);  mask = mask >> 1, sh++ )
      ;
    *shift = sh;
    *bits = (unsigned short)mask;
    for ( sh = 0;  mask;  mask = mask >> 1, sh++ )
      ;
    *nbits = sh;
    *a = (float)(pkv_rpower ( 2.0, sh ) - 1.0);
  }
} /*parse_colourmask*/

static boolean done = false;

void _xgle_MakePalette ( void )
{
  xgecolour_int *pal;
  int i, j;
  
  if ( !done ) {
    parse_colourmask ( xgevisual->red_mask, &xge_rgbmap.r_bits,
               &xge_rgbmap.nr_bits, &xge_rgbmap.r_shift, &xge_rgbmap.ar );
    parse_colourmask ( xgevisual->green_mask, &xge_rgbmap.g_bits,
               &xge_rgbmap.ng_bits, &xge_rgbmap.g_shift, &xge_rgbmap.ag );
    parse_colourmask ( xgevisual->blue_mask, &xge_rgbmap.b_bits,
               &xge_rgbmap.nb_bits, &xge_rgbmap.b_shift, &xge_rgbmap.ab );

    xge_palette = pal = (xgecolour_int*)(void*)&rgb_palette;
    for ( i = 0; i < XGLE_PALETTE_LENGTH; i++ ) {
      for ( j = 0; j < 3; j++ )
        xgle_palette[i][j] = (float)rgb_palette[i][j]/255.0;
      *pal++ = xge_PixelColour ( rgb_palette[i][0], rgb_palette[i][1],
                                 rgb_palette[i][2] );
    }
    done = true;
  }
} /*_xgle_MakePalette*/

