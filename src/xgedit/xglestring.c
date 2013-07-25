
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2008, 2010                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>

#include <GL/gl.h>
#include <GL/glx.h>

#include "xgledit.h"
#include "xgeprivate.h"

static boolean xgle_font_ready = false;
static Font font;
static GLuint listBase;

void xgleDrawString ( char *string, int x, int y )
{
  int lgt;

  if ( !xgle_font_ready ) {
    font = XLoadFont ( xgedisplay, "Fixed" );
    if ( !font )
      return;
    listBase = glGenLists ( 128 );
    glXUseXFont ( font, 0, 128, listBase );
    xgle_font_ready = true;
  }
  lgt = strlen ( string );
  glRasterPos2i ( x, y );
  glListBase ( listBase );
  glCallLists ( lgt, GL_BYTE, string );
} /*xgleDrawString*/

