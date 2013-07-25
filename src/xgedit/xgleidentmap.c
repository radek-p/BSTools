
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2008, 2009                            */
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

/* ///////////////////////////////////////////////////////////////////// */
/* This function sets up the OpenGL viewport and geometry transformation */
/* matrices so as to obtain the coordinate system coinciding with the    */
/* XWindow coordinate system for a window, i.e. with (0,0) at the upper  */
/* left corner of the window, y growing downwards. The viewport covers   */
/* the widget rectangle.                                                 */

void xgle_SetIdentMapping ( xge_widget *er )
{
  glViewport ( 0, 0, er->w, er->h );
  glMatrixMode ( GL_PROJECTION );
  glLoadIdentity ();
  glOrtho ( (double)er->x, (double)(er->x+er->w),
            (double)(er->y+er->h), (double)er->y, -1.0, 1.0 );
  glMatrixMode ( GL_MODELVIEW );
  glLoadIdentity ();
} /*xgle_SetIdentMapping*/

