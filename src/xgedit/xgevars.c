
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2014                            */
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

#include "pkgeom.h"
#include "multibs.h"

#include "xgedit.h"
#include "xgeprivate.h"

#ifdef USE_XEXT_SHAPE
#include <X11/extensions/shape.h>
#endif
 
int     xge_winnum = 0;                /* current number of windows */
int     xge_current_win = -1;          /* window currently active */
WinDesc xge_windesc[xge_MAX_WINDOWS];  /* window descriptors */

Display     *xgedisplay;
XVisualInfo *xgevisualinfo;
Colormap    xgecolormap;
int         xgescreen;
Window      xgeroot;
Window      xgewindow;
Pixmap      xgepixmap;
GC          xgegc;
Visual      *xgevisual;
XSizeHints  xgehints;
short       xge_cursnum = 0;            /* current number of cursors */
Cursor      xgecursor[xge_MAX_CURSORS];
KeySym      xgekeysym;

#ifdef USE_XEXT_SHAPE
/* ///////////////////////////////////////////////////////////////////////// */
/* If X Window shape extension is available, widgets may be drawn using      */
/* a special window, which may actually stick beyond the area of the regular */
/* application window. The window has no border and no title bar.            */
boolean    xge_try_ext = true,
           xge_use_specialwin = false,
           xge_specialwin_in_use = false;
SpecialWinDesc xge_specialwin[4] =
  {{None},{None},{None},{None}};
Pixmap     xge_specialpixmap = None;
XRectangle xge_specialpixmaprect;
GC         xge_specialwingc = None, xge_specialpixmapgc = None;
int        xge_specialevent_base, xge_specialerror_base;
#endif

unsigned int xge_current_width, xge_current_height;  /* window dimensions */
unsigned int xge_mouse_buttons = 0;
int        xge_mouse_x, xge_mouse_y;
short      xge_xx, xge_yy;

char       *xge_p_name;
char       xge_done = 0;
short      xge_prevx = -1, xge_prevy = -1;

Colormap   xgecolormap;
int        xge_nplanes;
xgecolour_int xge_foreground, xge_background;
xgecolour_int *xge_palette;
boolean    xge_notinfocus;

xge_widget    *xge_null_widget = NULL;

xge_rgbmap_bits xge_rgbmap;

char **xge_info_msgtext = NULL;

XEvent  xgeevent;

int         xge_win_rect_num = 0;

xge_widget *xge_lastwin = NULL;

int           xge_errmsg_win;
char          *xge_errmsg_msgtext;
xgecolour_int xge_msgbkcolour;

float xge_aspect = 1.0;

int (*xge_callback)(xge_widget*,int,int,short,short);

xge_widget *_xge_background_widget = NULL;
xge_widget *_xge_special_widget = NULL;
int _xge_argc;
char **_xge_argv;

