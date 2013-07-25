
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2009                            */
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

