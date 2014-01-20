
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

void xge_Init ( int argc, char *argv[],
                int (*callback)(xge_widget*,int,int,short,short),
                char *title )
{
  _xge_argc = argc;
  _xge_argv = argv;

  xge_p_name = argv[0];
  xgedisplay = XOpenDisplay ( "" );
  if ( !xgedisplay ) {
    printf ( "Error: Could not connect to the X server\n" );
    exit ( 1 );
  }
  xgescreen = DefaultScreen ( xgedisplay );
  xgeroot = RootWindow ( xgedisplay, xgescreen );
  xge_callback = callback;

  xgehints.flags = PPosition | PSize | PMinSize;
  xgehints.height = xge_HEIGHT;
  xgehints.width = xge_WIDTH;
  xgehints.x = 0;
  xgehints.y = 100;
  xgehints.min_height = xge_HEIGHT / 2;
  xgehints.min_width = xge_WIDTH / 2;
  xgevisual = XDefaultVisual ( xgedisplay, xgescreen );
  xge_nplanes = XDisplayPlanes ( xgedisplay, xgescreen );
  xgecolormap = XCreateColormap ( xgedisplay,
                      DefaultRootWindow(xgedisplay),
                      xgevisual, AllocNone );
  if ( !xgecolormap ) {
    printf ( "Error: Cannot create colormap\n" );
    exit ( 1 );
  }

  if ( !title )
    title = xge_p_name;
  xge_NewWindow ( title );
  xge_SetWindow ( 0 );

  xgegc = XCreateGC ( xgedisplay, xgewindow, 0, 0 );
  xge_NewCursor ( XC_crosshair );  /* 0 four windows widget */
  xge_NewCursor ( XC_hand1 );      /* 1 list box */
  xge_NewCursor ( XC_pencil );     /* 2 text editor */
  xge_NewCursor ( XC_fleur );      /* 3 */
  xge_NewCursor ( XC_draft_large /* XC_arrow */ );      /* 4 */
  xge_NewCursor ( XC_watch );      /* 5 */
  xge_NewCursor ( XC_circle );     /* 6 */
  xge_NewCursor ( XC_left_ptr );   /* 7 default arrow cursor */
  xge_NewCursor ( -1 );            /* 8 invisible (pixmap) cursor */
  xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
  xge_GetWindowSize ();
  _xge_FindAspect ();

  switch ( xgevisual->class )
  {
case TrueColor:
    _xge_MakePalette ();
    break;
default:
    printf ( "Error: True colour visual is needed." );
    exit ( 1 );
  }
  xge_InitRectAllocation ();
  xge_errmsg_edr = xge_NewEmptyWidget ( -1, NULL, -1, 1, 1, -1, -1 );
  _xge_background_widget = xge_NewWidget ( -1, NULL, -1, 1, 1, -1, -1, NULL, NULL,
                                           _xge_background_msg, xge_DrawEmpty );
  _xge_special_widget = xge_NewEmptyWidget ( -1, NULL, -1, 1, 1, -1, -1 );
  xge_null_widget  = xge_NewEmptyWidget ( -1, NULL, -1, 1, 1, -1, -1 );
} /*xge_Init*/

void xge_Cleanup ( void )
{
  int i;

#ifdef USE_XEXT_SHAPE
  _xge_DestroySpecialWin ();
#endif
  XFreeGC ( xgedisplay, xgegc );
  for ( i = xge_winnum-1; i >= 0; i-- ) {
    xge_SetWindow ( i );
    XFreePixmap ( xgedisplay, xgepixmap );
    XDestroyWindow ( xgedisplay, xgewindow );
  }
  XCloseDisplay ( xgedisplay );
} /*xge_Cleanup*/

