
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2008, 2014                            */
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

#define R  GLX_RED_SIZE
#define G  GLX_GREEN_SIZE
#define B  GLX_BLUE_SIZE
#define A  GLX_ALPHA_SIZE
#define AR GLX_ACCUM_RED_SIZE
#define AG GLX_ACCUM_GREEN_SIZE
#define AB GLX_ACCUM_BLUE_SIZE
#define AA GLX_ACCUM_ALPHA_SIZE
#define D  GLX_DEPTH_SIZE
#define DB GLX_DOUBLEBUFFER
#define ST GLX_STENCIL_SIZE


boolean _xgle_no_depth = false, _xgle_no_accum = false, _xgle_no_stencil = false;

static void _xgle_InitAttr ( boolean depth, boolean accum, boolean stencil,
                             int *nattr, int *attrs, int rgb_bits )
{
  int n;

  attrs[2] = attrs[4] = attrs[6] = attrs[8] = rgb_bits;
  n = 9;
  if ( depth ) {
    attrs[n++] = D;
    attrs[n++] = 8;
  }
  if ( accum ) {
    attrs[n++] = AR;
    attrs[n++] = 8;
    attrs[n++] = AG;
    attrs[n++] = 8;
    attrs[n++] = AB;
    attrs[n++] = 8;
    attrs[n++] = AA;
    attrs[n++] = 8;
  }
  if ( stencil ) {
    attrs[n++] = ST;
    attrs[n++] = 1;
  }
  attrs[n] = 0;
  *nattr = n;
} /*_xgle_InitAttr*/

boolean xgle_Init ( int argc, char *argv[],
                    int (*callback)(xge_widget*,int,int,short,short),
                    char *title,
                    boolean depth, boolean accum, boolean stencil )
{
  int     attrs[32] = { GLX_RGBA, R,8, G,8, B,8, A,8, 0 };
  int     rgb_bits[2] = { 8, 5 };
  int     errorBase, eventBase, nattr, i;
  boolean result;

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

  if ( !glXQueryExtension ( xgedisplay, &errorBase, &eventBase ) ) {
    printf ( "Error: GL X extension not supported\n" );
    exit ( 1 );
  }
  xgehints.flags = PPosition | PSize | PMinSize;
  xgehints.height = xge_HEIGHT;
  xgehints.width = xge_WIDTH;
  xgehints.x = 0;
  xgehints.y = 100;
  xgehints.min_height = xge_HEIGHT / 2;
  xgehints.min_width = xge_WIDTH / 2;
        /* negotiate the GL resources */
  result = true;
  _xgle_no_depth = _xgle_no_accum = _xgle_no_stencil = false;
  for ( i = 0; i < 2; i++ ) {
    _xgle_InitAttr ( depth != XGLE_NOT_NEEDED, accum != XGLE_NOT_NEEDED,
                     stencil != XGLE_NOT_NEEDED, &nattr, attrs, rgb_bits[i] );
    if ( (xgevisualinfo = glXChooseVisual ( xgedisplay, xgescreen, attrs )) )
      goto success;
  }
  if ( depth == XGLE_WOULD_BE_NICE ) {
    result = false;  _xgle_no_depth = true;
    for ( i = 0; i < 2; i++ ) {
      _xgle_InitAttr ( false, accum != XGLE_NOT_NEEDED,
                       stencil != XGLE_NOT_NEEDED, &nattr, attrs, rgb_bits[i] );
      if ( (xgevisualinfo = glXChooseVisual ( xgedisplay, xgescreen, attrs )) )
        goto success;
    }
  }
  if ( accum == XGLE_WOULD_BE_NICE ) {
    result = false;  _xgle_no_accum = true;
    for ( i = 0; i < 2; i++ ) {
      _xgle_InitAttr ( depth != XGLE_NOT_NEEDED, false,
                       stencil != XGLE_NOT_NEEDED, &nattr, attrs, rgb_bits[i] );
      if ( (xgevisualinfo = glXChooseVisual ( xgedisplay, xgescreen, attrs )) )
        goto success;
    }
  }
  if ( stencil == XGLE_WOULD_BE_NICE ) {
    result = false;  _xgle_no_stencil = true;
    for ( i = 0; i < 2; i++ ) {
      _xgle_InitAttr ( depth != XGLE_NOT_NEEDED, accum != XGLE_NOT_NEEDED,
                       false, &nattr, attrs, rgb_bits[i] );
      if ( (xgevisualinfo = glXChooseVisual ( xgedisplay, xgescreen, attrs )) )
        goto success;
    }
  }
        /* negotiations failed, no sense to continue */
  printf ( "Error: No matching visual\n" );
  exit ( 1 );

success:
  xgevisual = xgevisualinfo->visual;
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

  xglepixmap = XCreatePixmap ( xgedisplay, xgewindow,
                               xge_MAX_WIDTH, xge_MAX_HEIGHT, xge_nplanes );
  if ( !xglepixmap ) {
    printf ( "Error: Cannot create X pixmap\n" );
    exit ( 1 );
  }
  _xglepixmap = glXCreateGLXPixmap ( xgedisplay, xgevisualinfo, xglepixmap );
  if ( !_xglepixmap ) {
    printf ( "Error: Cannot create GLX pixmap\n" );
    exit ( 1 );
  }
  xglecontext = glXCreateContext ( xgedisplay, xgevisualinfo, NULL, GL_FALSE );
  if ( !xglecontext ) {
    printf ( "Error: Cannot create GL X context " );
    exit ( 1 );
  }
  if ( !glXMakeCurrent ( xgedisplay, _xglepixmap, (GLXContext)xglecontext ) ) {
    printf ( "Error: glXMakeCurrent failed\n" );
    exit ( 1 );
  }
  xgegc = XCreateGC ( xgedisplay, xgewindow, 0, NULL );
  xge_NewCursor ( XC_crosshair );  /* four windows widget */
  xge_NewCursor ( XC_hand1 );      /* list box */
  xge_NewCursor ( XC_pencil );     /* text editor */
  xge_NewCursor ( XC_fleur );
  xge_NewCursor ( XC_arrow );
  xge_NewCursor ( XC_watch );
  xge_NewCursor ( XC_circle );
  xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
  xge_GetWindowSize ();
  _xge_FindAspect ();

  switch ( xgevisual->class )
  {
case TrueColor:
    _xgle_MakePalette ();
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
  return result;
} /*xgle_Init*/

void xgle_Cleanup ( void )
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
} /*xgle_Cleanup*/

