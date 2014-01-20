
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
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

#include "xgedit.h"
#include "xgeprivate.h"

/* ///////////////////////////////////////////////////////////////////////// */
/* If X Window shape extension is available, the widget will be drawn using  */
/* a special window, which may actually stick beyond the area of the regular */
/* application window. The window has no border and no title bar.            */

#ifdef USE_XEXT_SHAPE
void (*_xge_compspecialwinsizes)( xge_widget *wdg ) = NULL;
void (*_xge_maskspecialwin)(xge_widget *wdg) = NULL;
void (*_xge_drawspecialwin)(int w, xge_widget *wdg) = NULL;
xge_widget *_xge_specialwdg = NULL;

boolean _xge_CreateSpecialWin ( void )
{
  XSetWindowAttributes attr;
  XGCValues            gcv;
  int                  i;

  if ( xge_try_ext ) {
    xge_try_ext = false;
    if ( !XShapeQueryExtension ( xgedisplay, &xge_specialevent_base,
                                 &xge_specialerror_base ) ) {
      xge_use_specialwin = false;
      return false;
    }
    xge_specialpixmap = None;
    xge_specialwingc = xge_specialpixmapgc = None;
    attr.background_pixmap = None;
    attr.save_under = True;
    attr.override_redirect = True;
    attr.event_mask = ExposureMask;
    for ( i = 0; i < 4; i++ ) {
      xge_specialwin[i].thewindow = XCreateWindow ( xgedisplay,
                       DefaultRootWindow(xgedisplay), 0, 0, 1, 1, 0,
                       CopyFromParent, InputOutput, CopyFromParent,
                       CWBackPixmap | CWSaveUnder | CWEventMask | CWOverrideRedirect,
                       &attr );
      if ( xge_specialwin[i].thewindow == None )
        goto failure;
      xge_specialwin[i].nonempty = xge_specialwin[i].mapped = false;
    }
    xge_use_specialwin = true;
    gcv.foreground = 1;
    gcv.line_width = 1;
    gcv.line_style = LineSolid;
    xge_specialwingc = XCreateGC ( xgedisplay, xge_specialwin[0].thewindow,
                                   GCForeground | GCLineWidth | GCLineStyle,
                                   &gcv );
    xge_specialpixmap = XCreatePixmap ( xgedisplay, xge_specialwin[0].thewindow,
                               xge_MAX_SPECIAL_WIDTH, xge_MAX_SPECIAL_HEIGHT, 1 );
    xge_specialpixmapgc = XCreateGC ( xgedisplay, xge_specialpixmap,
                                      GCForeground | GCLineWidth | GCLineStyle,
                                      &gcv );
    if ( xge_specialwingc == None || xge_specialpixmapgc == None ||
         xge_specialpixmapgc == None )
      goto failure;
    xge_specialwin_in_use = false;
  }
  return true;

failure:
  if ( xge_specialwingc != None ) {
    XFreeGC ( xgedisplay, xge_specialwingc );
    xge_specialwingc = None;
  }
  if ( xge_specialpixmap != None ) {
    XFreePixmap ( xgedisplay, xge_specialpixmap );
    xge_specialpixmap = None;
  }
  if ( xge_specialpixmapgc != None ) {
    XFreeGC ( xgedisplay, xge_specialpixmapgc );
    xge_specialpixmapgc = None;
  }
  for ( i--; i >= 0; i-- ) {
    XDestroyWindow ( xgedisplay, xge_specialwin[i].thewindow );
    xge_specialwin[i].thewindow = None;
  }
  return false;
} /*_xge_CreateSpecialWin*/

void _xge_UnmapSpecialWin ( void )
{
  int i;

  if ( xge_specialwin_in_use ) {
    for ( i = 0; i < 4; i++ )
      if ( xge_specialwin[i].mapped ) {
        XUnmapWindow ( xgedisplay, xge_specialwin[i].thewindow );
        xge_specialwin[i].mapped = false;
      }
    _xge_compspecialwinsizes = NULL;
    _xge_drawspecialwin = NULL;
    _xge_maskspecialwin = NULL;
    _xge_specialwdg = NULL;
    xge_specialwin_in_use = false;
  }
} /*_xge_UnmapSpecialWin*/

boolean _xge_MapSpecialWin ( int trwin,
                             void (*compspecialwinsizes)(xge_widget *wdg),
                             void (*maskspecialwin)(xge_widget *wdg),
                             void (*drawspecialwin)(int w, xge_widget *wdg),
                             xge_widget *wdg )
{
  int i;

  if ( xge_use_specialwin ) {
    _xge_UnmapSpecialWin ();
    _xge_compspecialwinsizes = compspecialwinsizes;
    _xge_maskspecialwin  = maskspecialwin;
    _xge_drawspecialwin  = drawspecialwin;
    _xge_specialwdg = wdg;
    for ( i = 0; i < 4; i++ )
      if ( xge_specialwin[i].nonempty ) {
        XSetTransientForHint ( xgedisplay, xge_specialwin[i].thewindow,
                               xge_windesc[trwin].thewindow );
        XMoveResizeWindow ( xgedisplay, xge_specialwin[i].thewindow,
                            xge_specialwin[i].xpos, xge_specialwin[i].ypos,
                            xge_specialwin[i].thewinrect.width,
                            xge_specialwin[i].thewinrect.height );
        XMapRaised ( xgedisplay, xge_specialwin[i].thewindow );
        xge_specialwin[i].mapped = true;
    }
    XSync ( xgedisplay, False );
    xge_specialwin_in_use = true;
    return true;
  }
  else
    return false;
} /*_xge_MapSpecialWin*/

void _xge_DestroySpecialWin ( void )
{
  int i;

  _xge_UnmapSpecialWin ();
  for ( i = 0; i < 4; i++ ) {
    XDestroyWindow ( xgedisplay, xge_specialwin[i].thewindow );
    xge_specialwin[i].thewindow = None;
  }
  if ( xge_specialwingc != None ) {
    XFreeGC ( xgedisplay, xge_specialwingc );
    xge_specialwingc = None;
  }
  if ( xge_specialpixmapgc != None ) {
    XFreeGC ( xgedisplay, xge_specialpixmapgc );
    xge_specialpixmapgc = None;
  }
  if ( xge_specialpixmap != None ) {
    XFreePixmap ( xgedisplay, xge_specialpixmap );
    xge_specialpixmap = None;
  }
  xge_try_ext = true;
} /*_xge_DestroySpecialWin*/

void _xge_RemaskSpecialWin ( void )
{
  int i;
        /* clear the mask */
  XSetForeground ( xgedisplay, xge_specialpixmapgc, 0 );
  XFillRectangle ( xgedisplay, xge_specialpixmap, xge_specialpixmapgc, 0, 0,
                   xge_MAX_SPECIAL_WIDTH, xge_MAX_SPECIAL_HEIGHT );
        /* create the new mask */
  if ( _xge_maskspecialwin )
    _xge_maskspecialwin ( _xge_specialwdg );
  for ( i = 0; i < 4; i++ ) {
    XShapeCombineMask ( xgedisplay, xge_specialwin[i].thewindow, ShapeClip,
                        xge_specialwin[i].xpix, xge_specialwin[i].ypix,
                        xge_specialpixmap, ShapeSet );
    XShapeCombineMask ( xgedisplay, xge_specialwin[i].thewindow, ShapeBounding,
                        xge_specialwin[i].xpix, xge_specialwin[i].ypix,
                        xge_specialpixmap, ShapeSet );
  }
} /*_xge_RemaskSpecialWin*/

void _xge_SpecialWinEvent ( void )
{
  int w;

  if ( !xge_specialwin_in_use )
    return;

  switch ( xgeevent.type ) {
case Expose:
    if ( _xge_drawspecialwin ) {
      for ( w = 0; w < 4; w++ )
        if ( xge_specialwin[w].thewindow == xgeevent.xany.window ) {
          _xge_drawspecialwin ( w, _xge_specialwdg );
          break;
        }
    }
    XSync ( xgedisplay, False );
    return;

default:
    return;
  }
} /*_xge_SpecialWinEvent*/

void _xge_OutSpecialMaskPixels ( const xpoint *buf, int n )
{
  XDrawPoints ( xgedisplay, xge_specialpixmap, xge_specialpixmapgc,
                (XPoint*)buf, n, CoordModeOrigin );
} /*_xge_OutSpecialPixels*/
#endif

