
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012                                  */
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

#include "pkgeom.h"
#include "multibs.h"

#include "xgedit.h"
#include "xgeprivate.h"

/* ///////////////////////////////////////////////////////////////////////// */
/* In a multithreaded application, if a thread sends an event to the         */
/* (main thread of) application using SendEvent (e.g. xge_PostIdleCommand),  */
/* and then terminates, the system is waiting for some other event (e.g.     */
/* mouse motion) before delivering this thread-sent event. The purpose of    */
/* xge_FlushFromAThread is to force XWindow to flush the event queue from a  */
/* thread other than main (i.e. the one executing xge_MessageLoop).          */
static Bool FlushPredicate ( Display *display, XEvent *event, XPointer arg )
{
  return True;
} /*FlushPredicate*/

void xge_FlushFromAThread ( void )
{
  XEvent event;

  XCheckIfEvent ( xgedisplay, &event, FlushPredicate, NULL );
} /*xge_FlushFromAThread*/

