
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camera.h"
#include "multibs.h"
#include "bsmesh.h"
#include "g2blendingd.h"
#include "egholed.h"
#include "xgedit.h"
#include "xgledit.h"

#include "widgets.h"
#include "editor.h"
#include "editor_bezc.h"
#include "editor_bsc.h"
#include "editor_bezp.h"
#include "editor_bsp.h"
#include "editor_bsm.h"
#include "editor_bsh.h"
#include "pozwalaj.h"
#include "calc.h"


xge_widget *InitBottom00Menu ( xge_widget *prev )
{
  xge_widget *menu;

  win0statl = xge_NewTextWidget ( win0, NULL, 0, xge_WIDTH-SIDEMENUWIDTH, 16, 0, 0,
                                  status0 );
  xge_SetWidgetPositioning ( win0statl, 0, 0, 0 );
  win0cmdl = xge_NewStringEd ( win0, win0statl, textedM02COMMAND,
                               xge_WIDTH-SIDEMENUWIDTH, 19, 0, 0,
                               MAX_COMMAND_LGT, command0, &command0_editor );
  xge_SetWidgetPositioning ( win0cmdl, 0, 0, BOTTOMMENUHEIGHT1-19 );

        /* initially this menu is empty, one of the above widgets is */
        /* placed in it when needed */
  menu = xge_NewMenu ( win0, prev, BOTTOMMENU0, xge_WIDTH-SIDEMENUWIDTH,
                       0, SIDEMENUWIDTH, xge_HEIGHT, win0cmdl );
  return menu;
} /*InitBottom00Menu*/

void ExecuteCommand0 ( void )
{
  double  value;
  boolean expr_ok;
  char    *p0, *p1, c;

/* at the moment the only command is to evaluate a simple expression, */
/* which is planned to change some time in future. */
  if ( command0[0] ) {
    expr_ok = EvaluateExpression ( command0, &value );
    if ( expr_ok ) {
      sprintf ( status0, "%20.15g", value );
        /* delete leading spaces */
      if ( status0[0] == ' ' ) {
        p0 = p1 = status0;
        while ( *p1 == ' ' )
          p1 ++;
        do {
          c = *p0 ++ = *p1 ++;
        } while ( c );
      }
      XStoreBytes ( xgedisplay, status0, strlen(status0)+1 );
    }
    else
      strcpy ( status0, "Error" );
    if ( win0commandline ) {
      xge_SetClipping ( win0statl );
      win0statl->redraw ( win0statl, true );
    }
    command0[0] = 0;
    command0_editor.start = command0_editor.pos = 0;
  }
} /*ExecuteCommand0*/

int Bottom00MenuCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_TEXT_EDIT_VERIFY:
    switch ( er->id ) {
  case textedM02COMMAND:
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_ENTER:
    switch ( er->id ) {
  case textedM02COMMAND:
      ExecuteCommand0 ();
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_ESCAPE:
    switch ( er->id ) {
  case textedM02COMMAND:
      return 1;
  default:
      return 0;
    }

default:
    return 0;
  }
} /*Bottom00MenuCallBack*/

