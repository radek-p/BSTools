
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


void InitNameEditor ( xge_string_ed *ed, char *name )
{
  ed->start = ed->pos = 0;
  ed->er->data0 = name;
} /*InitNameEditor*/

xge_widget *InitSide10Menu ( xge_widget *prev )
{
  xge_widget *w, *menu;

        /* widgets specific for no object */
  w = xge_NewSwitch ( win1, NULL, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_none = w;
        /* widgets specific for Bezier curves */
  InitSide10Menu_Bezc ();
        /* widgets specific for Bezier patches */
  InitSide10Menu_Bezp ();
        /* widgets specific for B-spline curves */
  InitSide10Menu_BSc ();
        /* widgets specific for B-spline patches */
  InitSide10Menu_BSp ();
        /* widgets specific for B-spline meshes */
  InitSide10Menu_BSm ();
        /* widgets specific for B-spline holes */
  InitSide10Menu_BSh ();

  menu = xge_NewMenu ( win1, prev, SIDEMENU1, SIDEMENUWIDTH,
                       xge_HEIGHT-TOPMENUHEIGHT, 0, TOPMENUHEIGHT,
                       side10wdg_none );
  return menu;
} /*InitSide10Menu*/

void SetupWin1StatusLine ( void )
{
  if ( win1statusline )
    command1[0] = 0;
  ResizeWindow1 ( xge_current_width, xge_current_height );
} /*SetupWin1StatusLine*/

void SetupWin1CommandLine ( void )
{
  if ( win1commandline ) {
    command1[0] = 0;
    command1_editor.start = command1_editor.pos = 0;
    command1_editor.er->state = xgestate_NOTHING;
  }
  ResizeWindow1 ( xge_current_width, xge_current_height );
} /*SetupWin1CommandLine*/

int Side10MenuCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_SWITCH_COMMAND:
    switch ( er->id ) {
  case swM11STATUS:
      SetupWin1StatusLine ();
      return 1;
  case swM11COMMAND:
      SetupWin1CommandLine ();
      return 1;
    }

default:
    return 0;
  }
} /*Side10MenuCallBack*/

