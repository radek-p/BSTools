
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
#include "bsfile.h"
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


xge_widget *InitTop10Menu ( xge_widget *prev )
{
  xge_widget *w, *menu;

  w = xge_NewButton ( win1, NULL, btnM10OBJECTS, 58, 19, 0, 0, txtObjects );
  w = xge_NewButton ( win1, w, btnM10EDIT, 58, 19, 60, 0, txtEdit );
  w = xge_NewButton ( win1, w, btnM10VIEW, 58, 19, 120, 0, txtView );
  w = xge_NewButton ( win1, w, btnM10DATA, 58, 19, 180, 0, txtData );
  w = xge_NewButton ( win1, w, btnM10OPTIONS, 58, 19, 240, 0, txtOptions );
  menu = xge_NewMenu ( win1, prev, TOPMENU1, xge_WIDTH, TOPMENUHEIGHT, 0, 0, w );
  return menu;
} /*InitTop10Menu*/

int Top10MenuCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnM10OBJECTS:
      if ( SetupObjectNameList () )
        OpenPopup ( popup12, false );
      return 1;
  case btnM10EDIT:
      SetObjectEditMenu ( current_go );
      goto resize_and_redraw;
  case btnM10VIEW:
      SetObjectViewMenu ( current_go );
      goto resize_and_redraw;
  case btnM10DATA:
      SetObjectDataMenu ( current_go );
      goto resize_and_redraw;
  case btnM10OPTIONS:
      SetObjectOptionsMenu ( current_go );
resize_and_redraw:
      if ( ChangeSide10MenuWidth ( side10menu->h ) )
        ResizeWindow1 ( xge_current_width, xge_current_height );
      else
        xge_Redraw ();
      return 1;
  default:
      return 0;
    }

default:
    return 0;
  }
} /*Top10MenuCallBack*/

