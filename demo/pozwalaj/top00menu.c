
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
#include "mengerc.h"
#include "bsfile.h"
#include "pkrender.h"
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


xge_widget *InitTop00Menu ( xge_widget *prev )
{
  xge_widget *w, *menu;

  w = xge_NewButton ( win0, NULL, btnM00FILE, 58, 19, 0, 0, txtFile );
  w = xge_NewButton ( win0, w, btnM00EDIT, 58, 19, 60, 0, txtEdit );
  w = xge_NewButton ( win0, w, btnM00TRANSFORM, 58, 19, 120, 0, txtTransform );
  w = xge_NewButton ( win0, w, btnM00PICTURE, 58, 19, 180, 0, txtPicture );
  w = xge_NewButton ( win0, w, btnM00ABOUT, 58, 19, xge_WIDTH-58, 0, txtAbout );
  xge_SetWidgetPositioning ( w, 1, -58, 0 );
  menu = xge_NewMenu ( win0, prev, TOPMENU0, xge_WIDTH, TOPMENUHEIGHT, 0, 0, w );
  return menu;
} /*InitTop00Menu*/

int Top00MenuCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnM00FILE:
      OpenPopup ( popup00, true );
      return 1;
  case btnM00EDIT:
      xge_SetMenuWidgets ( side00menu, side00widgets, false );
      editing_shapefunc = editing_lights = false;
      if ( ChangeSide00MenuWidth ( side00menu->h ) )
        ResizeWindow0 ( xge_current_width, xge_current_height );
      else
        xge_Redraw ();
      return 1;
  case btnM00TRANSFORM:
      xge_SetMenuWidgets ( side00menu, side00ewidgets, false );
      editing_shapefunc = editing_lights = false;
      Side00eMenuCallBack ( NULL, xgemsg_RESIZE, 0, side00menu->w, side00menu->h );
      if ( ChangeSide00MenuWidth ( side00menu->h ) )
        ResizeWindow0 ( xge_current_width, xge_current_height );
      else
        xge_Redraw ();
      return 1;
  case btnM00PICTURE:
      if ( picture_lights ) {
        xge_SetMenuWidgets ( side00menu, side00cwidgets, false );
        editing_lights = true;
        editing_shapefunc = false;
      }
      else {
        xge_SetMenuWidgets ( side00menu, side00bwidgets, false );
        editing_lights = false;
        editing_shapefunc = true;
      }
      if ( ChangeSide00MenuWidth ( side00menu->h ) )
        ResizeWindow0 ( xge_current_width, xge_current_height );
      else
        xge_Redraw ();
      return 1;
  case btnM00ABOUT:
      xge_DisplayInfoMessage ( InfoMsg, -1 );
      return 1;
  default:
      return 0;
    }

default:
    return 0;
  }
} /*Top00MenuCallBack*/

