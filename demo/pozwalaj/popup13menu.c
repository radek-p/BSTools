
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
#include "render.h"


xge_widget *InitPopup13 ( void )
{
  xge_widget *w, *menu;

  w = xge_NewSlidebardRGB ( win1, NULL, slP13_RED, 109, 10, 50, 40, &colour_rgb[0] );
  w = xge_NewSlidebardRGB ( win1, w, slP13_GREEN, 109, 10, 50, 60, &colour_rgb[1] );
  w = xge_NewSlidebardRGB ( win1, w, slP13_BLUE, 109, 10, 50, 80, &colour_rgb[2] );
  w = colour_box = xge_NewRGBSampled ( win1, w, spwP13_COLOUR, 50, 50, 180, 40,
                                   &colour_rgb[0] );
  w = xge_NewButton ( win1, w, btnP13_OK, 60, 18, 110, 100, txtOK );
  menu = xge_NewFMenu ( win1, NULL, POPUP13, 220, 100, 30, 30, w );
  menu->msgproc = xge_PopupMenuMsg;
  return menu;
} /*InitPopup13*/

int Popup13CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnP13_OK:
      xge_RemovePopup ( true );
      return 1;
  default:
      return 0;
    }

case xgemsg_SLIDEBAR_COMMAND:
    switch ( er->id ) {
  case slP13_RED:
  case slP13_GREEN:
  case slP13_BLUE:
      xge_SetClipping ( colour_box );
      colour_box->redraw ( colour_box, true );
      if ( current_go )
        memcpy ( current_go->colour, colour_rgb, 3*sizeof(double) );
      return 1;
  default:
      return 0;
    }

default:
    return 0;
  }
} /*Popup13CallBack*/

