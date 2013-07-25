
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


xge_widget *InitPopup03 ( void )
{
  xge_widget *w, *menu;

  w = xge_NewTextWidget ( win0, NULL, 0, 240, 16, 110, 70, MsgReconsider );
  w = xge_NewButton ( win0, w, btnP03EXIT, 58, 19, 110, 100, txtExit );
  w = xge_NewButton ( win0, w, btnP03SAVE, 58, 19, 210, 100, txtSave );
  w = xge_NewButton ( win0, w, btnP03CANCEL, 58, 19, 310, 100, txtCancel );
  menu = xge_NewFMenu ( win0, NULL, POPUP03, 300, 70, 90, 60, w );
  menu->msgproc = xge_PopupMenuMsg;
  return menu;
} /*InitPopup03*/

int Popup03CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnP03EXIT:
      xge_done = true;
      return 1;
  case btnP03SAVE:
      if ( FilenameCorrect ( filename ) ) {
        xge_RemovePopup ( true );
        Popup02SaveFile ();
      }
      else {
        xge_RemovePopup ( true );
        PreparePopup02 ();
        OpenPopup ( popup02, true );
      }
      return 1;
  case btnP03CANCEL:
      xge_RemovePopup ( true );
      return 1;
  default:
      return 0;
    }

default: 
    return 0;
  }
} /*Popup03CallBack*/

