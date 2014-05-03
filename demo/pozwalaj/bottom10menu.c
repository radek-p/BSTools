
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


xge_widget *InitBottom10Menu ( xge_widget *prev )
{
  xge_widget *menu;

  win1statl = xge_NewTextWidget ( win1, NULL, 0,
                                  xge_WIDTH-SIDEMENUWIDTH0, 16, 0, 0, status1 );
  xge_SetWidgetPositioning ( win1statl, 0, 0, 0 );
  win1cmdl = xge_NewStringEd ( win1, win1statl, textedM12COMMAND,
                               xge_WIDTH-SIDEMENUWIDTH0, 19, 0, 0,
                               MAX_COMMAND_LGT, command1, &command1_editor );
  xge_SetWidgetPositioning ( win1cmdl, 0, 0, BOTTOMMENUHEIGHT1-19 );

        /* initially this menu is empty, one of the above widgets is */
        /* placed in it when needed */
  menu = xge_NewMenu ( win1, prev, BOTTOMMENU1, xge_WIDTH-SIDEMENUWIDTH0,
                       0, SIDEMENUWIDTH0, xge_HEIGHT, win1cmdl );
  return menu;
} /*InitBottom10Menu*/

void ExecuteCommand1 ( void )
{
  if ( command1[0] ) {
printf ( "Command: %s\n", command1 );

    command1[0] = 0;
    command1_editor.start = command1_editor.pos = 0;
  }
} /*ExecuteCommand1*/

int Bottom10MenuCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_TEXT_EDIT_VERIFY:
    switch ( er->id ) {
  case textedM12COMMAND:
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_ENTER:
    switch ( er->id ) {
  case textedM12COMMAND:
      ExecuteCommand1 ();
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_ESCAPE:
    switch ( er->id ) {
  case textedM12COMMAND:
      return 1;
  default:
      return 0;
    }

default:
    return 1;
  }
} /*Bottom10MenuCallBack*/

