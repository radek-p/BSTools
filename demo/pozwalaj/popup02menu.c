
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


xge_widget *InitPopup02 ( void )
{
  xge_widget *w, *menu;

  w = xge_NewTextWidget ( win0, NULL, txtP02DIRSTR, 380, 16, 20+10, 40+10,
                          current_dir );
  w = xge_NewSwitch ( win0, w, swP02HIDDEN, 54, 16, 20+10, 40+30,
                      txtHidden, &showhiddenfiles );
  w = xge_NewButton ( win0, w, btnP02SAVE, 58, 19, 20+91, 40+232-30, txtSave );
  w = xge_NewButton ( win0, w, btnP02CANCEL, 58, 19, 20+251, 40+232-30, txtCancel );
  w = xge_NewListBox ( win0, w, lbP02DIRLIST, 180, 131, 20+10, 40+60, &dirlist2 );
  w = xge_NewListBox ( win0, w, lbP02FILELIST, 180, 131, 220+10, 40+60, &filelist2 );
  w = xge_NewTextWidget ( win0, w, 0, 42, 16, 220-44+10, 40+30, txtSaveAs );
  w = xge_NewStringEd ( win0, w, txtedP02FILENAME, 180, 19, 220+10, 40+30,
                        MAX_FILENAME_LGT, filename, &filename_editor );
  menu = xge_NewFMenu ( win0, NULL, POPUP02, 400, 232, 20, 40, w );
  menu->msgproc = xge_PopupMenuMsg;
  return menu;
} /*InitPopup02*/

void Popup02ChangeDir ( void )
{
  if ( !chdir ( &dirlist2.itemstr[dirlist2.itemind[dirlist2.currentitem]] ) ) {
    PreparePopup02 ();
    xge_SetClipping ( popup02 );
    popup02->redraw ( popup02, true );
  }
} /*Popup02ChangeDir*/

void Popup02ChangeDirAlt ( short x )
{
  short l1, l2;

  l1 = strlen ( current_directory );
  l2 = strlen ( current_dir );
  x += l2-l1;
  while ( x < l1 && current_directory[x] != '/' )
    x++;
  if ( x < l1 ) {
    current_directory[x] = 0;
    if ( chdir ( current_directory ) )
      current_directory[x] = '/';
    else {
      PreparePopup02 ();
      xge_SetClipping ( popup02 );
      popup02->redraw ( popup02, true );
    }
  }
} /*Popup02ChangeDirAlt*/

void Popup02SaveFile ( void )
{
        /* this procedure may be called also from a button */
        /* in popup00 and popup03 */
  if ( !GeomObjectWriteFile ( filename, false ) )
    xge_DisplayErrorMessage ( ErrorMsgCannotSave, 0 );
} /*Popup02SaveFile*/

int Popup02CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnP02SAVE:
      if ( FilenameCorrect ( filename ) ) {
        xge_RemovePopup ( true );
        Popup02SaveFile ();
      }
      else
        xge_DisplayErrorMessage ( ErrorMsgIncorrectFilename, -1 );
      return 1;
  case btnP02CANCEL:
      xge_RemovePopup ( true );
      return 1;
  default:
      return 0;
    }

case xgemsg_SWITCH_COMMAND:
    switch ( er->id ) {
  case swP02HIDDEN:
      PreparePopup02 ();
      xge_SetClipping ( popup02 );
      popup02->redraw ( popup02, true );
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_WIDGET_CLICK:
    switch ( er->id ) {
  case txtP02DIRSTR:
      Popup02ChangeDirAlt ( x );
      return 1;
  default:
      return 0;
    }

case xgemsg_LISTBOX_ITEM_PICK:
    switch ( er->id ) {
  case lbP02DIRLIST:
      Popup02ChangeDir ();
      return 1;
  case lbP02FILELIST:
      if ( xge_GetCurrentListBoxString ( &filelist2, filename ) ) {
        filename_editor.start = filename_editor.pos = 0;
        xge_SetClipping ( popup02 );
        popup02->redraw ( popup02, true );
      }
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_VERIFY:
    switch ( er->id ) {
  case txtedP02FILENAME:
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_ENTER:
    switch ( er->id ) {
  case txtedP02FILENAME:
      return 1;
  default:
      return 0;
    }

default:
    return 0;
  }
} /*Popup02CallBack*/

