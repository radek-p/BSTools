
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


xge_widget *InitPopup00 ( void )
{
  xge_widget *w, *menu;

  w = xge_NewButton ( win0, NULL, btnP00OPEN, 58, 19, 2, 22, txtOpen );
  w = xge_NewButton ( win0, w, btnP00SAVE, 58, 19, 2, 42, txtSave );
  w = xge_NewButton ( win0, w, btnP00SAVEAS, 58, 19, 2, 62, txtSaveAs );
  w = xge_NewButton ( win0, w, btnP00EXIT, 58, 19, 2, 82, txtExit );
  menu = xge_NewFMenu ( win0, NULL, POPUP00, 62, 83, 0, 20, w );
  menu->msgproc = xge_PopupMenuMsg;
  return menu;
} /*InitPopup00*/

void SetCurrentDir ( void )
{
  int l;

  getcwd ( current_directory, MAX_PATH_LGT+1 );
  l = strlen ( current_directory );
  if ( l <= MAX_PATH_SHRT )
    strcpy ( current_dir, current_directory );
  else {
    strcpy ( current_dir, "..." );
    strcpy ( &current_dir[3], &current_directory[l-MAX_PATH_SHRT+3] );
  }
} /*SetCurrentDir*/

void PreparePopup01 ( void )
{
  SetCurrentDir ();
  xge_SetupFileList ( &filelist1, ".", file_filter, showhiddenfiles );
  xge_SetupDirList ( &dirlist1, ".", NULL, showhiddenfiles, NULL );
} /*PreparePopup01*/

void PreparePopup02 ( void )
{
  SetCurrentDir ();
  xge_SetupFileList ( &filelist2, ".", file_filter, showhiddenfiles );
  xge_SetupDirList ( &dirlist2, ".", NULL, showhiddenfiles, NULL );
} /*PreparePopup02*/

boolean FilenameCorrect ( char *fn )
{
  return fn[0] != 0;
} /*FilenameCorrect*/

int Popup00CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnP00OPEN:
      xge_RemovePopup ( true );
      PreparePopup01 ();
      OpenPopup ( popup01, true );
      return 1;
  case btnP00SAVE:
      if ( FilenameCorrect ( filename ) ) {
        xge_RemovePopup ( true );
        Popup02SaveFile ();
      }
      else
        goto save_as;
      return 1;
  case btnP00SAVEAS:
save_as:
      xge_RemovePopup ( true );
      PreparePopup02 ();
      OpenPopup ( popup02, true );
      return 1;
  case btnP00EXIT:
      xge_RemovePopup ( true );
      OpenPopup ( popup03, true );
      return 1;
  default:
      return 0;
    }

default:
    return 0;
  }
} /*Popup00CallBack*/

