
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <sys/types.h> 
#include <signal.h>
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
#include "pkvthreads.h"
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

#define PARENT_SIDE
#include "pozwalajipc.h"
#undef PARENT_SIDE

/* this prototype seems to be missing from the my system header files */
int kill ( pid_t pid, int sig );

void RedrawGeom00Win ( void )
{
  int win;

  win = xge_CurrentWindow ();
  if ( win != win0 )
    xge_SetWindow ( win0 );
  xge_SetClipping ( geom00menu );
  geom00menu->redraw ( geom00menu, true );
  xge_ResetClipping ();
  if ( win != win0 )
    xge_SetWindow ( win );
} /*RedrawGeom00Win*/

void OpenPopup ( xge_widget *er, boolean allwin )
{
  xge_AddPopup ( er );
  xge_GrabFocus ( er, allwin );
  if ( allwin )
    xge_SetOtherWindowsCursor ( xgeCURSOR_CIRCLE );
} /*OpenPopup*/

void InitWindow0Widgets ( void )
{
  win0 = xge_CurrentWindow ();
  top00menu = InitTop00Menu ( NULL );
  side00menu = InitSide00Menu ( top00menu );
  bottom00menu = InitBottom00Menu ( side00menu );
  geom00menu = InitGeom00Menu ( bottom00menu );
  popup00 = InitPopup00 ();
  popup01 = InitPopup01 ();
  popup02 = InitPopup02 ();
  popup03 = InitPopup03 ();
  xge_SetWinEdRect ( geom00menu );
} /*InitWindow0Widgets*/

void InitWindow1Widgets ( void )
{
  win1 = xge_NewWindow ( "" );
  xge_SetWindow ( win1 );
  top10menu = InitTop10Menu ( NULL );
  side10menu = InitSide10Menu ( top10menu );
  bottom10menu = InitBottom10Menu ( side10menu );
  geom10menu = InitGeom10Menu ( bottom10menu );
  popup10 = InitPopup10 ();
  popup11 = InitPopup11 ();
  popup12 = InitPopup12 ();
  popup13 = InitPopup13 ();
  popup14 = InitPopup14 ();
  popup15 = InitPopup15 ();
  xge_SetWinEdRect ( geom10menu );
} /*InitWindow1Widgets*/

boolean ProcessCMDLineParameters ( int argc, char *argv[] )
{
  boolean result;
  int     i, j, l;
  char    *sstr, *s;

/* this program may be invoked with parameters - file names, with */
/* the suffix file_ext (.bs), containing data to be read in. */
  result = true;
  for ( i = 1; i < argc; i++ ) {
    sstr = strstr ( argv[i], file_ext );
    if ( sstr && !strcmp ( file_ext, sstr ) ) {
        /* remove the path, leaving only the file name */
      l = strlen ( argv[i] );
      s = argv[i] + l;
      for ( j = l; j > 0; j--, s-- ) {
        if ( *s == '/' ) {
          s++;
          break;
        }
      }
      strncpy ( filename, s, MAX_FILENAME_LGT+1 );
      if ( !GeomObjectReadFile ( filename, Popup01CameraReader ) ) {
        result = false;
        break;
      }
    }
  }
  current_go = first_go;
  if ( current_go )
    SetupObjectSpecificMenus ( current_go );
  return result;
} /*ProcessCMDLineParameters*/

void init_program ( int argc, char *argv[] )
{
  boolean files_ok;

  setvbuf ( stdout, NULL, _IONBF, 0 );
  if ( !pkv_InitScratchMem ( SCRATCHMEMSIZE ) ) {
    printf ( "Error: cannot allocate scratch memory stack\n" );
    exit ( 1 );
  }
  getcwd ( initial_directory, MAX_PATH_LGT+1 );
        /* pthreads so far are used by the renderer, */
        /* and also by pozwalaj_proc, which creates them for itself */
  ncpu = pkv_FindNCPU ();
  if ( ncpu > 1 ) {
    if ( !pkv_InitPThreads ( MAX_PTHREADS ) )
      ncpu = 1;
  }
  GeomObjectInitList ();
  InitWindow0Widgets ();
  InitWindow1Widgets ();
  RendInit ();
  files_ok = ProcessCMDLineParameters ( argc, argv );
  xge_RedrawAll ();
  xge_SetWindow ( win0 );
  if ( files_ok )
    xge_DisplayInfoMessage ( InfoMsg, -1 );
  else
    xge_DisplayErrorMessage ( ErrorMsgCannotOpen, -1 );
} /*init_program*/

void destroy_program ( void )
{
  RendDestroy ();
  printf ( "Scratch memory used: %d out of %d bytes\n",
           (int)pkv_MaxScratchTaken(), SCRATCHMEMSIZE );
  pkv_DestroyScratchMem ();
  if ( xge_ChildIsActive () ) {
    signal ( SIGCHLD, SIG_DFL );
    kill ( xge_child_pid, SIGKILL );
  }
} /*destroy_program*/

void ResizeWindow0 ( short x, short y )
{
  short w1, y1, h0, h1, smw;

  h1 = 0;
  if ( win0statusline )
    h1 = BOTTOMMENUHEIGHT0;
  if ( win0commandline )
    h1 += BOTTOMMENUHEIGHT1;
  h0 = y - TOPMENUHEIGHT - h1;
  y1 = TOPMENUHEIGHT + h0;
  ChangeSide00MenuWidth ( y-TOPMENUHEIGHT );
  smw = side00menu_wide ? SIDEMENUWIDTH1 : SIDEMENUWIDTH0;
  w1 = x-smw;
  top00menu->msgproc ( top00menu, xgemsg_RESIZE, 0, x, TOPMENUHEIGHT );
  side00menu->msgproc ( side00menu, xgemsg_RESIZE, 0, smw, y-TOPMENUHEIGHT );
  if ( side00menu->data1 == side00ewidgets )
    Side00eMenuCallBack ( NULL, xgemsg_RESIZE, 0, smw, y-TOPMENUHEIGHT );
  bottom00menu->y = y1;
  if ( win0statusline ) {
    win0statl->w = w1;
    win0statl->y = y - h1;
    win0statl->yofs = 0;
    win0cmdl->y = y - h1 + BOTTOMMENUHEIGHT0;
    win0cmdl->yofs = BOTTOMMENUHEIGHT0 + BOTTOMMENUHEIGHT1 - win0cmdl->h;
  }
  else {
    win0statl->y = y;
    win0statl->yofs = BOTTOMMENUHEIGHT1;
  }
  if ( win0commandline ) {
    win0cmdl->w = w1;
    command0_editor.chdisp = (short)((w1-4) / xge_CHAR_WIDTH);
    if ( !win0statusline ) {
      win0cmdl->y = y - win0cmdl->h;
      win0cmdl->yofs = BOTTOMMENUHEIGHT1 - win0cmdl->h;
    }
  }
  bottom00menu->x = geom00menu->x = geom00win->x = smw;
  bottom00menu->msgproc ( bottom00menu, xgemsg_RESIZE, 0, w1, h1 );
  geom00menu->msgproc ( geom00menu, xgemsg_RESIZE, 0, w1, h0 );
  geom00win->msgproc ( geom00win, xgemsg_RESIZE, 0, w1, h0 );
  xge_Redraw ();
} /*ResizeWindow0*/

int Window0CallBack ( int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_RESIZE:
    ResizeWindow0 ( x, y );
    return 1;

case xgemsg_KEY:
    switch ( key ) {
  case 'm':
      if ( xge_current_width != xge_WIDTH ||
           xge_current_height != xge_HEIGHT )
        xgeResizeWindow ( xge_WIDTH, xge_HEIGHT );
      else
        xgeResizeWindow ( 800, 600 );
      return 1;
  case 'M':
      if ( xge_current_width != xge_MAX_WIDTH ||
           xge_current_height != xge_MAX_HEIGHT )
        xgeResizeWindow ( xge_MAX_WIDTH, xge_MAX_HEIGHT );
      else
        xgeResizeWindow ( 1024, 768 );
      return 1;
  default:
      return 0;
    }

case xgemsg_IDLE_COMMAND:
    return ProcessIdleCommand ( key, x, y );

case xgemsg_CHILD_MESSAGE:
    ProcessChildMessage ( xgeevent.xclient.data.l[0],
                          xgeevent.xclient.data.l[1] );
    return 1;

default:
    return 0;
  }
} /*Window0CallBack*/

void ResizeWindow1 ( short x, short y )
{
  short w1, y1, h0, h1, smw;

  h1 = 0;
  if ( win1statusline )
    h1 = BOTTOMMENUHEIGHT0;
  if ( win1commandline )
    h1 += BOTTOMMENUHEIGHT1;
  h0 = y - TOPMENUHEIGHT - h1;
  y1 = TOPMENUHEIGHT + h0;
  ChangeSide10MenuWidth ( y-TOPMENUHEIGHT );
  smw = side10menu_wide ? SIDEMENUWIDTH1 : SIDEMENUWIDTH0;
  w1 = x-smw;
  top10menu->msgproc ( top10menu, xgemsg_RESIZE, 0, x, TOPMENUHEIGHT );
  side10menu->msgproc ( side10menu, xgemsg_RESIZE, 0, smw, y-TOPMENUHEIGHT );
  if ( side10menu->data1 == side10wdg_bsm_edit ||
       side10menu->data1 == side10wdg_bsm_opt3 )
    Side10MenuBsmCallBack ( side10menu, xgemsg_RESIZE, 0,
                            smw, y-TOPMENUHEIGHT );
  bottom10menu->y = y1;
  if ( win1statusline ) {
    win1statl->w = w1;
    win1statl->y = y - h1;
    win1statl->yofs = 0;
    win1cmdl->y = y - h1 + BOTTOMMENUHEIGHT0;
    win1cmdl->yofs = BOTTOMMENUHEIGHT0 + BOTTOMMENUHEIGHT1 - win1cmdl->h;
  }
  else {
    win1statl->y = y;
    win1statl->yofs = BOTTOMMENUHEIGHT1;
  }
  if ( win1commandline ) {
    win1cmdl->w = w1;
    command1_editor.chdisp = (short)((w1-4) / xge_CHAR_WIDTH);
    if ( !win1statusline ) {
      win1cmdl->y = y - win1cmdl->h;
      win1cmdl->yofs = BOTTOMMENUHEIGHT1 - win1cmdl->h;
    }
  }
  bottom10menu->x = geom10menu->x = smw;
  bottom10menu->msgproc ( bottom10menu, xgemsg_RESIZE, 0, x-smw, h1 );
  geom10menu->msgproc ( geom10menu, xgemsg_RESIZE, 0, x-smw, h0 );
  Geom10MenuResize ();
  xge_Redraw ();
} /*ResizeWindow1*/

int Window1CallBack ( int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_RESIZE:
    ResizeWindow1 ( x, y );
    return 1;

case xgemsg_KEY:
    switch ( key ) {
  case 'm':
      if ( xge_current_width != xge_WIDTH ||
           xge_current_height != xge_HEIGHT )
        xgeResizeWindow ( xge_WIDTH, xge_HEIGHT );
      else
        xgeResizeWindow ( 800, 600 );
      return 1;
  case 'M':
      if ( xge_current_width != xge_MAX_WIDTH ||
           xge_current_height != xge_MAX_HEIGHT )
        xgeResizeWindow ( xge_MAX_WIDTH, xge_MAX_HEIGHT );
      else
        xgeResizeWindow ( 1024, 768 );
      return 1;
  default:
      return 0;
    }

case xgemsg_IDLE_COMMAND:
    return ProcessIdleCommand ( key, x, y );

case xgemsg_CHILD_MESSAGE:
    ProcessChildMessage ( xgeevent.xclient.data.l[0],
                          xgeevent.xclient.data.l[1] );
    return 1;

default:
    return 0;
  }
} /*Window1CallBack*/

void HandleChildFailure ( void )
{
  xge_SetWindow ( win1 );
  if ( ipc_state != ipcstate_NO_CHILD ) {
    ipc_state = ipcstate_NO_CHILD;
    bsp_bl_optimizeG1->data0 = bsp_bl_optimizeG2->data0 =
    bsm_bl_optimize->data0 = txtOptimize;
    xge_DisplayErrorMessage ( ErrorMsgChildProcessFailure, 0 );
  }
} /*HandleChildFailure*/

int NonWidgetCallBack ( int msg, int key, short x, short y )
{
  int win;

  switch ( msg ) {
case xgemsg_POPUP_REMOVED:
    switch ( key ) {
  case POPUP12:
      CleanupPopup12 ();
      break;
  case POPUP14:
      CleanupPopup14 ();
      break;
  default:
      break;
    }
    return 1;
case xgemsg_CHILD_FAILURE:
    HandleChildFailure ();
    return 1;
default:
    win = xge_CurrentWindow ();
    if ( win == win0 )
      return Window0CallBack ( msg, key, x, y );
    else
      return Window1CallBack ( msg, key, x, y );
  }
} /*NonWidgetCallBack*/

int CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  if ( er ) {
    switch ( er->id & MENUMASK ) {
  case TOPMENU0:
      return Top00MenuCallBack ( er, msg, key, x, y );
  case SIDEMENU0:
      return Side00MenuCallBack ( er, msg, key, x, y );
  case SIDEMENU0e:
      return Side00eMenuCallBack ( er, msg, key, x, y );
  case SIDEMENU0b:
      return Side00bMenuCallBack ( er, msg, key, x, y );
  case BOTTOMMENU0:
      return Bottom00MenuCallBack ( er, msg, key, x, y );
  case GEOMMENU0:
      return Geom00MenuCallBack ( er, msg, key, x, y );
  case POPUP00:
      return Popup00CallBack ( er, msg, key, x, y );
  case POPUP01:
      return Popup01CallBack ( er, msg, key, x, y );
  case POPUP02:
      return Popup02CallBack ( er, msg, key, x, y );
  case POPUP03:
      return Popup03CallBack ( er, msg, key, x, y );
  case GEOMWIN2D0:
      return Geom00Win2DCallBack ( er, msg, key, x, y );
  case GEOMWIN3D0:
      return Geom00Win3DCallBack ( er, msg, key, x, y );
  case TOPMENU1:
      return Top10MenuCallBack ( er, msg, key, x, y );
  case SIDEMENU1:
      return Side10MenuCallBack ( er, msg, key, x, y );
  case SIDEMENU1_BEZC:
      return Side10MenuBezcCallBack ( er, msg, key, x, y );
  case SIDEMENU1_BEZP:
      return Side10MenuBezpCallBack ( er, msg, key, x, y );
  case SIDEMENU1_BSC:
      return Side10MenuBscCallBack ( er, msg, key, x, y );
  case SIDEMENU1_BSP:
      return Side10MenuBspCallBack ( er, msg, key, x, y );
  case SIDEMENU1_BSM:
      return Side10MenuBsmCallBack ( er, msg, key, x, y );
  case SIDEMENU1_BSH:
      return Side10MenuBshCallBack ( er, msg, key, x, y );
  case BOTTOMMENU1:
      return Bottom10MenuCallBack ( er, msg, key, x, y );
  case GEOMMENU1:
      return Geom10MenuCallBack ( er, msg, key, x, y );
  case POPUP10:
      return Popup10CallBack ( er, msg, key, x, y );
  case POPUP11:
      return Popup11CallBack ( er, msg, key, x, y );
  case POPUP12:
      return Popup12CallBack ( er, msg, key, x, y );
  case POPUP13:
      return Popup13CallBack ( er, msg, key, x, y );
  case POPUP14:
      return Popup14CallBack ( er, msg, key, x, y );
  case POPUP15:
      return Popup15CallBack ( er, msg, key, x, y );
  default:
      return 0;
    }
  }
  else
    return NonWidgetCallBack ( msg, key, x, y );
} /*CallBack*/

