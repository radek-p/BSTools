
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camerad.h"
#include "multibs.h"
#include "xgedit.h"
#include "xgeipc.h"

#include "render.h"
#include "spl3d.h"
#include "ed3ds.h"
#include "ed3dswidgets.h"
#include "pomnijipc.h"

int ipc_state = ipcSTATE_NOTHING;
boolean constraints_sent = false;

/* ///////////////////////////////////////////////////////////////////////// */
void BreakSprBind ( void )
{
  if ( bind_spr ) {
    bind_spr = false;
    xge_SetWindow ( win1 );
    xge_Redraw ();
  }
} /*BreakSprBind*/

void BreakNLBlending ( void )
{
  if ( ipc_state != ipcSTATE_NOTHING ) {
    xge_SignalTheChild (); 
    xge_ParentFlushPipe ();
    ipc_state = ipcSTATE_NOTHING;
    bl_optbtn->data0 = txtOptimize;
  }
} /*BreakNLBlending*/

void ProcessIdleCommand ( int key, short x, short y )
{
  switch ( key ) {
case IDLE_COMMAND_START_RENDERING:
    StartRendering ();
    break;

case IDLE_COMMAND_RENDER:
    ContRendering ();
    break;

case IDLE_COMMAND_SEND_OPTIONS:
    constraints_sent = false;
    SendOptionsToChild ();
    ipc_state = ipcSTATE_OPTIONS_SENT;
    break;

case IDLE_COMMAND_SEND_PATCH:
    SendPatchToChild ();
    ipc_state = ipcSTATE_PATCH_SENT;
    break;

case IDLE_COMMAND_SEND_CONSTRAINTS:
    SendConstraintsToChild ();
    ipc_state = ipcSTATE_CONSTRAINTS_SENT;
    constraints_sent = true;
    break;

case IDLE_COMMAND_RECEIVE_PATCH:
/*
printf ( "parent: getting the patch\n" );
*/
    GetPatchFromChild ();
    SetKWindNKnots ();
    ResizeObject ();
    switch ( x ) {
  case 0:
/*
printf ( "parent: continue\n" );
*/
      xge_CallTheChild ( cmdCONTINUE_OPT, 0 );
      break;
  case 1:
/*
printf ( "parent: finished\n" );
*/
      ipc_state = ipcSTATE_NOTHING;
      bl_optbtn->data0 = txtOptimize;
      break;
    }
    swind_picture = false;
    sw_nonlin_blending = true;
    sw_bind_blending = true;
    sw_triharmonic_blending = false;
    if ( RenderingIsOn )
      BreakRendering ( false );
    xge_RedrawAll ();
    break;

default:
    break;
  }
} /*ProcessIdleCommand*/

void ProcessChildMessage ( int msg, int size )
{
/*
printf ( "child msg %d, size %d, state = %d,\n", msg, size, ipc_state );
*/
  switch ( msg ) {
case cmdCHILD_READY:
    child_ready = true;
    break;

case cmdGOT_INTERRUPT:
    ipc_state = ipcSTATE_NOTHING;
    break;

case cmdSUCCESS:
/*
printf ( "parent: %d\n", ipc_state );
*/
    switch ( ipc_state ) {
  case ipcSTATE_NOTHING:
      break;
  case ipcSTATE_OPTIONS_SENT:
      xge_CallTheChild ( cmdGET_BSPATCH, PatchDataSize () );
      xge_PostIdleCommand ( IDLE_COMMAND_SEND_PATCH, 0, 0 );
      break;
  case ipcSTATE_PATCH_SENT:
      if ( sw_blending_constraints &&
           n_blending_constraints > 0 ) {
        xge_CallTheChild ( cmdGET_CONSTRAINTS, ConstraintDataSize () );
        xge_PostIdleCommand ( IDLE_COMMAND_SEND_CONSTRAINTS, 0, 0 );
      }
      else {
        ipc_state = ipcSTATE_G2BLOPT_LAUNCHED;
        xge_CallTheChild ( cmdG2BL_OPTIMIZE_LMT, 0 );
        bl_optbtn->data0 = txtInterrupt;
        xge_SetWindow ( win1 );
        xge_Redraw ();
      }
      break;
  case ipcSTATE_CONSTRAINTS_SENT:
      ipc_state = ipcSTATE_G2BLOPT_LAUNCHED;
      xge_CallTheChild ( cmdG2BL_OPTIMIZE_LMT, 0 );
      bl_optbtn->data0 = txtInterrupt;
      xge_SetWindow ( win1 );
      xge_Redraw ();
      break;
  case ipcSTATE_G2BLOPT_LAUNCHED:
                       /* the computation is complete, ask for the result */
/*
printf ( "parent: optimization finished\n" );
*/
      ipc_state = ipcSTATE_G2BLOPT_FINISHED;
      xge_CallTheChild ( cmdSEND_BSPATCH, 0 );
      xge_PostIdleCommand ( IDLE_COMMAND_RECEIVE_PATCH, 1, 0 );
      break;
  default:
      break;
    }
    break;

case cmdERROR:
    switch ( ipc_state ) {
  case ipcSTATE_G2BLOPT_LAUNCHED:
      break;
  default:
      break;
    }
    ipc_state = ipcSTATE_NOTHING;
    sw_nonlin_blending = false;
    bl_optbtn->data0 = txtOptimize;
    xge_SetWindow ( win1 );
    xge_Redraw ();
    break;

case cmdG2BL_PARTIAL_RESULT:
/*
printf ( "parent: partial\n" );
*/
    xge_CallTheChild ( cmdSEND_BSPATCH, 0 );
    xge_PostIdleCommand ( IDLE_COMMAND_RECEIVE_PATCH, 0, 0 );
    break;

case cmdG2BL_FINAL_RESULT:
    break;

default:
    break;
  }
/*
printf ( "%d\n", ipc_state );
*/
} /*ProcessChildMessage*/

int ProcessOtherMsg ( int msg, int key, short x, short y )
{
  int win;

  switch ( msg ) {
case xgemsg_KEY:
    switch ( key ) {
  case 'M':
      xgeResizeWindow ( xge_MAX_WIDTH, xge_MAX_HEIGHT );
      break;
  case 'm':
      xgeResizeWindow ( xge_WIDTH, xge_HEIGHT );
      break;
  case 'D': case 'd':
      xgeMoveWindow ( 0, 512 );
      break;
  case 'S': case 's':
      DumpData ();
      break;
/*
  case 'Q': case 'q':
      xge_done = 1;
      break;
*/
  default:
      break;
    }
    break;

case xgemsg_RESIZE:
    win = xge_CurrentWindow ();
    ResizeWinStatus ( win );
    if ( win == win0 ) {
      swind_picture = false;
      if ( RenderingIsOn )
        BreakRendering ( true );
    }
    xge_Redraw ();
    break;

case xgemsg_IDLE_COMMAND:
    ProcessIdleCommand ( key, x, y );
    break;

case xgemsg_CHILD_MESSAGE:
    ProcessChildMessage ( xgeevent.xclient.data.l[0],
                          xgeevent.xclient.data.l[1] );
    break;

case xgemsg_CHILD_FAILURE:
    if ( !MakeTheChildProcess () )
      xge_DisplayErrorMessage ( ErrorChildProcessTerminated, -1 );
    break;

default:
    break;
  }
  return 1;
} /*ProcessOtherMsg*/

void OpenFile ( void )
{
  int i;

  xge_RemovePopup ( true );
  if ( xge_GetCurrentListBoxString ( &filelist, filename ) ) {
    xge_ClearListBox ( &filelist );
    xge_ClearListBox ( &dirlist );
    ResetObject ();
    if ( ReadBSPatch ( filename ) ) {
      G2CheckIfClamped ();
      ResizeObject ();
      xge_T2KnotWindFindMapping ( &kwind );
/*      ClearPointMarking (
          (lastknot_u-degree_u)*(lastknot_v-degree_v),
          mkpoints ); */
    }
    else {
      ResetObject ();
      SetKWindNKnots ();
      xge_DisplayErrorMessage ( ErrorMsgCannotOpen, -1 );
    }
  }
  UpdateBlendingRangeWidgets ();
  FindBoundingBox ( &swind.RefBBox );
  xge_3DwindSetupParProj ( &swind, &swind.RefBBox );
  swind.PerspBBox = swind.RefBBox;
  xge_3DwindSetupPerspProj ( &swind, false );
  for ( i = 0; i < 4; i++ )
    ProjectSurface ( i );
  bind_spr = sw_bind_blending = blending_mat_valid =
  sw_triharmonic_blending = sw_nonlin_blending = false;
  if ( RenderingIsOn )
    BreakRendering ( true );
  swind_picture = false;
  xge_RedrawAll ();
} /*OpenFile*/

