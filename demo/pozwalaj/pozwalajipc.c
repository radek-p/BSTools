
/* ///////////////////////////////////////////////////  */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* ///////////////////////////////////////////////////  */

#include <sys/types.h>
#include <signal.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "eg2holed.h"  
#include "bsmesh.h"
#include "g2blendingd.h"
#include "bsfile.h"
#include "xgedit.h"
#include "xgeipc.h"

#include "widgets.h"
#include "editor.h"
#include "editor_bezc.h"
#include "editor_bsc.h"
#include "editor_bezp.h"
#include "editor_bsp.h"
#include "editor_bsm.h"
#include "editor_bsh.h"
#include "pozwalaj.h"
#define PARENT_SIDE
#include "pozwalajipc.h"
#undef PARENT_SIDE

int kill ( pid_t pid, int sig );


int ipc_state = ipcstate_NO_CHILD;

ipc_data_item ipc_buffer[IPC_BUFFER_LENGTH];
int           ipc_buf_count = 0, ipc_data_size = 0;

geom_object   *ipc_go_bound = NULL;

/* ////////////////////////////////////////////////////////////////////////// */
boolean LaunchAChildProcess ( void )
{
  if ( ipc_state == ipcstate_NO_CHILD ) {
    if ( chdir ( initial_directory ) )
      return false;
    if ( xge_MakeTheChild ( prog_argv[0], "_proc", POZWALAJ_IPC_MAGIC ) ) {
      ipc_state = ipcstate_CHILD_LAUNCHED;
      chdir ( current_directory );
      return true;
    }
  }
  chdir ( current_directory );
  return false;
} /*LaunchAChildProcess*/

void ResetIPCBuffer ( void )
{
  ipc_buf_count = ipc_data_size = 0;
} /*ResetIPCBuffer*/

boolean IPCAppendDataItem ( int desc, int size, void *ptr )
{
  if ( size >= 0 && ipc_buf_count < IPC_BUFFER_LENGTH ) {
    ipc_buffer[ipc_buf_count].desc = desc;
    ipc_buffer[ipc_buf_count].size = size;
    ipc_buffer[ipc_buf_count].ptr = ptr;
    ipc_data_size += size + 2*sizeof(int);
    ipc_buf_count ++;
    return true;
  }
  else
    return false;
} /*IPCAppendDataItem*/

void IPCSendData ( void )
{
  int i;

  for ( i = 0; i < ipc_buf_count; i++ ) {
    write ( xge_pipe_in[1], &ipc_buffer[i].desc, sizeof(int) );
    write ( xge_pipe_in[1], &ipc_buffer[i].size, sizeof(int) );
    if ( ipc_buffer[i].size > 0 )
    write ( xge_pipe_in[1], ipc_buffer[i].ptr, ipc_buffer[i].size );
  }
} /*IPCSendData*/

void IPCWakeUpChild ( void )
{
        /* make the child be ready for reading data of the total size given */
  if ( ipc_data_size ) {
    xge_CallTheChild ( ipccmd_GET_DATA, ipc_data_size );
    xge_PostIdleCommand ( IDLE_COMMAND_SEND_DATA_TO_CHILD, 0, 0 );
  }
} /*IPCWakeUpChild*/

void BindChildToGeomObject ( geom_object *go )
{
  ipc_go_bound = go;
  if ( go )
    go->bound_with_a_child = true;
} /*BindChildToGeomObject*/

void IPCInterruptTheChild ( void )
{
  if ( xge_ChildIsActive () ) {
printf ( "sending the kill signal to %d\n", xge_child_pid );
    kill ( xge_child_pid, SIGKILL );
    ipc_state = ipcstate_NO_CHILD;
    if ( GeomObjectIsInTheList ( ipc_go_bound ) ) {
      ipc_go_bound->bound_with_a_child = false;
      ipc_go_bound = NULL;
    }
  }
} /*IPCInterruptTheChild*/

static void *MallocCopy ( void *ptr, int size )
{
  void *_ptr;

  _ptr = malloc ( size );
  if ( _ptr )
    memcpy ( _ptr, ptr, size );
  return _ptr;
} /*MallocCopy*/

void AssignBSPatchData ( GO_BSplinePatch *obj, int size, char *buf )
{
  int          *ibuf, idesc, isize;
  ipc_blp_size *psize;
  double       *knotsu, *knotsv;
  point3d      *cp, *acp;
  int          lknu, lknv, ncp, i;
  byte         *mkcp, *amkcp;

  psize = NULL;
  cp    = NULL;
  knotsu = knotsv = NULL;
  mkcp  = NULL;
  do {
    ibuf = (int*)buf;
    idesc = ibuf[0];
    isize = ibuf[1];
    size -= 2*sizeof(int);
    buf = (char*)&ibuf[2];
    switch ( idesc ) {
  case ipcd_BLP_SIZE:
      psize = (ipc_blp_size*)buf;
      if ( psize->cpdimen != 3 )
        return;
      break;
  case ipcd_BLP_CPOINTS:
      cp = (point3d*)buf;
      break;
  case ipcd_BLP_MKCP:
      mkcp = (byte*)buf;
      break;
  default:
      break;
    }
    buf = &buf[isize];
    size -= isize;
  } while ( size > 0 );
  if ( psize && cp ) {
    lknu = psize->lastknotu;
    lknv = psize->lastknotv;
    knotsu = malloc ( (lknu+1)*sizeof(double) );
    knotsv = malloc ( (lknv+1)*sizeof(double) );
    ncp = (lknu-psize->degu)*(lknv-psize->degv);
    acp = malloc ( ncp*sizeof(point3d) );
    amkcp = malloc ( ncp*sizeof(char) );
    if ( knotsu && knotsv && acp && amkcp ) {
      for ( i = 0; i <= lknu; i++ ) knotsu[i] = (double)i;
      for ( i = 0; i <= lknv; i++ ) knotsv[i] = (double)i;
      memcpy ( acp, cp, ncp*sizeof(point3d) );
      if ( mkcp )
        memcpy ( amkcp, mkcp, ncp*sizeof(char) );
      else
        memset ( amkcp, MASK_CP_MOVEABLE, ncp );
      GeomObjectAssignBSPatch ( obj, 3, false, psize->degu, lknu, knotsu,
                                psize->degv, lknv, knotsv, (double*)acp, amkcp,
                                false, false );
    }
    else {
      if ( knotsu ) free ( knotsu );
      if ( knotsv ) free ( knotsv );
      if ( acp ) free ( acp );
      if ( amkcp ) free ( amkcp );
    }
  }
} /*AssignBSPatchData*/

void AssignBSMeshData ( GO_BSplineMesh *obj, int size, char *buf )
{
  int          *ibuf, idesc, isize;
  int          nv, nhe, nfac;
  ipc_bsm_size *msize;
  BSMvertex    *mv, *_mv;
  BSMhalfedge  *mhe, *_mhe;
  BSMfacet     *mfac, *_mfac;
  int          *mvhei, *mfhei, *_mvhei, *_mfhei;
  point3d      *mvpc, *_mvpc;
  byte         *mkcp, *_mkcp;

  msize = NULL;
  mv = NULL;  mhe = NULL;  mfac = NULL;  mvhei = mfhei = NULL;
  mvpc = NULL;
  mkcp = NULL;
  nv = nhe = nfac = 0;
  do {
    ibuf = (int*)buf;
    idesc = ibuf[0];
    isize = ibuf[1];
    buf = (char*)(&ibuf[2]);
    size -= 2*sizeof(int);
    switch ( idesc ) {
  case ipcd_BSM_SIZE:
      msize = (ipc_bsm_size*)buf;
      if ( msize->cpdimen != 3 )
        return;
      nv = msize->nv;  nhe = msize->nhe;  nfac = msize->nfac;
      break;
  case ipcd_BSM_VERT:
      mv = (BSMvertex*)buf;
      break;
  case ipcd_BSM_VHE:
      mvhei = (int*)buf;
      break;
  case ipcd_BSM_VERTC:
      mvpc = (point3d*)buf;
      break;
  case ipcd_BSM_VERTMK:  /* this item is optional */
      mkcp = (byte*)buf;
      break;
  case ipcd_BSM_HALFE:
      mhe = (BSMhalfedge*)buf;
      break;
  case ipcd_BSM_FAC:
      mfac = (BSMfacet*)buf;
      break;
  case ipcd_BSM_FHE:
      mfhei = (int*)buf;
      break;
  default:
      break;
    }
    buf = &buf[isize];
    size -= isize;
  } while ( size > 0 );
  if ( mv && mvhei && mvpc && mhe && mfac && mfhei ) {
    _mv = MallocCopy ( mv, nv*sizeof(BSMvertex) );
    _mvhei = MallocCopy ( mvhei, nhe*sizeof(int) );
    _mvpc = MallocCopy ( mvpc, nv*sizeof(point3d) );
    _mhe = MallocCopy ( mhe, nhe*sizeof(BSMhalfedge) );
    _mfac = MallocCopy ( mfac, nfac*sizeof(BSMfacet) );
    _mfhei = MallocCopy ( mfhei, nhe*sizeof(int) );
    if ( mkcp )
      _mkcp = MallocCopy ( mkcp, nv );
    else
      _mkcp = MallocCopy ( obj->mkcp, nv );
    if ( _mv && _mvhei && _mvpc && _mhe && _mfac && _mfhei && _mkcp ) {
      GeomObjectAssignBSplineMesh ( obj, 3, false, nv, _mv, _mvhei, (double*)_mvpc,
                                    nhe, _mhe, nfac, _mfac, _mfhei, _mkcp );
      if ( bsm_sw_log_it ) {
        if ( bsf_OpenOutputFile ( logfilename, true ) ) {
          GeomObjectWriteObj ( (geom_object*)obj );
          bsf_CloseOutputFile ();
        }
      }
    }
    else {
      if ( _mv ) free ( _mv );
      if ( _mvhei ) free ( _mvhei );
      if ( _mvpc ) free ( _mvpc );
      if ( _mhe ) free ( _mhe );
      if ( _mfac ) free ( _mfac );
      if ( _mfhei ) free ( _mfhei );
      if ( _mkcp ) free ( _mkcp );
    }
  }
} /*AssignBSMeshData*/

void ReadDataFromChild ( int size )
{
  void *sp;
  char *buf, b;
  int  *ibuf, s, i;

  sp = pkv_GetScratchMemTop ();
  buf = pkv_GetScratchMem ( size );
  if ( buf ) {
#ifdef DEBUG_IPC
printf ( "receiving result, size = %d\n", size );
#endif
    s = 0;
    do {
      i = read ( xge_pipe_out[0], &buf[s], size );
      s += i;
    } while ( s < size );
    if ( GeomObjectIsInTheList ( ipc_go_bound ) ) {
      if ( ipc_go_bound->bound_with_a_child ) {
        ibuf = (int*)buf;
        switch ( *ibuf ) {
      case ipcd_BLP_SIZE:
          AssignBSPatchData ( (GO_BSplinePatch*)ipc_go_bound, size, buf );
          if ( ipc_state == ipcstate_CHILD_READY &&
                ipc_go_bound == current_go ) {
            ipc_go_bound->bound_with_a_child = false;
            SetupBSplinePatchWidgets ( (GO_BSplinePatch*)ipc_go_bound );
          }
          break;

      case ipcd_BSM_SIZE:  /* mesh data */
          AssignBSMeshData ( (GO_BSplineMesh*)ipc_go_bound, size, buf );
          if ( ipc_state == ipcstate_CHILD_READY &&
               ipc_go_bound == current_go ) {
            ipc_go_bound->bound_with_a_child = false;
            SetupBSplineMeshWidgets ( (GO_BSplineMesh*)ipc_go_bound );
          }
          break;

      default:
          goto way_out;
        }
      }
      else goto way_out;
    }
#ifdef DEBUG_IPC
printf ( "read in %d bytes\n", s );
#endif
  }
  else {  /* reject the data */
    for ( i = 0; i < size; i++ )
      read ( xge_pipe_out[0], &b, 1 );
    goto way_out;
  }
  if ( ipc_go_bound )
    if ( ipc_go_bound == current_go || ipc_go_bound->active ) {
      rendered_picture = false;
      xge_SetWindow ( win0 );
      xge_Redraw ();
    }
way_out:
  pkv_SetScratchMemTop ( sp );
} /*ReadDataFromChild*/

void ProcessChildMessage ( int msg, int size )
{
#ifdef DEBUG_IPC
printf ( "message for PARENT %d, size %d, state %d\n", msg, size, ipc_state );
#endif
  switch ( msg ) {
case ipccmd_CHILD_READY:
    ipc_state = ipcstate_CHILD_READY;
    break;

case ipccmd_BEGIN_BLP:
case ipccmd_BEGIN_BSM:
    ipc_state = ipcstate_CHILD_BUSY;
    break;

case ipccmd_SUCCESS:
    break;

case ipccmd_ERROR:
    ipc_state = ipcstate_CHILD_READY;
    ipc_go_bound->bound_with_a_child = false;
    if ( ipc_go_bound == current_go ) {
      switch ( current_go->obj_type ) {
    case GO_BSPLINE_PATCH:
        SetupBSplinePatchWidgets ( (GO_BSplinePatch*)current_go );
        break;
    case GO_BSPLINE_MESH:
        SetupBSplineMeshWidgets ( (GO_BSplineMesh*)current_go );
        break;
    default:
        break;
      }
    }
    ipc_go_bound = NULL;
    xge_SetWindow ( win1 );
    xge_DisplayErrorMessage ( ErrorMsgOptimizationFailure, -1 );
    break;

case ipccmd_PARTIAL_RESULT:
    ReadDataFromChild ( size );
    break;

case ipccmd_FINAL_RESULT:
    ipc_state = ipcstate_CHILD_READY;
    ReadDataFromChild ( size );
    ipc_go_bound->bound_with_a_child = false;
    ipc_go_bound = NULL;
    xge_SetWindow ( win1 );
    xge_Redraw ();
    break;

default:
    break;
  }
} /*ProcessChildMessage*/

