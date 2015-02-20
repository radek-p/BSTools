
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

#define PARENT_SIDE
#include "pozwalajipc.h"
#undef PARENT_SIDE

static ipc_blp_size    psize;
static ipc_blp_options poptions;


boolean BlendingPatchOptimizationPrepareData ( GO_BSplinePatch *obj )
{
  int ncp, deg;
  char gcont;

  if ( obj->degree_u == 2 )
    gcont = 1;
  else if ( obj->degree_u == 3 )
    gcont = 2;
  else return false;
  if ( obj->me.cpdimen != 3 || obj->rational )
    return false;
  psize.gcont = gcont;
  psize.cpdimen = 3;
  obj->closed_u = psize.closed_u = false /*obj->closed_u*/;  /* for now */
  psize.clamped = false  /*obj->clamped*/;   /* for now */
  switch ( gcont ) {
case 1:
    if ( obj->degree_u != 2 || obj->degree_v != 2 )
      return false;
    deg = psize.degu = psize.degv = 2;
    break;
case 2:
    if ( obj->degree_u != 3 || obj->degree_v != 3 )
      return false;
    deg = psize.degu = psize.degv = 3;
    break;
default:
    return false;
  }
  BindChildToGeomObject ( (geom_object*)obj );
  ResetIPCBuffer ();
  psize.lastknotu = obj->lastknot_u;
  psize.lastknotv = obj->lastknot_v;
  psize.pitch = 3*(obj->lastknot_v-deg);
  IPCAppendDataItem ( ipcd_BLP_SIZE, sizeof(ipc_blp_size), &psize );
  ncp = (obj->lastknot_u-deg)*(obj->lastknot_v-deg);
  IPCAppendDataItem ( ipcd_BLP_CPOINTS, ncp*sizeof(point3d), obj->cpoints );
  IPCAppendDataItem ( ipcd_BLP_MKCP, ncp, obj->mkcp );
  memcpy ( poptions.bl_range, obj->blp_range, 4*sizeof(int) );
  poptions.nkn1 = obj->nkn1;
  poptions.nkn2 = obj->nkn2;
  poptions.maxit = obj->maxit;
  poptions.C = obj->blp_C;
  poptions.pretrans = obj->me.pretrans;
  IPCAppendDataItem ( ipcd_BLP_OPTIMIZE, sizeof(ipc_blp_options), &poptions );
  return true;
} /*BlendingPatchOptimizationPrepareData*/

void InitBlendingPatchOptimization ( void )
{
  switch ( ipc_state ) {
case ipcstate_NO_CHILD:
    if ( LaunchAChildProcess () )
      xge_PostIdleCommand ( IDLE_COMMAND_BSP_BLENDING_OPT_INIT, 0, 0 );
    else
      xge_DisplayErrorMessage ( ErrorMsgCannotLaunchAChild, 0 );
    break;
case ipcstate_CHILD_LAUNCHED:
    xge_PostIdleCommand ( IDLE_COMMAND_BSP_BLENDING_OPT_INIT, 0, 0 );
    break;
case ipcstate_CHILD_READY:
    IPCWakeUpChild ();
    break;
case ipcstate_CHILD_BUSY:
    xge_DisplayErrorMessage ( ErrorMsgChildProcessBusy, 0 );
    break;
default:
    break;
  }
} /*InitBlendingPatchOptimization*/

