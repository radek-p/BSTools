
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


/* data buffers for mesh optimization */
static ipc_bsm_size    meshsize, coarsesize;
static ipc_bsm_options meshoptions;

/* ////////////////////////////////////////////////////////////////////////// */
GO_BSplineMesh *BlendingMeshOptimizationFindCoarse ( GO_BSplineMesh *obj )
{
  geom_object *go;
  char        *coarse_name;

  coarse_name = obj->coarse_name;
  for ( go = first_go; go; go = go->next )
    if ( go->obj_type == GO_BSPLINE_MESH &&
         go != (geom_object*)obj &&
         ((GO_BSplineMesh*)go)->nv < obj->nv &&
         strcmp ( coarse_name, go->name ) == 0 )
      return (GO_BSplineMesh*)go;
  return NULL;
} /*BlendingMeshOptimizationFindCoarse*/

boolean BlendingMeshOptimizationPrepareData ( GO_BSplineMesh *obj )
{
  int            nv, nhe, nfac, cpdimen;
  GO_BSplineMesh *coarse;

  BindChildToGeomObject ( (geom_object*)obj );
  ResetIPCBuffer ();
  meshsize.nv = nv = obj->nv;
  meshsize.nhe = nhe = obj->nhe;
  meshsize.nfac = nfac = obj->nfac;
  meshsize.spdimen = obj->me.spdimen;
  meshsize.cpdimen = cpdimen = obj->me.cpdimen;
  meshsize.degree = obj->degree;
  IPCAppendDataItem ( ipcd_BSM_SIZE, sizeof(ipc_bsm_size), &meshsize );
  IPCAppendDataItem ( ipcd_BSM_VERT, nv*sizeof(BSMvertex), obj->meshv );
  IPCAppendDataItem ( ipcd_BSM_VHE, nhe*sizeof(int), obj->meshvhei );
  IPCAppendDataItem ( ipcd_BSM_VERTC, nv*cpdimen*sizeof(double), obj->meshvpc );
  IPCAppendDataItem ( ipcd_BSM_VERTMK, nv*sizeof(char), obj->mkcp );
  IPCAppendDataItem ( ipcd_BSM_HALFE, nhe*sizeof(BSMhalfedge), obj->meshhe );
  IPCAppendDataItem ( ipcd_BSM_FAC, nfac*sizeof(BSMfacet), obj->meshfac );
  IPCAppendDataItem ( ipcd_BSM_FHE, nhe*sizeof(int), obj->meshfhei );
  meshoptions.use_coarse = false;
  if ( obj->bl_use_coarse ) {
    coarse = BlendingMeshOptimizationFindCoarse ( obj );
    if ( coarse ) {
      memset ( &coarsesize, 0, sizeof(ipc_bsm_size) );
      coarsesize.nv = nv = coarse->nv;
      coarsesize.nhe = nhe = coarse->nhe;
      coarsesize.nfac = nfac = coarse->nfac;
      coarsesize.degree = obj->degree;
      IPCAppendDataItem ( ipcd_BSM_COARSE_SIZE, sizeof(ipc_bsm_size), &coarsesize );
      IPCAppendDataItem ( ipcd_BSM_COARSE_VERT, nv*sizeof(BSMvertex), coarse->meshv );
      IPCAppendDataItem ( ipcd_BSM_COARSE_VHE, nhe*sizeof(int), coarse->meshvhei );
      IPCAppendDataItem ( ipcd_BSM_COARSE_HALFE, nhe*sizeof(BSMhalfedge), coarse->meshhe );
      IPCAppendDataItem ( ipcd_BSM_COARSE_FAC, nfac*sizeof(BSMfacet), coarse->meshfac );
      IPCAppendDataItem ( ipcd_BSM_COARSE_FHE, nhe*sizeof(int), coarse->meshfhei );
      meshoptions.use_coarse = true;
    }
  }
  meshoptions.alt_multilevel = bsm_sw_alt_multilevel;
  meshoptions.nkn1 = obj->nkn1;
  meshoptions.nkn2 = obj->nkn2;
  meshoptions.nlevels = obj->nlevels;
  meshoptions.nblocks = obj->nblocks;
  meshoptions.maxit = obj->maxit;
  meshoptions.startbl = obj->startfrom;
  meshoptions.C = obj->bsm_bl_C;
  meshoptions.pretrans = obj->me.pretrans;
  meshoptions.use_constraints = obj->bl_constr;
  meshoptions.shape_only = obj->bl_shape_only;
  meshoptions.constr_mask = marking_mask;
  meshoptions.npthreads = (byte)bsm_npthreads;
  IPCAppendDataItem ( ipcd_BSM_OPTIMIZE, sizeof(ipc_bsm_options), &meshoptions );
  if ( bsm_sw_log_it ) {
    if ( bsf_OpenOutputFile ( logfilename, false ) ) {
      GeomObjectWriteObj ( (geom_object*)obj );
      bsf_CloseOutputFile ();
    }
  }
  return true;
} /*BlendingMeshOptimizationPrepareData*/

void InitBlendingMeshOptimization ( void )
{
  switch ( ipc_state ) {
case ipcstate_NO_CHILD:
    if ( LaunchAChildProcess () )
      xge_PostIdleCommand ( IDLE_COMMAND_BSM_BLENDING_OPT_INIT, 0, 0 );
    else
      xge_DisplayErrorMessage ( ErrorMsgCannotLaunchAChild, 0 );
    break;
case ipcstate_CHILD_LAUNCHED:
    xge_PostIdleCommand ( IDLE_COMMAND_BSM_BLENDING_OPT_INIT, 0, 0 );
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
} /*InitBlendingMeshOptimization*/

