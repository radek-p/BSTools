
/* ///////////////////////////////////////////////////  */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* ///////////////////////////////////////////////////  */

#include <sys/types.h>
#include <sys/times.h>
#include <signal.h>
#include <unistd.h>
#include <setjmp.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fpu_control.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "eg2holed.h"
#include "bsmesh.h"
#include "g1blendingd.h"
#include "g2blendingd.h"
#include "g2mblendingd.h"
#include "xgedit.h"
#include "xgeipc.h"

#define CHILD_SIDE
#include "pozwalajipc.h"
#undef CHILD_SIDE

#include "pozwalaj_proc.h"

/* /////////////////////////////////////////////////////////////////////////// */
void BeginBLPOptimization ( void )
{
  int    optlknu, optlknv, fcp;
  double *occp;

  finished = true;
  if ( blpsize.cpdimen == 3 ) {
    time0 = times ( &start );
    PKV_MALLOC ( _blpcp, blpcpsize );
    if ( !_blpcp )
      goto failure;
    memcpy ( _blpcp, blpcp, blpcpsize );
    if ( !InvertPretransformation ( &blpoptions.pretrans ) )
      goto failure;
    switch ( blpsize.gcont ) {
  case 1:
      optlknu = blpoptions.bl_range[1]-blpoptions.bl_range[0]+7;
      optlknv = blpoptions.bl_range[3]-blpoptions.bl_range[2]+7;
      fcp = (blpoptions.bl_range[0]-2)*blpsize.pitch +
            (blpoptions.bl_range[2]-2)*blpsize.cpdimen;
      occp = &_blpcp[fcp];      
      TransformCPoints ( &blpoptions.pretrans, optlknu-2, optlknv-2,
                         blpsize.pitch, (point3d*)occp );
      if ( g1bl_InitBlSurfaceOptLMTd ( optlknu, optlknv,
                    blpsize.pitch, (point3d*)occp, blpoptions.C,
                    0.0, 0.0, blpoptions.nkn1, blpoptions.nkn2, &optdata ) ) {
        xge_ChildCallYourself ( ipccmd_CONTINUE_BLP1 );
        finished = false;
        itn = 0;
      }
      return;
  case 2:
      optlknu = blpoptions.bl_range[1]-blpoptions.bl_range[0]+10;
      optlknv = blpoptions.bl_range[3]-blpoptions.bl_range[2]+10;
      fcp = (blpoptions.bl_range[0]-3)*blpsize.pitch +
            (blpoptions.bl_range[2]-3)*blpsize.cpdimen;
      occp = &_blpcp[fcp];      
      TransformCPoints ( &blpoptions.pretrans, optlknu-3, optlknv-3,
                         blpsize.pitch, (point3d*)occp );
      if ( g2bl_InitBlSurfaceOptLMTd ( optlknu, optlknv,
                    blpsize.pitch, (point3d*)occp, blpoptions.C,
                    0.0, 0.0, blpoptions.nkn1, blpoptions.nkn2, &optdata ) ) {
        xge_ChildCallYourself ( ipccmd_CONTINUE_BLP2 );
        finished = false;
        itn = 0;
      }
      return;
  default:
      break;
    }
  }
failure:
  xge_CallTheParent ( ipccmd_ERROR, 0 );
} /*BeginBLPOptimization*/

void PrepareBLPOutput ( void )
{
  int    optlknu, optlknv, fcp;
  double *occp;

  memcpy ( blpcp, _blpcp, blpcpsize );
  switch ( blpsize.gcont ) {
case 1:  /* quadratic */
    optlknu = blpoptions.bl_range[1]-blpoptions.bl_range[0]+7;
    optlknv = blpoptions.bl_range[3]-blpoptions.bl_range[2]+7;
    fcp = (blpoptions.bl_range[0]-2)*blpsize.pitch +
          (blpoptions.bl_range[2]-2)*blpsize.cpdimen;
    occp = &blpcp[fcp];      
    TransformCPoints ( &pretrans_inv, optlknu-2, optlknv-2,
                       blpsize.pitch, (point3d*)occp );
    break;
case 2:  /* cubic */
    optlknu = blpoptions.bl_range[1]-blpoptions.bl_range[0]+10;
    optlknv = blpoptions.bl_range[3]-blpoptions.bl_range[2]+10;
    fcp = (blpoptions.bl_range[0]-3)*blpsize.pitch +
          (blpoptions.bl_range[2]-3)*blpsize.cpdimen;
    occp = &blpcp[fcp];      
    TransformCPoints ( &pretrans_inv, optlknu-3, optlknv-3,
                       blpsize.pitch, (point3d*)occp );
    break;
default:
    return;
  }
  ResetIPCBuffer ();
  IPCAppendDataItem ( ipcd_BLP_SIZE, sizeof(ipc_blp_size), &blpsize );
  IPCAppendDataItem ( ipcd_BLP_CPOINTS, blpcpsize, blpcp );
  IPCAppendDataItem ( ipcd_BLP_MKCP, blpmkcpsize, blpmkcp );
} /*PrepareBLPOutput*/

void ContinueBLPOptimization1 ( void )
{
  if ( finished )
    return;
  itn ++;
  if ( g1bl_IterBlSurfaceOptLMTd ( optdata, &finished ) ) {
    if ( finished || itn >= blpoptions.maxit ) {
      g1bl_OptLMTDeallocated ( &optdata );
      PrepareBLPOutput ();
      xge_CallTheParent ( ipccmd_FINAL_RESULT, ipc_data_size );
      xge_ChildCallYourself ( ipccmd_SEND_RESULT );
      time1 = times ( &stop );
      printf ( "time: %6.2fs\n", ((double)(time1-time0))/sysconf(_SC_CLK_TCK) );
      finished = true;
    }
    else {
      PrepareBLPOutput ();
      xge_CallTheParent ( ipccmd_PARTIAL_RESULT, ipc_data_size );
      xge_ChildCallYourself ( ipccmd_SEND_RESULT );
      xge_ChildCallYourself ( ipccmd_CONTINUE_BLP1 );
    }
  }
  else {
    g1bl_OptLMTDeallocated ( &optdata );
    xge_CallTheParent ( ipccmd_ERROR, 0 );
    finished = true;
  }
} /*ContinueBLPOptimization1*/

void ContinueBLPOptimization2 ( void )
{
  if ( finished )
    return;
  itn ++;
  if ( g2bl_IterBlSurfaceOptLMTd ( optdata, &finished ) ) {
    if ( finished || itn >= blpoptions.maxit ) {
      g2bl_OptLMTDeallocated ( &optdata );
      PrepareBLPOutput ();
      xge_CallTheParent ( ipccmd_FINAL_RESULT, ipc_data_size );
      xge_ChildCallYourself ( ipccmd_SEND_RESULT );
      time1 = times ( &stop );
      printf ( "time: %6.2fs\n", ((double)(time1-time0))/sysconf(_SC_CLK_TCK) );
      finished = true;
    }
    else {
      PrepareBLPOutput ();
      xge_CallTheParent ( ipccmd_PARTIAL_RESULT, ipc_data_size );
      xge_ChildCallYourself ( ipccmd_SEND_RESULT );
      xge_ChildCallYourself ( ipccmd_CONTINUE_BLP2 );
    }
  }
  else {
    g2bl_OptLMTDeallocated ( &optdata );
    xge_CallTheParent ( ipccmd_ERROR, 0 );
    finished = true;
  }
} /*ContinueBLPOptimization2*/

