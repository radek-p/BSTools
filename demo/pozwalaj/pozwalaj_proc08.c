
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
#include "mengerc.h"
#include "xgedit.h"

#define CHILD_SIDE
#include "pozwalajipc.h"

#include "pozwalaj_proc.h"

/* /////////////////////////////////////////////////////////////////////////// */
void BeginBSCMCOptimization ( void )
{
  finished = true;
  if ( bscmcsize.cpdimen == 3 && bscmcsize.spdimen == 3 &&
       bscmcsize.degree >= 3 && bscmcsize.closed ) {
    printf ( "Minimization of integral Menger curvature\n" );
    printf ( "deg = %d, lkn = %d, ncp = %d, nvcp = %d, nvars = %d\n",
             bscmcsize.degree, bscmcsize.lkn,
             bscmcsize.lkn-bscmcsize.degree,
             bscmcsize.lkn-2*bscmcsize.degree,
             3*(bscmcsize.lkn-2*bscmcsize.degree) );
    printf ( "exponent = %6.2f\n", bscmcoptions.exponent );
    time0 = times ( &start );
    if ( !mengerc_InitMCOptimization ( bscmcsize.degree, bscmcsize.lkn,
                                       bscmcknots, (point3d*)bscmccp,
                                       bscmcoptions.exponent,
                                       bscmcoptions.pparam, bscmcoptions.nqkn,
                                       bscmcoptions.npthreads, bscmcoptions.ppopt,
                                       &bscmcdata ) ) {
      xge_CallTheParent ( ipccmd_ERROR, 0 );
      return;
    }
    xge_ChildCallYourself ( ipccmd_CONTINUE_BSCMC );
    finished = false;
    itn = 0;
  }
} /*BeginBSCMCOptimization*/

void PrepareBSCMCOutput ( mengerc_data *mcd )
{
  int ncp;

  ncp = bscmcsize.lkn-bscmcsize.degree;
  ResetIPCBuffer ();
  IPCAppendDataItem ( ipcd_BSC_SIZE, sizeof(ipc_bscmc_size), &bscmcsize );
  IPCAppendDataItem ( ipcd_BSC_KNOTS, (bscmcsize.lkn+1)*sizeof(double),
                      bscmcknots );
  IPCAppendDataItem ( ipcd_BSC_CPOINTS, ncp*bscmcsize.cpdimen*sizeof(double),
                      bscmccp );
  memset ( bscmcmkcp, 1, ncp*sizeof(byte) );
  if ( mcd->mdi < ncp )
    bscmcmkcp[mcd->mdi] = 9;
  IPCAppendDataItem ( ipcd_BSC_MKCP, ncp*sizeof(byte), bscmcmkcp );
  bscmcinfo.exponent = bscmcoptions.exponent;
  memcpy ( bscmcinfo.pparam, mcd->penalty_param,
           MENGERC_NPPARAM*sizeof(double) );
  bscmcinfo.imc = mcd->ffkM;
  bscmcinfo.imcp = mcd->ffkMe;
  bscmcinfo.length = mcd->lgt;
  memcpy ( bscmcinfo.pfv, mcd->ffR, MENGERC_NPPARAM*sizeof(double) );
  IPCAppendDataItem ( ipcd_BSC_MCINFO, sizeof(ipc_bscmc_info), &bscmcinfo );
} /*PrepareBSCMCOutput*/

void OutputBSCMCInfo ( mengerc_data *mcd )
{
  int  i;
  char code;

  if ( mcd->ppopt ) {
    printf ( "penalty parameters:\n" );
    for ( i = 0; i < MENGERC_NPPARAM; i++ )
      printf ( "%14.8e ", mcd->penalty_param[i] );
    printf ( "\n" );
  }
  switch ( mcd->itres ) {
case PKN_LMT_ERROR:            code = 'e';  break;
case PKN_LMT_CROSSED_LIMIT:    code = 'l';  break;
case PKN_LMT_CONTINUE_N:       code = '+';  break;
case PKN_LMT_CONTINUE_LM:      code = '-';  break;
case PKN_LMT_CONTINUE_LM_P:    code = '#';  break;
case PKN_LMT_FOUND_MINIMUM:    code = '.';  break;
case PKN_LMT_FOUND_ZEROGRAD:   code = 'z';  break;
case PKN_LMT_FOUND_ZEROGRAD_P: code = 'Z';  break;
case PKN_LMT_FOUND_BARRIER:    code = 'b';  break;
case PKN_LMT_NO_PROGRESS:      code = 's';  break;
default:                       code = ' ';  break;
  }
  printf ( "%3d (%c): %14.8g %14.8g\n", itn, code, mcd->lastf, mcd->gn );
} /*OutputBSCMCInfo*/

void ContinueBSCMCOptimization ( void )
{
  if ( finished )
    return;
  itn ++;
  if ( !mengerc_IterMCOptimization ( &bscmcdata, &finished ) )
    goto failure;
  OutputBSCMCInfo ( &bscmcdata );
  if ( finished || itn >= bscmcoptions.maxiter ) {
    mengerc_UntieTheCurve ( &bscmcdata );
    PrepareBSCMCOutput ( &bscmcdata );
    xge_CallTheParent ( ipccmd_FINAL_RESULT, ipc_data_size );
    xge_ChildCallYourself ( ipccmd_SEND_RESULT );
    time1 = times ( &stop );
    printf ( "time: %6.2fs\n", ((double)(time1-time0))/sysconf(_SC_CLK_TCK) );
    finished = true;
  }
  else {
    PrepareBSCMCOutput ( &bscmcdata );
    xge_CallTheParent ( ipccmd_PARTIAL_RESULT, ipc_data_size );
    xge_ChildCallYourself ( ipccmd_SEND_RESULT );
    xge_ChildCallYourself ( ipccmd_CONTINUE_BSCMC );
  }
  return;

failure:
  mengerc_UntieTheCurve ( &bscmcdata );
  xge_CallTheParent ( ipccmd_ERROR, 0 );
  finished = true;
} /*ContinueBSCMCOptimization*/

