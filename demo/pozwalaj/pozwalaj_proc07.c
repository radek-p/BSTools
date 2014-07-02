
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
void ReadBSCSize ( int isize )
{
  read ( xge_pipe_in[0], &bscmcsize, sizeof(ipc_bscmc_size) );
} /*ReadBSCSize*/

void ReadBSCKnots ( int isize )
{
  if ( bscmcknots ) PKV_FREE ( bscmcknots );
  PKV_MALLOC ( bscmcknots, isize );
  if ( bscmcknots ) {
    read ( xge_pipe_in[0], bscmcknots, isize );
    bscmcmkcpsize = isize;
  }
} /*ReadBSCKnots*/

void ReadBSCCPoints ( int isize )
{
  if ( bscmccp ) PKV_FREE ( bscmccp );
  PKV_MALLOC ( bscmccp, isize );
  if ( bscmccp ) {
    read ( xge_pipe_in[0], bscmccp, isize );
    bscmccpsize = isize;
  }
} /*ReadBSCCPoints*/

void ReadBSCMkcp ( int isize )
{
  if ( bscmcmkcp ) PKV_FREE ( bscmcmkcp );
  PKV_MALLOC ( bscmcmkcp, isize );
  if ( bscmcmkcp ) {
    read ( xge_pipe_in[0], bscmcmkcp, isize );
    bscmcmkcpsize = isize;
  }
} /*ReadBSCMkcp*/

void ReadBSCMCOptimizeOptions ( int isize )
{
  read ( xge_pipe_in[0], &bscmcoptions, sizeof(ipc_bscmc_options) );
  xge_CallTheParent ( ipccmd_BEGIN_BSCMC, 0 );
  xge_ChildCallYourself ( ipccmd_BEGIN_BSCMC );
} /*ReadBSCMCOptimizeOptions*/

