
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
void ReadBLPSize ( int isize )
{
  read ( xge_pipe_in[0], &blpsize, sizeof(ipc_blp_size) );
} /*ReadBLPSize*/

void ReadBLPCPoints ( int isize )
{
  if ( blpcp ) PKV_FREE ( blpcp );
  if ( _blpcp ) PKV_FREE ( _blpcp );
  PKV_MALLOC ( blpcp, isize );
  if ( blpcp ) {
    read ( xge_pipe_in[0], blpcp, isize );
    blpcpsize = isize;
  }
} /*ReadBLPCPoints*/

void ReadBLPMkcp ( int isize )
{
  if ( blpmkcp ) PKV_FREE ( blpmkcp );
  PKV_MALLOC ( blpmkcp, isize );
  if ( blpmkcp ) {
    read ( xge_pipe_in[0], blpmkcp, isize );
    blpmkcpsize = isize;
  }
} /*ReadBLPMkcp*/

void ReadBLPOptimizeOptions ( int isize )
{
  read ( xge_pipe_in[0], &blpoptions, sizeof(ipc_blp_options) );
  xge_CallTheParent ( ipccmd_BEGIN_BLP, 0 );
  xge_ChildCallYourself ( ipccmd_BEGIN_BLP );
} /*ReadBLPOptimizeOptions*/

