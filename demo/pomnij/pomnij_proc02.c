
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <sys/types.h>
#include <sys/times.h>
#include <signal.h>
#include <unistd.h>
#include <setjmp.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <fpu_control.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camerad.h"
#include "multibs.h"
#include "g2blendingd.h"
#include "xgedit.h"

#include "bsfile.h"

#include "pomnijipc.h"
#include "pomnij_proc.h"
#include "proc_regmem.h"

/* ////////////////////////////////////////////////////////////////////////// */
boolean InputBSPatch ( int size )
{
  int esize;

  if ( knotsu ) PKV_FREE ( knotsu )
  if ( knotsv ) PKV_FREE ( knotsv )
  if ( cpoints ) PKV_FREE ( cpoints )
  if ( size < 5*sizeof(int ) )
    return false;
        /* read the integer data from the pipe */
  read ( xge_pipe_in[0], &degu, sizeof(int) );
  read ( xge_pipe_in[0], &degv, sizeof(int) );
  read ( xge_pipe_in[0], &lknu, sizeof(int) );
  read ( xge_pipe_in[0], &lknv, sizeof(int) );
  read ( xge_pipe_in[0], &dim, sizeof(int) );
  if ( degu < 1 || degv < 1 || dim < 1 || dim > 4 ||
       lknu <= 2*degu || lknv <= 2*degv ) {
    return false;
  }
  size -= 5*sizeof(int);
        /* calculate the expected size of data in bytes */
  esize = (lknu+1 + lknv+1 + (lknu-degu)*(lknv-degv)*dim)*sizeof(double);
  if ( esize != size )
    return false;
  PKV_MALLOC ( knotsu, (lknu+1)*sizeof(double) )
  PKV_MALLOC ( knotsv, (lknv+1)*sizeof(double) )
  PKV_MALLOC ( cpoints, (lknu-degu)*(lknv-degv)*dim*sizeof(double) )
  if ( !knotsu || !knotsv || !cpoints )
    return false;
        /* read the data from the pipe */
  read ( xge_pipe_in[0], knotsu, (lknu+1)*sizeof(double) );
  read ( xge_pipe_in[0], knotsv, (lknv+1)*sizeof(double) );
  read ( xge_pipe_in[0], cpoints, (lknu-degu)*(lknv-degv)*dim*sizeof(double) );

/*
printf ( "Received\n" );
*/
/*
bsf_OpenOutputFile ( "proc_test.bs" );
bsf_WriteBSplinePatchd ( dim, degu, lknu, knotsu, degv, lknv, knotsv,
    dim*(lknv-degu), cpoints );
bsf_CloseOutputFile ();
*/
  return true;
} /*InputBSPatch*/

boolean InputOptions ( int size )
{
  int i, j;

        /* size is supposed to be sizeof(ipc_options) */
  if ( size != sizeof(ipc_options) )
    return false;
  read ( xge_pipe_in[0], &options, sizeof(ipc_options) );
/*
printf ( "Received options\n" );
*/
  if ( options.NLConst <= 0.0 || options.maxiter < 1 )
    return false;

  for ( i = 0; i < 3; i++ ) {
    for ( j = 0; j < 4; j++ )
      printf ( "%f ", options.trans.U1.a[i][j] );
    printf ( "\n" );
  }
  return true;
} /*InputOptions*/

boolean InputConstraints ( int size )
{
  void   *sp;
  char   *buf;
  int    *ip;
  double *dp;
  int    s;

  sp = pkv_GetScratchMemTop ();

  if ( ucknots ) PKV_FREE ( ucknots )
  if ( ucurvcp ) PKV_FREE ( ucurvcp )
  if ( cmat ) PKV_FREE ( cmat )
  if ( crhs ) PKV_FREE ( crhs )

  buf = pkv_GetScratchMem ( size );
  if ( !buf )
    goto failure;
  read ( xge_pipe_in[0], buf, size );
  ip = (int*)buf;
  nucknots = *ip;
  if ( nucknots < 1 || nucknots > MAX_BLENDING_CONSTRAINTS )
    goto failure;
        /* it is assumed that the patch has been received before */
  s = nucknots*(lknv-3)*sizeof(point4d);
  if ( sizeof(int) + nucknots*sizeof(double) + s != size )
    goto failure;
  PKV_MALLOC ( ucknots, nucknots*sizeof(double) )
  PKV_MALLOC ( ucurvcp, s )
  if ( !ucknots || !ucurvcp )
    goto failure;
  ip ++;
  dp = (double*)ip;
  memcpy ( ucknots, dp, nucknots*sizeof(double) );
  dp += nucknots;
  memcpy ( ucurvcp, dp, s );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*InputConstraints*/

