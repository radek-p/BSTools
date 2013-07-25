
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
#include "xgeipc.h"

#include "bsfile.h"

#include "pomnijipc.h"
#include "pomnij_proc.h"
#include "proc_regmem.h"

/* ////////////////////////////////////////////////////////////////////////// */
static void _PatchDataSize ( int *s0, int *s1, int *s2, int *s3 )
{
  *s0 = 5*sizeof(int);
  *s1 = (lknu+1)*sizeof(double);
  *s2 = (lknv+1)*sizeof(double);
  *s3 = ((lknu-degu)*(lknv-degv)*dim)*sizeof(double);
} /*_PatchDataSize*/

int PatchDataSize ( void )
{
  int s0, s1, s2, s3;

  _PatchDataSize ( &s0, &s1, &s2, &s3 );
  return s0 + s1 + s2 + s3;
} /*PatchDataSize*/

void OutputBSPatch ( void )
{
  int s0, s1, s2, s3;
  int intpar[5];

  _PatchDataSize ( &s0, &s1, &s2, &s3 );
  intpar[0] = degu;
  intpar[1] = degv;
  intpar[2] = lknu;
  intpar[3] = lknv;
  intpar[4] = dim;
  write ( xge_pipe_out[1], intpar, s0 ); 
  write ( xge_pipe_out[1], knotsu, s1 );
  write ( xge_pipe_out[1], knotsv, s2 );
  write ( xge_pipe_out[1], cpoints, s3 );
} /*OutputBSPatch*/

int ConstraintsSize ( void )
{
  return sizeof(int) + nucknots*sizeof(double) +
         nucknots*(lknv-degv)*sizeof(point4d);
} /*ConstraintsSize*/

void OutputConstraints ( void )
{
  write ( xge_pipe_out[1], &nucknots, sizeof(int) );
  write ( xge_pipe_out[1], ucknots, nucknots*sizeof(double) );
  write ( xge_pipe_out[1], ucurvcp, nucknots*(lknv-degv)*sizeof(point4d) );
} /*OutputConstraints*/

