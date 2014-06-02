
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

#define CHILD_SIDE
#include "pozwalajipc.h"
#undef CHILD_SIDE

#include "pozwalaj_proc.h"

ipc_data_item ipc_buffer[IPC_BUFFER_LENGTH];
int           ipc_buf_count = 0, ipc_data_size = 0;

void    *optdata = NULL;
int     itn;
struct  tms start, stop;
clock_t time0, time1;
boolean finished = true;

trans3d         pretrans_inv;

/* blending B-spline patch with optimized shape */
ipc_blp_size    blpsize;
ipc_blp_options blpoptions;
double          *blpcp     = NULL;
double          *_blpcp    = NULL;  /* pretransformed */
char            *blpmkcp   = NULL;
int             blpcpsize, blpmkcpsize;

/* blending mesh surface with optimized shape */
ipc_bsm_size    bsmsize;
ipc_bsm_options bsmoptions;
BSMvertex       *meshv     = NULL;
int             *meshvhei  = NULL;
double          *meshvpc   = NULL;
double          *_meshvpc  = NULL;  /* pretransformed */
char            *meshmkcp  = NULL;
BSMhalfedge     *meshhe    = NULL;
BSMfacet        *meshfac   = NULL;
int             *meshfhei  = NULL;
int             meshvpcsize;
/* coarse mesh to define a preconditioner */
ipc_bsm_size    cmeshsize;
BSMvertex       *cmeshv    = NULL;
int             *cmeshvhei = NULL;
BSMhalfedge     *cmeshhe   = NULL;
BSMfacet        *cmeshfac  = NULL;
int             *cmeshfhei = NULL;
/* refinement matrix - from coarse to fine mesh */
int             rmnnz;
index2          *rmnzi     = NULL;
double          *rmnzc     = NULL;

