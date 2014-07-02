
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

int ProcessIdleCommand ( int key, short x, short y )
{
  switch ( key ) {
case IDLE_COMMAND_RENDER_START:
    if ( InitRendering () )
      xge_PostIdleCommand ( IDLE_COMMAND_RENDER_CONT, 0, 0 );
    else
      StopRendering ();
    return 1;
case IDLE_COMMAND_RENDER_CONT:
    ContinueRendering ();
    return 1;
case IDLE_COMMAND_BSP_BLENDING_OPT_INIT:
    InitBlendingPatchOptimization ();
    return 1;
case IDLE_COMMAND_BSM_BLENDING_OPT_INIT:
    InitBlendingMeshOptimization ();
    return 1;
case IDLE_COMMAND_BSC_MENGERC_OPT_INIT:
    InitMengerCurvatureOptimization ();
    return 1;
case IDLE_COMMAND_SEND_DATA_TO_CHILD:
    IPCSendData ();
    return 1;
case IDLE_COMMAND_BSMESH_OPT_SPECIALS:
    BlendingMeshOptSpecialPatches ();
    return 1;
default:
    return 0;
  }
} /*ProcessIdleCommand*/

