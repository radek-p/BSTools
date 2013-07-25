
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

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
#include "bsfile.h"
#include "xgedit.h"
#include "xgledit.h"

/*#include "widgets.h"*/
#include "editor.h"  
#include "editor_bezc.h"
#include "editor_bsc.h"
#include "editor_bezp.h"
#include "editor_bsp.h"
#include "editor_bsm.h"
#include "editor_bsh.h"

#define DUMPFILENAME "dump.bs"

boolean GeomObjectOpenDumpFile ( void )
{
  return bsf_OpenOutputFile ( DUMPFILENAME, false );
} /*GeomObjectOpenDumpFile*/

void GeomObjectCloseDumpFile ( void )
{
  bsf_CloseOutputFile ();
  printf ( "%s\n", DUMPFILENAME );
} /*GeomObjectCloseDumpFile*/

void GeomObjectDumpBezierPatch ( int spdim, int cpdim, boolean rational,
                                 int udeg, int vdeg, const double *cp )
{
  bsf_WriteBezierPatchd ( spdim, cpdim, rational,
                          udeg, vdeg, cpdim*(vdeg+1), cp, NULL, "" );
} /*GeomObjectDumpBezierPatch*/

