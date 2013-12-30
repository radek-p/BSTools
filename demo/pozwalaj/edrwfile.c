
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

#include "editor.h"
#include "editor_bezc.h"
#include "editor_bsc.h"
#include "editor_bezp.h"
#include "editor_bsp.h"
#include "editor_bsm.h"
#include "editor_bsh.h"

/* ////////////////////////////////////////////////////////////////////////// */
void GeomObjectBeginReading ( void *usrdata, int obj_type )
{
  rw_object_attributes *attrib;

      /* initialise the attributes to default values */
  attrib = (rw_object_attributes*)usrdata;
        /* colour: */
  attrib->colour[0] = attrib->colour[1] = 1.0;
  attrib->colour[2] = 0.0;
        /* others - in future */
} /*GeomObjectBeginReading*/

void GeomObjectEndReading ( void *usrdata, boolean success )
{
  rw_object_attributes *attrib;

  if ( success ) {
      /* the object has been read in, it is the last in the list, */
      /* assign the attributes */
    attrib = (rw_object_attributes*)usrdata;
        /* colour */
    memcpy ( last_go->colour, attrib->colour, 3*sizeof(double) );
        /* others - in future */
  }
} /*GeomObjectEndReading*/

void GeomObjectReadColour ( void *usrdata, point3d *colour )
{
  rw_object_attributes *attrib;

  attrib = (rw_object_attributes*)usrdata;
  memcpy ( attrib->colour, colour, 3*sizeof(double) );
} /*GeomObjectReadColour*/

boolean GeomObjectReadFile ( char *filename, bsf_Camera_fptr CameraReader )
{
  bsf_UserReaders      readers;
  rw_object_attributes attrib;

        /* register procedures entering the objects read in */
  bsf_ClearReaders ( &readers );
  readers.userData = (void*)&attrib;
  bsf_BeginReadingFuncd ( &readers, GeomObjectBeginReading );
  bsf_EndReadingFuncd ( &readers, GeomObjectEndReading );
  bsf_BC4ReadFuncd  ( &readers, GeomObjectReadBezierCurve, MAX_DEGREE );
  bsf_BSC4ReadFuncd ( &readers, GeomObjectReadBSplineCurve,
                      MAX_DEGREE, MAX_BSC_KNOTS-1 );
  bsf_BP4ReadFuncd  ( &readers, GeomObjectReadBezierPatch, MAX_DEGREE );
  bsf_BSP4ReadFuncd ( &readers, GeomObjectReadBSplinePatch,
                      MAX_DEGREE, MAX_BSP_KNOTS-1 );
  bsf_BSM4ReadFuncd ( &readers, GeomObjectReadBSplineMesh, MAX_DEGREE,
                      MAX_BSM_NV, MAX_BSM_NHE, MAX_BSM_NFAC );
  bsf_BSH4ReadFuncd ( &readers, GeomObjectReadBSplineHole );
  bsf_CameraReadFuncd ( &readers, CameraReader );
  bsf_ColourReadFuncd ( &readers, GeomObjectReadColour );
        /* read the file */
  return bsf_ReadBSFiled ( filename, &readers );
} /*GeomObjectReadFile*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean GeomObjectWriteObj ( geom_object *go )
{
  switch ( go->obj_type ) {
case GO_BEZIER_CURVE:
    if ( !GeomObjectWriteBezierCurve ( (GO_BezierCurve*)go ) )
      return false;
    break;
case GO_BEZIER_PATCH:
    if ( !GeomObjectWriteBezierPatch ( (GO_BezierPatch*)go ) )
      return false;
    break;
case GO_BSPLINE_CURVE:
    if ( !GeomObjectWriteBSplineCurve ( (GO_BSplineCurve*)go ) )
      return false;
    break;
case GO_BSPLINE_PATCH:
    if ( !GeomObjectWriteBSplinePatch ( (GO_BSplinePatch*)go ) )
      return false;
    break;
case GO_BSPLINE_MESH:
    if ( !GeomObjectWriteBSplineMesh ( (GO_BSplineMesh*)go ) )
      return false;
    break;
case GO_BSPLINE_HOLE:
    if ( !GeomObjectWriteBSplineHole ( (GO_BSplineHole*)go ) )
      return false;
    break;
default:
    break;
  }
  return true;
} /*GeomObjectWriteObj*/

boolean GeomObjectWriteFile ( char *filename, char whattowrite,
                              boolean (*writeotherdata)( void *usrdata ),
                              void *usrdata,
                              boolean append )
{
  geom_object *go;

  if ( !bsf_OpenOutputFile ( filename, append ) )
    return false;

  switch ( whattowrite ) {
case GO_WRITE_CURRENT:
    if ( !GeomObjectWriteObj ( current_go ) )
      goto failure;
    break;

case GO_WRITE_ACTIVE:
    for ( go = first_go; go; go = go->next )
      if ( go->active || go == current_go ) {
        if ( !GeomObjectWriteObj ( go ) )
          goto failure;
      }
    break;

case GO_WRITE_ALL:
    for ( go = first_go; go; go = go->next )
      if ( !GeomObjectWriteObj ( go ) )
        goto failure;
    break;

default:
    goto failure;
  }

  if ( writeotherdata ) {
    if ( !writeotherdata ( usrdata ) )
      goto failure;
  }

  bsf_CloseOutputFile ();
  return true;

failure:
  bsf_CloseOutputFile ();
  return false;
} /*GeomObjectWriteFile*/

