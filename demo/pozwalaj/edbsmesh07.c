
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
#include "editor_bsm.h"


boolean GeomObjectWriteBSMAttributes ( GO_BSplineMesh *obj )
{
  return true;
} /*GeomObjectWriteBSMAttributes*/

boolean GeomObjectWriteBSplineMesh ( GO_BSplineMesh *obj )
{
  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return false;
  return bsf_WriteBSMeshd ( obj->me.spdimen, obj->me.cpdimen, obj->rational,
                  obj->degree, obj->nv, obj->meshv, obj->meshvhei, obj->meshvpc,
                  obj->nhe, obj->meshhe, obj->nfac, obj->meshfac, obj->meshfhei,
                  obj->me.name, obj->me.ident,
                  GeomObjectWriteAttributes,
                  (void*)&obj->me );
} /*GeomObjectWriteBSplineMesh*/

boolean GeomObjectBSMResolveDependencies ( GO_BSplineMesh *obj )
{
  return true;
} /*GeomObjectBSMResolveDependencies*/

void GeomObjectReadBSplineMesh ( void *usrdata,
                    const char *name, int ident, int degree,
                    int nv, const BSMvertex *mv, const int *mvhei,
                    const point4d *vc,
                    int nhe, const BSMhalfedge *mhe,
                    int nfac, const BSMfacet *mfac, const int *mfhei,
                    int spdimen, boolean rational )
{
  GO_BSplineMesh       *obj;
  BSMvertex            *vert;
  BSMhalfedge          *halfe;
  BSMfacet             *fac;
  double               *vpc;
  int                  *vhei, *fhei;
  byte                 *mkcp, *mkhe, *mkfac;
  int                  cpdimen;
  rw_object_attributes *attrib;

  attrib = (rw_object_attributes*)usrdata;
  obj = (GO_BSplineMesh*)attrib->go_being_read;
  if ( !bsm_CheckMeshIntegrity ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei ) ) {
    attrib->integrity_ok = false;
    return;
  }
  if ( obj ) {
    strncpy ( obj->me.name, name, MAX_NAME_LENGTH+1 );
    obj->me.spdimen = spdimen;
    obj->me.cpdimen = rational ? spdimen+1 : spdimen;
    obj->me.ident = ident;
    cpdimen = obj->me.cpdimen;
    vert = malloc ( nv*sizeof(BSMvertex) );
    halfe = malloc ( nhe*sizeof(BSMhalfedge) );
    fac  = malloc ( nfac*sizeof(BSMfacet) );
    vpc  = malloc ( nv*cpdimen*sizeof(double) );
    vhei = malloc ( nhe*sizeof(int) );
    fhei = malloc ( nhe*sizeof(int) );
    mkcp = malloc ( nv );
    mkhe = malloc ( nhe );
    mkfac = malloc ( nfac );
    if ( !vert || !halfe || !fac || !vpc || !vhei || !fhei ||
         !mkcp || !mkhe || !mkfac ) {
      if ( vert ) free ( vert );
      if ( halfe ) free ( halfe );
      if ( fac )  free ( fac );
      if ( vpc )  free ( vpc );
      if ( vhei ) free ( vhei );
      if ( fhei ) free ( fhei );
      if ( mkcp ) free ( mkcp );
      if ( mkhe ) free ( mkhe );
      if ( mkfac ) free ( mkfac );
      GeomObjectDeleteBSplineMesh ( obj );
      return;
    }
    GeomObjectAssignBSplineMesh ( obj, spdimen, rational,
                      nv, vert, vhei, vpc, nhe, halfe,
                      nfac, fac, fhei, mkcp, mkhe, mkfac );
    memcpy ( vert, mv, nv*sizeof(BSMvertex) );
    memcpy ( vhei, mvhei, nhe*sizeof(int) );
    GeomObjectSetupIniPoints ( spdimen, rational, &obj->me.cpdimen,
                               nv, (double*)vc, vpc );
    memcpy ( halfe, mhe, nhe*sizeof(BSMhalfedge) );
    memcpy ( fac, mfac, nfac*sizeof(BSMfacet) );
    memcpy ( fhei, mfhei, nhe*sizeof(int) );
  }
} /*GeomObjectReadBSplineMesh*/

void GeomObjectBSplineMeshDisplayInfoText ( GO_BSplineMesh *obj )
{
  char s[160];

  sprintf ( s, "nv = %d, nhe = %d, nfac = %d", obj->nv, obj->nhe, obj->nfac );
  SetStatusText ( s, true );
} /*GeomObjectBSplineMeshDisplayInfoText*/

