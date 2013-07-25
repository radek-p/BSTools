
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


boolean GeomObjectWriteBSplineMesh ( GO_BSplineMesh *obj )
{
  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return false;
  return bsf_WriteBSMeshd ( obj->me.spdimen, obj->me.cpdimen, obj->rational,
                  obj->degree, obj->nv, obj->meshv, obj->meshvhei, obj->meshvpc,
                  obj->nhe, obj->meshhe, obj->nfac, obj->meshfac, obj->meshfhei,
                  obj->mkcp, obj->me.name );
} /*GeomObjectWriteBSplineMesh*/

void GeomObjectReadBSplineMesh ( void *usrdata,
                    const char *name, int degree,
                    int nv, const BSMvertex *mv, const int *mvhei,
                    const point4d *vc,
                    int nhe, const BSMhalfedge *mhe,
                    int nfac, const BSMfacet *mfac, const int *mfhei,
                    int spdimen, boolean rational, byte *_mkcp )
{
  GO_BSplineMesh *obj;
  BSMvertex      *vert;
  BSMhalfedge    *halfe;
  BSMfacet       *fac;
  double         *vpc;
  int            *vhei, *fhei;
  byte           *mkcp;
  int            cpdimen, i;

  if ( !bsm_CheckMeshIntegrity ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei ) )
    return;
  obj = (GO_BSplineMesh*)GeomObjectAddBSplineMesh ( name, spdimen, rational );
  if ( obj ) {
    cpdimen = obj->me.cpdimen;
    vert = malloc ( nv*sizeof(BSMvertex) );
    halfe = malloc ( nhe*sizeof(BSMhalfedge) );
    fac  = malloc ( nfac*sizeof(BSMfacet) );
    vpc  = malloc ( nv*cpdimen*sizeof(double) );
    vhei = malloc ( nhe*sizeof(int) );
    fhei = malloc ( nhe*sizeof(int) );
    mkcp = malloc ( nv );
    if ( !vert || !halfe || !fac || !vpc || !vhei || !fhei || !mkcp ) {
      if ( vert ) free ( vert );
      if ( halfe ) free ( halfe );
      if ( fac )  free ( fac );
      if ( vpc )  free ( vpc );
      if ( vhei ) free ( vhei );
      if ( fhei ) free ( fhei );
      if ( mkcp ) free ( mkcp );
      GeomObjectDeleteBSplineMesh ( obj );
      return;
    }
    memcpy ( mkcp, _mkcp, nv );
    for ( i = 0; i < nv; i++ )
      mkcp[i] |= MASK_CP_MOVEABLE;
    GeomObjectAssignBSplineMesh ( obj, spdimen, rational,
                      nv, vert, vhei, vpc, nhe, halfe,
                      nfac, fac, fhei, mkcp );
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

