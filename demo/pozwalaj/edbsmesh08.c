
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

#include "render.h"
#include "editor.h"
#include "editor_bsm.h"


void GeomObjectOutputMeshBSPatchToRenderer3D ( int d, int *vertnum, int *mtab,
                                               void *usrptr )
{
  void           *sp;
  GO_BSplineMesh *obj;
  double         *mcp, *cp, *knu, *knv;
  int            dim, degree, lkn, i, j;

  sp = pkv_GetScratchMemTop ();
  obj = (GO_BSplineMesh*)usrptr;
  dim = obj->me.cpdimen;
  degree = obj->degree;
  if ( degree+1 != d )
    return;
  lkn = 2*degree+1;
  knu = pkv_GetScratchMemd ( lkn+1 );
  knv = pkv_GetScratchMemd ( lkn+1 );
  cp = pkv_GetScratchMemd ( d*d*dim );
  if ( knu && knv && cp ) {
    for ( i = 0; i <= lkn; i++ )
      knu[i] = knv[i] = (double)i;
    mcp = obj->meshvpc;
    for ( i = j = 0;  i < d*d;  i++, j += dim )
      memcpy ( &cp[j], &mcp[dim*vertnum[i]], dim*sizeof(double) );
    if ( dim == 3 )
      RendEnterBSPatch3d ( degree, lkn, knu, degree, lkn, knv, (point3d*)cp,
                           obj->me.colour );
    else if ( dim == 4 )
      RendEnterBSPatch3Rd ( degree, lkn, knu, degree, lkn, knv, (point4d*)cp,
                            obj->me.colour );
  }
  pkv_SetScratchMemTop ( sp );
} /*GeomObjectOutputMeshBSPatchToRenderer3D*/

void GeomObjectBSplineMeshOutputBezPatchToRenderer3D ( int n, int m, const double *cp,
                                                       void *usrptr )
{
  GO_BSplineMesh *obj;
  int            dim;

  obj = (GO_BSplineMesh*)usrptr;
  dim = obj->me.cpdimen;
  if ( dim == 3 )
    RendEnterBezPatch3d ( n, m, (point3d*)cp, obj->me.colour );
  else if ( dim == 4 )
    RendEnterBezPatch3Rd ( n, m, (point4d*)cp, obj->me.colour );
} /*GeomObjectBSplineMeshOutputBezPatchToRenderer3D*/

void GeomObjectBSplineMeshOutputHoleFillingToRenderer3D ( int d, int k,
               int *vertnum, int *mtab, void *usrptr )
{
  GO_BSplineMesh *obj;

  obj = (GO_BSplineMesh*)usrptr;
  GeomObjectBSplineMeshFillBicubicHole ( obj, k, vertnum,
            GeomObjectBSplineMeshOutputBezPatchToRenderer3D );
} /*GeomObjectBSplineMeshOutputHoleFillingToRenderer3D*/

void GeomObjectBSplineMeshOutputToRenderer3D ( GO_BSplineMesh *obj )
{
  bsm_special_elem_list *spvlist;
  bsm_special_el        *spel;
  int                   i, nsp, deg, d, *spvert;
  double                *spv;

  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return;
  if ( obj->me.spdimen != 3 )
    return;
  if ( obj->view_surf )
    bsm_FindRegularSubnets ( obj->nv, obj->meshv, obj->meshvhei,
                             obj->nhe, obj->meshhe,
                             obj->nfac, obj->meshfac, obj->meshfhei,
                             obj->degree+1, obj,
                             GeomObjectOutputMeshBSPatchToRenderer3D );
  if ( obj->subdivision ) {
    /* *************** */
  }
  else if ( obj->degree == 3 && obj->view_holefill ) {
    if ( obj->special_patches_ok ) {
      deg = obj->special_deg;
      d = 3*(deg+1)*(deg+1);
      for ( i = 0, spv = obj->special_patches;
            i < obj->nspecial_patches;
            i++, spv = &spv[d] )
        GeomObjectBSplineMeshOutputBezPatchToRenderer3D ( deg, deg, spv, obj );
    }
    else if ( GeomObjectBSplineMeshFindSpecialVertices ( obj, 2 ) ) {
      spvlist = &obj->spvlist;
      spel    = spvlist->spel;
      spvert  = spvlist->spvert;
      nsp     = spvlist->nspecials;
      for ( i = 0; i < nsp; i++ ) {
        if ( spel[i].el_type == 0 && spel[i].snet_rad == 2 )
          GeomObjectBSplineMeshFillBicubicHole ( obj, spel[i].degree,
              &spvert[spel[i].first_snet_vertex],
              GeomObjectBSplineMeshOutputBezPatchToRenderer3D );
      }
    }
  }
} /*GeomObjectBSplineMeshOutputToRenderer3D*/

