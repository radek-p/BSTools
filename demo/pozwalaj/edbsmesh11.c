
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
#include "eg1holed.h"
#include "eg2holed.h"
#include "bsfile.h"
#include "xgedit.h"
#include "xgledit.h"

#include "editor.h"
#include "editor_bsm.h"


static int _MyOptionProcG1 ( GHoleDomaind *domain, int query, int qn,
                             int *ndata, int **idata, double **fdata )
{
  switch ( query ) {
case G1HQUERY_QUADRATURE:
    return G1H_QUADRATURE_GAUSS_LEGENDRE;
default:
    return G1H_DEFAULT;
  }
} /*_MyOptionProcG1*/

static int _MyOptionProcG2 ( GHoleDomaind *domain, int query, int qn,
                             int *ndata, int **idata, double **fdata )
{
  switch ( query ) {
case G2HQUERY_QUADRATURE:
    return G2H_QUADRATURE_GAUSS_LEGENDRE;
default:
    return G2H_DEFAULT;
  }
} /*_MyOptionProcG2*/

static void _MyOutBezPatch ( int n, int m, const point3d *cp, void *usrptr )
{
  GO_BSplineMesh *obj;
  double         *thepatch;

  obj = (GO_BSplineMesh*)usrptr;
  thepatch = &obj->special_patches[obj->nspecial_patches*(n+1)*(n+1)*3];
  memcpy ( thepatch, cp, (n+1)*(n+1)*sizeof(point3d) );
  obj->nspecial_patches ++;
} /*_MyOutBezPatch*/

void GeomObjectBSplineMeshOptSpecialPatches ( GO_BSplineMesh *obj )
{
  void                  *sp;
  int                   i, j, np, deg, kmax, hole_k, lk;
  bsm_special_elem_list *spvlist;
  bsm_special_el        *spel;
  GHoleDomaind          *domain;
  int                   *vertnum;
  point3d               *hcp;
  int                   ii, l, m, n;
  point3d               *vpc;

printf ( "opt special patches\n" );

  sp = pkv_GetScratchMemTop ();
  domain = NULL;
  obj->special_patches_ok = false;
  if ( obj->me.cpdimen != 3 )
    goto failure;
  if ( !obj->spvlist_ok ) {
    if ( !GeomObjectBSplineMeshFindSpecialVertices ( obj, 2 ) )
      goto failure;
  }
  spvlist = &obj->spvlist;
  spel = obj->spvlist.spel;
  if ( spvlist->nspvert < 1 )
    goto failure;
        /* count the Bezier patches corresponding to the special elements */
  for ( i = np = kmax = 0;  i < spvlist->nspecials;  i++ ) {
    hole_k = spel[i].degree;
    kmax = max ( kmax, hole_k );
    np += hole_k;
  }
        /* what is their degree? */
  deg = obj->special_deg = obj->fill_G2 ? 9 : 5;
        /* allocate a suitable array for these patches */
  obj->nspecial_patches = 0;
  if ( obj->special_patches ) free ( obj->special_patches );
  obj->special_patches = malloc ( np*(deg+1)*(deg+1)*sizeof(point3d) );
  if ( !obj->special_patches )
    goto failure;
        /* optimize patches filling the k-gonal holes */
  hcp = pkv_GetScratchMem ( (12*kmax+1)*sizeof(point3d) );
  if ( !hcp )
    goto failure;
  lk = 0;
  domain = NULL;
  vpc = (point3d*)obj->meshvpc;
  for ( i = 0; i < spvlist->nspecials; i++ ) {
printf ( "<" );
    hole_k = spel[i].degree;
    if ( hole_k != lk ) {
      if ( domain ) gh_DestroyDomaind ( domain );
      j = hole_k == 3 ? 0 : hole_k-4;
      domain = gh_CreateDomaind ( hole_k, NULL, egh_eigendomcpd[j] );
      if ( !domain )
        goto failure;
      lk = hole_k;
      if ( obj->fill_G2 )
        g2h_SetOptionProcd ( domain, _MyOptionProcG2 );
      else
        g1h_SetOptionProcd ( domain, _MyOptionProcG1 );
    }
    vertnum = &spvlist->spvert[spvlist->spel[i].first_snet_vertex];
    memset ( hcp, 0, (12*hole_k+1)*sizeof(point3d) );
    memcpy ( hcp, &vpc[vertnum[0]], sizeof(point3d) );
    for ( ii = 0, m = n = 1;  ii < hole_k;  ii++, n += 3 )
      for ( j = 0;  j < 3;  j++, n++ )
        for ( l = 0;  l < 2;  l++, m++, n++ )
          memcpy ( &hcp[n], &vpc[vertnum[m]], sizeof(point3d) );
    if ( obj->fill_G2 ) {
      if ( !g2h_NLExtFillHoled ( domain, hcp, NULL, obj, _MyOutBezPatch ) )
        goto failure;
    }
    else if ( obj->fill_G1 ) {
      if ( !g1h_NLExtFillHoled ( domain, hcp, NULL, obj, _MyOutBezPatch ) )
        goto failure;
    }
    else if ( obj->fill_G1Q2 ) {
      if ( !g1h_Q2NLExtFillHoled ( domain, hcp, NULL, obj, _MyOutBezPatch ) )
        goto failure;
    }
printf ( ">\n" );
  }
  obj->special_patches_ok = true;
  obj->me.dlistmask &= ~BSM_DLM_HOLEFILL;

failure:
  if ( domain )
    gh_DestroyDomaind ( domain );
  pkv_SetScratchMemTop ( sp );
} /*GeomObjectBSplineMeshOptSpecialPatches*/

