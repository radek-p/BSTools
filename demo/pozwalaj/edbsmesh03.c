
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

/*#define NOEXTDEB*/

#define FILL_G1   0
#define FILL_G2   1
#define FILL_G1Q2 2

static GHoleDomaind *domains[GH_MAX_K-3] =
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
static double *g1patchmatrix[GH_MAX_K-3] =
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
static double *g1extpatchmatrix[GH_MAX_K-3] =
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
static double *g2patchmatrix[GH_MAX_K-3] =
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
static double *g2extpatchmatrix[GH_MAX_K-3] =
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
static double *g1q2patchmatrix[GH_MAX_K-3] =
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
static double g1q2const[GH_MAX_K-3] =
  {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
static double *g1q2extpatchmatrix[GH_MAX_K-3] =
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
static double g1q2extconst[GH_MAX_K-3] =
  {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

/* ///////////////////////////////////////////////////////////////////////// */
int GeomObjectBSplineMeshG1HoleOptionProc ( GHoleDomaind *domain,
                              int query, int qn,
                              int *ndata, int **idata, double **fdata )
{
  GO_BSplineMesh *obj;

  switch ( query ) {
case G1HQUERY_QUADRATURE:
    return G1H_QUADRATURE_GAUSS_LEGENDRE;
case G1HQUERY_Q2_FORM_CONSTANT:
    obj = (GO_BSplineMesh*)domain->usrptr;
    if ( obj ) {
      *fdata = &obj->g1q2param;
      return G1H_Q2_USE_SUPPLIED_CONSTANT;
    }
    else
      return G1H_DEFAULT;
default:
    return G1H_DEFAULT;
  }
} /*GeomObjectBSplineMeshG1HoleOptionProc*/

int GeomObjectBSplineMeshG2HoleOptionProc ( GHoleDomaind *domain,
                              int query, int qn,
                              int *ndata, int **idata, double **fdata )
{
  switch ( query ) {
case G2HQUERY_QUADRATURE:
    return G2H_QUADRATURE_GAUSS_LEGENDRE;
default:
    return G2H_DEFAULT;
  }
} /*GeomObjectBSplineMeshG2HoleOptionProc*/

/* ///////////////////////////////////////////////////////////////////////// */
void GeomObjectBSplineMeshSetHoleFillingOpt ( GO_BSplineMesh *obj,
                        boolean fillCoons, boolean fillBezier,
                        boolean fillG1, boolean fillG2, boolean fillG1Q2 )
{
  obj->fill_Coons = fillCoons;
  obj->fill_Bezier = !fillCoons;
  if ( fillG1 )
    obj->fill_G1 = true, obj->fill_G2 = false, obj->fill_G1Q2 = false;
  else if ( fillG2 )
    obj->fill_G1 = false, obj->fill_G2 = true, obj->fill_G1Q2 = false;
  else
    obj->fill_G1 = false, obj->fill_G2 = false, obj->fill_G1Q2 = true;
  obj->me.dlistmask &= ~BSM_DLM_HOLEFILL;
} /*GeomObjectBSplineMeshSetHoleFillingOpt*/

void GeomObjectBSplineMeshSetG1Q2param ( GO_BSplineMesh *obj, double param )
{
  obj->sl_g1q2param = param;
  obj->g1q2param = xge_LogSlidebarValued ( Q2COEFF_MIN, Q2COEFF_MAX, param );
  if ( obj->fill_G1Q2 )
    obj->me.dlistmask &= ~BSM_DLM_HOLEFILL;
} /*GeomObjectBSplineMeshSetG1Q2param*/

/* ///////////////////////////////////////////////////////////////////////// */
static int DomainNum ( int hole_k )
{
  if ( hole_k == 3 )
    return 0;
  else if ( hole_k >= 5 && hole_k <= GH_MAX_K )
    return hole_k-4;
  else
    return -1;
} /*DomainNum*/

GHoleDomaind *GeomObjectBSplineMeshSetupBicubicHoleDomain ( int hole_k )
{
  int          domnum;
  void         *sp;
  int          i;
  double       *domkn;
  GHoleDomaind *domp;

        /* valid numbers of hole sides are 3 and 5,6,...,16 */
  domnum = DomainNum ( hole_k );
  if ( domnum < 0 )
    return NULL;
        /* if the domain has been created before, just return the pointer */
  if ( domains[domnum] )
    return domains[domnum];
        /* else create it */
  sp = pkv_GetScratchMemTop ();
        /* setup equidistant knots */
  domp = NULL;
  domkn = pkv_GetScratchMemd ( hole_k*11 );
  if ( !domkn )
    goto way_out;
  domkn[0] = domkn[1] = 0.0;
  domkn[9] = domkn[10] = 1.0;
  for ( i = 2; i < 9; i++ )
    domkn[i] = (double)(i-1)/8.0;
  for ( i = 1; i < hole_k; i++ )
    memcpy ( &domkn[11*i], domkn, 11*sizeof(double) );
  domains[domnum] = domp = gh_CreateDomaind ( hole_k, domkn,
                                              egh_eigendomcpd[domnum] );
  if ( domp ) {
    g1h_SetOptionProcd ( domp, GeomObjectBSplineMeshG1HoleOptionProc );
    g2h_SetOptionProcd ( domp, GeomObjectBSplineMeshG2HoleOptionProc );
    domp->usrptr = NULL;
  }
way_out:
  pkv_SetScratchMemTop ( sp );
  return domp;
} /*GeomObjectBSplineMeshSetupBicubicHoleDomain*/

double *GeomObjectBSplineMeshSetupPatchMatrix ( int hole_k, char gcont, byte ext,
                                                double g1q2param )
{
  int          domnum;
  GHoleDomaind *dom;
  double       *pmat;

  domnum = DomainNum ( hole_k );
  dom = GeomObjectBSplineMeshSetupBicubicHoleDomain ( hole_k );
  if ( !dom )
    return NULL;
  switch ( gcont ) {
case FILL_G1:
    if ( ext ) {
      pmat = g1extpatchmatrix[domnum];
      if ( !pmat )
        g1extpatchmatrix[domnum] = pmat =
             malloc ( g1h_SymPatchMatrixSize ( hole_k )*sizeof(double) );
      if ( !pmat )
        return NULL;
      if ( !g1h_GetExtSymPatchMatrixd ( dom, pmat ) ) {
        free ( pmat );
        g1extpatchmatrix[domnum] = NULL;
        return NULL;
      }
    }
    else {
      pmat = g1patchmatrix[domnum];
      if ( !pmat )
        g1patchmatrix[domnum] = pmat =
             malloc ( g1h_SymPatchMatrixSize ( hole_k )*sizeof(double) );
      if ( !pmat )
        return NULL;
      if ( !g1h_GetSymPatchMatrixd ( dom, pmat ) ) {
        free ( pmat );
        g1patchmatrix[domnum] = NULL;
        return NULL;
      }
    }
    return pmat;
case FILL_G2:
    if ( ext ) {
      pmat = g2extpatchmatrix[domnum];
      if ( !pmat )
        g2extpatchmatrix[domnum] = pmat =
             malloc ( g2h_SymPatchMatrixSize ( hole_k )*sizeof(double) );
      if ( !pmat )
        return NULL;
#ifdef NOEXTDEB
      if ( !g2h_GetSymPatchMatrixd ( dom, pmat ) ) {
        free ( pmat );
        g2extpatchmatrix[domnum] = NULL;
        return NULL;
      }
#else
      if ( !g2h_GetExtSymPatchMatrixd ( dom, pmat ) ) {
        free ( pmat );
        g2extpatchmatrix[domnum] = NULL;
        return NULL;
      }
#endif
    }
    else {
      pmat = g2patchmatrix[domnum];
      if ( !pmat )
        g2patchmatrix[domnum] = pmat =
             malloc ( g2h_SymPatchMatrixSize ( hole_k )*sizeof(double) );
      if ( !pmat )
        return NULL;
      if ( !g2h_GetSymPatchMatrixd ( dom, pmat ) ) {
        free ( pmat );
        g2patchmatrix[domnum] = NULL;
        return NULL;
      }
    }
    return pmat;
case FILL_G1Q2:
    if ( ext ) {
      pmat = g1q2extpatchmatrix[domnum];
      if ( !pmat ) {
        g1q2extpatchmatrix[domnum] = pmat =
             malloc ( g1h_SymPatchMatrixSize ( hole_k )*sizeof(double) );
        g1q2extconst[domnum] = -1.0;
      }
      if ( !pmat )
        return NULL;
      if ( g1q2extconst[domnum] != g1q2param ) {
        g1h_DestroyQ2PrivateDatad ( dom );
        if ( !g1h_Q2GetExtSymPatchMatrixd ( dom, pmat ) ) {
          free ( pmat );
          g1q2extpatchmatrix[domnum] = NULL;
          return NULL;
        }
        g1q2extconst[domnum] = g1q2param;
      }
    }
    else {
      pmat = g1q2patchmatrix[domnum];
      if ( !pmat ) {
        g1q2extpatchmatrix[domnum] = pmat =
             malloc ( g1h_SymPatchMatrixSize ( hole_k )*sizeof(double) );
        g1q2const[domnum] = -1.0;
      }
      if ( !pmat )
        return NULL;
      if ( g1q2const[domnum] != g1q2param ) {
        g1h_DestroyQ2PrivateDatad ( dom );
        if ( !g1h_Q2GetSymPatchMatrixd ( dom, pmat ) ) {
          free ( pmat );
          g1q2patchmatrix[domnum] = NULL;
          return NULL;
        }
        g1q2const[domnum] = g1q2param;
      }
    }
    return pmat;
default:
    return NULL;
  }
} /*GeomObjectBSplineMeshSetupPatchMatrix*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean GeomObjectBSplineMeshFillBicubicHole ( GO_BSplineMesh *obj,
                          int hole_k, int *vertnum,
                          void (*outpatch)( int n, int m, const double *cp,
                                            void *usrptr ) )
{
  void         *sp;
  GHoleDomaind *domp;
  int          dim, ncp, i, j, l, m, n;
  double       *meshcp, *holecp, *patchmatrix;

  sp = pkv_GetScratchMemTop ();
  domp = GeomObjectBSplineMeshSetupBicubicHoleDomain ( hole_k );
  if ( !domp )
    goto failure;
  if ( domp->hole_k != hole_k )  /* checking just in case */
    goto failure;

        /* extract the control points from the mesh */
  dim = obj->me.cpdimen;
  meshcp = obj->meshvpc;
  ncp = 1+12*hole_k;
  holecp = pkv_GetScratchMemd ( ncp*dim );
  if ( !holecp )
    goto failure;
  memset ( holecp, 0, ncp*dim*sizeof(double) );
  memcpy ( holecp, &meshcp[vertnum[0]*dim], dim*sizeof(double) );
  for ( i = 0, m = n = 1;  i < hole_k;  i++, n += 3 )
    for ( j = 0;  j < 3;  j++, n++ )
      for ( l = 0;  l < 2;  l++, m++, n++ )
        memcpy ( &holecp[n*dim], &meshcp[vertnum[m]*dim], dim*sizeof(double) );

        /* construct the patches filling the hole */
  domp->usrptr = obj;
  if ( obj->fill_G2 ) {
    patchmatrix = GeomObjectBSplineMeshSetupPatchMatrix ( hole_k, FILL_G2,
                                                          obj->fill_Bezier, 0.0 );
    g2h_MatrixFillSymHoled ( hole_k, patchmatrix, dim, holecp, obj, outpatch );
  }
  else if ( obj->fill_G1Q2 ) {
    patchmatrix = GeomObjectBSplineMeshSetupPatchMatrix ( hole_k, FILL_G1Q2,
                                               obj->fill_Bezier, obj->g1q2param );
    g1h_MatrixFillSymHoled ( hole_k, patchmatrix, dim, holecp, obj, outpatch );
  }
  else {
    patchmatrix = GeomObjectBSplineMeshSetupPatchMatrix ( hole_k, FILL_G1,
                                                          obj->fill_Bezier, 0.0 );
    g1h_MatrixFillSymHoled ( hole_k, patchmatrix, dim, holecp, obj, outpatch );
  }
  domp->usrptr = NULL;
  obj->special_patches_ok = false;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*GeomObjectBSplineMeshFillBicubicHole*/

