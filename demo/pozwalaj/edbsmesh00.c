
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


/* tetrahedron mesh */
#define BSM_TETRAHEDRON_NV     4
#define BSM_TETRAHEDRON_NHE   12
#define BSM_TETRAHEDRON_NFAC   4
static const BSMvertex tetrahvert[BSM_TETRAHEDRON_NV] =
  {{3,0},{3,3},{3,6},{3,9}};
static const int tetrahvhei[BSM_TETRAHEDRON_NHE] =
  {1,4,10, 3,0,6, 5,2,8, 9,7,11};
static const point4d tetrahvc[BSM_TETRAHEDRON_NV] =
  {{-1.0,-1.0,-1.0},{1.0,1.0,-1.0},{-1.0,1.0,1.0},{1.0,-1.0,1.0}};
static const BSMhalfedge tetrahhe[BSM_TETRAHEDRON_NHE] =
  {{1,0,0,1},{0,1,1,0},{2,1,0,3},{1,2,2,2},{0,2,0,5},{2,0,3,4},
   {1,3,1,7},{3,1,2,6},{2,3,2,9},{3,2,3,8},{0,3,3,11},{3,0,1,10}};
static const BSMfacet tetrahfac[BSM_TETRAHEDRON_NFAC] =
  {{3,0},{3,3},{3,6},{3,9}};
static const int tetrahfhei[BSM_TETRAHEDRON_NHE] =
  {0,4,2, 1,6,11, 3,8,7, 5,10,9};

/* cube mesh */
#define BSM_CUBE_NV     8
#define BSM_CUBE_NHE   24
#define BSM_CUBE_NFAC   6
static const BSMvertex cubevert[BSM_CUBE_NV] =
  {{3,0},{3,3},{3,6},{3,9},
   {3,12},{3,15},{3,18},{3,21}};
static const int cubevhei[BSM_CUBE_NHE] =
  {1,9,4, 2,11,5, 3,13,6, 0,15,7,
   8,17,20, 10,18,21, 12,19,22, 14,16,23};
static const point4d cubevc[BSM_CUBE_NV] =
  {{-1.0,-1.0,-1.0,1.0},{1.0,-1.0,-1.0,1.0},{1.0,1.0,-1.0,1.0},{-1.0,1.0,-1.0,1.0},
    {-1.0,-1.0,1.0,1.0},{1.0,-1.0,1.0,1.0},{1.0,1.0,1.0,1.0},{-1.0,1.0,1.0,1.0}};
static const BSMhalfedge cubehe[BSM_CUBE_NHE] = {
  {3,0,0,4},{0,1,0,5},{1,2,0,6},{2,3,0,7},
  {0,3,4,0},{1,0,1,1},{2,1,2,2},{3,2,3,3},
  {4,0,4,9},{0,4,1,8},{5,1,1,11},{1,5,2,10},
  {6,2,2,13},{2,6,3,12},{7,3,3,15},{3,7,4,14},
  {7,4,4,20},{4,5,1,21},{5,6,2,22},{6,7,3,23},
  {4,7,5,16},{5,4,5,17},{6,5,5,18},{7,6,5,19}};
static const BSMfacet cubefac[BSM_CUBE_NFAC] =
  {{4,0},{4,4},{4,8},{4,12},{4,16},{4,20}};
static const int cubefachei[BSM_CUBE_NHE] =
  {0,1,2,3, 9,17,10,5, 11,18,12,6,
   13,19,14,7, 15,16,8,4, 20,23,22,21};


/* dodecahedron mesh */
#define BSM_DODECAHEDRON_NV    20
#define BSM_DODECAHEDRON_NHE   60
#define BSM_DODECAHEDRON_NFAC  12
#define DOD_D 1.309016994374947
#define DOD_E (DOD_D-0.5)
static const BSMvertex dodecvert[BSM_DODECAHEDRON_NV] =
  {{3,0},{3,3},{3,6},{3,9},{3,12},
   {3,15},{3,18},{3,21},{3,24},{3,27},
   {3,30},{3,33},{3,36},{3,39},{3,42},
   {3,45},{3,48},{3,51},{3,54},{3,57}};
static const int dodecvhei[BSM_DODECAHEDRON_NHE] =
  {0,6,25,   1,11,5,   2,16,10,  3,21,15,  4,26,20,
   30,36,55, 31,41,35, 32,46,40, 33,51,45, 34,56,50,
   8,49,52,  48,9,12,  13,44,47, 43,14,17, 18,39,42,
   38,19,22, 23,59,37, 58,24,27, 28,54,57, 53,29,7};
static const point4d dodecvc[BSM_DODECAHEDRON_NV] =
  {{0.0,0.5,-DOD_D,1.0}, /*  0 */
   {0.0,-0.5,-DOD_D,1.0}, /*  1 */
   {-DOD_E,-DOD_E,-DOD_E,1.0}, /*  2 */
   {-DOD_D,0.0,-0.5,1.0}, /*  3 */
   {-DOD_E,DOD_E,-DOD_E,1.0}, /*  4 */
   {0.0,0.5,DOD_D,1.0}, /*  5 */
   {0.0,-0.5,DOD_D,1.0}, /*  6 */
   {DOD_E,-DOD_E,DOD_E,1.0}, /*  7 */
   {DOD_D,0.0,0.5,1.0}, /*  8 */
   {DOD_E,DOD_E,DOD_E,1.0}, /*  9 */
   {DOD_D,0.0,-0.5,1.0}, /* 10 */
   {DOD_E,-DOD_E,-DOD_E,1.0}, /* 11 */
   {0.5,-DOD_D,0.0,1.0}, /* 12 */
   {-0.5,-DOD_D,0.0,1.0}, /* 13 */
   {-DOD_E,-DOD_E,DOD_E,1.0}, /* 14 */
   {-DOD_D,0.0,0.5,1.0}, /* 15 */
   {-DOD_E,DOD_E,DOD_E,1.0}, /* 16 */
   {-0.5,DOD_D,0.0,1.0}, /* 17 */
   {0.5,DOD_D,0.0,1.0}, /* 18 */
   {DOD_E,DOD_E,-DOD_E,1.0}};/* 19 */
static const BSMhalfedge dodeche[BSM_DODECAHEDRON_NHE] =
  {{0,1,0,5},{1,2,0,10},{2,3,0,15},{3,4,0,20},{4,0,0,25},
   {1,0,1,0},{0,19,1,29},{19,10,1,52},{10,11,1,48},{11,1,1,11},
   {2,1,2,1},{1,11,2,9},{11,12,2,47},{12,13,2,43},{13,2,2,16},
   {3,2,3,2},{2,13,3,14},{13,14,3,42},{14,15,3,38},{15,3,3,21},
   {4,3,4,3},{3,15,4,19},{15,16,4,37},{16,17,4,58},{17,4,4,26},
   {0,4,5,4},{4,17,5,24},{17,18,5,57},{18,19,5,53},{19,0,5,6},
   {5,6,6,35},{6,7,6,40},{7,8,6,45},{8,9,6,50},{9,5,6,55},
   {6,5,7,30},{5,16,7,59},{16,15,7,22},{15,14,7,18},{14,6,7,41},
   {7,6,8,31},{6,14,8,39},{14,13,8,17},{13,12,8,13},{12,7,8,46},
   {8,7,9,32},{7,12,9,44},{12,11,9,12},{11,10,9,8},{10,8,9,51},
   {9,8,10,33},{8,10,10,49},{10,19,10,7},{19,18,10,28},{18,9,10,56},
   {5,9,11,34},{9,18,11,54},{18,17,11,27},{17,16,11,23},{16,5,11,36}};
static const BSMfacet dodecfac[BSM_DODECAHEDRON_NFAC] =
  {{5,0},{5,5},{5,10},{5,15},{5,20},{5,25},
   {5,30},{5,35},{5,40},{5,45},{5,50},{5,55}};
static const int dodecfhei[BSM_DODECAHEDRON_NHE] =
  {0,1,2,3,4,      5,6,7,8,9,      10,11,12,13,14, 15,16,17,18,19,
   20,21,22,23,24, 25,26,27,28,29, 30,31,32,33,34, 35,36,37,38,39,
   40,41,42,43,44, 45,46,47,48,49, 50,51,52,53,54, 55,56,57,58,59};


boolean GeomObjectBSplineMeshInitCube ( GO_BSplineMesh *obj, boolean add )
{
  BSMvertex      *vert;
  BSMhalfedge    *halfe;
  BSMfacet       *fac;
  double         *vpc;
  int            *vhei, *fhei;
  byte           *mkcp, *mkhe, *mkfac;

  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return false;
  vert = malloc ( BSM_CUBE_NV*sizeof(BSMvertex) );
  halfe = malloc ( BSM_CUBE_NHE*sizeof(BSMhalfedge) );
  fac = malloc ( BSM_CUBE_NFAC*sizeof(BSMfacet) );
  vpc = malloc ( BSM_CUBE_NV*3*sizeof(double) );
  vhei = malloc ( BSM_CUBE_NHE*sizeof(int) );
  fhei = malloc ( BSM_CUBE_NHE*sizeof(int) );
  mkcp = malloc ( BSM_CUBE_NV );
  mkhe = malloc ( BSM_CUBE_NHE );
  mkfac = malloc ( BSM_CUBE_NFAC );
  if ( !vert || !halfe || !fac || !vpc || !vhei || !fhei ||
       !mkcp || !mkhe || !mkfac ) {
    if ( vert )  free ( vert );
    if ( halfe ) free ( halfe );
    if ( fac )   free ( fac );
    if ( vpc )   free ( vpc );
    if ( vhei )  free ( vhei );
    if ( fhei )  free ( fhei );
    if ( mkcp )  free ( mkcp );
    if ( mkhe )  free ( mkhe );
    if ( mkfac ) free ( mkfac );
    return false;
  }
  memcpy ( vert, cubevert, BSM_CUBE_NV*sizeof(BSMvertex) );
  memcpy ( halfe, cubehe, BSM_CUBE_NHE*sizeof(BSMhalfedge) );
  memcpy ( fac, cubefac, BSM_CUBE_NFAC*sizeof(BSMfacet) );
  memcpy ( vhei, cubevhei, BSM_CUBE_NHE*sizeof(int) );
  memcpy ( fhei, cubefachei, BSM_CUBE_NHE*sizeof(int) );
  memset ( mkcp, MASK_CP_MOVEABLE, BSM_CUBE_NV );
  memset ( mkhe, 0, BSM_CUBE_NHE );
  memset ( mkfac, 0, BSM_CUBE_NFAC );
  GeomObjectSetupIniPoints ( 3, false, &obj->me.cpdimen, BSM_CUBE_NV,
                             &cubevc[0].x, vpc );
  if ( add )
    GeomObjectExtendBSplineMesh ( obj, 3, false,
                      BSM_CUBE_NV, vert, vhei, vpc, BSM_CUBE_NHE, halfe,
                      BSM_CUBE_NFAC, fac, fhei, mkcp, mkhe, mkfac );
  else
    GeomObjectAssignBSplineMesh ( obj, 3, false,
                      BSM_CUBE_NV, vert, vhei, vpc, BSM_CUBE_NHE, halfe,
                      BSM_CUBE_NFAC, fac, fhei, mkcp, mkhe, mkfac );
  obj->me.dlistmask = 0;
  obj->spvlist_ok = false;
  return true;
} /*GeomObjectBSplineMeshInitCube*/

boolean GeomObjectBSplineMeshInitTetrahedron ( GO_BSplineMesh *obj, boolean add )
{
  BSMvertex      *vert;
  BSMhalfedge    *halfe;
  BSMfacet       *fac;
  double         *vpc;
  int            *vhei, *fhei;
  byte           *mkcp, *mkhe, *mkfac;

  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return false;
  vert = malloc ( BSM_TETRAHEDRON_NV*sizeof(BSMvertex) );
  halfe = malloc ( BSM_TETRAHEDRON_NHE*sizeof(BSMhalfedge) );
  fac = malloc ( BSM_TETRAHEDRON_NFAC*sizeof(BSMfacet) );
  vpc = malloc ( BSM_TETRAHEDRON_NV*3*sizeof(double) );
  vhei = malloc ( BSM_TETRAHEDRON_NHE*sizeof(int) );
  fhei = malloc ( BSM_TETRAHEDRON_NHE*sizeof(int) );
  mkcp = malloc ( BSM_TETRAHEDRON_NV );
  mkhe = malloc ( BSM_TETRAHEDRON_NHE );
  mkfac = malloc ( BSM_TETRAHEDRON_NFAC );
  if ( !vert || !halfe || !fac || !vpc || !vhei || !fhei ||
       !mkcp || !mkhe || !mkfac ) {
    if ( vert )  free ( vert );
    if ( halfe ) free ( halfe );
    if ( fac )   free ( fac );
    if ( vpc )   free ( vpc );
    if ( vhei )  free ( vhei );
    if ( fhei )  free ( fhei );
    if ( mkcp )  free ( mkcp );
    if ( mkhe )  free ( mkhe );
    if ( mkfac ) free ( mkfac );
    return false;
  }
  memcpy ( vert, tetrahvert, BSM_TETRAHEDRON_NV*sizeof(BSMvertex) );
  memcpy ( halfe, tetrahhe, BSM_TETRAHEDRON_NHE*sizeof(BSMhalfedge) );
  memcpy ( fac, tetrahfac, BSM_TETRAHEDRON_NFAC*sizeof(BSMfacet) );
  memcpy ( vhei, tetrahvhei, BSM_TETRAHEDRON_NHE*sizeof(int) );
  memcpy ( fhei, tetrahfhei, BSM_TETRAHEDRON_NHE*sizeof(int) );
  memset ( mkcp, MASK_CP_MOVEABLE, BSM_TETRAHEDRON_NV );
  memset ( mkhe, 0, BSM_TETRAHEDRON_NHE );
  memset ( mkfac, 0, BSM_TETRAHEDRON_NFAC );
  GeomObjectSetupIniPoints ( 3, false, &obj->me.cpdimen, BSM_TETRAHEDRON_NV,
                             &tetrahvc[0].x, vpc );
  if ( add )
    GeomObjectExtendBSplineMesh ( obj, 3, false,
                      BSM_TETRAHEDRON_NV, vert, vhei, vpc,
                      BSM_TETRAHEDRON_NHE, halfe,
                      BSM_TETRAHEDRON_NFAC, fac, fhei, mkcp, mkhe, mkfac );
  else
    GeomObjectAssignBSplineMesh ( obj, 3, false,
                      BSM_TETRAHEDRON_NV, vert, vhei, vpc,
                      BSM_TETRAHEDRON_NHE, halfe,
                      BSM_TETRAHEDRON_NFAC, fac, fhei, mkcp, mkhe, mkfac );
  obj->me.dlistmask = 0;  
  obj->spvlist_ok = false;
  return true;
} /*GeomObjectBSplineMeshInitTetrahedron*/

boolean GeomObjectBSplineMeshInitDodecahedron ( GO_BSplineMesh *obj, boolean add )
{
  BSMvertex      *vert;
  BSMhalfedge    *halfe;
  BSMfacet       *fac;
  double         *vpc;
  int            *vhei, *fhei;
  byte           *mkcp, *mkhe, *mkfac;

  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return false;
  vert = malloc ( BSM_DODECAHEDRON_NV*sizeof(BSMvertex) );
  halfe = malloc ( BSM_DODECAHEDRON_NHE*sizeof(BSMhalfedge) );
  fac = malloc ( BSM_DODECAHEDRON_NFAC*sizeof(BSMfacet) );
  vpc = malloc ( BSM_DODECAHEDRON_NV*3*sizeof(double) );
  vhei = malloc ( BSM_DODECAHEDRON_NHE*sizeof(int) );
  fhei = malloc ( BSM_DODECAHEDRON_NHE*sizeof(int) );
  mkcp = malloc ( BSM_DODECAHEDRON_NV );
  mkhe = malloc ( BSM_DODECAHEDRON_NHE );
  mkfac = malloc ( BSM_DODECAHEDRON_NFAC );
  if ( !vert || !halfe || !fac || !vpc || !vhei || !fhei ||
       !mkcp || !mkhe || !mkfac ) {
    if ( vert )  free ( vert );
    if ( halfe ) free ( halfe );
    if ( fac )   free ( fac );
    if ( vpc )   free ( vpc );
    if ( vhei )  free ( vhei );
    if ( fhei )  free ( fhei );
    if ( mkcp )  free ( mkcp );
    if ( mkhe )  free ( mkhe );
    if ( mkfac ) free ( mkfac );
    return false;
  }
  memcpy ( vert, dodecvert, BSM_DODECAHEDRON_NV*sizeof(BSMvertex) );
  memcpy ( halfe, dodeche, BSM_DODECAHEDRON_NHE*sizeof(BSMhalfedge) );
  memcpy ( fac, dodecfac, BSM_DODECAHEDRON_NFAC*sizeof(BSMfacet) );
  memcpy ( vhei, dodecvhei, BSM_DODECAHEDRON_NHE*sizeof(int) );
  memcpy ( fhei, dodecfhei, BSM_DODECAHEDRON_NHE*sizeof(int) );
  memset ( mkcp, MASK_CP_MOVEABLE, BSM_DODECAHEDRON_NV );
  memset ( mkhe, 0, BSM_DODECAHEDRON_NHE );
  memset ( mkfac, 0, BSM_DODECAHEDRON_NFAC );
  GeomObjectSetupIniPoints ( 3, false, &obj->me.cpdimen, BSM_DODECAHEDRON_NV,
                             &dodecvc[0].x, vpc );
  if ( add )
    GeomObjectExtendBSplineMesh ( obj, 3, false,
                      BSM_DODECAHEDRON_NV, vert, vhei, vpc,
                      BSM_DODECAHEDRON_NHE, halfe,
                      BSM_DODECAHEDRON_NFAC, fac, fhei, mkcp, mkhe, mkfac );
  else
    GeomObjectAssignBSplineMesh ( obj, 3, false,
                      BSM_DODECAHEDRON_NV, vert, vhei, vpc,
                      BSM_DODECAHEDRON_NHE, halfe,
                      BSM_DODECAHEDRON_NFAC, fac, fhei, mkcp, mkhe, mkfac );
  obj->me.dlistmask = 0;  
  obj->spvlist_ok = false;
  return true;
} /*GeomObjectBSplineMeshInitDodecahedron*/

boolean GeomObjectBSplineMeshInitKGon ( GO_BSplineMesh *obj, int k, boolean add )
{
  BSMvertex   *mv;
  BSMhalfedge *mhe;
  BSMfacet    *mfac;
  int         *mvhei, *mfhei;
  double      *mvpc;
  byte        *mkcp, *mkhe, *mkfac;
  int         i, cpdimen;
  double      phi;

  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return false;
  if ( k < 3 || k > MAX_BSM_DEGREE )
    return false;
  mv = malloc ( k*sizeof(BSMvertex) );
  mhe = malloc ( k*sizeof(BSMhalfedge) );
  mfac = malloc ( sizeof(BSMfacet) );
  mvhei = malloc ( k*sizeof(int) );
  mfhei = malloc ( k*sizeof(int) );
  mkcp = malloc ( k );
  mkhe = malloc ( k );
  mkfac = malloc ( 1 );
  cpdimen = obj->me.spdimen;
  mvpc = malloc ( k*cpdimen*sizeof(double) );
  if ( !mv || !mhe || !mfac || !mvpc || !mvhei || !mfhei ||
       !mkcp || !mkhe || !mkfac ) {
    if ( mv )    free ( mv );
    if ( mhe )   free ( mhe );
    if ( mfac )  free ( mfac );
    if ( mvpc )  free ( mvpc );
    if ( mvhei ) free ( mvhei );
    if ( mfhei ) free ( mfhei );
    if ( mkcp )  free ( mkcp );
    if ( mkhe )  free ( mkhe );
    if ( mkfac ) free ( mkfac );
    return false;
  }
  mfac[0].degree = k;
  mfac[0].firsthalfedge = 0;
  for ( i = 0; i < k; i++ ) {
    mfhei[i] = i;
    mvhei[i] = i;

    mhe[i].v0 = i;
    mhe[i].v1 = (i+1) % k;
    mhe[i].facetnum = 0;
    mhe[i].otherhalf = -1;

    mv[i].degree = 1;
    mv[i].firsthalfedge = i;

    phi = ((double)i-0.5)/(double)k*2.0*PI;
    mvpc[i*cpdimen]   = cos ( phi );
    mvpc[i*cpdimen+1] = sin ( phi );
    if ( cpdimen == 3 )
      mvpc[i*cpdimen+2] = 0.0;
  }
  memset ( mkcp, MASK_CP_MOVEABLE, k );
  memset ( mkhe, 0, k );
  mkfac[0] = 0;
  obj->rational = false;
  obj->me.cpdimen = cpdimen;
  if ( add )
    GeomObjectExtendBSplineMesh ( obj, cpdimen, false,
                        k, mv, mvhei, mvpc, k, mhe, 1, mfac, mfhei,
                        mkcp, mkhe, mkfac );
  else
    GeomObjectAssignBSplineMesh ( obj, cpdimen, false,
                        k, mv, mvhei, mvpc, k, mhe, 1, mfac, mfhei,
                        mkcp, mkhe, mkfac );
  obj->me.dlistmask = 0;
  obj->spvlist_ok = false;
  return true;
} /*GeomObjectBSplineMeshInitKGon*/

boolean GeomObjectBSplineMeshInitKPrism ( GO_BSplineMesh *obj, int k, boolean add )
{
  int         nv, nhe, nfac;
  BSMvertex   *mv;
  BSMhalfedge *mhe;
  BSMfacet    *mfac;
  int         *mvhei, *mfhei;
  double      *mvpc;
  byte        *mkcp, *mkhe, *mkfac;
  int         cpdimen;
  int         i, j, fhe;
  double      phi;

  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return false;
  if ( k < 3 || k > MAX_BSM_DEGREE )
    return false;
  nv   = 2*k;
  nhe  = 6*k;
  nfac = k+2;
  mv = malloc ( nv*sizeof(BSMvertex) );
  mhe = malloc ( nhe*sizeof(BSMhalfedge) );
  mfac = malloc ( nfac*sizeof(BSMfacet) );
  mvhei = malloc ( nhe*sizeof(int) );
  mfhei = malloc ( nhe*sizeof(int) );
  mkcp = malloc ( nv );
  mkhe = malloc ( nhe );
  mkfac = malloc ( nfac );
  cpdimen = obj->me.spdimen;
  mvpc = malloc ( nv*cpdimen*sizeof(double) );
  if ( !mv || !mhe || !mfac || !mvpc || !mvhei || !mfhei ||
       !mkcp || !mkhe || !mkfac ) {
    if ( mv )    free ( mv );
    if ( mhe )   free ( mhe );
    if ( mfac )  free ( mfac );
    if ( mvpc )  free ( mvpc );
    if ( mvhei ) free ( mvhei );
    if ( mfhei ) free ( mfhei );
    if ( mkcp )  free ( mkcp );
    if ( mkhe )  free ( mkhe );
    if ( mkfac ) free ( mkfac );
    return false;
  }
        /* setup the facets */
  mfac[0].degree = mfac[1].degree = k;
  mfac[0].firsthalfedge = 0;
  mfac[1].firsthalfedge = k;
  for ( i = 0; i < k; i++ ) {
    mfhei[i] = i;
    mfhei[k+i] = 4*k-1-i;
  }
  for ( i = 2; i < k+2; i++ ) {
    mfac[i].degree = 4;
    mfac[i].firsthalfedge = fhe = mfac[i-1].firsthalfedge + mfac[i-1].degree;
    mfhei[fhe] = k+i-2;
    mfhei[fhe+1] = 5*k+i-3;
    mfhei[fhe+2] = 2*k+i-2;
    mfhei[fhe+3] = 4*k+i-2;
  }
  mfhei[mfac[2].firsthalfedge+1] = 6*k-1;
        /* setup the vertices */
  for ( i = 0; i < k; i++ ) {
    mv[i].degree = 3;
    mv[i].firsthalfedge = fhe = 3*i;
    mvhei[fhe] = i;
    mvhei[fhe+1] = 5*k-1+i;
    mvhei[fhe+2] = k-1+i;

    mv[i+k].degree = 3;
    mv[i+k].firsthalfedge = fhe = 3*(i+k);
    mvhei[fhe] = 4*k-1+i;
    mvhei[fhe+1] = 2*k+i;
    mvhei[fhe+2] = 3*k-1+i;

    phi = ((double)i-0.5)/(double)k*2.0*PI;
    mvpc[i*cpdimen] = mvpc[(i+k)*cpdimen] =     cos ( phi );
    mvpc[i*cpdimen+1] = mvpc[(i+k)*cpdimen+1] = sin ( phi );
    if ( cpdimen == 3 ) {
      mvpc[i*cpdimen+2]   = -1.0;
      mvpc[(i+k)*cpdimen+2] =  1.0;
    }
  }
  fhe = mv[0].firsthalfedge;
  mvhei[fhe+1] = 6*k-1;
  mvhei[fhe+2] = 2*k-1;
  fhe = mv[k].firsthalfedge;
  mvhei[fhe] = 5*k-1;
  mvhei[fhe+2] = 4*k-1;
        /* setup the halfdedges */
  for ( i = 0; i < k; i++ ) {
    mhe[i].otherhalf = k+i;
    mhe[k+i].otherhalf = i;
    mhe[2*k+i].otherhalf = 3*k+i;
    mhe[3*k+i].otherhalf = 2*k+i;
    mhe[4*k+i].otherhalf = 5*k+i;
    mhe[5*k+i].otherhalf = 4*k+i;
  }
  for ( i = 0; i < nfac; i++ )
    for ( j = 0, fhe = mfac[i].firsthalfedge;  j < mfac[i].degree;  j++ )
      mhe[mfhei[fhe+j]].facetnum = i;
  for ( i = 0; i < nv; i++ )
    for ( j = 0, fhe = mv[i].firsthalfedge;  j < mv[i].degree;  j++ )
      mhe[mvhei[fhe+j]].v0 = i;
  for ( i = 0; i < nhe; i++ )
    mhe[i].v1 = mhe[mhe[i].otherhalf].v0;

  memset ( mkcp, MASK_CP_MOVEABLE, nv );
  memset ( mkhe, 0, nhe );
  memset ( mkfac, 0, nfac );
  obj->rational = false;
  obj->me.cpdimen = cpdimen;
  if ( add )
    GeomObjectExtendBSplineMesh ( obj, cpdimen, false,
                      nv, mv, mvhei, mvpc, nhe, mhe, nfac, mfac, mfhei,
                      mkcp, mkhe, mkfac );
  else
    GeomObjectAssignBSplineMesh ( obj, cpdimen, false,
                      nv, mv, mvhei, mvpc, nhe, mhe, nfac, mfac, mfhei,
                      mkcp, mkhe, mkfac );
  obj->me.dlistmask = 0;
  obj->spvlist_ok = false;
  return true;
} /*GeomObjectBSplineMeshInitKPrism*/

