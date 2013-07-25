
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "camera.h"
#include "raybezf.h"
#include "xgedit.h"

#define COPY_POSTSCRIPT
#ifdef COPY_POSTSCRIPT
#include "psout.h"
#endif

#include "render.h"


#define MIN_T 1.0e-5

typedef struct tnode {
    Box3f        bbox;
    struct tnode *left, *right;
    int          k;
  } tnode;

boolean RenderingIsOn = false;
XImage  *rendimage = NULL;

boolean swGaussian_c = false, swMean_c = false, swLambert_c = false,
        swReflection_c = false, swHighlight_c = false, swSections_c = false;
boolean swGaussian_d = false, swMean_d = false, swLambert_d = false,
        swReflection_d = false, swHighlight_d = false, swSections_d = false;

static vector3f lightdir[R_NLIGHTS];
static float    lightint[R_NLIGHTS+1];
static boolean Shadows, Antialias;

static char *imagedata = NULL;
static boolean RendererIsOk = false;

static CameraRecf _CPos;
static xge_widget *_xgw;

static int npatches;
static BezPatchTreefp *patch_tab = NULL;
static tnode *root = NULL;

static vector3f sectiondir = {0.0,0.0,2.0};
static point3f  rp0 = {0.0,0.0,-8.0};
static vector3f rv1 = {0.5,0.2,0.0};
static vector3f rv2 = {-0.4,1.0,0.0};
static float dfsf = 10.0;

static float minsf, maxsf;

static boolean (*dShapeFunc)( int n, int m, const point3f *cp, float u, float v,
                       point3f *p, vector3f *nv, vector3f *vv );
static float (*cShapeFunc)( int n, int m, const point3f *cp, float u, float v,
                     point3f *p, vector3f *nv, vector3f *vv );

static int y;
static byte *aabuf, *aaline[7];


/* ///////////////////////////////////////////////////////////////////////// */
/* discontinuous shape functions */
boolean dShapeFunc0 ( int n, int m, const point3f *cp, float u, float v,
                      point3f *p, vector3f *nv, vector3f *vv )
{
  return true;
} /*dShapeFunc0*/

boolean dHighlightLines ( int n, int m, const point3f *cp, float u, float v,
                           point3f *p, vector3f *nv, vector3f *vv )
{
  vector3f dp;
  double   s;

  SubtractPoints3f ( p, &rp0, &dp );
  s = dfsf*det3f ( nv, &dp, &rv2 ) / det3f ( nv, &rv1, &rv2 );
  if ( s < 0.0 )
    s = 1.0-s;
  return ((int)s) & 0x01;
} /*dHighlightLines*/

boolean dReflectionLines ( int n, int m, const point3f *cp, float u, float v,
                           point3f *p, vector3f *nv, vector3f *vv )
{
  vector3f r, dp;
  double   s;

  s = 2.0*DotProduct3f( nv, vv )/DotProduct3f ( nv, nv );
  SetVector3f ( &r, vv->x - s*nv->x, vv->y - s*nv->y, vv->z - s*nv->z );
  SubtractPoints3f ( p, &rp0, &dp );
  s = dfsf*det3f ( &r, &dp, &rv2 ) / det3f ( &r, &rv1, &rv2 );
  if ( s < 0.0 )
    s = 1.0-s;
  return !(((int)s) & 0x01);
} /*dReflectionLines*/

boolean dCrossSections ( int n, int m, const point3f *cp, float u, float v,
                         point3f *p, vector3f *nv, vector3f *vv )
{
  double s;

  s = dfsf*DotProduct3f ( p, &sectiondir );
  if ( s < 0.0 )
    s = 1.0-s;
  return !(((int)s) & 0x01);
} /*dCrossSections*/

boolean dLambertLight ( int n, int m, const point3f *cp, float u, float v,
                        point3f *p, vector3f *nv, vector3f *vv )
{
  double s;

  s = dfsf*DotProduct3f ( nv, &lightdir[0] );
  if ( s < 0.0 )
    s = 1.0-s;
  return !(((int)s) & 0x01);
} /*dLambertLight*/

boolean dMeanCurvature ( int n, int m, const point3f *cp, float u, float v,
                         point3f *p, vector3f *nv, vector3f *vv )
{
  float gauss, mean, s;

  mbs_GMCurvaturesBP3f ( n, m, cp, u, v, &gauss, &mean );
  s = dfsf*mean;
  if ( s < 0.0 )
    s = 1.0-s;
  return !(((int)s) & 0x01);
} /*dMeanCurvature*/

boolean dGaussianCurvature ( int n, int m, const point3f *cp, float u, float v,
                             point3f *p, vector3f *nv, vector3f *vv )
{
  float gauss, mean, s;

  mbs_GMCurvaturesBP3f ( n, m, cp, u, v, &gauss, &mean );
  s = dfsf*gauss;
  if ( s < 0.0 )
    s = 1.0-s;
  return !(((int)s) & 0x01);
} /*dMeanCurvature*/

/* ////////////////////////////////////////////////////////////////////////// */
/* continuous shape functions */
float cShapeFunc0 ( int n, int m, const point3f *cp, float u, float v,
                    point3f *p, vector3f *nv, vector3f *vv )
{
  return 0.0;
} /*cShapeFunc0*/

float cHighlightLines ( int n, int m, const point3f *cp, float u, float v,
                         point3f *p, vector3f *nv, vector3f *vv )
{
  vector3f dp;
  double   s;

  SubtractPoints3f ( p, &rp0, &dp );
  s = det3f ( nv, &dp, &rv2 ) / det3f ( nv, &rv1, &rv2 );
  return s;
} /*cHighlightLines*/

float cReflectionLines ( int n, int m, const point3f *cp, float u, float v,
                          point3f *p, vector3f *nv, vector3f *vv )
{
  vector3f r, dp;
  double   s;

  s = 2.0*DotProduct3f( nv, vv )/DotProduct3f ( nv, nv );
  SetVector3f ( &r, vv->x - s*nv->x, vv->y - s*nv->y, vv->z - s*nv->z );
  SubtractPoints3f ( p, &rp0, &dp );
  s = det3f ( &r, &dp, &rv2 ) / det3f ( &r, &rv1, &rv2 );
  return s;
} /*cReflectionLines*/

float cCrossSections ( int n, int m, const point3f *cp, float u, float v,
                       point3f *p, vector3f *nv, vector3f *vv )
{
  double s;

  s = DotProduct3f ( p, &sectiondir );
  return s;
} /*cCrossSections*/

float cLambertLight ( int n, int m, const point3f *cp, float u, float v,
                      point3f *p, vector3f *nv, vector3f *vv )
{
  double s;

  s = DotProduct3f ( nv, &lightdir[0] );
  return s;
} /*cLambertLight*/

float cMeanCurvature ( int n, int m, const point3f *cp, float u, float v,
                       point3f *p, vector3f *nv, vector3f *vv )
{
  float gauss, mean;

  mbs_GMCurvaturesBP3f ( n, m, cp, u, v, &gauss, &mean );
  return mean;
} /*cMeanCurvature*/

float cGaussianCurvature ( int n, int m, const point3f *cp, float u, float v,
                           point3f *p, vector3f *nv, vector3f *vv )
{
  float gauss, mean;

  mbs_GMCurvaturesBP3f ( n, m, cp, u, v, &gauss, &mean );
  return gauss;
} /*cMeanCurvature*/

/* ///////////////////////////////////////////////////////////////////////// */
void FindMinMaxShapeFunc ( void )
{
#define DENS 20
  int      i, j, k, m, n;
  float    u, v, f, fmin, fmax;
  point3f  *cp, p;
  vector3f nv, vv;
  boolean  ok;

  ok = false;
  fmin = 0.0;  fmax = 1.0;
  for ( i = 0; i < npatches; i++ ) {
    n = patch_tab[i]->n;
    m = patch_tab[i]->m;
    cp = patch_tab[i]->root->ctlpoints;
    if ( !ok ) {
      mbs_BCHornerNvP3f ( n, m, cp, 0.0, 0.0, &p, &nv );
      SubtractPoints3f ( &p, &_CPos.position, &vv );
      NormalizeVector3f ( &vv );
      fmin = fmax = cShapeFunc ( n, m, cp, 0.0, 0.0, &p, &nv, &vv );
      ok = true;
    }
    for ( j = 0; j <= DENS; j++ ) {
      u = (float)j/(float)DENS;
      for ( k = 0; k <= DENS; k++ ) {
        v = (float)k/(float)DENS;
        mbs_BCHornerNvP3f ( n, m, cp, u, v, &p, &nv );
        SubtractPoints3f ( &p, &_CPos.position, &vv ); 
        NormalizeVector3f ( &vv );
        f = cShapeFunc ( n, m, cp, u, v, &p, &nv, &vv );
        if ( f < fmin )
          fmin = f;
        else if ( f > fmax )
          fmax = f;
      }
    }
  }
  minsf = fmin;
  maxsf = fmax;
/*  printf ( "Shape function min = %f, max = %f\n", minsf, maxsf ); */
  if ( minsf >= maxsf ) {
    minsf -= 2.0;
    maxsf += 3.0;
  }
#undef DENS
} /*FindMinMaxShapeFunction*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean RendInit ( void )
{
  int nplanes;

  RendererIsOk = false;
  nplanes = XDisplayPlanes ( xgedisplay, xgescreen );
  rendimage = XCreateImage ( xgedisplay, xgevisual, nplanes, ZPixmap, 0,
                             NULL, xge_MAX_WIDTH, xge_MAX_HEIGHT, 8, 0 );
  imagedata = (char*)malloc ( xge_MAX_HEIGHT*rendimage->bytes_per_line );
  if ( rendimage && imagedata )
    rendimage->data = imagedata;
  else
    return false;
  patch_tab = malloc ( R_MAXPATCHES*sizeof(BezPatchTreefp) );
  if ( !patch_tab )
    return false;
  npatches = 0;
  cShapeFunc = cShapeFunc0;
  dShapeFunc = dShapeFunc0;
  if ( !(aabuf = malloc ( 7*3*(3*xge_MAX_WIDTH+4) )) )
    return false;
  RendererIsOk = true;
  return true;
} /*RendInit*/

void RendDestroy ( void )
{
  RendReset ();
  XDestroyImage ( rendimage );
  free ( patch_tab );
  free ( aabuf );
  RendererIsOk = false;
} /*RendDestroy*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean RendEnterBezPatchd ( int n, int m, const point3d *cp )
{
  void    *sp;
  int     size;
  point3f *cpf;

  sp = pkv_GetScratchMemTop ();
  if ( patch_tab && npatches < R_MAXPATCHES ) {
    size = (n+1)*(m+1);
    cpf = pkv_GetScratchMem ( size*sizeof(point3f) );
    if ( !cpf )
      goto failure;
    pkv_Selectdf ( 1, 3*size, 1, 1, (double*)cp, (float*)cpf );
    patch_tab[npatches] = rbez_NewBezPatchTreef ( 0, n, m, 0.0, 1.0, 0.0, 1.0, cpf );
    if ( patch_tab[npatches] )
      npatches ++;
    else
      goto failure;
    pkv_SetScratchMemTop ( sp );
    return true;
  }
failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*RendEnterBezPatchd*/

boolean RendEnterBSPatchd ( int n, int lknu, const double *knu,  
                            int m, int lknv, const double *knv,  
                            const point3d *cp )
{
  void   *sp;
  int    ku, kv, pitch, md, i, j;
  double *b, *c;

  sp = pkv_GetScratchMemTop ();
  if ( patch_tab ) {
    ku = mbs_NumKnotIntervalsd ( n, lknu, knu );
    kv = mbs_NumKnotIntervalsd ( m, lknv, knv );
    md = 3*(m+1);
    pitch = md*kv;
    b = pkv_GetScratchMemd ( pitch*(n+1)*ku );
    c = pkv_GetScratchMemd ( md*(n+1) );
    if ( b && c ) {
      mbs_BSPatchToBezd ( 3, n, lknu, knu, m, lknv, knv, 3*(lknv-m), (double*)cp,
                          &ku, NULL, NULL, &kv, NULL, NULL, pitch, b );
      for ( i = 0; i < ku; i++ )
        for ( j = 0; j < kv; j++ ) {
          pkv_Selectd ( n+1, md, pitch, md, &b[(n+1)*i*pitch+md*j], c );
          if ( !RendEnterBezPatchd ( n, m, (point3d*)c ) )
            goto failure;
        }
      pkv_SetScratchMemTop ( sp );
      return true;
    }
  }
failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*RendEnterBSPatchd*/

void FindSumBBox ( Box3f *box, Box3f *box1, Box3f *box2 )
{
  box->x0 = min ( box1->x0, box2->x0 );
  box->x1 = max ( box1->x1, box2->x1 );
  box->y0 = min ( box1->y0, box2->y0 );
  box->y1 = max ( box1->y1, box2->y1 );
  box->z0 = min ( box1->z0, box2->z0 );
  box->z1 = max ( box1->z1, box2->z1 );
} /*FindSumBBox*/

tnode* NewPatchNode ( int k )
{
  tnode *node;

  node = malloc ( sizeof(tnode) );
  if ( !node )
    exit ( 1 );
  node->left = node->right = NULL;
  node->k = k;
  memcpy ( &node->bbox, &patch_tab[k]->root->bbox, sizeof(Box3f) );
  return node;
} /*NewPatchNode*/

tnode* NewInnerNode ( tnode *left, tnode *right )
{
  tnode *node;

  node = malloc ( sizeof(tnode) );
  if ( !node )
    exit ( 1 );
  node->left = left;
  node->right = right;
  node->k = -1;
  FindSumBBox ( &node->bbox, &left->bbox, &right->bbox );
  return node;
} /*NewInnerNode*/

float BoxDiameter ( Box3f *box )
{
  float a, b, c;

  a = box->x1 - box->x0;
  b = box->y1 - box->y0;
  c = box->z1 - box->z0;
  return a*a + b*b + c*c;
} /*BoxDiameter*/

boolean CompNodePrio ( void *node1, void *node2 )
{
  float d1, d2;

  d1 = BoxDiameter ( &((tnode*)node1)->bbox );
  d2 = BoxDiameter ( &((tnode*)node2)->bbox );
  return d1 < d2;
} /*CompNodePrio*/

tnode* BuildPatchesTree ( int npatches )
{
  void  *sp;
  void  **pqueue;
  int   i, j, k, last;
  tnode *tn0, *tn1, *tn2;
  Box3f box;
  float d1, d2;

  sp = pkv_GetScratchMemTop ();
  pqueue = malloc ( npatches*sizeof(void*) );
  if ( !pqueue )
    exit ( 1 );
  for ( i = 0; i < npatches; i++ )
    pqueue[i] = NewPatchNode ( i );
  pkv_HeapOrder ( pqueue, npatches, CompNodePrio );
  for ( i = npatches-1; i > 0; i-- ) {
    tn0 = pqueue[0];
    last = i;
    pkv_HeapRemove ( pqueue, &last, 0, CompNodePrio );
    tn1 = pqueue[0];
    FindSumBBox ( &box, &tn0->bbox, &tn1->bbox );
    d1 = BoxDiameter ( &box );
    k = 0;
    for ( j = 1; j <= last; j++ ) {
      tn2 = pqueue[j];
      FindSumBBox ( &box, &tn0->bbox, &tn2->bbox );
      d2 = BoxDiameter ( &box );
      if ( d2 < d1 ) {
        k = j;
        d1 = d2;
        tn1 = tn2;
      }
    }
    pkv_HeapRemove ( pqueue, &last, k, CompNodePrio );
    tn2 = NewInnerNode ( tn0, tn1 );
    pkv_HeapInsert ( pqueue, &last, (void*)tn2, CompNodePrio );
  }
  tn0 = pqueue[0];
  pkv_SetScratchMemTop ( sp );
  return tn0;
} /*BuildPatchesTree*/

void DestroyTNodeTree ( tnode *node )
{
  if ( node ) {
    DestroyTNodeTree ( node->left );
    DestroyTNodeTree ( node->right );
    free ( node );
  }
} /*DestroyTNodeTree*/

#ifdef COPY_POSTSCRIPT
static void DumpPostscriptPicture ( void )
{
  void *sp;
  int  x, y;
  byte *buf;
  xgecolour_int pixel;

  sp = pkv_GetScratchMemTop ();
  if ( (buf = pkv_GetScratchMem ( 3*_xgw->w )) ) {
    ps_OpenFile ( "picture.ps", 600 );
    ps_Init_BitmapRGBP ( _xgw->w, _xgw->h, 0, 0 );
    for ( y = 0; y < _xgw->h; y++ ) {
      for ( x = 0; x < _xgw->w; x++ ) {
        pixel = XGetPixel ( rendimage, _xgw->x+x, _xgw->y+y );
        xge_GetPixelColour ( pixel, &buf[3*x], &buf[3*x+1], &buf[3*x+2] );
      }
      ps_Out_LineRGBP ( buf );
    }
    ps_CloseFile ();
  }
  pkv_SetScratchMemTop ( sp );
} /*DumpPostscriptPicture*/
#endif

boolean RendReset ( void )
{
  int i;

  if ( RenderingIsOn ) {
    RenderingIsOn = false;
    DestroyTNodeTree ( root );
    root = NULL;
    for ( i = 0; i < npatches; i++ ) {
      rbez_DestroyBezPatchTreef ( patch_tab[i] );
      patch_tab[i] = NULL;
    }
#ifdef COPY_POSTSCRIPT
    DumpPostscriptPicture ();
#endif
  }
  npatches = 0;
  return true;
} /*RendReset*/

static void Point3df ( const point3d *p, point3f *q )
{
  SetPoint3f ( q, (float)p->x, (float)p->y, (float)p->z );
} /*Point3df*/

boolean RendEnterCamerad ( CameraRecd *CPos, xge_widget *er )
{
  if ( CPos->parallel )
    return false;
      /* convert from double to single precision */
  _CPos.parallel = false;
  _CPos.upside   = CPos->upside;
  _CPos.c_fixed  = CPos->c_fixed;
  _CPos.magnification = CPos->magnification;
  _CPos.xmin   = CPos->xmin;
  _CPos.ymin   = CPos->ymin;
  _CPos.width  = CPos->width;
  _CPos.height = CPos->height;
  _CPos.aspect = (float)CPos->aspect;
  Point3df ( &CPos->g_centre, &_CPos.g_centre );
  Point3df ( &CPos->c_centre, &_CPos.c_centre );
  Point3df ( &CPos->position, &_CPos.position );
  _CPos.psi   = (float)CPos->psi;
  _CPos.theta = (float)CPos->theta;
  _CPos.phi   = (float)CPos->phi;
  _CPos.vd.persp.f     = (float)CPos->vd.persp.f;
  _CPos.vd.persp.xi0   = (float)CPos->vd.persp.xi0;
  _CPos.vd.persp.eta0  = (float)CPos->vd.persp.eta0;
  _CPos.vd.persp.dxi0  = (float)CPos->vd.persp.dxi0;
  _CPos.vd.persp.deta0 = (float)CPos->vd.persp.deta0;
  CameraSetMappingf ( &_CPos );
  _xgw = er;
  return true;
} /*RendEnterCamerad*/

void RendEnterLightsd ( int nlights, const vector3d *light_dir,
                        const double *light_int )
{
  int   i;
  float lint;

  for ( i = 0; i < nlights && i < R_NLIGHTS; i++ ) {
    SetVector3f ( &lightdir[i], light_dir[i].x, light_dir[i].y, light_dir[i].z );
    NormalizeVector3f ( &lightdir[i] );
    lint = min ( 1.0, light_int[i] );
    lightint[i] = max ( 0.0, lint );
  }
  lightint[R_NLIGHTS] = light_int[nlights];  /* ambient light intensity */
} /*RendEnterLightsd*/

boolean RendBegin ( boolean swShadows, boolean swAntialias )
{
  y = _xgw->y;

        /* build the scene hierarchy binary tree */
  root = BuildPatchesTree ( npatches );
        /* setup the shape function texturing procedures */
  if ( swGaussian_c )
    cShapeFunc = cGaussianCurvature;
  else if ( swMean_c )
    cShapeFunc = cMeanCurvature;
  else if ( swLambert_c )
    cShapeFunc = cLambertLight;
  else if ( swReflection_c )
    cShapeFunc = cReflectionLines;
  else if ( swHighlight_c )
    cShapeFunc = cHighlightLines;
  else if ( swSections_c )
    cShapeFunc = cCrossSections;
  else
    cShapeFunc = cShapeFunc0;
  if ( swGaussian_d )
    dShapeFunc = dGaussianCurvature;
  else if ( swMean_d )
    dShapeFunc = dMeanCurvature;
  else if ( swLambert_d )
    dShapeFunc = dLambertLight;
  else if ( swReflection_d )
    dShapeFunc = dReflectionLines;
  else if ( swHighlight_d )
    dShapeFunc = dHighlightLines;
  else if ( swSections_d )
    dShapeFunc = dCrossSections;
  else
    dShapeFunc = dShapeFunc0;
  FindMinMaxShapeFunc ();
  Shadows = swShadows;
  Antialias = swAntialias;
  RenderingIsOn = true;
  if ( Antialias )
    InitRenderingAA ();
  return true;
} /*RendBegin*/

/* ///////////////////////////////////////////////////////////////////////// */
static char ClipTestf ( float p, float q, float *t0, float *t1 )
{
#define EPS 5.0e-6
  float r;

  if ( p < 0.0 ) {
    r = q/p;
    if ( r > *t1+EPS ) return 0;
    else if ( r > *t0 ) *t0 = r;
  }
  else if ( p > 0.0 ) {
    r = q/p;
    if ( r < *t0-EPS ) return 0;
    else if ( r < *t1 ) *t1 = r;
  }
  else if ( q < -EPS ) return 0;
  return 1;
#undef EPS 
} /*ClipTestf*/

static char TestRayBBoxf ( ray3f *ray, Box3f *box )
{
  float t0, t1;

  t0 = 0.0;  t1 = 1.0e38;
  if ( ClipTestf ( -ray->v.x, ray->p.x - box->x0, &t0, &t1 ) )
    if ( ClipTestf (  ray->v.x, box->x1 - ray->p.x, &t0, &t1 ) )
      if ( ClipTestf ( -ray->v.y, ray->p.y - box->y0, &t0, &t1 ) )
        if ( ClipTestf (  ray->v.y, box->y1 - ray->p.y, &t0, &t1 ) )
          if ( ClipTestf ( -ray->v.z, ray->p.z - box->z0, &t0, &t1 ) )
            if ( ClipTestf (  ray->v.z, box->z1 - ray->p.z, &t0, &t1 ) )
              return 1;
  return 0;
} /*TestRayBBoxf*/

void rFindRayTrInters ( ray3f *ray, tnode *node,
                        int *ninters, RayObjectIntersf *inters, int *k )
{
  int _k, i, nint;

  if ( (_k = node->k) >= 0 ) {  /* it is a leaf, i.e. a Bezier patch */
           /* no ray/box intersection is tested here, as the ray/patch */
           /* intersection procedure has its own test */
    if ( rbez_FindRayBezPatchIntersf (
           patch_tab[_k], ray, MAXLEVEL, MAXINTERS, &nint, &inters[*ninters] ) ) {
           /* ignore all intersections too close to the ray origin */
      for ( i = 0; i < nint; )
        if ( inters[i].t < MIN_T )
          inters[i] = inters[--nint];
        else
          i ++;
      if ( nint ) {
           /* something was left */
        if ( !(*ninters) )
          *k = _k;
        nint += *ninters;
        for ( i = 1; i < nint; i++ )
          if ( inters[i].t < inters[0].t ) {
            inters[0] = inters[i];
            *k = _k;
          }
        *ninters = 1;
      }
    }
  }
  else if ( TestRayBBoxf ( ray, &node->bbox ) ) {
    rFindRayTrInters ( ray, node->left, ninters, inters, k );
    rFindRayTrInters ( ray, node->right, ninters, inters, k );
  }
} /*rFindRayTrInters*/

boolean FindRayTrInters ( ray3f *ray, RayObjectIntersf *inters, int *k )
{
  int              ninters;
  RayObjectIntersf _inters[MAXINTERS+1];

  *k = -1;
  ninters = 0;
  rFindRayTrInters ( ray, root, &ninters, _inters, k );
  if ( ninters ) {
    *inters = _inters[0];
    return true;
  }
  else
    return false;
} /*FindRayTrInters*/

void GetTexColour ( int n, int m, const point3f *cp, float u, float v,
                    point3f *p, vector3f *nv, vector3f *vv,
                    vector3f *texcolour )
{
  float t;

  t = cShapeFunc ( n, m, cp, u, v, p, nv, vv );
  t = (t-minsf)/(maxsf-minsf);
  t = max ( 0.0, t );
  t = min ( 1.0, t );
  /* "rainbow" texture */
  if      ( t < 0.2 ) SetPoint3f ( texcolour, 1.0, 0.0, 1.0-5.0*t );  
  else if ( t < 0.4 ) SetPoint3f ( texcolour, 1.0, 5.0*(t-0.2), 0.0 );
  else if ( t < 0.6 ) SetPoint3f ( texcolour, 1.0-5.0*(t-0.4), 1.0, 0.0 );
  else if ( t < 0.8 ) SetPoint3f ( texcolour, 0.0, 1.0, 5.0*(t-0.6) );
  else                SetPoint3f ( texcolour, 0.0, 1.0-5.0*(t-0.8), 1.0 );
  if ( !dShapeFunc ( n, m, cp, u, v, p, nv, vv ) )
    MultVector3f ( 0.85, texcolour, texcolour );
} /*GetTexColour*/

void AmbientLight ( vector3f *tex, float amb_intens, vector3f *rgb )
{
  MultVector3f ( amb_intens, tex, rgb );
} /*AmbientLight*/

boolean IsInShadow ( point3f *p, vector3f *light )
{
  ray3f ray;
  int   k;
  RayObjectIntersf inters;

  ray.p = *p;
  ray.v = *light;
  if ( FindRayTrInters ( &ray, &inters, &k ) )
    return true;
  else
    return false;
} /*IsInShadow*/

void PhongLight ( vector3f *tex, point3f *p, vector3f *nv,
                  vector3f *vv, vector3f *light, float light_intens,
                  vector3f *rgb )
{
  float    c, d;
  vector3f h;

  NormalizeVector3f ( nv );
  c = DotProduct3f ( nv, light );
  d = DotProduct3f ( nv, vv );
  if ( c*d <= 0.0 ) {  
    if ( Shadows ) {
      if ( IsInShadow ( p, light ) )
        return;
    }
    SubtractPoints3f ( vv, light, &h );
    NormalizeVector3f ( &h );
    d = fabs ( DotProduct3f ( nv, &h ) );
    d = d*d;  d = d*d;  d = d*d;  d = d*d;  d = d*d;  d = d*d;
    c = 0.8*light_intens*fabs(c);
    d *= 0.15*light_intens;
    rgb->x += (c*tex->x + d);
    rgb->y += (c*tex->y + d);
    rgb->z += (c*tex->z + d);
  }
} /*PhongLight*/

void GetPixelColour ( ray3f *ray, unsigned int *r, unsigned int *g,
                      unsigned int *b )
{
  int              k;
  RayObjectIntersf iint;
  vector3f         texcolour, rgb;
  unsigned int     a;

  *r = *g = *b = 220;
  if ( FindRayTrInters ( ray, &iint, &k ) ) {
    GetTexColour ( patch_tab[k]->n, patch_tab[k]->m, patch_tab[k]->root->ctlpoints,
                   iint.u, iint.v, &iint.p, &iint.nv, &ray->v, &texcolour );
    AmbientLight ( &texcolour, lightint[R_NLIGHTS], &rgb );
    for ( k = 0; k < R_NLIGHTS; k++ )
      if ( lightint[k] > 0.001 )
        PhongLight ( &texcolour, &iint.p, &iint.nv, &ray->v,
                     &lightdir[k], lightint[k], &rgb );
    a = (unsigned int)(rgb.x*255.0);  *r = min ( a, 255 );
    a = (unsigned int)(rgb.y*255.0);  *g = min ( a, 255 );
    a = (unsigned int)(rgb.z*255.0);  *b = min ( a, 255 );
  }
} /*GetPixelColour*/

int RenderLineA ( void )
{
  int   x;
  ray3f ray;
  unsigned int r, g, b;

  if ( !RenderingIsOn )
    return y;
  for ( x = _xgw->x; x < _xgw->x+_xgw->w; x++ ) {
        /* get the ray */
    CameraRayOfPixelf ( &_CPos, (float)x, (float)y, &ray );
        /* compute the pixel colour */
    GetPixelColour ( &ray, &r, &g, &b );
    XPutPixel ( rendimage, x, y, xge_PixelColour ( r, g, b ) );
  }
  y++;
  if ( y >= _xgw->y+_xgw->h )
    RendReset ();
  return y;
} /*RenderLineA*/

boolean OverThreshold1 ( byte *a, byte *b )
{
#define THR 5
  int  i;

  for ( i = 0; i < 3; i++ ) {
    if ( abs ( a[i]-b[i] ) > THR )
      return true;
  }
  return false;
#undef THR
} /*OverThreshold1*/

boolean OverThreshold2 ( byte *a, byte *b, byte *c, byte *d )
{
#define THR 5
  byte cmin, cmax;
  int  i;

  for ( i = 0; i < 3; i++ ) {
    if ( a[i] < b[i] ) { cmin = a[i];  cmax = b[i]; }
                  else { cmin = b[i];  cmax = a[i]; }
    if ( c[i] < d[i] ) { cmin = min ( cmin, c[i] );  cmax = max ( cmax, d[i] ); }
                  else { cmin = min ( cmin, d[i] );  cmax = max ( cmax, c[i] ); }
    if ( cmax-cmin > THR )
      return true;
  }
  return false;
#undef THR
} /*OverThreshold2*/

void Interpol2_3 ( byte *a, byte *b, byte *c )
{
  int i;

  for ( i = 0; i < 3; i++ )
    c[i] = (2*a[i]+b[i])/3;
} /*Interpol2_3*/

void GetSample ( float x, float y, byte *rgb )
{
  ray3f ray;
  unsigned int r, g, b;

  CameraRayOfPixelf ( &_CPos, (float)x, y, &ray );
  GetPixelColour ( &ray, &r, &g, &b );
  rgb[0] = (byte)r;
  rgb[1] = (byte)g;
  rgb[2] = (byte)b;
} /*GetSample*/

void RenderSubpixelLine ( float y, byte *aaline )
{
  int          x, i;

  if ( !RenderingIsOn )
    return;
  for ( i = 0, x = _xgw->x-1;  x <= _xgw->x+_xgw->w;  x++, i += 9 )
    GetSample ( (float)x, y, &aaline[i] );
  for ( i = 0, x = _xgw->x-1; x < _xgw->x+_xgw->w; x++, i += 9 ) {
    if ( OverThreshold1 ( &aaline[i], &aaline[i+9] ) ) {
      GetSample ( (float)x+(float)(1.0/3.0), y, &aaline[i+3] );
      GetSample ( (float)x+(float)(2.0/3.0), y, &aaline[i+6] );
    }
    else {
      Interpol2_3 ( &aaline[i], &aaline[i+9], &aaline[i+3] );
      Interpol2_3 ( &aaline[i+9], &aaline[i], &aaline[i+6] );
    }
  }
} /*RenderSubpixelLine*/

void GetSupersamples ( float y )
{
  int i, x;

  for ( i = 0, x = _xgw->x-1;  x <= _xgw->x+_xgw->w;  x++, i += 9 ) {
    if ( OverThreshold1 ( &aaline[3][i], &aaline[6][i] ) ) {
      GetSample ( (float)x, (float)y+(float)(1.0/3.0), &aaline[4][i] );
      GetSample ( (float)x, (float)y+(float)(2.0/3.0), &aaline[5][i] );
    }
    else {
      Interpol2_3 ( &aaline[3][i], &aaline[6][i], &aaline[4][i] );
      Interpol2_3 ( &aaline[6][i], &aaline[3][i], &aaline[5][i] );
    }
  }
  for ( i = 0, x = _xgw->x-1;  x < _xgw->x+_xgw->w;  x++, i += 9 ) {
    if ( OverThreshold2 ( &aaline[3][i], &aaline[3][i+9],
                          &aaline[6][i], &aaline[6][i+9] ) ) {
      GetSample ( (float)x+(float)(1.0/3.0), y+(float)(1.0/3.0), &aaline[4][i+3] );
      GetSample ( (float)x+(float)(2.0/3.0), y+(float)(1.0/3.0), &aaline[4][i+6] );
      GetSample ( (float)x+(float)(1.0/3.0), y+(float)(2.0/3.0), &aaline[5][i+3] );
      GetSample ( (float)x+(float)(2.0/3.0), y+(float)(2.0/3.0), &aaline[5][i+6] );
    }
    else {
      Interpol2_3 ( &aaline[4][i], &aaline[4][i+9], &aaline[4][i+3] );
      Interpol2_3 ( &aaline[4][i+0], &aaline[4][i], &aaline[4][i+6] );
      Interpol2_3 ( &aaline[5][i], &aaline[5][i+9], &aaline[5][i+3] );
      Interpol2_3 ( &aaline[5][i+0], &aaline[5][i], &aaline[5][i+6] );
    }
  }
} /*GetSupersamples*/

void InitRenderingAA ( void )
{
  int i, lgt;

  lgt = 3*(3*_CPos.width+4);
  for ( i = 0; i < 7; i++ )
    aaline[i] = &aabuf[i*lgt];
  RenderSubpixelLine ( (float)(y-1), aaline[3] );
  RenderSubpixelLine ( (float)y, aaline[6] );
  GetSupersamples ( (float)(y-1) );
} /*InitRenderingAA*/

int RenderLineAA ( void )
{
  int          i, j, k, l, x;
  byte         *a;
  unsigned int b[8][3];
      /* Gaussian filter coefficients for convolution */
  const static unsigned int fc[7] = {4,22,60,84,60,22,4};
/*  const static unsigned int fc[7] = {1,15,62,100,62,15,1};*/

  if ( RenderingIsOn ) {
        /* move previously computed samples up */
    a = aaline[0];  aaline[0] = aaline[3];  aaline[3] = aaline[6];
    aaline[6] = aaline[2];  aaline[2] = aaline[5];  aaline[5] = aaline[1];
    aaline[1] = aaline[4];  aaline[4] = a;
        /* get the samples for the current and next line */
    RenderSubpixelLine ( (float)(y+1), aaline[6] );
    GetSupersamples ( (float)y );
        /* filtering and writing the pixels to the image */
    for ( i = 0, x = _xgw->x;  x <= _xgw->x+_xgw->w;  x++, i += 9 ) {
      memset ( b, 0, 24*sizeof(unsigned int) );
      for ( j = 0; j < 7; j++ )
        for ( k = 0; k < 7; k++ )
          for ( l = 0; l < 3; l++ )
            b[j][l] += fc[k]*aaline[j][i+3*k+l];
      for ( j = 0; j < 7; j++ )
        for ( l = 0; l < 3; l++ )
          b[7][l] += fc[j]*b[j][l];
      XPutPixel ( rendimage, x, y,
                  xge_PixelColour ( b[7][0] >> 16, b[7][1] >> 16, b[7][2] >> 16 ) );
    }
    y ++;
    if ( y >= _xgw->y+_xgw->h )
      RendReset ();
  }
  return y;
} /*RenderLineAA*/

int RenderLine ( void )
{
  if ( Antialias )
    return RenderLineAA ();
  else
    return RenderLineA ();
} /*RenderLine*/

