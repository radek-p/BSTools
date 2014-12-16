
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
#include "pkvthreads.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "camera.h"
#include "raybez.h"
#include "xgedit.h"

#define _COPY_POSTSCRIPT
#ifdef COPY_POSTSCRIPT
#include "psout.h"
#endif

#include "pozwalaj0.h"
#include "render.h"

#define MIN_T 1.0e-5

/* object types, more to be added */
#define obj_TRIANGLE  1
#define obj_BEZPATCH  2
#define obj_RBEZPATCH 3
#define obj_BEZCURVE  4
#define obj_RBEZCURVE 5

typedef struct {
    Box3d           bbox;
    point3d         p0;
    vector3d        a1, a2, n;
  } triangle_data;

typedef struct {
    char            type;  /* == obj_TRIANGLE */
    double          colour[3];
    triangle_data   *trdata;
  } renderobj_triangle;

typedef struct {
    char            type;  /* == obj_BEZPATCH */
    double          colour[3];
    BezPatchTreedp  ptree;
  } renderobj_bezpatch;

typedef struct {
    char            type;  /* == obj_RBEZPATCH */
    double          colour[3];
    RBezPatchTreedp ptree;
  } renderobj_rbezpatch;

typedef struct {
    char            type;  /* == obj_BEZCURVE */
    double          colour[3];
    BezCurveTreedp  ctree;
  } renderobj_bezcurve;

typedef struct {
    char            type;  /* == obj_RBEZCURVE */
    double          colour[3];
    RBezCurveTreedp ctree;
  } renderobj_rbezcurve;

typedef union {
    char type;
    renderobj_triangle   triang;
    renderobj_bezpatch   bezp;
    renderobj_rbezpatch  rbezp;
    renderobj_bezcurve   bezc;
    renderobj_rbezcurve  rbezc;
  } renderobj;

/* binary tree of object hierarchy */
typedef struct tnode {
    Box3d        bbox;
    struct tnode *left, *right;
    int          k;
  } tnode;

/* this is to be replaced by a separate variable, dedicated for the renderer */
extern int rendering_npthreads;

static int renderer_npthreads = 1;
static clock_t tic;

boolean RenderingIsOn = false;
XImage  *rendimage = NULL;

boolean swGaussian_c = false, swMean_c = false, swLambert_c = false,
        swReflection_c = false, swHighlight_c = false, swSections_c = false;
boolean swGaussian_d = false, swMean_d = false, swLambert_d = false,
        swReflection_d = false, swHighlight_d = false, swSections_d = false;

static vector3d lightdir[R_NLIGHTS];
static double   lightint[R_NLIGHTS+1];
boolean swShadows = false, swAntialias = true;

static char *imagedata = NULL;
static boolean RendererIsOk = false;

static CameraRecd _CPos;
static xge_widget *_xgw;

static int        obj_tab_length;
static int        nobjects;
static renderobj  *obj_tab = NULL;
static tnode      *root = NULL;

static vector3d sectiondir = {0.0,0.0,2.0};
static point3d  rp0 = {0.0,0.0,-8.0};
static vector3d rv1 = {0.5,0.2,0.0};
static vector3d rv2 = {-0.4,1.0,0.0};
static point3d  hp0 = {0.0,0.0,-8.0};
static vector3d hv1 = {0.5,0.2,0.0};
static vector3d hv2 = {-0.4,1.0,0.0};
double          render_dfsf = 0.5;
static double   dfsf = 10.0;
double          render_cfrange[2] = {0.0,1.0};

static double minsf, maxsf; /* extremal values of the shape function */

static boolean (*dShapeFunc)( int n, int m, const point3d *cp, double u, double v,
                       point3d *p, vector3d *nv, vector3d *vv );
static double (*cShapeFunc)( int n, int m, const point3d *cp, double u, double v,
                       point3d *p, vector3d *nv, vector3d *vv );
static boolean (*dRShapeFunc)( int n, int m, const point4d *cp, double u, double v,
                       point3d *p, vector3d *nv, vector3d *vv );
static double (*cRShapeFunc)( int n, int m, const point4d *cp, double u, double v,
                     point3d *p, vector3d *nv, vector3d *vv );
static void (*GetTexColour) ( char type, double *colour,
                     int n, int m, const point4d *cp, double u, double v,
                     point3d *p, vector3d *nv, vector3d *vv,
                     vector3d *texcolour );

static int y;
static byte *aabuf, *aaline[7];


/* ///////////////////////////////////////////////////////////////////////// */
/* discontinuous shape functions */
boolean dShapeFunc0 ( int n, int m, const point3d *cp, double u, double v,
                      point3d *p, vector3d *nv, vector3d *vv )
{
  return true;
} /*dShapeFunc0*/

boolean dHighlightLines ( int n, int m, const point3d *cp, double u, double v,
                           point3d *p, vector3d *nv, vector3d *vv )
{
  vector3d dp;
  double   s;

  SubtractPoints3d ( p, &hp0, &dp );
  s = dfsf*det3d ( nv, &dp, &hv2 ) / det3d ( nv, &hv1, &hv2 );
  if ( s < 0.0 )
    s = 1.0-s;
  return ((int)s) & 0x01;
} /*dHighlightLines*/

boolean dReflectionLines ( int n, int m, const point3d *cp, double u, double v,
                           point3d *p, vector3d *nv, vector3d *vv )
{
  vector3d r, dp;
  double   s;

  s = 2.0*DotProduct3d( nv, vv )/DotProduct3d ( nv, nv );
  SetVector3d ( &r, vv->x - s*nv->x, vv->y - s*nv->y, vv->z - s*nv->z );
  SubtractPoints3d ( p, &rp0, &dp );
  s = dfsf*det3d ( &r, &dp, &rv2 ) / det3d ( &r, &rv1, &rv2 );
  if ( s < 0.0 )
    s = 1.0-s;
  return !(((int)s) & 0x01);
} /*dReflectionLines*/

boolean dCrossSections ( int n, int m, const point3d *cp, double u, double v,
                         point3d *p, vector3d *nv, vector3d *vv )
{
  double s;

  s = dfsf*DotProduct3d ( p, &sectiondir );
  if ( s < 0.0 )
    s = 1.0-s;
  return !(((int)s) & 0x01);
} /*dCrossSections*/

boolean dLambertLight ( int n, int m, const point3d *cp, double u, double v,
                        point3d *p, vector3d *nv, vector3d *vv )
{
  double s;

  s = dfsf*DotProduct3d ( nv, &lightdir[0] );
  if ( s < 0.0 )
    s = 1.0-s;
  return !(((int)s) & 0x01);
} /*dLambertLight*/

boolean dMeanCurvature ( int n, int m, const point3d *cp, double u, double v,
                          point3d *p, vector3d *nv, vector3d *vv )
{
  double gauss, mean, s;

  mbs_GMCurvaturesBP3d ( n, m, cp, u, v, &gauss, &mean );
  s = dfsf*mean;
  if ( s < 0.0 )
    s = 1.0-s;
  return !(((int)s) & 0x01);
} /*dMeanCurvature*/

boolean dGaussianCurvature ( int n, int m, const point3d *cp, double u, double v,
                             point3d *p, vector3d *nv, vector3d *vv )
{
  double gauss, mean, s;

  mbs_GMCurvaturesBP3d ( n, m, cp, u, v, &gauss, &mean );
  s = dfsf*gauss;
  if ( s < 0.0 )
    s = 1.0-s;
  return !(((int)s) & 0x01);
} /*dMeanCurvature*/

boolean dRShapeFunc0 ( int n, int m, const point4d *cp, double u, double v,
                       point3d *p, vector3d *nv, vector3d *vv )
{
  return true;
} /*dRShapeFunc0*/

boolean dRHighlightLines ( int n, int m, const point4d *cp, double u, double v,
                            point3d *p, vector3d *nv, vector3d *vv )
{
  vector3d dp;
  double   s;

  SubtractPoints3d ( p, &hp0, &dp );
  s = dfsf*det3d ( nv, &dp, &hv2 ) / det3d ( nv, &hv1, &hv2 );
  if ( s < 0.0 )
    s = 1.0-s;
  return ((int)s) & 0x01;
} /*dRHighlightLines*/

boolean dRReflectionLines ( int n, int m, const point4d *cp, double u, double v,
                            point3d *p, vector3d *nv, vector3d *vv )
{
  vector3d r, dp;
  double   s;

  s = 2.0*DotProduct3d( nv, vv )/DotProduct3d ( nv, nv );
  SetVector3d ( &r, vv->x - s*nv->x, vv->y - s*nv->y, vv->z - s*nv->z );
  SubtractPoints3d ( p, &rp0, &dp );
  s = dfsf*det3d ( &r, &dp, &rv2 ) / det3d ( &r, &rv1, &rv2 );
  if ( s < 0.0 )
    s = 1.0-s;
  return !(((int)s) & 0x01);
} /*dRReflectionLines*/

boolean dRCrossSections ( int n, int m, const point4d *cp, double u, double v,
                          point3d *p, vector3d *nv, vector3d *vv )
{
  double s;

  s = dfsf*DotProduct3d ( p, &sectiondir );
  if ( s < 0.0 )
    s = 1.0-s;
  return !(((int)s) & 0x01);
} /*dRCrossSections*/

boolean dRLambertLight ( int n, int m, const point4d *cp, double u, double v,
                         point3d *p, vector3d *nv, vector3d *vv )
{
  double s;

  s = dfsf*DotProduct3d ( nv, &lightdir[0] );
  if ( s < 0.0 )
    s = 1.0-s;
  return !(((int)s) & 0x01);
} /*dRLambertLight*/

boolean dRMeanCurvature ( int n, int m, const point4d *cp, double u, double v,
                          point3d *p, vector3d *nv, vector3d *vv )
{
  double gauss, mean, s;

  mbs_GMCurvaturesBP3Rd ( n, m, cp, u, v, &gauss, &mean );
  s = dfsf*mean;
  if ( s < 0.0 )
    s = 1.0-s;
  return !(((int)s) & 0x01);
} /*dRMeanCurvature*/

boolean dRGaussianCurvature ( int n, int m, const point4d *cp, double u, double v,
                              point3d *p, vector3d *nv, vector3d *vv )
{
  double gauss, mean, s;

  mbs_GMCurvaturesBP3Rd ( n, m, cp, u, v, &gauss, &mean );
  s = dfsf*gauss;
  if ( s < 0.0 )
    s = 1.0-s;
  return !(((int)s) & 0x01);
} /*dRMeanCurvature*/

/* ////////////////////////////////////////////////////////////////////////// */
/* continuous shape functions */
double cShapeFunc0 ( int n, int m, const point3d *cp, double u, double v,
                    point3d *p, vector3d *nv, vector3d *vv )
{
  return 0.0;
} /*cShapeFunc0*/

double cHighlightLines ( int n, int m, const point3d *cp, double u, double v,
                         point3d *p, vector3d *nv, vector3d *vv )
{
  vector3d dp;
  double   s;

  SubtractPoints3d ( p, &rp0, &dp );
  s = det3d ( nv, &dp, &rv2 ) / det3d ( nv, &rv1, &rv2 );
  return s;
} /*cHighlightLines*/

double cReflectionLines ( int n, int m, const point3d *cp, double u, double v,
                          point3d *p, vector3d *nv, vector3d *vv )
{
  vector3d r, dp;
  double   s;

  s = 2.0*DotProduct3d( nv, vv )/DotProduct3d ( nv, nv );
  SetVector3d ( &r, vv->x - s*nv->x, vv->y - s*nv->y, vv->z - s*nv->z );
  SubtractPoints3d ( p, &rp0, &dp );
  s = det3d ( &r, &dp, &rv2 ) / det3d ( &r, &rv1, &rv2 );
  return s;
} /*cReflectionLines*/

double cCrossSections ( int n, int m, const point3d *cp, double u, double v,
                       point3d *p, vector3d *nv, vector3d *vv )
{
  double s;

  s = DotProduct3d ( p, &sectiondir );
  return s;
} /*cCrossSections*/

double cLambertLight ( int n, int m, const point3d *cp, double u, double v,
                      point3d *p, vector3d *nv, vector3d *vv )
{
  double s;

  s = DotProduct3d ( nv, &lightdir[0] );
  return s;
} /*cLambertLight*/

double cMeanCurvature ( int n, int m, const point3d *cp, double u, double v,
                       point3d *p, vector3d *nv, vector3d *vv )
{
  double gauss, mean;

  mbs_GMCurvaturesBP3d ( n, m, cp, u, v, &gauss, &mean );
  return mean;
} /*cMeanCurvature*/

double cGaussianCurvature ( int n, int m, const point3d *cp, double u, double v,
                           point3d *p, vector3d *nv, vector3d *vv )
{
  double gauss, mean;

  mbs_GMCurvaturesBP3d ( n, m, cp, u, v, &gauss, &mean );
  return gauss;
} /*cMeanCurvature*/

double cRShapeFunc0 ( int n, int m, const point4d *cp, double u, double v,
                     point3d *p, vector3d *nv, vector3d *vv )
{
  return 0.0;
} /*cRShapeFunc0*/

double cRHighlightLines ( int n, int m, const point4d *cp, double u, double v,
                          point3d *p, vector3d *nv, vector3d *vv )
{
  vector3d dp;
  double   s;

  SubtractPoints3d ( p, &rp0, &dp );
  s = det3d ( nv, &dp, &rv2 ) / det3d ( nv, &rv1, &rv2 );
  return s;
} /*cRHighlightLines*/

double cRReflectionLines ( int n, int m, const point4d *cp, double u, double v,
                           point3d *p, vector3d *nv, vector3d *vv )
{
  vector3d r, dp;
  double   s;

  s = 2.0*DotProduct3d( nv, vv )/DotProduct3d ( nv, nv );
  SetVector3d ( &r, vv->x - s*nv->x, vv->y - s*nv->y, vv->z - s*nv->z );
  SubtractPoints3d ( p, &rp0, &dp );
  s = det3d ( &r, &dp, &rv2 ) / det3d ( &r, &rv1, &rv2 );
  return s;
} /*cRReflectionLines*/

double cRCrossSections ( int n, int m, const point4d *cp, double u, double v,
                        point3d *p, vector3d *nv, vector3d *vv )
{
  double s;

  s = DotProduct3d ( p, &sectiondir );
  return s;
} /*cRCrossSections*/

double cRLambertLight ( int n, int m, const point4d *cp, double u, double v,
                       point3d *p, vector3d *nv, vector3d *vv )
{
  double s;

  s = DotProduct3d ( nv, &lightdir[0] );
  return s;
} /*cRLambertLight*/

double cRMeanCurvature ( int n, int m, const point4d *cp, double u, double v,
                        point3d *p, vector3d *nv, vector3d *vv )
{
  double gauss, mean;

  mbs_GMCurvaturesBP3Rd ( n, m, cp, u, v, &gauss, &mean );
  return mean;
} /*cRMeanCurvature*/

double cRGaussianCurvature ( int n, int m, const point4d *cp, double u, double v,
                            point3d *p, vector3d *nv, vector3d *vv )
{
  double gauss, mean;

  mbs_GMCurvaturesBP3Rd ( n, m, cp, u, v, &gauss, &mean );
  return gauss;
} /*cRMeanCurvature*/

/* ///////////////////////////////////////////////////////////////////////// */
void OutputMinMaxShapeFunc ( double minsf, double maxsf )
{
  extern void SetStatusText ( char *text, boolean onscreen );
  char txt[256];

/*  printf ( "Shape function min = %f, max = %f\n", minsf, maxsf ); */
  sprintf ( txt, "Shape function min = %f, max = %f", minsf, maxsf );
  SetStatusText ( txt, true );
} /*OutputMinMaxShapeFunc*/

typedef struct {
    int     dens;
    int     *jobs;
    boolean ok;
    float   fmin, fmax;
  } minmax_struct;

static boolean _FindMinMaxShapeFunc ( void *usrdata, int3 *jobnum )
{
  minmax_struct *mms;
  int           i, i0, i1, j, k, n, m, dens;
  double        u, v, f, fmin, fmax;
  point3d       p, *cp3;
  point4d       *cp4;
  vector3d      nv, vv;
  boolean       ok;

  mms = (minmax_struct*)usrdata;
  dens = mms->dens;
  i = jobnum->x;
  i0 = mms->jobs[i];  i1 = mms->jobs[i+1];
  ok = false;
  fmin = -1.0;  fmax = 1.0;
  for ( i = i0; i < i1; i++ ) {
    switch ( obj_tab[i].type ) {
  case obj_BEZPATCH:
      n = obj_tab[i].bezp.ptree->n;
      m = obj_tab[i].bezp.ptree->m;
      cp3 = obj_tab[i].bezp.ptree->root->ctlpoints;
      if ( !ok ) {
        mbs_BCHornerNvP3d ( n, m, cp3, 0.0, 0.0, &p, &nv );
        SubtractPoints3d ( &p, &_CPos.position, &vv );
        NormalizeVector3d ( &vv );
        fmin = fmax = cShapeFunc ( n, m, cp3, 0.0, 0.0, &p, &nv, &vv );
        ok = true;
      }
      for ( j = 0; j <= dens; j++ ) {
        u = (double)j/(double)dens;
        for ( k = 0; k <= dens; k++ ) {
          v = (double)k/(double)dens;
          mbs_BCHornerNvP3d ( n, m, cp3, u, v, &p, &nv );
          SubtractPoints3d ( &p, &_CPos.position, &vv ); 
          NormalizeVector3d ( &vv );
          f = cShapeFunc ( n, m, cp3, u, v, &p, &nv, &vv );
          if ( f < fmin )
            fmin = f;
          else if ( f > fmax )
            fmax = f;
        }
      }
      break;
  case obj_RBEZPATCH:
      n = obj_tab[i].rbezp.ptree->n;
      m = obj_tab[i].rbezp.ptree->m;
      cp4 = obj_tab[i].rbezp.ptree->root->ctlpoints;
      if ( !ok ) {
        mbs_BCHornerNvP3Rd ( n, m, cp4, 0.0, 0.0, &p, &nv );
        SubtractPoints3d ( &p, &_CPos.position, &vv );
        NormalizeVector3d ( &vv );
        fmin = fmax = cRShapeFunc ( n, m, cp4, 0.0, 0.0, &p, &nv, &vv );
        ok = true;
      }
      for ( j = 0; j <= dens; j++ ) {
        u = (double)j/(double)dens;
        for ( k = 0; k <= dens; k++ ) {
          v = (double)k/(double)dens;
          mbs_BCHornerNvP3Rd ( n, m, cp4, u, v, &p, &nv );
          SubtractPoints3d ( &p, &_CPos.position, &vv ); 
          NormalizeVector3d ( &vv );
          f = cRShapeFunc ( n, m, cp4, u, v, &p, &nv, &vv );
          if ( f < fmin )
            fmin = f;
          else if ( f > fmax )
            fmax = f;
        }
      }
      break;
  default:
      break;
    }
  }
        /* store the minimal and maximal values found so far */
  if ( ncpu > 1 && rendering_npthreads > 1 ) {
    pthread_mutex_lock ( &raybez_mutex );
    if ( mms->ok ) {
      if ( fmin < mms->fmin ) mms->fmin = fmin;
      if ( fmax > mms->fmax ) mms->fmax = fmax;
    }
    else {
      mms->fmin = fmin;
      mms->fmax = fmax;
      mms->ok = true;
    }
    pthread_mutex_unlock ( &raybez_mutex );
  }
  else {
    if ( mms->ok ) {
      if ( fmin < mms->fmin ) mms->fmin = fmin;
      if ( fmax > mms->fmax ) mms->fmax = fmax;
    }
    else {
      mms->fmin = fmin;
      mms->fmax = fmax;
      mms->ok = true;
    }
  }
  return true;
} /*_FindMinMaxShapeFunc*/

void FindMinMaxShapeFunc ( int dens )
{
  void          *sp;
  minmax_struct mms;
  int           i, d, nthr;
  int3          jobsize;
  boolean       success;

  sp = pkv_GetScratchMemTop ();

pkv_Tic ( &tic );

  nthr = ncpu > 1 ? rendering_npthreads : 1;
  nthr = min ( nobjects, nthr );
  mms.dens = dens;
  mms.ok = false;
  mms.fmin = mms.fmax = 0.0;
  mms.jobs = pkv_GetScratchMemi ( nthr+1 );
  if ( !mms.jobs || nobjects <= 0 )
    goto failure;
  mms.jobs[0] = 0;
  mms.jobs[nthr] = nobjects;
  if ( nthr > 1 ) {

printf ( "A: ncpu = %d, nthr = %d, nobjects = %d\n", ncpu, nthr, nobjects );

    d = (nobjects+nthr-1)/nthr;
    for ( i = 1; i < nthr; i++ )
      mms.jobs[i] = mms.jobs[i-1]+d;
    jobsize.x = nthr;  jobsize.y = jobsize.z = 1;
    pkv_SetPThreadsToWork ( &jobsize, nthr, 4*1048576, 4*1048576,
                            (void*)&mms, _FindMinMaxShapeFunc, NULL, NULL,
                            &success );
  }
  else {

printf ( "B\n" );

    jobsize.x = jobsize.y = jobsize.z = 0;
    _FindMinMaxShapeFunc ( (void*)&mms, &jobsize );
  }

failure:
  minsf = mms.fmin + (mms.fmax-mms.fmin)*render_cfrange[0];
  maxsf = mms.fmin + (mms.fmax-mms.fmin)*render_cfrange[1];
  OutputMinMaxShapeFunc ( minsf, maxsf );
  if ( minsf >= maxsf ) {
    minsf -= 2.0;
    maxsf += 3.0;
  }

printf ( "min-max-sf time = %6.2f\n", pkv_Seconds ( pkv_Toc ( &tic ) ) );

  pkv_SetScratchMemTop ( sp );
} /*FindMinMaxShapeFunction*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean RendInit ( void )
{
  int nplanes;

  RendererIsOk = false;

  if ( ncpu > 1 ) {
    if ( !raybez_InitMutex () )
      ncpu = 1;
  }
  nplanes = XDisplayPlanes ( xgedisplay, xgescreen );
  rendimage = XCreateImage ( xgedisplay, xgevisual, nplanes, ZPixmap, 0,
                             NULL, xge_MAX_WIDTH, xge_MAX_HEIGHT, 8, 0 );
  imagedata = (char*)malloc ( xge_MAX_HEIGHT*rendimage->bytes_per_line );
  if ( rendimage && imagedata )
    rendimage->data = imagedata;
  else
    return false;
  obj_tab_length = 100;
  obj_tab = malloc ( obj_tab_length*sizeof(renderobj) );
  if ( !obj_tab )
    return false;
  nobjects = 0;
  cShapeFunc = cShapeFunc0;
  dShapeFunc = dShapeFunc0;
  cRShapeFunc = cRShapeFunc0;
  dRShapeFunc = dRShapeFunc0;
  if ( !(aabuf = malloc ( 7*3*(3*xge_MAX_WIDTH+4) )) )
    return false;
  RendererIsOk = true;
  return true;
} /*RendInit*/

void RendDestroy ( void )
{
  if ( ncpu > 1 )
    raybez_DestroyMutex ();
  RendReset ();
  XDestroyImage ( rendimage );
  free ( obj_tab );
  obj_tab = NULL;
  free ( aabuf );
  aabuf = NULL;
  RendererIsOk = false;
} /*RendDestroy*/

/* ////////////////////////////////////////////////////////////////////////// */
static boolean ReallocObjTab ( void )
{
  renderobj *newtab;
  int       newtablength;

  newtablength = 2*obj_tab_length;
  newtab = malloc ( newtablength*sizeof(renderobj) );
  if ( newtab ) {
    memcpy ( newtab, obj_tab, nobjects*sizeof(renderobj) );
    free ( obj_tab );
    obj_tab = newtab;
    obj_tab_length = newtablength;
    return true;
  }
  else
    return false;
} /*ReallocObjTab*/

boolean RendEnterTriangle3d ( point3d *p0, point3d *p1, point3d *p2,
                              double *colour )
{
#define TOL 1.0e-14
  renderobj_triangle *tr;
  triangle_data      *trd;
  double             m11, m21, m22, l11, l21, l22;
  vector3d           v1, v2;

  trd = NULL;
  if ( nobjects >= obj_tab_length ) {
    if ( !ReallocObjTab () )
      goto failure;
  }
  PKV_MALLOC ( trd, sizeof(triangle_data) );
  if ( trd ) {
    tr = &obj_tab[nobjects].triang;
    tr->type = obj_TRIANGLE;
    memcpy ( tr->colour, colour, 3*sizeof(double) );
    tr->trdata = trd;
        /* find the bounding box */
    trd->bbox.x0 = trd->bbox.x1 = p0->x;
    trd->bbox.y0 = trd->bbox.y1 = p0->y;
    trd->bbox.z0 = trd->bbox.z1 = p0->z;
    if ( p1->x < p2->x ) {
      if ( p1->x < trd->bbox.x0 ) trd->bbox.x0 = p1->x;
      if ( p2->x > trd->bbox.x1 ) trd->bbox.x1 = p2->x;
    }
    else {
      if ( p2->x < trd->bbox.x0 ) trd->bbox.x0 = p2->x;
      if ( p1->x > trd->bbox.x1 ) trd->bbox.x1 = p1->x;
    }
    if ( p1->y < p2->y ) {
      if ( p1->y < trd->bbox.y0 ) trd->bbox.y0 = p1->y;
      if ( p2->y > trd->bbox.y1 ) trd->bbox.y1 = p2->y;
    }
    else {
      if ( p2->y < trd->bbox.y0 ) trd->bbox.y0 = p2->y;
      if ( p1->y > trd->bbox.y1 ) trd->bbox.y1 = p1->y;
    }
    if ( p1->z < p2->z ) {
      if ( p1->z < trd->bbox.z0 ) trd->bbox.z0 = p1->z;
      if ( p2->z > trd->bbox.z1 ) trd->bbox.z1 = p2->z;
    }
    else {
      if ( p2->z < trd->bbox.z0 ) trd->bbox.z0 = p2->z;
      if ( p1->z > trd->bbox.z1 ) trd->bbox.z1 = p1->z;
    }
        /* compute the unit normal vector */
    trd->p0 = *p0;
    SubtractPoints3d ( p1, p0, &v1 );
    SubtractPoints3d ( p2, p0, &v2 );
    CrossProduct3d ( &v1, &v2, &trd->n );
    NormalizeVector3d ( &trd->n );
        /* find the pseudo-inversion of the 3x2 matrix A=[v1,v2] */
          /* compute A^T*A */
    m11 = DotProduct3d ( &v1, &v1 );
    if ( m11 <= TOL )
      goto failure;
    m21 = DotProduct3d ( &v2, &v1 );
    m22 = DotProduct3d ( &v2, &v2 );
          /* this is the Cholesky's decomposition of a 2x2 matrix A^T*A */
    l11 = sqrt ( m11 );
    l21 = m21/l11;
    l22 = m22-m21*m21/m11;
    if ( l22 <= TOL )
      goto failure;
          /* compute the rows a1^T, a2^T, of (A^T*A)^{-1}A^T */
    MultVector3d ( 1.0/l11, &v1, &v1 );
    AddVector3Md ( &v2, &v1, -l21, &v2 );
    MultVector3d ( 1.0/l22, &v2, &trd->a2 );
    AddVector3Md ( &v1, &trd->a2, -l21, &v1 );
    MultVector3d ( 1.0/l11, &v1, &trd->a1 );
    nobjects ++;
    return true;
  }
  else {
failure:
    if ( trd ) PKV_FREE ( trd );
    return false;
  }
#undef TOL
} /*RendEnterTriangle3d*/

boolean RendEnterBezPatch3d ( int n, int m, const point3d *cp, double *colour )
{
  if ( nobjects >= obj_tab_length ) {
    if ( !ReallocObjTab () )
      goto failure;
  }
  if ( nobjects < obj_tab_length ) {
    obj_tab[nobjects].type = obj_BEZPATCH;
    obj_tab[nobjects].bezp.ptree =
          rbez_NewBezPatchTreed ( nobjects, n, m, 0.0, 1.0, 0.0, 1.0, cp );
    if ( obj_tab[nobjects].bezp.ptree ) {
      memcpy ( &obj_tab[nobjects].bezp.colour[0], colour, 3*sizeof(double) );
      nobjects ++;
      return true;
    }
    else
      goto failure;
  }
failure:
  return false;
} /*RendEnterBezPatch3d*/

boolean RendEnterBSPatch3d ( int n, int lknu, const double *knu,  
                             int m, int lknv, const double *knv,  
                             const point3d *cp, double *colour )
{
  void   *sp;
  int    ku, kv, pitch1, pitch2, pitch3, i, j, start;
  double *b, *c;

  sp = pkv_GetScratchMemTop ();
  if ( obj_tab ) {
    ku = mbs_NumKnotIntervalsd ( n, lknu, knu );
    kv = mbs_NumKnotIntervalsd ( m, lknv, knv );
    pitch1 = (lknv-m)*3;
    pitch2 = (m+1)*3*kv;
    pitch3 = (m+1)*3;
    b = pkv_GetScratchMemd ( pitch2*ku*(n+1)*3 );
    c = pkv_GetScratchMemd ( (n+1)*pitch3 );
    if ( b && c ) {
      mbs_BSPatchToBezd ( 3, n, lknu, knu, m, lknv, knv, pitch1, (double*)cp,
                          &ku, NULL, NULL, &kv, NULL, NULL, pitch2, b );
      for ( i = 0; i < ku; i++ )
        for ( j = 0; j < kv; j++ ) {
          start = i*(n+1)*pitch2+pitch3*j;
          pkv_Selectd ( n+1, pitch3, pitch2, pitch3, &b[start], c );
          if ( !RendEnterBezPatch3d ( n, m, (point3d*)c, colour ) )
            goto failure;
        }
      pkv_SetScratchMemTop ( sp );
      return true;
    }
  }
failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*RendEnterBSPatch3d*/

boolean RendEnterBezPatch3Rd ( int n, int m, const point4d *cp, double *colour )
{
  if ( nobjects >= obj_tab_length ) {
    if ( !ReallocObjTab () )
      goto failure;
  }
  if ( nobjects < obj_tab_length ) {
    obj_tab[nobjects].type = obj_RBEZPATCH;
    obj_tab[nobjects].rbezp.ptree =
          rbez_NewRBezPatchTreed ( nobjects, n, m, 0.0, 1.0, 0.0, 1.0, cp );
    if ( obj_tab[nobjects].rbezp.ptree ) {
      memcpy ( &obj_tab[nobjects].rbezp.colour[0], colour, 3*sizeof(double) );
      nobjects ++;
      return true;
    }
    else
      goto failure;
  }
failure:
  return false;
} /*RendEnterBezPatch3Rd*/

boolean RendEnterBSPatch3Rd ( int n, int lknu, const double *knu,  
                              int m, int lknv, const double *knv,  
                              const point4d *cp, double *colour )
{
  void   *sp;
  int    ku, kv, pitch1, pitch2, pitch3, i, j, start;
  double *b, *c;

  sp = pkv_GetScratchMemTop ();
  if ( obj_tab ) {
    ku = mbs_NumKnotIntervalsd ( n, lknu, knu );
    kv = mbs_NumKnotIntervalsd ( m, lknv, knv );
    pitch1 = (lknv-m)*4;
    pitch2 = (m+1)*4*kv;
    pitch3 = (m+1)*4;
    b = pkv_GetScratchMemd ( pitch2*ku*(n+1)*4 );
    c = pkv_GetScratchMemd ( (n+1)*pitch3 );
    if ( b && c ) {
      mbs_BSPatchToBezd ( 4, n, lknu, knu, m, lknv, knv, pitch1, (double*)cp,
                          &ku, NULL, NULL, &kv, NULL, NULL, pitch2, b );
      for ( i = 0; i < ku; i++ )
        for ( j = 0; j < kv; j++ ) {
          start = i*(n+1)*pitch2+pitch3*j;
          pkv_Selectd ( n+1, pitch3, pitch2, pitch3, &b[start], c );
          if ( !RendEnterBezPatch3Rd ( n, m, (point4d*)c, colour ) )
            goto failure;
        }
      pkv_SetScratchMemTop ( sp );
      return true;
    }
  }
failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*RendEnterBSPatch3Rd*/

boolean RendEnterBezCurve3d ( int n, point3d *cp, double r, double *colour )
{
  if ( nobjects >= obj_tab_length ) {
    if ( !ReallocObjTab () )
      goto failure;
  }
  if ( nobjects < obj_tab_length ) {
    obj_tab[nobjects].type = obj_BEZCURVE;
    obj_tab[nobjects].bezc.ctree =
          rbez_NewBezCurveTreed ( nobjects, n, 0.0, 1.0, r, cp );
    if ( obj_tab[nobjects].bezc.ctree ) {
      memcpy ( &obj_tab[nobjects].bezc.colour[0], colour, 3*sizeof(double) );
      nobjects ++;
      return true;
    }
    else
      goto failure;
  }
failure:
  return false;
} /*RendEnterBezCurve3d*/

boolean RendEnterBSCurve3d ( int n, int lkn, const double *kn,
                             point3d *cp, double r, double *colour )
{
  void   *sp;
  int    i, k;
  double *b;

  sp = pkv_GetScratchMemTop ();
  if ( obj_tab ) {
    k = mbs_NumKnotIntervalsd ( n, lkn, kn );
    if ( k < 1 )
      goto failure;
    b = pkv_GetScratchMemd ( (n+1)*k*3 );
    if ( b ) {
      mbs_BSToBezC3d ( n, lkn, kn, cp, &k, NULL, NULL, b );
      for ( i = 0; i < k; i++ )
        if ( !RendEnterBezCurve3d ( n, (point3d*)&b[3*i*(n+1)], r, colour ) )
          goto failure;
/*
bsf_OpenOutputFile ( "mbc.bs" );
for ( i = 0; i < k; i++ )
  bsf_WriteBezierCurved ( 3, 3, false, n, &b[3*i*(n+1)], NULL, NULL );
bsf_CloseOutputFile ();
*/
    }
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*RendEnterBSCurve3d*/

boolean RendEnterBezCurve3Rd ( int n, point4d *cp, double r, double *colour )
{
  if ( nobjects >= obj_tab_length ) {
    if ( !ReallocObjTab () )
      goto failure;
  }
  if ( nobjects < obj_tab_length ) {
    obj_tab[nobjects].type = obj_RBEZCURVE;
    obj_tab[nobjects].rbezc.ctree =
          rbez_NewRBezCurveTreed ( nobjects, n, 0.0, 1.0, r, cp );
    if ( obj_tab[nobjects].rbezc.ctree ) {
      memcpy ( &obj_tab[nobjects].rbezc.colour[0], colour, 3*sizeof(double) );
      nobjects ++;
      return true;
    }
    else
      goto failure;
  }
failure:
  return false;
} /*RendEnterBezCurve3Rd*/

boolean RendEnterBSCurve3Rd ( int n, int lkn, const double *kn,
                              point4d *cp, double r, double *colour )
{
  void   *sp;
  int    i, k;
  double *b;

  sp = pkv_GetScratchMemTop ();
  if ( obj_tab ) {
    k = mbs_NumKnotIntervalsd ( n, lkn, kn );
    if ( k < 1 )
      goto failure;
    b = pkv_GetScratchMemd ( (n+1)*k*4 );
    if ( b ) {
      mbs_BSToBezC4d ( n, lkn, kn, cp, &k, NULL, NULL, b );
      for ( i = 0; i < k; i++ )
        if ( !RendEnterBezCurve3Rd ( n, (point4d*)&b[4*i*(n+1)], r, colour ) )
          goto failure;
    }
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*RendEnterBSCurve3Rd*/

tnode* NewPatchNode ( int k )
{
  tnode *node;

  node = malloc ( sizeof(tnode) );
  if ( !node )
    exit ( 1 );
  node->left = node->right = NULL;
  node->k = k;
  switch ( obj_tab[k].type ) {
case obj_TRIANGLE:
    memcpy ( &node->bbox, &obj_tab[k].triang.trdata->bbox, sizeof(Box3d) );
    break;
case obj_BEZPATCH:
    memcpy ( &node->bbox, &obj_tab[k].bezp.ptree->root->bbox, sizeof(Box3d) );
    break;
case obj_RBEZPATCH:
    memcpy ( &node->bbox, &obj_tab[k].rbezp.ptree->root->bbox, sizeof(Box3d) );
    break;
case obj_BEZCURVE:
    memcpy ( &node->bbox, &obj_tab[k].bezc.ctree->root->bbox, sizeof(Box3d) );
    break;
case obj_RBEZCURVE:
    memcpy ( &node->bbox, &obj_tab[k].rbezc.ctree->root->bbox, sizeof(Box3d) );
    break;
default:
    free ( node );
    return NULL;
  }
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
  rbez_FindSumBBoxd ( &left->bbox, &right->bbox, &node->bbox );
  return node;
} /*NewInnerNode*/

double BoxDiameter ( Box3d *box )
{
  double a, b, c;

  a = box->x1 - box->x0;
  b = box->y1 - box->y0;
  c = box->z1 - box->z0;
  return a*a + b*b + c*c;
} /*BoxDiameter*/

boolean CompNodePrio ( void *node1, void *node2 )
{
  double d1, d2;

  d1 = BoxDiameter ( &((tnode*)node1)->bbox );
  d2 = BoxDiameter ( &((tnode*)node2)->bbox );
  return d1 < d2;
} /*CompNodePrio*/

tnode* BuildObjectTree ( int nobjects )
{
  void   *sp;
  void   **pqueue;
  int    i, j, k, last;
  tnode  *tn0, *tn1, *tn2;
  Box3d  box;
  double d1, d2;

  sp = pkv_GetScratchMemTop ();
  if ( !nobjects )
    return NULL;
  pqueue = malloc ( nobjects*sizeof(void*) );
  if ( !pqueue )
    exit ( 1 );
  for ( i = 0; i < nobjects; i++ )
    pqueue[i] = NewPatchNode ( i );
  pkv_HeapOrder ( pqueue, nobjects, CompNodePrio );
  for ( i = nobjects-1; i > 0; i-- ) {
    tn0 = pqueue[0];
    last = i;
    pkv_HeapRemove ( pqueue, &last, 0, CompNodePrio );
    tn1 = pqueue[0];
    rbez_FindSumBBoxd ( &tn0->bbox, &tn1->bbox, &box );
    d1 = BoxDiameter ( &box );
    k = 0;
    for ( j = 1; j <= last; j++ ) {
      tn2 = pqueue[j];
      rbez_FindSumBBoxd ( &tn0->bbox, &tn2->bbox, &box );
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
} /*BuildObjectTree*/

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
    for ( i = 0; i < nobjects; i++ )
      switch ( obj_tab[i].type ) {
    case obj_TRIANGLE:
        PKV_FREE ( obj_tab[i].triang.trdata );
        break;
    case obj_BEZPATCH:
        rbez_DestroyBezPatchTreed ( obj_tab[i].bezp.ptree );
        break;
    case obj_RBEZPATCH:
        rbez_DestroyRBezPatchTreed ( obj_tab[i].rbezp.ptree );
        break;
    case obj_BEZCURVE:
        rbez_DestroyBezCurveTreed ( obj_tab[i].bezc.ctree );
        break;
    case obj_RBEZCURVE:
        rbez_DestroyRBezCurveTreed ( obj_tab[i].rbezc.ctree );
        break;
    default:
        break;
      }
    free ( obj_tab );
    obj_tab_length = 100;
    obj_tab = malloc ( obj_tab_length*sizeof(renderobj) );
    if ( !obj_tab )
      obj_tab_length = 0;
#ifdef COPY_POSTSCRIPT
    DumpPostscriptPicture ();
#endif
  }
  nobjects = 0;
  return obj_tab_length > 0;
} /*RendReset*/

boolean RendEnterCamerad ( CameraRecd *CPos, xge_widget *er )
{
  if ( CPos->parallel )
    return false;
      /* convert from double to single precision */
  _CPos = *CPos;
  CameraSetMappingd ( &_CPos );
  _xgw = er;
  return true;
} /*RendEnterCamerad*/

void RendEnterLightsd ( int nlights, const vector3d *light_dir,
                        const double *light_int )
{
  int   i;
  double lint;

  memset ( &lightint[0], 0, (R_NLIGHTS+1)*sizeof(double) );
  for ( i = 0; i < nlights && i < R_NLIGHTS; i++ ) {
    SetVector3d ( &lightdir[i], light_dir[i].x, light_dir[i].y, light_dir[i].z );
    NormalizeVector3d ( &lightdir[i] );
    lint = min ( 1.0, light_int[i] );
    lightint[i] = max ( 0.0, lint );
  }
  lightint[R_NLIGHTS] = light_int[nlights];  /* ambient light intensity */
} /*RendEnterLightsd*/

void RendEnterReflectionLinesFramed ( point3d rf[3] )
{
  rp0 = rf[0];
  SubtractPoints3d ( &rf[1], &rf[0], &rv1 );
  SubtractPoints3d ( &rf[2], &rf[0], &rv2 );
} /*RendEnterReflectionLinesFramed*/

void RendEnterHighlightLinesFramed ( point3d hf[3] )
{
  hp0 = hf[0];
  SubtractPoints3d ( &hf[1], &hf[0], &hv1 );
  SubtractPoints3d ( &hf[2], &hf[0], &hv2 );
} /*RendEnterHighlightLinesFramed*/

void RendEnterSectionPlanesNormald ( vector3d *spn )
{
  SetVector3d ( &sectiondir, (double)spn->x, (double)spn->y, (double)spn->z  );
} /*RendEnterSectionPlanesNormald*/

boolean RendBegin ( void )
{
  y = _xgw->y;

  dfsf = xge_LogSlidebarValued ( R_MINDFSF, R_MAXDFSF, render_dfsf );
        /* build the scene hierarchy binary tree */
  root = BuildObjectTree ( nobjects );
        /* setup the shape function texturing procedures */
  GetTexColour = GetTexColourSF;
  if ( swGaussian_c ) {
    cShapeFunc = cGaussianCurvature;
    cRShapeFunc = cRGaussianCurvature;
  }
  else if ( swMean_c ) {
    cShapeFunc = cMeanCurvature;
    cRShapeFunc = cRMeanCurvature;
  }
  else if ( swLambert_c ) {
    cShapeFunc = cLambertLight;
    cRShapeFunc = cRLambertLight;
  }
  else if ( swReflection_c ) {
    cShapeFunc = cReflectionLines;
    cRShapeFunc = cRReflectionLines;
  }
  else if ( swHighlight_c ) {
    cShapeFunc = cHighlightLines;
    cRShapeFunc = cRHighlightLines;
  }
  else if ( swSections_c ) {
    cShapeFunc = cCrossSections;
    cRShapeFunc = cRCrossSections;
  }
  else {
    GetTexColour = GetTexColourNoSF;
    cShapeFunc = cShapeFunc0;
    cRShapeFunc = cRShapeFunc0;
  }
  if ( swGaussian_d ) {
    dShapeFunc = dGaussianCurvature;
    dRShapeFunc = dRGaussianCurvature;
  }
  else if ( swMean_d ) {
    dShapeFunc = dMeanCurvature;
    dRShapeFunc = dRMeanCurvature;
  }
  else if ( swLambert_d ) {
    dShapeFunc = dLambertLight;
    dRShapeFunc = dRLambertLight;
  }
  else if ( swReflection_d ) {
    dShapeFunc = dReflectionLines;
    dRShapeFunc = dRReflectionLines;
  }
  else if ( swHighlight_d ) {
    dShapeFunc = dHighlightLines;
    dRShapeFunc = dRHighlightLines;
  }
  else if ( swSections_d ) {
    dShapeFunc = dCrossSections;
    dRShapeFunc = dRCrossSections;
  }
  else {
    dShapeFunc = dShapeFunc0;
    dRShapeFunc = dRShapeFunc0;
  }
  FindMinMaxShapeFunc ( 10 );

pkv_Tic ( &tic );

  renderer_npthreads = rendering_npthreads = ncpu > 1 ? rendering_npthreads : 1;
  RenderingIsOn = true;
  if ( swAntialias )
    InitRenderingAA ();
  return true;
} /*RendBegin*/

boolean RendRestart ( void )
{
  if ( RenderingIsOn ) {

pkv_Tic ( &tic );

    y = _xgw->y;
    if ( swAntialias )
      InitRenderingAA ();
    return true;
  }
  else
    return false;
} /*RendRestart*/

/* ///////////////////////////////////////////////////////////////////////// */
static void _rend_FindRayTriangleIntersd ( triangle_data *trd, ray3d *ray,
                                           int *nint, RayObjectIntersd *inters )
{
#define TOL 1.0e-10
  double   a, b, t;
  vector3d pq;
  int      i;

  SubtractPoints3d ( &trd->p0, &ray->p, &pq );
  a = DotProduct3d ( &ray->v, &trd->n );
  if ( fabs(a) < TOL )
    return;
  b = DotProduct3d ( &pq, &trd->n );
  t = b/a;
  if ( t <= 0.0 )
    return;
  AddVector3Md ( &pq, &ray->v, -t, &pq );
  a = DotProduct3d ( &trd->a1, &pq );
  if ( a >= 0.0 )
    return;
  b = DotProduct3d ( &trd->a2, &pq );
  if ( b >= 0.0 || (a+b) <= -1.0 )
    return;
  i = *nint;
  inters[i].object_id = 0;
  AddVector3Md ( &ray->p, &ray->v, t, &inters[i].p );
  inters[i].nv = trd->n;
  inters[i].u = a;  inters[i].v = b;  inters[i].t = t;
  *nint = i+1;
#undef TOL
} /*_rend_FindRayTriangleIntersd*/

boolean rFindRayTrInters ( ray3d *ray, tnode *node,
                           RayObjectIntersd *inters, int *k )
{
  void *sp;
  int  _k, nint;

  sp = pkv_GetScratchMemTop ();
  *k = -1;
  nint = 0;
  _k = node->k;
  if ( _k >= 0 ) {
    RayObjectIntersd *_inters, *_intp;
    int              i;
           /* it is a leaf, i.e. a Bezier patch or a Bezier curve; no ray/box */
           /* intersection is tested here, as the ray/patch or ray/curve offset */
           /* intersection procedures have their own tests */

    _inters = pkv_GetScratchMem ( MAXINTERS*sizeof(RayObjectIntersd) );
    if ( !_inters )
      return false;
    switch ( obj_tab[_k].type ) {
  case obj_TRIANGLE:
      _rend_FindRayTriangleIntersd ( obj_tab[_k].triang.trdata, ray,
                                     &nint, _inters );
      break;
  case obj_BEZPATCH:
      rbez_FindRayBezPatchIntersd ( obj_tab[_k].bezp.ptree, ray,
                    MAXPLEVEL, MAXINTERS, &nint, _inters );
      break;
  case obj_RBEZPATCH:
      rbez_FindRayRBezPatchIntersd ( obj_tab[_k].rbezp.ptree, ray,
                    MAXPLEVEL, MAXINTERS, &nint, _inters );
      break;
  case obj_BEZCURVE:
      rbez_FindRayBezcOffsetIntersd ( obj_tab[_k].bezc.ctree, ray,
                    MAXCLEVEL, MAXINTERS, &nint, _inters );
      break;
  case obj_RBEZCURVE:
      rbez_FindRayRBezcOffsetIntersd ( obj_tab[_k].rbezc.ctree, ray,
                    MAXCLEVEL, MAXINTERS, &nint, _inters );
      break;
  default:
      pkv_SetScratchMemTop ( sp );
      return false;
    }
    if ( nint > 0 ) {
           /* ignore all intersections too close to the ray origin */
      for ( i = 0; i < nint; )
        if ( _inters[i].t < MIN_T )
          _inters[i] = _inters[--nint];
        else
          i ++;
      if ( nint ) {
           /* something was left */
        _intp = _inters;
        for ( i = 1; i < nint; i++ )
          if ( _inters[i].t < _inters[0].t )
            _intp = &_inters[i];
        *inters = *_intp;
        *k = _k;
      }
    }
  }
  else if ( rbez_TestRayBBoxd ( ray, &node->bbox ) ) {
    int              k0, k1;
    boolean          r0, r1;
    RayObjectIntersd inters0, inters1;

    r0 = rFindRayTrInters ( ray, node->left, &inters0, &k0 );
    r1 = rFindRayTrInters ( ray, node->right, &inters1, &k1 );
    if ( r0 ) {
      if ( r1 ) {
        if ( inters0.t < inters1.t )
          { *inters = inters0;  *k = k0; }
        else
          { *inters = inters1;  *k = k1; }
      }
      else
        { *inters = inters0;  *k = k0; }
      nint = 1;
    }
    else {
      if ( r1 )
        { *inters = inters1;  *k = k1;  nint = 1; }
      else nint = 0;
    }
  }
  pkv_SetScratchMemTop ( sp );
  return nint > 0;
} /*rFindRayTrInters*/

boolean FindRayTrInters ( ray3d *ray, RayObjectIntersd *inters, int *k )
{
  if ( !root )
    return false;
  return rFindRayTrInters ( ray, root, inters, k );
} /*FindRayTrInters*/

boolean rFindShadowRayInters ( ray3d *ray, tnode *node )
{
  int              _k, nint;
  RayObjectIntersd inters;

  _k = node->k;
  nint = 0;
  if ( _k >= 0 ) {
           /* it is a leaf, i.e. a Bezier patch or a Bezier curve; no ray/box */
           /* intersection is tested here, as the ray/patch or ray/curve offset */
           /* intersection procedures have their own tests */
    switch ( obj_tab[_k].type ) {
  case obj_TRIANGLE:
      _rend_FindRayTriangleIntersd ( obj_tab[_k].triang.trdata, ray,
                                     &nint, &inters );
      break;
  case obj_BEZPATCH:
      rbez_FindRayBezPatchIntersd ( obj_tab[_k].bezp.ptree, ray,
                    MAXPLEVEL, 1, &nint, &inters );
      break;
  case obj_RBEZPATCH:
      rbez_FindRayRBezPatchIntersd ( obj_tab[_k].rbezp.ptree, ray,
                    MAXPLEVEL, 1, &nint, &inters );
      break;
  case obj_BEZCURVE:
      rbez_FindRayBezcOffsetIntersd ( obj_tab[_k].bezc.ctree, ray,
                    MAXCLEVEL, 1, &nint, &inters );
      break;
  case obj_RBEZCURVE:
      rbez_FindRayRBezcOffsetIntersd ( obj_tab[_k].rbezc.ctree, ray,
                    MAXCLEVEL, 1, &nint, &inters );
      break;
  default:
      return false;
    }
    return nint > 0;
  }
  else if ( rbez_TestRayBBoxd ( ray, &node->bbox ) ) {
    if ( rFindShadowRayInters ( ray, node->left ) )
      return true;
    return rFindShadowRayInters ( ray, node->right );
  }
  else
    return false;
} /*rFindShadowRayInters*/

boolean IsInShadow ( point3d *p, vector3d *light )
{
  ray3d ray;

  if ( !root )
    return false;
  ray.p = *p;
  ray.v = *light;
  AddVector3Md ( p, light, MIN_T, &ray.p );
  return rFindShadowRayInters ( &ray, root );
} /*IsInShadow*/

void GetTexColourSF ( char type, double *colour,
                      int n, int m, const point4d *cp, double u, double v,
                      point3d *p, vector3d *nv, vector3d *vv,
                      vector3d *texcolour )
{
  double t;

  switch ( type ) {
case obj_BEZPATCH:
    t = cShapeFunc ( n, m, (point3d*)cp, u, v, p, nv, vv );
    break;
case obj_RBEZPATCH:
    t = cRShapeFunc ( n, m, cp, u, v, p, nv, vv );
    break;
case obj_TRIANGLE:
case obj_BEZCURVE:
case obj_RBEZCURVE:
default:
    memcpy ( texcolour, colour, sizeof(vector3d) );
    return;
  }
  t = (t-minsf)/(maxsf-minsf);
  t = max ( 0.0, t );
  t = min ( 1.0, t );
  /* "rainbow" texture */
  if      ( t < 0.2 ) SetPoint3d ( texcolour, 1.0, 0.0, 1.0-5.0*t );  
  else if ( t < 0.4 ) SetPoint3d ( texcolour, 1.0, 5.0*(t-0.2), 0.0 );
  else if ( t < 0.6 ) SetPoint3d ( texcolour, 1.0-5.0*(t-0.4), 1.0, 0.0 );
  else if ( t < 0.8 ) SetPoint3d ( texcolour, 0.0, 1.0, 5.0*(t-0.6) );
  else                SetPoint3d ( texcolour, 0.0, 1.0-5.0*(t-0.8), 1.0 );
  switch ( type ) {
case obj_BEZPATCH:
    if ( !dShapeFunc ( n, m, (point3d*)cp, u, v, p, nv, vv ) )
      MultVector3d ( 0.85, texcolour, texcolour );
    break;
case obj_RBEZPATCH:
    if ( !dRShapeFunc ( n, m, cp, u, v, p, nv, vv ) )
      MultVector3d ( 0.85, texcolour, texcolour );
    break;
default:
    break;
  }
} /*GetTexColourSF*/

void GetTexColourNoSF ( char type, double *colour,
                        int n, int m, const point4d *cp, double u, double v,
                        point3d *p, vector3d *nv, vector3d *vv,
                        vector3d *texcolour )
{
  memcpy ( texcolour, colour, 3*sizeof(double) );
  switch ( type ) {
case obj_TRIANGLE:
    break;
case obj_BEZPATCH:
    if ( !dShapeFunc ( n, m, (point3d*)cp, u, v, p, nv, vv ) )
      MultVector3d ( 0.85, texcolour, texcolour );
    break;
case obj_RBEZPATCH:
    if ( !dRShapeFunc ( n, m, cp, u, v, p, nv, vv ) )
      MultVector3d ( 0.85, texcolour, texcolour );
    break;
case obj_BEZCURVE:
    break;
case obj_RBEZCURVE:
    break;
default:
    break;
  }
} /*GetTexColourNoSF*/

void AmbientLight ( vector3d *tex, double amb_intens, vector3d *rgb )
{
  MultVector3d ( amb_intens, tex, rgb );
} /*AmbientLight*/

void PhongLight ( vector3d *tex, point3d *p, vector3d *nv,
                  vector3d *vv, vector3d *light, double light_intens,
                  vector3d *rgb )
{
  double    c, d;
  vector3d h;

  NormalizeVector3d ( nv );
  c = DotProduct3d ( nv, light );
  d = DotProduct3d ( nv, vv );
  if ( c*d <= 0.0 ) {  
    if ( swShadows ) {
      if ( IsInShadow ( p, light ) )
        return;
    }
    SubtractPoints3d ( vv, light, &h );
    NormalizeVector3d ( &h );
    d = fabs ( DotProduct3d ( nv, &h ) );
    d = d*d;  d = d*d;  d = d*d;  d = d*d;  d = d*d;  d = d*d;
    c = 0.8*light_intens*fabs(c);
    d *= 0.15*light_intens;
    rgb->x += (c*tex->x + d);
    rgb->y += (c*tex->y + d);
    rgb->z += (c*tex->z + d);
  }
} /*PhongLight*/

void GetPixelColour ( ray3d *ray, unsigned int *r, unsigned int *g,
                      unsigned int *b )
{
  int              k;
  RayObjectIntersd iint;
  vector3d         texcolour, rgb;
  unsigned int     a;

  *r = *g = *b = 220;
  if ( FindRayTrInters ( ray, &iint, &k ) ) {
    switch ( obj_tab[k].type ) {
  case obj_TRIANGLE:
      GetTexColour ( obj_TRIANGLE, obj_tab[k].triang.colour, 0, 0,
                     (point4d*)&obj_tab[k].triang.trdata->p0,
                     iint.u, iint.v, &iint.p, &iint.nv, &ray->v, &texcolour );
      break;
  case obj_BEZPATCH:
      GetTexColour ( obj_BEZPATCH, obj_tab[k].bezp.colour,
                     obj_tab[k].bezp.ptree->n, obj_tab[k].bezp.ptree->m,
                     (point4d*)obj_tab[k].bezp.ptree->root->ctlpoints,
                     iint.u, iint.v, &iint.p, &iint.nv, &ray->v, &texcolour );
      break;
  case obj_RBEZPATCH:
      GetTexColour ( obj_RBEZPATCH, obj_tab[k].rbezp.colour,
                     obj_tab[k].rbezp.ptree->n, obj_tab[k].rbezp.ptree->m,
                     obj_tab[k].rbezp.ptree->root->ctlpoints,
                     iint.u, iint.v, &iint.p, &iint.nv, &ray->v, &texcolour );
      break;
  case obj_BEZCURVE:
      GetTexColour ( obj_BEZCURVE, obj_tab[k].bezc.colour,
                     obj_tab[k].bezc.ctree->degree, 0,
                     (point4d*)obj_tab[k].bezc.ctree->root->ctlpoints,
                     iint.u, iint.v, &iint.p, &iint.nv, &ray->v, &texcolour );
      break;
  case obj_RBEZCURVE:
      GetTexColour ( obj_RBEZCURVE, obj_tab[k].rbezc.colour,
                     obj_tab[k].rbezc.ctree->degree, 0,
                     (point4d*)obj_tab[k].rbezc.ctree->root->ctlpoints,
                     iint.u, iint.v, &iint.p, &iint.nv, &ray->v, &texcolour );
      break;
  default:
      break;
    }
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

/* ////////////////////////////////////////////////////////////////////////// */
#ifdef OLD_ATTEMPT
typedef struct {
    unsigned int *buf;
  } pixeldata;

static boolean RenderPixelA ( void *usrdata, int3 *jobnum )
{
  ray3d        ray;
  int          x, xx;
  unsigned int *buf;

  buf = ((pixeldata*)usrdata)->buf;
        /* get the ray */
  x = _xgw->x+jobnum->x;
  xx = 3*jobnum->x;
  CameraRayOfPixeld ( &_CPos, (double)x, (double)y, &ray );
        /* compute the pixel colour */
  GetPixelColour ( &ray, &buf[xx], &buf[xx+1], &buf[xx+2] );
  return true;
} /*RenderPixelA*/

int RenderLineA ( void )
{
  void         *sp;
  int          x, xx;
  ray3d        ray;
  unsigned int r, g, b;
  pixeldata    pd;
  int3         jobsize;
  boolean      success;

  if ( !RenderingIsOn )
    return y;

  sp = pkv_GetScratchMemTop ();
  if ( renderer_threads > 1 ) {
    pd.buf = (unsigned int*)pkv_GetScratchMemi ( 3*_xgw->w );
    if ( !pd.buf )
      goto sequential;
    jobsize.x = _xgw->w;  jobsize.y = jobsize.z = 1;
    if ( !pkv_SetPThreadsToWork ( &jobsize, renderer_threads, 4*1048576, 4*1048576,
                                  (void*)&pd, RenderPixelA, NULL, NULL,
                                  &success ) )
      goto sequential;
    for ( x = _xgw->x, xx = 0;  x < _xgw->x+_xgw->w;  x++, xx+=3 )
      XPutPixel ( rendimage, x, y,
                  xge_PixelColour ( pd.buf[xx], pd.buf[xx+1], pd.buf[xx+2] ) );
  }
  else {
sequential:
    for ( x = _xgw->x; x < _xgw->x+_xgw->w; x++ ) {
        /* get the ray */
      CameraRayOfPixeld ( &_CPos, (double)x, (double)y, &ray );
        /* compute the pixel colour */
      GetPixelColour ( &ray, &r, &g, &b );
      XPutPixel ( rendimage, x, y, xge_PixelColour ( r, g, b ) );
    }
  }
  y++;
  if ( y >= _xgw->y+_xgw->h ) {

printf ( "rendering time = %6.2f\n", pkv_Seconds ( pkv_Toc ( &tic ) ) );

    RenderingIsOn = false;
  }
  pkv_SetScratchMemTop ( sp );
  return y;
} /*RenderLineA*/

#else
typedef struct {
    int         *jobs;
    unsigned int *buf;
  } render_struct;

static boolean RenderSublineA ( void *usrdata, int3 *jobnum )
{
  render_struct *rs;
  ray3d        ray;
  int          x, xx, x0, x1;
  unsigned int *buf;

  rs = (render_struct*)usrdata;
  buf = rs->buf;
  x = jobnum->x;
  x0 = rs->jobs[x];
  x1 = rs->jobs[x+1];
  for ( x = x0, xx = 3*(x-rs->jobs[0]);  x < x1;  x++, xx += 3 ) {
        /* get the ray */
    CameraRayOfPixeld ( &_CPos, (double)x, (double)y, &ray );
        /* compute the pixel colour */
    GetPixelColour ( &ray, &buf[xx], &buf[xx+1], &buf[xx+2] );
  }
  return true;
} /*RenderSublineA*/

int RenderLineA ( void )
{
  void          *sp;
  render_struct rs;
  int           nthr;
  int           x, xx;
  ray3d         ray;
  unsigned int  r, g, b;
  int3          jobsize;
  boolean       success;

  if ( !RenderingIsOn )
    return y;

  sp = pkv_GetScratchMemTop ();
  nthr = ncpu > 1 ? rendering_npthreads : 1;
  nthr = min ( _xgw->w, nthr );
  rs.jobs = pkv_GetScratchMemi ( nthr+1 );
  if ( !rs.jobs ) {
    RenderingIsOn = false;
    return y;
  }
  rs.jobs[0] = _xgw->x;
  rs.jobs[nthr] = _xgw->x+_xgw->w;
  if ( nthr > 1 ) {
    rs.buf = (unsigned int*)pkv_GetScratchMemi ( 3*_xgw->w );
    if ( !rs.buf )
      goto sequential;
    xx = (_xgw->w+nthr-1)/nthr;
    for ( x = 1; x < nthr; x++ )
      rs.jobs[x] = rs.jobs[x-1]+xx;
    jobsize.x = nthr;  jobsize.y = jobsize.z = 1;
    if ( !pkv_SetPThreadsToWork ( &jobsize, nthr, 4*1048576, 4*1048576,
                                  (void*)&rs, RenderSublineA, NULL, NULL,
                                  &success ) )
      goto sequential;
    for ( x = _xgw->x, xx = 0;  x < _xgw->x+_xgw->w;  x++, xx+=3 )
      XPutPixel ( rendimage, x, y,
                  xge_PixelColour ( rs.buf[xx], rs.buf[xx+1], rs.buf[xx+2] ) );
  }
  else {
sequential:
    for ( x = _xgw->x; x < _xgw->x+_xgw->w; x++ ) {
        /* get the ray */
      CameraRayOfPixeld ( &_CPos, (double)x, (double)y, &ray );
        /* compute the pixel colour */
      GetPixelColour ( &ray, &r, &g, &b );
      XPutPixel ( rendimage, x, y, xge_PixelColour ( r, g, b ) );
    }
  }
  y++;
  if ( y >= _xgw->y+_xgw->h ) {

printf ( "rendering time = %6.2f\n", pkv_Seconds ( pkv_Toc ( &tic ) ) );

    RenderingIsOn = false;
  }
  pkv_SetScratchMemTop ( sp );
  return y;
} /*RenderLineA*/

#endif
/* ////////////////////////////////////////////////////////////////////////// */
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

void GetSample ( double x, double y, byte *rgb )
{
  ray3d ray;
  unsigned int r, g, b;

  CameraRayOfPixeld ( &_CPos, (double)x, y, &ray );
  GetPixelColour ( &ray, &r, &g, &b );
  rgb[0] = (byte)r;
  rgb[1] = (byte)g;
  rgb[2] = (byte)b;
} /*GetSample*/

void RenderSubpixelLine ( double y, byte *aaline )
{
  int          x, i;

  if ( !RenderingIsOn )
    return;
  for ( i = 0, x = _xgw->x-1;  x <= _xgw->x+_xgw->w;  x++, i += 9 )
    GetSample ( (double)x, y, &aaline[i] );
  for ( i = 0, x = _xgw->x-1; x < _xgw->x+_xgw->w; x++, i += 9 ) {
    if ( OverThreshold1 ( &aaline[i], &aaline[i+9] ) ) {
      GetSample ( (double)x+(double)(1.0/3.0), y, &aaline[i+3] );
      GetSample ( (double)x+(double)(2.0/3.0), y, &aaline[i+6] );
    }
    else {
      Interpol2_3 ( &aaline[i], &aaline[i+9], &aaline[i+3] );
      Interpol2_3 ( &aaline[i+9], &aaline[i], &aaline[i+6] );
    }
  }
} /*RenderSubpixelLine*/

void GetSupersamples ( double y )
{
  int i, x;

  for ( i = 0, x = _xgw->x-1;  x <= _xgw->x+_xgw->w;  x++, i += 9 ) {
    if ( OverThreshold1 ( &aaline[3][i], &aaline[6][i] ) ) {
      GetSample ( (double)x, (double)y+(double)(1.0/3.0), &aaline[4][i] );
      GetSample ( (double)x, (double)y+(double)(2.0/3.0), &aaline[5][i] );
    }
    else {
      Interpol2_3 ( &aaline[3][i], &aaline[6][i], &aaline[4][i] );
      Interpol2_3 ( &aaline[6][i], &aaline[3][i], &aaline[5][i] );
    }
  }
  for ( i = 0, x = _xgw->x-1;  x < _xgw->x+_xgw->w;  x++, i += 9 ) {
    if ( OverThreshold2 ( &aaline[3][i], &aaline[3][i+9],
                          &aaline[6][i], &aaline[6][i+9] ) ) {
      GetSample ( (double)x+(double)(1.0/3.0), y+(double)(1.0/3.0), &aaline[4][i+3] );
      GetSample ( (double)x+(double)(2.0/3.0), y+(double)(1.0/3.0), &aaline[4][i+6] );
      GetSample ( (double)x+(double)(1.0/3.0), y+(double)(2.0/3.0), &aaline[5][i+3] );
      GetSample ( (double)x+(double)(2.0/3.0), y+(double)(2.0/3.0), &aaline[5][i+6] );
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
  RenderSubpixelLine ( (double)(y-1), aaline[3] );
  RenderSubpixelLine ( (double)y, aaline[6] );
  GetSupersamples ( (double)(y-1) );
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
    RenderSubpixelLine ( (double)(y+1), aaline[6] );
    GetSupersamples ( (double)y );
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
      RenderingIsOn = false;
  }
  return y;
} /*RenderLineAA*/

int RenderLine ( void )
{
  if ( swAntialias )
    return RenderLineAA ();
  else
    return RenderLineA ();
} /*RenderLine*/

