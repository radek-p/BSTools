
#include <stdio.h>
#include <math.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camera.h"
#include "multibs.h"
#include "raybez.h"
#include "eg1holef.h"

#include "oldxgedit.h"
#include "g1ekernel.h"
#include "edg1hole.h"
#include "render.h"


point3f RendPoints[10] =
  {{0.0,0.0,-1.0},
   {0.0,0.0,-4.0},{0.5,0.0,-4.0},{0.0,0.5,-4.0},
   {0.0,0.0,-4.0},{0.5,0.0,-4.0},{0.0,0.5,-4.0},
   {0.0,0.0,0.0},{1.0,0.0,0.0},{0.0,1.0,0.0}};

vector3f lightdir;
point3f  rp0;
vector3f rv1, rv2;

XImage  *theimage;
boolean RenderingIsOn = false;

/* rendering switches */
boolean rswGaussian = false;
boolean rswMean     = false;
boolean rswSections = false;
boolean rswVDepRefl = false;
boolean rswVIndRefl = false;
boolean rswChess    = false;
boolean eswSections = false;
boolean eswLines1   = false;
boolean eswLines2   = false;
boolean eswLight    = false;
boolean eswEdRendering = false;

static BezPatchTreefp patchtree[5*GH_MAX_K];
static int nptrees;
static int lastline;

static char *imagedata;
static boolean (*RefLines)( vector3f *v, point3f *p, vector3f *n );

static struct {
         unsigned short r_bits,  g_bits,  b_bits;
         unsigned char  r_shift, g_shift, b_shift;
       } rgbmap;

/* ///////////////////////////////////////////////////////////////////////// */
static void parse_colourmask ( int mask, unsigned short *bits, unsigned char *shift )
{
  unsigned char sh;

  if ( !mask )
    *bits = *shift = 0;
  else {
    for ( sh = 0;  !(mask & 0x01);  mask = mask >> 1, sh++ )
      ;
    *shift = sh;
    *bits = mask;
  }
} /*parse_colourmask*/

void alloc_image ()
{
  int nplanes;

  nplanes = XDisplayPlanes ( thedisplay, thescreen );
  parse_colourmask ( thevisual->red_mask,   &rgbmap.r_bits, &rgbmap.r_shift );
  parse_colourmask ( thevisual->green_mask, &rgbmap.g_bits, &rgbmap.g_shift );
  parse_colourmask ( thevisual->blue_mask,  &rgbmap.b_bits, &rgbmap.b_shift );
/*
  printf ( "  red bits = %d, shift = %d\n", rgbmap.r_bits, rgbmap.r_shift );
  printf ( "green bits = %d, shift = %d\n", rgbmap.g_bits, rgbmap.g_shift );
  printf ( " blue bits = %d, shift = %d\n", rgbmap.b_bits, rgbmap.b_shift );
  printf ( "depth = %d\n", nplanes );
*/
  theimage = XCreateImage ( thedisplay, thevisual, nplanes, ZPixmap, 0,
                            NULL, MAX_WIDTH, MAX_HEIGHT, 8, 0 );
  if ( !theimage )
    exit ( 1 );   
  imagedata = (char*)malloc ( MAX_HEIGHT*theimage->bytes_per_line );
  theimage->data = imagedata;
/*
  printf ( "theimage = %x\n", theimage );
  printf ( "r: %x,  g: %x,  b:%x\n",
           thevisual->red_mask, thevisual->green_mask, thevisual->blue_mask );
*/
} /*alloc_image*/

int PixelColour ( float r, float g, float b )
{
  int pix;

  r = r < 0.0 ? 0.0 : ( r > 1.0 ? 1.0 : r );
  g = g < 0.0 ? 0.0 : ( g > 1.0 ? 1.0 : g );
  b = b < 0.0 ? 0.0 : ( b > 1.0 ? 1.0 : b );
  pix = (((int)(r*(float)rgbmap.r_bits)) << rgbmap.r_shift) +
        (((int)(g*(float)rgbmap.g_bits)) << rgbmap.g_shift) +
        (((int)(b*(float)rgbmap.b_bits)) << rgbmap.b_shift); 
  return pix;
} /*PixelColour*/

boolean InitRenderer ()
{
  alloc_image ();
  memset ( patchtree, 0, 4*GH_MAX_K*sizeof(BezPatchTreefp) );
  nptrees = 0;
  return true;
} /*InitRenderer*/

void DestroyRenderer ()
{
} /*DestroyRenderer*/

boolean ResetRenderer ()
{
  int i;

  for ( i = 0; i < 4*GH_MAX_K; i++ )
    if ( patchtree[i] ) {
      rbez_DestroyBezPatchTreef ( patchtree[i] );
      patchtree[i] = NULL;
    }
  nptrees = 0;
  return true;
} /*ResetRenderer*/

boolean RendEnterBezierPatch ( int n, int m, point3f *cp )
{
  if ( nptrees < 4*GH_MAX_K ) {
    patchtree[nptrees] = rbez_NewBezPatchTreef ( 0, n, m, 0.0, 1.0, 0.0, 1.0, cp );
    if ( !patchtree[nptrees] )
      return false;
    nptrees ++;
    return true;
  }
  else
    return false;
} /*RendEnterBezierPatch*/

static float *GetKnotSequencef ( int i )
{
  if ( i < 0 ) i += hole_k;
  else if ( i >= hole_k ) i -= hole_k;

  return knots[i];
} /*GetKnotSequencef*/

/* ////////////////////////////////////////////////////////////////////////// */
static float minGauss, maxGauss, minMean, maxMean;

void FindCurvatures ()
{
#define DENS 20
  int     i, j, k, n, m;
  float   u, v, Mean, Gauss;
  point3f *cp;
  boolean ok;

  ok = false;
  for ( i = 0; i < nptrees; i++ ) {
    n = patchtree[i]->n;
    m = patchtree[i]->m;
    cp = patchtree[i]->root->ctlpoints;
    if ( !ok ) {
      mbs_GMCurvaturesBP3f ( n, m, cp, 0.0, 0.0, &minGauss, &minMean );
      maxGauss = minGauss;
      maxMean = minMean;
      ok = true;
    }
    for ( j = 0; j <= DENS; j++ ) {
      u = (float)j/(float)DENS;
      for ( k = 0; k <= DENS; k++ ) {
        v = (float)k/(float)DENS;
        mbs_GMCurvaturesBP3f ( n, m, cp, u, v, &Gauss, &Mean );
        if ( Gauss < minGauss ) minGauss = Gauss;
        else if ( Gauss > maxGauss ) maxGauss = Gauss;
        if ( Mean < minMean ) minMean = Mean;
        else if ( Mean > maxMean ) maxMean = Mean;
      }
    }
  }
  printf ( "Gaussian: %f, %f\n", minGauss, maxGauss );
  printf ( "    Mean: %f, %f\n", minMean, maxMean );
#undef DENS
} /*FindCurvatures*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean RefLines0 ( vector3f *v, point3f *p, vector3f *n )
{
  return true;
} /*RefLines0*/

boolean RefLines1 ( vector3f *v, point3f *p, vector3f *n )
{
/* viewer-independent */
  vector3f dp;
  double   s;

  SubtractPoints3f ( p, &rp0, &dp );
  s = det3f ( n, &dp, &rv2 ) / det3f ( n, &rv1, &rv2 );
  if ( s < 0.0 )
    s = 1.0-s;
  return ((int)s) & 0x01;
} /*RefLines1*/

boolean RefLines2 ( vector3f *v, point3f *p, vector3f *n )
{
/* viewer-dependent */
  vector3f r, dp;
  double   s;

  s = 2.0*DotProduct3f( n, v )/DotProduct3f ( n, n );
  SetVector3f ( &r, v->x - s*n->x, v->y - s*n->y, v->z - s*n->z );
  SubtractPoints3f ( p, &rp0, &dp );
  s = det3f ( &r, &dp, &rv2 ) / det3f ( &r, &rv1, &rv2 );
  if ( s < 0.0 )
    s = 1.0-s;
  return !(((int)s) & 0x01);
} /*RefLines2*/

boolean RefLines3 ( vector3f *v, point3f *p, vector3f *n )
{
/* sections */
  vector3f r;
  double   s;

  SubtractPoints3f ( p, &RendPoints[7], &r );
  s = DotProduct3f( &r, &rv1 );
  if ( s < 0.0 )
    s = 1.0-s;
  return (((int)s) & 0x01);
} /*RefLines3*/

boolean RefLines1c ( vector3f *v, point3f *p, vector3f *n )
{
/* viewer-independent, chess */
  vector3f dp;
  double   d, d1, d2;

  SubtractPoints3f ( p, &rp0, &dp );
  d  = det3f ( n, &rv1, &rv2 );
  d1 = det3f ( n, &dp, &rv2 )/d;
  d2 = det3f ( n, &rv1, &dp )/d;
  if ( d1 < 0.0 )
    d1 = 1.0-d1;
  if ( d2 < 0.0 )
    d2 = 1.0-d2;

  return (((int)d1) & 0x01) ^ (((int)d2) & 0x01);
} /*RefLines1c*/

boolean RefLines2c ( vector3f *v, point3f *p, vector3f *n )
{
/* viewer-dependent, chess */
  vector3f r, dp;
  double   s, d1, d2, d;

  s = 2.0*DotProduct3f( n, v )/DotProduct3f ( n, n );
  SetVector3f ( &r, v->x - s*n->x, v->y - s*n->y, v->z - s*n->z );
  SubtractPoints3f ( p, &rp0, &dp );
  d  = det3f ( &r, &rv1, &rv2 );
  d1 = det3f ( &r, &dp, &rv2 )/d;
  d2 = det3f ( &r, &rv1, &dp )/d;
  if ( d1 < 0.0 )
    d1 = 1.0-d1;
  if ( d2 < 0.0 )
    d2 = 1.0-d2;

  return (((int)d1) & 0x01) ^ (((int)d2) & 0x01);
} /*RefLines2c*/

boolean RefLines3c ( vector3f *v, point3f *p, vector3f *n )
{
/* sections, chess */
  vector3f r;
  double   s1, s2;

  SubtractPoints3f ( p, &RendPoints[7], &r );
  s1 = DotProduct3f ( &r, &rv1 );
  if ( s1 < 0.0 )
    s1 = 1.0-s1;
  s2 = DotProduct3f ( &r, &rv2 );
  if ( s2 < 0.0 )
    s2 = 1.0-s2;

  return (((int)s1) & 0x01) ^ (((int)s2) & 0x01);
} /*RefLines3c*/

void PhongLight ( vector3f *tex, point3f *p, vector3f *n, vector3f *v, vector3f *light,
                  float *r, float *g, float *b )
{
  float    a, c, d, f;
  vector3f h;

  NormalizeVector3f ( n );
  c = -DotProduct3f ( n, light );
  d = DotProduct3f ( n, v );
  if ( c*d < 0.0 )
    c = d = 0.0;
  else {
    MidPoint3f ( v, light, &h );
    NormalizeVector3f ( &h );
    d = fabs ( DotProduct3f ( n, &h ) );
    d = d*d;  d = d*d;  d = d*d;  d = d*d;  d = d*d;  d = d*d;
    c = fabs(c);
  }
  f = 0.3;
  if ( !RefLines ( v, p, n ) )
    { c *= 0.8;  f *= 0.8; }
  a = ((f+0.65*c)*tex->x + 0.15*d);  *r = min ( a, 1.0 );
  a = ((f+0.65*c)*tex->y + 0.15*d);  *g = min ( a, 1.0 );
  a = ((f+0.65*c)*tex->z + 0.15*d);  *b = min ( a, 1.0 );
} /*PhongLight*/

void GetTexColour ( int n, int m, const point3f *cp, float u, float v,
                    vector3f *texcolour )
{
  float Mean, Gauss, t, s;

  if ( rswMean ) {
    mbs_GMCurvaturesBP3f ( n, m, cp, u, v, &Gauss, &Mean );
    t = (Mean-minMean)/(maxMean-minMean);
    t = max ( 0.0, t );
    t = min ( 1.0, t );
    s = 1.0-t;
    texcolour->x = t+t*s;
    texcolour->y = s+t*s;
    texcolour->z = 0.0;
  }
  else if ( rswGaussian ) {
    mbs_GMCurvaturesBP3f ( n, m, cp, u, v, &Gauss, &Mean );
    t = (Gauss-minGauss)/(maxGauss-minGauss);
    t = max ( 0.0, t );
    t = min ( 1.0, t );
    s = 1.0-t;
    texcolour->x = texcolour->y = t+t*s;
    texcolour->z = s+t*s;
  }
  else {
    texcolour->x = texcolour->y = 1.0;
    texcolour->z = 0.25;
  }
} /*GetTexColour*/

void GetPixelColour ( ray3f *ray, float *r, float *g, float *b )
{
#define MAXLEVEL 15
#define MAXINTERS 10
  int      i, j, k = 0, ninters, ii;
  RayObjectIntersf inters[MAXINTERS], iint;
  vector3f texcolour;

  *r = *g = *b = 0.67;
  ii = 0;
  for ( i = 0; i < nptrees; i++ ) {
    ninters = rbez_FindRayBezPatchIntersf ( patchtree[i], ray,
       MAXLEVEL, MAXINTERS, &ninters, inters );
    if ( ninters ) {
      if ( !ii ) {
        iint = inters[0];
        ii = 1;  k = i;
      }
      for ( j = 0; j < ninters; j++ )
        if ( inters[j].t < iint.t ) {
          iint = inters[j];
          k = i;
        }
    }
  }
  if ( ii ) {
    GetTexColour ( patchtree[k]->n, patchtree[k]->m,
                   patchtree[k]->root->ctlpoints,
                   iint.u, iint.v, &texcolour );
    PhongLight ( &texcolour, &iint.p, &iint.nv, &ray->v, &lightdir, r, g, b );
  }

#undef MAXLEVEL
#undef MAXINTERS
} /*GetPixelColour*/

void SetupTextureProc ()
{
  RefLines = RefLines0;
  if ( rswChess ) {
    if ( rswVDepRefl )
      RefLines = RefLines2c;
    else if ( rswVIndRefl )
      RefLines = RefLines1c;
    else if ( rswSections )
      RefLines = RefLines3c;
  }
  else {
    if ( rswVDepRefl )
      RefLines = RefLines2;
    else if ( rswVIndRefl )
      RefLines = RefLines1;
    else if ( rswSections )
      RefLines = RefLines3;
  }
} /*SetupTextureProc*/

void InitRendVectors ()
{
  int i;

  lightdir = RendPoints[0];
  if ( rswSections ) i = 7;
  else if ( rswVDepRefl ) i = 1;
  else if ( rswVIndRefl ) i = 4;
  else return;

  rp0 = RendPoints[i];
  SubtractPoints3f ( &RendPoints[i+1], &RendPoints[i], &rv1 );
  SubtractPoints3f ( &RendPoints[i+2], &RendPoints[i], &rv2 );
} /*InitRendVectors*/

boolean BeginRendering ( ed_rect *er )
{
  void    *sp;
  int     i, j, k;
  int     *ind;
  point3f *cp, *cq;
  float   *ukn, *vkn;

  sp = pkv_GetScratchMemTop ();
  if ( swDisplaySurfPatches ) {
    ind = pkv_GetScratchMem ( 16*sizeof(int) );
    cp  = pkv_GetScratchMem ( 32*sizeof(point3f) );
    if ( ind && cp ) {
      cq = &cp[16];
      for ( i = 0; i < hole_k; i++ )
        for ( j = 0; j < 3; j++ ) {
          ukn = GetKnotSequencef ( i-1 );  ukn += 3;
          vkn = GetKnotSequencef ( i );    vkn += j;
          gh_GetBspInd ( hole_k, i, j, ind );
          for ( k = 0; k < 16; k++ )
            cp[k] = surfcp[ind[k]];
          mbs_BSPatchToBezf ( 3, 3, 7, ukn, 3, 7, vkn, 12, (float*)cp,
                        NULL, NULL, NULL, NULL, NULL, NULL, 12, (float*)cq );
          RendEnterBezierPatch ( 3, 3, cq );
        }
    }
    else
      goto failure;
  }
  if ( swDisplayFinalPatches && FinalSurfValid ) {
    for ( i = 0; i < hole_k; i++ )
      RendEnterBezierPatch ( G1H_FINALDEG, G1H_FINALDEG, FinalCP[i] );
  }
  if ( swDisplayNLFinalPatches && NLFinalSurfValid ) {
    for ( i = hole_k; i < 2*hole_k; i++ )
      RendEnterBezierPatch ( G1H_FINALDEG, G1H_FINALDEG, FinalCP[i] );
  }

  if ( nptrees > 0 ) {
    XGetSubImage ( thedisplay, thepixmap, er->x, er->y, er->w, er->h,
                   0xFFFFFFFF, ZPixmap, theimage, 0, 0 );
    FindCurvatures ();
    SetupTextureProc ();
    InitRendVectors ();
    lastline = -1;
    RenderingIsOn = true;
    pkv_SetScratchMemTop ( sp );
    return true;
  }

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*BeginRendering*/

void StopRendering ()
{
  RenderingIsOn = false;
} /*StopRendering*/

void RenderLine ()
{
  int   x;
  ray3f ray;
  float r, g, b;

  if ( RenderingIsOn ) {
    lastline ++;
    for ( x = 0; x < CPos.width; x++ ) {
      CameraRayOfPixelf ( &CPos, CPos.xmin+x, CPos.ymin+lastline, &ray );
      GetPixelColour ( &ray, &r, &g, &b );
      XPutPixel ( theimage, x, lastline, PixelColour ( r, g, b ) );
    }
    if ( lastline == CPos.height-1 )
      RenderingIsOn = false;
  }
/*
printf ( "%3d", lastline );
*/
} /*RenderLine*/

/* ///////////////////////////////////////////////////////////////////////// */
void RedrawRendPoints ( int id )
{
  int i0, i1, i;
  point3f p, q, r;

  if ( eswLight ) {
    SetPoint3f ( &p, 0.0, 0.0, 0.0 );
    if ( id < 3 ) {
      CameraProjectPoint3f ( &PPos[id], &p, &q );
      CameraProjectPoint3f ( &PPos[id], &RendPoints[0], &r );
    }
    else {
      CameraProjectPoint3f ( &CPos, &p, &q );
      CameraProjectPoint3f ( &CPos, &RendPoints[0], &r );
    }
    XSetForeground ( thedisplay, thegc, c_lt_green );
    MyDrawLine ( (point2f*)(void*)&q, (point2f*)(void*)&r );
    XSetForeground ( thedisplay, thegc, c_white );
    MyMarkRect ( (point2f*)(void*)&q );
    XSetForeground ( thedisplay, thegc, c_lt_red );
    MyMarkRect ( (point2f*)(void*)&r );
    return;
  }
  else if ( eswSections )
    { i0 = 7;  i1 = 9; }
  else if ( eswLines1 )
    { i0 = 1;  i1 = 3; }
  else if ( eswLines2 )
    { i0 = 4;  i1 = 6; }
  else return;

  if ( id < 3 )
    CameraProjectPoint3f ( &PPos[id], &RendPoints[i0], &q );
  else
    CameraProjectPoint3f ( &CPos, &RendPoints[i0], &q );
  for ( i = i0; i <= i1; i++ ) {
    if ( id < 3 )
      CameraProjectPoint3f ( &PPos[id], &RendPoints[i], &r );
    else
      CameraProjectPoint3f ( &CPos, &RendPoints[i], &r );
    XSetForeground ( thedisplay, thegc, c_lt_green );
    MyDrawLine ( (point2f*)(void*)&q, (point2f*)(void*)&r );
    XSetForeground ( thedisplay, thegc, c_lt_red );
    MyMarkRect ( (point2f*)(void*)&r );
  }
  MyMarkRect ( (point2f*)(void*)&q );
} /*RedrawRendPoints*/

int FindNearestRendPoint ( int id, int x, int y )
{
#define TOL 10.0
  int     i0, i1, i, j;
  float   d, e;
  point3f q;

  if ( eswLight )         { i0 = i1 = 0; }
  else if ( eswSections ) { i0 = 7;  i1 = 9; }
  else if ( eswLines1 )   { i0 = 1;  i1 = 3; }
  else if ( eswLines2 )   { i0 = 4;  i1 = 6; }
  else return -1;

  d = TOL+1.0;
  j = -1;
  for ( i = i0; i <= i1; i++ ) {
    CameraProjectPoint3f ( &PPos[id], &RendPoints[i], &q );
    e = fabs((float)x-q.x) + fabs((float)y-q.y);
    if ( e < d ) {
      j = i;
      d = e;
    }
  }
 
  if ( d > TOL )
    j = -1;
  return j;
#undef TOL
} /*FindNearestRendPoint*/

void SetRendPoint ( int id, int np, int x, int y )
{
  point3f  q, r;
  vector3f v;

  CameraProjectPoint3f ( &PPos[id], &RendPoints[np], &q );
  q.x = (float)x;
  q.y = (float)y;
  CameraUnProjectPoint3f ( &PPos[id], &q, &r );
  SubtractPoints3f ( &r, &RendPoints[np], &v );
  RendPoints[np] = r;

  if ( np == 0 )
    NormalizeVector3f ( &RendPoints[0] );
  else if ( eswSections && np == 7 ) {
    AddVector3f ( &RendPoints[8], &v, &RendPoints[8] );
    AddVector3f ( &RendPoints[9], &v, &RendPoints[9] );
  }
  else if ( eswLines1 && np == 1 ) {
    AddVector3f ( &RendPoints[2], &v, &RendPoints[2] );
    AddVector3f ( &RendPoints[3], &v, &RendPoints[3] );
  }
  else if ( eswLines2 && np == 4 ) {
    AddVector3f ( &RendPoints[5], &v, &RendPoints[5] );
    AddVector3f ( &RendPoints[6], &v, &RendPoints[6] );
  }
} /*SetRendPoint*/

void ResetRendLines ()
{
  point3f DefaultRendPoints[10] =
    {{0.0,0.0,-1.0},
     {0.0,0.0,-4.0},{0.5,0.0,-4.0},{0.0,0.5,-4.0},
     {0.0,0.0,-4.0},{0.5,0.0,-4.0},{0.0,0.5,-4.0},
     {0.0,0.0,0.0},{1.0,0.0,0.0},{0.0,1.0,0.0}};

  memcpy ( RendPoints, DefaultRendPoints, 10*sizeof(point3f) );
} /*ResetRendLines*/

