
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "camera.h"
#include "multibs.h"

#include "xgedit.h"
#include "bookg1holef.h"
#include "datagen.h"
#include "drawg1hole.h"
#include "drawbez.h"


xge_3Dwinf swind;

point2f rBpt[4][MAX_BPTS];   /* BS ctl points projections in windows */

point3f savedcp[MAX_BPTS];
boolean cpmark[MAX_BPTS];

/* items to display */

point3f Bezpt[MAX_K][3][16]; /* Bezier patches surrounding the hole, bicubic */
point3f TBezpt[MAX_K][36];   /* target Bezier patches, biquintic */

int     nauxp;
point3f AuxBezpt[3*MAX_K][8];/* auxiliary patches, of degree (1,3) */
int     nauxc;
point3f AuxCurvpt[MAX_K][4]; /* auxiliary curves, of degree 3 */
float   auxct;
#define central_point Bpt[hole_np]
point3f StarCurvept[MAX_K][4]; /* final filling patches boundary curves */

float beta1=0.0, beta2=0.0;    /* filling procedure parameters */
boolean autocpoint = true;


/* editor */
boolean cpoint_ok = false;     /* manually set central point initialized */


/* displaying switches */
boolean bsnet = true;
boolean bezp = true;
boolean auxcurv = false;
boolean starcurv = false;
boolean auxpatch = false;
boolean targpatch = true;
boolean filling_ok, auxcurv_ok, centp_ok, starcurv_ok, auxpatch_ok;


/* ///////////////////////////////////////////////////////////////////////// */
static void GetSurroundingPatch ( int i, int j, point3f *bcp )
{
  int     k;
  int     *ind;
  point3f *q;
  void    *sp;
  float   kn[8] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 6.0};

  sp  = pkv_GetScratchMemTop ();
  ind = pkv_GetScratchMem ( 16*sizeof(int) );
  q   = pkv_GetScratchMem ( 16*sizeof(point3f) );
  if ( !ind || !q )
    exit ( 1 );

  GetBspInd ( i, j, ind );
  for ( k = 0; k < 16; k++ )
    q[k] = Bpt[ind[k]];
  mbs_BSPatchToBezf ( 3, 3, 7, kn, 3, 7, kn, 12, (float*)q,
                      NULL, NULL, NULL, NULL, NULL, NULL, 12, (float*)bcp );

  pkv_SetScratchMemTop ( sp );
} /*GetSurroundingPatch*/

static point3f* GetSurrPatchCP ( int i, int j )
{
  return &Bezpt[i][j][0];
} /*GetSurrPatchCP*/

static void OutAuxPatches ( int npatches, int degu, int degv,
                            const point3f *apcp )
{
  if ( degu != 1 || degv != 3 )
    auxpatch_ok = false;
  else {
    if ( !auxpatch_ok )
      nauxp = 0;
    memcpy ( &AuxBezpt[nauxp], apcp, npatches*8*sizeof(point3f) );
    nauxp += npatches;
    auxpatch_ok = true;
  }
} /*OutAuxPatches*/

static void OutAuxCurves ( int ncurves, int degree, const point3f *accp,
                           float t )
{
  if ( degree != 3 )
    auxcurv_ok = false;
  else {
    nauxc = ncurves;
    memcpy ( AuxCurvpt, accp, ncurves*4*sizeof(point3f) );
    auxct = t;
    auxcurv_ok = true;
  }
} /*OutAuxCurves*/

static void OutCentralPoint ( point3f *p )
{
  if ( autocpoint || !cpoint_ok ) {
    central_point = *p;
    cpoint_ok = true;
  }
  else if ( cpoint_ok )
    *p = central_point;
  centp_ok = true;
} /*OutCentralPoint*/

static void OutStarCurves ( int ncurves, int degree, const point3f *sccp )
{
  if ( degree != 3 )
    starcurv_ok = false;
  else {
    memcpy ( StarCurvept, sccp, ncurves*4*sizeof(point3f) );
    starcurv_ok = true;
  }
} /*OutStarCurves*/

void UpdateSurface ( void )
{
  int i, j;

  for ( i = 0; i < hole_k; i++ )
    for ( j = 0; j < 3; j++ )
      GetSurroundingPatch ( i, j, &Bezpt[i][j][0] );

  G1OutCentralPointf = OutCentralPoint;
  G1OutAuxCurvesf    = OutAuxCurves;
  G1OutStarCurvesf   = OutStarCurves;
  G1OutAuxPatchesf   = OutAuxPatches;
  auxcurv_ok = centp_ok = starcurv_ok = auxpatch_ok = false;
  filling_ok = FillG1Holef ( hole_k, GetSurrPatchCP, (float)(1.0+beta1),
                             (float)(0.5+beta2), &TBezpt[0][0] );
} /*UpdateSurface*/

/* ///////////////////////////////////////////////////////////////////////// */
void SelectCPoints ( int id, short x0, short x1, short y0, short y1 )
{
  int i;

  id &= 0x03;
  if ( x1-x0 >= 2 || y1-y0 >= 2 ) {
    for ( i = 0; i <= 12*hole_k; i++ )
      if ( rBpt[id][i].x >= x0 && rBpt[id][i].x <= x1 &&
           rBpt[id][i].y >= y0 && rBpt[id][i].y <= y1 )
        cpmark[i] = true;
  }
  else {
    if ( FindNearestPoint ( id, x0, y0, xge_MINDIST ) )
      cpmark[swind.current_point] = true;
  }
} /*SelectCPoints*/

void UnSelectCPoints ( int id, short x0, short x1, short y0, short y1 )
{
  int i;

  id &= 0x03;
  if ( x1-x0 >= 2 || y1-y0 >= 2 ) {
    for ( i = 0; i <= 12*hole_k; i++ )
      if ( rBpt[id][i].x >= x0 && rBpt[id][i].x <= x1 &&
           rBpt[id][i].y >= y0 && rBpt[id][i].y <= y1 )
        cpmark[i] = false;
  }
  else {
    if ( FindNearestPoint ( id, x0, y0, xge_MINDIST ) )
      cpmark[swind.current_point] = false;
  }
} /*UnSelectCPoints*/

void SaveCPoints ( void )
{
  memcpy ( savedcp, Bpt, (12*hole_k+1)*sizeof(point3f) );
} /*SaveCPoints*/

void TransformCPoints ( trans3f *tr )
{
  int i;

  for ( i = 0; i <= 12*hole_k; i++ )
    if ( cpmark[i] )
      TransPoint3f ( tr, &savedcp[i], &Bpt[i] );
  UpdateSurface ();
  for ( i = 0; i < 4; i++ )
    ProjectScene ( i );
} /*TransformCPoints*/

/* ///////////////////////////////////////////////////////////////////////// */
void FindBoundingBox ( Box3f *box )
{
  int i;

  box->x0 = box->x1 = Bpt[0].x;
  box->y0 = box->y1 = Bpt[0].y;
  box->z0 = box->z1 = Bpt[0].z;
  for ( i = 1; i <= 12*hole_k; i++ ) {
    if ( Bpt[i].x < box->x0 )      box->x0 = Bpt[i].x;
    else if ( Bpt[i].x > box->x1 ) box->x1 = Bpt[i].x;
    if ( Bpt[i].y < box->y0 )      box->y0 = Bpt[i].y;
    else if ( Bpt[i].y > box->y1 ) box->y1 = Bpt[i].y;
    if ( Bpt[i].z < box->z0 )      box->z0 = Bpt[i].z;
    else if ( Bpt[i].z > box->z1 ) box->z1 = Bpt[i].z;
  }
} /*FindBoundingBox*/

/* ///////////////////////////////////////////////////////////////////////// */
void ResizeObject ( void )
{
  int i;

  for ( i = 0; i < 4; i++ )
    ProjectScene ( i );
} /*ResizeObject*/

void ResetObject ( void )
{
  dataparam[0] = dataparam[1] = dataparam[2] = 0.0;
  beta1 = 0.0;
  beta2 = 0.5;
  InitHole ( hole_k, dataparam[0], dataparam[1], dataparam[2] );
  memset ( cpmark, 0, (12*hole_k+1)*sizeof(boolean) );
  UpdateSurface ();
  ResizeObject ();
} /*ResetObject*/

boolean FindNearestPoint ( int id, int x, int y, int mindist )
{
  int   i, n0, n1;
  float d, e;

  id &= 0x03;
  n0 = 0;
  if ( !bsnet ) n0 = hole_np;
  n1 = hole_np;
  if ( centp_ok && !autocpoint ) n1++;
  e = (float)mindist;
  for ( i = n0; i < n1; i++ ) {
    d = (float)(fabs((float)x-rBpt[id][i].x) + fabs((float)y-rBpt[id][i].y));
    if ( d < e ) {
      swind.current_point = i;
      e = d;
    }
  }
  return (boolean)(e < mindist);
} /*FindNearestPoint*/

void SetCPoint ( int id, int x, int y )
{
  int      j, n;
  point3f  r;

  id &= 0x03;
  j = swind.current_point;
  n = hole_np;
  if ( cpoint_ok && !autocpoint )
    n++;
  if ( swind.current_point < n ) {
    CameraProjectPoint3f ( &swind.CPos[id], &Bpt[j], &r );
    r.x = (float)x;  r.y = (float)y;
    CameraUnProjectPoint3f ( &swind.CPos[id], &r, &Bpt[j] );
  }
  UpdateSurface ();
  ResizeObject ();
} /*SetCPoint*/

static void RzutujPunkt ( int id, point3f *p, point2f* q )
{
  point3f qq;

  id &= 0x03;
  CameraProjectPoint3f ( &swind.CPos[id], p, &qq );
  q->x = qq.x;  q->y = qq.y;
} /*RzutujPunkt*/

void RzutujXPunkt ( int id, point3f *p, XPoint* q )
{
  point3f qq;

  id &= 0x03;
  CameraProjectPoint3f ( &swind.CPos[id], p, &qq );
  q->x = (short)qq.x;  q->y = (short)qq.y;
} /*RzutujXPunkt*/

void RzutujPK ( int id, const point3f *p, point3f *q )
{
  id &= 0x03;
  TransPoint3f ( &swind.CPos[id].CTr, p, q );
  if ( id < 3 )
    q->z = 1.0;
  else
    SetPoint3f ( q, q->x+swind.CPos[id].vd.persp.xi0*q->z,
                 q->y+swind.CPos[id].vd.persp.eta0*q->z, q->z );
} /*RzutujPK*/

void ProjectScene ( int id )
{
  int j, n;

  id &= 0x03;
  n = hole_np;
  if ( cpoint_ok )
    n++;
  for ( j = 0; j < n; j++ )
    RzutujPunkt ( id, &Bpt[j], &rBpt[id][j] );
} /*ProjectScene*/

static void RectMark ( int x, int y )
{
  xgeFillRectangle ( 3, 3, x-1, y-1 );
} /*RectMark*/

static void DrawBSNet ( int i, point2f *pts )
{
  int     ind[16];
  int     j;
  int     x, y;

  xgeSetForeground ( xgec_Green );
  GetBspInd ( i, 0, ind );
  for ( i = 0; i <= 3; i++ )
    for ( j = 0; j < 3; j++ )
      xgeDrawLine ( pts[ind[4*i+j]].x, pts[ind[4*i+j]].y,
                    pts[ind[4*i+j+1]].x, pts[ind[4*i+j+1]].y );
  for ( j = 1; j <= 3; j++ )
    for ( i = 0; i < 3; i++ )
      xgeDrawLine ( pts[ind[4*i+j]].x, pts[ind[4*i+j]].y,
                    pts[ind[4*i+j+4]].x, pts[ind[4*i+j+4]].y );
  xgeSetForeground ( xgec_Yellow );
  for ( j = 0; j < hole_np; j++ )
    if ( !cpmark[j] ) {
      x = (int)pts[j].x; 
      y = (int)pts[j].y;
      RectMark ( x, y );
    }
  xgeSetForeground ( xgec_OrangeRed );
  for ( j = 0; j < hole_np; j++ )
    if ( cpmark[j] ) {
      x = (int)pts[j].x; 
      y = (int)pts[j].y;
      RectMark ( x, y );
    }
} /*DrawBSNet*/

void DrawBSNets ( int id )
{
  int i;

  id &= 0x03;
  for ( i = 0; i < hole_k; i++ )
    DrawBSNet ( i, rBpt[id] );
} /*DrawBSNets*/

void DrawBezPatches ( int id )
{
  int i, j;

  id &= 0x03;
  if ( id == 3 ) {
    for ( i = 0; i < hole_k; i++ )
      for ( j = 0; j < 3; j++ )
        DrawBezPatch3a ( &swind.CPos[id], 3, 3, Bezpt[i][j], 6, 6,
                         xgec_White, xgec_Grey3 );
  }
  else {
    for ( i = 0; i < hole_k; i++ )
      for ( j = 0; j < 3; j++ )
        DrawBezPatch3b ( &swind.CPos[id], 3, 3, Bezpt[i][j], 6, 6,
                         xgec_White, xgec_Grey3 );
  }
} /*DrawBezPatches*/

void DrawAuxCurves ( int id )
{
  int     i;
  point3f p;
  XPoint  q;

  id &= 0x03;
  if ( auxcurv_ok ) {
    if ( id == 3) {
      for ( i = 0; i < nauxc; i++ ) {
        DrawBezCurve3a ( &swind.CPos[id], 3, AuxCurvpt[i], xgec_DodgerBlue );
        mbs_BCHornerC3f ( 3, AuxCurvpt[i], auxct, &p );
        RzutujXPunkt ( id, &p, &q );
        RectMark ( q.x, q.y );
      }
    }
    else {
      for ( i = 0; i < nauxc; i++ ) {
        DrawBezCurve3b ( &swind.CPos[id], 3, AuxCurvpt[i], xgec_DodgerBlue );
        mbs_BCHornerC3f ( 3, AuxCurvpt[i], auxct, &p );
        RzutujXPunkt ( id, &p, &q );
        RectMark ( q.x, q.y );
      }
    }
  }
} /*DrawAuxCurves*/

void DrawStarCurves ( int id )
{
  int i;

  id &= 0x03;
  if ( starcurv_ok ) {
    if ( id == 3 ) {
      for ( i = 0; i < hole_k; i++ )
        DrawBezCurve3a ( &swind.CPos[id], 3, StarCurvept[i], xgec_Yellow );
    }
    else {
      for ( i = 0; i < hole_k; i++ )
        DrawBezCurve3b ( &swind.CPos[id], 3, StarCurvept[i], xgec_Yellow );
    }
  }
} /*DrawStarCurves*/

void DrawAuxPatches ( int id )
{
  int i;

  id &= 0x03;
  if ( auxpatch_ok ) {
    if ( id == 3 ) {
      for ( i = 0; i < nauxp; i++ )
        DrawBezPatch3a ( &swind.CPos[id], 1, 3, AuxBezpt[i], 1, 6,
                         xgec_Yellow, xgec_Red );
    }
    else {
      for ( i = 0; i < nauxp; i++ )
        DrawBezPatch3b ( &swind.CPos[id], 1, 3, AuxBezpt[i], 1, 6,
                         xgec_Yellow, xgec_Red );
    }
  }
} /*DrawAuxPatches*/

void DrawTargetPatches ( int id )
{
  int i;

  id &= 0x03;
  if ( filling_ok ) {
    if ( id == 3 ) {
      for ( i = 0; i < hole_k; i++ )
        DrawBezPatch3a ( &swind.CPos[id], 5, 5, TBezpt[i], 6, 6,
                         xgec_White, xgec_Grey1 );
    }
    else {
      for ( i = 0; i < hole_k; i++ )
        DrawBezPatch3b ( &swind.CPos[id], 5, 5, TBezpt[i], 6, 6,
                         xgec_White, xgec_Grey1 );
    }
  }
} /*DrawTargetPatches*/

static void DrawCentralPoint ( int id )
{
  XPoint q;

  id &= 0x03;
  if ( centp_ok ) {
    RzutujXPunkt ( id, &central_point, &q );
    xgeSetForeground ( xgec_Red );
    RectMark ( q.x, q.y );
  }
} /*DrawCentralPoint*/

void DrawSurface ( int id )
{
  id &= 0x03;
  if ( bsnet )
    DrawBSNets ( id );
  if ( bezp )
    DrawBezPatches ( id );
  if ( auxcurv )
    DrawAuxCurves ( id );
  if ( targpatch )
    DrawTargetPatches ( id );
  if ( starcurv )
    DrawStarCurves ( id );
  if ( auxpatch )
    DrawAuxPatches ( id );
  if ( !autocpoint )
    DrawCentralPoint ( id );
} /*DrawSurface*/

/* ///////////////////////////////////////////////////////////////////////// */
void UstawK3 ( void )
{
  netk3 = true;  netk5 = netk6 = netk8 = false;
  cpoint_ok = false;  autocpoint = true;
  InitHole ( hole_k = 3, dataparam[0], dataparam[1], dataparam[2] );
  memset ( cpmark, 0, (12*hole_k+1)*sizeof(boolean) );
  UpdateSurface ();
  ResizeObject ();
} /*UstawK3*/

void UstawK5 ( void )
{
  netk5 = true;  netk3 = netk6 = netk8 = false;
  cpoint_ok = false;  autocpoint = true;
  InitHole ( hole_k = 5, dataparam[0], dataparam[1], dataparam[2] );
  memset ( cpmark, 0, (12*hole_k+1)*sizeof(boolean) );
  UpdateSurface ();
  ResizeObject ();
} /*UstawK5*/

void UstawK6 ( void )
{
  netk6 = true;  netk3 = netk5 = netk8 = false;
  cpoint_ok = false;  autocpoint = true;
  InitHole ( hole_k = 6, dataparam[0], dataparam[1], dataparam[2] );
  memset ( cpmark, 0, (12*hole_k+1)*sizeof(boolean) );
  UpdateSurface ();
  ResizeObject ();
} /*UstawK6*/

void UstawK8 ( void )
{
  netk8 = true;  netk3 = netk5 = netk6 = false;
  cpoint_ok = false;  autocpoint = true;
  InitHole ( hole_k = 8, dataparam[0], dataparam[1], dataparam[2] );
  memset ( cpmark, 0, (12*hole_k+1)*sizeof(boolean) );
  UpdateSurface ();
  ResizeObject ();
} /*UstawK8*/

void UstawAutoCPoint ( void )
{
  cpoint_ok = false;
  UpdateSurface ();
  ResizeObject ();
} /*UstawAutoCPoint*/

