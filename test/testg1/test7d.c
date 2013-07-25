
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/times.h>
#include <unistd.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camera.h"
#include "psout.h"
#include "multibs.h"  
#include "eg1holed.h"

#include "datagend.h"
#include "bslapd.h"
#include "drawitd.h"

#define _DRAW_BASIS
#define _DRAW_LAPLACIAN
#define DRAW_LAPLACIAN_JUMP
#define _TEST_MATRIX1
#define _TEST_MATRIX2

#define QUAD_FACTOR 10  /* 10 */

#define NK 1 /*3*/
#define M1 2
#define M2 4

#define HOLE_K 5
#define FUNCA21

/*
#define HOLE_K 6
#define FUNCA25
*/
char fn1temp[] = "sbasD%03dd.ps";   /* this is a template for file names! */
char fn2temp[] = "sbasD%03dLd.ps";
char fn3temp[] = "sbasD%03dJd.ps";
char fn3[] = "smatrixd.ps";

GHoleDomaind *domain;
int   nfunc_a, nfunc_b, nfunc_c, nfunc_d;

/* ///////////////////////////////////////////////////////////////////////// */
void SetupPProj ( double w, double h, double x, double y, double diag )
{
  CameraInitFramed ( &PPos, true, true, w, h, x, y, 1.0, 4 );
  CameraInitPosd ( &PPos );
  PPos.vd.para.diag = diag;
  CameraSetMappingd ( &PPos );
} /*SetupPProj*/

void SetupCPos ( double w, double h, double x, double y )
{
  vector3d v;

  CameraInitFramed ( &CPos, false, true, w, h, x, y, 1.0, 4 );
  CameraInitPosd ( &CPos );
  SetVector3d ( &v, 0.0, 0.0, -46.0 );
  CameraMoveGd ( &CPos, &v );
  CameraRotXGd ( &CPos, -0.65*PI );
  CameraRotZGd ( &CPos, 0.05*PI );
  CameraZoomd ( &CPos, 11.0 );
} /*SetupCPos*/

/* ///////////////////////////////////////////////////////////////////////// */
void DrawFN ( double x, double y, int fn )
{
  char s[50];

  ps_Set_Gray ( 0.0 );
  sprintf ( s, "(%d) show", fn );
  ps_MoveTo ( x, y );
  ps_Write_Command ( s );
} /*DrawFN*/

double scf = 1.0;
int   cnt;

void SetScf ( int fn )
{
  if ( fn > nfunc_a+nfunc_b+nfunc_c ) {
     switch ( (fn - (nfunc_a+nfunc_b+nfunc_c)) % 2 ) {
  case 0: scf = 1.0;  break;
  case 1: scf = 2.0;  break;
     }
  }
  else
    scf = 1.0;
} /*SetScf*/

void DrawAuxBSPatch ( int n, int lknu, const double *knu,
                      int m, int lknv, const double *knv, const point3d *cp )
{
  void    *sp;
  point3d *p;
  int     i, s;

  sp = pkv_GetScratchMemTop ();
  s = (lknu-n)*(lknv-m);
  p = pkv_GetScratchMem ( s*sizeof(point3d) );
  if ( !p )
    exit ( 1 );
  memcpy ( p, cp, s*sizeof(point3d) );
  for ( i = 0; i < s; i++ )
    p[i].z *= scf;
  DrawBSPatch3Rd ( n, lknu, knu, m, lknv, knv, p, 0.0, 0.35, 0.0, 1.0, 6, 4 );
  pkv_SetScratchMemTop ( sp );
} /*DrawAuxBSPatch*/

void AltDrawAuxBSPatch ( int n, int lknu, const double *knu,
                         int m, int lknv, const double *knv, const point3d *cp )
{
  void    *sp;
  point3d *p;
  int     i, s;

  if ( !cnt ) {
    sp = pkv_GetScratchMemTop ();
    s = (lknu-n)*(lknv-m);
    p = pkv_GetScratchMem ( s*sizeof(point3d) );
    if ( !p )
      exit ( 1 );
    memcpy ( p, cp, s*sizeof(point3d) );
    for ( i = 0; i < s; i++ )
      p[i].z *= scf;
    DrawBSPatch3Rd ( n, lknu, knu, m, lknv, knv, p, 0.0, 0.35, 0.0, 1.0, 6, 6 );
    pkv_SetScratchMemTop ( sp );
  }
  cnt ++;
} /*AltDrawAuxBSPatch*/

void DrawAuxBasisFuncBSPatch ( int fn, double x, double y )
{
  char s[50];

  SetScf ( fn );
  ps_GSave ();
  sprintf ( s, "%6.2f %6.2f translate", x, y );
  ps_Write_Command ( s );
  ps_Set_Gray ( 0.5 );
  g1h_DrawSplBasAuxPatchesd ( domain, fn, DrawAuxBSPatch );
  ps_GRestore ();
} /*DrawAuxBasisFuncBSPatch*/

void DrawBSPatch ( int n, int lknu, const double *knu,
                   int m, int lknv, const double *knv, const point3d *cp )
{
  void    *sp;
  point3d *p;
  int     i, s;

  sp = pkv_GetScratchMemTop ();
  s = (lknu-n)*(lknv-m);
  p = pkv_GetScratchMem ( s*sizeof(point3d) );
  if ( !p )
    exit ( 1 );
  memcpy ( p, cp, s*sizeof(point3d) );
  for ( i = 0; i < s; i++ )
    p[i].z *= scf;
  DrawBSPatch3d ( n, lknu, knu, m, lknv, knv, p, 6, 6 );
  pkv_SetScratchMemTop ( sp );
} /*DrawBSPatch*/

void EmptyDrawBSPatch ( int n, int lknu, const double *knu,
                        int m, int lknv, const double *knv, const point3d *cp )
{
} /*EmptyDrawBSPatch*/

void DrawBasisFuncBSPatch ( int fn, double x, double y )
{
  char s[50];

  SetScf ( fn );
  ps_GSave ();
  sprintf ( s, "%6.2f %6.2f translate", x, y );
  ps_Write_Command ( s );
  DrawFN ( 0.0, 1100.0, fn );
  ps_Set_Gray ( 0.5 );
  g1h_DrawSplBasFunctiond ( domain, fn, DrawBSPatch );
  if ( fn >= nfunc_a+nfunc_b+nfunc_c ) {
    cnt = -((fn-(nfunc_a+nfunc_b+nfunc_c)) / (nfunc_d/HOLE_K));
    ps_Set_RGB ( 1.0, 0.0, 0.0 );
    g1h_DrawSplBasAuxPatchesd ( domain, fn, AltDrawAuxBSPatch );
  }
  ps_GRestore ();
} /*DrawBasisFuncBSPatch*/

void MakePicture1x ( const char *fn, int ffn, int nf )
{
  ps_OpenFile ( fn, 600 );
  ps_DenseScreen ();
  ps_Write_Command ( "1 setlinecap" );
  ps_Write_Command ( "/Times-Roman findfont 60 scalefont setfont" );
  SetupCPos ( 1600.0, 1160.0, 0.0, 0.0 );

  DrawBasisFuncBSPatch ( ffn+0, 20.0, 3610.0 );
  if ( nf > 1 )
    DrawBasisFuncBSPatch ( ffn+1, 1660.0, 3610.0 );
  if ( nf > 2 )
    DrawBasisFuncBSPatch ( ffn+2, 20.0, 2410.0 );
  if ( nf > 3 )
    DrawBasisFuncBSPatch ( ffn+3, 1660.0, 2410.0 );
  if ( nf > 4 )
    DrawBasisFuncBSPatch ( ffn+4, 20.0, 1210.0 );
  if ( nf > 5 )
    DrawBasisFuncBSPatch ( ffn+5, 1660.0, 1210.0 );
  if ( nf > 6 )
    DrawBasisFuncBSPatch ( ffn+6, 20.0, 10.0 );
  if ( nf > 7 )
    DrawBasisFuncBSPatch ( ffn+7, 1660.0, 10.0 );
  ps_CloseFile ();
  printf ( "%s\n", fn );
} /*MakePicture1x*/

void MakePicture1 ( void )
{
  char fn[40];
  int i, j, n;

  n = nfunc_a+nfunc_b+nfunc_c+nfunc_d;
  for ( i = j = 0; j < n; i++, j += 8 ) {
    sprintf ( fn, fn1temp, i );
    MakePicture1x ( fn, j, min (8, n-j) );
  }
} /*MakePicture1*/

/* ///////////////////////////////////////////////////////////////////////// */
double maxlap;

void DrawBSPatchL1 ( int n, int lknu, const double *knu,
                     int m, int lknv, const double *knv, const point3d *cp )
{
  double lap;

  lap = FindBSPatchMaxLaplaciand ( n, lknu, knu, m, lknv, knv, cp, 20, 20 );
  maxlap = max ( maxlap, lap );
} /*DrawBSPatchL1*/

void DrawBSPatchL2 ( int n, int lknu, const double *knu,
                     int m, int lknv, const double *knv, const point3d *cp )
{
  void    *sp;
  point3d *p;
  int     i, s;

  sp = pkv_GetScratchMemTop ();
  s = (lknu-n)*(lknv-m);
  p = pkv_GetScratchMem ( s*sizeof(point3d) );
  if ( !p )
    exit ( 1 );
  memcpy ( p, cp, s*sizeof(point3d) );
  for ( i = 0; i < s; i++ )
    p[i].z *= scf;
  DrawBSPatchLaplaciand ( n, lknu, knu, m, lknv, knv, p, 6, 6 );
  pkv_SetScratchMemTop ( sp );
} /*DrawBSPatchL2*/

void DrawBezPatchL1 ( int n, int m, const point3d *cp )
{
  double lap;

  lap = FindBezPatchMaxLaplaciand ( n, m, cp, 20, 20 );
  maxlap = max ( maxlap, lap );
} /*DrawBezPatchL1*/

void DrawBezPatchL2 ( int n, int m, const point3d *cp )
{
  void    *sp;
  point3d *p;
  int     i, s;

  sp = pkv_GetScratchMemTop ();
  s = (n+1)*(m+1);
  p = pkv_GetScratchMem ( s*sizeof(point3d) );
  if ( !p )
    exit ( 1 );
  memcpy ( p, cp, s*sizeof(point3d) );
  for ( i = 0; i < s; i++ )
    p[i].z *= scf;
  DrawBezPatchLaplaciand ( n, m, p, 6, 6 );
  pkv_SetScratchMemTop ( sp );
} /*DrawBezPatchL2*/

void DrawBasisFuncBSPatchL ( int fn, double x, double y )
{
  char s[50];

  ps_GSave ();
  sprintf ( s, "%6.2f %6.2f translate", x, y );
  ps_Write_Command ( s );
  ps_Set_Line_Width ( 1.0 );
  ps_Draw_Rect ( CPos.width, CPos.height, CPos.xmin, CPos.ymin );
  ps_Set_Clip_Rect ( CPos.width, CPos.height, CPos.xmin, CPos.ymin );
  DrawFN ( 0, 1100.0, fn );
  ps_Set_Gray ( 0.5 );
  if ( fn < nfunc_a || fn >= nfunc_a+nfunc_b ) {
    maxlap = 0.0;
    g1h_DrawSplBasFunctiond ( domain, fn, DrawBSPatchL1 );
    scf = 0.5/maxlap;
    g1h_DrawSplBasFunctiond ( domain, fn, DrawBSPatchL2 );
  }
  else {
    maxlap = 0.0;
    g1h_DrawBasBFunctiond ( domain, fn-nfunc_a, DrawBezPatchL1 );
    scf = 0.5/maxlap;
    g1h_DrawBasBFunctiond ( domain, fn-nfunc_a, DrawBezPatchL2 );
  }
  ps_GRestore ();
} /*DrawBasisFuncBSPatchL*/

void MakePicture2x ( const char *fn, int ffn, int nf )
{
  ps_OpenFile ( fn, 600 );
  ps_DenseScreen ();
  ps_Write_Command ( "1 setlinecap" );
  ps_Write_Command ( "/Times-Roman findfont 60 scalefont setfont" );
  SetupCPos ( 1600.0, 1160.0, 0.0, 0.0 );

  DrawBasisFuncBSPatchL ( ffn+0, 20.0, 3610.0 );
  if ( nf > 1 )
    DrawBasisFuncBSPatchL ( ffn+1, 1660.0, 3610 );
  if ( nf > 2 )
    DrawBasisFuncBSPatchL ( ffn+2, 20.0, 2410.0 );
  if ( nf > 3 )
    DrawBasisFuncBSPatchL ( ffn+3, 1660.0, 2410.0 );
  if ( nf > 4 )
    DrawBasisFuncBSPatchL ( ffn+4, 20.0, 1210.0 );
  if ( nf > 5 )
    DrawBasisFuncBSPatchL ( ffn+5, 1660.0, 1210.0 );
  if ( nf > 6 )
    DrawBasisFuncBSPatchL ( ffn+6, 20.0, 10.0 );
  if ( nf > 7 )
    DrawBasisFuncBSPatchL ( ffn+7, 1660.0, 10.0 );
  ps_CloseFile ();
  printf ( "%s\n", fn );
} /*MakePicture2x*/

void MakePicture2 ( void )
{
  char fn[40];
  int i, j, n;

  n = nfunc_a+nfunc_b+nfunc_c+nfunc_d;
  for ( i = j = 0; j < n; i++, j += 8 ) {
    sprintf ( fn, fn2temp, i );
    MakePicture2x ( fn, j, min (8, n-j) );
  }
} /*MakePicture2*/

/* ///////////////////////////////////////////////////////////////////////// */
int nbfp;
struct {
    int   n, lkn;
    double knots[40];
    point3d cp[900];
  } bfptab[64];

void DrawBPatchJ1 ( int n, int m, const point3d *cp )
{
  if ( nbfp < 64 ) {
    bfptab[nbfp].n = n;
    bfptab[nbfp].lkn = 0;
    memcpy ( bfptab[nbfp].cp, cp, (n+1)*(m+1)*sizeof(point3d) );
    nbfp ++;
  }
} /*DrawBPatchJ1*/

void DrawBSPatchJ1 ( int n, int lknu, const double *knu,
                     int m, int lknv, const double *knv, const point3d *cp )
{
  if ( nbfp < 64 ) {
    bfptab[nbfp].n = n;
    bfptab[nbfp].lkn = lknu;
    memcpy ( bfptab[nbfp].knots, knu, (lknu+1)*sizeof(double) );
    memcpy ( bfptab[nbfp].cp, cp, (lknu-n)*(lknv-m)*sizeof(point3d) );
    nbfp ++;
  }
} /*DrawBSPatchJ1*/

void DrawBasAJump ( void )
{
#define DD 32
  void  *sp;
  int   i, j, k;
  double u;
  point3d p1, p2;
  point2d *c, *d;

  sp = pkv_GetScratchMemTop ();
  c = pkv_GetScratchMem ( (DD+1)*sizeof(point2d) );
  d = pkv_GetScratchMem ( (DD+1)*sizeof(point2d) );
  if ( !c || !d )
    exit ( 1 );
  for ( i = 0, j = nbfp-1;  i < nbfp;  j = i++ ) {
    for ( k = 0; k <= DD; k++ ) {
      u = (double)k/(double)DD;
      mbs_BCHornerPd ( bfptab[j].n, bfptab[j].n, 3, &bfptab[j].cp[0].x,
                       0.0, u, &p1.x );
      p1.z = scf*ComputeBezLaplaciand ( bfptab[j].n, bfptab[j].n,
                       &bfptab[j].cp[0], 0.0, u );
      mbs_BCHornerPd ( bfptab[i].n, bfptab[i].n, 3, &bfptab[i].cp[0].x,
                       u, 0.0, &p2.x );
      p2.z = scf*ComputeBezLaplaciand ( bfptab[i].n, bfptab[i].n,
                       &bfptab[i].cp[0], u, 0.0 );
      p1.z = p2.z-p1.z;
      CameraProjectPoint3d ( &CPos, &p1, &p2 );
      SetPoint2d ( &d[k], p2.x, p2.y );
      p1.z = 0.0;
      CameraProjectPoint3d ( &CPos, &p1, &p2 );
      SetPoint2d ( &c[k], p2.x, p2.y );
    }
    ps_Set_Gray ( 0.0 );
    ps_Draw_Polyline2d ( DD+1, c );
    ps_Set_RGB ( 0.0, 0.25, 1.0 );
    ps_Draw_Polyline2d ( DD+1, d );

    for ( k = 0; k <= DD; k++ ) {
      u = (double)k/(double)DD;
      mbs_BCHornerPd ( bfptab[i].n, bfptab[i].n, 3, &bfptab[i].cp[0].x,
                       u, 1.0, &p1.x );
      p1.z = -scf*ComputeBezLaplaciand ( bfptab[i].n, bfptab[i].n,
                       &bfptab[i].cp[0], u, 1.0 );
      CameraProjectPoint3d ( &CPos, &p1, &p2 );
      SetPoint2d ( &d[k], p2.x, p2.y );
      p1.z = 0.0;
      CameraProjectPoint3d ( &CPos, &p1, &p2 );
      SetPoint2d ( &c[k], p2.x, p2.y );
    }
    ps_Set_Gray ( 0.0 );
    ps_Draw_Polyline2d ( DD+1, c );
    ps_Set_RGB ( 0.0, 0.25, 1.0 );
    ps_Draw_Polyline2d ( DD+1, d );

    for ( k = 0; k <= DD; k++ ) {
      u = (double)k/(double)DD;
      mbs_BCHornerPd ( bfptab[i].n, bfptab[i].n, 3, &bfptab[i].cp[0].x,
                       1.0, u, &p1.x );
      p1.z = -scf*ComputeBezLaplaciand ( bfptab[i].n, bfptab[i].n,
                       &bfptab[i].cp[0], 1.0, u );
      CameraProjectPoint3d ( &CPos, &p1, &p2 );
      SetPoint2d ( &d[k], p2.x, p2.y );
      p1.z = 0.0;
      CameraProjectPoint3d ( &CPos, &p1, &p2 );
      SetPoint2d ( &c[k], p2.x, p2.y );
    }
    ps_Set_Gray ( 0.0 );
    ps_Draw_Polyline2d ( DD+1, c );
    ps_Set_RGB ( 0.0, 0.25, 1.0 );
    ps_Draw_Polyline2d ( DD+1, d );
  }
  pkv_SetScratchMemTop ( sp );
#undef DD
} /*DrawBasAJump*/

void DrawBasBJump ( void )
{
#define DD 32
  void  *sp;
  int   i, j, jj, k;
  double u;
  point3d p1, p2;
  point2d *c, *d;

  sp = pkv_GetScratchMemTop ();
  c = pkv_GetScratchMem ( (DD+1)*sizeof(point2d) );
  d = pkv_GetScratchMem ( (DD+1)*sizeof(point2d) );
  if ( !c || !d )
    exit ( 1 );
  for ( i = 0, j = HOLE_K-1;  i < HOLE_K;  j = i++ ) {
    for ( k = 0; k <= DD; k++ ) {
      u = (double)k/(double)DD;
      mbs_BCHornerPd ( bfptab[3*HOLE_K+j].n, bfptab[3*HOLE_K+j].n, 3,
                       &bfptab[3*HOLE_K+j].cp[0].x, 0.0, u, &p1.x );
      p1.z = scf*ComputeBezLaplaciand ( bfptab[3*HOLE_K+j].n, bfptab[3*HOLE_K+j].n,
                       &bfptab[3*HOLE_K+j].cp[0], 0.0, u );
      mbs_BCHornerPd ( bfptab[3*HOLE_K+i].n, bfptab[3*HOLE_K+i].n, 3,
                       &bfptab[3*HOLE_K+i].cp[0].x, u, 0.0, &p2.x );
      p2.z = scf*ComputeBezLaplaciand ( bfptab[3*HOLE_K+i].n, bfptab[3*HOLE_K+i].n,
                       &bfptab[3*HOLE_K+i].cp[0], u, 0.0 );
      p1.z = p2.z-p1.z;
      CameraProjectPoint3d ( &CPos, &p1, &p2 );
      SetPoint2d ( &d[k], p2.x, p2.y );
      p1.z = 0.0;
      CameraProjectPoint3d ( &CPos, &p1, &p2 );
      SetPoint2d ( &c[k], p2.x, p2.y );
    }
    ps_Set_Gray ( 0.0 );
    ps_Draw_Polyline2d ( DD+1, c );
    ps_Set_RGB ( 0.0, 0.25, 1.0 );
    ps_Draw_Polyline2d ( DD+1, d );

    jj = (i+1) % HOLE_K;
    for ( k = 0; k <= DD; k++ ) {
      u = (double)k/(double)DD;
      mbs_BCHornerPd ( bfptab[3*HOLE_K+i].n, bfptab[3*HOLE_K+i].n, 3,
                       &bfptab[3*HOLE_K+i].cp[0].x, u, 1.0, &p1.x );
      p1.z = scf*ComputeBezLaplaciand ( bfptab[3*HOLE_K+i].n, bfptab[3*HOLE_K+i].n,
                       &bfptab[3*HOLE_K+i].cp[0], u, 1.0 );
      mbs_BCHornerPd ( bfptab[3*jj+1].n, bfptab[3*jj+1].n, 3,
                       &bfptab[3*jj+1].cp[0].x, 0.0, u, &p2.x );
      p2.z = scf*ComputeBezLaplaciand ( bfptab[3*jj+1].n, bfptab[3*jj+1].n,
                       &bfptab[3*jj+1].cp[0], 0.0, u );
      p1.z = p2.z-p1.z;
      CameraProjectPoint3d ( &CPos, &p1, &p2 );
      SetPoint2d ( &d[k], p2.x, p2.y );
      p1.z = 0.0;
      CameraProjectPoint3d ( &CPos, &p1, &p2 );
      SetPoint2d ( &c[k], p2.x, p2.y );
    }
    ps_Set_Gray ( 0.0 );
    ps_Draw_Polyline2d ( DD+1, c );
    ps_Set_RGB ( 0.0, 0.25, 1.0 );
    ps_Draw_Polyline2d ( DD+1, d );

    for ( k = 0; k <= DD; k++ ) {
      u = (double)k/(double)DD;
      mbs_BCHornerPd ( bfptab[3*HOLE_K+i].n, bfptab[3*HOLE_K+i].n, 3,
                       &bfptab[3*HOLE_K+i].cp[0].x, 1.0, u, &p1.x );
      p1.z = scf*ComputeBezLaplaciand ( bfptab[3*HOLE_K+i].n, bfptab[3*HOLE_K+i].n,
                       &bfptab[3*HOLE_K+i].cp[0], 1.0, u );
      mbs_BCHornerPd ( bfptab[3*i+2].n, bfptab[3*i+2].n, 3,
                       &bfptab[3*i+2].cp[0].x, 0.0, u, &p2.x );
      p2.z = scf*ComputeBezLaplaciand ( bfptab[3*i+2].n, bfptab[3*i+2].n,
                             &bfptab[3*i+2].cp[0], 0.0, u );
      p1.z = p2.z-p1.z;
      CameraProjectPoint3d ( &CPos, &p1, &p2 );
      SetPoint2d ( &d[k], p2.x, p2.y );
      p1.z = 0.0;
      CameraProjectPoint3d ( &CPos, &p1, &p2 );
      SetPoint2d ( &c[k], p2.x, p2.y );
    }
    ps_Set_Gray ( 0.0 );
    ps_Draw_Polyline2d ( DD+1, c );
    ps_Set_RGB ( 0.0, 0.25, 1.0 );
    ps_Draw_Polyline2d ( DD+1, d );
  }
  pkv_SetScratchMemTop ( sp );
#undef DD
} /*DrawBasBJump*/

void DrawBasCDJump ( void )
{
#define DD 32
  void    *sp;
  int     ku, kv, pitch1, pitch2, pitch3;
  double   *b1, *b2, *c1, *c2, u, *knu, *knv;
  int     i, j, k, l, ll, m, n, lknu, lknv, start1, start2;
  point3d *cp1, *cp2, p1, p2;
  point2d *cc, *dd;

  sp = pkv_GetScratchMemTop ();
  n = m = bfptab[0].n;
  lknu = lknv = bfptab[0].lkn;
  knu = knv = bfptab[0].knots;
  ku = kv = mbs_NumKnotIntervalsd ( n, lknu, knu );
  pitch1 = (lknv-m)*3;
  pitch2 = (m+1)*3*kv;
  pitch3 = (m+1)*3;
  b1 = pkv_GetScratchMemd ( pitch2*ku*(n+1) );
  b2 = pkv_GetScratchMemd ( pitch2*ku*(n+1) );
  c1 = pkv_GetScratchMemd ( (n+1)*pitch3 );
  c2 = pkv_GetScratchMemd ( (n+1)*pitch3 );
  cc = pkv_GetScratchMem ( (DD+1)*sizeof(point2d) );
  dd = pkv_GetScratchMem ( (DD+1)*sizeof(point2d) );
  if ( b1 && b2 && c1 && c2 && cc && dd ) {
    for ( i = 0, j = HOLE_K-1;  i < HOLE_K;  j = i++ ) {
      cp1 = bfptab[j].cp;
      cp2 = bfptab[i].cp;
      mbs_BSPatchToBezd ( 3, n, lknu, knu, m, lknv, knv, pitch1, (double*)cp1,
                          &ku, NULL, NULL, &kv, NULL, NULL, pitch2, b1 );
      mbs_BSPatchToBezd ( 3, n, lknu, knu, m, lknv, knv, pitch1, (double*)cp2,
                          &ku, NULL, NULL, &kv, NULL, NULL, pitch2, b2 );
      for ( k = 0; k < ku; k++ ) {
        start1 = k*pitch3;
        pkv_Selectd ( n+1, pitch3, pitch2, pitch3, &b1[start1], c1 );
        start2 = k*(n+1)*pitch2;
        pkv_Selectd ( n+1, pitch3, pitch2, pitch3, &b2[start2], c2 );
        for ( ll = 0; ll <= DD; ll++ ) {
          u = (double)ll/(double)DD;
          mbs_BCHornerPd ( n, n, 3, c1, 0.0, u, &p1.x );
          p1.z = scf*ComputeBezLaplaciand ( n, n, (point3d*)c1, 0.0, u );
          mbs_BCHornerPd ( n, n, 3, c2, u, 0.0, &p2.x );
          p2.z = scf*ComputeBezLaplaciand ( n, n, (point3d*)c2, u, 0.0 );
          p1.z = p2.z-p1.z;
          CameraProjectPoint3d ( &CPos, &p1, &p2 );
          SetPoint2d ( &dd[ll], p2.x, p2.y );
          p1.z = 0.0;
          CameraProjectPoint3d ( &CPos, &p1, &p2 );
          SetPoint2d ( &cc[ll], p2.x, p2.y );
        }
        ps_Set_Gray ( 0.0 );
        ps_Draw_Polyline2d ( DD+1, cc );
        ps_Set_RGB ( 0.0, 0.25, 1.0 );
        ps_Draw_Polyline2d ( DD+1, dd );
      }

      for ( k = 1; k < ku; k++ ) {
        for ( l = 0; l < kv; l++ ) {
          start1 = (k-1)*(n+1)*pitch2+l*pitch3;
          pkv_Selectd ( n+1, pitch3, pitch2, pitch3, &b2[start1], c1 );
          start2 = k*(n+1)*pitch2+l*pitch3;
          pkv_Selectd ( n+1, pitch3, pitch2, pitch3, &b2[start2], c2 );
          for ( ll = 0; ll <= DD; ll++ ) {
            u = (double)ll/(double)DD;
            mbs_BCHornerPd ( n, n, 3, c1, 1.0, u, &p1.x );
            p1.z = scf*ComputeBezLaplaciand ( n, n, (point3d*)c1, 1.0, u );
            mbs_BCHornerPd ( n, n, 3, c2, 0.0, u, &p2.x );
            p2.z = scf*ComputeBezLaplaciand ( n, n, (point3d*)c2, 0.0, u );
            p1.z = p2.z-p1.z;
            CameraProjectPoint3d ( &CPos, &p1, &p2 );
            SetPoint2d ( &dd[ll], p2.x, p2.y );
            p1.z = 0.0;
            CameraProjectPoint3d ( &CPos, &p1, &p2 );
            SetPoint2d ( &cc[ll], p2.x, p2.y );
          }
          ps_Set_Gray ( 0.0 );
          ps_Draw_Polyline2d ( DD+1, cc );
          ps_Set_RGB ( 0.0, 0.75, 0.25 );
          ps_Draw_Polyline2d ( DD+1, dd );
        }
      }
      for ( l = 0; l < kv; l++ ) {
        start1 = (ku-1)*(n+1)*pitch2+l*pitch3;
        pkv_Selectd ( n+1, pitch3, pitch2, pitch3, &b2[start1], c1 );
        for ( ll = 0; ll <= DD; ll++ ) {
          u = (double)ll/(double)DD;
          mbs_BCHornerPd ( n, n, 3, c1, 1.0, u, &p1.x );
          p1.z = -scf*ComputeBezLaplaciand ( n, n, (point3d*)c1, 1.0, u );
          CameraProjectPoint3d ( &CPos, &p1, &p2 );
          SetPoint2d ( &dd[ll], p2.x, p2.y );
          p1.z = 0.0;
          CameraProjectPoint3d ( &CPos, &p1, &p2 );
          SetPoint2d ( &cc[ll], p2.x, p2.y );
        }
        ps_Set_Gray ( 0.0 );
        ps_Draw_Polyline2d ( DD+1, cc );
        ps_Set_RGB ( 0.0, 0.75, 0.25 );
        ps_Draw_Polyline2d ( DD+1, dd );
      }

      for ( k = 0; k < ku; k++ )
        for ( l = 1; l < kv; l++ ) {
          start1 = k*(n+1)*pitch2+(l-1)*pitch3;
          pkv_Selectd ( n+1, pitch3, pitch2, pitch3, &b2[start1], c1 );
          start2 = k*(n+1)*pitch2+l*pitch3;
          pkv_Selectd ( n+1, pitch3, pitch2, pitch3, &b2[start2], c2 );
          for ( ll = 0; ll <= DD; ll++ ) {
            u = (double)ll/(double)DD;
            mbs_BCHornerPd ( n, n, 3, c1, u, 1.0, &p1.x );
            p1.z = scf*ComputeBezLaplaciand ( n, n, (point3d*)c1, u, 1.0 );
            mbs_BCHornerPd ( n, n, 3, c2, u, 0.0, &p2.x );
            p2.z = scf*ComputeBezLaplaciand ( n, n, (point3d*)c2, u, 0.0 );
            p1.z = p2.z-p1.z;
            CameraProjectPoint3d ( &CPos, &p1, &p2 );
            SetPoint2d ( &dd[ll], p2.x, p2.y );
            p1.z = 0.0;
            CameraProjectPoint3d ( &CPos, &p1, &p2 );
            SetPoint2d ( &cc[ll], p2.x, p2.y );
          }
          ps_Set_Gray ( 0.0 );
          ps_Draw_Polyline2d ( DD+1, cc );
          ps_Set_RGB ( 0.25, 0.0, 0.75 );
          ps_Draw_Polyline2d ( DD+1, dd );
        }
      for ( k = 0; k < ku; k++ ) {
        start1 = k*(n+1)*pitch2+(kv-1)*pitch3;
        pkv_Selectd ( n+1, pitch3, pitch2, pitch3, &b2[start1], c1 );
        for ( ll = 0; ll <= DD; ll++ ) {
          u = (double)ll/(double)DD;
          mbs_BCHornerPd ( n, n, 3, c1, u, 1.0, &p1.x );
          p1.z = -scf*ComputeBezLaplaciand ( n, n, (point3d*)c1, u, 1.0 );
          CameraProjectPoint3d ( &CPos, &p1, &p2 );
          SetPoint2d ( &dd[ll], p2.x, p2.y );
          p1.z = 0.0;
          CameraProjectPoint3d ( &CPos, &p1, &p2 );
          SetPoint2d ( &cc[ll], p2.x, p2.y );
        }
        ps_Set_Gray ( 0.0 );
        ps_Draw_Polyline2d ( DD+1, cc );
        ps_Set_RGB ( 0.25, 0.0, 0.75 );
        ps_Draw_Polyline2d ( DD+1, dd );
      }
    }
  }
  pkv_SetScratchMemTop ( sp );
#undef DD
} /*DrawBasCDJump*/

void DrawBasisFuncBSPatchLJ ( int fn, double x, double y )
{
  char s[50];

  ps_GSave ();
  sprintf ( s, "%6.2f %6.2f translate", x, y );
  ps_Write_Command ( s );
  DrawFN ( 0, 1100.0, fn );
  ps_Set_Gray ( 0.5 );
  maxlap = 0.0;
  g1h_DrawSplBasFunctiond ( domain, fn, DrawBSPatchL1 );
  scf = 0.5/maxlap;
  nbfp = 0;
  if ( fn < nfunc_a ) {
    g1h_DrawBasAFunctiond ( domain, fn, DrawBPatchJ1 );
    DrawBasAJump ();
  }
  else if ( fn < nfunc_a+nfunc_b ) {
    g1h_DrawBasBFunctiond ( domain, fn-nfunc_a, DrawBPatchJ1 );
    DrawBasBJump ();
  }
  else {
    g1h_DrawSplBasFunctiond ( domain, fn, DrawBSPatchJ1 );
    DrawBasCDJump ();
  }
  ps_GRestore ();
} /*DrawBasisFuncBSPatchLJ*/

void MakePicture3x ( const char *fn, int ffn, int nf )
{
  ps_OpenFile ( fn, 600 );
  ps_DenseScreen ();
  ps_Write_Command ( "1 setlinecap" );
  ps_Write_Command ( "/Times-Roman findfont 60 scalefont setfont" );
  SetupCPos ( 1600.0, 1160.0, 0.0, 0.0 );

  DrawBasisFuncBSPatchLJ ( ffn+0, 20.0, 3610.0 );
  if ( nf > 1 )
    DrawBasisFuncBSPatchLJ ( ffn+1, 1660.0, 3610 );
  if ( nf > 2 )
    DrawBasisFuncBSPatchLJ ( ffn+2, 20.0, 2410.0 );
  if ( nf > 3 )
    DrawBasisFuncBSPatchLJ ( ffn+3, 1660.0, 2410.0 );
  if ( nf > 4 )
    DrawBasisFuncBSPatchLJ ( ffn+4, 20.0, 1210.0 );
  if ( nf > 5 )
    DrawBasisFuncBSPatchLJ ( ffn+5, 1660.0, 1210.0 );
  if ( nf > 6 )
    DrawBasisFuncBSPatchLJ ( ffn+6, 20.0, 10.0 );
  if ( nf > 7 )
    DrawBasisFuncBSPatchLJ ( ffn+7, 1660.0, 10.0 );
  ps_CloseFile ();
  printf ( "%s\n", fn );
} /*MakePicture3x*/

void MakePicture3 ( void )
{
  char fn[40];
  int i, j, n;

  n = nfunc_a+nfunc_b+nfunc_c+nfunc_d;
  for ( i = j = 0; j < n; i++, j += 8 ) {
    sprintf ( fn, fn3temp, i );
    MakePicture3x ( fn, j, min (8, n-j) );
  }
} /*MakePicture3*/

/* ///////////////////////////////////////////////////////////////////////// */
void DrawSplMatrix ( int k, int r, int s, int t, double *A, double *B )
{
  void *sp;
  int  n, m, i, j, p;
  byte *buf;

  sp = pkv_GetScratchMemTop ();
  n = k*(r+s)+t;
  m = 6*HOLE_K+1;
  buf = pkv_GetScratchMem ( max(n,m) );
  if ( buf ) {
    ps_GSave ();
    ps_Write_Command ( "100 100 translate 10 dup scale" );
    ps_Init_BitmapP ( n, n, 0, 0 );
    for ( i = 0; i < n; i++ ) {
      memset ( buf, 0xDF, n );
      for ( j = 0; j < n; j++ ) {
        p = pkn_Block2FindElemPos ( k, r, s, t, i, j );
        if ( p >= 0 ) {
          if ( A[p] )
            buf[j] = 0;
        }
        else buf[j] = 0xAF;
      }
      ps_Out_LineP ( buf );
    }
    ps_Init_BitmapP ( m, n, n+10, 0 );
    for ( i = p = 0; i < n; i++ ) {
      memset ( buf, 0xDF, m );
      for ( j = 0; j < m; j++, p++ )
        if ( B[p] ) buf[j] = 0;
      ps_Out_LineP ( buf );
    }
    ps_GRestore ();
  }
  pkv_SetScratchMemTop ( sp );
} /*DrawSplMatrix*/

void MakePicture4 ( void )
{
  ps_OpenFile ( fn3, 600 );
  g1h_DrawSplMatricesd ( domain, DrawSplMatrix );
  ps_CloseFile ();
  printf ( "%s\n", fn3 );
} /*MakePicture4*/

/* ///////////////////////////////////////////////////////////////////////// */
void MakePictures ( void )
{
#ifdef DRAW_BASIS
  MakePicture1 ();
#endif
#ifdef DRAW_LAPLACIAN
  MakePicture2 ();
#endif
#ifdef DRAW_LAPLACIAN_JUMP
  MakePicture3 ();
#endif
  MakePicture4 ();
} /*MakePictures*/

/* ///////////////////////////////////////////////////////////////////////// */
int   mk, mr, ms, mt, ms1, ms2;
double *Amat = NULL, *Bmat = NULL;
double *AAmat = NULL, *BBmat = NULL;
double *Lmat = NULL;

int     cnt;
int     sdegu1[HOLE_K], lknu1[HOLE_K], sdegv1[HOLE_K], lknv1[HOLE_K];
double   *knu1[HOLE_K], *knv1[HOLE_K];
point3d *scp1[HOLE_K];
int     sdegu2[HOLE_K], lknu2[HOLE_K], sdegv2[HOLE_K], lknv2[HOLE_K];
double   *knu2[HOLE_K], *knv2[HOLE_K];
point3d *scp2[HOLE_K];

void AllocPatchStorage ( void )
{
  int i, n;

  n = G1H_FINALDEG+1+NK*max(M1+4,M2);
  for ( i = 0; i < HOLE_K; i++ ) {
    knu1[i] = malloc ( (n+G1H_FINALDEG+1)*sizeof(double) );
    knv1[i] = malloc ( (n+G1H_FINALDEG+1)*sizeof(double) );
    scp1[i] = malloc ( n*n*sizeof(point3d) );
    knu2[i] = malloc ( (n+G1H_FINALDEG+1)*sizeof(double) );
    knv2[i] = malloc ( (n+G1H_FINALDEG+1)*sizeof(double) );
    scp2[i] = malloc ( n*n*sizeof(point3d) );
    if ( !knu1[i] || !knv1[i] || !scp1[i] ||
         !knu2[i] || !knv2[i] || !scp2[i] )
      exit ( 1 );
  }
} /*AllocPatchStorage*/

void SavePatch1 ( int n, int lknu, const double *knu,
                  int m, int lknv, const double *knv, const point3d *cp )
{
  if ( cnt < HOLE_K ) {
    sdegu1[cnt] = n;
    lknu1[cnt] = lknu;
    memcpy ( knu1[cnt], knu, (lknu+1)*sizeof(double) );
    sdegv1[cnt] = m;
    lknv1[cnt] = lknv;
    memcpy ( knv1[cnt], knv, (lknv+1)*sizeof(double) );
    memcpy ( scp1[cnt], cp, (lknu-n)*(lknv-m)*sizeof(point3d) );
    cnt++;
  }
  else
    exit ( 1 );
} /*SavePatch1*/

void SavePatch2 ( int n, int lknu, const double *knu,
                  int m, int lknv, const double *knv, const point3d *cp )
{
  if ( cnt < HOLE_K ) {
    sdegu2[cnt] = n;
    lknu2[cnt] = lknu;
    memcpy ( knu2[cnt], knu, (lknu+1)*sizeof(double) );
    sdegv2[cnt] = m;
    lknv2[cnt] = lknv;
    memcpy ( knv2[cnt], knv, (lknv+1)*sizeof(double) );
    memcpy ( scp2[cnt], cp, (lknu-n)*(lknv-m)*sizeof(point3d) );
    cnt++;
  }
  else
    exit ( 1 );
} /*SavePatch2*/

void SaveMatrices ( int k, int r, int s, int t, double *A, double *B )
{
  mk = k;  mr = r;  ms = s;  mt = t;
  ms1 = pkn_Block2ArraySize ( k, r, s, t );
  ms2 = (k*(r+s)+t)*(6*k+1);
  printf ( "s1=%d, s2=%d\n", ms1, ms2 );
  Amat = malloc ( ms1*sizeof(double) );
  Bmat = malloc ( ms2*sizeof(double) );
  if ( !Amat || !Bmat )
    exit ( 1 );
  memcpy ( Amat, A, ms1*sizeof(double) );
  memcpy ( Bmat, B, ms2*sizeof(double) );
} /*SaveMatrices*/

boolean IsZnonZero ( int n, int lknu, const double *knu,
                     int m, int lknv, const double *knv, const point3d *cp )
{
  int ncp, i;

  ncp = (lknu-n)*(lknv-m);
  for ( i = 0; i < ncp; i++ )
    if ( cp[i].z )
      return true;
  return false;
} /*IsZnonZero*/

static boolean dumpspl = false;

double CompIntegral ( void )
{
  void     *sp;
  int      i, j, k, sz, d1, d2;
  double   s, ss, sss;
  double    u, v;
  point3d  *ap1, *ap2, p1[16], p2[16];
  double    lap1, lap2;
  double    jac1, jac2, gx, gy, gxx, gxy, gyy;
  FILE     *f;

  if ( dumpspl )
    f = fopen ( "intspl1.txt", "w+" );
  sp = pkv_GetScratchMemTop ();
  sz = 0;
  for ( i = 0; i < HOLE_K; i++ ) {
    sz = max ( sz, lknu1[i]-sdegu1[i] );
    sz = max ( sz, lknv1[i]-sdegv1[i] );
    sz = max ( sz, lknu2[i]-sdegu2[i] );
    sz = max ( sz, lknv2[i]-sdegv2[i] );
  }
  ap1 = pkv_GetScratchMem ( 3*sz*sizeof(point3d) );
  ap2 = pkv_GetScratchMem ( 3*sz*sizeof(point3d) );
  if ( !ap1 || !ap2 )
    exit ( 1 );
  sss = 0.0;
  for ( i = 0; i < HOLE_K; i++ ) {
    if ( IsZnonZero ( sdegu1[i], lknu1[i], knu1[i],
                      sdegv1[i], lknv1[i], knv1[i], scp1[i] ) &&
         IsZnonZero ( sdegu2[i], lknu2[i], knu2[i],
                      sdegv2[i], lknv2[i], knv2[i], scp2[i] ) ) {
      d1 = lknv1[i]-sdegv1[i];
      d2 = lknv2[i]-sdegv2[i];
      ss = 0.0;
      for ( j = 0; j < QUAD_FACTOR*(NK+1); j++ ) {
        u = (double)(j+j+1)/(double)(2*QUAD_FACTOR*(NK+1));
        mbs_multideBoorDer2d ( sdegu1[i], lknu1[i], knu1[i], 1, 3*d1, 0,
            (double*)scp1[i], u, &ap1[0].x, &ap1[d1].x, &ap1[2*d1].x );
        mbs_multideBoorDer2d ( sdegu2[i], lknu2[i], knu2[i], 1, 3*d2, 0,
            (double*)scp2[i], u, &ap2[0].x, &ap2[d2].x, &ap2[2*d2].x );
        s = 0.0;
        for ( k = 0; k < QUAD_FACTOR*(NK+1); k++ ) {
          v = (double)(k+k+1)/(double)(2*QUAD_FACTOR*(NK+1));
          mbs_multideBoorDer2d ( sdegv1[i], lknv1[i], knv1[i], 3, 3, 3*d1,
              (double*)ap1, v, &p1[0].x, &p1[4].x, &p1[8].x );
          mbs_multideBoorDer2d ( sdegv2[i], lknv2[i], knv2[i], 3, 3, 3*d2,
              (double*)ap2, v, &p2[0].x, &p2[4].x, &p2[8].x );
          pkn_Comp2iDerivatives2d ( p1[1].x, p1[1].y, p1[4].x, p1[4].y,
              p1[2].x, p1[2].y, p1[5].x, p1[5].y, p1[8].x, p1[8].y,
              1, &p1[1].z, &p1[4].z, &p1[2].z, &p1[5].z, &p1[8].z,
              &gx, &gy, &gxx, &gxy, &gyy );
          lap1 = gxx+gyy;
          jac1 = fabs(p1[4].x*p1[1].y-p1[4].y*p1[1].x);
          pkn_Comp2iDerivatives2d ( p2[1].x, p2[1].y, p2[4].x, p2[4].y,
              p2[2].x, p2[2].y, p2[5].x, p2[5].y, p2[8].x, p2[8].y,
              1, &p2[1].z, &p2[4].z, &p2[2].z, &p2[5].z, &p2[8].z,
              &gx, &gy, &gxx, &gxy, &gyy );
          lap2 = gxx+gyy;
          jac2 = fabs(p2[4].x*p2[1].y-p2[4].y*p2[1].x);
          s += lap1*lap2*jac1;
          if ( dumpspl )
            fprintf ( f, "%d, %f, %f, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e\n",
                      j*QUAD_FACTOR*(NK+1)+k, u, v,
                      p2[1].x, p2[1].y, p2[4].x, p2[4].y,
                      p2[2].x, p2[2].y, p2[5].x, p2[5].y, p2[8].x, p2[8].y,
                      p2[1].z, p2[4].z, p2[2].z, p2[5].z, p2[8].z,
                      lap1, lap2, jac1 );
        }
        ss += s;
      }
      sss += ss;
    }
  }
  pkv_SetScratchMemTop ( sp );
  sss *= 1.0/(double)(QUAD_FACTOR*QUAD_FACTOR*(NK+1)*(NK+1));
  if ( dumpspl ) {
    fprintf ( f, "\n" );
    fclose ( f );
  }
  return sss;
} /*CompIntegral*/

void ComputeMatrices ( void )
{
/* for verification - this computation is slow, but it has */
/* a better chance of being programmed without bugs */
  int i, j, fni, fnj, k;

  memset ( AAmat, 0, ms1*sizeof(double) );
  memset ( BBmat, 0, ms2*sizeof(double) );
  for ( i = 0; i < nfunc_a+nfunc_c+nfunc_d; i++ ) {
    if ( i < nfunc_c+nfunc_d ) fni = nfunc_a+nfunc_b+i;
    else fni = i-nfunc_c-nfunc_d;
    printf ( "%4d %4d\b\b\b\b\b\b\b\b\b", i, fni );
    cnt = 0;
    g1h_DrawSplBasFunctiond ( domain, fni, SavePatch1 );
    for ( j = 0; j <= i; j++ ) {

dumpspl = (boolean)(i == 0 && j == 0);

      if ( j < nfunc_c+nfunc_d ) fnj = nfunc_a+nfunc_b+j;
      else fnj = j-nfunc_c-nfunc_d;
      k = pkn_Block2FindElemPos ( mk, mr, ms, mt, i, j );
      if ( k >= 0 ) {  /* this element belongs to a nonzero block of A */
        cnt = 0;
        g1h_DrawSplBasFunctiond ( domain, fnj, SavePatch2 );
        AAmat[k] = CompIntegral ();
      }
    }
dumpspl = false;
    for ( j = 0; j < nfunc_b; j++ ) {
      cnt = 0;
      g1h_DrawSplBasFunctiond ( domain, nfunc_a+j, SavePatch2 );
      BBmat[nfunc_b*i+j] = CompIntegral ();
    }
  }
  printf ( "\n" );
} /*ComputeMatrices*/

void CompareMatrices ( void )
{
  FILE *cmp;
  int i;
  double e;

  cmp = fopen ( "matr.txt", "w+" );
  for ( i = 0; i < ms1; i++ ) {
    e = (AAmat[i]-Amat[i])/Amat[i];
    fprintf ( cmp, "%6d: %8g, %8g, %8g, %8g ", i,
             Amat[i], AAmat[i], AAmat[i]-Amat[i], e );
    if ( fabs(e) > 1.0e-3 )
      fprintf ( cmp, "*" );
    fprintf ( cmp, "\n" );
  }
  fprintf ( cmp, "\n" );
  for ( i = 0; i < ms2; i++ ) {
    e = (BBmat[i]-Bmat[i])/Bmat[i];
    fprintf ( cmp, "%6d: %8g, %8g, %8g, %8g ", i,
             Bmat[i], BBmat[i], BBmat[i]-Bmat[i], e );
    if ( fabs(e) > 1.0e-3 )
      fprintf ( cmp, "*" );
    fprintf ( cmp, "\n" );
  }
  fclose ( cmp );
} /*CompareMatrices*/

void WriteBlockAddr ( void )
{
  int i, j;

  for ( i = 0; i <= 2*HOLE_K; i++ ) {
    for ( j = 0; j <= i; j++ )
      printf ( "%5d ",
               pkn_Block2FindBlockPos ( mk, mr, ms, mt, i, j ) );
    printf ( "\n" );
  }
  printf ( "\n" );
} /*WriteBlockAddr*/

void TestMatrix1 ( void )
{
  AllocPatchStorage ();
  g1h_DrawSplMatricesd ( domain, SaveMatrices );
  WriteBlockAddr ();
  AAmat = malloc ( ms1*sizeof(double) );
  BBmat = malloc ( ms2*sizeof(double) );
  if ( !AAmat || !BBmat )
    exit ( 1 );
  ComputeMatrices ();
  CompareMatrices ();
  free ( Amat );
  free ( Bmat );
  free ( AAmat );
  free ( BBmat );
} /*TestMatrix1*/

/* ///////////////////////////////////////////////////////////////////////// */
void FindCondNumber ( void )
{
  void  *sp;
  int   n;
  double *x, *y, ca, pr, yl, lambda1, lambda2;
  int   i;

  sp = pkv_GetScratchMemTop ();
  n = mk*(mr+ms)+mt;
  x = pkv_GetScratchMemd ( n );
  y = pkv_GetScratchMemd ( n );
  if ( !x || !y )
    exit ( 1 );

  for ( x[0] = 1.0, i = 1; i < n; i++ )
    x[i] = x[i-1]*-1.0001;
  pr = pkn_SecondNormd ( n, x );
  pkn_MultMatrixNumd ( 1, n, 0, x, 1.0/pr, 0, x );
  for ( i = 0; i < 1000; i++ ) {
    pkn_Block2SymMatrixMultd ( mk, mr, ms, mt, Amat, 1, 1, x, 1, y );
    pr = pkn_ScalarProductd ( n, x, y );
    yl = pkn_SecondNormd ( n, y );
    ca = pr/yl;
    pkn_MultMatrixNumd ( 1, n, 0, y, 1.0/yl, 0, x );
    if ( ca > 0.9999999 )
      break;
  }
  printf ( "i = %d, ca = %f, yl = %f\n", i, ca, yl );
  lambda1 = yl;

  for ( x[0] = 1.0, i = 1; i < n; i++ )
    x[i] = x[i-1]*-1.0001;
  pr = pkn_SecondNormd ( n, x );
  pkn_MultMatrixNumd ( 1, n, 0, x, 1.0/pr, 0, x );
  for ( i = 0; i < 1000; i++ ) {
    memcpy ( y, x, n*sizeof(double) );
    pkn_Block2LowerTrMSolved ( mk, mr, ms, mt, Lmat, 1, 1, y );
    pkn_Block2UpperTrMSolved ( mk, mr, ms, mt, Lmat, 1, 1, y );
    pr = pkn_ScalarProductd ( n, x, y );
    yl = pkn_SecondNormd ( n, y );
    ca = pr/yl;
    pkn_MultMatrixNumd ( 1, n, 0, y, 1.0/yl, 0, x );
    if ( ca > 0.9999999 )
      break;
  }
  printf ( "i = %d, ca = %f, yl = %f\n", i, ca, yl );
  lambda2 = 1.0/yl;
  printf ( "cond A = %e\n", lambda1/lambda2 );

  pkv_SetScratchMemTop ( sp );
} /*FindCondNumber*/

void TestMatrix2 ( void )
{
  AllocPatchStorage ();
  g1h_DrawSplMatricesd ( domain, SaveMatrices );
  WriteBlockAddr ();
  Lmat = malloc ( ms1*sizeof(double) );
  if ( !Lmat )
    exit ( 1 );
  memcpy ( Lmat, Amat, ms1*sizeof(double) );
  if ( pkn_Block2CholeskyDecompMd ( mk, mr, ms, mt, Lmat ) ) {
    FindCondNumber ();
  }
  else
    printf ( "Cannot decompose.\n" );
  free ( Lmat );
  free ( Amat );
  free ( Bmat );
} /*TestMatrix2*/

/* ///////////////////////////////////////////////////////////////////////// */
int main ()
{
  struct tms start, stop;
  double time;
  double param[2] = {0.0,0.0};

  setvbuf ( stdout, NULL, _IONBF, 0 );
  pkv_InitScratchMem ( 16777216 ); /* 16MB */
  InitKnots ( HOLE_K );
  ModifyKnots ( HOLE_K );
  InitDomain ( HOLE_K, param );
  times ( &start );
  if ( (domain = gh_CreateDomaind ( HOLE_K, knots, Dompt )) ) {
    if ( g1h_ComputeSplBasisd ( domain, NK, M1, M2 ) ) {
      if ( g1h_ComputeSplFormMatrixd ( domain ) ) {
        times ( &stop );
        time = (double)(stop.tms_utime-start.tms_utime)/
               (double)(sysconf(_SC_CLK_TCK));
        printf ( "time = %8.3f\n", time );
        g1h_DrawSplBasFuncNumd ( domain, &nfunc_a, &nfunc_b, &nfunc_c, &nfunc_d );
        printf ( "k = %d, fa = %d, fb = %d, fc = %d, fd = %d\n",
                 HOLE_K, nfunc_a, nfunc_b, nfunc_c, nfunc_d );
        MakePictures ();
#ifdef TEST_MATRIX1
        times ( &start );
        TestMatrix1 ();
        times ( &stop );
        time = (double)(stop.tms_utime-start.tms_utime)/
               (double)(sysconf(_SC_CLK_TCK));
        printf ( "time = %8.3f\n", time );
#endif
#ifdef TEST_MATRIX2
        times ( &start );
        TestMatrix2 ();
        times ( &stop );
        time = (double)(stop.tms_utime-start.tms_utime)/
               (double)(sysconf(_SC_CLK_TCK));
        printf ( "time = %8.3f\n", time );
#endif
      }
    }
    gh_DestroyDomaind ( domain );
  }

  printf ( "Scratch memory used: %d bytes\a\n", pkv_MaxScratchTaken() );
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

