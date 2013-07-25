
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

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

#include "eg1holef.h"
#include "datagenf.h"
#include "testgraphf.h"

#define HOLE_K 6

char fn1[] = "g1holeextf.ps";
char fn2[] = "g1extmatrixf.ps";

CameraRecf   CPos;

GHoleDomainf *domain;

int nfinalp;
point3f FinalCPoints[HOLE_K][121];

float iparams[3] = {0.5,0.0,0.5};

/* ///////////////////////////////////////////////////////////////////////// */
void SetupCamera ( float w, float h, float x, float y )
{
  CameraInitFramef ( &CPos, false, true, w, h, x, y, 1.0, 4 );
  CameraInitPosf ( &CPos );
  SetPoint3f ( &CPos.position, 20.402914, -7.444557, -26.986362 );
  CPos.psi = -0.173305;  CPos.theta = 0.677664;  CPos.phi = -1.920662;
  CPos.vd.persp.f = 4.0;
  CameraSetMappingf ( &CPos );
} /*SetupCamera*/

void DrawBezCurve3f ( int n, const point3f *cp, float t0, float t1 )
{
#define DD 32
  void    *sp;
  int     i;
  point3f p, q;
  point2f *cc;

  sp = pkv_GetScratchMemTop ();
  cc = (point2f*)pkv_GetScratchMem ( (DD+1)*sizeof(point2f) );
  if ( !cc )
    exit ( 1 );

  for ( i = 0; i <= DD; i++ ) {
    mbs_BCHornerC3f ( n, cp, (float)i/(float)DD, &p );
    CameraProjectPoint3f ( &CPos, &p, &q );
    SetPoint2f ( &cc[i], q.x, q.y );
  }
  ps_Draw_Polyline2f ( DD+1, cc );

  pkv_SetScratchMemTop ( sp );
#undef DD
} /*DrawBezCurve3f*/

void SetLW ( int i, int i0, int i1 )
{
    if ( i == i0 || i == i1 )
      ps_Set_Line_Width ( 4.0 );
    else
      ps_Set_Line_Width ( 2.0 );
} /*SetLW*/

void DrawBezPatch3f ( int n, int m, const point3f *cp )
{
#define D 8
  void    *sp;
  int     i;
  point3f *cc;

  sp = pkv_GetScratchMemTop ();
  cc = (point3f*)pkv_GetScratchMem ( (max(n,m)+1)*sizeof(point3f) );
  if ( !cc )
    exit ( 1 );

  for ( i = 0; i <= D; i++ ) {
    mbs_multiBCHornerf ( n, 1, 3*(m+1), 0, (float*)cp, (float)i/(float)D,
                         (float*)cc );
    SetLW ( i, 0, D );
    DrawBezCurve3f ( m, cc, 0.0, 1.0 );
  }
  for ( i = 0; i <= D; i++ ) {
    mbs_multiBCHornerf ( m, n+1, 3, 3*(m+1), (float*)cp, (float)i/(float)D,
                         (float*)cc );
    SetLW ( i, 0, D );
    DrawBezCurve3f ( n, cc, 0.0, 1.0 );
  }

  pkv_SetScratchMemTop ( sp );
#undef D
} /*DrawBezPatch3f*/

void DrawPatchNet3f ( int n, int m, const point3f *cp )
{
  void    *sp;
  int     i, j;
  point3f q;
  point2f *r;

  sp = pkv_GetScratchMemTop ();
  r = pkv_GetScratchMem ( (n+1)*(m+1)*sizeof(point2f) );
  if ( !r )
    exit ( 1 );
  for ( i = 0; i < (n+1)*(m+1); i++ ) {
    CameraProjectPoint3f ( &CPos, &cp[i], &q );
    SetPoint2f ( &r[i], q.x, q.y );
  }
  ps_Set_Line_Width ( 1.0 );
  for ( i = 0; i <= n; i++ )
    ps_Draw_Polyline2f ( m+1, &r[i*(m+1)] );
  for ( j = 0; j <= m; j++ )
    for ( i = 0; i < n; i++ )
      ps_Draw_Line ( r[i*(m+1)+j].x, r[i*(m+1)+j].y,
                     r[(i+1)*(m+1)+j].x, r[(i+1)*(m+1)+j].y );
  for ( i = 0; i < (n+1)*(m+1); i++ )
    ps_Mark_Circle ( r[i].x, r[i].y );

  pkv_SetScratchMemTop ( sp );
} /*DrawPatchNet3f*/

/* ///////////////////////////////////////////////////////////////////////// */
void GetHoleSurrndPatch ( int i, int j, point3f *dombcp )
{
  int     k;
  int     *ind;
  point3f *q;
  void    *sp;
  float   *ukn, *vkn;

  sp  = pkv_GetScratchMemTop ();
  ind = pkv_GetScratchMem ( 16*sizeof(int) );
  q   = pkv_GetScratchMem ( 16*sizeof(point3f) );
  if ( !ind || !q )
    exit ( 1 );

  gh_GetBspInd ( HOLE_K, i, j, ind );
  for ( k = 0; k < 16; k++ )
    q[k] = Bpt[ind[k]];

  ukn = &knots[11*((i+HOLE_K-1) % HOLE_K)+3];
  vkn = &knots[11*i+j];
  mbs_BSPatchToBezf ( 3, 3, 7, ukn, 3, 7, vkn, 12, (float*)q,
                      NULL, NULL, NULL, NULL, NULL, NULL, 12, (float*)dombcp );

  pkv_SetScratchMemTop ( sp );
} /*GetHoleSurrndPatch*/

void DrawBezPatches ()
{
  int     i, j;
  point3f cp[16];

  ps_Set_Gray ( 0.72 );
  for ( i = 0; i < HOLE_K; i++ )
    for ( j = 0; j < 3; j++ ) {
      GetHoleSurrndPatch ( i, j, cp );
      DrawBezPatch3f ( 3, 3, cp );
    }
} /*DrawBezPatches*/

void DrawFinalPatches ()
{
  int i;

  ps_Set_Gray ( 0.25 );
  for ( i = 0; i < nfinalp; i++ )
    DrawBezPatch3f ( G1H_FINALDEG, G1H_FINALDEG, FinalCPoints[i] );
} /*DrawFinalPatches*/

void MakePicture1 ()
{
  ps_OpenFile ( fn1, 600 );
  ps_Write_Command ( "1 setlinecap" );
  SetupCamera ( 1800, 1460, 0.0, 0.0 );

  DrawBezPatches ();
  DrawFinalPatches ();

  ps_CloseFile ();
  printf ( "%s\n", fn1 );
} /*MakePicture1*/

/* ///////////////////////////////////////////////////////////////////////// */
float amax;

float CRadius ( float a, float scf )
{
  if ( a < 0.0 )
    ps_Set_Gray ( 0.68 );
  else
    ps_Set_Gray ( 0.0 );
  return 0.66*scf*pow( fabs(a)/amax, 0.2 );
} /*CRadius*/

void DrawSymBlock ( int size, float scf, float posx, float posy, float *coeff )
{
  int   i, j, k;
  float r;

  for ( i = k = 0;  i < size; i++ ) {
    for ( j = 0;  j < i;  j++, k++ ) {
      r = CRadius ( coeff[k], scf );
      if ( r > 0.05 ) {
        ps_Fill_Circle ( posx + i*scf + 0.5*scf, posy - j*scf - 0.5*scf, r );
        ps_Fill_Circle ( posx + j*scf + 0.5*scf, posy - i*scf - 0.5*scf, r );
      }
    }
    r = CRadius ( coeff[k], scf );
    if ( r > 0.05 )
      ps_Fill_Circle ( posx + i*scf + 0.5*scf, posy - i*scf - 0.5*scf, r );
    k++;
  }
} /*DrawSymBlock*/

void DrawFullBlock ( int nrows, int ncols, float scf,
                     float posx, float posy, float *coeff )
{
  int   i, j, k;
  float r;

  for ( i = k = 0;  i < nrows;  i++ )
    for ( j = 0;  j < ncols;  j++, k++ ) {
      r = CRadius ( coeff[k], scf );
      if ( r > 0.05 )
        ps_Fill_Circle ( posx + j*scf + 0.5*scf, posy - i*scf - 0.5*scf, r );
    }
} /*DrawFullBlock*/

void DrawFullBlockT ( int nrows, int ncols, float scf,
                      float posx, float posy, float *coeff )
{
  int   i, j, k;
  float r;

  for ( i = k = 0;  i < nrows;  i++ )
    for ( j = 0;  j < ncols;  j++, k++ ) {
      r = CRadius ( coeff[k], scf );
      if ( r > 0.05 )
        ps_Fill_Circle ( posx + i*scf + 0.5*scf, posy - j*scf - 0.5*scf, r );
    }
} /*DrawFullBlockT*/

void DrawMatrix ( int k, int r, int s, float *Aii, float *Bi )
{
  int   i, n;
  float scf;
  float *Aki, *Akk;

  Aki = &Aii[pkn_Block1FindBlockPos ( k, r, s, k, 0 )];
  Akk = &Aii[pkn_Block1FindBlockPos ( k, r, s, k, k )];
  n = s+k*r;
  printf ( "n = %d\n", n );
  scf = 25.0;

  amax = fabs(Aii[0]);
  for ( i = 1; i < k*r*(r+1)/2; i++ )
    amax = max ( amax, fabs(Aii[i]) );
  for ( i = 0; i < s*(s+1)/2; i++ )
    amax = max ( amax, fabs(Akk[i]) );
  for ( i = 0; i < k*r*s; i++ )
    amax = max ( amax, fabs(Aki[i]) );

  for ( i = 0; i < k; i++ ) {
    DrawSymBlock ( r, scf, i*r*scf, (n-i*r)*scf, &Aii[r*(r+1)/2] );
    DrawFullBlock ( s, r, scf, i*r*scf, s*scf, &Aki[i*r*s] );
    DrawFullBlockT ( s, r, scf, k*r*scf, (n-i*r)*scf, &Aki[i*r*s] );
  }
  DrawSymBlock ( s, scf, k*r*scf, s*scf, Akk );
} /*DrawMatrix*/

void ExtCondNumber ( int k, int r, int s, float *Aii, float *Bi )
{
#define ITER 50
  void  *sp;
  int   n, i;
  float *Lii, *x, *y;
  float ca, pr, yl, lmin, lmax;

  sp = pkv_GetScratchMemTop ();
  n = k*r+s;
  Lii = pkv_GetScratchMemf ( pkn_Block1ArraySize ( k, r, s ) );
  x = pkv_GetScratchMemf ( n );
  y = pkv_GetScratchMemf ( n );
  if ( Lii && x && y ) {
    memcpy ( Lii, Aii, pkn_Block1ArraySize ( k, r, s )*sizeof(float) );
    if ( !pkn_Block1CholeskyDecompMf ( k, r, s, Lii ) )
      exit ( 1 );

    memset ( x, 0, n*sizeof(float) );
    x[0] = 1.0;
    for ( i = 0; i < ITER; i++ ) {
      pkn_Block1SymMatrixMultf ( k, r, s, Aii, 1, 1, x, 1, y );
      pr = pkn_ScalarProductf ( n, x, y );
      yl = pkn_SecondNormf ( n, y );
      ca = pr/yl;pkn_MultMatrixNumf ( 1, n, 0, y, 1.0/yl, 0, x );
      if  (ca > 0.999999 )
        break;
    }
    lmax = yl;

    memset ( x, 0, n*sizeof(float) );
    x[0] = 1.0;
    for ( i = 0; i < ITER; i++ ) {
      memcpy ( y, x, n*sizeof(float) );
      pkn_Block1LowerTrMSolvef ( k, r, s, Lii, 1, 1, y );
      pkn_Block1UpperTrMSolvef ( k, r, s, Lii, 1, 1, y );
      pr = pkn_ScalarProductf ( n, x, y );
      yl = pkn_SecondNormf ( n, y );
      ca = pr/yl;pkn_MultMatrixNumf ( 1, n, 0, y, 1.0/yl, 0, x );
      if  (ca > 0.999999 )
        break;
    }
    lmin = 1.0/yl;
    printf ( "lambda_max = %f,  lambda_min = %f,  cond_2(A) = %f\n",
             lmax, lmin, lmax*yl );
  }

  pkv_SetScratchMemTop ( sp );
#undef ITER
} /*ExtCondNumber*/

void MakePicture2 ()
{
  ps_OpenFile ( fn2, 600 );
  g1h_DrawExtMatricesf ( domain, DrawMatrix );
  g1h_DrawExtMatricesf ( domain, ExtCondNumber );
  ps_CloseFile ();
  printf ( "%s\n", fn2 );
} /*MakePicture2*/

/* ///////////////////////////////////////////////////////////////////////// */
void MakePictures ()
{
  MakePicture1 ();
  MakePicture2 ();
} /*MakePictures*/

/* ///////////////////////////////////////////////////////////////////////// */
void OutPatch ( int n, int m, const float *cp, void *usrptr )
{
  if ( n != G1H_FINALDEG || m != G1H_FINALDEG || nfinalp >= HOLE_K )
    exit ( 1 );
  memcpy ( FinalCPoints[nfinalp], cp, 121*sizeof(point3f) );
  nfinalp ++;
} /*OutPatch*/

int main ()
{
  struct tms start, stop;
  float time;
  float param[2] = {0.0,0.0};
  float funcval, *acoeff;

  pkv_InitScratchMem ( 4194304 );
  acoeff = pkv_GetScratchMemf ( 474 );
  if ( !acoeff )
    goto failure;
  InitKnots ( HOLE_K );
  InitDomain ( HOLE_K, param );
  times ( &start );
  if ( (domain = gh_CreateDomainf ( HOLE_K, knots, Dompt )) ) {
    if ( g1h_ComputeBasisf ( domain ) ) {
      if ( g1h_ComputeExtFormMatrixf ( domain ) ) {
        nfinalp = 0;
        InitHole ( HOLE_K, iparams );
        if ( !g1h_ExtFillHolef ( domain, 3, (float*)Bpt, acoeff, NULL, OutPatch ) )
          goto failure;
        times ( &stop );
        funcval = g1h_ExtFunctionalValuef ( domain, 3, (float*)Bpt, acoeff );
        printf ( "Functional value = %f\n", funcval );
        MakePictures ();
        time = (float)(stop.tms_utime-start.tms_utime)/
               (float)(sysconf(_SC_CLK_TCK));
        printf ( "time = %8.3f\n", time );
      }
    }
    gh_DestroyDomainf ( domain );
  }

  printf ( "Scratch memory used: %d bytes\n", pkv_MaxScratchTaken () );
  pkv_DestroyScratchMem ();
  exit ( 0 );

failure:
  exit ( 1 );
} /*main*/

