
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

#define NEWDATA

#include "eg1holef.h"
#include "datagenf.h"
#ifdef NEWDATA
#include "datagenf.new.h"
#endif
#include "testgraphf.h"
#include "drawitf.h"

#define HOLE_K 5

char fn1[]  = "g1hq2nonlinf.ps";

CameraRecf   CPos;

GHoleDomainf *domain;

int nfinalp;
point3f FinalCPoints[2*HOLE_K][(G1H_FINALDEG+1)*(G1H_FINALDEG+1)];

float dparam[2] = {0.0,0.0};
float iparams[3] = {0.0,0.0,0.0} /*{0.5,0.0,0.5}*/;

#ifdef NEWDATA
float newparams[DATAGEN_SURF_PARAMS] = {0.0,0.0,0.5,0.0,0.0};
#endif

float q2const = 10.0;

vector3f nlnvect;

/* ///////////////////////////////////////////////////////////////////////// */
void SetupCamera ( float w, float h, float x, float y )
{
  CameraInitFramef ( &CPos, false, true, w, h, x, y, 1.0, 4 );
  CameraInitPosf ( &CPos );
  SetPoint3f ( &CPos.position, 20.402914, -7.444557, -26.986362 );
  CPos.psi = -0.173305;  CPos.theta = 0.677664;  CPos.phi = -1.920662;
  CPos.vd.persp.f = 5.0;
  CameraSetMappingf ( &CPos );
#ifdef NEWDATA
  CameraRotXCf ( &CPos, PI );
#endif
} /*SetupCamera*/

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
      DrawBezPatch3f ( 3, 3, cp, 8, 8 );
    }
} /*DrawBezPatches*/

void DrawFinalPatches ( int first )
{
  int i;

  ps_Set_Gray ( 0.25 );
  for ( i = first; i < first+HOLE_K; i++ )
    DrawBezPatch3f ( G1H_FINALDEG, G1H_FINALDEG, FinalCPoints[i], 8, 8 );
} /*DrawFinalPatches*/

void DrawNLNV ()
{
  point3f p0, p1;

  SetPoint3f ( &p1, 0.0, 0.0, 0.0 );
  CameraProjectPoint3f ( &CPos, &p1, &p0 );
  CameraProjectPoint3f ( &CPos, &nlnvect, &p1 );
  ps_Set_Gray ( 0.0 );
  psl_SetLine ( p0.x, p0.y, p1.x, p1.y, 0.0, 1.0 );
  psl_Draw ( 0.0, 0.9, 6.0 );
  psl_Arrow ( 1.0, true );
} /*DrawNLNV*/

void MakePicture1 ()
{
  ps_OpenFile ( fn1, 600 );
  ps_Write_Command ( "1 setlinecap" );

  SetupCamera ( 2400, 1800, 0.0, 0.0 );
  DrawBezPatches ();
  DrawFinalPatches ( 0 );

  SetupCamera ( 2400, 1800, 2400.0, 0.0 );
  DrawBezPatches ();
  DrawFinalPatches ( HOLE_K );

  ps_CloseFile ();
  printf ( "%s\n", fn1 );
} /*MakePicture1*/

void MakePictures ()
{
  MakePicture1 ();
} /*MakePictures*/

/* ///////////////////////////////////////////////////////////////////////// */
float amax;

float CRadius ( float a, float scf )
{
  if ( a < 0.0 )
    ps_Set_RGB ( 1.0, 0.0, 0.0 );
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

void DrawMatrix ( int k, int r, int s, float *Aii, float *Aki, float *Akk )
{
  int   i, n;
  float scf;

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

void ExtCondNumber ( int k, int r, int s, float *Aii )
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
      ca = pr/yl;
      pkn_MultMatrixNumf ( 1, n, 0, y, 1.0/yl, 0, x );
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

void DrawExtMatrix ( int k, int r, int s, float *Aii, float *Aki, float *Akk )
{
  ps_OpenFile ( "nlextmatrix.ps", 600 );
  DrawMatrix ( k, r, s, Aii, Aki, Akk );
  ps_CloseFile ();
  exit ( 0 );
} /*DrawExtMatrix*/

/* ///////////////////////////////////////////////////////////////////////// */
void OutPatch ( int n, int m, const point3f *cp, void *usrptr )
{
  if ( n != G1H_FINALDEG || m != G1H_FINALDEG || nfinalp >= 2*HOLE_K )
    exit ( 1 );
  memcpy ( FinalCPoints[nfinalp], cp,
          (G1H_FINALDEG+1)*(G1H_FINALDEG+1)*sizeof(point3f) );
  nfinalp ++;
} /*OutPatch*/

/* ///////////////////////////////////////////////////////////////////////// */
void OutputNLNV ( vector3f *nv )
{
  nlnvect = *nv;
} /*OutputNLNV*/

void SetupCameraG ( float w, float h, float x, float y )
{
  vector3f v;

  CameraInitFramef ( &CPos, false, true, w, h, x, y, 1.0, 4 );
  CameraInitPosf ( &CPos );
  SetVector3f ( &v, 0.0, 0.0, -46.0 );
  CameraMoveGf ( &CPos, &v );
  CameraRotXGf ( &CPos, -0.65*PI );
  CameraRotZGf ( &CPos, 0.05*PI );
  CameraZoomf ( &CPos, 10.0 );
} /*SetupCameraG*/

void DrawFGraph ( int k, int n, point2f *xy, float *z )
{
/*#define SF 0.03125*/
/*#define SF 0.0625*/
#define SF 0.025
  int i, j, l;
  point3f p, q, r;

  ps_OpenFile ( "grnl.ps", 600 );
  SetupCameraG ( 1800, 1350, 100, 100 );
  SetPoint3f ( &p, 0.0, 0.0, 0.0 );
  CameraProjectPoint3f ( &CPos, &p, &r );
  SetPoint3f ( &p, 1.0, 0.0, 0.0 );
  CameraProjectPoint3f ( &CPos, &p, &q );
  psl_SetLine ( r.x, r.y, q.x, q.y, 0.0, 1.0  );
  ps_Set_RGB ( 1.0, 0.0, 0.0 );
  psl_Draw ( 0.0, 0.95, 4.0 );
  psl_Arrow ( 1.0, true );
  ps_Set_RGB ( 0.0, 1.0, 0.0 );
  SetPoint3f ( &p, 0.0, 1.0, 0.0 );
  CameraProjectPoint3f ( &CPos, &p, &q );
  psl_SetLine ( r.x, r.y, q.x, q.y, 0.0, 1.0  );
  psl_Draw ( 0.0, 0.95, 4.0 );
  psl_Arrow ( 1.0, true );

/*
  ps_Set_Gray ( 0.68 );
  for ( i = 0; i < n; i++ ) {
    SetPoint3f ( &p, xy[i].x, xy[i].y, 0.0 );
    CameraProjectPoint3f ( &CPos, &p, &q );
    ps_Fill_Circle ( q.x, q.y, 3.0 );
  }
*/
  ps_Set_Line_Width ( 2.0 );
  for ( i = 0; i < k; i++ ) {
    ps_Set_Gray ( 0.3+0.1*i );
    for ( j = 0; j < n; j++ ) {
      SetPoint3f ( &p, xy[i*n*n+j*n].x, xy[i*n*n+j*n].y, z[i*n*n+j*n]*SF );
      CameraProjectPoint3f ( &CPos, &p, &q );
      for ( l = 1; l < n; l++ ) {
        r = q;
        SetPoint3f ( &p, xy[i*n*n+j*n+l].x, xy[i*n*n+j*n+l].y, z[i*n*n+j*n+l]*SF );
        CameraProjectPoint3f ( &CPos, &p, &q );
        ps_Draw_Line ( r.x, r.y, q.x, q.y );
      }
    }
    for ( j = 0; j < n; j++ ) {
      SetPoint3f ( &p, xy[i*n*n+j].x, xy[i*n*n+j].y, z[i*n*n+j]*SF );
      CameraProjectPoint3f ( &CPos, &p, &q );
      for ( l = 1; l < n; l++ ) {
        r = q;
        SetPoint3f ( &p, xy[i*n*n+l*n+j].x, xy[i*n*n+l*n+j].y, z[i*n*n+l*n+j]*SF );
        CameraProjectPoint3f ( &CPos, &p, &q );
        ps_Draw_Line ( r.x, r.y, q.x, q.y );
      }
    }
  }
  ps_CloseFile ();
  printf ( "grnl.ps\n" );
} /*DrawFGraph*/

/* ///////////////////////////////////////////////////////////////////////// */
int GetOption ( GHoleDomainf *domain, int query, int qn,
                int *ndata, int **idata, float **fdata )
{
  switch ( query ) {
/*
case G1HQUERY_BASIS:
    return G1H_USE_RESTRICTED_BASIS;
*/
case G1HQUERY_Q2_FORM_CONSTANT:
    *ndata = 1;
    *fdata = &q2const;
    return G1H_Q2_USE_SUPPLIED_CONSTANT;
default:
    return G1H_DEFAULT;
  }
} /*GetOption*/

int main ( void )
{
  struct tms start, stop;
  float time;
/*  float funcval[3]; */
  vector3f *acoeff;
#ifdef NEWDATA
  vector3f hvect[HOLE_K];
#endif
/*  int i; */

  pkv_InitScratchMem ( 2097152 );
  acoeff = (vector3f*)pkv_GetScratchMem ( 158*sizeof(vector3f) );
  if ( !acoeff )
    goto failure;
  InitKnots ( HOLE_K );
  InitDomain ( HOLE_K, dparam );
#ifdef NEWDATA
  InitGHVectors3f ( HOLE_K, hvect );
#endif
  times ( &start );
  if ( (domain = gh_CreateDomainf ( HOLE_K, knots, Dompt )) ) {
    g1h_SetOptionProcf ( domain, GetOption );
    if ( g1h_ComputeBasisf ( domain ) ) {
      InitHole ( HOLE_K, iparams );
#ifdef NEWDATA
      InitGHSurfNetf ( HOLE_K, hvect, 5, newparams, Bpt );
#endif
      nfinalp = 0;
      if ( !g1h_Q2NLFillHolef ( domain, Bpt, (float*)acoeff, NULL, OutPatch ) )
        goto failure;
/*
      if ( !g1h_NLFunctionalValuef ( domain, Bpt, acoeff, &funcval[1] ) )
        goto failure;
      printf ( "E_a = %f, E_b = %f\n", funcval[1], funcval[2] );
*/
      if ( !g1h_Q2NLExtFillHolef ( domain, Bpt, (float*)acoeff, NULL, OutPatch ) )
        goto failure;
/*
      if ( !g1h_NLExtFunctionalValuef ( domain, Bpt, acoeff, &funcval[1] ) )
        goto failure;
      printf ( "E_a = %f, E_b = %f\n", funcval[1], funcval[2] );
*/
      times ( &stop );
      MakePictures ();
      time = (float)(stop.tms_utime-start.tms_utime)/
             (float)(sysconf(_SC_CLK_TCK));
      printf ( "time = %8.3f\n", time );
    }
/*
for ( i = 0; i <= 12*HOLE_K; i++ )      
  MultVector3f ( 10.0, &Bpt[i], &Bpt[i] );      
nfinalp = 0;
if ( !g1h_NLExtFillHolef ( domain, Bpt, (float*)acoeff, OutPatch ) )
  goto failure;
if ( !g1h_NLExtFunctionalValuef ( domain, Bpt, acoeff, &funcval[1] ) )   
  goto failure;
printf ( "E_a = %f, E_b = %f\n", funcval[1], funcval[2] );
*/
    gh_DestroyDomainf ( domain );
  }
/*
  pkn_MultMatrixNumf ( 1, 2*(12*HOLE_K+1), 0, &Dompt[0].x, 2.0,
                       0, &Dompt[0].x );
  if ( (domain = gh_CreateDomainf ( HOLE_K, knots, Dompt )) ) {
    g1h_SetOptionProcf ( domain, GetOption );
    if ( g1h_ComputeBasisf ( domain ) ) {
      InitHole ( HOLE_K, iparams );
#ifdef NEWDATA
      InitGHSurfNetf ( HOLE_K, hvect, 5, newparams, Bpt );
#endif
      nfinalp = 0;
      if ( !g1h_Q2NLFillHolef ( domain, Bpt, (float*)acoeff, OutPatch ) )
        goto failure;
      if ( !g1h_Q2NLExtFillHolef ( domain, Bpt, (float*)acoeff, OutPatch ) )
        goto failure;
      times ( &stop );
      time = (float)(stop.tms_utime-start.tms_utime)/
             (float)(sysconf(_SC_CLK_TCK));
      printf ( "time = %8.3f\n", time );
    }
    gh_DestroyDomainf ( domain );
  }
*/

  printf ( "Scratch memory used: %d bytes\n", pkv_MaxScratchTaken() );
  pkv_DestroyScratchMem ();
  exit ( 0 );

failure:
  exit ( 1 );
} /*main*/

