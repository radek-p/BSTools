
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

#include "eg2holef.h"
#include "datagenf.h"
#include "drawitf.h"

#define _DRAW_BASIS
#define _DRAW_LAPLACIAN

#define QUAD_FACTOR 10  /* 10 */

#define NK 3
#define M1 3
#define M2 5

#define HOLE_K 5
#define FUNCA21

/*
#define HOLE_K 6
#define FUNCA25
*/
char fn1temp[] = "sbasD%03df.ps";   /* this is a template for file names! */
char fn2temp[] = "sbasD%03dLf.ps";
char fn3[] = "smatrixf.ps";

GHoleDomainf *domain;
int   nfunc_a, nfunc_b, nfunc_c, nfunc_d;

/* ///////////////////////////////////////////////////////////////////////// */
void SetupPProj ( float w, float h, float x, float y, float diag )
{
  CameraInitFramef ( &PPos, true, true, w, h, x, y, 1.0, 4 );
  CameraInitPosf ( &PPos );
  PPos.vd.para.diag = diag;
  CameraSetMappingf ( &PPos );
} /*SetupPProj*/

void SetupCPos ( float w, float h, float x, float y )
{
  vector3f v;

  CameraInitFramef ( &CPos, false, true, w, h, x, y, 1.0, 4 );
  CameraInitPosf ( &CPos );
  SetVector3f ( &v, 0.0, 0.0, -46.0 );
  CameraMoveGf ( &CPos, &v );
  CameraRotXGf ( &CPos, -0.65*PI );
  CameraRotZGf ( &CPos, 0.05*PI );
  CameraZoomf ( &CPos, 11.0 );
} /*SetupCPos*/

/* ///////////////////////////////////////////////////////////////////////// */
void DrawFN ( float x, float y, int fn )
{
  char s[50];

  ps_Set_Gray ( 0.0 );
  sprintf ( s, "(%d) show", fn );
  ps_MoveTo ( x, y );
  ps_Write_Command ( s );
} /*DrawFN*/

float scf = 1.0;
int   cnt;

void SetScf ( int fn )
{
  if ( fn > nfunc_a+nfunc_b+nfunc_c ) {
     switch ( (fn - (nfunc_a+nfunc_b+nfunc_c)) % 3 ) {
  case 0: scf = 1.0;  break;
  case 1: scf = 2.0;  break;
  case 2: scf = 8.0;  break;
     }
  }
  else
    scf = 1.0;
} /*SetScf*/

void DrawAuxBSPatch ( int n, int lknu, const float *knu,
                      int m, int lknv, const float *knv, const point3f *cp )
{
  void    *sp;
  point3f *p;
  int     i, s;

  sp = pkv_GetScratchMemTop ();
  s = (lknu-n)*(lknv-m);
  p = pkv_GetScratchMem ( s*sizeof(point3f) );
  if ( !p )
    exit ( 1 );
  memcpy ( p, cp, s*sizeof(point3f) );
  for ( i = 0; i < s; i++ )
    p[i].z *= scf;
  DrawBSPatch3Rf ( n, lknu, knu, m, lknv, knv, p, 0.0, 0.35, 0.0, 1.0, 6, 4 );
  pkv_SetScratchMemTop ( sp );
} /*DrawAuxBSPatch*/

void AltDrawAuxBSPatch ( int n, int lknu, const float *knu,
                         int m, int lknv, const float *knv, const point3f *cp )
{
  void    *sp;
  point3f *p;
  int     i, s;

  if ( !cnt ) {
    sp = pkv_GetScratchMemTop ();
    s = (lknu-n)*(lknv-m);
    p = pkv_GetScratchMem ( s*sizeof(point3f) );
    if ( !p )
      exit ( 1 );
    memcpy ( p, cp, s*sizeof(point3f) );
    for ( i = 0; i < s; i++ )
      p[i].z *= scf;
    DrawBSPatch3Rf ( n, lknu, knu, m, lknv, knv, p, 0.0, 0.35, 0.0, 1.0, 6, 6 );
    pkv_SetScratchMemTop ( sp );
  }
  cnt ++;
} /*AltDrawAuxBSPatch*/

void DrawAuxBasisFuncBSPatch ( int fn, float x, float y )
{
  char s[50];

  SetScf ( fn );
  ps_GSave ();
  sprintf ( s, "%6.2f %6.2f translate", x, y );
  ps_Write_Command ( s );
  ps_Set_Gray ( 0.5 );
  g2h_DrawSplBasAuxPatchesf ( domain, fn, DrawAuxBSPatch );
  ps_GRestore ();
} /*DrawAuxBasisFuncBSPatch*/

void DrawBSPatch ( int n, int lknu, const float *knu,
                   int m, int lknv, const float *knv, const point3f *cp )
{
  void    *sp;
  point3f *p;
  int     i, s;

  sp = pkv_GetScratchMemTop ();
  s = (lknu-n)*(lknv-m);
  p = pkv_GetScratchMem ( s*sizeof(point3f) );
  if ( !p )
    exit ( 1 );
  memcpy ( p, cp, s*sizeof(point3f) );
  for ( i = 0; i < s; i++ )
    p[i].z *= scf;
  DrawBSPatch3f ( n, lknu, knu, m, lknv, knv, p, 6, 6 );
  pkv_SetScratchMemTop ( sp );
} /*DrawBSPatch*/

void EmptyDrawBSPatch ( int n, int lknu, const float *knu,
                        int m, int lknv, const float *knv, const point3f *cp )
{
} /*EmptyDrawBSPatch*/

void DrawBasisFuncBSPatch ( int fn, float x, float y )
{
  char s[50];

  SetScf ( fn );
  ps_GSave ();
  sprintf ( s, "%6.2f %6.2f translate", x, y );
  ps_Write_Command ( s );
  DrawFN ( 0.0, 1100.0, fn );
  ps_Set_Gray ( 0.5 );
  g2h_DrawSplBasFunctionf ( domain, fn, DrawBSPatch );
  if ( fn >= nfunc_a+nfunc_b+nfunc_c ) {
    cnt = -((fn-(nfunc_a+nfunc_b+nfunc_c)) / (nfunc_d/HOLE_K));
    ps_Set_RGB ( 1.0, 0.0, 0.0 );
    g2h_DrawSplBasAuxPatchesf ( domain, fn, AltDrawAuxBSPatch );
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
float maxlap;

void DrawBSPatchL1 ( int n, int lknu, const float *knu,
                     int m, int lknv, const float *knv, const point3f *cp )
{
  float lap;

  lap = FindBSPatchMaxLaplacianf ( n, lknu, knu, m, lknv, knv, cp, 20, 20 );
  maxlap = max ( maxlap, lap );
} /*DrawBSPatchL1*/

void DrawBSPatchL2 ( int n, int lknu, const float *knu,
                     int m, int lknv, const float *knv, const point3f *cp )
{
  void    *sp;
  point3f *p;
  int     i, s;

  sp = pkv_GetScratchMemTop ();
  s = (lknu-n)*(lknv-m);
  p = pkv_GetScratchMem ( s*sizeof(point3f) );
  if ( !p )
    exit ( 1 );
  memcpy ( p, cp, s*sizeof(point3f) );
  for ( i = 0; i < s; i++ )
    p[i].z *= scf;
  DrawBSPatchLaplacianf ( n, lknu, knu, m, lknv, knv, p, 6, 6 );
  pkv_SetScratchMemTop ( sp );
} /*DrawBSPatchL2*/

void DrawBasisFuncBSPatchL ( int fn, float x, float y )
{
  char s[50];

  ps_GSave ();
  sprintf ( s, "%6.2f %6.2f translate", x, y );
  ps_Write_Command ( s );
  DrawFN ( 0, 1100.0, fn );
  ps_Set_Gray ( 0.5 );
  maxlap = 0.0;
  g2h_DrawSplBasFunctionf ( domain, fn, DrawBSPatchL1 );
  scf = 0.5/maxlap;
  g2h_DrawSplBasFunctionf ( domain, fn, DrawBSPatchL2 );
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
  int i, j;

  for ( i = j = 0; j < nfunc_d; i++, j += 8 ) {
    sprintf ( fn, fn2temp, i );
    MakePicture2x ( fn, j+nfunc_a+nfunc_b+nfunc_c, min (8, nfunc_d-j) );
  }
} /*MakePicture2*/

/* ///////////////////////////////////////////////////////////////////////// */
void DrawSplMatrix ( int k, int r, int s, int t, float *A, float *B )
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
    ps_Write_Command ( "100 100 translate 1 dup scale" );
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

void MakePicture3 ( void )
{
  ps_OpenFile ( fn3, 600 );
  g2h_DrawSplMatricesf ( domain, DrawSplMatrix );
  ps_CloseFile ();
  printf ( "%s\n", fn3 );
} /*MakePicture3*/

/* ///////////////////////////////////////////////////////////////////////// */
void MakePictures ( void )
{
#ifdef DRAW_BASIS
  MakePicture1 ();
#endif
#ifdef DRAW_LAPLACIAN
  MakePicture2 ();
#endif
  MakePicture3 ();
} /*MakePictures*/

/* ///////////////////////////////////////////////////////////////////////// */
int   mk, mr, ms, mt, ms1, ms2;
float *Amat = NULL, *Bmat = NULL;
float *AAmat = NULL, *BBmat = NULL;
float *Lmat = NULL;

int     cnt;
int     sdegu1[HOLE_K], lknu1[HOLE_K], sdegv1[HOLE_K], lknv1[HOLE_K];
float   *knu1[HOLE_K], *knv1[HOLE_K];
point3f *scp1[HOLE_K];
int     sdegu2[HOLE_K], lknu2[HOLE_K], sdegv2[HOLE_K], lknv2[HOLE_K];
float   *knu2[HOLE_K], *knv2[HOLE_K];
point3f *scp2[HOLE_K];

void AllocPatchStorage ( void )
{
  int i, n;

  n = G2H_FINALDEG+1+NK*max(M1+4,M2);
  for ( i = 0; i < HOLE_K; i++ ) {
    knu1[i] = malloc ( (n+G2H_FINALDEG+1)*sizeof(float) );
    knv1[i] = malloc ( (n+G2H_FINALDEG+1)*sizeof(float) );
    scp1[i] = malloc ( n*n*sizeof(point3f) );
    knu2[i] = malloc ( (n+G2H_FINALDEG+1)*sizeof(float) );
    knv2[i] = malloc ( (n+G2H_FINALDEG+1)*sizeof(float) );
    scp2[i] = malloc ( n*n*sizeof(point3f) );
    if ( !knu1[i] || !knv1[i] || !scp1[i] ||
         !knu2[i] || !knv2[i] || !scp2[i] )
      exit ( 1 );
  }
} /*AllocPatchStorage*/

void SavePatch1 ( int n, int lknu, const float *knu,
                  int m, int lknv, const float *knv, const point3f *cp )
{
  if ( cnt < HOLE_K ) {
    sdegu1[cnt] = n;
    lknu1[cnt] = lknu;
    memcpy ( knu1[cnt], knu, (lknu+1)*sizeof(float) );
    sdegv1[cnt] = m;
    lknv1[cnt] = lknv;
    memcpy ( knv1[cnt], knv, (lknv+1)*sizeof(float) );
    memcpy ( scp1[cnt], cp, (lknu-n)*(lknv-m)*sizeof(point3f) );
    cnt++;
  }
  else
    exit ( 1 );
} /*SavePatch1*/

void SavePatch2 ( int n, int lknu, const float *knu,
                  int m, int lknv, const float *knv, const point3f *cp )
{
  if ( cnt < HOLE_K ) {
    sdegu2[cnt] = n;
    lknu2[cnt] = lknu;
    memcpy ( knu2[cnt], knu, (lknu+1)*sizeof(float) );
    sdegv2[cnt] = m;
    lknv2[cnt] = lknv;
    memcpy ( knv2[cnt], knv, (lknv+1)*sizeof(float) );
    memcpy ( scp2[cnt], cp, (lknu-n)*(lknv-m)*sizeof(point3f) );
    cnt++;
  }
  else
    exit ( 1 );
} /*SavePatch2*/

void SaveMatrices ( int k, int r, int s, int t, float *A, float *B )
{
  mk = k;  mr = r;  ms = s;  mt = t;
  ms1 = pkn_Block2ArraySize ( k, r, s, t );
  ms2 = (k*(r+s)+t)*(6*k+1);
  printf ( "s1=%d, s2=%d\n", ms1, ms2 );
  Amat = malloc ( ms1*sizeof(float) );
  Bmat = malloc ( ms2*sizeof(float) );
  if ( !Amat || !Bmat )
    exit ( 1 );
  memcpy ( Amat, A, ms1*sizeof(float) );
  memcpy ( Bmat, B, ms2*sizeof(float) );
} /*SaveMatrices*/

boolean IsZnonZero ( int n, int lknu, const float *knu,
                     int m, int lknv, const float *knv, const point3f *cp )
{
  int ncp, i;

  ncp = (lknu-n)*(lknv-m);
  for ( i = 0; i < ncp; i++ )
    if ( cp[i].z )
      return true;
  return false;
} /*IsZnonZero*/

static boolean dumpspl = false;

float CompIntegral ( void )
{
  void     *sp;
  int      i, j, k, sz, d1, d2;
  double   s, ss, sss;
  float    u, v;
  point3f  *ap1, *ap2, p1[16], p2[16];
  vector2f lgr1, lgr2;
  float    jac1, jac2, gx, gy, gxx, gxy, gyy, gxxx, gxxy, gxyy, gyyy;
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
  ap1 = pkv_GetScratchMem ( 4*sz*sizeof(point3f) );
  ap2 = pkv_GetScratchMem ( 4*sz*sizeof(point3f) );
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
        u = (float)(j+j+1)/(float)(2*QUAD_FACTOR*(NK+1));
        mbs_multideBoorDer3f ( sdegu1[i], lknu1[i], knu1[i], 1, 3*d1, 0,
            (float*)scp1[i], u, &ap1[0].x, &ap1[d1].x, &ap1[2*d1].x, &ap1[3*d1].x );
        mbs_multideBoorDer3f ( sdegu2[i], lknu2[i], knu2[i], 1, 3*d2, 0,
            (float*)scp2[i], u, &ap2[0].x, &ap2[d2].x, &ap2[2*d2].x, &ap2[3*d2].x );
        s = 0.0;
        for ( k = 0; k < QUAD_FACTOR*(NK+1); k++ ) {
          v = (float)(k+k+1)/(float)(2*QUAD_FACTOR*(NK+1));
          mbs_multideBoorDer3f ( sdegv1[i], lknv1[i], knv1[i], 4, 3, 3*d1,
              (float*)ap1, v, &p1[0].x, &p1[4].x, &p1[8].x, &p1[12].x );
          mbs_multideBoorDer3f ( sdegv2[i], lknv2[i], knv2[i], 4, 3, 3*d2,
              (float*)ap2, v, &p2[0].x, &p2[4].x, &p2[8].x, &p2[12].x );
          pkn_Comp2iDerivatives3f ( p1[1].x, p1[1].y, p1[4].x, p1[4].y,
              p1[2].x, p1[2].y, p1[5].x, p1[5].y, p1[8].x, p1[8].y,
              p1[3].x, p1[3].y, p1[6].x, p1[6].y, p1[9].x, p1[9].y,
              p1[12].x, p1[12].y, 1, &p1[1].z, &p1[4].z, &p1[2].z, &p1[5].z,
              &p1[8].z, &p1[3].z, &p1[6].z, &p1[9].z, &p1[12].z,
              &gx, &gy, &gxx, &gxy, &gyy, &gxxx, &gxxy, &gxyy, &gyyy );
          lgr1.x = gxxx+gxyy;
          lgr1.y = gxxy+gyyy;
          jac1 = fabs(p1[4].x*p1[1].y-p1[4].y*p1[1].x);
          pkn_Comp2iDerivatives3f ( p2[1].x, p2[1].y, p2[4].x, p2[4].y,
              p2[2].x, p2[2].y, p2[5].x, p2[5].y, p2[8].x, p2[8].y,
              p2[3].x, p2[3].y, p2[6].x, p2[6].y, p2[9].x, p2[9].y,
              p2[12].x, p2[12].y, 1, &p2[1].z, &p2[4].z, &p2[2].z, &p2[5].z,
              &p2[8].z, &p2[3].z, &p2[6].z, &p2[9].z, &p2[12].z,
              &gx, &gy, &gxx, &gxy, &gyy, &gxxx, &gxxy, &gxyy, &gyyy );
          lgr2.x = gxxx+gxyy;
          lgr2.y = gxxy+gyyy;
          jac2 = fabs(p2[4].x*p2[1].y-p2[4].y*p2[1].x);
          s += (lgr1.x*lgr2.x+lgr1.y*lgr2.y)*jac1;
          if ( dumpspl )
            fprintf ( f, "%d, %e, %e, %e, %e, %e\n",
                      j*QUAD_FACTOR*(NK+1)+k,
                      lgr1.x, lgr1.y, lgr2.x, lgr2.y, jac1 );
        }
        ss += s;
      }
      sss += ss;
    }
  }
  pkv_SetScratchMemTop ( sp );
  sss *= 1.0/(float)(QUAD_FACTOR*QUAD_FACTOR*(NK+1)*(NK+1));
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

  memset ( AAmat, 0, ms1*sizeof(float) );
  memset ( BBmat, 0, ms2*sizeof(float) );
  for ( i = 0; i < nfunc_a+nfunc_c+nfunc_d; i++ ) {
    if ( i < nfunc_c+nfunc_d ) fni = nfunc_a+nfunc_b+i;
    else fni = i-nfunc_c-nfunc_d;
    printf ( "%4d %4d\b\b\b\b\b\b\b\b\b", i, fni );
    cnt = 0;
    g2h_DrawSplBasFunctionf ( domain, fni, SavePatch1 );
    for ( j = 0; j <= i; j++ ) {
    if ( j < nfunc_c+nfunc_d ) fnj = nfunc_a+nfunc_b+j;
    else fnj = j-nfunc_c-nfunc_d;
      k = pkn_Block2FindElemPos ( mk, mr, ms, mt, i, j );
      if ( k >= 0 ) {  /* this element belongs to a nonzero block of A */
        cnt = 0;
        g2h_DrawSplBasFunctionf ( domain, fnj, SavePatch2 );
        AAmat[k] = CompIntegral ();
      }
    }
    for ( j = 0; j < nfunc_b; j++ ) {
      cnt = 0;
      g2h_DrawSplBasFunctionf ( domain, nfunc_a+j, SavePatch2 );
      BBmat[nfunc_b*i+j] = CompIntegral ();
    }
  }
  printf ( "\n" );
} /*ComputeMatrices*/

void CompareMatrices ( void )
{
  FILE *cmp;
  int i;
  float e;

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
  g2h_DrawSplMatricesf ( domain, SaveMatrices );
  WriteBlockAddr ();
  AAmat = malloc ( ms1*sizeof(float) );
  BBmat = malloc ( ms2*sizeof(float) );
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
  float *x, *y, ca, pr, yl, lambda1, lambda2;
  int   i;

  sp = pkv_GetScratchMemTop ();
  n = mk*(mr+ms)+mt;
  x = pkv_GetScratchMemf ( n );
  y = pkv_GetScratchMemf ( n );
  if ( !x || !y )
    exit ( 1 );

  for ( x[0] = 1.0, i = 1; i < n; i++ )
    x[i] = x[i-1]*-1.0001;
  pr = pkn_SecondNormf ( n, x );
  pkn_MultMatrixNumf ( 1, n, 0, x, 1.0/pr, 0, x );
  for ( i = 0; i < 1000; i++ ) {
    pkn_Block2SymMatrixMultf ( mk, mr, ms, mt, Amat, 1, 1, x, 1, y );
    pr = pkn_ScalarProductf ( n, x, y );
    yl = pkn_SecondNormf ( n, y );
    ca = pr/yl;
    pkn_MultMatrixNumf ( 1, n, 0, y, 1.0/yl, 0, x );
    if ( ca > 0.9999999 )
      break;
  }
  printf ( "i = %d, ca = %f, yl = %f\n", i, ca, yl );
  lambda1 = yl;

  for ( x[0] = 1.0, i = 1; i < n; i++ )
    x[i] = x[i-1]*-1.0001;
  pr = pkn_SecondNormf ( n, x );
  pkn_MultMatrixNumf ( 1, n, 0, x, 1.0/pr, 0, x );
  for ( i = 0; i < 1000; i++ ) {
    memcpy ( y, x, n*sizeof(float) );
    pkn_Block2LowerTrMSolvef ( mk, mr, ms, mt, Lmat, 1, 1, y );
    pkn_Block2UpperTrMSolvef ( mk, mr, ms, mt, Lmat, 1, 1, y );
    pr = pkn_ScalarProductf ( n, x, y );
    yl = pkn_SecondNormf ( n, y );
    ca = pr/yl;
    pkn_MultMatrixNumf ( 1, n, 0, y, 1.0/yl, 0, x );
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
  g2h_DrawSplMatricesf ( domain, SaveMatrices );
  WriteBlockAddr ();
  Lmat = malloc ( ms1*sizeof(float) );
  if ( !Lmat )
    exit ( 1 );
  memcpy ( Lmat, Amat, ms1*sizeof(float) );
  if ( pkn_Block2CholeskyDecompMf ( mk, mr, ms, mt, Lmat ) ) {
    FindCondNumber ();
  }
  else
    printf ( "Cannot decompose." );
  free ( Lmat );
  free ( Amat );
  free ( Bmat );
} /*TestMatrix2*/

/* ///////////////////////////////////////////////////////////////////////// */
int main ()
{
  struct tms start, stop;
  float time;
  float param[2] = {0.0,0.0};

  setvbuf ( stdout, NULL, _IONBF, 0 );
  pkv_InitScratchMem ( 16777216 ); /* 16MB */
  InitKnots ( HOLE_K );
  ModifyKnots ( HOLE_K );
  InitDomain ( HOLE_K, param );
  times ( &start );
  if ( (domain = gh_CreateDomainf ( HOLE_K, knots, Dompt )) ) {
    if ( g2h_ComputeSplBasisf ( domain, NK, M1, M2 ) ) {
      if ( g2h_ComputeSplFormMatrixf ( domain ) ) {
        times ( &stop );
        time = (float)(stop.tms_utime-start.tms_utime)/
               (float)(sysconf(_SC_CLK_TCK));
        printf ( "time = %8.3f\n", time );

g2h_DrawSplBasFunctionf ( domain, 552, EmptyDrawBSPatch );

        g2h_DrawSplBasFuncNumf ( domain, &nfunc_a, &nfunc_b, &nfunc_c, &nfunc_d );
        printf ( "k = %d, fa = %d, fb = %d, fc = %d, fd = %d\n",
                 HOLE_K, nfunc_a, nfunc_b, nfunc_c, nfunc_d );
        MakePictures ();
        times ( &start );
        TestMatrix2 ();
        times ( &stop );
        time = (float)(stop.tms_utime-start.tms_utime)/
               (float)(sysconf(_SC_CLK_TCK));
        printf ( "time = %8.3f\n", time );
      }
    }
    gh_DestroyDomainf ( domain );
  }

  printf ( "Scratch memory used: %d bytes\a\n", pkv_MaxScratchTaken() );
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

