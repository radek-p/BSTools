
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
#include "writdatf.h"

#define _DRAW_BASIS
#define _DRAW_LAPLACIAN

#define WIDTH  920.0
#define HEIGHT 690.0

#define QUAD_FACTOR 10  /* 10 */

#define NK 1
#define M1 1
#define M2 1

#define HOLE_K 5
#define FUNCA21

char fn[] = "eg2hnlsplf.ps";
char datfn1[] = "eg2hnlsplf.dat";

GHoleDomainf *domain;

int nfpatches = 0;
int fpdeg, fplkn;
float *fpknots;
point3f *fpcp[HOLE_K];

int nfbpatches = 0;
int bfpdeg;
point3f *bfpcp[2*HOLE_K];

/* ///////////////////////////////////////////////////////////////////////// */
void SetupCPos ( float w, float h, float x, float y )
{
  vector3f v;

  CameraInitFramef ( &CPos, false, true, w, h, x, y, 1.0, 4 );
  CameraInitPosf ( &CPos );

  SetVector3f ( &v, 0.0, 0.0, -60.0 );
  CameraMoveGf ( &CPos, &v );
  CameraRotXGf ( &CPos, 0.25*PI );
  CameraRotZGf ( &CPos, 0.2*PI );
  CameraZoomf ( &CPos, 12.0 );
/*
  SetVector3f ( &v, 0.0, 0.0, -46.0 );
  CameraMoveGf ( &CPos, &v );
  CameraRotXGf ( &CPos, -0.65*PI );
  CameraRotZGf ( &CPos, 0.05*PI ); 
  CameraZoomf ( &CPos, 11.0 );
*/
} /*SetupCPos*/

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

void DrawBezPatches ( void )
{
  int     i, j;
  point3f cp[16];

  ps_Set_Gray ( 0.5 );
  for ( i = 0; i < HOLE_K; i++ )
    for ( j = 0; j < 3; j++ ) {
      GetHoleSurrndPatch ( i, j, cp );
      DrawBezPatch3f ( 3, 3, cp, 6, 6 );
    }
} /*DrawBezPatches*/

void DrawFinalBezPatches ( void )
{
  int i;

  ps_Set_RGB ( 1.0, 0.0, 0.0 );
  for ( i = 0; i < nfpatches; i++ )
    DrawBezPatch3f ( fpdeg, fpdeg, bfpcp[i], 2*(NK+1), 2*(NK+1) );
} /*DrawFinalBezPatches*/

void DrawFinalPatches ( void )
{
  int i;

  ps_Set_Gray ( 0.25 );
  for ( i = 0; i < nfpatches; i++ )
    DrawBSPatch3f ( fpdeg, fplkn, fpknots, fpdeg, fplkn, fpknots,
                    fpcp[i], 4, 4 );
} /*DrawFinalPatches*/

void MakePicture1 ( void )
{
  ps_OpenFile ( fn, 600 );
  ps_Write_Command ( "1 setlinecap" );
  SetupCPos ( 2*WIDTH, 2*HEIGHT, 0.0, 0.0 );
  ps_Set_Line_Width ( 1.0 );
  ps_Draw_Rect ( CPos.width, CPos.height, 0, 0 );
  ps_Set_Clip_Rect ( CPos.width, CPos.height, 0, 0 );
  DrawBezPatches ();
  DrawFinalBezPatches ();
  DrawFinalPatches ();

  ps_CloseFile ();
  printf ( "%s\n", fn );
} /*MakePicture1*/

/* ///////////////////////////////////////////////////////////////////////// */
void MakePictures ( void )
{
  MakePicture1 ();
} /*MakePictures*/

/* ///////////////////////////////////////////////////////////////////////// */
void WriteData ( void )
{
  int i;
  vector3f ld = {0.2,-0.3,1.0};

  OpenDataFile ( datfn1 );
  SetupCPos ( WIDTH, HEIGHT, 0.0, 0.0 );
  WriteCamera ( &CPos );
  NormalizeVector3f ( &ld );
  WriteLightDir ( &ld );
  WriteKnots ( HOLE_K, (void*)knots );
  WriteDomCP ( 12*HOLE_K+1, Dompt );
  WriteSurfCP ( 12*HOLE_K+1, Bpt );
  WriteSurfPatches ( HOLE_K, (void*)knots, Bpt );
  for ( i = 0; i < nfpatches; i++ )
    WriteBSPatch (  fpdeg, fplkn, fpknots, fpdeg, fplkn, fpknots, fpcp[i] );
  CloseDataFile ();
  printf ( "%s\n", datfn1 );
/*
  OpenDataFile ( datfn2 );
  SetupCPos ( WIDTH, HEIGHT, 0.0, 0.0 );
  WriteCamera ( &CPos );
  WriteLightDir ( &ld );
  WriteKnots ( HOLE_K, (void*)knots );
  WriteDomCP ( 12*HOLE_K+1, Dompt );
  WriteSurfCP ( 12*HOLE_K+1, Bpt );
  WriteSurfPatches ( HOLE_K, (void*)knots, Bpt );
  for ( i = 0; i < HOLE_K; i++ )
    WriteBezPatch (  fpdeg, fpdeg, &bfpcp[HOLE_K+i]->x );
  CloseDataFile ();
  printf ( "%s\n", datfn2 );

  OpenDataFile ( datfn3 );
  SetupCPos ( WIDTH, HEIGHT, 0.0, 0.0 );
  WriteCamera ( &CPos );
  WriteLightDir ( &ld );
  WriteKnots ( HOLE_K, (void*)knots );
  WriteDomCP ( 12*HOLE_K+1, Dompt );
  WriteSurfCP ( 12*HOLE_K+1, Bpt );
  WriteSurfPatches ( HOLE_K, (void*)knots, Bpt );
  for ( i = 0; i < HOLE_K; i++ )
    WriteBezPatch (  fpdeg, fpdeg, &bfpcp[i]->x );
  CloseDataFile ();
  printf ( "%s\n", datfn3 );
*/
} /*WriteData*/

/* ///////////////////////////////////////////////////////////////////////// */
void OutputSplPatch ( int degu, int lknu, const float *knotsu,
                      int degv, int lknv, const float *knotsv,
                      const point3f *cp, void *usrptr )
{
  int s;

  if ( nfpatches == 0 ) {
    fpknots = malloc ( (lknu+1)*sizeof(float) );
    if ( !fpknots )
      exit ( 1 );
    memcpy ( fpknots, knotsu, (lknu+1)*sizeof(float) );
    fpdeg = degu;
    fplkn = lknu;
    printf ( "%d %d %d %d\n", degu, lknu, degv, lknv );
  }
  s = (lknu-degu)*(lknu-degu);
  fpcp[nfpatches] = malloc ( s*sizeof(point3f) );
  if ( !fpcp[nfpatches] )
    exit ( 1 );
  memcpy ( fpcp[nfpatches], cp, s*sizeof(point3f) );
  nfpatches ++;
} /*OutputSplPatch*/

void OutputBezPatch ( int degu, int degv, const float *cp, void *usrptr )
{
  bfpdeg = degu;
  bfpcp[nfbpatches] = malloc ( (degu+1)*(degu+1)*sizeof(point3f) );
  if ( !bfpcp[nfbpatches] )
    exit ( 1 );
  memcpy ( bfpcp[nfbpatches], cp, (degu+1)*(degu+1)*sizeof(point3f) );
  nfbpatches ++;
} /*OutputBezPatch*/

int main ( void )
{
  struct tms start, stop;
  float time;
  float param[2] = { 0.0, 0.0 };
  float iparams[3] = /*{ 0.0, 0.0, 0.0 }*/ { 0.7, 0.1, 1.0 };
  float *acoeff, fval;
  int nfunc_a, nfunc_b, nfunc_c, nfunc_d;

  setvbuf ( stdout, NULL, _IONBF, 0 );
  pkv_InitScratchMem ( 67108864 ); /* 64MB */
  acoeff = pkv_GetScratchMemf ( 262144 );
  if ( !acoeff )
    exit ( 1 );
  InitKnots ( HOLE_K );
  ModifyKnots ( HOLE_K );
  InitDomain ( HOLE_K, param );
  InitHole ( HOLE_K, iparams );
  times ( &start );
  if ( (domain = gh_CreateDomainf ( HOLE_K, knots, Dompt )) ) {
    if ( g2h_ComputeSplBasisf ( domain, NK, M1, M2 ) ) {
      g2h_DrawSplBasFuncNumf ( domain, &nfunc_a, &nfunc_b, &nfunc_c, &nfunc_d );
      printf ( "a=%d, b=%d, c=%d, d=%d\n", nfunc_a, nfunc_b, nfunc_c, nfunc_d );
      printf ( "    A    B    C    D\n" );
      printf ( "%5d%5d%5d%5d%5d\n", 0, nfunc_a, nfunc_a+nfunc_b,
               nfunc_a+nfunc_b+nfunc_c, nfunc_a+nfunc_b+nfunc_c+nfunc_d );
      printf ( "%5d     %5d%5d\n", nfunc_c+nfunc_d, 0, nfunc_c );
      nfbpatches = 0;
      if ( !g2h_FillHolef ( domain, 3, (float*)Bpt, acoeff, NULL, OutputBezPatch ) )
        exit ( 1 );
      fval = g2h_FunctionalValuef ( domain, 3, (float*)Bpt, acoeff );
      printf ( "fval = %f\n", fval );
      if ( !g2h_ExtFillHolef ( domain, 3, (float*)Bpt, acoeff, NULL, OutputBezPatch ) )
        exit ( 1 );
      fval = g2h_ExtFunctionalValuef ( domain, 3, (float*)Bpt, acoeff );
      printf ( "fval = %f\n", fval );
      nfpatches = 0;
      if ( g2h_NLSplFillHolef ( domain, Bpt, NULL, NULL, OutputSplPatch ) ) {
        times ( &stop );
        time = (float)(stop.tms_utime-start.tms_utime)/
               (float)(sysconf(_SC_CLK_TCK));
        printf ( "time = %8.3f\n", time );
        MakePictures ();
        WriteData ();
      }
    }
    gh_DestroyDomainf ( domain );
  }

  printf ( "Scratch memory used: %d bytes\a\n", pkv_MaxScratchTaken() );
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

