
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
#include "g2ps.h"

char fn[] = "g2hdomain.ps";

int hole_k;

GHoleDomainf *domain;
CameraRecf PPos;

static void SetupCamera ( void )
{
  CameraInitFramef ( &PPos, true, false, 1400, 1400, 0, 0, 1.0, 4 );
  CameraInitPosf ( &PPos );
  PPos.vd.para.diag = 14.0;
  CameraSetMappingf ( &PPos );
} /*SetupCamera*/

static void MakePicture ( char *fn, boolean net )
{
  int i;

  if ( net )
    ps_WriteBBox ( 4, 16, 333, 169 );
  else
    ps_WriteBBox ( 0, 7, 181, 204 );
  ps_OpenFile ( fn, 600 );
  ps_DenseScreen ();
  ps_Write_Command ( "1 setlinecap" );
  for ( i = 0; i < hole_k; i++ )
    DrawBSDomPatch ( hole_k, knots, Dompt, &PPos, i, true, false );
  for ( i = 0; i < hole_k; i++ )
    DrawBSDomPatch ( hole_k, knots, Dompt, &PPos, i, false, false );
  if ( net )
    for ( i = 0; i < hole_k; i++ )
      DrawBSDomNet ( hole_k, Dompt, &PPos, i, false );
  ps_GSave ();
  ps_Write_Command ( "1400 0 translate" );
  DrawBSDomNetNum ( hole_k, Dompt, &PPos );
  ps_GRestore ();

  ps_CloseFile ();
  printf ( "%s\n", fn );
} /*MakePicture*/

int main ( void )
{
  float dparams[3] = { 0.0, 0.0, 0.0 };

  pkv_InitScratchMem ( 2097152 );
  SetupCamera ();
  InitKnots ( hole_k = 5 );
  knots[2] = 0.23;
  knots[43-2] = (float)(1.0-knots[2]);
  knots[8] = 0.96;
  knots[9] = 0.98;
  knots[32-8] = (float)(1.0-knots[8]);
  knots[32-9] = (float)(1.0-knots[9]);
  InitDomain ( hole_k, dparams );
  hole_np = 12*hole_k+1;
  if ( (domain = gh_CreateDomainf ( hole_k, knots, Dompt )) ) {
    if ( !g2h_ComputeBasisf ( domain ) )
      exit ( 1 );
    MakePicture ( fn, true );
    gh_DestroyDomainf ( domain );
  }
  else
    exit ( 1 );

  printf ( "Scratch memory used: %d bytes\n", (int)pkv_MaxScratchTaken() );
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

