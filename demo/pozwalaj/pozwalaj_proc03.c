
/* ///////////////////////////////////////////////////  */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* ///////////////////////////////////////////////////  */

#include <sys/types.h>
#include <sys/times.h>
#include <signal.h>
#include <unistd.h>
#include <setjmp.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fpu_control.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "eg2holed.h"
#include "bsmesh.h"
#include "g1blendingd.h"
#include "g2blendingd.h"
#include "g2mblendingd.h"
#include "xgedit.h"
#include "xgeipc.h"

#define CHILD_SIDE
#include "pozwalajipc.h"
#undef CHILD_SIDE

#include "pozwalaj_proc.h"

/* /////////////////////////////////////////////////////////////////////////// */
void ReadBSMSize ( int isize )
{
  read ( xge_pipe_in[0], &bsmsize, sizeof(ipc_bsm_size) );
} /*ReadBSMSize*/

void ReadBSMVert ( int isize )
{
  if ( meshv ) PKV_FREE ( meshv );
  PKV_MALLOC ( meshv, isize );
  if ( meshv )
    read ( xge_pipe_in[0], meshv, isize );
} /*ReadBSMVert*/

void ReadBSMVertHE ( int isize )
{
  if ( meshvhei ) PKV_FREE ( meshvhei );
  PKV_MALLOC ( meshvhei, isize );
  if ( meshvhei )
    read ( xge_pipe_in[0], meshvhei, isize );
} /*ReadBSMVertHE*/

void ReadBSMVertPC ( int isize )
{
  if ( meshvpc ) PKV_FREE ( meshvpc );
  if ( _meshvpc ) PKV_FREE ( _meshvpc );
  PKV_MALLOC ( meshvpc, isize );
  if ( meshvpc ) {
    read ( xge_pipe_in[0], meshvpc, isize );
    meshvpcsize = isize;
  }
} /*ReadBSMVertPC*/

void ReadBSMVertMK ( int isize )
{
  if ( meshmkcp ) PKV_FREE ( meshmkcp );
  PKV_MALLOC ( meshmkcp, isize );
  if ( meshmkcp )
    read ( xge_pipe_in[0], meshmkcp, isize );
} /*ReadBSMVertMK*/

void ReadBSMHalfedges ( int isize )
{
  if ( meshhe ) PKV_FREE ( meshhe );
  PKV_MALLOC ( meshhe, isize );
  if ( meshhe )
    read ( xge_pipe_in[0], meshhe, isize );
} /*ReadBSMHalfedges*/

void ReadBSMFacets ( int isize )
{
  if ( meshfac ) PKV_FREE ( meshfac );
  PKV_MALLOC ( meshfac, isize );
  if ( meshfac )
    read ( xge_pipe_in[0], meshfac, isize );
} /*ReadBSMFacets*/

void ReadBSMFacetHE ( int isize )
{
  if ( meshfhei ) PKV_FREE ( meshfhei );
  PKV_MALLOC ( meshfhei, isize );
  if ( meshfhei )
    read ( xge_pipe_in[0], meshfhei, isize );
} /*ReadBSMFacetHE*/

void ReadBSMOptimizeOptions ( int isize )
{
  read ( xge_pipe_in[0], &bsmoptions, sizeof(ipc_bsm_options) );
    /* at this point everything ought to be read in */
  xge_CallTheParent ( ipccmd_BEGIN_BSM, 0 );
  xge_ChildCallYourself ( ipccmd_BEGIN_BSM );
} /*ReadBSMOptimizeOptions*/

/* /////////////////////////////////////////////////////////////////////////// */
void ReadBSMCSize ( int isize )
{
  read ( xge_pipe_in[0], &cmeshsize, sizeof(ipc_bsm_size) );
} /*ReadBSMCSize*/

void ReadBSMCVert ( int isize )
{
  if ( cmeshv ) PKV_FREE ( cmeshv );
  PKV_MALLOC ( cmeshv, isize );
  if ( cmeshv )
    read ( xge_pipe_in[0], cmeshv, isize );
} /*ReadBSMCVert*/

void ReadBSMCVertHE ( int isize )
{
  if ( cmeshvhei ) PKV_FREE ( cmeshvhei );
  PKV_MALLOC ( cmeshvhei, isize );
  if ( cmeshvhei )
    read ( xge_pipe_in[0], cmeshvhei, isize );
} /*ReadBSMCVertHE*/

void ReadBSMCHalfedges ( int isize )
{
  if ( cmeshhe ) PKV_FREE ( cmeshhe );
  PKV_MALLOC ( cmeshhe, isize );
  if ( cmeshhe )
    read ( xge_pipe_in[0], cmeshhe, isize );
} /*ReadBSMCHalfedges*/

void ReadBSMCFacets ( int isize )
{
  if ( cmeshfac ) PKV_FREE ( cmeshfac );
  PKV_MALLOC ( cmeshfac, isize );
  if ( cmeshfac )
    read ( xge_pipe_in[0], cmeshfac, isize );
} /*ReadBSMCFacets*/

void ReadBSMCFacetHE ( int isize )
{
  if ( cmeshfhei ) PKV_FREE ( cmeshfhei );
  PKV_MALLOC ( cmeshfhei, isize );
  if ( cmeshfhei )
    read ( xge_pipe_in[0], cmeshfhei, isize );
} /*ReadBSMCFacetHE*/

