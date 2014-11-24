
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
#include "mengerc.h"
#include "xgedit.h"

#define CHILD_SIDE
#include "pozwalajipc.h"

#include "pozwalaj_proc.h"

/* /////////////////////////////////////////////////////////////////////////// */
void SetupRefinementMatrix ( void )
{
  void         *sp;
  int          nv, nhe, nfac, *mvhei, *mfhei;
  BSMvertex    *mv;
  BSMhalfedge  *mhe;
  BSMfacet     *mfac;
  int          nv1, nhe1, nfac1, *mvhei1, *mfhei1;
  BSMvertex    *mv1;
  BSMhalfedge  *mhe1;
  BSMfacet     *mfac1;
  unsigned int rmnnz1, rmnnz2, nmult;
  index2       *rmnzi1, *rmnzi2;
  double       *rmnzc1, *rmnzc2;
  unsigned int *permut1, *permut2;
  int          *cols1, *cols2, *abpos;
  index2       *aikbkj;

  sp = pkv_GetScratchMemTop ();
printf ( "setting up the refinement matrix\n" );
  mv = mv1 = NULL;  mhe = mhe1 = NULL;  mfac = mfac1 = NULL;
  mvhei = mvhei1 = mfhei = mfhei1 = NULL;
  rmnzi1 = rmnzi2 = NULL;
  rmnzc1 = rmnzc2 = NULL;
  permut1 = permut2 = NULL;

  if ( cmeshsize.nv >= bsmsize.nv )
    goto failure;
  if ( rmnzi ) PKV_FREE ( rmnzi );
  if ( rmnzc ) PKV_FREE ( rmnzc );
  if ( !bsm_RefinementMatd ( 3, cmeshsize.nv, cmeshv, cmeshvhei,
                             cmeshsize.nhe, cmeshhe,
                             cmeshsize.nfac, cmeshfac, cmeshfhei,
                             &nv, &mv, &mvhei, &nhe, &mhe,
                             &nfac, &mfac, &mfhei,
                             &rmnnz, &rmnzi, &rmnzc ) )
    goto failure;
        /* the coarse mesh is not needed any more */
  PKV_FREE ( cmeshv );
  PKV_FREE ( cmeshvhei );
  PKV_FREE ( cmeshhe );
  PKV_FREE ( cmeshfac );
  PKV_FREE ( cmeshfhei );
printf ( "nv = %d, rmnnz = %d\n", nv, rmnnz );

  while ( nv < bsmsize.nv ) {
    if ( !bsm_RefinementMatd ( 3, nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                               &nv1, &mv1, &mvhei1, &nhe1, &mhe1,
                               &nfac1, &mfac1, &mfhei1,
                               (int*)&rmnnz1, &rmnzi1, &rmnzc1 ) )
      goto failure;
printf ( "nv1 = %d, rmnnz1 = %d\n", nv1, rmnnz1 );
    permut1 = (unsigned int*)pkv_GetScratchMemi ( rmnnz+rmnnz1+cmeshsize.nv+nv+2 );
    if ( !permut1 )
      goto failure;
    permut2 = &permut1[rmnnz];
    cols1 = (int*)&permut2[rmnnz1];
    cols2 = &cols1[cmeshsize.nv+1];
/* old and slow method */
/*
    rmnnz2 = pkn_SPMSizeMM ( nv1, nv, rmnnz1, rmnzi1, permut2, false,
                             cmeshsize.nv, rmnnz, rmnzi, permut1, false );
printf ( "rmnnz2 = %d\n", rmnnz2 );
    PKV_MALLOC ( rmnzi2, rmnnz2*sizeof(index2) );
    PKV_MALLOC ( rmnzc2, rmnnz2*sizeof(double) );
    if ( !rmnzi2 || !rmnzc2 )
      goto failure;
    if ( !pkn_SPMMultMMd ( nv1, nv, rmnnz1, rmnzi1, rmnzc1, permut2, true,
                           cmeshsize.nv, rmnnz, rmnzi, rmnzc, permut1, true,
                           &rmnnz2, rmnzi2, rmnzc2 ) )
      goto failure;
*/
/* new method */
    if ( !pkn_SPMCountMMnnzC ( nv1, nv, cmeshsize.nv,
                               rmnnz1, rmnzi1, permut2, cols2, false,
                               rmnnz, rmnzi, permut1, cols1, false,
                               &rmnnz2, &nmult ) )
      goto failure;
printf ( "rmnnz2 = %d\n", rmnnz2 );
    abpos = pkv_GetScratchMemi ( rmnnz2+1 );
    aikbkj = (index2*)pkv_GetScratchMem ( nmult*sizeof(index2) );
    if ( !abpos || !aikbkj )
      goto failure;
    PKV_MALLOC ( rmnzi2, rmnnz2*sizeof(index2) );
    PKV_MALLOC ( rmnzc2, rmnnz2*sizeof(double) );
    if ( !rmnzi2 || !rmnzc2 )
      goto failure;
    if ( !pkn_SPMFindMMnnzC ( nv1, nv, cmeshsize.nv,
                              rmnnz1, rmnzi1, permut2, cols2,
                              rmnnz, rmnzi, permut1, cols1,
                              rmnzi2, abpos, aikbkj ) )
      goto failure;
    pkn_SPMFastMultMMd ( rmnzc1, rmnzc, rmnnz2, abpos, aikbkj, rmnzc2 );

    PKV_FREE ( mv );
    PKV_FREE ( mvhei );
    PKV_FREE ( mhe );
    PKV_FREE ( mfac );
    PKV_FREE ( mfhei );
    mv    = mv1;     mv1    = NULL;
    mvhei = mvhei1;  mvhei1 = NULL;
    mhe   = mhe1;    mhe1   = NULL;
    mfac  = mfac1;   mfac1  = NULL;
    mfhei = mfhei1;  mfhei1 = NULL;
    nv = nv1;
    nhe = nhe1;
    nfac = nfac1;

    PKV_FREE ( rmnzi );
    PKV_FREE ( rmnzc );
    PKV_FREE ( rmnzi1 );
    PKV_FREE ( rmnzc1 );
    rmnnz = rmnnz2;
    rmnzi = rmnzi2;  rmnzi2 = NULL;
    rmnzc = rmnzc2;  rmnzc2 = NULL;

    pkv_SetScratchMemTop ( sp );
  }
        /* verify, whether the mesh topology produced by the refinement  */
        /* procedures here is identical with that of the fine mesh to be */
        /* optimized */
  if ( nv > bsmsize.nv || nhe != bsmsize.nhe || nfac != bsmsize.nfac )
    goto failure;
/*
  if ( !memcmp ( mv, meshv, nv*sizeof(BSMvertex) ) )       goto failure;
  if ( !memcmp ( mvhei, meshvhei, nhe*sizeof(int) ) )      goto failure;
  if ( !memcmp ( mhe, meshhe, nhe*sizeof(BSMhalfedge) ) )  goto failure;
  if ( !memcmp ( mfac, meshfac, nfac*sizeof(BSMfacet) ) )  goto failure;
  if ( !memcmp ( mfhei, meshfhei, nhe*sizeof(int) ) )      goto failure;
*/
        /* deallocate the arrays unnecessary from now on */
  if ( mv )    PKV_FREE ( mv );
  if ( mhe )   PKV_FREE ( mhe );
  if ( mfac )  PKV_FREE ( mfac );
  if ( mvhei ) PKV_FREE ( mvhei );
  if ( mfhei ) PKV_FREE ( mfhei );

  pkv_SetScratchMemTop ( sp );
  return;

failure:
  if ( mv1 )    PKV_FREE ( mv1 );
  if ( mhe1 )   PKV_FREE ( mhe1 );
  if ( mfac1 )  PKV_FREE ( mfac1 );
  if ( mvhei1 ) PKV_FREE ( mvhei1 );
  if ( mfhei1 ) PKV_FREE ( mfhei1 );
  if ( mv )    PKV_FREE ( mv );
  if ( mhe )   PKV_FREE ( mhe );
  if ( mfac )  PKV_FREE ( mfac );
  if ( mvhei ) PKV_FREE ( mvhei );
  if ( mfhei ) PKV_FREE ( mfhei );
  if ( rmnzi ) PKV_FREE ( rmnzi );
  if ( rmnzc ) PKV_FREE ( rmnzc );
  if ( rmnzi1 ) PKV_FREE ( rmnzi1 );
  if ( rmnzc1 ) PKV_FREE ( rmnzc1 );
  if ( rmnzi2 ) PKV_FREE ( rmnzi2 );
  if ( rmnzc2 ) PKV_FREE ( rmnzc2 );
  bsmoptions.use_coarse = false;
printf ( "Rejecting the coarse mesh\n" );
  pkv_SetScratchMemTop ( sp );
} /*SetupRefinementMatrix*/

/* /////////////////////////////////////////////////////////////////////////// */
void DumpRefinementMatrix ( int m, int n, int nnz, index2 *nzi, double *nzc )
{
  FILE *f;
  int  i;

  f = fopen ( "rmat.txt", "w+" );
  fprintf ( f, "m = %d, n = %d, nnz = %d\n", m, n, nnz );
  for ( i = 0; i < nnz; i++ )
    fprintf ( f, "%4d,%4d:%12.9f\n", nzi[i].i, nzi[i].j, nzc[i] );
  fclose ( f );
} /*DumpRefinementMatrix*/

/* /////////////////////////////////////////////////////////////////////////// */
void BeginBSMOptimization ( void )
{
  void    *sp;
  boolean result;
  int     i;
  byte    *ccp;

  sp = pkv_GetScratchMemTop ();
  finished = true;
  if ( bsmsize.cpdimen == 3 &&
       bsm_CheckMeshIntegrity ( bsmsize.nv, meshv, meshvhei,
                                bsmsize.nhe, meshhe,
                                bsmsize.nfac, meshfac, meshfhei,
                                NULL, NULL ) ) {
    time0 = times ( &start );
    PKV_MALLOC ( _meshvpc, meshvpcsize );
    if ( !meshvpc )
      goto failure;
    memcpy ( _meshvpc, meshvpc, meshvpcsize );
    if ( !InvertPretransformation ( &bsmoptions.pretrans ) )
      goto failure;
    TransformCPoints ( &bsmoptions.pretrans, 1, bsmsize.nv, 0, (point3d*)_meshvpc );
    if ( bsmoptions.use_constraints ) {
      if ( !meshmkcp )
        goto failure;
      ccp = pkv_GetScratchMem ( bsmsize.nv );
      if ( !ccp )
        goto failure;
      for ( i = 0; i < bsmsize.nv; i++ )
        ccp[i] = meshmkcp[i] & bsmoptions.constr_mask;
    }
    else
      ccp = NULL;
    _g2mbl_npthreads = bsmoptions.npthreads;
    if ( bsmoptions.nlevels > 0 ) {
      if ( bsmoptions.nlevels == 1 )
        bsmoptions.use_coarse = false;
      if ( bsmoptions.use_coarse )
        SetupRefinementMatrix ();
printf ( "Using the multilevel method" );
      if ( bsmoptions.shape_only ) {
printf ( ", shape only" );
        if ( bsmoptions.use_coarse ) {
printf ( ", with coarse mesh preconditioner\n" );
          result = g2mbl_MLSCMPOptInitd ( bsmsize.nv, meshv, meshvhei,
                      (point3d*)_meshvpc,
                      bsmsize.nhe, meshhe, bsmsize.nfac, meshfac, meshfhei, ccp,
                      cmeshsize.nv, rmnnz, rmnzi, rmnzc,
                      bsmoptions.nkn1, bsmoptions.nkn2, max(1,bsmoptions.nlevels),
                      &optdata );
        }
        else {
printf ( "\n" );
          result = g2mbl_MLSOptInitd ( bsmsize.nv, meshv, meshvhei,
                      (point3d*)_meshvpc, bsmsize.nhe, meshhe,
                      bsmsize.nfac, meshfac, meshfhei, ccp,
                      bsmoptions.nkn1, bsmoptions.nkn2, max(1,bsmoptions.nlevels),
                      &optdata );
        }
      }
      else {
        if ( bsmoptions.use_coarse ) {

printf ( ", with coarse mesh preconditioner\n" );
/*
DumpRefinementMatrix ( bsmsize.nv, cmeshsize.nv, rmnnz, rmnzi, rmnzc );
*/
          result = g2mbl_MLCMPOptInitd ( bsmsize.nv, meshv, meshvhei,
                      (point3d*)_meshvpc,
                      bsmsize.nhe, meshhe, bsmsize.nfac, meshfac, meshfhei, ccp,
                      cmeshsize.nv, rmnnz, rmnzi, rmnzc,
                      bsmoptions.C, 0.0, 0.0,
                      bsmoptions.nkn1, bsmoptions.nkn2, max(1,bsmoptions.nlevels),
                      &optdata );
        }
        else {

printf ( "\n" );

          result = g2mbl_MLOptInitd ( bsmsize.nv, meshv, meshvhei,
                      (point3d*)_meshvpc,
                      bsmsize.nhe, meshhe, bsmsize.nfac, meshfac, meshfhei,
                      ccp, bsmoptions.C, 0.0, 0.0,
                      bsmoptions.nkn1, bsmoptions.nkn2, max(1,bsmoptions.nlevels),
                      &optdata );
        }
      }
      if ( result ) {
        g2mbl_MLSetNextBlock ( optdata, bsmoptions.startbl );
        g2mbl_MLSetLogLeveld ( optdata, 2 );
      }
    }
    else if ( bsmoptions.nblocks == 1 ) {
printf ( "Using the nonblock algorithm\n" );
      result = g2mbl_InitBlSurfaceOptLMTd ( bsmsize.nv, meshv, meshvhei,
                    (point3d*)_meshvpc,
                    bsmsize.nhe, meshhe, bsmsize.nfac, meshfac,
                    meshfhei, ccp, bsmoptions.C, 0.0, 0.0,
                    bsmoptions.nkn1, bsmoptions.nkn2, &optdata );
    }
    else {
printf ( "Using the two-level procedure" );
      if ( bsmoptions.use_coarse ) {
printf ( ", with coarse mesh preconditioner\n" );
        SetupRefinementMatrix ();
        result = g2mbl_InitBlCMPSurfaceOptd ( bsmsize.nv, meshv, meshvhei,
                      (point3d*)_meshvpc, bsmsize.nhe, meshhe, bsmsize.nfac,
                      meshfac, meshfhei, ccp,
                      cmeshsize.nv, rmnnz, rmnzi, rmnzc,
                      bsmoptions.C, 0.0, 0.0,
                      bsmoptions.nkn1, bsmoptions.nkn2,
                      bsmoptions.nblocks, &optdata );
      }
      else {
printf ( "\n" );
        result = g2mbl_InitBlSurfaceOptAltBLMTd ( bsmsize.nv, meshv, meshvhei,
                      (point3d*)_meshvpc,
                      bsmsize.nhe, meshhe, bsmsize.nfac, meshfac,
                      meshfhei, ccp, bsmoptions.C, 0.0, 0.0,
                      bsmoptions.nkn1, bsmoptions.nkn2,
                      bsmoptions.nblocks, &optdata );
      }
    }
    if ( result ) {
      xge_ChildCallYourself ( ipccmd_CONTINUE_BSM );
      finished = false;
      itn = 0;
      pkv_SetScratchMemTop ( sp );
      return;
    }
  }
failure:
  xge_CallTheParent ( ipccmd_ERROR, 0 );
  pkv_SetScratchMemTop ( sp );
  return;
} /*BeginBSMOptimization*/

void PrepareBSMOutput ( void )
{
  memcpy ( meshvpc, _meshvpc, meshvpcsize );
  TransformCPoints ( &pretrans_inv, 1, bsmsize.nv, 0, (point3d*)meshvpc );
  ResetIPCBuffer ();
  IPCAppendDataItem ( ipcd_BSM_SIZE, sizeof(ipc_bsm_size), &bsmsize );
  IPCAppendDataItem ( ipcd_BSM_VERT, bsmsize.nv*sizeof(BSMvertex), meshv );
  IPCAppendDataItem ( ipcd_BSM_VHE, bsmsize.nhe*sizeof(int), meshvhei );
  IPCAppendDataItem ( ipcd_BSM_VERTC, bsmsize.nv*sizeof(point3d), meshvpc );
  if ( meshmkcp )
    IPCAppendDataItem ( ipcd_BSM_VERTMK, bsmsize.nv*sizeof(char), meshmkcp );
  IPCAppendDataItem ( ipcd_BSM_HALFE, bsmsize.nhe*sizeof(BSMhalfedge), meshhe );
  IPCAppendDataItem ( ipcd_BSM_FAC, bsmsize.nfac*sizeof(BSMfacet), meshfac );
  IPCAppendDataItem ( ipcd_BSM_FHE, bsmsize.nhe*sizeof(int), meshfhei );
} /*PrepareBSMOutput*/

void MarkBlockPoints ( boolean mark )
{
  int bln, nv, nvcp, c0, bvncp;
  int *nncpi, *vncpi, *bvncpi, *vpermut;
  int i;

  if ( !meshmkcp )
    return;
  g2mbl_GetBLMBlockNumd ( optdata, &bln );
  g2mbl_GetBLMTBlockInfod ( optdata, bln, &nv, &nvcp, &nncpi, &c0,
                            &bvncp, &vncpi, &bvncpi, &vpermut );
  for ( i = 0; i < nv; i++ )
    if ( nncpi[i] >= 0 )
      meshmkcp[i] &= ~MASK_CP_SPECIAL;
  if ( mark && bvncpi )
    for ( i = 0; i < bvncp; i++ )
      meshmkcp[bvncpi[i]] |= MASK_CP_SPECIAL;
} /*MarkBlockPoints*/

void MarkMLBlockPoints ( void )
{
  int bl, nvcp, *vncpi, seed;
  int i;

  g2mbl_MLGetBlockVCPNumbersd ( optdata, 0, &nvcp, &vncpi, &seed );
  for ( i = 0; i < nvcp; i++ )
    meshmkcp[vncpi[i]] &= ~MASK_CP_SPECIAL;
  bl = g2mbl_MLGetLastBlockd ( optdata );
  if ( bl > 0 ) {
    g2mbl_MLGetBlockVCPNumbersd ( optdata, bl, &nvcp, &vncpi, &seed );
    for ( i = 0; i < nvcp; i++ )
      meshmkcp[vncpi[i]] |= MASK_CP_SPECIAL;
  }
} /*MarkMLBlockPoints*/

void ContinueBSMOptimization ( void )
{
  boolean result;
  float   time_prep, time_h, time_cg;

  if ( finished )
    return;
  if ( bsmoptions.nlevels > 0 ) {
    if ( bsmoptions.shape_only ) {
      if ( bsmoptions.alt_multilevel )
        result = g2mbl_MLSCOptIterd ( optdata, &finished );
      else
        result = g2mbl_MLSOptIterd ( optdata, &finished );
    }
    else {
      if ( bsmoptions.alt_multilevel )
        result = g2mbl_MLCOptIterd ( optdata, &finished );
      else
        result = g2mbl_MLOptIterd ( optdata, &finished );
    }
  }
  else if ( bsmoptions.nblocks == 1 )
    result = g2mbl_IterBlSurfaceOptLMTd ( optdata, &finished );
  else
    result = g2mbl_IterBlSurfaceOptAltBLMTd ( optdata, &finished );
  if ( result ) {
    if ( bsmoptions.nlevels == 0 )
      itn ++;
    if ( finished || itn >= bsmoptions.maxit ) {
      if ( bsmoptions.nlevels > 0 ) {
#ifdef G2MBL_TIME_IT
        g2mbl_MLGetTimes ( optdata, &time_prep, &time_h, &time_cg );
        printf ( "times: prep %6.2f, h %6.2f, cg %6.2f, ",
                 time_prep, time_h, time_cg );
#endif
        MarkMLBlockPoints ();
        g2mbl_MLOptDeallocated ( &optdata );
      }
      else {
        if ( bsmoptions.nblocks > 1 )
          MarkBlockPoints ( false );
        g2mbl_OptLMTDeallocated ( &optdata );
      }
      PrepareBSMOutput ();
      xge_CallTheParent ( ipccmd_FINAL_RESULT, ipc_data_size );
      xge_ChildCallYourself ( ipccmd_SEND_RESULT );
#ifdef G2MBL_TIME_IT
      time1 = times ( &stop );
      printf ( "total %6.2fs\n",
               ((double)(time1-time0))/sysconf(_SC_CLK_TCK) );
#endif
      finished = true;
    }
    else {
      if ( bsmoptions.nlevels > 0 )
        MarkMLBlockPoints ();
      else if ( bsmoptions.nblocks > 1 )
        MarkBlockPoints ( true );
      PrepareBSMOutput ();
      xge_CallTheParent ( ipccmd_PARTIAL_RESULT, ipc_data_size );
      xge_ChildCallYourself ( ipccmd_SEND_RESULT );
      xge_ChildCallYourself ( ipccmd_CONTINUE_BSM );
    }
  }
  else {
    if ( bsmoptions.nlevels > 0 )
      g2mbl_MLOptDeallocated ( &optdata );
    else
      g2mbl_OptLMTDeallocated ( &optdata );
    xge_CallTheParent ( ipccmd_ERROR, 0 );
    finished = true;
  }
} /*ContinueBSMOptimization*/

