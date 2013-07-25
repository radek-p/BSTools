
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/times.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkgeom.h"
#include "multibs.h"
#include "g1blendingd.h"
#include "convh.h"
#include "camerad.h"
#include "xgedit.h"
#include "xgeipc.h"

#include "render.h"
#include "spl3d.h"
#include "ed3ds.h"

/* ////////////////////////////////////////////////////////////////////////// */
void InitG1BlendingSurface ( void )
{
  int    i, j, k;
  double xh, yh;

  sw_clamped_blending = false;
  sw_blending_constraints = false;
        /* nonclosed at the moment */
  kwind.closed_u = kwind.closed_v = false;
        /* the surface is supposed to be bicubic */
  degree_u = degree_v = 2;
        /* set equidistant knots */
  lastknot_u = max ( lastknot_u, 7 );
  SetEquidistantU ();
  lastknot_v = max ( lastknot_v, 7 );
  SetEquidistantV ();
        /* the control net is initially quite regular and flat */
  xh = 2.0/(double)(lastknot_u-3);
  yh = 2.0/(double)(lastknot_v-3);
  for ( i = k = 0;  i < lastknot_u-2;  i++ )
    for ( j = 0;  j < lastknot_v-2;  j++, k++ )
      SetPoint4d ( &cpoints[k], (double)i*xh-1.0, (double)j*yh-1.0, 0.0, 1.0 );

        /* now setup the projections and project the data */
  FindBoundingBox ( &swind.RefBBox );
  swind.PerspBBox = swind.RefBBox;
  xge_3DwindSetupParProj ( &swind, &swind.RefBBox );
  xge_3DwindSetupPerspProj ( &swind, false );
  for ( i = 0; i < 4; i++ )
    ProjectSurface ( i );
  xge_T2KnotWindFindMapping ( &kwind );
  ClearPointMarking ( (lastknot_u-2)*(lastknot_v-2), mkpoints );
} /*InitG1BlendingSurface*/

void RefineG1BlendingConstrKnots ( void )
{
/*
  int i;

  blending_constr_knots[0] = 0.0;
  for ( i = 1; i <= n_blending_constraints; i++ ) {
    blending_constr_knots[i] = 2.0*blending_constr_knots[i] - 3.0;
    FindG1BlendingConstrCP ( i );
  }
  for ( i = 0; i < 4; i++ )
    ProjectG1BlendingConstrCPoly ( i );
*/
} /*RefineG1BlendingConstrKnots*/

boolean RefineG1BlendingSurface ( void )
{
  void    *sp;
  int     nlknu, nlknv, pitch1, pitch2;
  double  *auxcp1, *auxcp2;
  boolean clamped;

  sp = pkv_GetScratchMemTop ();
  if ( degree_u != 2 || degree_v != 2 )
    goto failure;

        /* compute the final numbers of knots */
  nlknu = 2*(lastknot_u - degree_u);
  nlknv = 2*(lastknot_v - degree_v);
  if ( nlknu >= MAX_KNOTS || nlknv >= MAX_KNOTS )
    goto failure;
  pitch1 = 4*(lastknot_v-degree_v);
  pitch2 = 4*(nlknv-degree_v);
  auxcp1 = pkv_GetScratchMemd ( pitch1*(nlknu-degree_u) );
  auxcp2 = pkv_GetScratchMemd ( pitch2*(nlknu-degree_u) );
  if ( !auxcp1 || !auxcp2 )
    goto failure;
  clamped = sw_clamped_blending;
  G1ClampedBoundaryToFree ();
  sw_clamped_blending = false;
        /* use the Lane-Riesenfeld algorithm */
  if ( !mbs_multiLaneRiesenfeldd ( 4, lastknot_u-degree_u,
                             degree_v, lastknot_v,
                             pitch1, (double*)cpoints, &nlknv, pitch2, auxcp1 ) )
    goto failure;
  if ( !mbs_multiLaneRiesenfeldd ( pitch2, 1, degree_u, lastknot_u,
                             0, auxcp1, &nlknu, 0, auxcp2 ) )
    goto failure;
        /* success; store the result in the arrays */
  memset ( mkpoints, 0, pitch2*(nlknu-degree_u) );
  memcpy ( cpoints, auxcp2, pitch2*(nlknu-degree_u)*sizeof(double) );
  lastknot_u = nlknu;
  lastknot_v = nlknv;
  SetEquidistantU ();
  SetEquidistantV ();
  sw_triharmonic_blending = blending_mat_valid = false;
  sw_bind_blending = sw_nonlin_blending = false;
  if ( clamped ) {
    G1FreeBoundaryToClamped ();
    sw_clamped_blending = true;
  }
  RefineG1BlendingConstrKnots ();
  ResizeObject ();
  SetKWindNKnots ();
  xge_T2KnotWindFindMapping ( &kwind );

  blending_opt_part[0] = 2*blending_opt_part[0]-3;
  blending_opt_part[1] = 2*blending_opt_part[1]+1;
  blending_opt_part[2] = 2*blending_opt_part[2]-3;
  blending_opt_part[3] = 2*blending_opt_part[3]+1;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*RefineG1BlendingSurface*/

void G1CheckIfClamped ( void )
{
  int i;

  sw_clamped_blending = false;

  for ( i = 2; i < degree_u; i++ )
    if ( knots_u[i-1] < knots_u[i] )
      return;
  for ( i = lastknot_u-degree_u+1; i < lastknot_u; i++ )
    if ( knots_u[i-1] < knots_u[i] )
      return;
  for ( i = 2; i < degree_v; i++ )
    if ( knots_v[i-1] < knots_v[i] )
      return;
  for ( i = lastknot_v-degree_v+1; i < lastknot_v; i++ )
    if ( knots_v[i-1] < knots_v[i] )
      return;

  sw_clamped_blending = true;
} /*G1CheckIfClamped*/

void SetupG1BLKnots ( void )
{
  double  *knots;
  int     lkn, i;
  boolean clu, clv;

  knots = knots_u;
  lkn = lastknot_u;
  clu = knots[1] >= knots[2] && knots[lkn-2] >= knots[lkn-1];
  mbs_TransformAffKnotsd ( 2, lkn, knots, (double)2, (double)(lkn-2), knots );
  for ( i = 0; i <= lkn; i++ )
    knots[i] = (double)i;

  knots = knots_v;
  lkn = lastknot_v;
  clv = knots[1] >= knots[2] && knots[lkn-2] >= knots[lkn-1];
  mbs_TransformAffKnotsd ( 2, lkn, knots, (double)2, (double)(lkn-2), knots );
  for ( i = 0; i <= lkn; i++ )
    knots[i] = (double)i;

  if ( clu && clv ) {
    knots = knots_u;
    lkn = lastknot_u;
    knots[0] = knots[1] = knots[2];
    knots[lkn] = knots[lkn-1] = knots[lkn-2];
    knots = knots_v;
    lkn = lastknot_v;
    knots[0] = knots[1] = knots[2];
    knots[lkn] = knots[lkn-1] = knots[lkn-2];
    sw_clamped_blending = true;
  }
  else
    sw_clamped_blending = false;
} /*SetupG1BLKnots*/

boolean SetupG1BlendingMatrix ( void )
{
  if ( blending_Amat ) { free ( blending_Amat );  blending_Amat = NULL; }
  if ( blending_Arow ) { free ( blending_Arow );  blending_Arow = NULL; }
  if ( blending_prof ) { free ( blending_prof );  blending_prof = NULL; }
  blending_mat_valid = false;
  if ( degree_u != 2 || degree_v != 2 ||
       lastknot_u < 7 || lastknot_v < 7 )
    return false;
  blending_lknu = blending_opt_part[1]-blending_opt_part[0]+7;
  blending_lknv = blending_opt_part[3]-blending_opt_part[2]+7;
  if ( blending_lknu < 7 || blending_lknv < 7 )
    return false;
  if ( g1bl_SetupBiharmAMatrixd ( blending_lknu, blending_lknv,
                               &blending_n, &blending_prof, &blending_Amat,
                               &blending_Arow ) ) {
    blending_mat_valid = pkn_NRBSymCholeskyDecompd ( blending_n,
                               blending_prof, blending_Amat, blending_Arow, NULL );
  }
  return blending_mat_valid;
} /*SetupG1BlendingMatrix*/

boolean ConstructBiharmG1BlendingSurface ( void )
{
  void   *sp;
  double *rhs, *cp;
  int    fcp, pitch;

  if ( !blending_mat_valid )
    return false;
  sp = pkv_GetScratchMemTop ();
  rhs = pkv_GetScratchMemd ( 4*blending_n );
  if ( !rhs )
    goto failure;
printf ( ".1\n" );
  pitch = 4*(lastknot_v-2);
  fcp = (blending_opt_part[0]-2)*(lastknot_v-2)+(blending_opt_part[2]-2);
  cp = &cpoints[fcp].x;
  if ( !g1bl_SetupBiharmRHSd ( blending_lknu, blending_lknv,
                               4, pitch, cp, rhs ) )
    goto failure;
  if ( !pkn_NRBLowerTrSolved ( blending_n, blending_prof,
                               blending_Amat, blending_Arow,
                               4, 4, rhs, 4, rhs ) )
    goto failure;
  if ( !pkn_NRBUpperTrSolved ( blending_n, blending_prof,
                               blending_Amat, blending_Arow,
                               4, 4, rhs, 4, rhs ) )
    goto failure;
  pkv_Selectd ( blending_lknu-6, 4*(blending_lknv-6),
                4*(blending_lknv-6), pitch, rhs, &cp[2*pitch+8] );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*ConstructBiharmG1BlendingSurface*/

boolean ConstructBiharmG1BlendingConstrSurface ( void )
{
#ifdef IT_HAS_FINALLY_BEEN_IMPLEMENTED
  void   *sp;
  int    w, bls, i, nconstr;
  double *w1T, *w1, *w2;
  double *rhs, *cp, *rhsc, *arhs, *ww;
  int    fcp, pitch;

FILE *f;

  if ( !blending_mat_valid )
    return false;
  sp = pkv_GetScratchMemTop ();
  ww = pkv_GetScratchMemd ( n_blending_constraints );
  if ( !ww )
    goto failure;
  ValidateG1ConstrKnots ();
  for ( i = 1, w = 0;  i <= n_blending_constraints;  i++ ) {
    if ( blending_constr_poly_valid[i] )
      w++;
    ww[i-1] = blending_constr_knots[i]-(double)(blending_opt_part[0]-3);
  }

  if ( !w ) {  /* no valid constraints */
    pkv_SetScratchMemTop ( sp );
    return ConstructBiharmG1BlendingSurface ();
  }

  pitch = 4*(lastknot_v-3);
  fcp = (blending_opt_part[0]-3)*(lastknot_v-3)+(blending_opt_part[2]-3);
  cp = &cpoints[fcp].x;

        /* setup the constraint equations matrix */
  bls = blending_lknv-9;
  w1  = pkv_GetScratchMemd ( w*bls*blending_n );
  w1T = pkv_GetScratchMemd ( w*bls*blending_n );
  w2  = pkv_GetScratchMemd ( w*bls*max(bls*6,blending_n) );
  rhsc = pkv_GetScratchMemd ( 4*w*bls );
  if ( !w1 || !w1T || !w2 || !rhsc )
    goto failure;

printf ( ",1\n" );
  if ( !g2bl_SetupULConstraintsd ( blending_lknu, blending_lknv, 4,
                                   pitch, cp, n_blending_constraints, ww,
                                   pitch, (double*)blending_constr_cp,
                                   &nconstr, w1, rhsc ) )
    goto failure;

        /* solve the system L*E = w_1^T */
  pkv_TransposeMatrixd ( nconstr, blending_n, blending_n, w1, nconstr, w1T );
  if ( !pkn_NRBLowerTrSolved ( blending_n, blending_prof,
                               blending_Amat, blending_Arow,
                               nconstr, nconstr, w1T, nconstr, w2 ) )
    goto failure;

        /* find the QR decomposition of the matrix E */
  arhs = pkv_GetScratchMemd ( max(2*nconstr, 4*blending_n) );
  if ( !arhs )
    goto failure;
  if ( !pkn_QRDecomposeMatrixd ( blending_n, nconstr, w2, arhs ) )
    goto failure;

        /* setup the right hand side, -B*b */
  rhs = pkv_GetScratchMemd ( 4*blending_n );
  if ( !rhs )
    goto failure;
  if ( !g2bl_SetupBiharmRHSd ( blending_lknu, blending_lknv,
                                4, pitch, cp, rhs ) )
    goto failure;

        /* solve the systems of equations L*f = -B*b, L^T*e = f */
  if ( !pkn_NRBLowerTrSolved ( blending_n, blending_prof,
                               blending_Amat, blending_Arow,
                               4, 4, rhs, 4, rhs ) )
    goto failure;
  if ( !pkn_NRBUpperTrSolved ( blending_n, blending_prof,
                               blending_Amat, blending_Arow,
                               4, 4, rhs, 4, rhs ) )
    goto failure;

        /* solve the systems of equations R^T*g = C-W*e, R*c = g */
  pkn_MultMatrixSubd ( nconstr, blending_n, blending_n, w1, 4, 4, rhs, 4, rhsc );
  pkn_multiMultInvTrUTVectord ( nconstr, w2, 4, 4, rhsc, 4, rhsc );
  pkn_multiMultInvUTVectord ( nconstr, w2, 4, 4, rhsc, 4, rhsc );

        /* solve the systems of equations L*h = W^T*c, L^T*k = h */
  pkn_MultMatrixd ( blending_n, nconstr, nconstr, w1T, 4, 4, rhsc, 4, arhs );
  if ( !pkn_NRBLowerTrSolved ( blending_n, blending_prof,
                               blending_Amat, blending_Arow,
                               4, 4, arhs, 4, arhs ) )
    goto failure;
  if ( !pkn_NRBUpperTrSolved ( blending_n, blending_prof,
                               blending_Amat, blending_Arow,
                               4, 4, arhs, 4, arhs ) )
    goto failure;

        /* compute a = e + h and insert the computed points to the control net */
  pkn_AddMatrixd ( 1, 4*blending_n, 0, rhs, 0, arhs, 0, rhs );
  pkv_Selectd ( blending_lknu-9, 4*(blending_lknv-9),
                4*(blending_lknv-9), pitch, rhs, &cp[3*pitch+12] );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
#endif
  return false;
} /*ConstructBiharmG1BlendingConstrSurface*/

boolean ConstructG1BlendingSurface ( void )
{
  boolean result;

  SetupG1BLKnots ();
  sw_nonlin_blending = false;
  if ( !blending_mat_valid )
    SetupG1BlendingMatrix ();
  if ( !blending_mat_valid )
    return false;

  if ( sw_clamped_blending )
    G1ClampedBoundaryToFree ();
  if ( sw_blending_constraints && n_blending_constraints ) {
    result = ConstructBiharmG1BlendingConstrSurface ();
  }
  else {
    result = ConstructBiharmG1BlendingSurface ();
  }

  if ( sw_clamped_blending )
    G1FreeBoundaryToClamped ();
  ResizeObject ();
  return result;
} /*ConstructG1BlendingSurface*/

/* ////////////////////////////////////////////////////////////////////////// */
/* the two procedures below are intended to be used only for quadratic        */
/* blending surfaces with equidistant (subsequent integers) internal knots    */
void G1FreeBoundaryToClamped ( void )
{
  void   *sp;
  double newkn[4];
  int    i, pitch;

  if ( degree_u != 2 || degree_v != 2 )
    return;
  sp = pkv_GetScratchMemTop ();
  pitch = 4*(lastknot_v-degree_v);
  newkn[0] = 0.0;
  for ( i = 1; i <= 2; i++ )
    newkn[i] = 2.0;
  mbs_multiBSChangeLeftKnotsd ( lastknot_u-degree_u, 4,
                                degree_v, knots_v,
                                pitch, (double*)cpoints, newkn );
  for ( i = 0; i < 2; i++ )
    newkn[i] = (double)(lastknot_v-2);
  newkn[2] = (double)lastknot_v;
  mbs_multiBSChangeRightKnotsd ( lastknot_u-degree_u, 4,
                                 degree_v, lastknot_v, knots_v,
                                 pitch, (double*)cpoints, newkn );
  newkn[0] = 0.0;
  for ( i = 1; i <= 2; i++ )
    newkn[i] = 2.0;
  mbs_multiBSChangeLeftKnotsd ( 1, pitch, degree_u,
                                knots_u, 0, (double*)cpoints, newkn );
  for ( i = 0; i < 2; i++ )
    newkn[i] = (double)(lastknot_u-2);
  newkn[2] = (double)lastknot_u;
  mbs_multiBSChangeRightKnotsd ( 1, pitch, degree_u, lastknot_u,
                                 knots_u, 0, (double*)cpoints, newkn );

  pkv_SetScratchMemTop ( sp );
} /*G1FreeBoundaryToClamped*/

void G1ClampedBoundaryToFree ( void )
{
  void   *sp;
  double newkn[4];
  int    i, pitch;

  if ( degree_u != 2 || degree_v != 2 )
    return;
  sp = pkv_GetScratchMemTop ();
  pitch = 4*(lastknot_v-degree_v);
  for ( i = 0; i <= 2; i++ )
    newkn[i] = (double)i;
  mbs_multiBSChangeLeftKnotsd ( 1, pitch, degree_u,
                                knots_u, 0, (double*)cpoints, newkn );
  for ( i = 0; i <= 2; i++ )
    newkn[i] = (double)(lastknot_u-2+i);
  mbs_multiBSChangeRightKnotsd ( 1, pitch, degree_u, lastknot_u,
                                 knots_u, 0, (double*)cpoints, newkn );
  for ( i = 0; i <= 2; i++ )
    newkn[i] = (double)i;
  mbs_multiBSChangeLeftKnotsd ( lastknot_u-degree_u, 4,
                                degree_v, knots_v,
                                pitch, (double*)cpoints, newkn );
  for ( i = 0; i <= 2; i++ )
    newkn[i] = (double)(lastknot_v-2+i);
  mbs_multiBSChangeRightKnotsd ( lastknot_u-degree_u, 4,
                                 degree_v, lastknot_v, knots_v,
                                 pitch, (double*)cpoints, newkn );

  pkv_SetScratchMemTop ( sp );
} /*G1ClampedBoundaryToFree*/

