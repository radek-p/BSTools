
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
#include "g2blendingd.h"
#include "convh.h"
#include "camerad.h"
#include "xgedit.h"

#include "render.h"
#include "spl3d.h"
#include "ed3ds.h"
#include "ed3dswidgets.h"

/* ////////////////////////////////////////////////////////////////////////// */
void InitG2BlendingSurface ( void )
{
  int    i, j, k;
  double xh, yh;

  sw_clamped_blending = false;
  sw_blending_constraints = false;
        /* nonclosed at the moment */
  kwind.closed_u = kwind.closed_v = false;
        /* the surface is supposed to be bicubic */
  degree_u = degree_v = 3;
        /* set equidistant knots */
  lastknot_u = max ( lastknot_u, 10 );
  SetEquidistantU ();
  lastknot_v = max ( lastknot_v, 10 );
  SetEquidistantV ();
        /* the control net is initially quite regular and flat */
  xh = 2.0/(double)(lastknot_u-4);
  yh = 2.0/(double)(lastknot_v-4);
  for ( i = k = 0;  i < lastknot_u-3;  i++ )
    for ( j = 0;  j < lastknot_v-3;  j++, k++ )
      SetPoint4d ( &cpoints[k], (double)i*xh-1.0, (double)j*yh-1.0, 0.0, 1.0 );

        /* now setup the projections and project the data */
  FindBoundingBox ( &swind.RefBBox );
  swind.PerspBBox = swind.RefBBox;
  xge_3DwindSetupParProj ( &swind, &swind.RefBBox );
  xge_3DwindSetupPerspProj ( &swind, false );
  for ( i = 0; i < 4; i++ )
    ProjectSurface ( i );
  xge_T2KnotWindFindMapping ( &kwind );
  ClearPointMarking ( (lastknot_u-3)*(lastknot_v-3), mkpoints );
} /*InitG2BlendingSurface*/

void RefineG2BlendingConstrKnots ( void )
{
  int i;

  blending_constr_knots[0] = 0.0;
  for ( i = 1; i <= n_blending_constraints; i++ ) {
    blending_constr_knots[i] = 2.0*blending_constr_knots[i] - 3.0;
    FindG2BlendingConstrCP ( i );
  }
  for ( i = 0; i < 4; i++ )
    ProjectG2BlendingConstrCPoly ( i );
} /*RefineG2BlendingConstrKnots*/

boolean RefineG2BlendingSurface ( void )
{
  void    *sp;
  int     nlknu, nlknv, pitch1, pitch2;
  double  *auxcp1, *auxcp2;
  boolean clamped;

  sp = pkv_GetScratchMemTop ();
  if ( degree_u != 3 || degree_v != 3 )
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
  G2ClampedBoundaryToFree ();
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
  if ( kwind.closed_u ) {
    if ( (blending_opt_part[0] == 0 &&
          blending_opt_part[1] == lastknot_u-2*degree_u-1) ||
         (blending_opt_part[0] == blending_opt_part[1]+1) ) {
      blending_opt_part[0] = 0;
      blending_opt_part[1] = nlknu-2*degree_u-1;
    }
    else {
      blending_opt_part[0] = 2*blending_opt_part[0]-3;
      blending_opt_part[1] = 2*blending_opt_part[1]+1;
      if ( blending_opt_part[0] < 0 )
        blending_opt_part[0] += nlknu-2*degree_u;
      else if ( blending_opt_part[0] > nlknu-2*degree_u )
        blending_opt_part[0] -= nlknu-2*degree_u;
      if ( blending_opt_part[1] < 0 )
        blending_opt_part[1] += nlknu-2*degree_u;
      else if ( blending_opt_part[1] > nlknu-2*degree_u )
        blending_opt_part[1] -= nlknu-2*degree_u;
    }
  }
  else {
    blending_opt_part[0] = 2*blending_opt_part[0]-3;
    blending_opt_part[1] = 2*blending_opt_part[1]+1;
  }
  blending_opt_part[2] = 2*blending_opt_part[2]-3;
  blending_opt_part[3] = 2*blending_opt_part[3]+1;

  memset ( mkpoints, 0, pitch2*(nlknu-degree_u) );
  memcpy ( cpoints, auxcp2, pitch2*(nlknu-degree_u)*sizeof(double) );
  lastknot_u = nlknu;
  lastknot_v = nlknv;
  SetEquidistantU ();
  SetEquidistantV ();
  sw_triharmonic_blending = blending_mat_valid = false;
  sw_bind_blending = sw_nonlin_blending = false;
  if ( clamped ) {
    G2FreeBoundaryToClamped ();
    sw_clamped_blending = true;
  }
  if ( kwind.closed_u ) {
    kwind.clcKu = lastknot_u-2*degree_u;
    kwind.clcTu = (double)(lastknot_u-2*degree_u);
  }
  RefineG2BlendingConstrKnots ();
  ResizeObject ();
  SetKWindNKnots ();
  xge_T2KnotWindFindMapping ( &kwind );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*RefineG2BlendingSurface*/

void G2CheckIfClamped ( void )
{
  int i;

  sw_clamped_blending = false;

  for ( i = 2; i < degree_v; i++ )
    if ( knots_v[i-1] < knots_v[i] )
      return;
  for ( i = lastknot_v-degree_v+1; i < lastknot_v; i++ )
    if ( knots_v[i-1] < knots_v[i] )
      return;
  if ( kwind.closed_u ) {
    sw_clamped_blending = true;
    return ;
  }
  for ( i = 2; i < degree_u; i++ )
    if ( knots_u[i-1] < knots_u[i] )
      return;
  for ( i = lastknot_u-degree_u+1; i < lastknot_u; i++ )
    if ( knots_u[i-1] < knots_u[i] )
      return;
  sw_clamped_blending = true;
} /*G2CheckIfClamped*/

void SetupG2BLKnots ( void )
{
  double  *knots;
  int     lkn, i;
  boolean clu, clv;

  knots = knots_u;
  lkn = lastknot_u;
  if ( kwind.closed_u )
    clu = true;
  else
    clu = knots[1] >= knots[2] && knots[2] >= knots[3] &&
          knots[lkn-3] >= knots[lkn-2] && knots[lkn-2] >= knots[lkn-1];
  mbs_TransformAffKnotsd ( 3, lkn, knots, (double)3, (double)(lkn-3), knots );
  for ( i = 0; i <= lkn; i++ )
    knots[i] = (double)i;

  kwind.closed_v = false;
  knots = knots_v;
  lkn = lastknot_v;
  clv = knots[1] >= knots[2] && knots[2] >= knots[3] &&
         knots[lkn-3] >= knots[lkn-2] && knots[lkn-2] >= knots[lkn-1];
  mbs_TransformAffKnotsd ( 3, lkn, knots, (double)3, (double)(lkn-3), knots );
  for ( i = 0; i <= lkn; i++ )
    knots[i] = (double)i;

  if ( clu && clv ) {
    if ( !kwind.closed_u ) {
      knots = knots_u;
      lkn = lastknot_u;
      knots[0] = knots[1] = knots[2] = knots[3];
      knots[lkn] = knots[lkn-1] = knots[lkn-2] = knots[lkn-3];
    }
    knots = knots_v;
    lkn = lastknot_v;
    knots[0] = knots[1] = knots[2] = knots[3];
    knots[lkn] = knots[lkn-1] = knots[lkn-2] = knots[lkn-3];
    sw_clamped_blending = true;
  }
  else
    sw_clamped_blending = false;
} /*SetupG2BLKnots*/

boolean SetupG2BlendingMatrix ( void )
{
  if ( blending_Amat ) { free ( blending_Amat );  blending_Amat = NULL; }
  if ( blending_Arow ) { free ( blending_Arow );  blending_Arow = NULL; }
  if ( blending_prof ) { free ( blending_prof );  blending_prof = NULL; }
  blending_mat_valid = false;
  if ( degree_u != 3 || degree_v != 3 ||
       lastknot_u < 10 || lastknot_v < 10 )
    return false;
  if ( kwind.closed_u ) {
/* for now only the entire patch */
    blending_lknu = lastknot_u;
    blending_lknv = blending_opt_part[3]-blending_opt_part[2]+10;
    if ( blending_lknu < 10 || blending_lknv < 10 )
      return false;
    if ( g2bl_SetupClosedTriharmAMatrixd ( blending_lknu, blending_lknv,
                                 &blending_n, &blending_prof, &blending_Amat,
                                 &blending_Arow ) ) {
      blending_mat_valid = pkn_NRBSymCholeskyDecompd ( blending_n,
                                 blending_prof, blending_Amat, blending_Arow, NULL );
    }
  }
  else {
    blending_lknu = blending_opt_part[1]-blending_opt_part[0]+10;
    blending_lknv = blending_opt_part[3]-blending_opt_part[2]+10;
    if ( blending_lknu < 10 || blending_lknv < 10 )
      return false;
    if ( g2bl_SetupTriharmAMatrixd ( blending_lknu, blending_lknv,
                                 &blending_n, &blending_prof, &blending_Amat,
                                 &blending_Arow ) ) {
      blending_mat_valid = pkn_NRBSymCholeskyDecompd ( blending_n,
                                 blending_prof, blending_Amat, blending_Arow, NULL );
    }
  }
  return blending_mat_valid;
} /*SetupG2BlendingMatrix*/

boolean ConstructTriharmG2BlendingSurface ( void )
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
/*printf ( ".1\n" ); */
  pitch = 4*(lastknot_v-3);
  if ( kwind.closed_u ) {
/* for now only the entire patch */
    fcp = blending_opt_part[2]-3;
    cp = &cpoints[fcp].x;
    if ( !g2bl_SetupClosedTriharmRHSd ( blending_lknu, blending_lknv,
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
    pkv_Selectd ( lastknot_u-6, 4*(blending_lknv-9), 4*(blending_lknv-9), pitch,
                  rhs, &cp[12] );
    pkv_Selectd ( 3, pitch, pitch, pitch, &cpoints[0].x,
                  &cpoints[(lastknot_u-6)*(lastknot_v-3)].x );
  }
  else {
    fcp = (blending_opt_part[0]-3)*(lastknot_v-3)+(blending_opt_part[2]-3);
    cp = &cpoints[fcp].x;
    if ( !g2bl_SetupTriharmRHSd ( blending_lknu, blending_lknv,
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
    pkv_Selectd ( blending_lknu-9, 4*(blending_lknv-9),
                  4*(blending_lknv-9), pitch, rhs, &cp[3*pitch+12] );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*ConstructTriharmG2BlendingSurface*/

boolean ConstructTriharmG2BlendingConstrSurface ( void )
{
  void   *sp;
  int    w, bls, i, nconstr;
  double *w1T, *w1, *w2;
  double *rhs, *cp, *rhsc, *arhs, *ww;
  int    fcp, pitch;

  if ( !blending_mat_valid )
    return false;
  sp = pkv_GetScratchMemTop ();
  ww = pkv_GetScratchMemd ( n_blending_constraints );
  if ( !ww )
    goto failure;
  ValidateG2ConstrKnots ();
  for ( i = 1, w = 0;  i <= n_blending_constraints;  i++ ) {
    if ( blending_constr_poly_valid[i] )
      w++;
    ww[i-1] = blending_constr_knots[i]-(double)(blending_opt_part[0]-3);
  }
  if ( !w ) {  /* no valid constraints */
    pkv_SetScratchMemTop ( sp );
    return ConstructTriharmG2BlendingSurface ();
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
  if ( !g2bl_SetupTriharmRHSd ( blending_lknu, blending_lknv,
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
  return false;
} /*ConstructTriharmG2BlendingConstrSurface*/

boolean ConstructTriharmG2BlendingConstrClosedSurface ( void )
{
  void   *sp;
  int    w, bls, i, nconstr;
  double *w1T, *w1, *w2;
  double *rhs, *cp, *rhsc, *arhs, *ww;
  int    fcp, pitch;

  if ( !blending_mat_valid )
    return false;
  sp = pkv_GetScratchMemTop ();
  ww = pkv_GetScratchMemd ( n_blending_constraints );
  if ( !ww )
    goto failure;
  ValidateG2ConstrKnots ();
  for ( i = 1, w = 0;  i <= n_blending_constraints; i++ ) {
    if ( blending_constr_poly_valid[i] )
      w++;
    ww[i-1] = blending_constr_knots[i];
  }
  if ( !w ) { /* no valid constraints */
    pkv_SetScratchMemTop ( sp );
    return ConstructTriharmG2BlendingSurface ();
  }

  pitch = 4*(lastknot_v-3);
  fcp = blending_opt_part[2]-3; /* for now */
  cp = &cpoints[fcp].x;
        /* setup the constraint equations matrix */
  bls = blending_lknv-9;
  w1 = pkv_GetScratchMemd ( w*bls*blending_n );
  w1T = pkv_GetScratchMemd ( w*bls*blending_n );
  w2 = pkv_GetScratchMemd ( w*bls*max(bls*6,blending_n) );
  rhsc = pkv_GetScratchMemd ( 4*w*bls );
  if ( !w1 || !w1T || !w2 || !rhsc )
    goto failure;

  if ( !g2bl_SetupClosedULConstraintsd ( blending_lknu, blending_lknv, 4,
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
  if ( !g2bl_SetupClosedTriharmRHSd ( blending_lknu, blending_lknv,
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
  pkv_Selectd ( blending_lknu-6, 4*(blending_lknv-9),
                4*(blending_lknv-9), pitch, rhs, &cp[12] );
  memcpy ( &cpoints[(lastknot_u-6)*(lastknot_v-3)], cpoints, 3*pitch*sizeof(double) );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*ConstructTriharmG2BlendingConstrClosedSurface*/

boolean ConstructG2BlendingSurface ( void )
{
  boolean result;

  SetupG2BLKnots ();
  sw_nonlin_blending = false;
  if ( !blending_mat_valid )
    SetupG2BlendingMatrix ();
  if ( !blending_mat_valid )
    return false;

  if ( sw_clamped_blending )
    G2ClampedBoundaryToFree ();
  if ( sw_blending_constraints && n_blending_constraints ) {
    if ( kwind.closed_u )
      result = ConstructTriharmG2BlendingConstrClosedSurface ();
    else
      result = ConstructTriharmG2BlendingConstrSurface ();
  }
  else {
    result = ConstructTriharmG2BlendingSurface ();
  }

  if ( sw_clamped_blending )
    G2FreeBoundaryToClamped ();
  ResizeObject ();
  return result;
} /*ConstructG2BlendingSurface*/

/* ////////////////////////////////////////////////////////////////////////// */
/* the two procedures below are intended to be used only for cubic blending */
/* surfaces with equidistant (subsequent integers) internal knots */
void G2FreeBoundaryToClamped ( void )
{
  void   *sp;
  double newkn[4];
  int    i, pitch;

  if ( degree_u != 3 || degree_v != 3 )
    return;
  sp = pkv_GetScratchMemTop ();
  pitch = 4*(lastknot_v-degree_v);
  newkn[0] = 0.0;
  for ( i = 1; i <= 3; i++ )
    newkn[i] = 3.0;
  mbs_multiBSChangeLeftKnotsd ( lastknot_u-degree_u, 4,
                                degree_v, knots_v,
                                pitch, (double*)cpoints, newkn );
  for ( i = 0; i < 3; i++ )
    newkn[i] = (double)(lastknot_v-3);
  newkn[3] = (double)lastknot_v;
  mbs_multiBSChangeRightKnotsd ( lastknot_u-degree_u, 4,
                                 degree_v, lastknot_v, knots_v,
                                 pitch, (double*)cpoints, newkn );
  if ( !kwind.closed_u ) {
    newkn[0] = 0.0;
    for ( i = 1; i <= 3; i++ )
      newkn[i] = 3.0;
    mbs_multiBSChangeLeftKnotsd ( 1, pitch, degree_u,
                                  knots_u, 0, (double*)cpoints, newkn );
    for ( i = 0; i < 3; i++ )
      newkn[i] = (double)(lastknot_u-3);
    newkn[3] = (double)lastknot_u;
    mbs_multiBSChangeRightKnotsd ( 1, pitch, degree_u, lastknot_u,
                                   knots_u, 0, (double*)cpoints, newkn );
  }
  pkv_SetScratchMemTop ( sp );
} /*G2FreeBoundaryToClamped*/

void G2ClampedBoundaryToFree ( void )
{
  void   *sp;
  double newkn[4];
  int    i, pitch;

  if ( degree_u != 3 || degree_v != 3 )
    return;
  sp = pkv_GetScratchMemTop ();
  pitch = 4*(lastknot_v-degree_v);
  for ( i = 0; i <= 3; i++ )
    newkn[i] = (double)i;
  mbs_multiBSChangeLeftKnotsd ( 1, pitch, degree_u,
                                knots_u, 0, (double*)cpoints, newkn );
  for ( i = 0; i <= 3; i++ )
    newkn[i] = (double)(lastknot_u-3+i);
  mbs_multiBSChangeRightKnotsd ( 1, pitch, degree_u, lastknot_u,
                                 knots_u, 0, (double*)cpoints, newkn );
  for ( i = 0; i <= 3; i++ )
    newkn[i] = (double)i;
  mbs_multiBSChangeLeftKnotsd ( lastknot_u-degree_u, 4,
                                degree_v, knots_v,
                                pitch, (double*)cpoints, newkn );
  for ( i = 0; i <= 3; i++ )
    newkn[i] = (double)(lastknot_v-3+i);
  mbs_multiBSChangeRightKnotsd ( lastknot_u-degree_u, 4,
                                 degree_v, lastknot_v, knots_v,
                                 pitch, (double*)cpoints, newkn );

  pkv_SetScratchMemTop ( sp );
} /*G2ClampedBoundaryToFree*/

