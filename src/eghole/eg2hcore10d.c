
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2010                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

boolean g2h_DecomposeMatrixd ( GHoleDomaind *domain )
{
  void   *sp;
  G2HolePrivateRecd *privateG2;
  double *lmat;
  int    nfunc_a;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  if ( !privateG2->Lmat ) {
    if ( !privateG2->Amat )
      if ( !g2h_ComputeFormMatrixd ( domain ) )
        goto failure;
    nfunc_a = privateG2->nfunc_a;
    lmat = privateG2->Lmat = malloc ( (nfunc_a*(nfunc_a+1)/2)*sizeof(double) );
    if ( !lmat )
      goto failure;
    memcpy ( lmat, privateG2->Amat, (nfunc_a*(nfunc_a+1)/2)*sizeof(double) );
    if ( !pkn_CholeskyDecompd ( nfunc_a, lmat ) ) {
      domain->error_code = G2H_ERROR_NONPOSITIVE_MATRIX;
      goto failure;
    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_DecomposeMatrixd*/

boolean _g2h_SetRightSided ( GHoleDomaind *domain,
                             int spdimen, CONST_ double *hole_cp,
                             double *fc00, double *b )
{
  void   *sp;
  GHolePrivateRecd  *privateG;
  G2HolePrivateRecd *privateG2;
  unsigned char *bfcpn;
  double *bc00, *bc01, *bc02, *bc10, *bc11, *bc12,
         *bd00, *bd01, *bd02, *bd10, *bd11, *bd12;
  double        *fc01, *fc02, *fc10, *fc11, *fc12,
         *fd00, *fd01, *fd02, *fd10, *fd11, *fd12;
  double *x, *y, *cp;
  int    hole_k, nfunc_a, nfunc_b;
  int    i, j, k, l;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG  = domain->privateG;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  bfcpn   = privateG->bfcpn;

  x = pkv_GetScratchMemd ( spdimen*nfunc_b );
  y = pkv_GetScratchMemd ( spdimen*(G2H_FINALDEG+1)*(G2H_FINALDEG+1) );
  if ( !x || !y )
    goto failure;

  memset ( fc00, 0, (G2_CROSSDEGSUM+6)*2*hole_k*spdimen*sizeof(double) );
  G2GetFCAddresses ();

  for ( j = 0; j < nfunc_b; j++ ) {
    l = bfcpn[j];
    cp = &hole_cp[l*spdimen];
    memcpy ( &x[j*spdimen], cp, spdimen*sizeof(double) );
    for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
      _g2h_GetBFBPatchCurvesd ( domain, j, i,
                                &bc00, &bc01, &bc02, &bc10, &bc11, &bc12,
                                &bd00, &bd01, &bd02, &bd10, &bd11, &bd12 );

      pkn_MultMatrixd ( G2_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS00DEG+1), 0, &fc00[k*(G2_CROSS00DEG+1)],
                       0, y, 0, &fc00[k*(G2_CROSS00DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS01DEG+1, 1, 1, bc01, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS01DEG+1), 0, &fc01[k*(G2_CROSS01DEG+1)],
                       0, y, 0, &fc01[k*(G2_CROSS01DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS02DEG+1, 1, 1, bc02, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS02DEG+1), 0, &fc02[k*(G2_CROSS02DEG+1)],
                       0, y, 0, &fc02[k*(G2_CROSS02DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS10DEG+1, 1, 1, bc10, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS10DEG+1), 0, &fc10[k*(G2_CROSS10DEG+1)],
                       0, y, 0, &fc10[k*(G2_CROSS10DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS11DEG+1, 1, 1, bc11, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS11DEG+1), 0, &fc11[k*(G2_CROSS11DEG+1)],
                       0, y, 0, &fc11[k*(G2_CROSS11DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS12DEG+1, 1, 1, bc12, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS12DEG+1), 0, &fc12[k*(G2_CROSS12DEG+1)],
                       0, y, 0, &fc12[k*(G2_CROSS12DEG+1)] );

      pkn_MultMatrixd ( G2_CROSS00DEG+1, 1, 1, bd00, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS00DEG+1), 0, &fd00[k*(G2_CROSS00DEG+1)],
                       0, y, 0, &fd00[k*(G2_CROSS00DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS01DEG+1, 1, 1, bd01, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS01DEG+1), 0, &fd01[k*(G2_CROSS01DEG+1)],
                       0, y, 0, &fd01[k*(G2_CROSS01DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS02DEG+1, 1, 1, bd02, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS02DEG+1), 0, &fd02[k*(G2_CROSS02DEG+1)],
                       0, y, 0, &fd02[k*(G2_CROSS02DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS10DEG+1, 1, 1, bd10, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS10DEG+1), 0, &fd10[k*(G2_CROSS10DEG+1)],
                       0, y, 0, &fd10[k*(G2_CROSS10DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS11DEG+1, 1, 1, bd11, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS11DEG+1), 0, &fd11[k*(G2_CROSS11DEG+1)],
                       0, y, 0, &fd11[k*(G2_CROSS11DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS12DEG+1, 1, 1, bd12, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS12DEG+1), 0, &fd12[k*(G2_CROSS12DEG+1)],
                       0, y, 0, &fd12[k*(G2_CROSS12DEG+1)] );
    }
  }

  if ( b )
    pkn_MultMatrixd ( nfunc_a, nfunc_b, nfunc_b, privateG2->Bmat,
                      spdimen, spdimen, x, spdimen, b );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_SetRightSided*/

boolean _g2h_OutputPatchesd ( GHoleDomaind *domain,
                int spdimen, CONST_ double *x, double *fc00, void *usrptr,
                void (*outpatch) ( int n, int m, const double *cp, void *usrptr ) )
{
  void   *sp;
  G2HolePrivateRecd *privateG2;
  double *bc00, *bc01, *bc02, *bd00, *bd01, *bd02;
  double        *fc01, *fc02, *fc10, *fc11, *fc12,
         *fd00, *fd01, *fd02, *fd10, *fd11, *fd12;
  double *y, *cp;
  int    hole_k, nfunc_a, degu, degv;
  int    i, j, k;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;

  y = pkv_GetScratchMemd ( spdimen*(G2H_FINALDEG+1)*(G2H_FINALDEG+1) );
  if ( !y )
    goto failure;

  G2GetFCAddresses ();

  for ( j = 0; j < nfunc_a; j++ ) {
    cp = &x[j*spdimen];
    for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
      _g2h_GetBFAPatchCurvesd ( domain, j, i,
                                &bc00, &bc01, &bc02, &bd00, &bd01, &bd02 );

      pkn_MultMatrixd ( G2_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixd ( 1, spdimen*(G2_CROSS00DEG+1), 0, &fc00[k*(G2_CROSS00DEG+1)],
                            0, y, 0, &fc00[k*(G2_CROSS00DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS01DEG+1, 1, 1, bc01, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixd ( 1, spdimen*(G2_CROSS01DEG+1), 0, &fc01[k*(G2_CROSS01DEG+1)],
                            0, y, 0, &fc01[k*(G2_CROSS01DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS02DEG+1, 1, 1, bc02, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixd ( 1, spdimen*(G2_CROSS02DEG+1), 0, &fc02[k*(G2_CROSS02DEG+1)],
                            0, y, 0, &fc02[k*(G2_CROSS02DEG+1)] );

      pkn_MultMatrixd ( G2_CROSS00DEG+1, 1, 1, bd00, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixd ( 1, spdimen*(G2_CROSS00DEG+1), 0, &fd00[k*(G2_CROSS00DEG+1)],
                            0, y, 0, &fd00[k*(G2_CROSS00DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS01DEG+1, 1, 1, bd01, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixd ( 1, spdimen*(G2_CROSS01DEG+1), 0, &fd01[k*(G2_CROSS01DEG+1)],
                            0, y, 0, &fd01[k*(G2_CROSS01DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS02DEG+1, 1, 1, bd02, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixd ( 1, spdimen*(G2_CROSS02DEG+1), 0, &fd02[k*(G2_CROSS02DEG+1)],
                            0, y, 0, &fd02[k*(G2_CROSS02DEG+1)] );
    }
  }

  for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
    mbs_BezC2CoonsToBezd ( spdimen,
                   G2_CROSS00DEG, &fc00[k*(G2_CROSS00DEG+1)],
                   G2_CROSS01DEG, &fc01[k*(G2_CROSS01DEG+1)],
                   G2_CROSS02DEG, &fc02[k*(G2_CROSS02DEG+1)],
                   G2_CROSS10DEG, &fc10[k*(G2_CROSS10DEG+1)],
                   G2_CROSS11DEG, &fc11[k*(G2_CROSS11DEG+1)],
                   G2_CROSS12DEG, &fc12[k*(G2_CROSS12DEG+1)],
                   G2_CROSS00DEG, &fd00[k*(G2_CROSS00DEG+1)],
                   G2_CROSS01DEG, &fd01[k*(G2_CROSS01DEG+1)],
                   G2_CROSS02DEG, &fd02[k*(G2_CROSS02DEG+1)],
                   G2_CROSS10DEG, &fd10[k*(G2_CROSS10DEG+1)],
                   G2_CROSS11DEG, &fd11[k*(G2_CROSS11DEG+1)],
                   G2_CROSS12DEG, &fd12[k*(G2_CROSS12DEG+1)],
                   &degu, &degv, y );
    outpatch ( degu, degv, y, usrptr );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_OutputPatchesd*/

boolean g2h_FillHoled ( GHoleDomaind *domain,
                        int spdimen, CONST_ double *hole_cp,
                        double *acoeff, void *usrptr,
                        void (*outpatch) ( int n, int m, const double *cp,
                                           void *usrptr ) )
{
  void   *sp;
  G2HolePrivateRecd *privateG2;
  int    hole_k, nfunc_a, nfunc_b;
  double *lmat, *b, *x, *fc00;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  if ( !privateG2->Lmat )
    if ( !g2h_DecomposeMatrixd ( domain ) )
      goto failure;
    
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  lmat = privateG2->Lmat;
  x = pkv_GetScratchMemd ( spdimen*max(nfunc_a,nfunc_b) );
  b = pkv_GetScratchMemd ( spdimen*nfunc_a );
  fc00 = pkv_GetScratchMemd ( (G2_CROSSDEGSUM+6)*2*hole_k*spdimen );
  if ( !b || !x || !fc00 )
    goto failure;

  if ( !_g2h_SetRightSided ( domain, spdimen, hole_cp, fc00, b ) )
    goto failure;

  pkn_LowerTrMatrixSolved ( nfunc_a, lmat, spdimen, spdimen, b, spdimen, x );
  pkn_UpperTrMatrixSolved ( nfunc_a, lmat, spdimen, spdimen, x, spdimen, x );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nfunc_a*sizeof(double) );

  if ( !_g2h_OutputPatchesd ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_FillHoled*/

boolean g2h_FillHoleConstrd ( GHoleDomaind *domain,
                              int spdimen, CONST_ double *hole_cp,
                              int nconstr, CONST_ double *constr,
                              double *acoeff, void *usrptr,
                              void (*outpatch) ( int n, int m, const double *cp,
                                                 void *usrptr ) )
{
  void   *sp;
  G2HolePrivateRecd *privateG2;
  int    hole_k, nfunc_a, nfunc_b;
  double *lmat, *cmat, *rcmat, *b, *x, *y;
  double *fc00;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  if ( !privateG2->Lmat || !privateG2->Cmat || !privateG2->RCmat ||
       privateG2->nconstr != nconstr ) {
    domain->error_code = G2H_ERROR_UNDEFINED_CONSTR;
    goto failure;
  }
    
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  lmat = privateG2->Lmat;
  cmat = privateG2->Cmat;
  rcmat = privateG2->RCmat;
  b = pkv_GetScratchMemd ( spdimen*nfunc_a );
  x = pkv_GetScratchMemd ( spdimen*max(nfunc_a,nfunc_b) );
  fc00 = pkv_GetScratchMemd ( (G2_CROSSDEGSUM+6)*2*hole_k*spdimen );
  if ( !b || !x || !fc00 ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  if ( !_g2h_SetRightSided ( domain, spdimen, hole_cp, fc00, b ) )
    goto failure;

  pkn_LowerTrMatrixSolved ( nfunc_a, lmat, spdimen, spdimen, b, spdimen, x );
  pkn_UpperTrMatrixSolved ( nfunc_a, lmat, spdimen, spdimen, x, spdimen, x );
              /* now x is the solution without the constraints and it is time */
              /* to employ the constraints that have been previously set */
  pkn_MultMatrixd ( nconstr, nfunc_a, nfunc_a, cmat,
                    spdimen, spdimen, x, spdimen, b );
  pkn_AddMatrixd ( 1, spdimen*nconstr, 0, b, 0, constr, 0, b );
              /* solve the system with the Schur matrix */
  pkn_LowerTrMatrixSolved ( nconstr, rcmat, spdimen, spdimen, b, spdimen, b );
  pkn_UpperTrMatrixSolved ( nconstr, rcmat, spdimen, spdimen, b, spdimen, b );
              /* compute the constraints correction */
  y = pkv_GetScratchMemd ( spdimen*nfunc_a );
  if ( !y ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  pkn_MultTMatrixd ( nconstr, nfunc_a, nfunc_a, cmat,
                     spdimen, spdimen, b, spdimen, y );
  pkn_LowerTrMatrixSolved ( nfunc_a, lmat, spdimen, spdimen, y, spdimen, y );
  pkn_UpperTrMatrixSolved ( nfunc_a, lmat, spdimen, spdimen, y, spdimen, y );
  pkn_SubtractMatrixd ( 1, spdimen*nfunc_a, 0, x, 0, y, 0, x );
  pkv_FreeScratchMemd ( spdimen*nfunc_a );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nfunc_a*sizeof(double) );

  if ( !_g2h_OutputPatchesd ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_FillHoleConstrd*/

boolean g2h_FillHoleAltConstrd ( GHoleDomaind *domain,
                              int spdimen, CONST_ double *hole_cp,
                              int naconstr, CONST_ double *constr,
                              double *acoeff, void *usrptr,
                              void (*outpatch) ( int n, int m, const double *cp,
                                                 void *usrptr ) )
{
  void   *sp;
  G2HolePrivateRecd *privateG2;
  int    hole_k, nfunc_a, nfunc_b;
  double *lmat, *acmat, *arcmat, *b, *x, *y;
  double *fc00;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  if ( !privateG2->Lmat || !privateG2->ACmat || !privateG2->ARCmat ||
       privateG2->naconstr != naconstr || privateG2->acdim != spdimen ) {
    domain->error_code = G2H_ERROR_UNDEFINED_CONSTR;
    goto failure;
  }
    
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  lmat = privateG2->Lmat;
  acmat = privateG2->ACmat;
  arcmat = privateG2->ARCmat;
  b = pkv_GetScratchMemd ( spdimen*nfunc_a );
  x = pkv_GetScratchMemd ( spdimen*max(nfunc_a,nfunc_b) );
  y = pkv_GetScratchMemd ( spdimen*nfunc_a );
  fc00 = pkv_GetScratchMemd ( (G2_CROSSDEGSUM+6)*2*hole_k*spdimen );
  if ( !b || !x || !y || !fc00 ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  if ( !_g2h_SetRightSided ( domain, spdimen, hole_cp, fc00, b ) )
    goto failure;

  pkn_LowerTrMatrixSolved ( nfunc_a, lmat, spdimen, spdimen, b, spdimen, x );
  pkn_UpperTrMatrixSolved ( nfunc_a, lmat, spdimen, spdimen, x, spdimen, x );
              /* now x is the solution without the constraints and it is time */
              /* to employ the constraints that have been previously set */
  pkv_TransposeMatrixd ( nfunc_a, spdimen, spdimen, x, nfunc_a, y );
  pkn_MultMatrixd ( naconstr, spdimen*nfunc_a, spdimen*nfunc_a, acmat,
                    1, 1, y, 1, b );
  pkn_AddMatrixd ( 1, naconstr, 0, b, 0, constr, 0, b );
              /* solve the system with the Schur matrix */
  pkn_LowerTrMatrixSolved ( naconstr, arcmat, 1, 1, b, 1, b );
  pkn_UpperTrMatrixSolved ( naconstr, arcmat, 1, 1, b, 1, b );
              /* compute the constraints correction */
  pkn_MultTMatrixd ( naconstr, nfunc_a*spdimen, nfunc_a*spdimen, acmat,
                     1, 1, b, 1, y );
  pkv_TransposeMatrixd ( spdimen, nfunc_a, nfunc_a, y, spdimen, b );
  pkn_LowerTrMatrixSolved ( nfunc_a, lmat, spdimen, spdimen, b, spdimen, y );
  pkn_UpperTrMatrixSolved ( nfunc_a, lmat, spdimen, spdimen, y, spdimen, y );
  pkn_SubtractMatrixd ( 1, spdimen*nfunc_a, 0, x, 0, y, 0, x );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nfunc_a*sizeof(double) );

  if ( !_g2h_OutputPatchesd ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_FillHoleAltConstrd*/

