
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2010                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

boolean g1h_DecomposeMatrixd ( GHoleDomaind *domain )
{
  void   *sp;
  G1HolePrivateRecd *privateG1;
  double *lmat;
  int    nfunc_a;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  if ( !privateG1->Lmat ) {
    if ( !privateG1->Amat )
      if ( !g1h_ComputeFormMatrixd ( domain ) )
        goto failure;
    nfunc_a = privateG1->nfunc_a;
    lmat = privateG1->Lmat = malloc ( (nfunc_a*(nfunc_a+1)/2)*sizeof(double) );
    if ( !lmat )
      goto failure;
    memcpy ( lmat, privateG1->Amat, (nfunc_a*(nfunc_a+1)/2)*sizeof(double) );
    if ( !pkn_CholeskyDecompd ( nfunc_a, lmat ) ) {
      domain->error_code = G1H_ERROR_NONPOSITIVE_MATRIX;
      goto failure;
    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_DecomposeMatrixd*/

boolean _g1h_SetRightSided ( GHoleDomaind *domain, const double *Bmat,
                             int spdimen, CONST_ double *hole_cp,
                             double *fc00, double *b )
{
  void   *sp;
  GHolePrivateRecd  *privateG;
  G1HolePrivateRecd *privateG1;
  unsigned char *bfcpn;
  double *bc00, *bc01, *bc10, *bc11, *bd00, *bd01, *bd10, *bd11;
  double        *fc01, *fc10, *fc11, *fd00, *fd01, *fd10, *fd11;
  double *x, *y, *cp;
  int    hole_k, nfunc_a, nfunc_b;
  int    i, j, k, l;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG  = domain->privateG;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  bfcpn   = privateG->bfcpn;

  x = pkv_GetScratchMemd ( spdimen*nfunc_b );
  y = pkv_GetScratchMemd ( spdimen*(G1H_FINALDEG+1)*(G1H_FINALDEG+1) );
  if ( !x || !y )
    goto failure;

  memset ( fc00, 0, (G1_CROSSDEGSUM+4)*2*hole_k*spdimen*sizeof(double) );
  G1GetFCAddresses ();

  for ( j = 0; j < nfunc_b; j++ ) {
    l = bfcpn[j];
    cp = &hole_cp[l*spdimen];
    memcpy ( &x[j*spdimen], cp, spdimen*sizeof(double) );
    for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
      _g1h_GetBFBPatchCurvesd ( domain, j, i, &bc00, &bc01, &bc10, &bc11,
                                &bd00, &bd01, &bd10, &bd11 );

      pkn_MultMatrixd ( G1_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G1_CROSS00DEG+1), 0, &fc00[k*(G1_CROSS00DEG+1)],
                       0, y, 0, &fc00[k*(G1_CROSS00DEG+1)] );
      pkn_MultMatrixd ( G1_CROSS01DEG+1, 1, 1, bc01, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G1_CROSS01DEG+1), 0, &fc01[k*(G1_CROSS01DEG+1)],
                       0, y, 0, &fc01[k*(G1_CROSS01DEG+1)] );
      pkn_MultMatrixd ( G1_CROSS10DEG+1, 1, 1, bc10, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G1_CROSS10DEG+1), 0, &fc10[k*(G1_CROSS10DEG+1)],
                       0, y, 0, &fc10[k*(G1_CROSS10DEG+1)] );
      pkn_MultMatrixd ( G1_CROSS11DEG+1, 1, 1, bc11, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G1_CROSS11DEG+1), 0, &fc11[k*(G1_CROSS11DEG+1)],
                       0, y, 0, &fc11[k*(G1_CROSS11DEG+1)] );

      pkn_MultMatrixd ( G1_CROSS00DEG+1, 1, 1, bd00, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G1_CROSS00DEG+1), 0, &fd00[k*(G1_CROSS00DEG+1)],
                       0, y, 0, &fd00[k*(G1_CROSS00DEG+1)] );
      pkn_MultMatrixd ( G1_CROSS01DEG+1, 1, 1, bd01, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G1_CROSS01DEG+1), 0, &fd01[k*(G1_CROSS01DEG+1)],
                       0, y, 0, &fd01[k*(G1_CROSS01DEG+1)] );
      pkn_MultMatrixd ( G1_CROSS10DEG+1, 1, 1, bd10, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G1_CROSS10DEG+1), 0, &fd10[k*(G1_CROSS10DEG+1)],
                       0, y, 0, &fd10[k*(G1_CROSS10DEG+1)] );
      pkn_MultMatrixd ( G1_CROSS11DEG+1, 1, 1, bd11, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G1_CROSS11DEG+1), 0, &fd11[k*(G1_CROSS11DEG+1)],
                       0, y, 0, &fd11[k*(G1_CROSS11DEG+1)] );
    }
  }

  if ( b )
    pkn_MultMatrixd ( nfunc_a, nfunc_b, nfunc_b, Bmat,
                      spdimen, spdimen, x, spdimen, b );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_SetRightSided*/

boolean _g1h_OutputPatchesd ( GHoleDomaind *domain, int spdimen,
                       CONST_ double *x, double *fc00, void *usrptr,
                       void (*outpatch) ( int n, int m, const double *cp,
                                          void *usrptr ) )
{
  void   *sp;
  G1HolePrivateRecd *privateG1;
  double *bc00, *bc01, *bd00, *bd01;
  double *fc01, *fc10, *fc11, *fd00, *fd01, *fd10, *fd11;
  double *y, *cp;
  int    hole_k, nfunc_a, degu, degv;
  int    i, j, k;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;

  y = pkv_GetScratchMemd ( spdimen*(G1H_FINALDEG+1)*(G1H_FINALDEG+1) );
  if ( !y )
    goto failure;

  G1GetFCAddresses ();

  for ( j = 0; j < nfunc_a; j++ ) {
    cp = &x[j*spdimen];
    for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
      _g1h_GetBFAPatchCurvesd ( domain, j, i, &bc00, &bc01, &bd00, &bd01 );

      pkn_MultMatrixd ( G1_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixd ( 1, spdimen*(G1_CROSS00DEG+1), 0, &fc00[k*(G1_CROSS00DEG+1)],
                            0, y, 0, &fc00[k*(G1_CROSS00DEG+1)] );
      pkn_MultMatrixd ( G1_CROSS01DEG+1, 1, 1, bc01, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixd ( 1, spdimen*(G1_CROSS01DEG+1), 0, &fc01[k*(G1_CROSS01DEG+1)],
                            0, y, 0, &fc01[k*(G1_CROSS01DEG+1)] );

      pkn_MultMatrixd ( G1_CROSS00DEG+1, 1, 1, bd00, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixd ( 1, spdimen*(G1_CROSS00DEG+1), 0, &fd00[k*(G1_CROSS00DEG+1)],
                            0, y, 0, &fd00[k*(G1_CROSS00DEG+1)] );
      pkn_MultMatrixd ( G1_CROSS01DEG+1, 1, 1, bd01, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixd ( 1, spdimen*(G1_CROSS01DEG+1), 0, &fd01[k*(G1_CROSS01DEG+1)],
                            0, y, 0, &fd01[k*(G1_CROSS01DEG+1)] );
    }
  }

  for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
    mbs_BezC1CoonsToBezd ( spdimen,
                   G1_CROSS00DEG, &fc00[k*(G1_CROSS00DEG+1)],
                   G1_CROSS01DEG, &fc01[k*(G1_CROSS01DEG+1)],
                   G1_CROSS10DEG, &fc10[k*(G1_CROSS10DEG+1)],
                   G1_CROSS11DEG, &fc11[k*(G1_CROSS11DEG+1)],
                   G1_CROSS00DEG, &fd00[k*(G1_CROSS00DEG+1)],
                   G1_CROSS01DEG, &fd01[k*(G1_CROSS01DEG+1)],
                   G1_CROSS10DEG, &fd10[k*(G1_CROSS10DEG+1)],
                   G1_CROSS11DEG, &fd11[k*(G1_CROSS11DEG+1)],
                   &degu, &degv, y );
    outpatch ( degu, degv, y, usrptr );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_OutputPatchesd*/

boolean g1h_FillHoled ( GHoleDomaind *domain,
                        int spdimen, CONST_ double *hole_cp,
                        double *acoeff, void *usrptr,
                        void (*outpatch) ( int n, int m, const double *cp,
                                           void *usrptr ) )
{
  void   *sp;
  G1HolePrivateRecd *privateG1;
  int    hole_k, nfunc_a, nfunc_b;
  double *lmat, *b, *x, *fc00;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  if ( !privateG1->Lmat )
    if ( !g1h_DecomposeMatrixd ( domain ) )
      goto failure;
    
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  lmat = privateG1->Lmat;
  x = pkv_GetScratchMemd ( spdimen*max(nfunc_a,nfunc_b) );
  b = pkv_GetScratchMemd ( spdimen*nfunc_a );
  fc00 = pkv_GetScratchMemd ( (G1_CROSSDEGSUM+4)*2*hole_k*spdimen );
  if ( !b || !x || !fc00 )
    goto failure;

  if ( !_g1h_SetRightSided ( domain, privateG1->Bmat,
                             spdimen, hole_cp, fc00, b ) )
    goto failure;

  pkn_LowerTrMatrixSolved ( nfunc_a, lmat, spdimen, spdimen, b, spdimen, x );
  pkn_UpperTrMatrixSolved ( nfunc_a, lmat, spdimen, spdimen, x, spdimen, x );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nfunc_a*sizeof(double) );

  if ( !_g1h_OutputPatchesd ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_FillHoled*/

boolean g1h_FillHoleConstrd ( GHoleDomaind *domain,
                              int spdimen, CONST_ double *hole_cp,
                              int nconstr, CONST_ double *constr,
                              double *acoeff, void *usrptr,
                              void (*outpatch) ( int n, int m, const double *cp,
                                                 void *usrptr ) )
{
  void   *sp;
  G1HolePrivateRecd *privateG1;
  int    hole_k, nfunc_a, nfunc_b;
  double *lmat, *cmat, *rcmat, *b, *x, *y;
  double *fc00;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  if ( !privateG1->Lmat || !privateG1->Cmat || !privateG1->RCmat ||
       privateG1->nconstr != nconstr ) {
    domain->error_code = G1H_ERROR_UNDEFINED_CONSTR;
    goto failure;
  }
    
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  lmat = privateG1->Lmat;
  cmat = privateG1->Cmat;
  rcmat = privateG1->RCmat;
  b = pkv_GetScratchMemd ( spdimen*nfunc_a );
  x = pkv_GetScratchMemd ( spdimen*max(nfunc_a,nfunc_b) );
  fc00 = pkv_GetScratchMemd ( (G1_CROSSDEGSUM+4)*2*hole_k*spdimen );
  if ( !b || !x || !fc00 ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  if ( !_g1h_SetRightSided ( domain, privateG1->Bmat,
                             spdimen, hole_cp, fc00, b ) )
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
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
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

  if ( !_g1h_OutputPatchesd ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_FillHoleConstrd*/

boolean g1h_FillHoleAltConstrd ( GHoleDomaind *domain,
                              int spdimen, CONST_ double *hole_cp,
                              int naconstr, CONST_ double *constr,
                              double *acoeff, void *usrptr,
                              void (*outpatch) ( int n, int m, const double *cp,
                                                 void *usrptr ) )
{
  void   *sp;
  G1HolePrivateRecd *privateG1;
  int    hole_k, nfunc_a, nfunc_b;
  double *lmat, *acmat, *arcmat, *b, *x, *y;
  double *fc00;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  if ( !privateG1->Lmat || !privateG1->ACmat || !privateG1->ARCmat ||
       privateG1->naconstr != naconstr || privateG1->acdim != spdimen ) {
    domain->error_code = G1H_ERROR_UNDEFINED_CONSTR;
    goto failure;
  }
    
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  lmat = privateG1->Lmat;
  acmat = privateG1->ACmat;
  arcmat = privateG1->ARCmat;
  b = pkv_GetScratchMemd ( spdimen*nfunc_a );
  x = pkv_GetScratchMemd ( spdimen*max(nfunc_a,nfunc_b) );
  y = pkv_GetScratchMemd ( spdimen*nfunc_a );
  fc00 = pkv_GetScratchMemd ( (G1_CROSSDEGSUM+4)*2*hole_k*spdimen );
  if ( !b || !x || !y || !fc00 ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  if ( !_g1h_SetRightSided ( domain, privateG1->Bmat,
                             spdimen, hole_cp, fc00, b ) )
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

  if ( !_g1h_OutputPatchesd ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_FillHoleAltConstrd*/

