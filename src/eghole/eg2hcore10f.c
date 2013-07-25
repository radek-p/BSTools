
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2010                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

boolean g2h_DecomposeMatrixf ( GHoleDomainf *domain )
{
  void   *sp;
  G2HolePrivateRecf *privateG2;
  float  *lmat;
  int    nfunc_a;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  if ( !privateG2->Lmat ) {
    if ( !privateG2->Amat )
      if ( !g2h_ComputeFormMatrixf ( domain ) )
        goto failure;
    nfunc_a = privateG2->nfunc_a;
    lmat = privateG2->Lmat = malloc ( (nfunc_a*(nfunc_a+1)/2)*sizeof(float) );
    if ( !lmat )
      goto failure;
    memcpy ( lmat, privateG2->Amat, (nfunc_a*(nfunc_a+1)/2)*sizeof(float) );
    if ( !pkn_CholeskyDecompf ( nfunc_a, lmat ) ) {
      domain->error_code = G2H_ERROR_NONPOSITIVE_MATRIX;
      goto failure;
    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_DecomposeMatrixf*/

boolean _g2h_SetRightSidef ( GHoleDomainf *domain,
                             int spdimen, CONST_ float *hole_cp,
                             float *fc00, float *b )
{
  void   *sp;
  GHolePrivateRecf  *privateG;
  G2HolePrivateRecf *privateG2;
  unsigned char *bfcpn;
  float  *bc00, *bc01, *bc02, *bc10, *bc11, *bc12,
         *bd00, *bd01, *bd02, *bd10, *bd11, *bd12;
  float         *fc01, *fc02, *fc10, *fc11, *fc12,
         *fd00, *fd01, *fd02, *fd10, *fd11, *fd12;
  float  *x, *y, *cp;
  int    hole_k, nfunc_a, nfunc_b;
  int    i, j, k, l;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG  = domain->privateG;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  bfcpn   = privateG->bfcpn;

  x = pkv_GetScratchMemf ( spdimen*nfunc_b );
  y = pkv_GetScratchMemf ( spdimen*(G2H_FINALDEG+1)*(G2H_FINALDEG+1) );
  if ( !x || !y )
    goto failure;

  memset ( fc00, 0, (G2_CROSSDEGSUM+6)*2*hole_k*spdimen*sizeof(float) );
  G2GetFCAddresses ();

  for ( j = 0; j < nfunc_b; j++ ) {
    l = bfcpn[j];
    cp = &hole_cp[l*spdimen];
    memcpy ( &x[j*spdimen], cp, spdimen*sizeof(float) );
    for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
      _g2h_GetBFBPatchCurvesf ( domain, j, i,
                                &bc00, &bc01, &bc02, &bc10, &bc11, &bc12,
                                &bd00, &bd01, &bd02, &bd10, &bd11, &bd12 );

      pkn_MultMatrixf ( G2_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS00DEG+1), 0, &fc00[k*(G2_CROSS00DEG+1)],
                       0, y, 0, &fc00[k*(G2_CROSS00DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS01DEG+1, 1, 1, bc01, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS01DEG+1), 0, &fc01[k*(G2_CROSS01DEG+1)],
                       0, y, 0, &fc01[k*(G2_CROSS01DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS02DEG+1, 1, 1, bc02, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS02DEG+1), 0, &fc02[k*(G2_CROSS02DEG+1)],
                       0, y, 0, &fc02[k*(G2_CROSS02DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS10DEG+1, 1, 1, bc10, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS10DEG+1), 0, &fc10[k*(G2_CROSS10DEG+1)],
                       0, y, 0, &fc10[k*(G2_CROSS10DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS11DEG+1, 1, 1, bc11, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS11DEG+1), 0, &fc11[k*(G2_CROSS11DEG+1)],
                       0, y, 0, &fc11[k*(G2_CROSS11DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS12DEG+1, 1, 1, bc12, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS12DEG+1), 0, &fc12[k*(G2_CROSS12DEG+1)],
                       0, y, 0, &fc12[k*(G2_CROSS12DEG+1)] );

      pkn_MultMatrixf ( G2_CROSS00DEG+1, 1, 1, bd00, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS00DEG+1), 0, &fd00[k*(G2_CROSS00DEG+1)],
                       0, y, 0, &fd00[k*(G2_CROSS00DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS01DEG+1, 1, 1, bd01, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS01DEG+1), 0, &fd01[k*(G2_CROSS01DEG+1)],
                       0, y, 0, &fd01[k*(G2_CROSS01DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS02DEG+1, 1, 1, bd02, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS02DEG+1), 0, &fd02[k*(G2_CROSS02DEG+1)],
                       0, y, 0, &fd02[k*(G2_CROSS02DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS10DEG+1, 1, 1, bd10, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS10DEG+1), 0, &fd10[k*(G2_CROSS10DEG+1)],
                       0, y, 0, &fd10[k*(G2_CROSS10DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS11DEG+1, 1, 1, bd11, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS11DEG+1), 0, &fd11[k*(G2_CROSS11DEG+1)],
                       0, y, 0, &fd11[k*(G2_CROSS11DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS12DEG+1, 1, 1, bd12, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS12DEG+1), 0, &fd12[k*(G2_CROSS12DEG+1)],
                       0, y, 0, &fd12[k*(G2_CROSS12DEG+1)] );
    }
  }

  if ( b )
    pkn_MultMatrixf ( nfunc_a, nfunc_b, nfunc_b, privateG2->Bmat,
                      spdimen, spdimen, x, spdimen, b );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_SetRightSidef*/

boolean _g2h_OutputPatchesf ( GHoleDomainf *domain,
                int spdimen, CONST_ float *x, float *fc00, void *usrptr,
                void (*outpatch) ( int n, int m, const float *cp,
                                   void *usrptr ) )
{
  void   *sp;
  G2HolePrivateRecf *privateG2;
  float  *bc00, *bc01, *bc02, *bd00, *bd01, *bd02;
  float         *fc01, *fc02, *fc10, *fc11, *fc12,
         *fd00, *fd01, *fd02, *fd10, *fd11, *fd12;
  float  *y, *cp;
  int    hole_k, nfunc_a, degu, degv;
  int    i, j, k;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;

  y = pkv_GetScratchMemf ( spdimen*(G2H_FINALDEG+1)*(G2H_FINALDEG+1) );
  if ( !y )
    goto failure;

  G2GetFCAddresses ();

  for ( j = 0; j < nfunc_a; j++ ) {
    cp = &x[j*spdimen];
    for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
      _g2h_GetBFAPatchCurvesf ( domain, j, i,
                                &bc00, &bc01, &bc02, &bd00, &bd01, &bd02 );

      pkn_MultMatrixf ( G2_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixf ( 1, spdimen*(G2_CROSS00DEG+1), 0, &fc00[k*(G2_CROSS00DEG+1)],
                            0, y, 0, &fc00[k*(G2_CROSS00DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS01DEG+1, 1, 1, bc01, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixf ( 1, spdimen*(G2_CROSS01DEG+1), 0, &fc01[k*(G2_CROSS01DEG+1)],
                            0, y, 0, &fc01[k*(G2_CROSS01DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS02DEG+1, 1, 1, bc02, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixf ( 1, spdimen*(G2_CROSS02DEG+1), 0, &fc02[k*(G2_CROSS02DEG+1)],
                            0, y, 0, &fc02[k*(G2_CROSS02DEG+1)] );

      pkn_MultMatrixf ( G2_CROSS00DEG+1, 1, 1, bd00, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixf ( 1, spdimen*(G2_CROSS00DEG+1), 0, &fd00[k*(G2_CROSS00DEG+1)],
                            0, y, 0, &fd00[k*(G2_CROSS00DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS01DEG+1, 1, 1, bd01, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixf ( 1, spdimen*(G2_CROSS01DEG+1), 0, &fd01[k*(G2_CROSS01DEG+1)],
                            0, y, 0, &fd01[k*(G2_CROSS01DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS02DEG+1, 1, 1, bd02, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixf ( 1, spdimen*(G2_CROSS02DEG+1), 0, &fd02[k*(G2_CROSS02DEG+1)],
                            0, y, 0, &fd02[k*(G2_CROSS02DEG+1)] );
    }
  }

  for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
    mbs_BezC2CoonsToBezf ( spdimen,
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
} /*_g2h_OutputPatchesf*/

boolean g2h_FillHolef ( GHoleDomainf *domain,
                        int spdimen, CONST_ float *hole_cp,
                        float *acoeff, void *usrptr,
                        void (*outpatch) ( int n, int m, const float *cp,
                                           void *usrptr ) )
{
  void   *sp;
  G2HolePrivateRecf *privateG2;
  int    hole_k, nfunc_a, nfunc_b;
  float  *lmat, *b, *x, *fc00;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  if ( !privateG2->Lmat )
    if ( !g2h_DecomposeMatrixf ( domain ) )
      goto failure;
    
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  lmat = privateG2->Lmat;
  x = pkv_GetScratchMemf ( spdimen*max(nfunc_a,nfunc_b) );
  b = pkv_GetScratchMemf ( spdimen*nfunc_a );
  fc00 = pkv_GetScratchMemf ( (G2_CROSSDEGSUM+6)*2*hole_k*spdimen );
  if ( !b || !x || !fc00 )
    goto failure;

  if ( !_g2h_SetRightSidef ( domain, spdimen, hole_cp, fc00, b ) )
    goto failure;

  pkn_LowerTrMatrixSolvef ( nfunc_a, lmat, spdimen, spdimen, b, spdimen, x );
  pkn_UpperTrMatrixSolvef ( nfunc_a, lmat, spdimen, spdimen, x, spdimen, x );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nfunc_a*sizeof(float) );

  if ( !_g2h_OutputPatchesf ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_FillHolef*/

boolean g2h_FillHoleConstrf ( GHoleDomainf *domain,
                              int spdimen, CONST_ float *hole_cp,
                              int nconstr, CONST_ float *constr,
                              float *acoeff, void *usrptr,
                              void (*outpatch) ( int n, int m, const float *cp,
                                                 void *usrptr ) )
{
  void   *sp;
  G2HolePrivateRecf *privateG2;
  int    hole_k, nfunc_a, nfunc_b;
  float  *lmat, *cmat, *rcmat, *b, *x, *y;
  float  *fc00;

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
  b = pkv_GetScratchMemf ( spdimen*nfunc_a );
  x = pkv_GetScratchMemf ( spdimen*max(nfunc_a,nfunc_b) );
  fc00 = pkv_GetScratchMemf ( (G2_CROSSDEGSUM+6)*2*hole_k*spdimen );
  if ( !b || !x || !fc00 ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  if ( !_g2h_SetRightSidef ( domain, spdimen, hole_cp, fc00, b ) )
    goto failure;

  pkn_LowerTrMatrixSolvef ( nfunc_a, lmat, spdimen, spdimen, b, spdimen, x );
  pkn_UpperTrMatrixSolvef ( nfunc_a, lmat, spdimen, spdimen, x, spdimen, x );
              /* now x is the solution without the constraints and it is time */
              /* to employ the constraints that have been previously set */
  pkn_MultMatrixf ( nconstr, nfunc_a, nfunc_a, cmat,
                    spdimen, spdimen, x, spdimen, b );
  pkn_AddMatrixf ( 1, spdimen*nconstr, 0, b, 0, constr, 0, b );
              /* solve the system with the Schur matrix */
  pkn_LowerTrMatrixSolvef ( nconstr, rcmat, spdimen, spdimen, b, spdimen, b );
  pkn_UpperTrMatrixSolvef ( nconstr, rcmat, spdimen, spdimen, b, spdimen, b );
              /* compute the constraints correction */
  y = pkv_GetScratchMemf ( spdimen*nfunc_a );
  if ( !y ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  pkn_MultTMatrixf ( nconstr, nfunc_a, nfunc_a, cmat,
                     spdimen, spdimen, b, spdimen, y );
  pkn_LowerTrMatrixSolvef ( nfunc_a, lmat, spdimen, spdimen, y, spdimen, y );
  pkn_UpperTrMatrixSolvef ( nfunc_a, lmat, spdimen, spdimen, y, spdimen, y );
  pkn_SubtractMatrixf ( 1, spdimen*nfunc_a, 0, x, 0, y, 0, x );
  pkv_FreeScratchMemf ( spdimen*nfunc_a );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nfunc_a*sizeof(float) );

  if ( !_g2h_OutputPatchesf ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_FillHoleConstrf*/

boolean g2h_FillHoleAltConstrf ( GHoleDomainf *domain,
                              int spdimen, CONST_ float *hole_cp,
                              int naconstr, CONST_ float *constr,
                              float *acoeff, void *usrptr,
                              void (*outpatch) ( int n, int m, const float *cp,
                                                 void *usrptr ) )
{
  void   *sp;
  G2HolePrivateRecf *privateG2;
  int    hole_k, nfunc_a, nfunc_b;
  float  *lmat, *acmat, *arcmat, *b, *x, *y;
  float  *fc00;

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
  b = pkv_GetScratchMemf ( spdimen*nfunc_a );
  x = pkv_GetScratchMemf ( spdimen*max(nfunc_a,nfunc_b) );
  y = pkv_GetScratchMemf ( spdimen*nfunc_a );
  fc00 = pkv_GetScratchMemf ( (G2_CROSSDEGSUM+6)*2*hole_k*spdimen );
  if ( !b || !x || !y || !fc00 ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  if ( !_g2h_SetRightSidef ( domain, spdimen, hole_cp, fc00, b ) )
    goto failure;

  pkn_LowerTrMatrixSolvef ( nfunc_a, lmat, spdimen, spdimen, b, spdimen, x );
  pkn_UpperTrMatrixSolvef ( nfunc_a, lmat, spdimen, spdimen, x, spdimen, x );
              /* now x is the solution without the constraints and it is time */
              /* to employ the constraints that have been previously set */
  pkv_TransposeMatrixf ( nfunc_a, spdimen, spdimen, x, nfunc_a, y );
  pkn_MultMatrixf ( naconstr, spdimen*nfunc_a, spdimen*nfunc_a, acmat,
                    1, 1, y, 1, b );
  pkn_AddMatrixf ( 1, naconstr, 0, b, 0, constr, 0, b );
              /* solve the system with the Schur matrix */
  pkn_LowerTrMatrixSolvef ( naconstr, arcmat, 1, 1, b, 1, b );
  pkn_UpperTrMatrixSolvef ( naconstr, arcmat, 1, 1, b, 1, b );
              /* compute the constraints correction */
  pkn_MultTMatrixf ( naconstr, nfunc_a*spdimen, nfunc_a*spdimen, acmat,
                     1, 1, b, 1, y );
  pkv_TransposeMatrixf ( spdimen, nfunc_a, nfunc_a, y, spdimen, b );
  pkn_LowerTrMatrixSolvef ( nfunc_a, lmat, spdimen, spdimen, b, spdimen, y );
  pkn_UpperTrMatrixSolvef ( nfunc_a, lmat, spdimen, spdimen, y, spdimen, y );
  pkn_SubtractMatrixf ( 1, spdimen*nfunc_a, 0, x, 0, y, 0, x );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nfunc_a*sizeof(float) );

  if ( !_g2h_OutputPatchesf ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_FillHoleAltConstrf*/

