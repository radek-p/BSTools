
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

boolean _g2h_TabDiPatchJac3d ( int nkn, const double *kn, const double *hfunc,
        const double *dhfunc, const double *ddhfunc, const double *dddhfunc,
        const vector2d *c00, const vector2d *c01, const vector2d *c02,
        const vector2d *c10, const vector2d *c11, const vector2d *c12,
        const vector2d *d00, const vector2d *d01, const vector2d *d02,
        const vector2d *d10, const vector2d *d11, const vector2d *d12,
        double *jac, double *trd )
{
  void     *sp;
  int      i, k;
  vector2d *tabpu, *tabpv, *tabpuu, *tabpuv, *tabpvv,
           *tabpuuu, *tabpuuv, *tabpuvv, *tabpvvv;

  sp = pkv_GetScratchMemTop ();
  k = nkn*nkn;
  tabpu = (vector2d*)pkv_GetScratchMem ( 9*k*sizeof(vector2d) );
  if ( !tabpu )
    goto failure;
  tabpv = &tabpu[k];      tabpuu = &tabpv[k];    tabpuv = &tabpuu[k];
  tabpvv = &tabpuv[k];    tabpuuu = &tabpvv[k];  tabpuuv = &tabpuuu[k];
  tabpuvv = &tabpuuv[k];  tabpvvv = &tabpuvv[k];

  if ( !mbs_TabBezC2CoonsDer3d ( 2, nkn, kn, hfunc, dhfunc, ddhfunc, dddhfunc,
           nkn, kn, hfunc, dhfunc, ddhfunc, dddhfunc,
           G2_CROSS00DEG, (double*)c00, G2_CROSS01DEG, (double*)c01,
           G2_CROSS02DEG, (double*)c02,
           G2_CROSS10DEG, (double*)c10, G2_CROSS11DEG, (double*)c11,
           G2_CROSS12DEG, (double*)c12,
           G2_CROSS00DEG, (double*)d00, G2_CROSS01DEG, (double*)d01,
           G2_CROSS02DEG, (double*)d02,
           G2_CROSS10DEG, (double*)d10, G2_CROSS11DEG, (double*)d11,
           G2_CROSS12DEG, (double*)d12,
           NULL, (double*)tabpu, (double*)tabpv, (double*)tabpuu,
           (double*)tabpuv, (double*)tabpvv, (double*)tabpuuu,
           (double*)tabpuuv, (double*)tabpuvv, (double*)tabpvvv ) )
    goto failure;

  for ( i = 0; i < k; i++ )
    if ( !_g2h_DiJacobian3d ( &tabpu[i], &tabpv[i], &tabpuu[i], &tabpuv[i], &tabpvv[i],
                              &tabpuuu[i], &tabpuuv[i], &tabpuvv[i], &tabpvvv[i],
                              &jac[i], &trd[18*i] ) )
      goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_TabDiPatchJac3d*/

boolean _g2h_TabLaplacianGradd ( int nkn, const double *tkn,
        const double *hfunc, const double *dhfunc, const double *ddhfunc,
        const double *dddhfunc,
        const double *fc00, const double *fc01, const double *fc02,
        const double *fc10, const double *fc11, const double *fc12,
        const double *fd00, const double *fd01, const double *fd02,
        const double *fd10, const double *fd11, const double *fd12,
        const double *trd,
        vector2d *lapgrad )
{
  void   *sp;
  int    i, k;
  double *tabfu, *tabfv, *tabfuu, *tabfuv, *tabfvv,
         *tabfuuu, *tabfuuv, *tabfuvv, *tabfvvv;

  sp = pkv_GetScratchMemTop ();
  k = nkn*nkn;

  tabfu = pkv_GetScratchMemd ( 9*k );
  if ( !tabfu )
    goto failure;
  tabfv = &tabfu[k];  tabfuu = &tabfv[k];  tabfuv = &tabfuu[k];
  tabfvv = &tabfuv[k];  tabfuuu = &tabfvv[k];  tabfuuv = &tabfuuu[k];
  tabfuvv = &tabfuuv[k];  tabfvvv = &tabfuvv[k];

  if ( !mbs_TabBezC2CoonsDer3d ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                     nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                     G2_CROSS00DEG, fc00, G2_CROSS01DEG, fc01, G2_CROSS02DEG, fc02,
                     G2_CROSS10DEG, fc10, G2_CROSS11DEG, fc11, G2_CROSS12DEG, fc12,
                     G2_CROSS00DEG, fd00, G2_CROSS01DEG, fd01, G2_CROSS02DEG, fd02,
                     G2_CROSS10DEG, fd10, G2_CROSS11DEG, fd11, G2_CROSS12DEG, fd12,
                     NULL, tabfu, tabfv, tabfuu, tabfuv, tabfvv,
                     tabfuuu, tabfuuv, tabfuvv, tabfvvv ) )
    goto failure;

  for ( i = k = 0;  i < nkn*nkn;  i++, k += 18 ) {
    lapgrad[i].x = trd[k+0]*tabfu[i] + trd[k+1]*tabfv[i] +
                   trd[k+2]*tabfuu[i] + trd[k+3]*tabfuv[i] +
                   trd[k+4]*tabfvv[i] +
                   trd[k+5]*tabfuuu[i] + trd[k+6]*tabfuuv[i] +
                   trd[k+7]*tabfuvv[i] + trd[k+8]*tabfvvv[i];
    lapgrad[i].y = trd[k+9]*tabfu[i] + trd[k+10]*tabfv[i] +
                   trd[k+11]*tabfuu[i] + trd[k+12]*tabfuv[i] +
                   trd[k+13]*tabfvv[i] +
                   trd[k+14]*tabfuuu[i] + trd[k+15]*tabfuuv[i] +
                   trd[k+16]*tabfuvv[i] + trd[k+17]*tabfvvv[i];
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_TabLaplacianGradd*/

boolean _g2h_TabLaplacianGrad0d ( int nkn, const double *tkn,
        const double *hfunc, const double *dhfunc, const double *ddhfunc,
        const double *dddhfunc,
        const double *fc00, const double *fc01, const double *fc02,
        const double *fd00, const double *fd01, const double *fd02,
        const double *trd,
        vector2d *lapgrad )
{
  void   *sp;
  int    i, k;
  double *tabfu, *tabfv, *tabfuu, *tabfuv, *tabfvv,
         *tabfuuu, *tabfuuv, *tabfuvv, *tabfvvv;

  sp = pkv_GetScratchMemTop ();
  k = nkn*nkn;

  tabfu = pkv_GetScratchMemd ( 9*k );
  if ( !tabfu )
    goto failure;
  tabfv = &tabfu[k];  tabfuu = &tabfv[k];  tabfuv = &tabfuu[k];
  tabfvv = &tabfuv[k];  tabfuuu = &tabfvv[k];  tabfuuv = &tabfuuu[k];
  tabfuvv = &tabfuuv[k];  tabfvvv = &tabfuvv[k];

  if ( !mbs_TabBezC2Coons0Der3d ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                     nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                     G2_CROSS00DEG, fc00, G2_CROSS01DEG, fc01, G2_CROSS02DEG, fc02,
                     G2_CROSS00DEG, fd00, G2_CROSS01DEG, fd01, G2_CROSS02DEG, fd02,
                     NULL, tabfu, tabfv, tabfuu, tabfuv, tabfvv,
                     tabfuuu, tabfuuv, tabfuvv, tabfvvv ) )
    goto failure;

  for ( i = k = 0;  i < nkn*nkn;  i++, k += 18 ) {
    lapgrad[i].x = trd[k+0]*tabfu[i] + trd[k+1]*tabfv[i] +
                   trd[k+2]*tabfuu[i] + trd[k+3]*tabfuv[i] +
                   trd[k+4]*tabfvv[i] +
                   trd[k+5]*tabfuuu[i] + trd[k+6]*tabfuuv[i] +
                   trd[k+7]*tabfuvv[i] + trd[k+8]*tabfvvv[i];
    lapgrad[i].y = trd[k+9]*tabfu[i] + trd[k+10]*tabfv[i] +
                   trd[k+11]*tabfuu[i] + trd[k+12]*tabfuv[i] +
                   trd[k+13]*tabfvv[i] +
                   trd[k+14]*tabfuuu[i] + trd[k+15]*tabfuuv[i] +
                   trd[k+16]*tabfuvv[i] + trd[k+17]*tabfvvv[i];
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_TabLaplacianGrad0d*/

boolean g2h_ComputeFormMatrixd ( GHoleDomaind *domain )
{
  void     *sp;
  GHolePrivateRecd  *privateG;
  G2HolePrivateRecd *privateG2;
  int      hole_k, nfunc_a, nfunc_b;
  int      i, j, l, n;
  double   *tkn, *hfunc, *dhfunc, *ddhfunc, *dddhfunc;
  vector2d *c00, *c01, *c02, *c10, *c11, *c12,
           *d00, *d01, *d02, *d10, *d11, *d12;
  double   *fc00, *fc01, *fc02, *fc10, *fc11, *fc12,
           *fd00, *fd01, *fd02, *fd10, *fd11, *fd12;
  double   *trd, *jac, *Amat, *Bmat;
  vector2d *lgr;
  unsigned short *support_b;

  sp = pkv_GetScratchMemTop ();
  if ( !domain->basisG2 ) {
    if ( !g2h_ComputeBasisd ( domain ) )
      goto failure;
  }
  hole_k  = domain->hole_k;
  privateG  = domain->privateG;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  support_b = privateG->support_b;
  Amat = privateG2->Amat = malloc ( (nfunc_a*(nfunc_a+1)/2)*sizeof(double) );
  Bmat = privateG2->Bmat = malloc ( nfunc_a*nfunc_b*sizeof(double) );
  if ( !Amat || !Bmat )
    goto failure;

  n = hole_k*G2_NQUADSQ;

  tkn = pkv_GetScratchMemd ( 25*G2_NQUAD );
  jac = pkv_GetScratchMemd ( n );
  trd = pkv_GetScratchMemd ( 18*n );
  lgr = (vector2d*)pkv_GetScratchMem ( (nfunc_a+1)*n*sizeof(vector2d) );
  if ( !tkn || !jac || !trd || !lgr )
    goto failure;
  hfunc = &tkn[G2_NQUAD];         dhfunc = &hfunc[6*G2_NQUAD];
  ddhfunc = &dhfunc[6*G2_NQUAD];  dddhfunc = &ddhfunc[6*G2_NQUAD];
  _gh_PrepareTabKnotsd ( G2_NQUAD, privateG2->opt_quad, tkn );
  if ( !mbs_TabQuinticHFuncDer3d ( 0.0, 1.0, G2_NQUAD, tkn,
                                   hfunc, dhfunc, ddhfunc, dddhfunc ) )
    goto failure;

  for ( i = 0; i < hole_k; i++ ) {
    _g2h_GetDiPatchCurvesd ( domain, i,
                             &c00, &c01, &c02, &c10, &c11, &c12,
                             &d00, &d01, &d02, &d10, &d11, &d12 );
    _g2h_TabDiPatchJac3d ( G2_NQUAD, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                           c00, c01, c02, c10, c11, c12,
                           d00, d01, d02, d10, d11, d12,
                           &jac[i*G2_NQUADSQ], &trd[i*G2_NQUADSQ*18] );
  }

  for ( j = 0; j < nfunc_a; j++ )
    for ( i = 0; i < hole_k; i++ ) {
      _g2h_GetBFAPatchCurvesd ( domain, j, i,
                                &fc00, &fc01, &fc02, &fd00, &fd01, &fd02 );
      _g2h_TabLaplacianGrad0d ( G2_NQUAD, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                                fc00, fc01, fc02, fd00, fd01, fd02,
                                &trd[i*G2_NQUADSQ*18],
                                &lgr[(j*hole_k+i)*G2_NQUADSQ] );
    }

  for ( i = l = 0;  i < nfunc_a;  i++ ) {
    for ( j = 0;  j <= i;  j++, l++ )
      Amat[l] = _g2h_Integrald ( hole_k, G2_NQUADSQ, jac,
                                 0xFFFF, &lgr[i*n], 0xFFFF, &lgr[j*n] );
  }
  
  for ( j = 0; j < nfunc_b; j++ ) {
    for ( i = 0; i < hole_k; i++ )
      if ( support_b[j] & (0x0001 << i) ) {
        _g2h_GetBFBPatchCurvesd ( domain, j, i,
                                  &fc00, &fc01, &fc02, &fc10, &fc11, &fc12,
                                  &fd00, &fd01, &fd02, &fd10, &fd11, &fd12 );
        _g2h_TabLaplacianGradd ( G2_NQUAD, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                                 fc00, fc01, fc02, fc10, fc11, fc12,
                                 fd00, fd01, fd02, fd10, fd11, fd12,
                                 &trd[i*G2_NQUADSQ*18],
                                 &lgr[(nfunc_a*hole_k+i)*G2_NQUADSQ] );
      }
    for ( i = 0; i < nfunc_a; i++ )
      Bmat[i*nfunc_b+j] = _g2h_Integrald ( hole_k, G2_NQUADSQ, jac,
                            0xFFFF, &lgr[i*n], support_b[j], &lgr[nfunc_a*n] );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_ComputeFormMatrixd*/

