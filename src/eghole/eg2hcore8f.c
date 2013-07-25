
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2010                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

boolean _g2h_TabDiPatchJac3f ( int nkn, const float *kn, const float *hfunc,
        const float *dhfunc, const float *ddhfunc, const float *dddhfunc,
        const vector2f *c00, const vector2f *c01, const vector2f *c02,
        const vector2f *c10, const vector2f *c11, const vector2f *c12,
        const vector2f *d00, const vector2f *d01, const vector2f *d02,
        const vector2f *d10, const vector2f *d11, const vector2f *d12,
        float *jac, float *trd )
{
  void     *sp;
  int      i, k;
  vector2f *tabpu, *tabpv, *tabpuu, *tabpuv, *tabpvv,
           *tabpuuu, *tabpuuv, *tabpuvv, *tabpvvv;

  sp = pkv_GetScratchMemTop ();
  k = nkn*nkn;
  tabpu = (vector2f*)pkv_GetScratchMem ( 9*k*sizeof(vector2f) );
  if ( !tabpu )
    goto failure;
  tabpv = &tabpu[k];      tabpuu = &tabpv[k];    tabpuv = &tabpuu[k];
  tabpvv = &tabpuv[k];    tabpuuu = &tabpvv[k];  tabpuuv = &tabpuuu[k];
  tabpuvv = &tabpuuv[k];  tabpvvv = &tabpuvv[k];

  if ( !mbs_TabBezC2CoonsDer3f ( 2, nkn, kn, hfunc, dhfunc, ddhfunc, dddhfunc,
           nkn, kn, hfunc, dhfunc, ddhfunc, dddhfunc,
           G2_CROSS00DEG, (float*)c00, G2_CROSS01DEG, (float*)c01,
           G2_CROSS02DEG, (float*)c02,
           G2_CROSS10DEG, (float*)c10, G2_CROSS11DEG, (float*)c11,
           G2_CROSS12DEG, (float*)c12,
           G2_CROSS00DEG, (float*)d00, G2_CROSS01DEG, (float*)d01,
           G2_CROSS02DEG, (float*)d02,
           G2_CROSS10DEG, (float*)d10, G2_CROSS11DEG, (float*)d11,
           G2_CROSS12DEG, (float*)d12,
           NULL, (float*)tabpu, (float*)tabpv, (float*)tabpuu,
           (float*)tabpuv, (float*)tabpvv, (float*)tabpuuu,
           (float*)tabpuuv, (float*)tabpuvv, (float*)tabpvvv ) )
    goto failure;

  for ( i = 0; i < k; i++ )
    _g2h_DiJacobian3f ( &tabpu[i], &tabpv[i], &tabpuu[i], &tabpuv[i], &tabpvv[i],
                        &tabpuuu[i], &tabpuuv[i], &tabpuvv[i], &tabpvvv[i],
                        &jac[i], &trd[18*i] );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_TabDiPatchJac3f*/

boolean _g2h_TabLaplacianGradf ( int nkn, const float *tkn,
        const float *hfunc, const float *dhfunc, const float *ddhfunc,
        const float *dddhfunc,
        const float *fc00, const float *fc01, const float *fc02,
        const float *fc10, const float *fc11, const float *fc12,
        const float *fd00, const float *fd01, const float *fd02,
        const float *fd10, const float *fd11, const float *fd12,
        const float *trd,
        vector2f *lapgrad )
{
  void  *sp;
  int   i, k;
  float *tabfu, *tabfv, *tabfuu, *tabfuv, *tabfvv,
        *tabfuuu, *tabfuuv, *tabfuvv, *tabfvvv;

  sp = pkv_GetScratchMemTop ();
  k = nkn*nkn;

  tabfu = pkv_GetScratchMemf ( 9*k );
  if ( !tabfu )
    goto failure;
  tabfv = &tabfu[k];  tabfuu = &tabfv[k];  tabfuv = &tabfuu[k];
  tabfvv = &tabfuv[k];  tabfuuu = &tabfvv[k];  tabfuuv = &tabfuuu[k];
  tabfuvv = &tabfuuv[k];  tabfvvv = &tabfuvv[k];

  if ( !mbs_TabBezC2CoonsDer3f ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
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
} /*_g2h_TabLaplacianGradf*/

boolean _g2h_TabLaplacianGrad0f ( int nkn, const float *tkn,
        const float *hfunc, const float *dhfunc, const float *ddhfunc,
        const float *dddhfunc,
        const float *fc00, const float *fc01, const float *fc02,
        const float *fd00, const float *fd01, const float *fd02,
        const float *trd,
        vector2f *lapgrad )
{
  void  *sp;
  int   i, k;
  float *tabfu, *tabfv, *tabfuu, *tabfuv, *tabfvv,
        *tabfuuu, *tabfuuv, *tabfuvv, *tabfvvv;

  sp = pkv_GetScratchMemTop ();
  k = nkn*nkn;

  tabfu = pkv_GetScratchMemf ( 9*k );
  if ( !tabfu )
    goto failure;
  tabfv = &tabfu[k];  tabfuu = &tabfv[k];  tabfuv = &tabfuu[k];
  tabfvv = &tabfuv[k];  tabfuuu = &tabfvv[k];  tabfuuv = &tabfuuu[k];
  tabfuvv = &tabfuuv[k];  tabfvvv = &tabfuvv[k];

  if ( !mbs_TabBezC2Coons0Der3f ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
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
} /*_g2h_TabLaplacianGrad0f*/

boolean g2h_ComputeFormMatrixf ( GHoleDomainf *domain )
{
  void     *sp;
  GHolePrivateRecf  *privateG;
  G2HolePrivateRecf *privateG2;
  int      hole_k, nfunc_a, nfunc_b;
  int      i, j, l, n;
  float    *tkn, *hfunc, *dhfunc, *ddhfunc, *dddhfunc;
  vector2f *c00, *c01, *c02, *c10, *c11, *c12,
           *d00, *d01, *d02, *d10, *d11, *d12;
  float    *fc00, *fc01, *fc02, *fc10, *fc11, *fc12,
           *fd00, *fd01, *fd02, *fd10, *fd11, *fd12;
  float    *trd, *jac, *Amat, *Bmat;
  vector2f *lgr;
  unsigned short *support_b;

  sp = pkv_GetScratchMemTop ();
  if ( !domain->basisG2 ) {
    if ( !g2h_ComputeBasisf ( domain ) )
      goto failure;
  }
  hole_k  = domain->hole_k;
  privateG  = domain->privateG;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  support_b = privateG->support_b;
  Amat = privateG2->Amat = malloc ( (nfunc_a*(nfunc_a+1)/2)*sizeof(float) );
  Bmat = privateG2->Bmat = malloc ( nfunc_a*nfunc_b*sizeof(float) );
  if ( !Amat || !Bmat )
    goto failure;

  n = hole_k*G2_NQUADSQ;

  tkn = pkv_GetScratchMemf ( 25*G2_NQUAD );
  jac = pkv_GetScratchMemf ( n );
  trd = pkv_GetScratchMemf ( 18*n );
  lgr = (vector2f*)pkv_GetScratchMem ( (nfunc_a+1)*n*sizeof(vector2f) );
  if ( !tkn || !jac || !trd || !lgr )
    goto failure;
  hfunc = &tkn[G2_NQUAD];         dhfunc = &hfunc[6*G2_NQUAD];
  ddhfunc = &dhfunc[6*G2_NQUAD];  dddhfunc = &ddhfunc[6*G2_NQUAD];
  _gh_PrepareTabKnotsf ( G2_NQUAD, privateG2->opt_quad, tkn );
  mbs_TabQuinticHFuncDer3f ( 0.0, 1.0, G2_NQUAD, tkn,
                             hfunc, dhfunc, ddhfunc, dddhfunc );

  for ( i = 0; i < hole_k; i++ ) {
    _g2h_GetDiPatchCurvesf ( domain, i,
                             &c00, &c01, &c02, &c10, &c11, &c12,
                             &d00, &d01, &d02, &d10, &d11, &d12 );
    _g2h_TabDiPatchJac3f ( G2_NQUAD, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                           c00, c01, c02, c10, c11, c12,
                           d00, d01, d02, d10, d11, d12,
                           &jac[i*G2_NQUADSQ], &trd[i*G2_NQUADSQ*18] );
  }

  for ( j = 0; j < nfunc_a; j++ )
    for ( i = 0; i < hole_k; i++ ) {
      _g2h_GetBFAPatchCurvesf ( domain, j, i,
                                &fc00, &fc01, &fc02, &fd00, &fd01, &fd02 );
      _g2h_TabLaplacianGrad0f ( G2_NQUAD, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                                fc00, fc01, fc02, fd00, fd01, fd02,
                                &trd[i*G2_NQUADSQ*18],
                                &lgr[(j*hole_k+i)*G2_NQUADSQ] );
    }

  for ( i = l = 0;  i < nfunc_a;  i++ ) {
    for ( j = 0;  j <= i;  j++, l++ )
      Amat[l] = _g2h_Integralf ( hole_k, G2_NQUADSQ, jac,
                                 0xFFFF, &lgr[i*n], 0xFFFF, &lgr[j*n] );
  }
  
  for ( j = 0; j < nfunc_b; j++ ) {
    for ( i = 0; i < hole_k; i++ )
      if ( support_b[j] & (0x0001 << i) ) {
        _g2h_GetBFBPatchCurvesf ( domain, j, i,
                                  &fc00, &fc01, &fc02, &fc10, &fc11, &fc12,
                                  &fd00, &fd01, &fd02, &fd10, &fd11, &fd12 );
        _g2h_TabLaplacianGradf ( G2_NQUAD, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                                 fc00, fc01, fc02, fc10, fc11, fc12,
                                 fd00, fd01, fd02, fd10, fd11, fd12,
                                 &trd[i*G2_NQUADSQ*18],
                                 &lgr[(nfunc_a*hole_k+i)*G2_NQUADSQ] );
      }
    for ( i = 0; i < nfunc_a; i++ )
      Bmat[i*nfunc_b+j] = _g2h_Integralf ( hole_k, G2_NQUADSQ, jac,
                            0xFFFF, &lgr[i*n], support_b[j], &lgr[nfunc_a*n] );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_ComputeFormMatrixf*/

