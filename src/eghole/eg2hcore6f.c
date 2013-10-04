
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

static boolean FindDiCrossDeraf ( int spdimen,
               const float *b1, const float *c1,
               const float *b2, const float *c2,
               const float *b1b1, const float *twob1c1, const float *c1c1,
               const float *r0, const float *r0s, const float *r0ss,
               float *pv, float *pvv )
{
  void  *sp;
  float *r0u, *r0uu, *r0us;
  float *aux1, *aux2, *aux3;
  int   deg1, deg2;

  sp = pkv_GetScratchMemTop ();

  r0u = pkv_GetScratchMemf ( (3*G2_CROSS00DEG-2)*spdimen );
  aux1 = pkv_GetScratchMemf ( 3*(G2H_FINALDEG+1)*spdimen );
  if ( !r0u || !aux1 )
    goto failure;
  r0uu = &r0u[G2_CROSS00DEG*spdimen];
  r0us = &r0uu[(G2_CROSS00DEG-1)*spdimen];
  aux2 = &aux1[(G2H_FINALDEG+1)*spdimen];
  aux3 = &aux2[(G2H_FINALDEG+1)*spdimen];

  mbs_multiFindBezDerivativef ( G2_CROSS00DEG, 1, spdimen, 0, r0, 0, r0u );
  mbs_multiFindBezDerivativef ( G2_CROSS00DEG-1, 1, spdimen, 0, r0u, 0, r0uu );
  mbs_multiFindBezDerivativef ( G2_CROSS00DEG-1, 1, spdimen, 0, r0s, 0, r0us );

        /* compute pv = b1*ru + c1*rs */
  if ( !mbs_multiMultBezCf ( 1, G2_BF01DEG, 0, b1, spdimen, 1, G2_CROSS00DEG-1, 0, r0u,
                       &deg1, 0, aux1 ) )
    goto failure;
  if ( !mbs_multiMultBezCf ( 1, G2_CG01DEG, 0, c1, spdimen, 1, G2_CROSS00DEG-1, 0, r0s,
                       &deg2, 0, aux2 ) )
    goto failure;

            /* at this point deg1 = 8, deg2 = 7 */
  mbs_multiBCDegElevf ( 1, spdimen, 0, deg2, aux2, deg1-deg2, 0, &deg2, pv );
  pkn_AddMatrixf ( 1, spdimen*(deg1+1), 0, pv, 0, aux1, 0, pv );

        /* compute pvv = b2*ru + c2*rs + b1^2ruu + 2*b1*c1*rus + c1^2*rss */
  if ( !mbs_multiMultBezCf ( 1, G2_BF02DEG, 0, b2, spdimen, 1, G2_CROSS00DEG-1, 0, r0u,
                       &deg1, 0, aux1 ) )
    goto failure;
  if ( !mbs_multiMultBezCf ( 1, G2_CG02DEG, 0, c2, spdimen, 1, G2_CROSS00DEG-1, 0, r0s,
                       &deg2, 0, aux2 ) )
    goto failure;

            /* at this point deg1 = 9, deg2 = 8 */
  mbs_multiBCDegElevf ( 1, spdimen, 0, deg2, aux2, deg1-deg2, 0, &deg2, pvv );
  pkn_AddMatrixf ( 1, spdimen*(deg1+1), 0, pvv, 0, aux1, 0, pvv );
  if ( !mbs_multiMultBezCf ( 1, 2*G2_BF01DEG, 0, b1b1, spdimen, 1, G2_CROSS00DEG-2, 0,
                       r0uu, &deg1, 0, aux1 ) )
    goto failure;
  if ( !mbs_multiMultBezCf ( 1, G2_BF01DEG+G2_CG01DEG, 0, twob1c1, spdimen, 1, G2_CROSS00DEG-2, 0,
                       r0us, &deg2, 0, aux2 ) )
    goto failure;
            /* at this point deg1 = 9, deg2 = 8 */
  mbs_multiBCDegElevf ( 1, spdimen, 0, deg2, aux2, deg1-deg2, 0, &deg2, aux3 );
  pkn_AddMatrixf ( 1, spdimen*(deg1+1), 0, aux1, 0, aux3, 0, aux1 );
  if ( !mbs_multiMultBezCf ( 1, 2*G2_CG01DEG, 0, c1c1, spdimen, 1, G2_CROSS00DEG-2, 0,
                       r0ss, &deg2, 0, aux2 ) )
    goto failure;
            /* at this point deg1 = 9, deg2 = 7 */
  mbs_multiBCDegElevf ( 1, spdimen, 0, deg2, aux2, deg1-deg2, 0, &deg2, aux3 );
  pkn_AddMatrixf ( 1, spdimen*(deg1+1), 0, aux1, 0, aux3, 0, aux1 );
  pkn_AddMatrixf ( 1, spdimen*(G2_CROSS02DEG+1), 0, pvv, 0, aux1, 0, pvv );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*FindDiCrossDeraf*/

static boolean FindDiCrossDerbf ( int spdimen,
               const float *b1, const float *c1,
               const float *b2, const float *c2,
               const float *b1b1, const float *twob1c1, const float *c1c1,
               const float *r0, const float *r0s, const float *r0ss,
               float *pv, float *pvv )
{
  void  *sp;
  float *r0u, *r0uu, *r0us;
  float *aux1, *aux2, *aux3;
  int   deg1, deg2;

  sp = pkv_GetScratchMemTop ();

  r0u = pkv_GetScratchMemf ( 8*spdimen );
  aux1 = pkv_GetScratchMemf ( 3*(G2H_FINALDEG+1)*spdimen );
  if ( !r0u || !aux1 )
    goto failure;
  r0uu = &r0u[3*spdimen];
  r0us = &r0uu[2*spdimen];
  aux2 = &aux1[(G2H_FINALDEG+1)*spdimen];
  aux3 = &aux2[(G2H_FINALDEG+1)*spdimen];

  mbs_multiFindBezDerivativef ( 3, 1, spdimen, 0, r0, 0, r0u );
  mbs_multiFindBezDerivativef ( 2, 1, spdimen, 0, r0u, 0, r0uu );
  mbs_multiFindBezDerivativef ( 3, 1, spdimen, 0, r0s, 0, r0us );

        /* compute pv = b1*ru + c1*rs */
  if ( !mbs_multiMultBezCf ( 1, G2_BF11DEG, 0, b1, spdimen, 1, 2, 0, r0u,
                       &deg1, 0, aux1 ) )
    goto failure;
  if ( !mbs_multiMultBezCf ( 1, G2_CG11DEG, 0, c1, spdimen, 1, 3, 0, r0s,
                       &deg2, 0, aux2 ) )
    goto failure;

            /* at this point deg1 = deg2 = 6 = G2_CROSS11DEG */
  pkn_AddMatrixf ( 1, spdimen*(deg1+1), 0, aux1, 0, aux2, 0, pv );

        /* compute pvv = b2*ru + c2*rs + b1^2ruu + 2*b1*c1*rus + c1^2*rss */
  if ( !mbs_multiMultBezCf ( 1, G2_BF12DEG, 0, b2, spdimen, 1, 2, 0, r0u,
                       &deg1, 0, aux1 ) )
    goto failure;
  if ( !mbs_multiMultBezCf ( 1, G2_CG12DEG, 0, c2, spdimen, 1, 3, 0, r0s,
                       &deg2, 0, aux2 ) )
    goto failure;
            /* at this point deg1 = deg2 = 7 = G2_CROSS12DEG-2 */
  pkn_AddMatrixf ( 1, spdimen*(deg1+1), 0, aux1, 0, aux2, 0, aux3 );
  mbs_multiBCDegElevf ( 1, spdimen, 0, deg1, aux3, G2_CROSS12DEG-deg1, 0, &deg1, pvv );

  if ( !mbs_multiMultBezCf ( 1, 2*G2_BF11DEG, 0, b1b1, spdimen, 1, 1, 0,
                       r0uu, &deg1, 0, aux1 ) )
    goto failure;
  if ( !mbs_multiMultBezCf ( 1, G2_BF11DEG+G2_CG11DEG, 0, twob1c1, spdimen, 1, 2, 0,
                       r0us, &deg2, 0, aux2 ) )
    goto failure;
            /* at this point deg1 = deg2 = 9 = G2_CROSS12DEG */
  pkn_AddMatrixf ( 1, spdimen*(deg1+1), 0, aux1, 0, aux2, 0, aux1 );
  if ( !mbs_multiMultBezCf ( 1, 2*G2_CG11DEG, 0, c1c1, spdimen, 1, 3, 0,
                       r0ss, &deg2, 0, aux2 ) )
    goto failure;
            /* at this point deg1 = deg2 = 9 = G2_CROSS12DEG */
  pkn_AddMatrixf ( 1, spdimen*(deg1+1), 0, aux1, 0, aux2, 0, aux1 );
  pkn_AddMatrixf ( 1, spdimen*(G2_CROSS12DEG+1), 0, pvv, 0, aux1, 0, pvv );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*FindDiCrossDerbf*/

static boolean FindDiPatchesf ( GHoleDomainf *domain )
{
  void     *sp;
  GHolePrivateRecf *privateG;
  G2HolePrivateRecf *privateG2;
  int      hole_k, i, j, n;
  float    *b01, *c01, *f01, *g01, *b11, *c11, *f11, *g11,
           *b02, *c02, *f02, *g02, *b12, *c12, *f12, *g12,
           *b01b01, *twob01c01, *c01c01, *f01f01, *twof01g01, *g01g01,
           *b11b11, *twob11c11, *c11c11, *f11f11, *twof11g11, *g11g11;
  vector2f *surrpc;
  vector2f *omc, *omcd, *omcdd;
  vector2f *dir0cr1, *diq0cr1, *dir1cr1, *diq1cr1,
           *dir0cr2, *diq0cr2, *dir1cr2, *diq1cr2;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG = domain->privateG;
  privateG2 = domain->privateG2;
  surrpc = privateG->surrpc;
  n = G2_CROSS01DEG+G2_CROSS02DEG+G2_CROSS11DEG+G2_CROSS12DEG;
  privateG2->dicross = malloc ( hole_k*2*(n+4)*sizeof(vector2f) );
  if ( !privateG2->dicross )
    goto failure;

  G2GetDiCrossAddresses ();

  G2GetPolynomialAddresses ( privateG2->jfunc,
      b01, c01, f01, g01, b11, c11, f11, g11,
      b02, c02, f02, g02, b12, c12, f12, g12, b01b01, twob01c01, c01c01,
      f01f01, twof01g01, g01g01, b11b11, twob11c11, c11c11,
      f11f11, twof11g11, g11g11 );

  omc = privateG2->omc;
  omcd = &omc[(G2H_OMCDEG+1)*hole_k];
  omcdd = &omcd[G2H_OMCDEG*hole_k];

  for ( i = 0; i < hole_k; i++ ) {
    j = (i+1) % hole_k;
    FindDiCrossDeraf ( 2,
        &b01[i*(G2_BF01DEG+1)], &c01[i*(G2_CG01DEG+1)],
        &b02[i*(G2_BF02DEG+1)], &c02[i*(G2_CG02DEG+1)],
        &b01b01[i*(2*G2_BF01DEG+1)], &twob01c01[i*(G2_BF01DEG+G2_CG01DEG+1)],
        &c01c01[i*(2*G2_CG01DEG+1)],
        (float*)&omc[i*(G2H_OMCDEG+1)], (float*)&omcd[i*G2H_OMCDEG],
        (float*)&omcdd[i*(G2H_OMCDEG-1)],
        (float*)&dir0cr1[i*(G2_CROSS01DEG+1)],
        (float*)&dir0cr2[i*(G2_CROSS02DEG+1)] );
    FindDiCrossDeraf ( 2,
        &f01[i*(G2_BF01DEG+1)], &g01[i*(G2_CG01DEG+1)],
        &f02[i*(G2_BF02DEG+1)], &g02[i*(G2_CG02DEG+1)],
        &f01f01[i*(2*G2_BF01DEG+1)], &twof01g01[i*(G2_BF01DEG+G2_CG01DEG+1)],
        &g01g01[i*(2*G2_CG01DEG+1)],
        (float*)&omc[j*(G2H_OMCDEG+1)], (float*)&omcd[j*G2H_OMCDEG],
        (float*)&omcdd[j*(G2H_OMCDEG-1)],
        (float*)&diq0cr1[i*(G2_CROSS01DEG+1)],
        (float*)&diq0cr2[i*(G2_CROSS02DEG+1)] );
    FindDiCrossDerbf ( 2,
        &b11[i*(G2_BF11DEG+1)], &c11[i*(G2_CG11DEG+1)],
        &b12[i*(G2_BF12DEG+1)], &c12[i*(G2_CG12DEG+1)],
        &b11b11[i*(2*G2_BF11DEG+1)], &twob11c11[i*(G2_BF11DEG+G2_CG11DEG+1)],
        &c11c11[i*(2*G2_CG11DEG+1)],
        (float*)&surrpc[24*j], (float*)&surrpc[24*j+4],
        (float*)&surrpc[24*j+8],
        (float*)&dir1cr1[i*(G2_CROSS11DEG+1)],
        (float*)&dir1cr2[i*(G2_CROSS12DEG+1)] );
    FindDiCrossDerbf ( 2,
        &f11[i*(G2_BF11DEG+1)], &g11[i*(G2_CG11DEG+1)],
        &f12[i*(G2_BF12DEG+1)], &g12[i*(G2_CG12DEG+1)],
        &f11f11[i*(2*G2_BF11DEG+1)], &twof11g11[i*(G2_BF11DEG+G2_CG11DEG+1)],
        &g11g11[i*(2*G2_CG11DEG+1)],
        (float*)&surrpc[12*(2*i+1)], (float*)&surrpc[12*(2*i+1)+4],
        (float*)&surrpc[12*(2*i+1)+8],
        (float*)&diq1cr1[i*(G2_CROSS11DEG+1)],
        (float*)&diq1cr2[i*(G2_CROSS12DEG+1)] );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*FindDiPatchesf*/

