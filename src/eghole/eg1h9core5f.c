
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
               const float *r0, const float *r0s,
               float *pv )
{
  void  *sp;
  float *r0u, *aux1, *aux2;
  int   deg1, deg2;

  sp = pkv_GetScratchMemTop ();

  r0u = pkv_GetScratchMemf ( G1_CROSS00DEG*spdimen );
  aux1 = pkv_GetScratchMemf ( 2*(G1H_FINALDEG+1)*spdimen );
  if ( !r0u || !aux1 )
    goto failure;
  aux2 = &aux1[(G1H_FINALDEG+1)*spdimen];

  mbs_multiFindBezDerivativef ( G1_CROSS00DEG, 1, spdimen, 0, r0, 0, r0u );

        /* compute pv = b1*ru + c1*rs */
  if ( !mbs_multiMultBezCf ( 1, G1_BF01DEG, 0, b1, spdimen, 1, G1_CROSS00DEG-1, 0, r0u,
                             &deg1, 0, aux1 ) )
    goto failure;
  if ( !mbs_multiMultBezCf ( 1, G1_CG01DEG, 0, c1, spdimen, 1, G1_CROSS00DEG-1, 0, r0s,
                             &deg2, 0, aux2 ) )
    goto failure;

            /* at this point deg1 = 5, deg2 = 4 */
  mbs_multiBCDegElevf ( 1, spdimen, 0, deg2, aux2, deg1-deg2, 0, &deg2, pv );
  pkn_AddMatrixf ( 1, spdimen*(deg1+1), 0, pv, 0, aux1, 0, pv );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*FindDiCrossDeraf*/

static boolean FindDiCrossDerbf ( int spdimen,
               const float *b1, const float *c1,
               const float *r0, const float *r0s,
               float *pv )
{
  void  *sp;
  float *r0u;
  float *aux1;
  int   deg1, deg2;

  sp = pkv_GetScratchMemTop ();

  r0u = pkv_GetScratchMemf ( 8*spdimen );
  aux1 = pkv_GetScratchMemf ( 3*(G1H_FINALDEG+1)*spdimen );
  if ( !r0u || !aux1 )
    goto failure;

  mbs_multiFindBezDerivativef ( 3, 1, spdimen, 0, r0, 0, r0u );

        /* compute pv = b1*ru + c1*rs */
  if ( !mbs_multiMultBezCf ( 1, G1_BF11DEG, 0, b1, spdimen, 1, 2, 0, r0u,
                             &deg1, 0, aux1 ) )
    goto failure;
  if ( !mbs_multiMultBezCf ( 1, G1_CG11DEG, 0, c1, spdimen, 1, 3, 0, r0s,
                             &deg2, 0, pv ) )
    goto failure;

            /* at this point deg1 = deg2 = 5 = G1_CROSS11DEG */
  pkn_AddMatrixf ( 1, spdimen*(deg1+1), 0, aux1, 0, pv, 0, pv );

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
  G1HolePrivateRecf *privateG1;
  int      hole_k, i, j, n;
  float    *b01, *c01, *f01, *g01, *b11, *c11, *f11, *g11;
  vector2f *surrpc;
  vector2f *omc, *omcd;
  vector2f *dir0cr1, *diq0cr1, *dir1cr1, *diq1cr1;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG = domain->privateG;
  privateG1 = domain->privateG1;
  surrpc = privateG->surrpc;
  n = G1_CROSS01DEG+G1_CROSS11DEG;
  privateG1->dicross = malloc ( hole_k*2*(n+2)*sizeof(vector2f) );
  if ( !privateG1->dicross )
    goto failure;

  G1GetDiCrossAddresses ();

  G1GetPolyAddr ( privateG1->jfunc, b01, c01, f01, g01, b11, c11, f11, g11 );

  omc = privateG1->omc;
  omcd = &omc[(G1H_OMCDEG+1)*hole_k];

  for ( i = 0; i < hole_k; i++ ) {
    j = (i+1) % hole_k;
    FindDiCrossDeraf ( 2,
        &b01[i*(G1_BF01DEG+1)], &c01[i*(G1_CG01DEG+1)],
        (float*)&omc[i*(G1H_OMCDEG+1)], (float*)&omcd[i*G1H_OMCDEG],
        (float*)&dir0cr1[i*(G1_CROSS01DEG+1)] );
    FindDiCrossDeraf ( 2,
        &f01[i*(G1_BF01DEG+1)], &g01[i*(G1_CG01DEG+1)],
        (float*)&omc[j*(G1H_OMCDEG+1)], (float*)&omcd[j*G1H_OMCDEG],
        (float*)&diq0cr1[i*(G1_CROSS01DEG+1)] );
    FindDiCrossDerbf ( 2,
        &b11[i*(G1_BF11DEG+1)], &c11[i*(G1_CG11DEG+1)],
        (float*)&surrpc[24*j], (float*)&surrpc[24*j+4],
        (float*)&dir1cr1[i*(G1_CROSS11DEG+1)] );
    FindDiCrossDerbf ( 2,
        &f11[i*(G1_BF11DEG+1)], &g11[i*(G1_CG11DEG+1)],
        (float*)&surrpc[12*(2*i+1)], (float*)&surrpc[12*(2*i+1)+4],
        (float*)&diq1cr1[i*(G1_CROSS11DEG+1)] );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*FindDiPatchesf*/

