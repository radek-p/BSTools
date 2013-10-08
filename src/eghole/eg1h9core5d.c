
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

static boolean FindDiCrossDerad ( int spdimen,
               const double *b1, const double *c1,
               const double *r0, const double *r0s,
               double *pv )
{
  void  *sp;
  double *r0u, *aux1, *aux2;
  int   deg1, deg2;

  sp = pkv_GetScratchMemTop ();

  r0u = pkv_GetScratchMemd ( G1_CROSS00DEG*spdimen );
  aux1 = pkv_GetScratchMemd ( 2*(G1H_FINALDEG+1)*spdimen );
  if ( !r0u || !aux1 )
    goto failure;
  aux2 = &aux1[(G1H_FINALDEG+1)*spdimen];

  mbs_multiFindBezDerivatived ( G1_CROSS00DEG, 1, spdimen, 0, r0, 0, r0u );

        /* compute pv = b1*ru + c1*rs */
  if ( !mbs_multiMultBezCd ( 1, G1_BF01DEG, 0, b1, spdimen, 1, G1_CROSS00DEG-1, 0, r0u,
                             &deg1, 0, aux1 ) )
    goto failure;
  if ( !mbs_multiMultBezCd ( 1, G1_CG01DEG, 0, c1, spdimen, 1, G1_CROSS00DEG-1, 0, r0s,
                             &deg2, 0, aux2 ) )
    goto failure;

            /* at this point deg1 = 5, deg2 = 4 */
  if ( !mbs_multiBCDegElevd ( 1, spdimen, 0, deg2, aux2, deg1-deg2, 0, &deg2, pv ) )
    goto failure;
  pkn_AddMatrixd ( 1, spdimen*(deg1+1), 0, pv, 0, aux1, 0, pv );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*FindDiCrossDerad*/

static boolean FindDiCrossDerbd ( int spdimen,
               const double *b1, const double *c1,
               const double *r0, const double *r0s,
               double *pv )
{
  void  *sp;
  double *r0u;
  double *aux1;
  int   deg1, deg2;

  sp = pkv_GetScratchMemTop ();

  r0u = pkv_GetScratchMemd ( 8*spdimen );
  aux1 = pkv_GetScratchMemd ( 3*(G1H_FINALDEG+1)*spdimen );
  if ( !r0u || !aux1 )
    goto failure;

  mbs_multiFindBezDerivatived ( 3, 1, spdimen, 0, r0, 0, r0u );

        /* compute pv = b1*ru + c1*rs */
  if ( !mbs_multiMultBezCd ( 1, G1_BF11DEG, 0, b1, spdimen, 1, 2, 0, r0u,
                             &deg1, 0, aux1 ) )
    goto failure;
  if ( !mbs_multiMultBezCd ( 1, G1_CG11DEG, 0, c1, spdimen, 1, 3, 0, r0s,
                             &deg2, 0, pv ) )
    goto failure;

            /* at this point deg1 = deg2 = 5 = G1_CROSS11DEG */
  pkn_AddMatrixd ( 1, spdimen*(deg1+1), 0, aux1, 0, pv, 0, pv );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*FindDiCrossDerbd*/

static boolean FindDiPatchesd ( GHoleDomaind *domain )
{
  void     *sp;
  GHolePrivateRecd *privateG;
  G1HolePrivateRecd *privateG1;
  int      hole_k, i, j, n;
  double    *b01, *c01, *f01, *g01, *b11, *c11, *f11, *g11;
  vector2d *surrpc;
  vector2d *omc, *omcd;
  vector2d *dir0cr1, *diq0cr1, *dir1cr1, *diq1cr1;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG = domain->privateG;
  privateG1 = domain->privateG1;
  surrpc = privateG->surrpc;
  n = G1_CROSS01DEG+G1_CROSS11DEG;
  privateG1->dicross = malloc ( hole_k*2*(n+2)*sizeof(vector2d) );
  if ( !privateG1->dicross )
    goto failure;

  G1GetDiCrossAddresses ();

  G1GetPolyAddr ( privateG1->jfunc, b01, c01, f01, g01, b11, c11, f11, g11 );

  omc = privateG1->omc;
  omcd = &omc[(G1H_OMCDEG+1)*hole_k];

  for ( i = 0; i < hole_k; i++ ) {
    j = (i+1) % hole_k;
    FindDiCrossDerad ( 2,
        &b01[i*(G1_BF01DEG+1)], &c01[i*(G1_CG01DEG+1)],
        (double*)&omc[i*(G1H_OMCDEG+1)], (double*)&omcd[i*G1H_OMCDEG],
        (double*)&dir0cr1[i*(G1_CROSS01DEG+1)] );
    FindDiCrossDerad ( 2,
        &f01[i*(G1_BF01DEG+1)], &g01[i*(G1_CG01DEG+1)],
        (double*)&omc[j*(G1H_OMCDEG+1)], (double*)&omcd[j*G1H_OMCDEG],
        (double*)&diq0cr1[i*(G1_CROSS01DEG+1)] );
    FindDiCrossDerbd ( 2,
        &b11[i*(G1_BF11DEG+1)], &c11[i*(G1_CG11DEG+1)],
        (double*)&surrpc[24*j], (double*)&surrpc[24*j+4],
        (double*)&dir1cr1[i*(G1_CROSS11DEG+1)] );
    FindDiCrossDerbd ( 2,
        &f11[i*(G1_BF11DEG+1)], &g11[i*(G1_CG11DEG+1)],
        (double*)&surrpc[12*(2*i+1)], (double*)&surrpc[12*(2*i+1)+4],
        (double*)&diq1cr1[i*(G1_CROSS11DEG+1)] );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*FindDiPatchesd*/

