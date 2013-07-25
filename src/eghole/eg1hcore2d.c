
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2008                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

void _g1h_GetDiPatchCurvesd ( GHoleDomaind *domain, int i,
                point2d **c00, vector2d **c01, point2d **c10, vector2d **c11,
                point2d **d00, vector2d **d01, point2d **d10, vector2d **d11 )
{
  int       hole_k, j;
  GHolePrivateRecd *privateG;
  G1HolePrivateRecd *privateG1;
  vector2d *omc, *surrpc;
  vector2d *dir0cr1, *diq0cr1, *dir1cr1, *diq1cr1;
            
  hole_k = domain->hole_k;
  privateG = domain->privateG;
  privateG1 = domain->privateG1;
  omc = privateG1->omc;
  surrpc = privateG->surrpc;
  G1GetDiCrossAddresses ();

  j = (i+1) % hole_k;
  *c00 = &omc[i*(G1H_OMCDEG+1)];
  *c01 = &dir0cr1[i*(G1_CROSS01DEG+1)];
  *c10 = &surrpc[24*j];
  *c11 = &dir1cr1[i*(G1_CROSS11DEG+1)];
  *d00 = &omc[j*(G1H_OMCDEG+1)];
  *d01 = &diq0cr1[i*(G1_CROSS01DEG+1)];
  *d10 = &surrpc[12*(2*i+1)];
  *d11 = &diq1cr1[i*(G1_CROSS11DEG+1)];
} /*_g1h_GetDiPatchCurvesd*/

void _g1h_GetBFAPatchCurvesd ( GHoleDomaind *domain, int fn, int i,
                    double **c00, double **c01, double **d00, double **d01 )
{
  int    hole_k, nfunc_a, j;
  G1HolePrivateRecd *privateG1;
  double *bbr0, *bbr0cr1, *bbq0cr1;

  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;

  G1GetBFuncACrossAddresses();

  j = (i+1) % hole_k;
  *c00 = &bbr0[(fn*hole_k+i)*(G1H_OMCDEG+1)];
  *c01 = &bbr0cr1[(fn*hole_k+i)*(G1_CROSS01DEG+1)];
  *d00 = &bbr0[(fn*hole_k+j)*(G1H_OMCDEG+1)];
  *d01 = &bbq0cr1[(fn*hole_k+i)*(G1_CROSS01DEG+1)];
} /*_g1h_GetBFAPatchCurvesd*/

void _g1h_GetBFBPatchCurvesd ( GHoleDomaind *domain, int fn, int i,
                    double **c00, double **c01, double **c10, double **c11,
                    double **d00, double **d01, double **d10, double **d11 )
{
  int    hole_k, nfunc_b, j;
  G1HolePrivateRecd *privateG1;
  double *bbr0, *bbr1, *bbq1, *bbr0cr1, *bbr1cr1, *bbq0cr1, *bbq1cr1;

  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  nfunc_b = privateG1->nfunc_b;

  G1GetBFuncBCrossAddresses();

  j = (i+1) % hole_k;
  *c00 = &bbr0[(fn*hole_k+i)*(G1H_OMCDEG+1)];
  *c01 = &bbr0cr1[(fn*hole_k+i)*(G1_CROSS01DEG+1)];
  *c10 = &bbr1[(fn*hole_k+i)*(G1_CROSS10DEG+1)];
  *c11 = &bbr1cr1[(fn*hole_k+i)*(G1_CROSS11DEG+1)];
  *d00 = &bbr0[(fn*hole_k+j)*(G1H_OMCDEG+1)];
  *d01 = &bbq0cr1[(fn*hole_k+i)*(G1_CROSS01DEG+1)];
  *d10 = &bbq1[(fn*hole_k+i)*(G1_CROSS10DEG+1)];
  *d11 = &bbq1cr1[(fn*hole_k+i)*(G1_CROSS11DEG+1)];
} /*_g1h_GetBFBPatchCurvesd*/

/* ////////////////////////////////////////////////////////////////////////// */
#ifdef NICO
static void _g1h_GetDomSurrndBFuncd ( GHoleDomaind *domain,
                                      int i, int j, const double *bfc, double *bf )
{
  int     hole_k, k;
  int     *ind;
  double  *q;
  void   *sp;
  double *ukn, *vkn;

  sp  = pkv_GetScratchMemTop ();
  ind = pkv_GetScratchMem ( 16*sizeof(int) );
  q   = pkv_GetScratchMemd ( 16 );
  if ( !ind || !q )
    exit ( 1 );

  hole_k = domain->hole_k;
  gh_GetBspInd ( hole_k, i, j, ind );
  for ( k = 0; k < 16; k++ )
    q[k] = bfc[ind[k]];

  ukn = &domain->hole_knots[11*((i+hole_k-1) % hole_k)+3];
  vkn = &domain->hole_knots[11*i+j];
  mbs_BSPatchToBezd ( 1, 3, 7, ukn, 3, 7, vkn, 4, q,
                      NULL, NULL, NULL, NULL, NULL, NULL, 4, bf );

  pkv_SetScratchMemTop ( sp );
} /*_g1h_GetDomSurrndBFuncd*/
#endif
