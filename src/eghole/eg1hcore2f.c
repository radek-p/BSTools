
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2008                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

void _g1h_GetDiPatchCurvesf ( GHoleDomainf *domain, int i,
                point2f **c00, vector2f **c01, point2f **c10, vector2f **c11,
                point2f **d00, vector2f **d01, point2f **d10, vector2f **d11 )
{
  int       hole_k, j;
  GHolePrivateRecf *privateG;
  G1HolePrivateRecf *privateG1;
  vector2f *omc, *surrpc;
  vector2f *dir0cr1, *diq0cr1, *dir1cr1, *diq1cr1;
            
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
} /*_g1h_GetDiPatchCurvesf*/

void _g1h_GetBFAPatchCurvesf ( GHoleDomainf *domain, int fn, int i,
                    float **c00, float **c01, float **d00, float **d01 )
{
  int   hole_k, nfunc_a, j;
  G1HolePrivateRecf *privateG1;
  float *bbr0, *bbr0cr1, *bbq0cr1;

  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;

  G1GetBFuncACrossAddresses();

  j = (i+1) % hole_k;
  *c00 = &bbr0[(fn*hole_k+i)*(G1H_OMCDEG+1)];
  *c01 = &bbr0cr1[(fn*hole_k+i)*(G1_CROSS01DEG+1)];
  *d00 = &bbr0[(fn*hole_k+j)*(G1H_OMCDEG+1)];
  *d01 = &bbq0cr1[(fn*hole_k+i)*(G1_CROSS01DEG+1)];
} /*_g1h_GetBFAPatchCurvesf*/

void _g1h_GetBFBPatchCurvesf ( GHoleDomainf *domain, int fn, int i,
                    float **c00, float **c01, float **c10, float **c11,
                    float **d00, float **d01, float **d10, float **d11 )
{
  int   hole_k, nfunc_b, j;
  G1HolePrivateRecf *privateG1;
  float *bbr0, *bbr1, *bbq1, *bbr0cr1, *bbr1cr1, *bbq0cr1, *bbq1cr1;

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
} /*_g1h_GetBFBPatchCurvesf*/

