
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2008                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

void _g2h_GetDiPatchCurvesf ( GHoleDomainf *domain, int i,
                    point2f **c00, vector2f **c01, vector2f **c02,
                    point2f **c10, vector2f **c11, vector2f **c12,
                    point2f **d00, vector2f **d01, vector2f **d02,
                    point2f **d10, vector2f **d11, vector2f **d12 )
{
  int       hole_k, j;
  GHolePrivateRecf *privateG;
  G2HolePrivateRecf *privateG2;
  vector2f *omc, *surrpc;
  vector2f *dir0cr1, *diq0cr1, *dir1cr1, *diq1cr1,
           *dir0cr2, *diq0cr2, *dir1cr2, *diq1cr2;
            
  hole_k = domain->hole_k;
  privateG = domain->privateG;
  privateG2 = domain->privateG2;
  omc = privateG2->omc;
  surrpc = privateG->surrpc;
  G2GetDiCrossAddresses ();

  j = (i+1) % hole_k;
  *c00 = &omc[i*(G2H_OMCDEG+1)];
  *c01 = &dir0cr1[i*(G2_CROSS01DEG+1)];
  *c02 = &dir0cr2[i*(G2_CROSS02DEG+1)];
  *c10 = &surrpc[24*j];
  *c11 = &dir1cr1[i*(G2_CROSS11DEG+1)];
  *c12 = &dir1cr2[i*(G2_CROSS12DEG+1)];
  *d00 = &omc[j*(G2H_OMCDEG+1)];
  *d01 = &diq0cr1[i*(G2_CROSS01DEG+1)];
  *d02 = &diq0cr2[i*(G2_CROSS02DEG+1)];
  *d10 = &surrpc[12*(2*i+1)];
  *d11 = &diq1cr1[i*(G2_CROSS11DEG+1)];
  *d12 = &diq1cr2[i*(G2_CROSS12DEG+1)];
} /*_g2h_GetDiPatchCurvesf*/

void _g2h_GetBFAPatchCurvesf ( GHoleDomainf *domain, int fn, int i,
                    float **c00, float **c01, float **c02,
                    float **d00, float **d01, float **d02 )
{
  int   hole_k, nfunc_a, j;
  G2HolePrivateRecf *privateG2;
  float *bbr0, *bbr0cr1, *bbq0cr1, *bbr0cr2, *bbq0cr2;

  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;

  G2GetBFuncACrossAddresses();

  j = (i+1) % hole_k;
  *c00 = &bbr0[(fn*hole_k+i)*(G2_CROSS00DEG+1)];
  *c01 = &bbr0cr1[(fn*hole_k+i)*(G2_CROSS01DEG+1)];
  *c02 = &bbr0cr2[(fn*hole_k+i)*(G2_CROSS02DEG+1)];
  *d00 = &bbr0[(fn*hole_k+j)*(G2_CROSS00DEG+1)];
  *d01 = &bbq0cr1[(fn*hole_k+i)*(G2_CROSS01DEG+1)];
  *d02 = &bbq0cr2[(fn*hole_k+i)*(G2_CROSS02DEG+1)];
} /*_g2h_GetBFAPatchCurvesf*/

void _g2h_GetBFBPatchCurvesf ( GHoleDomainf *domain, int fn, int i,
                    float **c00, float **c01, float **c02,
                    float **c10, float **c11, float **c12,
                    float **d00, float **d01, float **d02,
                    float **d10, float **d11, float **d12 )
{
  int   hole_k, nfunc_b, j;
  G2HolePrivateRecf *privateG2;
  float *bbr0, *bbr1, *bbq1, *bbr0cr1, *bbr1cr1, *bbq0cr1, *bbq1cr1,
        *bbr0cr2, *bbr1cr2, *bbq0cr2, *bbq1cr2;

  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  nfunc_b = privateG2->nfunc_b;

  G2GetBFuncBCrossAddresses();

  j = (i+1) % hole_k;
  *c00 = &bbr0[(fn*hole_k+i)*(G2_CROSS00DEG+1)];
  *c01 = &bbr0cr1[(fn*hole_k+i)*(G2_CROSS01DEG+1)];
  *c02 = &bbr0cr2[(fn*hole_k+i)*(G2_CROSS02DEG+1)];
  *c10 = &bbr1[(fn*hole_k+i)*(G2_CROSS10DEG+1)];
  *c11 = &bbr1cr1[(fn*hole_k+i)*(G2_CROSS11DEG+1)];
  *c12 = &bbr1cr2[(fn*hole_k+i)*(G2_CROSS12DEG+1)];
  *d00 = &bbr0[(fn*hole_k+j)*(G2_CROSS00DEG+1)];
  *d01 = &bbq0cr1[(fn*hole_k+i)*(G2_CROSS01DEG+1)];
  *d02 = &bbq0cr2[(fn*hole_k+i)*(G2_CROSS02DEG+1)];
  *d10 = &bbq1[(fn*hole_k+i)*(G2_CROSS10DEG+1)];
  *d11 = &bbq1cr1[(fn*hole_k+i)*(G2_CROSS11DEG+1)];
  *d12 = &bbq1cr2[(fn*hole_k+i)*(G2_CROSS12DEG+1)];
} /*_g2h_GetBFBPatchCurvesf*/

