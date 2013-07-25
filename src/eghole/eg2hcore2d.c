
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2008                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

void _g2h_GetDiPatchCurvesd ( GHoleDomaind *domain, int i,
                    point2d **c00, vector2d **c01, vector2d **c02,
                    point2d **c10, vector2d **c11, vector2d **c12,
                    point2d **d00, vector2d **d01, vector2d **d02,
                    point2d **d10, vector2d **d11, vector2d **d12 )
{
  int       hole_k, j;
  GHolePrivateRecd *privateG;
  G2HolePrivateRecd *privateG2;
  vector2d *omc, *surrpc;
  vector2d *dir0cr1, *diq0cr1, *dir1cr1, *diq1cr1,
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
} /*_g2h_GetDiPatchCurvesd*/

void _g2h_GetBFAPatchCurvesd ( GHoleDomaind *domain, int fn, int i,
                    double **c00, double **c01, double **c02,
                    double **d00, double **d01, double **d02 )
{
  int    hole_k, nfunc_a, j;
  G2HolePrivateRecd *privateG2;
  double *bbr0, *bbr0cr1, *bbq0cr1, *bbr0cr2, *bbq0cr2;

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
} /*_g2h_GetBFAPatchCurvesd*/

void _g2h_GetBFBPatchCurvesd ( GHoleDomaind *domain, int fn, int i,
                    double **c00, double **c01, double **c02,
                    double **c10, double **c11, double **c12,
                    double **d00, double **d01, double **d02,
                    double **d10, double **d11, double **d12 )
{
  int    hole_k, nfunc_b, j;
  G2HolePrivateRecd *privateG2;
  double *bbr0, *bbr1, *bbq1, *bbr0cr1, *bbr1cr1, *bbq0cr1, *bbq1cr1,
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
} /*_g2h_GetBFBPatchCurvesd*/

