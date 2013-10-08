
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

#include "eg1holef.h"
#include "eg1hprivatef.h"
#include "eg1herror.h"

/* ////////////////////////////////////////////////////////////////////////// */
void g1h_DrawDomAuxPatchesf ( GHoleDomainf *domain,
               void (*drawpatch) ( int n, int m, const point2f *cp ) )
{
  void    *sp;
  G1HolePrivateRecf *privateG1;
  int     i, j, hole_k;
  point2f *omc, *omcd;
  point2f *auxp;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  hole_k = domain->hole_k;
  omc = privateG1->omc;
  omcd = &omc[(G1H_OMCDEG+1)*hole_k];

  auxp = (point2f*)pkv_GetScratchMem ( 10*sizeof(point2f) );
  if ( !auxp )
    goto finish;

  for ( i = 0; i < hole_k; i++ ) {
    memcpy ( auxp, omc, (G1H_OMCDEG+1)*sizeof(point2f) );
    if ( !mbs_BCDegElevC2f ( 3, omcd, 1, &j, &auxp[(G1H_OMCDEG+1)] ) )
      goto finish;
    for ( j = 0; j <= G1H_OMCDEG; j++ )
      AddVector2f ( &auxp[j], &auxp[(G1H_OMCDEG+1)+j], &auxp[(G1H_OMCDEG+1)+j] );
    drawpatch ( 1, G1H_OMCDEG, auxp );
    omc += (G1H_OMCDEG+1);
    omcd += G1H_OMCDEG;
  }

finish:
  pkv_SetScratchMemTop ( sp );
} /*g1h_DrawDomAuxPatchesf*/

void g1h_DrawBasAuxPatchesf ( GHoleDomainf *domain, int fn,
               void (*drawpatch) ( int n, int m, const float *cp ) )
{
  void    *sp;
  G1HolePrivateRecf *privateG1;
  int     i, j, hole_k;
  float   *bezfc, *omc, *omcd;
  float   *auxp;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  hole_k = domain->hole_k;

  auxp = pkv_GetScratchMemf ( 24+(2*G1H_OMCDEG+1)*hole_k );
  if ( !auxp )
    goto finish;

  omc = &auxp[24];
  omcd = &omc[(G1H_OMCDEG+1)*hole_k];
  if ( fn < privateG1->nfunc_a )
    _g1h_GetABasisAuxpf ( domain, fn, omc, omcd );
  else {
    bezfc = pkv_GetScratchMemf ( 32*hole_k );
    if ( !bezfc )
      goto finish;
    _g1h_GetBBasisAuxpf ( domain, fn-privateG1->nfunc_a, bezfc, omc, omcd );
  }

  for ( i = 0; i < hole_k; i++ ) {
    memcpy ( auxp, omc, (G1H_OMCDEG+1)*sizeof(float) );
    if ( !mbs_BCDegElevC1f ( G1H_OMCDEG-1, omcd, 1, &j, &auxp[G1H_OMCDEG+1] ) )
      goto finish;
    for ( j = 0; j <= G1H_OMCDEG; j++ )
      auxp[G1H_OMCDEG+1+j] += auxp[j];
    drawpatch ( 1, G1H_OMCDEG, auxp );
    omc += (G1H_OMCDEG+1);
    omcd += G1H_OMCDEG;
  }

finish:
  pkv_SetScratchMemTop ( sp );
} /*g1h_DrawBasAuxPatchesf*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean g1h_DrawJFunctionf ( GHoleDomainf *domain, int k, int l,  
                             void (*drawpoly) ( int deg, const float *f ) )
{
  void  *sp;
  G1HolePrivateRecf *privateG1;
  int   hole_k, deg;
  float *b01, *c01, *f01, *g01, *b11, *c11, *f11, *g11, *f;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  if ( !privateG1->jfunc )
    goto failure;
  G1GetPolyAddr ( privateG1->jfunc, b01, c01, f01, g01, b11, c11, f11, g11 );

  switch ( l ) {
 case 0: f = b01;     deg = G1_BF01DEG;    break;
 case 1: f = c01;     deg = G1_CG01DEG;    break;
 case 2: f = f01;     deg = G1_BF01DEG;    break;
 case 3: f = g01;     deg = G1_CG01DEG;    break;
 case 4: f = b11;     deg = G1_BF11DEG;    break;
 case 5: f = c11;     deg = G1_CG11DEG;    break;
 case 6: f = f11;     deg = G1_BF11DEG;    break;
 case 7: f = g11;     deg = G1_CG11DEG;    break;
default:
    PKV_SIGNALERROR ( LIB_EGHOLE, ERRCODE_3, ERRMSG_3 );
    goto failure;
  }
  drawpoly ( deg, &f[k*(deg+1)] );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_DrawJFunctionf*/

/* ////////////////////////////////////////////////////////////////////////// */
void g1h_DrawDiPatchesf ( GHoleDomainf *domain,
                      void (*drawpatch) ( int n, int m, const point2f *cp ) )
{
  void     *sp;
  G1HolePrivateRecf *privateG1;
  int      hole_k, i, degu, degv;
  vector2f *c00, *c01, *c10, *c11, *d00, *d01, *d10, *d11;
  point2f  *di;

  privateG1 = domain->privateG1;
  if ( !privateG1->omc || !privateG1->dicross )
    return;
  sp = pkv_GetScratchMemTop ();
  di = (point2f*)pkv_GetScratchMem (
         (G1H_FINALDEG+1)*(G1H_FINALDEG+1)*sizeof(point2f) );
  if ( di ) {
    hole_k = domain->hole_k;

    for ( i = 0; i < hole_k; i++ ) {
      _g1h_GetDiPatchCurvesf ( domain, i,
                               &c00, &c01, &c10, &c11, &d00, &d01, &d10, &d11 );
      mbs_BezC1CoonsToBezf ( 2,
         G1_CROSS00DEG, (float*)c00, G1_CROSS01DEG, (float*)c01,
         3, (float*)c10, G1_CROSS11DEG, (float*)c11,
         G1_CROSS00DEG, (float*)d00, G1_CROSS01DEG, (float*)d01,
         3, (float*)d10, G1_CROSS11DEG, (float*)d11,
         &degu, &degv, (float*)di );
      drawpatch ( degu, degv, di );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g1h_DrawDiPatchesf*/

/* ////////////////////////////////////////////////////////////////////////// */
void g1h_ExtractPartitionf ( GHoleDomainf *domain,
                             int *hole_k, int *hole_m,
                             float *partition,
                             float *part_delta,
                             float *spart_alpha,
                             float *spart_malpha,
                             float *spart_salpha,
                             float *spart_knot,
                             float *alpha0,
                             boolean *spart_sgn,
                             boolean *spart_both )
{
  G1HolePrivateRecf *privateG1;
  int i, hk, hm;

  privateG1 = domain->privateG1;
  *hole_k = hk = domain->hole_k;
  *hole_m = hm = privateG1->hole_m;

  if ( partition )
    for ( i = 0; i < hk; i++ )
      partition[i] = privateG1->partition[i];
  if ( part_delta )
    for ( i = 0; i < hk; i++ )
      part_delta[i] = privateG1->partition[hk+i];
  if ( spart_alpha )
    for ( i = 0; i < hk-hm; i++ )
      spart_alpha[i] = privateG1->spartition[i].alpha;
  if ( spart_malpha )
    for ( i = 0; i < hk-hm; i++ )
      spart_malpha[i] = privateG1->spartition[i].malpha;
  if ( spart_salpha )
    for ( i = 0; i < hk-hm; i++ )
      spart_salpha[i] = privateG1->spartition[i].salpha;
  if ( spart_knot )
    for ( i = 0; i < hk-hm; i++ )
      spart_knot[i] = privateG1->spartition[i].knot;
  if ( spart_sgn )
    for ( i = 0; i < hk-hm; i++ )
      spart_sgn[i] = privateG1->spartition[i].sgn;
  if ( spart_both )
    for ( i = 0; i < hk-hm; i++ )
      spart_both[i] = privateG1->spartition[i].both;
  if ( alpha0 )
    *alpha0 = privateG1->spart_alpha0;
} /*g1h_ExtractPartitionf*/

/* ////////////////////////////////////////////////////////////////////////// */
void g1h_ExtractCentralPointf ( GHoleDomainf *domain, 
                                point2f *centp, vector2f *centder )
{
  G1HolePrivateRecf *privateG1;
  vector2f          *omcbc;
  int               hole_k, i;

  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  if ( privateG1->omcbc ) {
    omcbc = privateG1->omcbc;
    *centp = omcbc[0];
    if ( centder )
      for ( i = 0; i < hole_k; i++ )
        centder[i] = omcbc[(G1H_OMCDEG+1)*i+1];
  }
} /*g1h_ExtractCentralPointf*/

/* ////////////////////////////////////////////////////////////////////////// */
void g1h_DrawBasAFunctionf ( GHoleDomainf *domain, int fn,
               void (*drawpatch) ( int n, int m, const point3f *cp ) )
{
  void    *sp;
  int     hole_k, ncp, i, degu, degv;
  float   *c00, *c01, *c10, *c11, *d00, *d01, *d10, *d11;
  point2f *dicp;
  float   *fcp;
  point3f *dficp;
  float   zero[4] = {0.0,0.0,0.0,0.0};

  sp = pkv_GetScratchMemTop ();
  ncp = (G1H_FINALDEG+1)*(G1H_FINALDEG+1);
  dicp = (point2f*)pkv_GetScratchMem ( ncp*sizeof(point2f) );
  fcp = pkv_GetScratchMemf ( ncp );
  dficp = (point3f*)pkv_GetScratchMem ( ncp*sizeof(point3f) );
  if ( dicp && fcp && dficp ) {
    hole_k = domain->hole_k;

    for ( i = 0; i < hole_k; i++ ) {
      _g1h_GetDiPatchCurvesf ( domain, i,
                     (point2f**)(void*)&c00, (point2f**)(void*)&c01,
                     (point2f**)(void*)&c10, (point2f**)(void*)&c11,
                     (point2f**)(void*)&d00, (point2f**)(void*)&d01,
                     (point2f**)(void*)&d10, (point2f**)(void*)&d11 );
      mbs_BezC1CoonsToBezf ( 2,
                     G1_CROSS00DEG, c00, G1_CROSS01DEG, c01,
                     3, c10, G1_CROSS11DEG, c11,
                     G1_CROSS00DEG, d00, G1_CROSS01DEG, d01,
                     3, d10, G1_CROSS11DEG, d11,
                     &degu, &degv, (float*)dicp );
      _g1h_GetBFAPatchCurvesf ( domain, fn, i, &c00, &c01, &d00, &d01 );
      mbs_BezC1CoonsToBezf ( 1,
                     G1_CROSS00DEG, c00, G1_CROSS01DEG, c01, 2, zero, 2, zero,
                     G1_CROSS00DEG, d00, G1_CROSS01DEG, d01, 2, zero, 2, zero,
                     &degu, &degv, fcp );
      pkv_Selectf ( ncp, 2, 2, 3, (float*)dicp, (float*)dficp );
      pkv_Selectf ( ncp, 1, 1, 3, (float*)fcp, (float*)&dficp[0].z );
      drawpatch ( degu, degv, dficp );
    }
  }

  pkv_SetScratchMemTop ( sp );
} /*g1h_DrawBasAFunctionf*/

void g1h_DrawBasBFunctionf ( GHoleDomainf *domain, int fn,
               void (*drawpatch) ( int n, int m, const point3f *cp ) )
{
  void    *sp;
  int     hole_k, ncp, i, j, degu, degv;
  float   *c00, *c01, *c10, *c11, *d00, *d01, *d10, *d11;
  point2f *dicp, *spcp;
  float   *fcp;
  point3f *dficp;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;

  ncp = (G1H_FINALDEG+1)*(G1H_FINALDEG+1);
  dicp = (point2f*)pkv_GetScratchMem ( ncp*sizeof(point2f) );
  fcp = pkv_GetScratchMemf ( ncp );
  dficp = (point3f*)pkv_GetScratchMem ( ncp*sizeof(point3f) );
  if ( dicp && fcp && dficp ) {
    for ( i = 0; i < hole_k; i++ )
      for ( j = 0; j < 3; j++ ) {
        if ( j == 0 ) {
          _gh_FindDomSurrndPatchf ( domain, i, j, dicp );
          spcp = dicp;
        }
        else
          spcp = _gh_GetDomSurrndPatchf ( domain, i, j );
        gh_GetDomSurrndBFuncf ( domain, fn, i, j, fcp );
        pkv_Selectf ( 16, 2, 2, 3, (float*)spcp, (float*)dficp );
        pkv_Selectf ( 16, 1, 1, 3, fcp, (float*)&dficp[0].z );
        drawpatch ( 3, 3, dficp );
      }

    for ( i = 0; i < hole_k; i++ ) {
      _g1h_GetDiPatchCurvesf ( domain, i,
                     (point2f**)(void*)&c00, (point2f**)(void*)&c01,
                     (point2f**)(void*)&c10, (point2f**)(void*)&c11,
                     (point2f**)(void*)&d00, (point2f**)(void*)&d01,
                     (point2f**)(void*)&d10, (point2f**)(void*)&d11 );
      mbs_BezC1CoonsToBezf ( 2,
                     G1_CROSS00DEG, c00, G1_CROSS01DEG, c01, 3, c10, G1_CROSS11DEG, c11,
                     G1_CROSS00DEG, d00, G1_CROSS01DEG, d01, 3, d10, G1_CROSS11DEG, d11,
                     &degu, &degv, (float*)dicp );
      _g1h_GetBFBPatchCurvesf ( domain, fn, i,
                                &c00, &c01, &c10, &c11, &d00, &d01, &d10, &d11 );
      mbs_BezC1CoonsToBezf ( 1,
                     G1_CROSS00DEG, c00, G1_CROSS01DEG, c01,
                     G1_CROSS10DEG, c10, G1_CROSS11DEG, c11,
                     G1_CROSS00DEG, d00, G1_CROSS01DEG, d01,
                     G1_CROSS10DEG, d10, G1_CROSS11DEG, d11,
                     &degu, &degv, fcp );
      pkv_Selectf ( ncp, 2, 2, 3, (float*)dicp, (float*)dficp );
      pkv_Selectf ( ncp, 1, 1, 3, fcp, (float*)&dficp[0].z );
      drawpatch ( degu, degv, dficp );
    }
  }

  pkv_SetScratchMemTop ( sp );
} /*g1h_DrawBasBFunctionf*/

boolean g1h_DrawBasCNetf ( GHoleDomainf *domain, int fn,   
               void (*drawnet) ( int n, int m, const point3f *cp ) )
{
  void    *sp;
  GHolePrivateRecf  *privateG;
  int     hole_k, i, j;
  int     *ind;
  point2f *domain_cp;
  point3f *cp;

  sp = pkv_GetScratchMemTop ();
  ind = pkv_GetScratchMem ( 16*sizeof(int) );
  cp  = pkv_GetScratchMem ( 16*sizeof(point3f) );
  if ( !ind || !cp ) {
    PKV_SIGNALERROR ( LIB_EGHOLE, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }

  privateG  = domain->privateG;
  domain_cp = domain->domain_cp;
  hole_k = domain->hole_k;
  for ( i = 0; i < hole_k; i++ ) {
    gh_GetBspInd ( hole_k, i, 0, ind );
    for ( j = 0; j < 16; j++ ) {
      SetPoint3f ( &cp[j], domain_cp[ind[j]].x, domain_cp[ind[j]].y, 0.0 );
      if ( fn >= 0 && privateG->bfcpn[fn] == ind[j] )
        cp[j].z = 1.0;
    }
    drawnet ( 3, 3, cp );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_DrawBasCNetf*/

/* ////////////////////////////////////////////////////////////////////////// */
void g1h_DrawBFAomcf ( GHoleDomainf *domain, int fn,
                       void (*drawpoly)(int degree, const float *coeff) )
{
  void  *sp;
  int   hole_k, i;
  float *c00, *c01, *d00, *d01, *c;

  sp = pkv_GetScratchMemTop ();
  if ( (c = pkv_GetScratchMemf ( G1H_OMCDEG+1 )) ) {
    hole_k = domain->hole_k;
    for ( i = 0; i < hole_k; i++ ) {
      _g1h_GetBFAPatchCurvesf ( domain, fn, i, &c00, &c01, &d00, &d01 );
      memcpy ( c, c00, (G1H_OMCDEG+1)*sizeof(float) );
      drawpoly ( G1H_OMCDEG, c );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g1h_DrawBFAomcf*/

void g1h_DrawBFBomcf ( GHoleDomainf *domain, int fn,
                       void (*drawpoly)(int degree, const float *coeff) )
{
  void  *sp;
  int   hole_k, i;
  float *c00, *c01, *c10, *c11, *d00, *d01, *d10, *d11, *c;

  sp = pkv_GetScratchMemTop ();
  if ( (c = pkv_GetScratchMemf ( G1H_OMCDEG+1 )) ) {
    hole_k = domain->hole_k;
    for ( i = 0; i < hole_k; i++ ) {
      _g1h_GetBFBPatchCurvesf ( domain, fn, i, &c00, &c01, &c10, &c11,
                                &d00, &d01, &d10, &d11 );
      memcpy ( c, c00, (G1H_OMCDEG+1)*sizeof(float) );
      drawpoly ( G1H_OMCDEG, c );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g1h_DrawBFBomcf*/

void g1h_DrawFinalSurfBCf ( GHoleDomainf *domain,
                            int spdimen, const float *hole_cp,
                            const float *acoeff,
                            void (*drawcurve)(int degree, int spdimen,
                                              const float *cp) )
{
  void *sp;
  int  hole_k, nfunc_a, nfunc_b;
  GHolePrivateRecf  *privateG;
  G1HolePrivateRecf *privateG1;
  unsigned char *bfcpn;
  float *x, *y;
  float *c00, *c01, *c10, *c11, *d00, *d01, *d10, *d11;
  int   i, j;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG = domain->privateG;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  bfcpn = privateG->bfcpn;
  if ( (x = pkv_GetScratchMemf ( 2*spdimen*(G1H_OMCDEG+1) )) ) {
    y = &x[spdimen*(G1H_OMCDEG+1)];
    for ( i = 0; i < hole_k; i++ ) {
      memset ( x, 0, spdimen*(G1H_OMCDEG+1)*sizeof(float) );
      for ( j = 0; j < nfunc_a; j++ ) {
        _g1h_GetBFAPatchCurvesf ( domain, j, i, &c00, &c01, &d00, &d01 );
        pkn_MultMatrixf ( G1H_OMCDEG+1, 1, 1, c00, spdimen, 0,
                          &acoeff[j*spdimen], spdimen, y );
        pkn_AddMatrixf ( 1, spdimen*(G1H_OMCDEG+1), 0, x, 0, y, 0, x );
      }
      for ( j = 0; j < nfunc_b; j++ ) {
        _g1h_GetBFBPatchCurvesf ( domain, j, i, &c00, &c01, &c10, &c11,
                                  &d00, &d01, &d10, &d11 );
        pkn_MultMatrixf ( G1H_OMCDEG+1, 1, 1, c00, spdimen, 0,
                          &hole_cp[bfcpn[j]], spdimen, y );
        pkn_AddMatrixf ( 1, spdimen*(G1H_OMCDEG+1), 0, x, 0, y, 0, x );
      }
      drawcurve ( G1H_OMCDEG, spdimen, x );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g1h_DrawFinalSurfBCf*/

void g1h_ExtDrawFinalSurfBCf ( GHoleDomainf *domain,
                               int spdimen, const float *hole_cp,
                               const float *acoeff,
                               void (*drawcurve)(int degree, int spdimen,
                                                 const float *cp) )
{
  void *sp;
  int  hole_k, nfunc_a, nfunc_b, nfunc_c;
  GHolePrivateRecf  *privateG;
  G1HolePrivateRecf *privateG1;
  unsigned char *bfcpn;
  float *x, *y;
  float *c00, *c01, *c10, *c11, *d00, *d01, *d10, *d11;
  int   i, j;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG = domain->privateG;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = hole_k*G1_DBDIM;
  bfcpn = privateG->bfcpn;
  if ( (x = pkv_GetScratchMemf ( 2*spdimen*(G1H_OMCDEG+1) )) ) {
    y = &x[spdimen*(G1H_OMCDEG+1)];
    for ( i = 0; i < hole_k; i++ ) {
      memset ( x, 0, spdimen*(G1H_OMCDEG+1)*sizeof(float) );
      for ( j = 0; j < nfunc_a; j++ ) {
        _g1h_GetBFAPatchCurvesf ( domain, j, i, &c00, &c01, &d00, &d01 );
        pkn_MultMatrixf ( G1H_OMCDEG+1, 1, 1, c00, spdimen, 0,
                          &acoeff[(nfunc_c+j)*spdimen], spdimen, y );
        pkn_AddMatrixf ( 1, spdimen*(G1H_OMCDEG+1), 0, x, 0, y, 0, x );
      }
      for ( j = 0; j < nfunc_b; j++ ) {
        _g1h_GetBFBPatchCurvesf ( domain, j, i, &c00, &c01, &c10, &c11,
                                  &d00, &d01, &d10, &d11 );
        pkn_MultMatrixf ( G1H_OMCDEG+1, 1, 1, c00, spdimen, 0,
                          &hole_cp[bfcpn[j]], spdimen, y );
        pkn_AddMatrixf ( 1, spdimen*(G1H_OMCDEG+1), 0, x, 0, y, 0, x );
      }
      drawcurve ( G1H_OMCDEG, spdimen, x );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g1h_ExtDrawFinalSurfBCf*/

/* ////////////////////////////////////////////////////////////////////////// */
void g1h_DrawMatricesf ( GHoleDomainf *domain,
                         void (*drawmatrix)(int nfa, int nfb,
                                            float *amat, float *bmat) )
{
  void              *sp;
  float             *amat, *bmat;
  int               na, nb;
  G1HolePrivateRecf *privateG1;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  if ( !(privateG1->Amat && privateG1->Bmat) )
    goto wayout;

  na = privateG1->nfunc_a;
  na = na*(na+1)/2;
  nb = privateG1->nfunc_a*privateG1->nfunc_b;
  amat = pkv_GetScratchMemf ( na );
  bmat = pkv_GetScratchMemf ( nb );
  if ( !amat || !bmat )
    goto wayout;

  memcpy ( amat, privateG1->Amat, na*sizeof(float) );
  memcpy ( bmat, privateG1->Bmat, nb*sizeof(float) );
  drawmatrix ( privateG1->nfunc_a, privateG1->nfunc_b, amat, bmat );

wayout:
  pkv_SetScratchMemTop ( sp );
} /*g1h_DrawMatricesf*/

/* ////////////////////////////////////////////////////////////////////////// */
void g1h_DrawExtMatricesf ( GHoleDomainf *domain,
                            void (*drawmatrix)(int k, int r, int s, float *Aii,
                                               float *Bi) )
{
  void  *sp;
  float *aii, *bi, *Aii, *Aki, *Akk, *Bi, *Bk, *Lii;
  G1HolePrivateRecf *privateG1;
  int   hole_k, nfunc_a, nfunc_b, nfunc_c, size;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  if ( !(privateG1->EAmat && privateG1->EBmat) )
    goto wayout;

  _g1h_GetExtBlockAddressesf ( domain, &Aii, &Aki, &Akk, &Bi, &Bk, &Lii );
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = hole_k*G1_DBDIM;

  aii = pkv_GetScratchMemf ( (size = pkn_Block1ArraySize ( hole_k, G1_DBDIM, nfunc_a)) );
  bi  = pkv_GetScratchMemf ( (nfunc_c+nfunc_a)*nfunc_b );
  if ( aii && bi ) {
    memcpy ( aii, Aii, size*sizeof(float) );
    memcpy ( bi, Bi, (nfunc_c+nfunc_a)*nfunc_b*sizeof(float) );
    drawmatrix ( hole_k, G1_DBDIM, nfunc_a, aii, bi );
  }
wayout:
  pkv_SetScratchMemTop ( sp );
} /*g1h_DrawExtMatricesf*/

