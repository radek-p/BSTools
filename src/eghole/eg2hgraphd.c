
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
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
#include "eg2holed.h"

#include "eg2hprivated.h"
#include "eg2herror.h"

/* ////////////////////////////////////////////////////////////////////////// */
void g2h_DrawDomAuxPatchesd ( GHoleDomaind *domain,
               void (*drawpatch) ( int n, int m, const point2d *cp ) )
{
  void    *sp;
  G2HolePrivateRecd *privateG2;
  int     i, j, hole_k;
  point2d *omc, *omcd, *omcdd;
  point2d *auxp;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  hole_k = domain->hole_k;
  omc = privateG2->omc;
  omcd = &omc[(G2H_OMCDEG+1)*hole_k];
  omcdd = &omcd[G2H_OMCDEG*hole_k];

  auxp = (point2d*)pkv_GetScratchMem ( 24*sizeof(point2d) );
  if ( !auxp )
    goto finish;

  for ( i = 0; i < hole_k; i++ ) {
    memcpy ( auxp, omc, (G2H_OMCDEG+1)*sizeof(point2d) );
    mbs_BCDegElevC2d ( 6, omcd, 1, &j, &auxp[8] );
    mbs_BCDegElevC2d ( 5, omcdd, 2, &j, &auxp[16] );
    for ( j = 0; j <= G2H_OMCDEG; j++ ) {
      MultVector2d ( 0.5, &auxp[8+j], &auxp[8+j] );
      AddVector2Md ( &auxp[8+j], &auxp[16+j], 0.5, &auxp[16+j] );
      AddVector2d ( &auxp[j], &auxp[8+j], &auxp[8+j] );
      AddVector2d ( &auxp[8+j], &auxp[16+j], &auxp[16+j] );
    }
    drawpatch ( 2, G2H_OMCDEG, auxp );
    omc += (G2H_OMCDEG+1);
    omcd += G2H_OMCDEG;
    omcdd += (G2H_OMCDEG-1);
  }

finish:
  pkv_SetScratchMemTop ( sp );
} /*g2h_DrawDomAuxPatchesd*/

void g2h_DrawBasAuxPatchesd ( GHoleDomaind *domain, int fn,
               void (*drawpatch) ( int n, int m, const double *cp ) )
{
  void    *sp;
  G2HolePrivateRecd *privateG2;
  int     i, j, hole_k;
  double  *bezfc, *omc, *omcd, *omcdd;
  double  *auxp;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  hole_k = domain->hole_k;

  auxp = pkv_GetScratchMemd ( 24+3*G2H_OMCDEG*hole_k );
  if ( !auxp )
    goto finish;

  omc = &auxp[24];
  omcd = &omc[(G2H_OMCDEG+1)*hole_k];
  omcdd = &omcd[G2H_OMCDEG*hole_k];
  if ( fn < privateG2->nfunc_a )
    _g2h_GetABasisAuxpd ( domain, fn, omc, omcd, omcdd );
  else {
    bezfc = pkv_GetScratchMemd ( 32*hole_k );
    if ( !bezfc )
      goto finish;
    _g2h_GetBBasisAuxpd ( domain, fn-privateG2->nfunc_a, bezfc, omc, omcd, omcdd );
  }

  for ( i = 0; i < hole_k; i++ ) {
    memcpy ( auxp, omc, (G2H_OMCDEG+1)*sizeof(double) );
    mbs_BCDegElevC1d ( 6, omcd, 1, &j, &auxp[8] );
    mbs_BCDegElevC1d ( 5, omcdd, 2, &j, &auxp[16] );
    for ( j = 0; j <= G2H_OMCDEG; j++ ) {
      auxp[8+j] *= 0.5;
      auxp[16+j] = (double)(auxp[8+j] + 0.5*auxp[16+j]);
      auxp[8+j] += auxp[j];
      auxp[16+j] += auxp[8+j];
    }
    drawpatch ( 2, G2H_OMCDEG, auxp );
    omc += (G2H_OMCDEG+1);
    omcd += G2H_OMCDEG;
    omcdd += (G2H_OMCDEG-1);
  }

finish:
  pkv_SetScratchMemTop ( sp );
} /*g2h_DrawBasAuxPatchesd*/

/* ////////////////////////////////////////////////////////////////////////// */
void g2h_DrawJFunctiond ( GHoleDomaind *domain, int k, int l,  
                          void (*drawpoly) ( int deg, const double *f ) )
{
  void   *sp;
  G2HolePrivateRecd *privateG2;
  int    hole_k, deg;
  double *b01, *c01, *f01, *g01, *b11, *c11, *f11, *g11,
         *b02, *c02, *f02, *g02, *b12, *c12, *f12, *g12,
         *b01b01, *twob01c01, *c01c01, *f01f01, *twof01g01, *g01g01,
         *b11b11, *twob11c11, *c11c11, *f11f11, *twof11g11, *g11g11, *f;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  if ( !privateG2->jfunc )
    return;
  G2GetPolynomialAddresses ( privateG2->jfunc,
        b01, c01, f01, g01, b11, c11, f11, g11,
        b02, c02, f02, g02, b12, c12, f12, g12, b01b01, twob01c01, c01c01,
        f01f01, twof01g01, g01g01, b11b11, twob11c11, c11c11,
        f11f11, twof11g11, g11g11 );

  switch ( l ) {
 case 0: f = b01;     deg = G2_BF01DEG;    break;
 case 1: f = c01;     deg = G2_CG01DEG;    break;
 case 2: f = f01;     deg = G2_BF01DEG;    break;
 case 3: f = g01;     deg = G2_CG01DEG;    break;
 case 4: f = b02;     deg = G2_BF02DEG;    break;
 case 5: f = c02;     deg = G2_CG02DEG;    break;
 case 6: f = f02;     deg = G2_BF02DEG;    break;
 case 7: f = g02;     deg = G2_CG02DEG;    break;
 case 8: f = b11;     deg = G2_BF11DEG;    break;
 case 9: f = c11;     deg = G2_CG11DEG;    break;
case 10: f = f11;     deg = G2_BF11DEG;    break;
case 11: f = g11;     deg = G2_CG11DEG;    break;
case 12: f = b12;     deg = G2_BF12DEG;    break;
case 13: f = c12;     deg = G2_CG12DEG;    break;
case 14: f = f12;     deg = G2_BF12DEG;    break;
case 15: f = g12;     deg = G2_CG12DEG;    break;
case 16: f = b01b01;  deg = 2*G2_BF01DEG;  break;
case 17: f = twob01c01;  deg = G2_BF01DEG+G2_CG01DEG;  break;
case 18: f = c01c01;  deg = 2*G2_CG01DEG;  break;
case 19: f = f01f01;  deg = 2*G2_BF01DEG;  break;
case 20: f = twof01g01;  deg = G2_BF01DEG+G2_CG01DEG;  break;
case 21: f = g01g01;  deg = 2*G2_CG01DEG;  break;
case 22: f = b11b11;  deg = 2*G2_BF11DEG;  break;
case 23: f = twob11c11;  deg = G2_BF11DEG+G2_CG11DEG;  break;
case 24: f = c11c11;  deg = 2*G2_CG11DEG;  break;
case 25: f = f11f11;  deg = 2*G2_BF11DEG;  break;
case 26: f = twof11g11;  deg = G2_BF11DEG+G2_CG11DEG;  break;
case 27: f = g11g11;  deg = 2*G2_CG11DEG;  break;
default: exit ( 1 );
  }
  drawpoly ( deg, &f[k*(deg+1)] );
  pkv_SetScratchMemTop ( sp );
} /*g2h_DrawJFunctiond*/

/* ////////////////////////////////////////////////////////////////////////// */
void g2h_DrawDiPatchesd ( GHoleDomaind *domain,
                      void (*drawpatch) ( int n, int m, const point2d *cp ) )
{
  void     *sp;
  G2HolePrivateRecd *privateG2;
  int      hole_k, i, degu, degv;
  vector2d *c00, *c01, *c02, *c10, *c11, *c12,
           *d00, *d01, *d02, *d10, *d11, *d12;
  point2d  *di;

  privateG2 = domain->privateG2;
  if ( !privateG2->omc || !privateG2->dicross )
    return;
  sp = pkv_GetScratchMemTop ();
  di = (point2d*)pkv_GetScratchMem (
         (G2H_FINALDEG+1)*(G2H_FINALDEG+1)*sizeof(point2d) );
  if ( di ) {
    hole_k = domain->hole_k;

    for ( i = 0; i < hole_k; i++ ) {
      _g2h_GetDiPatchCurvesd ( domain, i,
                               &c00, &c01, &c02, &c10, &c11, &c12,
                               &d00, &d01, &d02, &d10, &d11, &d12 );
      mbs_BezC2CoonsToBezd ( 2,
         G2_CROSS00DEG, (double*)c00, G2_CROSS01DEG, (double*)c01, G2_CROSS02DEG, (double*)c02,
         3, (double*)c10, G2_CROSS11DEG, (double*)c11, G2_CROSS12DEG, (double*)c12,
         G2_CROSS00DEG, (double*)d00, G2_CROSS01DEG, (double*)d01, G2_CROSS02DEG, (double*)d02,
         3, (double*)d10, G2_CROSS11DEG, (double*)d11, G2_CROSS12DEG, (double*)d12,
         &degu, &degv, (double*)di );
      drawpatch ( degu, degv, di );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g2h_DrawDiPatchesd*/

/* ////////////////////////////////////////////////////////////////////////// */
void g2h_ExtractPartitiond ( GHoleDomaind *domain,
                             int *hole_k, int *hole_m,
                             double *partition,
                             double *part_delta,
                             double *spart_alpha,
                             double *spart_malpha,
                             double *spart_salpha,
                             double *spart_knot,
                             double *alpha0,
                             boolean *spart_sgn,
                             boolean *spart_both )
{
  G2HolePrivateRecd *privateG2;
  int i, hk, hm;

  privateG2 = domain->privateG2;
  *hole_k = hk = domain->hole_k;
  *hole_m = hm = privateG2->hole_m;

  if ( partition )
    for ( i = 0; i < hk; i++ )
      partition[i] = privateG2->partition[i];
  if ( part_delta )
    for ( i = 0; i < hk; i++ )
      part_delta[i] = privateG2->partition[hk+i];
  if ( spart_alpha )
    for ( i = 0; i < hk-hm; i++ )
      spart_alpha[i] = privateG2->spartition[i].alpha;
  if ( spart_malpha )
    for ( i = 0; i < hk-hm; i++ )
      spart_malpha[i] = privateG2->spartition[i].malpha;
  if ( spart_salpha )
    for ( i = 0; i < hk-hm; i++ )
      spart_salpha[i] = privateG2->spartition[i].salpha;
  if ( spart_knot )
    for ( i = 0; i < hk-hm; i++ )
      spart_knot[i] = privateG2->spartition[i].knot;
  if ( spart_sgn )
    for ( i = 0; i < hk-hm; i++ )
      spart_sgn[i] = privateG2->spartition[i].sgn;
  if ( spart_both )
    for ( i = 0; i < hk-hm; i++ )
      spart_both[i] = privateG2->spartition[i].both;
  if ( alpha0 )
    *alpha0 = privateG2->spart_alpha0;
} /*g2h_ExtractPartitiond*/

/* ////////////////////////////////////////////////////////////////////////// */
void g2h_ExtractCentralPointd ( GHoleDomaind *domain, 
                                point2d *centp, vector2d *centder )
{
  G2HolePrivateRecd *privateG2;
  vector2d          *omcbc;
  int               hole_k, i;

  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  if ( privateG2->omcbc ) {
    omcbc = privateG2->omcbc;
    *centp = omcbc[0];
    if ( centder )
      for ( i = 0; i < hole_k; i++ )
        centder[i] = omcbc[(G2H_OMCDEG+1)*i+1];
  }
} /*g2h_ExtractCentralPointd*/

/* ////////////////////////////////////////////////////////////////////////// */
void g2h_DrawBasAFunctiond ( GHoleDomaind *domain, int fn,
               void (*drawpatch) ( int n, int m, const point3d *cp ) )
{
  void    *sp;
  int     hole_k, ncp, i, degu, degv;
  double  *c00, *c01, *c02, *c10, *c11, *c12,
          *d00, *d01, *d02, *d10, *d11, *d12;
  point2d *dicp;
  double  *fcp;
  point3d *dficp;
  double  zero[4] = {0.0,0.0,0.0,0.0};

  sp = pkv_GetScratchMemTop ();
  ncp = (G2H_FINALDEG+1)*(G2H_FINALDEG+1);
  dicp = (point2d*)pkv_GetScratchMem ( ncp*sizeof(point2d) );
  fcp = pkv_GetScratchMemd ( ncp );
  dficp = (point3d*)pkv_GetScratchMem ( ncp*sizeof(point3d) );
  if ( dicp && fcp && dficp ) {
    hole_k = domain->hole_k;

    for ( i = 0; i < hole_k; i++ ) {
      _g2h_GetDiPatchCurvesd ( domain, i,
           (point2d**)(void*)&c00, (point2d**)(void*)&c01, (point2d**)(void*)&c02,
           (point2d**)(void*)&c10, (point2d**)(void*)&c11, (point2d**)(void*)&c12,
           (point2d**)(void*)&d00, (point2d**)(void*)&d01, (point2d**)(void*)&d02,
           (point2d**)(void*)&d10, (point2d**)(void*)&d11, (point2d**)(void*)&d12 );
      mbs_BezC2CoonsToBezd ( 2,
                     G2_CROSS00DEG, c00, G2_CROSS01DEG, c01, G2_CROSS02DEG, c02,
                     3, c10, G2_CROSS11DEG, c11, G2_CROSS12DEG, c12,
                     G2_CROSS00DEG, d00, G2_CROSS01DEG, d01, G2_CROSS02DEG, d02,
                     3, d10, G2_CROSS11DEG, d11, G2_CROSS12DEG, d12,
                     &degu, &degv, (double*)dicp );
      _g2h_GetBFAPatchCurvesd ( domain, fn, i,
                                &c00, &c01, &c02, &d00, &d01, &d02 );
      mbs_BezC2CoonsToBezd ( 1,
                     G2_CROSS00DEG, c00, G2_CROSS01DEG, c01, G2_CROSS02DEG, c02,
                     2, zero, 2, zero, 2, zero,
                     G2_CROSS00DEG, d00, G2_CROSS01DEG, d01, G2_CROSS02DEG, d02,
                     2, zero, 2, zero, 2, zero,
                     &degu, &degv, fcp );
      pkv_Selectd ( ncp, 2, 2, 3, (double*)dicp, (double*)dficp );
      pkv_Selectd ( ncp, 1, 1, 3, (double*)fcp, (double*)&dficp[0].z );
      drawpatch ( degu, degv, dficp );
    }
  }

  pkv_SetScratchMemTop ( sp );
} /*g2h_DrawBasAFunctiond*/

void g2h_DrawBasBFunctiond ( GHoleDomaind *domain, int fn,
               void (*drawpatch) ( int n, int m, const point3d *cp ) )
{
  void    *sp;
  int     hole_k, ncp, i, j, degu, degv;
  double  *c00, *c01, *c02, *c10, *c11, *c12,
          *d00, *d01, *d02, *d10, *d11, *d12;
  point2d *dicp, *spcp;
  double  *fcp;
  point3d *dficp;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  ncp = (G2H_FINALDEG+1)*(G2H_FINALDEG+1);
  dicp = (point2d*)pkv_GetScratchMem ( ncp*sizeof(point2d) );
  fcp = pkv_GetScratchMemd ( ncp );
  dficp = (point3d*)pkv_GetScratchMem ( ncp*sizeof(point3d) );
  if ( dicp && fcp && dficp ) {
    for ( i = 0; i < hole_k; i++ )
      for ( j = 0; j < 3; j++ ) {
        if ( j == 0 ) {
          _gh_FindDomSurrndPatchd ( domain, i, j, dicp );
          spcp = dicp;
        }
        else
          spcp = _gh_GetDomSurrndPatchd ( domain, i, j );
        gh_GetDomSurrndBFuncd ( domain, fn, i, j, fcp );
        pkv_Selectd ( 16, 2, 2, 3, (double*)spcp, (double*)dficp );
        pkv_Selectd ( 16, 1, 1, 3, fcp, (double*)&dficp[0].z );
        drawpatch ( 3, 3, dficp );
      }

    for ( i = 0; i < hole_k; i++ ) {
      _g2h_GetDiPatchCurvesd ( domain, i,
           (point2d**)(void*)&c00, (point2d**)(void*)&c01, (point2d**)(void*)&c02,
           (point2d**)(void*)&c10, (point2d**)(void*)&c11, (point2d**)(void*)&c12,
           (point2d**)(void*)&d00, (point2d**)(void*)&d01, (point2d**)(void*)&d02,
           (point2d**)(void*)&d10, (point2d**)(void*)&d11, (point2d**)(void*)&d12 );
      mbs_BezC2CoonsToBezd ( 2,
                     G2_CROSS00DEG, c00, G2_CROSS01DEG, c01, G2_CROSS02DEG, c02,
                     3, c10, G2_CROSS11DEG, c11, G2_CROSS12DEG, c12,
                     G2_CROSS00DEG, d00, G2_CROSS01DEG, d01, G2_CROSS02DEG, d02,
                     3, d10, G2_CROSS11DEG, d11, G2_CROSS12DEG, d12,
                     &degu, &degv, (double*)dicp );
      _g2h_GetBFBPatchCurvesd ( domain, fn, i,
                                &c00, &c01, &c02, &c10, &c11, &c12,
                                &d00, &d01, &d02, &d10, &d11, &d12 );
      mbs_BezC2CoonsToBezd ( 1,
                     G2_CROSS00DEG, c00, G2_CROSS01DEG, c01, G2_CROSS02DEG, c02,
                     G2_CROSS10DEG, c10, G2_CROSS11DEG, c11, G2_CROSS12DEG, c12,
                     G2_CROSS00DEG, d00, G2_CROSS01DEG, d01, G2_CROSS02DEG, d02,
                     G2_CROSS10DEG, d10, G2_CROSS11DEG, d11, G2_CROSS12DEG, d12,
                     &degu, &degv, fcp );
      pkv_Selectd ( ncp, 2, 2, 3, (double*)dicp, (double*)dficp );
      pkv_Selectd ( ncp, 1, 1, 3, fcp, (double*)&dficp[0].z );
      drawpatch ( degu, degv, dficp );
    }
  }

  pkv_SetScratchMemTop ( sp );
} /*g2h_DrawBasBFunctiond*/

void g2h_DrawBasCNetd ( GHoleDomaind *domain, int fn,   
               void (*drawnet) ( int n, int m, const point3d *cp ) )
{
  void    *sp;
  GHolePrivateRecd  *privateG;
  int     hole_k, i, j;
  int     *ind;
  point2d *domain_cp;
  point3d *cp;

  sp = pkv_GetScratchMemTop ();
  ind = pkv_GetScratchMem ( 16*sizeof(int) );
  cp  = pkv_GetScratchMem ( 16*sizeof(point3d) );
  if ( !ind || !cp )
    exit ( 1 );

  privateG  = domain->privateG;
  domain_cp = domain->domain_cp;
  hole_k = domain->hole_k;
  for ( i = 0; i < hole_k; i++ ) {
    gh_GetBspInd ( hole_k, i, 0, ind );
    for ( j = 0; j < 16; j++ ) {
      SetPoint3d ( &cp[j], domain_cp[ind[j]].x, domain_cp[ind[j]].y, 0.0 );
      if ( fn >= 0 && privateG->bfcpn[fn] == ind[j] )
        cp[j].z = 1.0;
    }
    drawnet ( 3, 3, cp );
  }
  pkv_SetScratchMemTop ( sp );
} /*g2h_DrawBasCNetd*/

/* ////////////////////////////////////////////////////////////////////////// */
void g2h_DrawBFAomcd ( GHoleDomaind *domain, int fn,
                       void (*drawpoly)(int degree, const double *coeff) )
{
  void  *sp;
  int   hole_k, i;
  double *c00, *c01, *c02, *d00, *d01, *d02, *c;

  sp = pkv_GetScratchMemTop ();
  if ( (c = pkv_GetScratchMemd ( G2H_OMCDEG+1 )) ) {
    hole_k = domain->hole_k;
    for ( i = 0; i < hole_k; i++ ) {
      _g2h_GetBFAPatchCurvesd ( domain, fn, i,
                                &c00, &c01, &c02, &d00, &d01, &d02 );
      memcpy ( c, c00, (G2H_OMCDEG+1)*sizeof(double) );
      drawpoly ( G2H_OMCDEG, c );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g2h_DrawBFAomcd*/

void g2h_DrawBFBomcd ( GHoleDomaind *domain, int fn,
                       void (*drawpoly)(int degree, const double *coeff) )
{
  void  *sp;
  int   hole_k, i;
  double *c00, *c01, *c02, *c10, *c11, *c12,
        *d00, *d01, *d02, *d10, *d11, *d12, *c;

  sp = pkv_GetScratchMemTop ();
  if ( (c = pkv_GetScratchMemd ( G2H_OMCDEG+1 )) ) {
    hole_k = domain->hole_k;
    for ( i = 0; i < hole_k; i++ ) {
      _g2h_GetBFBPatchCurvesd ( domain, fn, i,
                                &c00, &c01, &c02, &c10, &c11, &c12,
                                &d00, &d01, &d02, &d10, &d11, &d12 );
      memcpy ( c, c00, (G2H_OMCDEG+1)*sizeof(double) );
      drawpoly ( G2H_OMCDEG, c );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g2h_DrawBFBomcd*/

void g2h_DrawFinalSurfBCd ( GHoleDomaind *domain,
                            int spdimen, const double *hole_cp,
                            const double *acoeff,
                            void (*drawcurve)(int degree, int spdimen,
                                              const double *cp) )
{
  void *sp;
  int  hole_k, nfunc_a, nfunc_b;
  GHolePrivateRecd  *privateG;
  G2HolePrivateRecd *privateG2;
  unsigned char *bfcpn;
  double *x, *y;
  double *c00, *c01, *c02, *c10, *c11, *c12, *d00, *d01, *d02, *d10, *d11, *d12;
  int   i, j;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG = domain->privateG;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  bfcpn = privateG->bfcpn;
  if ( (x = pkv_GetScratchMemd ( 2*spdimen*(G2H_OMCDEG+1) )) ) {
    y = &x[spdimen*(G2H_OMCDEG+1)];
    for ( i = 0; i < hole_k; i++ ) {
      memset ( x, 0, spdimen*(G2H_OMCDEG+1)*sizeof(double) );
      for ( j = 0; j < nfunc_a; j++ ) {
        _g2h_GetBFAPatchCurvesd ( domain, j, i, &c00, &c01, &c02,
                                  &d00, &d01, &d02 );
        pkn_MultMatrixd ( G2H_OMCDEG+1, 1, 1, c00, spdimen, 0,
                          &acoeff[j*spdimen], spdimen, y );
        pkn_AddMatrixd ( 1, spdimen*(G2H_OMCDEG+1), 0, x, 0, y, 0, x );
      }
      for ( j = 0; j < nfunc_b; j++ ) {
        _g2h_GetBFBPatchCurvesd ( domain, j, i,
                                  &c00, &c01, &c02, &c10, &c11, &c12,
                                  &d00, &d01, &d02, &d10, &d11, &d12 );
        pkn_MultMatrixd ( G2H_OMCDEG+1, 1, 1, c00, spdimen, 0,
                          &hole_cp[bfcpn[j]], spdimen, y );
        pkn_AddMatrixd ( 1, spdimen*(G2H_OMCDEG+1), 0, x, 0, y, 0, x );
      }
      drawcurve ( G2H_OMCDEG, spdimen, x );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g2h_DrawFinalSurfBCd*/

void g2h_ExtDrawFinalSurfBCd ( GHoleDomaind *domain,
                               int spdimen, const double *hole_cp,
                               const double *acoeff,
                               void (*drawcurve)(int degree, int spdimen,
                                                 const double *cp) )
{
  void *sp;
  int  hole_k, nfunc_a, nfunc_b, nfunc_c;
  GHolePrivateRecd  *privateG;
  G2HolePrivateRecd *privateG2;
  unsigned char *bfcpn;
  double *x, *y;
  double *c00, *c01, *c02, *c10, *c11, *c12, *d00, *d01, *d02, *d10, *d11, *d12;
  int   i, j;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG = domain->privateG;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  nfunc_c = hole_k*G2_DBDIM;
  bfcpn = privateG->bfcpn;
  if ( (x = pkv_GetScratchMemd ( 2*spdimen*(G2H_OMCDEG+1) )) ) {
    y = &x[spdimen*(G2H_OMCDEG+1)];
    for ( i = 0; i < hole_k; i++ ) {
      memset ( x, 0, spdimen*(G2H_OMCDEG+1)*sizeof(double) );
      for ( j = 0; j < nfunc_a; j++ ) {
        _g2h_GetBFAPatchCurvesd ( domain, j, i, &c00, &c01, &c02,
                                  &d00, &d01, &d02 );
        pkn_MultMatrixd ( G2H_OMCDEG+1, 1, 1, c00, spdimen, 0,
                          &acoeff[(nfunc_c+j)*spdimen], spdimen, y );
        pkn_AddMatrixd ( 1, spdimen*(G2H_OMCDEG+1), 0, x, 0, y, 0, x );
      }
      for ( j = 0; j < nfunc_b; j++ ) {
        _g2h_GetBFBPatchCurvesd ( domain, j, i,
                                  &c00, &c01, &c02, &c10, &c11, &c12,
                                  &d00, &d01, &d02, &d10, &d11, &d12 );
        pkn_MultMatrixd ( G2H_OMCDEG+1, 1, 1, c00, spdimen, 0,
                          &hole_cp[bfcpn[j]], spdimen, y );
        pkn_AddMatrixd ( 1, spdimen*(G2H_OMCDEG+1), 0, x, 0, y, 0, x );
      }
      drawcurve ( G2H_OMCDEG, spdimen, x );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g2h_ExtDrawFinalSurfBCd*/

/* ////////////////////////////////////////////////////////////////////////// */
void g2h_DrawMatricesd ( GHoleDomaind *domain,
                         void (*drawmatrix)(int nfa, int nfb,
                                            double *amat, double *bmat) )
{
  void              *sp;
  double            *amat, *bmat;
  int               na, nb;
  G2HolePrivateRecd *privateG2;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  if ( !(privateG2->Amat && privateG2->Bmat) )
    goto wayout;

  na = privateG2->nfunc_a;
  na = na*(na+1)/2;
  nb = privateG2->nfunc_a*privateG2->nfunc_b;
  amat = pkv_GetScratchMemd ( na );
  bmat = pkv_GetScratchMemd ( nb );
  if ( !amat || !bmat )
    goto wayout;

  memcpy ( amat, privateG2->Amat, na*sizeof(double) );
  memcpy ( bmat, privateG2->Bmat, nb*sizeof(double) );
  drawmatrix ( privateG2->nfunc_a, privateG2->nfunc_b, amat, bmat );

wayout:
  pkv_SetScratchMemTop ( sp );
} /*g2h_DrawMatricesd*/

/* ////////////////////////////////////////////////////////////////////////// */
void g2h_DrawExtMatricesd ( GHoleDomaind *domain,
                            void (*drawmatrix)(int k, int r, int s,
                                               double *Aii, double *Bi) )
{
  void   *sp;
  double *aii, *bi, *Aii, *Aki, *Akk, *Bi, *Bk, *Lii;
  G2HolePrivateRecd *privateG2;
  int    hole_k, nfunc_a, nfunc_b, nfunc_c, size;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  if ( !(privateG2->EAmat && privateG2->EBmat) )
    goto wayout;

  _g2h_GetExtBlockAddressesd ( domain, &Aii, &Aki, &Akk, &Bi, &Bk, &Lii );
  hole_k = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  nfunc_c = hole_k*G2_DBDIM;

  aii = pkv_GetScratchMemd ( (size = pkn_Block1ArraySize ( hole_k, G2_DBDIM, nfunc_a)) );
  bi  = pkv_GetScratchMemd ( (nfunc_c+nfunc_a)*nfunc_b );
  if ( aii && bi ) {
    memcpy ( aii, Aii, size*sizeof(double) );
    memcpy ( bi, Bi, (nfunc_c+nfunc_a)*nfunc_b*sizeof(double) );
    drawmatrix ( hole_k, G2_DBDIM, nfunc_a, aii, bi );
  }
wayout:
  pkv_SetScratchMemTop ( sp );
} /*g2h_DrawExtMatricesd*/

