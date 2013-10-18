
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

#undef CONST_
#define CONST_

#include "eg1holed.h"
#include "eg1hprivated.h"
#include "eg1herror.h"
#include "eg2holed.h"
#include "eg2hprivated.h"
#include "eg2herror.h"


int gh_GetDefaultOptiond ( GHoleDomaind *domain, int query, int qn,
                          int *ndata, int **idata, double **fdata )
{
  return G1H_DEFAULT;
} /*gh_GetDefaultOptiond*/

GHoleDomaind* gh_CreateDomaind ( int     hole_k,
                                 double  *hole_knots,         
                                 point2d *domain_cp )
{
  int               i;
  GHoleDomaind      *domain;
  GHolePrivateRecd  *privateG;
  G1HolePrivateRecd *privateG1;
  G2HolePrivateRecd *privateG2;
  double            *hkn;
  point2d           *domcp;

  domain = (GHoleDomaind*)malloc ( sizeof(GHoleDomaind) );
  if ( domain ) {
    domain->basisG1 = domain->basisG2 = false;
    privateG = (GHolePrivateRecd*)malloc ( sizeof(GHolePrivateRecd) );
    privateG1 = (G1HolePrivateRecd*)malloc ( sizeof(G1HolePrivateRecd) );
    privateG2 = (G2HolePrivateRecd*)malloc ( sizeof(G2HolePrivateRecd) );
    hkn = (double*)malloc ( hole_k*11*sizeof(double) );
    domcp = (point2d*)malloc ( (1+12*hole_k)*sizeof(point2d) );
    if ( privateG && privateG1 && privateG2 && hkn && domcp ) {
      privateG->dombezpt = privateG->surrpc = privateG->surrpcbc = NULL;
      privateG->domsurrbf = NULL;
      privateG->support_b = NULL;
      privateG->bfcpn = NULL;
      privateG->diam = 0.0;
      domain->privateG = privateG;

      privateG1->omcbc = privateG1->omc = privateG1->dicross = NULL;
      privateG1->jfunc = privateG1->basis_a = privateG1->basis_b =
      privateG1->AuxiMat = privateG1->partition = NULL;
      privateG1->spartition = NULL;
      privateG1->Amat = privateG1->Bmat = privateG1->BBmat = privateG1->Lmat = NULL;
      privateG1->EAmat = privateG1->EBmat = privateG1->ELmat = NULL;
      privateG1->Cmat = privateG1->RCmat =
      privateG1->ACmat = privateG1->ARCmat =
      privateG1->ECmat = privateG1->ERCmat =
      privateG1->AECmat = privateG1->AERCmat = NULL;
      privateG1->nconstr = privateG1->extnconstr = 0;
      privateG1->Q2AMat = privateG1->Q2BMat = privateG1->Q2LMat = NULL;
      privateG1->Q2RCMat = privateG1->Q2ARCMat = NULL;
      privateG1->Q2EAMat = privateG1->Q2EBMat = privateG1->Q2ELMat = NULL;
      privateG1->Q2ERCMat = privateG1->Q2AERCMat = NULL;
      domain->privateG1 = privateG1;
      domain->SprivateG1 = NULL;

      privateG2->omcbc = privateG2->omc = privateG2->dicross = NULL; 
      privateG2->jfunc = privateG2->basis_a = privateG2->basis_b =   
      privateG2->AuxiMat = privateG2->partition = NULL;
      privateG2->spartition = NULL;
      privateG2->Amat = privateG2->Bmat = privateG2->BBmat = privateG2->Lmat = NULL;
      privateG2->EAmat = privateG2->EBmat = privateG2->ELmat = NULL;
      privateG2->Cmat = privateG2->RCmat =
      privateG2->ACmat = privateG2->ARCmat =
      privateG2->ECmat = privateG2->ERCmat =
      privateG2->AECmat = privateG2->AERCmat = NULL;
      privateG2->nconstr = privateG2->extnconstr = 0;
      domain->privateG2 = privateG2;
      domain->SprivateG2 = NULL;

      domain->hole_k = hole_k;
      domain->hole_knots = hkn;
      domain->domain_cp = domcp;
      if ( hole_knots )
        memcpy ( hkn, hole_knots, hole_k*11*sizeof(double) );
      else { /* generate default equidistant knots */
        hkn[0] = hkn[1] = 0.0;
        for ( i = 1; i < 8; i++ ) hkn[i+1] = (double)i/8.0;
        hkn[9] = hkn[10] = 1.0;
        for ( i = 1; i < hole_k; i++ )
          memcpy ( &hkn[11*i], hkn, 11*sizeof(double) );
      }
      memcpy ( domcp, domain_cp, (1+12*hole_k)*sizeof(point2d) );
      domain->error_code = 0;
      privateG1->GetOption = gh_GetDefaultOptiond;
      privateG2->GetOption = gh_GetDefaultOptiond;
    }
    else {
      if ( privateG )  free ( privateG );
      if ( privateG1 ) free ( privateG1 );
      if ( privateG2 ) free ( privateG2 );
      if ( domcp ) free ( domcp );
      if ( hkn ) free ( hkn );
      free ( domain );
      domain = NULL;
    }
  }
  return domain;
} /*gh_CreateDomaind*/

void gh_DestroyDomaind ( GHoleDomaind *domain )
{
  GHolePrivateRecd  *privateG;
  G1HolePrivateRecd *privateG1;
  G2HolePrivateRecd *privateG2;

  privateG = domain->privateG;
  if ( privateG->dombezpt )   free ( privateG->dombezpt );
  if ( privateG->surrpc )     free ( privateG->surrpc );
  if ( privateG->surrpcbc )   free ( privateG->surrpcbc );
  if ( privateG->domsurrbf )  free ( privateG->domsurrbf );
  if ( privateG->bfcpn )      free ( privateG->bfcpn );
  if ( privateG->support_b )  free ( privateG->support_b );
  free ( privateG );

  g1h_DestroySPrivateDatad ( domain );
  g1h_DestroyQ2PrivateDatad ( domain );
  privateG1 = domain->privateG1;
  if ( privateG1->omcbc )      free ( privateG1->omcbc );
  if ( privateG1->omc )        free ( privateG1->omc );
  if ( privateG1->jfunc )      free ( privateG1->jfunc );
  if ( privateG1->dicross )    free ( privateG1->dicross );
  if ( privateG1->AuxiMat )    free ( privateG1->AuxiMat );
  if ( privateG1->basis_a )    free ( privateG1->basis_a );
  if ( privateG1->basis_b )    free ( privateG1->basis_b );
  if ( privateG1->partition )  free ( privateG1->partition );
  if ( privateG1->spartition ) free ( privateG1->spartition );
  if ( privateG1->Amat )       free ( privateG1->Amat );
  if ( privateG1->Bmat )       free ( privateG1->Bmat );
  if ( privateG1->Lmat )       free ( privateG1->Lmat );
  if ( privateG1->EAmat )      free ( privateG1->EAmat );
  if ( privateG1->EBmat )      free ( privateG1->EBmat );
  if ( privateG1->ELmat )      free ( privateG1->ELmat );
  if ( privateG1->Cmat )       free ( privateG1->Cmat );
  if ( privateG1->RCmat )      free ( privateG1->RCmat );
  if ( privateG1->ACmat )      free ( privateG1->ACmat );
  if ( privateG1->ARCmat )     free ( privateG1->ARCmat );
  if ( privateG1->ECmat )      free ( privateG1->ECmat );
  if ( privateG1->ERCmat )     free ( privateG1->ERCmat );
  if ( privateG1->AECmat )     free ( privateG1->AECmat );
  if ( privateG1->AERCmat )    free ( privateG1->AERCmat );
  if ( privateG1->BBmat )      free ( privateG1->BBmat );
  if ( privateG1->Q2RCMat )    free ( privateG1->Q2RCMat );
  if ( privateG1->Q2ERCMat )   free ( privateG1->Q2ERCMat );
  if ( privateG1->Q2ARCMat )   free ( privateG1->Q2ARCMat );
  if ( privateG1->Q2AERCMat )  free ( privateG1->Q2AERCMat );
  free ( privateG1 );

  g2h_DestroySPrivateDatad ( domain );
  privateG2 = domain->privateG2;
  if ( privateG2->omcbc )      free ( privateG2->omcbc );
  if ( privateG2->omc )        free ( privateG2->omc );
  if ( privateG2->jfunc )      free ( privateG2->jfunc );
  if ( privateG2->dicross )    free ( privateG2->dicross );
  if ( privateG2->AuxiMat )    free ( privateG2->AuxiMat );
  if ( privateG2->basis_a )    free ( privateG2->basis_a );
  if ( privateG2->basis_b )    free ( privateG2->basis_b );
  if ( privateG2->partition )  free ( privateG2->partition );
  if ( privateG2->spartition ) free ( privateG2->spartition );
  if ( privateG2->Amat )       free ( privateG2->Amat );
  if ( privateG2->Bmat )       free ( privateG2->Bmat );
  if ( privateG2->Lmat )       free ( privateG2->Lmat );
  if ( privateG2->EAmat )      free ( privateG2->EAmat );
  if ( privateG2->EBmat )      free ( privateG2->EBmat );
  if ( privateG2->ELmat )      free ( privateG2->ELmat );
  if ( privateG2->Cmat )       free ( privateG2->Cmat );
  if ( privateG2->RCmat )      free ( privateG2->RCmat );
  if ( privateG2->ACmat )      free ( privateG2->ACmat );
  if ( privateG2->ARCmat )     free ( privateG2->ARCmat );
  if ( privateG2->ECmat )      free ( privateG2->ECmat );
  if ( privateG2->ERCmat )     free ( privateG2->ERCmat );
  if ( privateG2->AECmat )     free ( privateG2->AECmat );
  if ( privateG2->AERCmat )    free ( privateG2->AERCmat );
  if ( privateG2->BBmat )      free ( privateG2->BBmat );
  free ( privateG2 );

  free ( domain->hole_knots );
  free ( domain->domain_cp );
  free ( domain );
} /*gh_DestroyDomaind*/

/* ///////////////////////////////////////////////////////////////////////// */
double *_gh_GetKnotSequenced ( GHoleDomaind *domain, int i )
{
  int hole_k;

  hole_k = domain->hole_k; 
  if ( i < 0 ) i += hole_k;
  else if ( i >= hole_k ) i -= hole_k;

  return &(domain->hole_knots[11*i]);
} /*_gh_GetKnotSequenced*/

point2d *_gh_GetDomSurrndPatchd ( GHoleDomaind *domain, int i, int j )
{
  int hole_k;
  GHolePrivateRecd *privateG;

  hole_k = domain->hole_k;
  if ( i < 0 ) i += hole_k;
  else if ( i >= hole_k ) i -= hole_k;

  privateG = domain->privateG;
  if ( !privateG )
    return NULL; 
  if ( !privateG->dombezpt )
    return NULL;
  return &(privateG->dombezpt[(2*i+j-1)*16]);
} /*_gh_GetDomSurrndPatchd*/

boolean _gh_FindDomSurrndPatchd ( GHoleDomaind *domain,
                                  int i, int j, point2d *bezcp )
{
  void    *sp;
  int     hole_k, k;
  point2d *domain_cp;
  int     *ind;
  point2d *q;
  double  *ukn, *vkn;

  sp = pkv_GetScratchMemTop ();
  ind = (int*)pkv_GetScratchMem ( 16*sizeof(int) );
  q = (point2d*)pkv_GetScratchMem ( 16*sizeof(point2d) );
  if ( !ind || !q ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  hole_k = domain->hole_k;
  domain_cp = domain->domain_cp;

  ukn = _gh_GetKnotSequenced ( domain, i-1 );  ukn += 3;
  vkn = _gh_GetKnotSequenced ( domain, i );    vkn += j;
  gh_GetBspInd ( hole_k, i, j, ind );
  for ( k = 0; k < 16; k++ )
    q[k] = domain_cp[ind[k]];
  if ( !mbs_BSPatchToBezd ( 2, 3, 7, ukn, 3, 7, vkn, 8, (double*)q,
                            NULL, NULL, NULL, NULL, NULL, NULL, 8, (double*)bezcp ) )
    goto failure;
  if ( j == 1 )
    mbs_multiReverseBSCurved ( 3, 0, NULL, 4, 2, 8, (double*)bezcp );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_gh_FindDomSurrndPatchd*/

boolean gh_FindDomSurrndBezPatchesd ( GHoleDomaind *domain )
{
  GHolePrivateRecd *privateG;
  int      hole_k;
  point2d  *dombezpt;
  vector2d *surrpc, *surrpcbc;
  int      i, j;

  hole_k = domain->hole_k;
  privateG = domain->privateG;
        /* at the second call return immediately */
  if ( privateG->dombezpt && privateG->surrpc && privateG->surrpcbc )
    return true;
        /* at the first call do the work */
  privateG->dombezpt = dombezpt =
    (point2d*)malloc ( hole_k*32*sizeof(point2d) );
  privateG->surrpc = surrpc =
    (vector2d*)malloc ( hole_k*24*sizeof(vector2d) );
  privateG->surrpcbc = surrpcbc =
    (vector2d*)malloc ( hole_k*36*sizeof(vector2d) );
  if ( !dombezpt || !surrpc || !surrpcbc ) {
    domain->error_code = G1H_ERROR_CANNOT_MALLOC;
    return false;
  }

  for ( i = 0; i < hole_k; i++ )
    for ( j = 0; j < 2; j++ ) {
      if ( !_gh_FindDomSurrndPatchd ( domain, i, j+1, dombezpt ) )
        return false;
      if ( !mbs_multiBCHornerDer2d ( 3, 1, 8, 0, (double*)dombezpt, 0.0,
              (double*)surrpc, (double*)&surrpc[4], (double*)&surrpc[8] ) )
        return false;
      if ( !mbs_multiBCHornerDer2d ( 3, 3, 2, 8, (double*)surrpc, 0.0,
              (double*)surrpcbc, (double*)&surrpcbc[3], (double*)&surrpcbc[6] ) )
        return false;
      if ( !mbs_multiBCHornerDer2d ( 3, 3, 2, 8, (double*)surrpc, 1.0,
              (double*)&surrpcbc[9], (double*)&surrpcbc[12], (double*)&surrpcbc[15] ) )
        return false;
      dombezpt += 16;
      surrpc += 12;
      surrpcbc += 18;
    }

  return true;
} /*gh_FindDomSurrndBezPatchesd*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean _gh_FindDomSurrndBFuncPatchesd ( GHoleDomaind *domain )
{
  void             *sp;
  GHolePrivateRecd *privateG;
  int              hole_k, k, i, j, fn, nfunc_b;
  int              *ind;
  double           *q, *bfc, *domsurrbf, *fcp;
  unsigned char    *bfcpn;
  unsigned short   *support_b;
  double           *ukn, *vkn;

  sp  = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG = (GHolePrivateRecd*)domain->privateG;
        /* if it has been called before, return at once */
  if ( privateG->domsurrbf && privateG->bfcpn && privateG->support_b )
    return true;

  nfunc_b = 6*hole_k+1;
  bfcpn = privateG->bfcpn = malloc ( (12*hole_k+1)*sizeof(unsigned char) );
  domsurrbf = privateG->domsurrbf = malloc ( nfunc_b*hole_k*(16*3*sizeof(double)) );
  support_b = privateG->support_b = malloc ( nfunc_b*sizeof(unsigned short) );
  if ( !bfcpn || !domsurrbf || !support_b )
    goto failure;

  if ( !_gh_FindBasisFuncbSupport ( hole_k, nfunc_b, support_b, bfcpn ) )
    goto failure;

  ind = pkv_GetScratchMem ( 16*sizeof(int) );
  q   = pkv_GetScratchMemd ( 16 );
  bfc = pkv_GetScratchMemd ( 12*hole_k+1 );
  if ( !ind || !q || !bfc )
    goto failure;

  memset ( bfc, 0, (12*hole_k+1)*sizeof(double) );
  for ( fn = 0; fn < nfunc_b; fn++ ) {
    bfc[bfcpn[fn]] = 1.0;
    for ( i = 0; i < hole_k; i++ ) {
      ukn = &domain->hole_knots[11*((i+hole_k-1) % hole_k)+3];
      for ( j = 0; j < 3; j++ ) {
        vkn = &domain->hole_knots[11*i+j];
        gh_GetBspInd ( hole_k, i, j, ind );
        for ( k = 0; k < 16; k++ )
          q[k] = bfc[ind[k]];
        fcp = &domsurrbf[((fn*hole_k+i)*3+j)*16];
        if ( !mbs_BSPatchToBezd ( 1, 3, 7, ukn, 3, 7, vkn, 4, q,
                                  NULL, NULL, NULL, NULL, NULL, NULL, 4, fcp ) )
          goto failure;
        if ( j == 1 )
          mbs_multiReverseBSCurved ( 3, 0, NULL, 4, 1, 4, fcp );
      }
    }
    bfc[bfcpn[fn]] = 0.0;
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  if ( bfcpn )     { free ( bfcpn );      privateG->bfcpn = NULL; }
  if ( domsurrbf ) { free ( domsurrbf );  privateG->domsurrbf = NULL; }
  if ( support_b ) { free ( support_b );  privateG->support_b = NULL; }
  return false;
} /*_gh_FindDomSurrndBFuncPatchesd*/

void gh_GetDomSurrndBFuncd ( GHoleDomaind *domain, int fn, int i, int j,
                             double *bf )
{
  GHolePrivateRecd *privateG;
  int              hole_k;
  
  hole_k = domain->hole_k;
  privateG = domain->privateG;
  memcpy ( bf, &privateG->domsurrbf[((fn*hole_k+i)*3+j)*16], 16*sizeof(double) );
} /*gh_GetDomSurrndBFuncd*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean _gh_AnalyzePartitiond ( GHoleDomaind *domain, int omcdeg,
                                vector2d *omcbc, int *_hole_m, int *_auxp_k,
                                double **_partition, GHoleSgnPartd **_spartition,
                                double *_spart_alpha0 )
{
  int            hole_k, i, auxp_k;
  double         alpha, alpha0, delta;
  double         *partition, *part_delta;
  GHoleSgnPartd  *spartition;

  hole_k = domain->hole_k;
  *_partition = partition = (double*)malloc ( 2*hole_k*sizeof(double) );
  *_spartition = spartition =
      (GHoleSgnPartd*)malloc ( hole_k*sizeof(GHoleSgnPartd) );
  if ( !partition || !spartition ) {
    domain->error_code = G2H_ERROR_CANNOT_MALLOC;
    return false;
  }
  memset ( spartition, 0, hole_k*sizeof(GHoleSgnPartd) );
  part_delta = &partition[hole_k];

  for ( i = 0; i < hole_k; i++ ) {
    alpha = atan2 ( omcbc[(omcdeg+1)*i+1].y, omcbc[(omcdeg+1)*i+1].x );
    if ( alpha < 0.0 )
      alpha += 2.0*PI;
    spartition[i].alpha = partition[i] = (double)alpha;
  }
  alpha0 = partition[0];
  for ( i = 0; i < hole_k; i++ ) {
    alpha = partition[i];
    delta = partition[(i+1) % hole_k]-alpha;
    if ( delta < 0.0 )         delta += 2.0*PI;
    else if ( delta > 2.0*PI ) delta -= 2.0*PI;
    if ( delta <= MIN_DELTA || delta >= MAX_DELTA ) {
      domain->error_code = G2H_ERROR_INVALID_PARTITION;
      return false;
    }
    part_delta[i] = (double)delta;
    delta = alpha-alpha0;
    if ( delta < -ANGLETOL )
      delta += 2.0*PI;
    else if ( delta < 0.0 )
      delta = 0.0;
    if ( delta < PI-ANGLETOL ) {
      spartition[i].malpha = (double)alpha;
      spartition[i].sgn = false;
    }
    else {
      if ( alpha >= PI )
        spartition[i].malpha = (double)(alpha-PI);
      else
        spartition[i].malpha = (double)(alpha+PI);
      spartition[i].sgn = true;
    }
    alpha = spartition[i].malpha-alpha0;
    if ( alpha < -ANGLETOL )
      alpha += 2.0*PI;
    else if ( alpha < 0.0 )
      alpha = 0.0;
    spartition[i].salpha = (double)alpha;
    spartition[i].both = false;
  }
       /* sort by the values of the salpha fields - the third parameter */
  if ( pkv_SortFast ( sizeof(double), ID_IEEE754_DOUBLE, sizeof(GHoleSgnPartd),
                      2*sizeof(double), hole_k, spartition ) != SORT_OK ) {
    domain->error_code = G2H_ERROR_CANNOT_MALLOC;
    return false;
  }

  auxp_k = 0;
  for ( i = 0; i < hole_k; i++ ) {
    delta = spartition[i].malpha-spartition[auxp_k].malpha;
    if ( fabs(delta-2.0*PI) < ANGLETOL )
      delta = 0.0;
    if ( delta < -ANGLETOL )
      delta += 2.0*PI;
    if ( delta >= MIN_DELTA ) {
      auxp_k++;
      spartition[auxp_k] = spartition[i];
    }
    else if ( delta < ANGLETOL &&
              spartition[i].sgn != spartition[auxp_k].sgn ) {
      spartition[auxp_k].both = true;
    }
    else if ( i != auxp_k ) {
      domain->error_code = G2H_ERROR_INVALID_PARTITION;
      return false;
    }
  }
  auxp_k++;
  *_hole_m = hole_k-auxp_k;

  if ( spartition[auxp_k-1].malpha > spartition[0].malpha )
    alpha0 = 0.5*(spartition[0].malpha+spartition[auxp_k-1].malpha);
  else
    alpha0 = 0.5*((spartition[0].malpha-2.0*PI)+spartition[auxp_k-1].malpha);
  *_spart_alpha0 = (double)alpha0;
  for ( i = 0; i < auxp_k; i++ ) {
    alpha = alpha0-spartition[i].malpha;
    spartition[i].knot = tan(alpha);
  }
  *_auxp_k = auxp_k;
  return true;
} /*_gh_AnalyzePartitiond*/

/* ///////////////////////////////////////////////////////////////////////// */
void _gh_PrepareTabKnotsd ( int nquad, int opt, double *knots )
{
  double a, b, h;
  int    i;

  switch ( opt ) {
case G1H_QUADRATURE_GAUSS_LEGENDRE:
         /* composite Gauss-Legendre of order 4 formula, nquad must be even */
    nquad /= 2;
    h = 1.0/nquad;
    a = (1.0-1.0/SQRT3)/2.0;
    b = (1.0+1.0/SQRT3)/2.0;
    for ( i = 0; i < nquad; i++ ) {
      knots[2*i] = h*((double)i+a);
      knots[2*i+1] = h*((double)i+b);
    }
    break;
default:      /* composite rectangles formula */
    for ( i = 0; i < nquad; i++ )
      knots[i] = ((double)(i+i+1))/((double)(2*nquad));
    break;
  }
} /*_gh_PrepareTabKnotsd*/

boolean _g2h_DiJacobian3d ( const vector2d *du, const vector2d *dv,
                            const vector2d *duu, const vector2d *duv,
                            const vector2d *dvv,
                            const vector2d *duuu, const vector2d *duuv,
                            const vector2d *duvv, const vector2d *dvvv,
                            double *jac, double *trd )
{
  vector2d gx, gy, gxx, gxy, gyy, gxxx, gxxy, gxyy, gyyy;
  double   A31[8], A32[12], A33[16];

  *jac = (double)fabs ( du->x*dv->y - du->y*dv->x );

  if ( !pkn_f2iDerivatives3d ( du->x, du->y, dv->x, dv->y,
          duu->x, duu->y, duv->x, duv->y, dvv->x, dvv->y,
          duuu->x, duuu->y, duuv->x, duuv->y, duvv->x, duvv->y, dvvv->x, dvvv->y,
          (double*)&gx, (double*)&gy, (double*)&gxx, (double*)&gxy, (double*)&gyy,
          (double*)&gxxx, (double*)&gxxy, (double*)&gxyy, (double*)&gyyy ) )
    return false;

  pkn_Setup2DerA31Matrixd ( gxxx.x, gxxx.y, gxxy.x, gxxy.y,
                            gxyy.x, gxyy.y, gyyy.x, gyyy.y, A31 );
  pkn_Setup2DerA32Matrixd ( gx.x, gx.y, gy.x, gy.y,
                            gxx.x, gxx.y, gxy.x, gxy.y, gyy.x, gyy.y, A32 );
  pkn_Setup2DerA33Matrixd ( gx.x, gx.y, gy.x, gy.y, A33 );

  trd[0]  = A31[0]+A31[4];   trd[1]  = A31[1]+A31[5];
  trd[2]  = A32[0]+A32[6];   trd[3]  = A32[1]+A32[7];   trd[4] = A32[2]+A32[8];
  trd[5]  = A33[0]+A33[8];   trd[6]  = A33[1]+A33[9];
  trd[7]  = A33[2]+A33[10];  trd[8]  = A33[3]+A33[11];

  trd[9]  = A31[2]+A31[6];   trd[10] = A31[3]+A31[7];
  trd[11] = A32[3]+A32[9];   trd[12] = A32[4]+A32[10];  trd[13] = A32[5]+A32[11];
  trd[14] = A33[4]+A33[12];  trd[15] = A33[5]+A33[13];
  trd[16] = A33[6]+A33[14];  trd[17] = A33[7]+A33[15];
  return true;
} /*_g2h_DiJacobian3d*/

double _g2h_Integrald ( int hole_k, int nknsq, double *jac,
                        unsigned short supp1, vector2d *func1,
                        unsigned short supp2, vector2d *func2 )
{
  int            i, k;
  long double    s;
  unsigned short supp;
  vector2d       *f1, *f2;
  double         *jc;

  supp = (unsigned short)(supp1 & supp2);
  s = 0.0;
  for ( k = 0; k < hole_k; k++ )
    if ( supp & (0x0001 << k) ) {
      jc = &jac[k*nknsq];
      f1 = &func1[k*nknsq];
      f2 = &func2[k*nknsq];
      for ( i = 0; i < nknsq; i++ )
        s += (f1[i].x*f2[i].x + f1[i].y*f2[i].y)*jc[i];
    }
  s /= (double)(nknsq);
  return (double)s;
} /*_g2h_Integrald*/

