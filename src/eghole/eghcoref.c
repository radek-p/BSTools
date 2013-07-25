
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

#include "eg1holef.h"
#include "eg1hprivatef.h"
#include "eg1herror.h"
#include "eg2holef.h"
#include "eg2hprivatef.h"
#include "eg2herror.h"


int gh_GetDefaultOptionf ( GHoleDomainf *domain, int query, int qn,
                           int *ndata, int **idata, float **fdata )
{
  return G1H_DEFAULT;
} /*gh_GetDefaultOptionf*/

GHoleDomainf* gh_CreateDomainf ( int     hole_k,
                                 float   *hole_knots,         
                                 point2f *domain_cp )
{
  int               i;
  GHoleDomainf      *domain;
  GHolePrivateRecf  *privateG;
  G1HolePrivateRecf *privateG1;
  G2HolePrivateRecf *privateG2;
  float             *hkn;
  point2f           *domcp;

  domain = (GHoleDomainf*)malloc ( sizeof(GHoleDomainf) );
  if ( domain ) {
    domain->basisG1 = domain->basisG2 = false;
    privateG = (GHolePrivateRecf*)malloc ( sizeof(GHolePrivateRecf) );
    privateG1 = (G1HolePrivateRecf*)malloc ( sizeof(G1HolePrivateRecf) );
    privateG2 = (G2HolePrivateRecf*)malloc ( sizeof(G2HolePrivateRecf) );
    hkn = (float*)malloc ( hole_k*11*sizeof(float) );
    domcp = (point2f*)malloc ( (1+12*hole_k)*sizeof(point2f) );
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
        memcpy ( hkn, hole_knots, hole_k*11*sizeof(float) );
      else { /* generate default equidistant knots */
        hkn[0] = hkn[1] = 0.0;
        for ( i = 1; i < 8; i++ ) hkn[i+1] = (float)i/8.0;
        hkn[9] = hkn[10] = 1.0;
        for ( i = 1; i < hole_k; i++ )
          memcpy ( &hkn[11*i], hkn, 11*sizeof(float) );
      }
      memcpy ( domcp, domain_cp, (1+12*hole_k)*sizeof(point2f) );
      domain->error_code = 0;
      privateG1->GetOption = gh_GetDefaultOptionf;
      privateG2->GetOption = gh_GetDefaultOptionf;
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
} /*gh_CreateDomainf*/

void gh_DestroyDomainf ( GHoleDomainf *domain )
{
  GHolePrivateRecf  *privateG;
  G1HolePrivateRecf *privateG1;
  G2HolePrivateRecf *privateG2;

  privateG = domain->privateG;
  if ( privateG->dombezpt )   free ( privateG->dombezpt );
  if ( privateG->surrpc )     free ( privateG->surrpc );
  if ( privateG->surrpcbc )   free ( privateG->surrpcbc );
  if ( privateG->domsurrbf )  free ( privateG->domsurrbf );
  if ( privateG->bfcpn )      free ( privateG->bfcpn );
  if ( privateG->support_b )  free ( privateG->support_b );
  free ( privateG );

  g1h_DestroySPrivateDataf ( domain );
  g1h_DestroyQ2PrivateDataf ( domain );
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

  g2h_DestroySPrivateDataf ( domain );
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
} /*gh_DestroyDomainf*/

/* ///////////////////////////////////////////////////////////////////////// */
float *_gh_GetKnotSequencef ( GHoleDomainf *domain, int i )
{
  int hole_k;

  hole_k = domain->hole_k; 
  if ( i < 0 ) i += hole_k;
  else if ( i >= hole_k ) i -= hole_k;

  return &(domain->hole_knots[11*i]);
} /*_gh_GetKnotSequencef*/

point2f *_gh_GetDomSurrndPatchf ( GHoleDomainf *domain, int i, int j )
{
  int hole_k;
  GHolePrivateRecf *privateG;

  hole_k = domain->hole_k;
  if ( i < 0 ) i += hole_k;
  else if ( i >= hole_k ) i -= hole_k;

  privateG = domain->privateG;
  if ( !privateG )
    return NULL; 
  if ( !privateG->dombezpt )
    return NULL;
  return &(privateG->dombezpt[(2*i+j-1)*16]);
} /*_gh_GetDomSurrndPatchf*/

boolean _gh_FindDomSurrndPatchf ( GHoleDomainf *domain,
                                  int i, int j, point2f *bezcp )
{
  void    *sp;
  int     hole_k, k;
  point2f *domain_cp;
  int     *ind;
  point2f *q;
  float   *ukn, *vkn;

  sp = pkv_GetScratchMemTop ();
  ind = (int*)pkv_GetScratchMem ( 16*sizeof(int) );
  q = (point2f*)pkv_GetScratchMem ( 16*sizeof(point2f) );
  if ( !ind || !q ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  hole_k = domain->hole_k;
  domain_cp = domain->domain_cp;

  ukn = _gh_GetKnotSequencef ( domain, i-1 );  ukn += 3;
  vkn = _gh_GetKnotSequencef ( domain, i );    vkn += j;
  gh_GetBspInd ( hole_k, i, j, ind );
  for ( k = 0; k < 16; k++ )
    q[k] = domain_cp[ind[k]];
  mbs_BSPatchToBezf ( 2, 3, 7, ukn, 3, 7, vkn, 8, (float*)q,
                      NULL, NULL, NULL, NULL, NULL, NULL, 8, (float*)bezcp );
  if ( j == 1 )
    mbs_multiReverseBSCurvef ( 3, 0, NULL, 4, 2, 8, (float*)bezcp );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_gh_FindDomSurrndPatchf*/

boolean gh_FindDomSurrndBezPatchesf ( GHoleDomainf *domain )
{
  GHolePrivateRecf *privateG;
  int      hole_k;
  point2f  *dombezpt;
  vector2f *surrpc, *surrpcbc;
  int      i, j;

  hole_k = domain->hole_k;
  privateG = domain->privateG;
        /* at the second call return immediately */
  if ( privateG->dombezpt && privateG->surrpc && privateG->surrpcbc )
    return true;
        /* at the first call do the work */
  privateG->dombezpt = dombezpt =
    (point2f*)malloc ( hole_k*32*sizeof(point2f) );
  privateG->surrpc = surrpc =
    (vector2f*)malloc ( hole_k*24*sizeof(vector2f) );
  privateG->surrpcbc = surrpcbc =
    (vector2f*)malloc ( hole_k*36*sizeof(vector2f) );
  if ( !dombezpt || !surrpc || !surrpcbc ) {
    domain->error_code = G1H_ERROR_CANNOT_MALLOC;
    return false;
  }

  for ( i = 0; i < hole_k; i++ )
    for ( j = 0; j < 2; j++ ) {
      if ( !_gh_FindDomSurrndPatchf ( domain, i, j+1, dombezpt ) )
        return false;
      mbs_multiBCHornerDer2f ( 3, 1, 8, 0, (float*)dombezpt, 0.0,
          (float*)surrpc, (float*)&surrpc[4], (float*)&surrpc[8] );
      mbs_multiBCHornerDer2f ( 3, 3, 2, 8, (float*)surrpc, 0.0,
          (float*)surrpcbc, (float*)&surrpcbc[3], (float*)&surrpcbc[6] );
      mbs_multiBCHornerDer2f ( 3, 3, 2, 8, (float*)surrpc, 1.0,
          (float*)&surrpcbc[9], (float*)&surrpcbc[12], (float*)&surrpcbc[15] );
      dombezpt += 16;
      surrpc += 12;
      surrpcbc += 18;
    }

  return true;
} /*gh_FindDomSurrndBezPatchesf*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean _gh_FindDomSurrndBFuncPatchesf ( GHoleDomainf *domain )
{
  void             *sp;
  GHolePrivateRecf *privateG;
  int              hole_k, k, i, j, fn, nfunc_b;
  int              *ind;
  float            *q, *bfc, *domsurrbf, *fcp;
  unsigned char    *bfcpn;
  unsigned short   *support_b;
  float            *ukn, *vkn;

  sp  = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG = (GHolePrivateRecf*)domain->privateG;
        /* if it has been called before, return at once */
  if ( privateG->domsurrbf && privateG->bfcpn && privateG->support_b )
    return true;

  nfunc_b = 6*hole_k+1;
  bfcpn = privateG->bfcpn = malloc ( (12*hole_k+1)*sizeof(unsigned char) );
  domsurrbf = privateG->domsurrbf = malloc ( nfunc_b*hole_k*(16*3*sizeof(float)) );
  support_b = privateG->support_b = malloc ( nfunc_b*sizeof(unsigned short) );
  if ( !bfcpn || !domsurrbf || !support_b )
    goto failure;

  if ( !_gh_FindBasisFuncbSupport ( hole_k, nfunc_b, support_b, bfcpn ) )
    goto failure;

  ind = pkv_GetScratchMem ( 16*sizeof(int) );
  q   = pkv_GetScratchMemf ( 16 );
  bfc = pkv_GetScratchMemf ( 12*hole_k+1 );
  if ( !ind || !q || !bfc )
    goto failure;

  memset ( bfc, 0, (12*hole_k+1)*sizeof(float) );
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
        mbs_BSPatchToBezf ( 1, 3, 7, ukn, 3, 7, vkn, 4, q,
                            NULL, NULL, NULL, NULL, NULL, NULL, 4, fcp );
        if ( j == 1 )
          mbs_multiReverseBSCurvef ( 3, 0, NULL, 4, 1, 4, fcp );
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
} /*_gh_FindDomSurrndBFuncPatchesf*/

void gh_GetDomSurrndBFuncf ( GHoleDomainf *domain, int fn, int i, int j,
                             float *bf )
{
  GHolePrivateRecf *privateG;
  int              hole_k;
  
  hole_k = domain->hole_k;
  privateG = domain->privateG;
  memcpy ( bf, &privateG->domsurrbf[((fn*hole_k+i)*3+j)*16], 16*sizeof(float) );
} /*gh_GetDomSurrndBFuncf*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean _gh_AnalyzePartitionf ( GHoleDomainf *domain, int omcdeg,
                                vector2f *omcbc, int *_hole_m, int *_auxp_k,
                                float **_partition, GHoleSgnPartf **_spartition,
                                float *_spart_alpha0 )
{
  int            hole_k, i, auxp_k;
  double         alpha, alpha0, delta;
  float          *partition, *part_delta;
  GHoleSgnPartf  *spartition;

  hole_k = domain->hole_k;
  *_partition = partition = (float*)malloc ( 2*hole_k*sizeof(float) );
  *_spartition = spartition =
      (GHoleSgnPartf*)malloc ( hole_k*sizeof(GHoleSgnPartf) );
  if ( !partition || !spartition ) {
    domain->error_code = G2H_ERROR_CANNOT_MALLOC;
    return false;
  }
  memset ( spartition, 0, hole_k*sizeof(GHoleSgnPartf) );
  part_delta = &partition[hole_k];

  for ( i = 0; i < hole_k; i++ ) {
    alpha = atan2 ( omcbc[(omcdeg+1)*i+1].y, omcbc[(omcdeg+1)*i+1].x );
    if ( alpha < 0.0 )
      alpha += 2.0*PI;
    spartition[i].alpha = partition[i] = (float)alpha;
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
    part_delta[i] = (float)delta;
    delta = alpha-alpha0;
    if ( delta < -ANGLETOL )
      delta += 2.0*PI;
    else if ( delta < 0.0 )
      delta = 0.0;
    if ( delta < PI-ANGLETOL ) {
      spartition[i].malpha = (float)alpha;
      spartition[i].sgn = false;
    }
    else {
      if ( alpha >= PI )
        spartition[i].malpha = (float)(alpha-PI);
      else
        spartition[i].malpha = (float)(alpha+PI);
      spartition[i].sgn = true;
    }
    alpha = spartition[i].malpha-alpha0;
    if ( alpha < -ANGLETOL )
      alpha += 2.0*PI;
    else if ( alpha < 0.0 )
      alpha = 0.0;
    spartition[i].salpha = (float)alpha;
    spartition[i].both = false;
  }
       /* sort by the values of the salpha fields - the third parameter */
  if ( pkv_SortFast ( sizeof(float), ID_IEEE754_FLOAT, sizeof(GHoleSgnPartf),
                      2*sizeof(float), hole_k, spartition ) != SORT_OK ) {
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
  *_spart_alpha0 = (float)alpha0;
  for ( i = 0; i < auxp_k; i++ ) {
    alpha = alpha0-spartition[i].malpha;
    spartition[i].knot = (float)tan(alpha);
  }
  *_auxp_k = auxp_k;
  return true;
} /*_gh_AnalyzePartitionf*/

/* ///////////////////////////////////////////////////////////////////////// */
void _gh_PrepareTabKnotsf ( int nquad, int opt, float *knots )
{
  float a, b, h;
  int   i;

  switch ( opt ) {
case G1H_QUADRATURE_GAUSS_LEGENDRE:
         /* composite Gauss-Legendre of order 4 formula, nquad must be even */
    nquad /= 2;
    h = (float)(1.0/nquad);
    a = (float)((1.0-1.0/SQRT3)/2.0);
    b = (float)((1.0+1.0/SQRT3)/2.0);
    for ( i = 0; i < nquad; i++ ) {
      knots[2*i] = h*((float)i+a);
      knots[2*i+1] = h*((float)i+b);
    }
    break;
default:      /* composite rectangles formula */
    for ( i = 0; i < nquad; i++ )
      knots[i] = ((float)(i+i+1))/((float)(2*nquad));
    break;
  }
} /*_gh_PrepareTabKnotsf*/

void _g2h_DiJacobian3f ( const vector2f *du, const vector2f *dv,
                         const vector2f *duu, const vector2f *duv,
                         const vector2f *dvv,
                         const vector2f *duuu, const vector2f *duuv,
                         const vector2f *duvv, const vector2f *dvvv,
                         float *jac, float *trd )
{
  vector2f gx, gy, gxx, gxy, gyy, gxxx, gxxy, gxyy, gyyy;
  float    A31[8], A32[12], A33[16];

  *jac = (float)fabs ( du->x*dv->y - du->y*dv->x );

  pkn_f2iDerivatives3f ( du->x, du->y, dv->x, dv->y,
      duu->x, duu->y, duv->x, duv->y, dvv->x, dvv->y,
      duuu->x, duuu->y, duuv->x, duuv->y, duvv->x, duvv->y, dvvv->x, dvvv->y,
      (float*)&gx, (float*)&gy, (float*)&gxx, (float*)&gxy, (float*)&gyy,
      (float*)&gxxx, (float*)&gxxy, (float*)&gxyy, (float*)&gyyy );

  pkn_Setup2DerA31Matrixf ( gxxx.x, gxxx.y, gxxy.x, gxxy.y,
                            gxyy.x, gxyy.y, gyyy.x, gyyy.y, A31 );
  pkn_Setup2DerA32Matrixf ( gx.x, gx.y, gy.x, gy.y,
                            gxx.x, gxx.y, gxy.x, gxy.y, gyy.x, gyy.y, A32 );
  pkn_Setup2DerA33Matrixf ( gx.x, gx.y, gy.x, gy.y, A33 );

  trd[0]  = A31[0]+A31[4];   trd[1]  = A31[1]+A31[5];
  trd[2]  = A32[0]+A32[6];   trd[3]  = A32[1]+A32[7];   trd[4] = A32[2]+A32[8];
  trd[5]  = A33[0]+A33[8];   trd[6]  = A33[1]+A33[9];
  trd[7]  = A33[2]+A33[10];  trd[8]  = A33[3]+A33[11];

  trd[9]  = A31[2]+A31[6];   trd[10] = A31[3]+A31[7];
  trd[11] = A32[3]+A32[9];   trd[12] = A32[4]+A32[10];  trd[13] = A32[5]+A32[11];
  trd[14] = A33[4]+A33[12];  trd[15] = A33[5]+A33[13];
  trd[16] = A33[6]+A33[14];  trd[17] = A33[7]+A33[15];
} /*_g2h_DiJacobian3f*/

float _g2h_Integralf ( int hole_k, int nknsq, float *jac,
                       unsigned short supp1, vector2f *func1,
                       unsigned short supp2, vector2f *func2 )
{
  int            i, k;
  double         s;
  unsigned short supp;
  vector2f       *f1, *f2;
  float          *jc;

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
  return (float)s;
} /*_g2h_Integralf*/

