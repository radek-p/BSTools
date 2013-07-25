
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#include <stdio.h>    
#include <stdlib.h>   
#include <math.h>
#include <malloc.h>  
#include <string.h>     

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "pknum.h"  
#include "pkgeom.h"   
#include "camerad.h"
#include "multibs.h"
#include "eg1holed.h"
#include "eg2holed.h"
#include "xgedit.h"   

#include "render.h"
#include "edpolicz.h"
#include "edwidgets.h"
#include "drawbezd.h"
#include "splhole.h"
#include "datagend.h"
  

/* ////////////////////////////////////////////////////////////////////////// */
boolean Option1SetOrder ( int k )
{
  if ( k != options1.order ) {
    switch ( k )
    {
  case 1:
      options1.order = 1;
      if ( options1.nk > G1H_S_MAX_NK )
        options1.nk = G1H_S_MAX_NK;
      if ( options1.m1 > G1H_S_MAX_M1 )
        options1.m1 = G1H_S_MAX_M1;
      if ( options1.m2 > G1H_S_MAX_M2 )
        options1.m2 = G1H_S_MAX_M2;
      if ( !domain1->basisG1 )
        g1h_ComputeBasisd ( domain1 );
      if ( domain_bcp1 )
        free ( domain_bcp1 );
      domain_bcp1 = NULL;
      g1h_DrawDiPatchesd ( domain1, OutputDomainPatch1 );
      break;
  case 2:
      options1.order = 2;
      options1.quasiG2 = false;
      if ( !domain1->basisG2 )
        g2h_ComputeBasisd ( domain1 );
      if ( domain_bcp1 )
        free ( domain_bcp1 );
      domain_bcp1 = NULL;
      g2h_DrawDiPatchesd ( domain1, OutputDomainPatch1 );
      break;
 default:
      return false;
    }
    options1.spl_basis_valid = false;
    InvalFinalSurface1 ();
    InitConstraintFrame ( 1 );
    view_surf_1 = false;
    swind_picture = false;
  }
  return true;
} /*Option1SetOrder*/

boolean Option1SetNK ( int k )
{
  if ( options1.spline && k >= 1 ) {
    switch ( options1.order ) {
  case 1:
      if ( k <= G1H_S_MAX_NK )
        options1.nk = k;
      break;
  case 2:
      if ( k <= G2H_S_MAX_NK )
        options1.nk = k;
      break;
 default:
      return false;
    }
    options1.spl_basis_valid = false;
    InvalFinalSurface1 ();
    InitConstraintFrame ( 1 );
    view_surf_1 = false;
    swind_picture = false;
  }
  return true;
} /*Option1SetNK*/

boolean Option1SetM1 ( int k )
{
  if ( options1.spline && k >= 1 ) {
    switch ( options1.order ) {
  case 1:
      if ( k <= G1H_S_MAX_M1 )
        options1.m1 = k;
      break;
  case 2:
      if ( k <= G2H_S_MAX_M1 )
        options1.m1 = k;
      break;
 default:
      return false;
    }
    options1.spl_basis_valid = false;
    InvalFinalSurface1 ();
    InitConstraintFrame ( 1 );
    view_surf_1 = false;
    swind_picture = false;
  }
  return true;
} /*Option1SetM1*/

boolean Option1SetM2 ( int k )
{
  if ( options1.spline && k >= 1 ) {
    switch ( options1.order ) {
  case 1:
      if ( k <= G1H_S_MAX_M2 )
        options1.m2 = k;
      break;
  case 2:
      if ( k <= G2H_S_MAX_M2 )
        options1.m2 = k;
      break;
 default:
      return false;
    }
    options1.spl_basis_valid = false;
    InvalFinalSurface1 ();
    view_surf_1 = false;
    swind_picture = false;
  }
  return true;
} /*Option1SetM2*/

boolean Option2SetOrder ( int k )
{
  if ( k != options2.order ) {
    switch ( k )
    {
  case 1:
      options2.order = 1;
      if ( options2.nk > G1H_S_MAX_NK )
        options2.nk = G1H_S_MAX_NK;
      if ( options2.m1 > G1H_S_MAX_M1 )
        options2.m1 = G1H_S_MAX_M1;
      if ( options2.m2 > G1H_S_MAX_M2 )
        options2.m2 = G1H_S_MAX_M2;
      if ( !domain2->basisG1 )
        g1h_ComputeBasisd ( domain2 );
      if ( domain_bcp2 )
        free ( domain_bcp2 );
      domain_bcp2 = NULL;
      g1h_DrawDiPatchesd ( domain1, OutputDomainPatch2 );
      break;
  case 2:
      options2.order = 2;
      options2.quasiG2 = false;
      if ( !domain2->basisG2 )
        g2h_ComputeBasisd ( domain2 );
      if ( domain_bcp2 )
        free ( domain_bcp2 );
      domain_bcp2 = NULL;
      g2h_DrawDiPatchesd ( domain1, OutputDomainPatch2 );
      break;
 default:
      return false;
    }
    options2.spl_basis_valid = false;
    InvalFinalSurface2 ();
    InitConstraintFrame ( 2 );
    view_surf_2 = false;
    swind_picture = false;
  }
  return true;
} /*Option2SetOrder*/

boolean Option2SetNK ( int k )
{
  if ( options2.spline && k >= 1 ) {
    switch ( options2.order ) {
  case 1:
      if ( k <= G1H_S_MAX_NK )
        options2.nk = k;
      break;
  case 2:
      if ( k <= G2H_S_MAX_NK )
        options2.nk = k;
      break;
 default:
      return false;
    }
    options2.spl_basis_valid = false;
    InvalFinalSurface2 ();
    InitConstraintFrame ( 2 );
    view_surf_2 = false;
    swind_picture = false;
  }
  return true;
} /*Option2SetNK*/

boolean Option2SetM1 ( int k )
{
  if ( options2.spline && k >= 1 ) {
    switch ( options2.order ) {
  case 1:
      if ( k <= G1H_S_MAX_M1 )
        options2.m1 = k;
      break;
  case 2:
      if ( k <= G2H_S_MAX_M1 )
        options2.m1 = k;
      break;
 default:
      return false;
    }
    options2.spl_basis_valid = false;
    InvalFinalSurface2 ();
    InitConstraintFrame ( 2 );
    view_surf_2 = false;
    swind_picture = false;
  }
  return true;
} /*Option2SetM1*/

boolean Option2SetM2 ( int k )
{
  if ( options2.spline && k >= 1 ) {
    switch ( options2.order ) {
  case 1:
      if ( k <= G1H_S_MAX_M2 )
        options2.m2 = k;
      break;
  case 2:
      if ( k <= G2H_S_MAX_M2 )
        options2.m2 = k;
      break;
 default:
      return false;
    }
    options2.spl_basis_valid = false;
    InvalFinalSurface2 ();
    view_surf_2 = false;
    swind_picture = false;
  }
  return true;
} /*Option2SetM2*/

/* ////////////////////////////////////////////////////////////////////////// */
int MyG1OptionProc1 ( GHoleDomaind *domain, int query, int qn,
                      int *ndata, int **idata, double **fdata )
{
  switch ( query ) {
case G1HQUERY_CENTRAL_POINT:
    if ( options1.altcentre )
      return G1H_CENTRAL_POINT_ALT;
    else
      return G1H_DEFAULT;
case G1HQUERY_BASIS:
    if ( options1.restricted )
      return G1H_USE_RESTRICTED_BASIS;
    else
      return G1H_DEFAULT;
case G1HQUERY_QUADRATURE:
    if ( options1.gausslegendre )
      return G1H_QUADRATURE_GAUSS_LEGENDRE;
    else
      return G1H_DEFAULT;
case G1HQUERY_Q2_FORM_CONSTANT:
    *ndata = 1;
    *fdata = &options1.q2coeff;
    return G1H_Q2_USE_SUPPLIED_CONSTANT;
default:
    return G1H_DEFAULT;
  }
} /*MyG1OptionProc1*/

int MyG2OptionProc1 ( GHoleDomaind *domain, int query, int qn,
                      int *ndata, int **idata, double **fdata )
{
  switch ( query ) {
case G2HQUERY_CENTRAL_POINT:
    if ( options1.altcentre )
      return G2H_CENTRAL_POINT_ALT;
    else
      return G2H_DEFAULT;
case G2HQUERY_BASIS:
    if ( options1.restricted )
      return G2H_USE_RESTRICTED_BASIS;
    else
      return G2H_DEFAULT;
case G2HQUERY_QUADRATURE:
    if ( options1.gausslegendre )
      return G2H_QUADRATURE_GAUSS_LEGENDRE;
    else
      return G2H_DEFAULT;
default:
    return G2H_DEFAULT;
  }
} /*MyG2OptionProc1*/

int MyG1OptionProc2 ( GHoleDomaind *domain, int query, int qn,
                      int *ndata, int **idata, double **fdata )
{
  switch ( query ) {
case G1HQUERY_CENTRAL_POINT:
    if ( options1.altcentre )
      return G1H_CENTRAL_POINT_ALT;
    else
      return G1H_DEFAULT;
case G1HQUERY_BASIS:
    if ( options2.restricted )
      return G1H_USE_RESTRICTED_BASIS;
    else
      return G1H_DEFAULT;
case G1HQUERY_QUADRATURE:
    if ( options2.gausslegendre )
      return G1H_QUADRATURE_GAUSS_LEGENDRE;
    else
      return G1H_DEFAULT;
case G1HQUERY_Q2_FORM_CONSTANT:
    *ndata = 1;
    *fdata = &options2.q2coeff;
    return G1H_Q2_USE_SUPPLIED_CONSTANT;
default:
    return G1H_DEFAULT;
  }
} /*MyG1OptionProc2*/

int MyG2OptionProc2 ( GHoleDomaind *domain, int query, int qn,
                      int *ndata, int **idata, double **fdata )
{
  switch ( query ) {
case G2HQUERY_CENTRAL_POINT:
    if ( options2.altcentre )
      return G2H_CENTRAL_POINT_ALT;
    else
      return G2H_DEFAULT;
case G2HQUERY_BASIS:
    if ( options2.restricted )
      return G2H_USE_RESTRICTED_BASIS;
    else
      return G2H_DEFAULT;
case G2HQUERY_QUADRATURE:
    if ( options2.gausslegendre )
      return G2H_QUADRATURE_GAUSS_LEGENDRE;
    else
      return G2H_DEFAULT;
default:
    return G2H_DEFAULT;
  }
} /*MyG2OptionProc2*/

/* ////////////////////////////////////////////////////////////////////////// */
void OutputDomainPatch1 ( int n, int m, const point2d *cp )
{
  int size;

  size = (n+1)*(m+1);
  if ( !domain_bcp1 ) {
    domain_bcp1 = malloc ( hole_k*size*sizeof(point2d) );
    domain_np1 = 0;
    domain_deg1 = n;
  }
  if ( domain_bcp1 && domain_np1 < hole_k ) {
    memcpy ( &domain_bcp1[size*domain_np1], cp, size*sizeof(point2d) );
    domain_np1 ++;
  }
} /*OutputDomainPatch1*/

void OutputDomainPatch2 ( int n, int m, const point2d *cp )
{
  int size;

  size = (n+1)*(m+1);
  if ( !domain_bcp2 ) {
    domain_bcp2 = malloc ( hole_k*size*sizeof(point2d) );
    domain_np2 = 0;
    domain_deg2 = n;
  }
  if ( domain_bcp2 && domain_np2 < hole_k ) {
    memcpy ( &domain_bcp2[size*domain_np2], cp, size*sizeof(point2d) );
    domain_np2 ++;
  }
} /*OutputDomainPatch1*/

boolean UpdateDomain1 ( void )
{
        /* destroy the old one */
  if ( domain1 ) {
    gh_DestroyDomaind ( domain1 );
    domain1 = NULL;
    options1.spl_basis_valid = false;
  }
  if ( domain_bcp1 )
    free ( domain_bcp1 );
  domain_bcp1 = NULL;
        /* create a new domain */
  domain1 = gh_CreateDomaind ( hole_k, knots, domain_cp );
  if ( domain1 ) {
    g1h_SetOptionProcd ( domain1, MyG1OptionProc1 );
    g2h_SetOptionProcd ( domain1, MyG2OptionProc1 );
    if ( options1.order == 1 ) {
      if ( g1h_ComputeBasisd ( domain1 ) )
        g1h_DrawDiPatchesd ( domain1, OutputDomainPatch1 );
    }
    else if ( options1.order == 2 ) {
      if ( g2h_ComputeBasisd ( domain1 ) )
        g2h_DrawDiPatchesd ( domain1, OutputDomainPatch1 );
    }
  }
  options1.constr_matrix_valid = false;
  InvalFinalSurface1 ();
  return domain1 != NULL;
} /*UpdateDomain1*/

boolean UpdateDomain2 ( void )
{
        /* destroy the old one */
  if ( domain2 ) {
    gh_DestroyDomaind ( domain2 );
    domain2 = NULL;
    options2.spl_basis_valid = false;
  }
  if ( domain_bcp2 )
    free ( domain_bcp2 );
  domain_bcp2 = NULL;
        /* create a new domain */
  domain2 = gh_CreateDomaind ( hole_k, knots, domain_cp );
  if ( domain2 ) {
    g1h_SetOptionProcd ( domain2, MyG1OptionProc2 );
    g2h_SetOptionProcd ( domain2, MyG2OptionProc2 );
    if ( options2.order == 1 ) {
      if ( g1h_ComputeBasisd ( domain2 ) )
        g1h_DrawDiPatchesd ( domain2, OutputDomainPatch2 );
    }
    else if ( options2.order == 2 ) {
      if ( g2h_ComputeBasisd ( domain2 ) )
        g2h_DrawDiPatchesd ( domain2, OutputDomainPatch2 );
    }
  }
  options2.constr_matrix_valid = false;
  InvalFinalSurface2 ();
  return domain2 != NULL;
} /*UpdateDomain2*/

boolean UpdateDomains ( void )
{
  UpdateDomain1 ();
  UpdateDomain2 ();
  return domain1 && domain2;
} /*UpdateDomains*/

boolean InitGHObject ( int k )
{
  if ( k >= 3 && k <= GH_MAX_K ) {
    hole_k = k;
    InitGHKnotsd ( hole_k, knots );
    InitGHVectors2d ( hole_k, domcvect );
    InitGHVectors3d ( hole_k, surfcvect );
    nctrlp = InitGHDomainNetd ( hole_k, domcvect, 2, domcparam, domain_cp );
    InitGHSurfNetd ( hole_k, surfcvect, 5, surfcparam, hole_cp );
    FindBoundingBox ( &swind.DefBBox );
    xge_3DwindInitProjections ( &swind, swind.DefBBox.x0, swind.DefBBox.x1,
        swind.DefBBox.y0, swind.DefBBox.y1, swind.DefBBox.z0, swind.DefBBox.z1 );
    FindDomainBoundingBox ( &domwind.DefBBox );
    xge_2DwindInitProjection ( &domwind, domwind.DefBBox.x0, domwind.DefBBox.x1,
        domwind.DefBBox.y0, domwind.DefBBox.y1 );
        /* clear control point selection */
    memset ( mkdcp, 0, nctrlp );
    memset ( mkhcp, 0, nctrlp );

    options1.bcmat = options2.bcmat = NULL;
    InvalFinalSurfaces ();
    ProjectSurfaceNet ();
    ProjectDomainNet ();
    if ( !UpdateDomains () )
      goto failure;
    InitConstraintFrame ( 1 );
    InitConstraintFrame ( 2 );
    return true;
  }
  else {
failure:
    return false;
  }
} /*InitGHObject*/

