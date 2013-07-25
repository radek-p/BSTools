
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <unistd.h>
#include <sys/times.h>

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

#define SHOW_TIME

/* ////////////////////////////////////////////////////////////////////////// */
void InvalFinalSurface1 ( void )
{
  if ( acoeff1 )      { free ( acoeff1 );       acoeff1 = NULL; }
  if ( final_cp1 )    { free ( final_cp1 );     final_cp1 = NULL; }
  if ( final_knots1 ) { free ( final_knots1 );  final_knots1 = NULL; }
  final_np1 = 0;
  options1.surf_valid = false;
} /*InvalFinalSurface1*/

void InvalFinalSurface2 ( void )
{
  if ( acoeff2 )      { free ( acoeff2 );       acoeff2 = NULL; }
  if ( final_cp2 )    { free ( final_cp2 );     final_cp2 = NULL; }
  if ( final_knots2 ) { free ( final_knots2 );  final_knots2 = NULL; }
  final_np2 = 0;
  options2.surf_valid = false;
} /*InvalFinalSurface2*/

void InvalFinalSurfaces ( void )
{
  InvalFinalSurface1 ();
  InvalFinalSurface2 ();
} /*InvalFinalSurfaces*/

void OutputFinalBezPatch1 ( int n, int m, const double *cp, void *usrptr )
{
  int size;

  size = (n+1)*(m+1);
  if ( !final_cp1 ) {
    final_cp1 = malloc ( hole_k*size*sizeof(point3d) );
    final_np1 = 0;
    final_deg1 = n;
  }
  if ( final_cp1 && final_np1 < hole_k ) {
    memcpy ( &final_cp1[final_np1*size], cp, size*sizeof(point3d) );
    final_np1 ++;
  }
} /*OutputFinalBezPatch1*/

void OutputFinalBezPatch2 ( int n, int m, const double *cp, void *usrptr )
{
  int size;

  size = (n+1)*(m+1);
  if ( !final_cp2 ) {
    final_cp2 = malloc ( hole_k*size*sizeof(point3d) );
    final_np2 = 0;
    final_deg2 = n;
  }
  if ( final_cp2 && final_np2 < hole_k ) {
    memcpy ( &final_cp2[final_np2*size], cp, size*sizeof(point3d) );
    final_np2 ++;
  }
} /*OutputFinalBezPatch2*/

void OutputFinalBSPatch1 ( int n, int lknu, const double *knu,
                           int m, int lknv, const double *knv,
                           const double *cp, void *usrptr )
{
  int size;

  size = (lknu-n)*(lknv-m);
  if ( !final_cp1 ) {
    final_cp1 = malloc ( hole_k*size*sizeof(point3d) );
    if ( (final_knots1 = malloc ( (lknu+1)*sizeof(double) )) )
      memcpy ( final_knots1, knu, (lknu+1)*sizeof(double) );
    final_np1 = 0;
    final_deg1 = n;
    final_lkn1 = lknu;
  }
  if ( final_cp1 && final_np1 < hole_k ) {
    memcpy ( &final_cp1[final_np1*size], cp, size*sizeof(point3d) );
    final_np1 ++;
  }
} /*OutputFinalBSPatch1*/

void OutputFinalBSPatch2 ( int n, int lknu, const double *knu,
                           int m, int lknv, const double *knv,
                           const double *cp, void *usrptr )
{
  int size;

  size = (lknu-n)*(lknv-m);
  if ( !final_cp2 ) {
    final_cp2 = malloc ( hole_k*size*sizeof(point3d) );
    if ( (final_knots2 = malloc ( (lknu+1)*sizeof(double) )) )
      memcpy ( final_knots2, knu, (lknu+1)*sizeof(double) );
    final_np2 = 0;
    final_deg2 = n;
    final_lkn2 = lknu;
  }
  if ( final_cp2 && final_np2 < hole_k ) {
    memcpy ( &final_cp2[final_np2*size], cp, size*sizeof(point3d) );
    final_np2 ++;
  }
} /*OutputFinalBSPatch2*/

boolean UpdateFinalSurfaceNC ( int surfno, GHoleDomaind *domain, GHoptions *options,
                  double **acoeff,
                  void (*OutputBezPatch) ( int n, int m, const double *cp,
                                           void *usrptr ),
                  void (*OutputBSPatch) ( int n, int lknu, const double *knu,
                                          int m, int lknv, const double *knv,
                                          const double *cp, void *usrptr ),
                  int *final_np )
{
  char       *errstr;
  int        dimV0;
#ifdef SHOW_TIME
  struct tms t0, t1;
  clock_t    clk0, clk1;
  long       cps;
#endif

#ifdef SHOW_TIME
  clk0 = times ( &t0 );
#endif
  switch ( options->order ) {
case 1:
    if ( options->bezier ) {
      dimV0 = g1h_ExtV0SpaceDimd ( domain );
      *acoeff = malloc ( dimV0*3*sizeof(double) );
      if ( options->quasiG2 ) {
        if ( options->lin )
          g1h_Q2ExtFillHoled ( domain, 3, (double*)hole_cp, *acoeff,
                               NULL, OutputBezPatch );
        else
          g1h_Q2NLExtFillHoled ( domain, hole_cp, *acoeff,
                                 NULL, (OutputBezPatch3)OutputBezPatch );
      }
      else {
        if ( options->lin )
          g1h_ExtFillHoled ( domain, 3, (double*)hole_cp, *acoeff,
                             NULL, OutputBezPatch );
        else 
          g1h_NLExtFillHoled ( domain, hole_cp, *acoeff,
                               NULL, (OutputBezPatch3)OutputBezPatch );
      }
    }
    else if ( options->spline ) {
      if ( !options->spl_basis_valid )
        options->spl_basis_valid = g1h_ComputeSplBasisd ( domain,
                                  options->nk, options->m1, options->m2 );
      if ( options->spl_basis_valid ) {
        dimV0 = g1h_SplV0SpaceDimd ( domain );
        *acoeff = malloc ( dimV0*3*sizeof(double) );
        if ( options->quasiG2 ) {
          if ( options->lin )
            g1h_Q2SplFillHoled ( domain, 3, (double*)hole_cp, *acoeff,
                                 NULL, OutputBSPatch );
          else
            g1h_Q2NLSplFillHoled ( domain, hole_cp, *acoeff,
                                   NULL, (OutputBSPatch3)OutputBSPatch );
        }
        else {
          if ( options->lin )
            g1h_SplFillHoled ( domain, 3, (double*)hole_cp, *acoeff,
                               NULL, OutputBSPatch );
          else
            g1h_NLSplFillHoled ( domain, hole_cp, *acoeff,
                                 NULL, (OutputBSPatch3)OutputBSPatch );
        }
      }
    }
    else {
      if ( options->quasiG2 ) {
        dimV0 = g1h_V0SpaceDimd ( domain );
        *acoeff = malloc ( dimV0*3*sizeof(double) );
        if ( options->lin )
          g1h_Q2FillHoled ( domain, 3, (double*)hole_cp, *acoeff,
                            NULL, OutputBezPatch );
        else
          g1h_Q2NLFillHoled ( domain, hole_cp, *acoeff,
                              NULL, (OutputBezPatch3)OutputBezPatch );
      }
      else {
        if ( options->lin )
          g1h_FillHoled ( domain, 3, (double*)hole_cp, *acoeff,
                          NULL, OutputBezPatch );
        else
          g1h_NLFillHoled ( domain, hole_cp, *acoeff,
                            NULL, (OutputBezPatch3)OutputBezPatch );
      }
    }
    break;
case 2:
    if ( options->bezier ) {
      dimV0 = g2h_ExtV0SpaceDimd ( domain );
      *acoeff = malloc ( dimV0*3*sizeof(double) );
      if ( options->lin )
        g2h_ExtFillHoled ( domain, 3, (double*)hole_cp, *acoeff,
                           NULL, OutputBezPatch );
      else
        g2h_NLExtFillHoled ( domain, hole_cp, *acoeff,
                             NULL, (OutputBezPatch3)OutputBezPatch );
    }
    else if ( options->spline ) {
      if ( !options->spl_basis_valid )
        options->spl_basis_valid = g2h_ComputeSplBasisd ( domain,
                                  options->nk, options->m1, options->m2 );
      if ( options->spl_basis_valid ) {
        dimV0 = g2h_SplV0SpaceDimd ( domain );
        *acoeff = malloc ( dimV0*3*sizeof(double) );
        if ( options->lin )
          g2h_SplFillHoled ( domain, 3, (double*)hole_cp, *acoeff,
                             NULL, OutputBSPatch );
        else
          g2h_NLSplFillHoled ( domain, hole_cp, *acoeff,
                               NULL, (OutputBSPatch3)OutputBSPatch );
      }
    }
    else {
      dimV0 = g2h_V0SpaceDimd ( domain );
      *acoeff = malloc ( dimV0*3*sizeof(double) );
      if ( options->lin )
        g2h_FillHoled ( domain, 3, (double*)hole_cp, *acoeff,
                        NULL, OutputBezPatch );
      else
        g2h_NLFillHoled ( domain, hole_cp, *acoeff,
                          NULL, (OutputBezPatch3)OutputBezPatch );
    }
    break;
default:
    break;
  }
#ifdef SHOW_TIME
  clk1 = times ( &t1 );
  cps = sysconf ( _SC_CLK_TCK );
  printf ( "time = %f\n", (double)(clk1-clk0)/(double)cps );
#endif
  if ( *final_np == hole_k ) {
    options->surf_valid = true;
    return true;
  }
  else {
    switch ( options->order ) {
  case 1:
      g1h_GetErrorCoded ( domain, &errstr );
      break;
  case 2:
      g2h_GetErrorCoded ( domain, &errstr );
      break;
    }
    xge_DisplayErrorMessage ( errstr, -1 );
    options->surf_valid = false;
    return false;
  }
} /*UpdateFinalSurfaceNC*/

boolean UpdateFinalSurfaceC1 ( int surfno, GHoleDomaind *domain, GHoptions *options,
                  double **acoeff,
                  void (*OutputBezPatch) ( int n, int m, const double *cp,
                                           void *usrptr ),
                  void (*OutputBSPatch) ( int n, int lknu, const double *knu,
                                          int m, int lknv, const double *knv,
                                          const double *cp, void *usrptr ),
                  int *final_np )
{
  void       *sp;
  char       *errstr;
  int        dimV0, constr_no;
  vector3d   *rhs;
#ifdef SHOW_TIME
  struct tms t0, t1;
  clock_t    clk0, clk1;
  long       cps;
#endif

#ifdef SHOW_TIME
  clk0 = times ( &t0 );
#endif

  sp = pkv_GetScratchMemTop ();
  if ( options->spline && !options->spl_basis_valid ) {
    options->constr_matrix_valid = false;
    if ( options->order == 1 )
      options->spl_basis_valid = g1h_ComputeSplBasisd ( domain,
                                 options->nk, options->m1, options->m2 );
    else
      options->spl_basis_valid = g2h_ComputeSplBasisd ( domain,
                                 options->nk, options->m1, options->m2 );
    if ( !options->spl_basis_valid )
      goto failure;
  }
  if ( !options->constr_matrix_valid ) {
    if ( !ComputeConstraintMatrices ( surfno ) )
      goto failure;
  }
  constr_no = options->constr_no;
  rhs = pkv_GetScratchMem ( constr_no*sizeof(vector3d) );
  if ( !rhs )
    goto failure;
  if ( !ComputeConstraintRightSide ( surfno, rhs ) )
    goto failure;

  switch ( options->order ) {
case 1:
    if ( options->bezier ) {
      dimV0 = g1h_ExtV0SpaceDimd ( domain );
      *acoeff = malloc ( dimV0*3*sizeof(double) );
      if ( options->quasiG2 ) {
        if ( options->lin )
          g1h_Q2ExtFillHoleConstrd ( domain, 3, (double*)hole_cp,
                                     constr_no, &rhs[0].x,
                                     *acoeff, NULL, OutputBezPatch );
        else
          g1h_Q2NLExtFillHoleConstrd ( domain, hole_cp, constr_no, rhs, *acoeff,
                                       NULL, (OutputBezPatch3)OutputBezPatch );
      }
      else {
        if ( options->lin )
          g1h_ExtFillHoleConstrd ( domain, 3, (double*)hole_cp,
                         constr_no, &rhs[0].x, *acoeff, NULL, OutputBezPatch );
        else
          g1h_NLExtFillHoleConstrd ( domain, hole_cp, constr_no, rhs, *acoeff,
                                     NULL, (OutputBezPatch3)OutputBezPatch );
      }
    }
    else if ( options->spline ) {
      if ( options->spl_basis_valid ) {
        dimV0 = g1h_SplV0SpaceDimd ( domain );
        *acoeff = malloc ( dimV0*3*sizeof(double) );
        if ( options->quasiG2 ) {
          if ( options->lin )
            g1h_Q2SplFillHoleConstrd ( domain, 3, (double*)hole_cp,
                                       constr_no, &rhs[0].x,
                                       *acoeff, NULL, OutputBSPatch );
          else
            g1h_Q2NLSplFillHoleConstrd ( domain, hole_cp, constr_no, rhs,
                                         *acoeff, NULL,
                                         (OutputBSPatch3)OutputBSPatch );
        }
        else {
          if ( options->lin )
            g1h_SplFillHoleConstrd ( domain, 3, (double*)hole_cp,
                                     constr_no, &rhs[0].x, *acoeff,
                                     NULL, OutputBSPatch );
          else
            g1h_NLSplFillHoleConstrd ( domain, hole_cp, constr_no, rhs,
                                       *acoeff, NULL,
                                       (OutputBSPatch3)OutputBSPatch );
        }
      }
    }
    else {
      if ( options->quasiG2 ) {
        dimV0 = g1h_V0SpaceDimd ( domain );
        *acoeff = malloc ( dimV0*3*sizeof(double) );
        if ( options->lin )
          g1h_Q2FillHoleConstrd ( domain, 3, (double*)hole_cp,
                                  constr_no, &rhs[0].x, *acoeff,
                                  NULL, OutputBezPatch );
        else
          g1h_Q2NLFillHoleConstrd ( domain, hole_cp, constr_no, rhs, *acoeff,
                                    NULL, (OutputBezPatch3)OutputBezPatch );
      }
      else {
        if ( options->lin )
          g1h_FillHoleConstrd ( domain, 3, (double*)hole_cp,
                                constr_no, &rhs[0].x, *acoeff, NULL,
                                OutputBezPatch );
        else
          g1h_NLFillHoleConstrd ( domain, hole_cp, constr_no, rhs, *acoeff,
                                  NULL, (OutputBezPatch3)OutputBezPatch );
      }
    }
    break;
case 2:
    if ( options->bezier ) {
      dimV0 = g2h_ExtV0SpaceDimd ( domain );
      *acoeff = malloc ( dimV0*3*sizeof(double) );
      if ( options->lin )
        g2h_ExtFillHoleConstrd ( domain, 3, (double*)hole_cp,
                       constr_no, &rhs[0].x, *acoeff, NULL,
                       OutputBezPatch );
      else
        g2h_NLExtFillHoleConstrd ( domain, hole_cp, constr_no, rhs, *acoeff,
                                   NULL, (OutputBezPatch3)OutputBezPatch );
    }
    else if ( options->spline ) {
      if ( options->spl_basis_valid ) {
        dimV0 = g2h_SplV0SpaceDimd ( domain );
        *acoeff = malloc ( dimV0*3*sizeof(double) );
        if ( options->lin )
          g2h_SplFillHoleConstrd ( domain, 3, (double*)hole_cp,
                                   constr_no, &rhs[0].x, *acoeff,
                                   NULL, OutputBSPatch );
        else
          g2h_NLSplFillHoleConstrd ( domain, hole_cp, constr_no, rhs, *acoeff,
                                     NULL, (OutputBSPatch3)OutputBSPatch );
      }
    }
    else {
      dimV0 = g2h_V0SpaceDimd ( domain );
      *acoeff = malloc ( dimV0*3*sizeof(double) );
      if ( options->lin )
        g2h_FillHoleConstrd ( domain, 3, (double*)hole_cp,
                              constr_no, &rhs[0].x, *acoeff, NULL,
                              OutputBezPatch );
      else
        g2h_NLFillHoleConstrd ( domain, hole_cp, constr_no, rhs, *acoeff,
                                NULL, (OutputBezPatch3)OutputBezPatch );
    }
    break;
default:
    break;
  }
#ifdef SHOW_TIME
  clk1 = times ( &t1 );
  cps = sysconf ( _SC_CLK_TCK );
  printf ( "time = %f\n", (double)(clk1-clk0)/(double)cps );
#endif
  if ( *final_np == hole_k )
    goto success;
  else {
    switch ( options->order ) {
  case 1:
      g1h_GetErrorCoded ( domain, &errstr );
      break;
  case 2:
      g2h_GetErrorCoded ( domain, &errstr );
      break;
    }
    xge_DisplayErrorMessage ( errstr, -1 );
    goto failure;
  }

success:
  pkv_SetScratchMemTop ( sp );
  options->surf_valid = true;
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  options->surf_valid = false;
  return false;
} /*UpdateFinalSurfaceC1*/

boolean UpdateFinalSurfaceC2 ( int surfno, GHoleDomaind *domain, GHoptions *options,
                  double **acoeff,
                  void (*OutputBezPatch) ( int n, int m, const double *cp,
                                          void *usrptr ),
                  void (*OutputBSPatch) ( int n, int lknu, const double *knu,
                                          int m, int lknv, const double *knv,
                                          const double *cp, void *usrptr ),
                  int *final_np )
{
  void       *sp;
  char       *errstr;
  int        dimV0, constr_no;
  double     *rhs;
#ifdef SHOW_TIME
  struct tms t0, t1;
  clock_t    clk0, clk1;
  long       cps;
#endif

#ifdef SHOW_TIME
  clk0 = times ( &t0 );
#endif

  sp = pkv_GetScratchMemTop ();
  if ( options->spline && !options->spl_basis_valid ) {
    options->constr_matrix_valid = false;
    if ( options->order == 1 )
      options->spl_basis_valid = g1h_ComputeSplBasisd ( domain,
                                 options->nk, options->m1, options->m2 );
    else
      options->spl_basis_valid = g2h_ComputeSplBasisd ( domain,
                                 options->nk, options->m1, options->m2 );
    if ( !options->spl_basis_valid )
      goto failure;
  }
  if ( !options->constr_matrix_valid ) {
    if ( !ComputeConstraintMatricesAlt ( surfno ) )
      goto failure;
  }
  constr_no = options->constr_no;
  rhs = pkv_GetScratchMemd ( constr_no );
  if ( !rhs )
    goto failure;
  if ( !ComputeConstraintRightSideAlt ( surfno, rhs ) )
    goto failure;

  switch ( options->order ) {
case 1:
    if ( options->bezier ) {
      dimV0 = g1h_ExtV0SpaceDimd ( domain );
      *acoeff = malloc ( dimV0*3*sizeof(double) );
      if ( options->quasiG2 ) {
        if ( options->lin )
          g1h_Q2ExtFillHoleAltConstrd ( domain, 3, (double*)hole_cp,
                                  constr_no, rhs, *acoeff, NULL,
                                  OutputBezPatch );
        else
          g1h_Q2NLExtFillHoleAltConstrd ( domain, hole_cp, constr_no, rhs,
                                  *acoeff, NULL,
                                  (OutputBezPatch3)OutputBezPatch );
      }
      else {
        if ( options->lin )
          g1h_ExtFillHoleAltConstrd ( domain, 3, (double*)hole_cp,
                                  constr_no, rhs, *acoeff, NULL,
                                  OutputBezPatch );
        else
          g1h_NLExtFillHoleAltConstrd ( domain, hole_cp, constr_no, rhs,
                                  *acoeff, NULL,
                                  (OutputBezPatch3)OutputBezPatch );
      }
    }
    else if ( options->spline ) {
      if ( options->spl_basis_valid ) {
        dimV0 = g1h_SplV0SpaceDimd ( domain );
        *acoeff = malloc ( dimV0*3*sizeof(double) );
        if ( options->quasiG2 ) {
          if ( options->lin )
            g1h_Q2SplFillHoleAltConstrd ( domain, 3, (double*)hole_cp,
                              constr_no, rhs, *acoeff, NULL,
                              OutputBSPatch );
          else
            g1h_Q2NLSplFillHoleAltConstrd ( domain, hole_cp, constr_no, rhs,
                              *acoeff, NULL,
                              (OutputBSPatch3)OutputBSPatch );
        }
        else {
          if ( options->lin )
            g1h_SplFillHoleAltConstrd ( domain, 3, (double*)hole_cp,
                              constr_no, rhs, *acoeff, NULL,
                              OutputBSPatch );
          else
            g1h_NLSplFillHoleAltConstrd ( domain, hole_cp, constr_no, rhs,
                              *acoeff, NULL,
                              (OutputBSPatch3)OutputBSPatch );
        }
      }
    }
    else {
      dimV0 = g1h_V0SpaceDimd ( domain );
      *acoeff = malloc ( dimV0*3*sizeof(double) );
      if ( options->quasiG2 ) {
        if ( options->lin )
          g1h_Q2FillHoleAltConstrd ( domain, 3, (double*)hole_cp,
                                  constr_no, rhs, *acoeff, NULL,
                                  OutputBezPatch );
        else
          g1h_Q2NLFillHoleAltConstrd ( domain, hole_cp, constr_no, rhs,
                                  *acoeff, NULL,
                                  (OutputBezPatch3)OutputBezPatch );
      }
      else {
        if ( options->lin )
          g1h_FillHoleAltConstrd ( domain, 3, (double*)hole_cp,
                                   constr_no, rhs, *acoeff, NULL,
                                   OutputBezPatch );
        else
          g1h_NLFillHoleAltConstrd ( domain, hole_cp, constr_no, rhs,
                                     *acoeff, NULL,
                                     (OutputBezPatch3)OutputBezPatch );
      }
    }
    break;

case 2:
    if ( options->bezier ) {
      dimV0 = g2h_ExtV0SpaceDimd ( domain );
      *acoeff = malloc ( dimV0*3*sizeof(double) );
      if ( options->lin )
        g2h_ExtFillHoleAltConstrd ( domain, 3, (double*)hole_cp,
                                constr_no, rhs, *acoeff, NULL,
                                OutputBezPatch );
      else
        g2h_NLExtFillHoleAltConstrd ( domain, hole_cp, constr_no, rhs,
                                *acoeff, NULL,
                                (OutputBezPatch3)OutputBezPatch );
    }
    else if ( options->spline ) {
      if ( options->spl_basis_valid ) {
        dimV0 = g2h_SplV0SpaceDimd ( domain );
        *acoeff = malloc ( dimV0*3*sizeof(double) );
        if ( options->lin )
          g2h_SplFillHoleAltConstrd ( domain, 3, (double*)hole_cp,
                            constr_no, rhs, *acoeff, NULL,
                            OutputBSPatch );
        else
          g2h_NLSplFillHoleAltConstrd ( domain, hole_cp, constr_no, rhs,
                            *acoeff, NULL,
                            (OutputBSPatch3)OutputBSPatch );
      }
    }
    else {
      dimV0 = g2h_V0SpaceDimd ( domain );
      *acoeff = malloc ( dimV0*3*sizeof(double) );
      if ( options->lin )
        g2h_FillHoleAltConstrd ( domain, 3, (double*)hole_cp,
                                 constr_no, rhs, *acoeff, NULL,
                                 OutputBezPatch );
      else
        g2h_NLFillHoleAltConstrd ( domain, hole_cp, constr_no, rhs,
                                   *acoeff, NULL,
                                   (OutputBezPatch3)OutputBezPatch );
    }
    break;

default:
    break;
  }
#ifdef SHOW_TIME
  clk1 = times ( &t1 );
  cps = sysconf ( _SC_CLK_TCK );
  printf ( "time = %f\n", (double)(clk1-clk0)/(double)cps );
#endif
  if ( *final_np == hole_k )
    goto success;
  else {
    switch ( options->order ) {
  case 1:
      g1h_GetErrorCoded ( domain, &errstr );
      break;
  case 2:
      g2h_GetErrorCoded ( domain, &errstr );
      break;
    }
    xge_DisplayErrorMessage ( errstr, -1 );
    goto failure;
  }

success:
  pkv_SetScratchMemTop ( sp );
  options->surf_valid = true;
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  options->surf_valid = false;
  return false;
} /*UpdateFinalSurfaceC2*/

boolean UpdateFinalSurface1 ( void )
{
  if ( options1.constr_type == 0 || options1.constr_no == 0 ) {
    return UpdateFinalSurfaceNC ( 1, domain1, &options1, &acoeff1,
                                  OutputFinalBezPatch1, OutputFinalBSPatch1,
                                  &final_np1 );
  }
  else if ( options1.constr_type == 1 ) {
    return UpdateFinalSurfaceC1 ( 1, domain1, &options1, &acoeff1,
                                  OutputFinalBezPatch1, OutputFinalBSPatch1,
                                  &final_np1 );
  }
  else {
    return UpdateFinalSurfaceC2 ( 1, domain1, &options1, &acoeff1,
                                  OutputFinalBezPatch1, OutputFinalBSPatch1,
                                  &final_np1 );
  }
} /*UpdateFinalSurface1*/

boolean UpdateFinalSurface2 ( void )
{
  if ( options2.constr_type == 0 || options2.constr_no == 0 ) {
    return UpdateFinalSurfaceNC ( 2, domain2, &options2, &acoeff2,
                                  OutputFinalBezPatch2, OutputFinalBSPatch2,
                                  &final_np2 );
  }
  else if ( options2.constr_type == 1 ) {
    return UpdateFinalSurfaceC1 ( 2, domain2, &options2, &acoeff2,
                                  OutputFinalBezPatch2, OutputFinalBSPatch2,
                                  &final_np2 );
  }
  else {
    return UpdateFinalSurfaceC2 ( 2, domain2, &options2, &acoeff2,
                                  OutputFinalBezPatch2, OutputFinalBSPatch2,
                                  &final_np2 );
  }
} /*UpdateFinalSurface2*/

boolean SurfFast ( GHoptions *opt )
{
  return ( (opt->fast = (!opt->spline && opt->lin)) );
} /*SurfFast*/

