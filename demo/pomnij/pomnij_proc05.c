
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <sys/types.h>
#include <sys/times.h>
#include <signal.h>
#include <unistd.h>
#include <setjmp.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <fpu_control.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camerad.h"
#include "multibs.h"
#include "g1blendingd.h"
#include "g2blendingd.h"
#include "xgedit.h"
#include "xgeipc.h"

#include "bsfile.h"

#include "pomnijipc.h"
#include "pomnij_proc.h"
#include "proc_regmem.h"

/* ////////////////////////////////////////////////////////////////////////// */
int     ncp, optlknu, optlknv, fcp;
point3d *ccp = NULL, *occp;
void    *optdata = NULL;
int     itn;
boolean finished = true;
struct tms start;

/* ////////////////////////////////////////////////////////////////////////// */
boolean G1OptimizeLMTInit ( void )
{
  void   *sp;
  int    i;
  double dO, dM;

  sp = pkv_GetScratchMemTop ();
  if ( !knotsu || !knotsv || !cpoints || !(dim == 3 || dim == 4) ||
       degu != 2 || degv != 2 || lknu <= 5 || lknv <= 5 ||
       options.NLConst <= 0.0 || options.maxiter < 1 || options.gcont != 1 )
    goto failure;

  if ( !G1AdjustOptRange () )
    goto failure;
  ncp = (lknu-degu)*(lknv-degv);
  PKV_MALLOC ( ccp, ncp*sizeof(point3d) );
  if ( !ccp )
    goto failure;
  if ( dim == 3 )
    memcpy ( ccp, cpoints, ncp*sizeof(point3d) );
  else {
    for ( i = 0; i < ncp; i++ )
      Point4to3d ( (point4d*)&cpoints[4*i], &ccp[i] );
  }
  SetupG1BLKnots ();
  G1ClampedBoundaryToFree ( ccp );
  PreTransformation ( ncp, ccp );

  if ( options.closed ) {
        /* ************* */
    goto failure;
  }
  else {
    optlknu = options.opt_range[1]-options.opt_range[0] + 7;
    optlknv = options.opt_range[3]-options.opt_range[2] + 7;
    fcp = (options.opt_range[0]-2)*(lknv-2) + (options.opt_range[2]-2);
  }
  occp = &ccp[fcp];
  if ( options.dumpdata )
    WriteThePatch ( 3, 2, optlknu, knotsu, 2, optlknv, knotsv, 3*(lknv-2),
                    (double*)occp );
  dO = (lknu-4)*(lknu-4) + (lknv-4)*(lknv-4);
  dM = g1bl_SurfNetDiameterSqd ( lknu, lknv, 3*(lknv-2), ccp );
  if ( options.nconstr ) {
        /* ************* */
    goto failure;
  }
  else {
/*no_constr:*/
    times ( &start );
    if ( options.closed ) {
        /* ************* */
      goto failure;
    }
    else {
      if ( !g1bl_InitBlSurfaceOptLMTd ( optlknu, optlknv, 3*(lknv-2), occp,
                     options.NLConst, dO, dM,
                     options.nkn1, options.nkn2, &optdata ) )
        goto failure;
    }
  }
  finished = false;
  itn = 0;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*G1OptimizeLMTInit*/

boolean G1OptimizeLMTIter ( void )
{
  void       *sp;
  int        i;
  struct tms stop;
  double     time;
  point3d    *acp;

  sp = pkv_GetScratchMemTop ();
  if ( optdata ) {
    if ( options.nconstr ) {
        /* ************* */
      goto failure;
    }
    else {
      if ( options.closed ) {
        /* ************* */
        goto failure;
      }
      else {
        if ( !g1bl_IterBlSurfaceOptLMTd ( optdata, &finished ) )
          goto failure;
      }
    }
    itn ++;
    if ( itn >= options.maxiter || finished || options.send_partial ) {
        /* convert the partial or final result */
      acp = (point3d*)pkv_GetScratchMem ( ncp*sizeof(point3d) );
      if ( !acp )
        goto failure;

      memcpy ( acp, ccp, ncp*sizeof(point3d) );
      if ( !PostTransformation ( ncp, acp ) )
        goto failure;
      if ( clu || clv ) {
        for ( i = 0; i <= lknu; i++ )
          knotsu[i] = (double)i;
        for ( i = 0; i <= lknv; i++ )
          knotsv[i] = (double)i;
        G1FreeBoundaryToClamped ( acp );
      }
      if ( dim == 3 )
        memcpy ( cpoints, acp, ncp*sizeof(point3d) );
      else {
        pkv_Selectd ( ncp, 3, 3, 4, acp, cpoints );
        for ( i = 0; i < ncp; i++ )
          cpoints[4*i+3] = 1.0;
      }
    }
    if ( itn >= options.maxiter || finished ) {
      finished = true;
      if ( options.nconstr ) {
        /* ************* */
        goto failure;
      }
      else {
        if ( options.closed ) {
        /* ************* */
          goto failure;
        }
        else
          g1bl_OptLMTDeallocated ( &optdata );
      }
      PKV_FREE ( ccp );
      times ( &stop );
      time = (double)(stop.tms_utime-start.tms_utime)/(double)(sysconf(_SC_CLK_TCK));
      printf ( "time = %7.2f\n", time );
    }
  }
  else
    goto failure;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( optdata ) {
    finished = true;
    if ( options.nconstr ) {
        /* ************* */
    }
    else {
      if ( options.closed )
        ;
      else
        g1bl_OptLMTDeallocated ( &optdata );
    }
  }
  if ( ccp ) PKV_FREE ( ccp );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*G1OptimizeLMTIter*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean G2OptimizeLMTInit ( void )
{
  void    *sp, *sp1;
  int     i, bls, nvars, nuccp;
  double  dO, dM;
  point3d *uccp, *ouccp;

  sp = pkv_GetScratchMemTop ();
  if ( !knotsu || !knotsv || !cpoints || !(dim == 3 || dim == 4) ||
       degu != 3 || degv != 3 || lknu <= 9 || lknv <= 9 ||
       options.NLConst <= 0.0 || options.maxiter < 1 || options.gcont != 2 )
    goto failure;

  if ( !G2AdjustOptRange () )
    goto failure;
  ncp = (lknu-degu)*(lknv-degv);
  PKV_MALLOC ( ccp, ncp*sizeof(point3d) );
  if ( !ccp )
    goto failure;
  if ( dim == 3 )
    memcpy ( ccp, cpoints, ncp*sizeof(point3d) );
  else {
    for ( i = 0; i < ncp; i++ )
      Point4to3d ( (point4d*)&cpoints[4*i], &ccp[i] );
  }
  SetupG2BLKnots ();
  G2ClampedBoundaryToFree ( ccp );
  PreTransformation ( ncp, ccp );

  if ( options.closed ) {
    options.opt_range[0] = 0;
    options.opt_range[1] = lknu-10;
    optlknu = lknu;
    optlknv = options.opt_range[3]-options.opt_range[2] + 10;
    fcp = options.opt_range[0]*(lknv-3) + (options.opt_range[2]-3);
  }
  else {
    optlknu = options.opt_range[1]-options.opt_range[0] + 10;
    optlknv = options.opt_range[3]-options.opt_range[2] + 10;
    fcp = (options.opt_range[0]-3)*(lknv-3) + (options.opt_range[2]-3);
  }
  occp = &ccp[fcp];
  if ( options.dumpdata )
    WriteThePatch ( 3, 3, optlknu, knotsu, 3, optlknv, knotsv, 3*(lknv-3),
                    (double*)occp );
  dO = (lknu-6)*(lknu-6) + (lknv-6)*(lknv-6);
  dM = g2bl_SurfNetDiameterSqd ( lknu, lknv, 3*(lknv-3), ccp );
  if ( options.nconstr ) {
    sp1 = pkv_GetScratchMemTop ();
    bls = lknv-9;
    if ( options.closed )
      nvars = bls*(lknu-6);
    else
      nvars = bls*(lknu-9);
    PKV_MALLOC ( cmat, 9*nvars*bls*nucknots*sizeof(double) );
    PKV_MALLOC ( crhs, 3*bls*nucknots*sizeof(double) );
    if ( !cmat || !crhs )
      goto failure;
    nuccp = nucknots*(lknv-3);
    uccp = pkv_GetScratchMem ( nuccp*sizeof(point3d) );
    for ( i = 0; i < nuccp; i++ )
      Point4to3d ( (point4d*)&ucurvcp[4*i], &uccp[i] );
    for ( i = 0; i < nucknots; ) {
      if ( !options.closed )
        ucknots[i] -= (double)(options.opt_range[0]-3);
      if ( ucknots[i] > 3.0 && ucknots[i] < (double)(optlknu-3) )
        i ++;
      else {
        nucknots --;
        ucknots[i] = ucknots[nucknots];
        memcpy ( &uccp[i*(lknv-3)], &uccp[nucknots*(lknv-3)],
                 (lknv-3)*sizeof(point3d) );
      }
    }
printf ( "nucknots = %d\n", nucknots );
    if ( !nucknots ) {
      options.nconstr = 0;
      goto no_constr;
    }
    if ( options.dumpdata && options.nconstr )
      WriteTheConstraints ( 3, nucknots, ucknots, 3, lknv, knotsv,
                            3*(lknv-3), (double*)uccp );
    PreTransformation ( nuccp, uccp );
    if ( options.closed ) {
      ouccp = &uccp[options.opt_range[2]-3];
      if ( !g2bl_SetupClosedUNLConstraintsd ( optlknu, optlknv, 3*(lknv-3), occp,
                                              nucknots, ucknots, 3*(lknv-3), ouccp,
                                              &options.nconstr, cmat, crhs ) )
        goto failure;
        pkv_SetScratchMemTop ( sp1 );
        times ( &start );
        if ( !g2bl_ClosedInitBlSurfaceConstrOptLMTd ( optlknu, optlknv, 3*(lknv-3),
                                  occp, options.nconstr, cmat, crhs,
                                  options.NLConst, dO, dM,
                                  options.nkn1, options.nkn2, &optdata ) )
          goto failure;
    }
    else {
      ouccp = &uccp[options.opt_range[2]-3];
      if ( !g2bl_SetupUNLConstraintsd ( optlknu, optlknv, 3*(lknv-3), occp,
                                        nucknots, ucknots, 3*(lknv-3), ouccp,
                                        &options.nconstr, cmat, crhs ) )
        goto failure;
      pkv_SetScratchMemTop ( sp1 );
      times ( &start );
      if ( !g2bl_InitBlSurfaceConstrOptLMTd ( optlknu, optlknv, 3*(lknv-3), occp,
                                     options.nconstr, cmat, crhs,
                                     options.NLConst, dO, dM,
                                     options.nkn1, options.nkn2, &optdata ) )
        goto failure;
    }
  }
  else {
no_constr:
    times ( &start );
    if ( options.closed ) {
      if ( !g2bl_ClosedInitBlSurfaceOptLMTd ( optlknu, optlknv, 3*(lknv-3), occp,
                     options.NLConst, dO, dM,
                     options.nkn1, options.nkn2, &optdata ) )
        goto failure;
    }
    else {
      if ( !g2bl_InitBlSurfaceOptLMTd ( optlknu, optlknv, 3*(lknv-3), occp,
                     options.NLConst, dO, dM,
                     options.nkn1, options.nkn2, &optdata ) )
        goto failure;
    }
  }
  finished = false;
  itn = 0;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*G2OptimizeLMTInit*/

boolean G2OptimizeLMTIter ( void )
{
  void       *sp;
  int        i;
  struct tms stop;
  double     time;
  point3d    *acp;

  sp = pkv_GetScratchMemTop ();
  if ( optdata ) {
    if ( options.nconstr ) {
      if ( options.closed ) {
        if ( !g2bl_ClosedIterBlSurfaceConstrOptLMTd ( optdata, &finished ) )
          goto failure;
      }
      else {
        if ( !g2bl_IterBlSurfaceConstrOptLMTd ( optdata, &finished ) )
          goto failure;
      }
    }
    else {
      if ( options.closed ) {
        if ( !g2bl_ClosedIterBlSurfaceOptLMTd ( optdata, &finished ) )
          goto failure;
      }
      else {
        if ( !g2bl_IterBlSurfaceOptLMTd ( optdata, &finished ) )
          goto failure;
      }
    }
    itn ++;
    if ( itn >= options.maxiter || finished || options.send_partial ) {
        /* convert the partial or final result */
      acp = (point3d*)pkv_GetScratchMem ( ncp*sizeof(point3d) );
      if ( !acp )
        goto failure;

      memcpy ( acp, ccp, ncp*sizeof(point3d) );
      if ( !PostTransformation ( ncp, acp ) )
        goto failure;
      if ( clu || clv ) {
        for ( i = 0; i <= lknu; i++ )
          knotsu[i] = (double)i;
        for ( i = 0; i <= lknv; i++ )
          knotsv[i] = (double)i;
        G2FreeBoundaryToClamped ( acp );
      }
      if ( dim == 3 )
        memcpy ( cpoints, acp, ncp*sizeof(point3d) );
      else {
        pkv_Selectd ( ncp, 3, 3, 4, acp, cpoints );
        for ( i = 0; i < ncp; i++ )
          cpoints[4*i+3] = 1.0;
      }
    }
    if ( itn >= options.maxiter || finished ) {
      finished = true;
      if ( options.nconstr ) {
        if ( options.closed )
          g2bl_ClosedConstrOptLMTDeallocated ( &optdata );
        else
          g2bl_ConstrOptLMTDeallocated ( &optdata );
      }
      else {
        if ( options.closed )
          g2bl_ClosedOptLMTDeallocated ( &optdata );
        else
          g2bl_OptLMTDeallocated ( &optdata );
      }
      PKV_FREE ( ccp );
      times ( &stop );
      time = (double)(stop.tms_utime-start.tms_utime)/(double)(sysconf(_SC_CLK_TCK));
      printf ( "time = %7.2f\n", time );
    }
  }
  else
    goto failure;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( optdata ) {
    finished = true;
    if ( options.nconstr ) {
      if ( options.closed )
        g2bl_ClosedConstrOptLMTDeallocated ( &optdata );
      else
        g2bl_ConstrOptLMTDeallocated ( &optdata );
    }
    else {
      if ( options.closed )
        g2bl_ClosedOptLMTDeallocated ( &optdata );
      else
        g2bl_OptLMTDeallocated ( &optdata );
    }
  }
  if ( ccp ) PKV_FREE ( ccp );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*G2OptimizeLMTIter*/

