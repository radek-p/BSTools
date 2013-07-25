
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "anima.h"


double ic_knots[MAX_KEY_POSES+3];
boolean sw_periodic = false;
boolean spline_ok = false;

static int naparams, nkposes, lastbsknot;
static double *keypose = NULL;
static double *acpoints = NULL;
static double *pos_knots;
static double bs_knots[MAX_KEY_POSES+2*IC_DEGREE+1];


int GetNPoses ( void )
{
  return nkposes;
} /*GetNPoses*/

boolean SetNArtParams ( int nparams )
{
  if ( nparams < 1 || nparams > MAX_APARAMS )
    return false;
  pos_knots = &ic_knots[1];
  naparams = nparams;
  if ( keypose )  free ( keypose );
  if ( acpoints ) free ( acpoints );
  keypose  = (double*)malloc ( naparams*MAX_KEY_POSES*sizeof(double) );
  acpoints = (double*)malloc ( naparams*(MAX_KEY_POSES+2*IC_DEGREE)*sizeof(double) );
  nkposes = 0;
  spline_ok = false;
  return keypose && acpoints;
} /*SetNArtParams*/

boolean EnterKeyPose ( double time, double *artp )
{
  int     i;
  boolean found;

  if ( nkposes < MAX_KEY_POSES ) {
    for ( found = false, i = 0;  i < nkposes;  i++ )
      if ( fabs(time-pos_knots[i]) < TIME_MIN_STEP ) {
        found = true;
        break;
      }
    if ( found )
      memcpy ( &keypose[i*naparams], artp, naparams*sizeof(double) );
    else {
      for ( i = nkposes; i > 0 && time < pos_knots[i-1]; i-- ) {
        pos_knots[i] = pos_knots[i-1];
        memcpy ( &keypose[i*naparams], &keypose[(i-1)*naparams],
                 naparams*sizeof(double) );
      }
      pos_knots[i] = time;
      memcpy ( &keypose[i*naparams], artp, naparams*sizeof(double) );
      nkposes ++;
      ic_knots[0] = ic_knots[1];
      ic_knots[nkposes+1] = ic_knots[nkposes];
    }
    spline_ok = false;
    return true;
  }
  else
    return false;
} /*EnterKeyPose*/

boolean DeleteKeyPose ( int k )
{
  if ( nkposes > 2 && k >= 0 && k < nkposes ) {
    if ( k < nkposes-1 ) {
      memmove ( &pos_knots[k-1], &pos_knots[k], (nkposes-k)*sizeof(double) );
      memmove ( &keypose[(k-1)*naparams], &keypose[k*naparams],
                (nkposes-k)*naparams*sizeof(double) );
    }
    nkposes --;
    spline_ok = false;
    return true;
  }
  else
    return false;
} /*DeleteKeyPose*/

boolean MoveKeyKnot ( int k )
{
  spline_ok = false;
  return true;
} /*MoveKeyKnot*/

boolean ConstructInterpCurve ( void )
{
  void *sp;

  sp = pkv_GetScratchMemTop ();
  if ( sw_periodic ) {
    mbs_multiBSCubicClosedInterpd ( nkposes-1, pos_knots, 1, naparams,
                              0, keypose, &lastbsknot, bs_knots, 0, acpoints );
  }
  else {
    mbs_multiBSCubicInterpd ( nkposes-1, pos_knots, 1, naparams,
                              0, keypose, 0,
                              BS3_BC_SECOND_DER0, NULL,
                              BS3_BC_SECOND_DER0, NULL,
                              &lastbsknot, bs_knots, 0, acpoints );
  }
  spline_ok = true;
  pkv_SetScratchMemTop ( sp );
  return true;
} /*ConstructInterpCurve*/

void ClampIt ( double *artp )
{
  int i;

  for ( i = 0; i < naparams; i++ ) {
    artp[i] = max ( artp[i], 0.0 );
    artp[i] = min ( artp[i], 1.0 );
  }
} /*ClampIt*/

boolean GetPose ( double time, double *artp )
{
  void   *sp;
  double t, *ap;
  int    i;

  sp = pkv_GetScratchMemTop ();
  if ( nkposes < 1 )
    goto failure;
  for ( i = 1; i < nkposes; i++ )
    if ( pos_knots[i-1] >= pos_knots[i] )
      goto failure;

        /* clamp the time */
  time = max ( time, pos_knots[0] );
  time = min ( time, pos_knots[nkposes-1] );

  switch ( nkposes ) {
case 1:
    memcpy ( artp, keypose, naparams*sizeof(double) );
    break;

case 2:   /* linear interpolation */
    t = (time-pos_knots[0])/(pos_knots[1]-pos_knots[0]);
    pkn_MatrixLinCombd ( 1, naparams, 0, &keypose[0], 1.0-t,
                         0, &keypose[naparams], t, 0, artp );
    break;

case 3:   /* quadratic interpolation - Aitken algorithm */
    ap = pkv_GetScratchMemd ( naparams );
    if ( !ap )
      goto failure;
    t = (time-pos_knots[0])/(pos_knots[1]-pos_knots[0]);
    pkn_MatrixLinCombd ( 1, naparams, 0, &keypose[0], 1.0-t,
                         0, &keypose[naparams], t, 0, artp );
    t = (time-pos_knots[1])/(pos_knots[2]-pos_knots[1]);
    pkn_MatrixLinCombd ( 1, naparams, 0, &keypose[naparams], 1.0-t,
                         0, &keypose[2*naparams], t, 0, ap );
    t = (time-pos_knots[0])/(pos_knots[2]-pos_knots[0]);
    pkn_MatrixLinCombd ( 1, naparams, 0, artp, 1.0-t, 0, ap, t, 0, artp );
    pkv_SetScratchMemTop ( sp );
    break;

case 4:  /* cubic interpolation - Aitken algorithm */
    ap = pkv_GetScratchMemd ( 2*naparams );
    if ( !ap )
      goto failure;
    t = (time-pos_knots[0])/(pos_knots[1]-pos_knots[0]);
    pkn_MatrixLinCombd ( 1, naparams, 0, &keypose[0], 1.0-t,
                         0, &keypose[naparams], t, 0, artp );
    t = (time-pos_knots[1])/(pos_knots[2]-pos_knots[1]);
    pkn_MatrixLinCombd ( 1, naparams, 0, &keypose[naparams], 1.0-t,
                         0, &keypose[2*naparams], t, 0, ap );
    t = (time-pos_knots[2])/(pos_knots[3]-pos_knots[2]);
    pkn_MatrixLinCombd ( 1, naparams, 0, &keypose[2*naparams], 1.0-t,
                         0, &keypose[3*naparams], t, 0, &ap[naparams] );
    t = (time-pos_knots[0])/(pos_knots[2]-pos_knots[0]);
    pkn_MatrixLinCombd ( 1, naparams, 0, artp, 1.0-t, 0, ap, t, 0, artp );
    t = (time-pos_knots[1])/(pos_knots[3]-pos_knots[1]);
    pkn_MatrixLinCombd ( 1, naparams, 0, ap, 1.0-t, 0, &ap[naparams], t, 0, ap );
    t = (time-pos_knots[0])/(pos_knots[3]-pos_knots[0]);
    pkn_MatrixLinCombd ( 1, naparams, 0, artp, 1.0-t, 0, ap, t, 0, artp );
    break;

default:  /* spline interpolation */
    if ( !spline_ok ) {
      if ( !ConstructInterpCurve () )
        goto failure;
    }
    mbs_multideBoord ( 3, lastbsknot, bs_knots, 1, naparams,
                       0, acpoints, time, artp );
    break;
  }
  ClampIt ( artp );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*GetPose*/

