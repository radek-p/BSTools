
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkgeom.h"
#include "multibs.h"
#include "convh.h"
#include "camerad.h"

#include "xgedit.h"
#include "spl3d.h"

/* ///////////////////////////////////////////////////////////////////////// */
void DisplayKnots ( boolean high )
{
  xge_widget *er;
  int        i;
  int        x, y, uxa, uxb, vya, vyb;
  double     uu;

  er = kwind.er;
  xgeSetForeground ( xgec_Green4 );
  x = kwind.CPos.xmin;
  y = kwind.CPos.ymin + kwind.CPos.height;
  xgeDrawLine ( kwind.CPos.xmin-KNOT_MARGIN+5, y,
                kwind.CPos.xmin+kwind.CPos.width-5, y );
  xgeDrawLine ( x, kwind.CPos.ymin+5,
                x, kwind.CPos.ymin+kwind.CPos.height+KNOT_MARGIN-5 );
  xgeSetForeground ( xgec_White );

  uxa = xge_T2KnotWindMapKnotU ( &kwind, knots_u[degree_u] );
  uxa = max ( er->x+1, uxa );
  uxb = xge_T2KnotWindMapKnotU ( &kwind, knots_u[lastknot_u-degree_u] );
  uxb = min ( er->x+er->w-2, uxb );
  xgeDrawLine ( uxa, y, uxb, y );
  if ( !high )
    xgeSetForeground ( xgec_Grey4 );
  for ( uu = (double)(knots_u[1]-1.0), uxb = 0, i = 1;
        i < lastknot_u;
        i++ ) {
    if ( knots_u[i] > uu ) {
      uxb = y + 5;
      uxa = xge_T2KnotWindMapKnotU ( &kwind, (uu = knots_u[i]) );
    }
    else
      uxb += 5;
    xgeFillRectangle ( 3, 3, uxa-1, uxb-1 );
  }
  
  xgeSetForeground ( xgec_White );
  vya = xge_T2KnotWindMapKnotV ( &kwind, knots_v[degree_v] );
  vya = max ( er->y+1, vya );
  vyb = xge_T2KnotWindMapKnotV ( &kwind, knots_v[lastknot_v-degree_v] );
  vyb = min ( er->y+er->h-2, vyb );
  xgeDrawLine ( x, vya, x, vyb );
  for ( uu = (double)(knots_v[1]-1.0), vyb = 0, i = 1;
        i < lastknot_v;
        i++ ) {
    if ( knots_v[i] > uu ) {
      vyb = x - 5;
      vya = xge_T2KnotWindMapKnotV ( &kwind, (uu = knots_v[i]) );
    }
    else
      vyb -= 5;
    xgeFillRectangle ( 3, 3, vyb-1, vya-1 );
  }
} /*DisplayKnots*/

void DisplayKnotLines ( void )
{
  int   i;
  int   uxa, uxb, vya, vyb;
  double uu;

  xgeSetForeground ( xgec_Grey2 );

  vya = xge_T2KnotWindMapKnotV ( &kwind, knots_v[degree_v] );
  vyb = xge_T2KnotWindMapKnotV ( &kwind, knots_v[lastknot_v-degree_v] );
  for ( uu = (double)(knots_u[1]-1.0), i = degree_u;
        i <= lastknot_u-degree_u;
        i++ ) {
    if ( knots_u[i] > uu ) {
      uxa = xge_T2KnotWindMapKnotU ( &kwind, (uu = knots_u[i]) );
      xgeDrawLine ( uxa, vya, uxa, vyb );
    }
  }
  
  uxa = xge_T2KnotWindMapKnotU ( &kwind, knots_u[degree_u] );
  uxb = xge_T2KnotWindMapKnotU ( &kwind, knots_u[lastknot_u-degree_u] );
  for ( uu = (double)(knots_v[1]-1.0), i = degree_v;
        i <= lastknot_v-degree_v;
        i++ ) {
    if ( knots_v[i] > uu ) {
      vya = xge_T2KnotWindMapKnotV ( &kwind, (uu = knots_v[i]) );
      xgeDrawLine ( uxa, vya, uxb, vya );
    }
  }
} /*DisplayKnotLines*/

void TriangleMark ( int x, int y )
{
  xgeDrawPoint ( x, y-1 );
  xgeDrawLine ( x-1, y, x+1, y );
  xgeDrawLine ( x-1, y+1, x+1, y+1 );
  xgeDrawLine ( x-2, y+2, x+2, y+2 );
} /*TriangleMark*/

void DisplayBlendingKnots ( void )
{
  int    i, uxa, uxb, y, vya, vyb;
  double uu;

  if ( n_blending_constraints <= 0 )
    return;
  y = kwind.CPos.ymin + kwind.CPos.height;
  vya = xge_T2KnotWindMapKnotV ( &kwind, knots_v[degree_v] );
  vyb = xge_T2KnotWindMapKnotV ( &kwind, knots_v[lastknot_v-degree_v] );
  xgeSetForeground ( xgec_Orchid1 );
  for ( uu = blending_constr_knots[1]-1.0, uxa = uxb = 0, i = 1;
        i <= n_blending_constraints;
        i++ ) {
    if ( blending_constr_knots[i] > uu ) {
      uxb = y + 5;
      uxa = xge_T2KnotWindMapKnotU ( &kwind, (uu = blending_constr_knots[i]) );
    }
    else
      uxb += 5;
    xgeDrawLine ( uxa, vya, uxa, vyb );
    TriangleMark ( uxa, uxb );
  }
} /*DisplayBlendingKnots*/

void DisplayDomain ( void )
{
  DisplayKnots ( !sw_blending_constraints );
  DisplayKnotLines ();
  if ( sw_blending_constraints )
    DisplayBlendingKnots ();
} /*DisplayDomain*/

double GrevilleAbscissa ( int degree, double *knots, int i )
{
  int j;
  double xi;

  xi = knots[i+1];
  for ( j = 2; j <= degree; j++ )
    xi += knots[i+j];
  return xi/(double)degree;
} /*GrevilleAbscissa*/

void DisplayDomainNet ( void )
{
  void *sp;
  int  i, j, k, pitch;
  int  xy, x0, x1, y0, y1;
  int  *xg, *yg;

  x0 = kwind.CPos.xmin-KNOT_MARGIN;
  xy = xge_T2KnotWindMapKnotU ( &kwind, GrevilleAbscissa ( degree_u, knots_u, 0 ) );
  x0 = max ( x0, xy );
  x1 = kwind.CPos.xmin+kwind.CPos.width;
  xy = xge_T2KnotWindMapKnotU ( &kwind, GrevilleAbscissa ( degree_u, knots_u,
                                lastknot_u-degree_u-1 ) );
  x1 = min ( x1, xy );
  if ( x1 < x0 )
    return;
  y0 = kwind.CPos.ymin;
  xy = xge_T2KnotWindMapKnotV ( &kwind, GrevilleAbscissa ( degree_v, knots_v,
                                lastknot_v-degree_v-1 ) );
  y0 = max ( y0, xy );
  y1 = kwind.CPos.ymin+kwind.CPos.height+KNOT_MARGIN;
  xy = xge_T2KnotWindMapKnotV ( &kwind, GrevilleAbscissa ( degree_v, knots_v, 0 ) );
  y1 = min ( y1, xy );
  if ( y1 < y0 )
    return;

  sp = pkv_GetScratchMemTop ();
  xg = pkv_GetScratchMem ( (lastknot_u-degree_u)*sizeof(int) );
  yg = pkv_GetScratchMem ( (lastknot_v-degree_v)*sizeof(int) );
  if ( !xg || !yg )
    goto leave_it;

  xgeSetForeground ( xgec_Green1 );
  for ( i = 0; i < lastknot_u-degree_u; i++ ) {
    xg[i] = xy = xge_T2KnotWindMapKnotU ( &kwind, GrevilleAbscissa ( degree_u, knots_u, i ) );
    if ( xy >= x0 && xy <= x1 )
      xgeDrawLine ( xy, y0, xy, y1 );
  }
  for ( i = 0; i < lastknot_v-degree_v; i++ ) {
    yg[i] = xy = xge_T2KnotWindMapKnotV ( &kwind, GrevilleAbscissa ( degree_v, knots_v, i ) );
    if ( xy >= y0 && xy <= y1 )
      xgeDrawLine ( x0, xy, x1, xy );
  }

  pitch = lastknot_v-degree_v;
  xgeSetForeground ( xgec_Yellow );
  if ( win1_contents == WIN1_BLENDING ) {
    for ( i = k = 0; i < lastknot_u-degree_u;  i++ )
      for ( j = 0;  j < lastknot_v-degree_v;  j++, k++ ) {
        if ( !(mkpoints[k] & 0x01) )
          xgeFillRectangle ( 3, 3, xg[i]-1, yg[j]-1 );
      }
    xgeSetForeground ( xgec_OrangeRed );
    for ( i = k = 0; i < lastknot_u-degree_u;  i++ )
      for ( j = 0;  j < lastknot_v-degree_v;  j++, k++ ) {
        if ( (mkpoints[k] & 0x01) )
          xgeFillRectangle ( 3, 3, xg[i]-1, yg[j]-1 );
      }
    xgeSetForeground ( xgec_Green );
    if ( kwind.closed_u ) {
      if ( blending_opt_part[0] <= blending_opt_part[1] ) {
        for ( i = k = 0; i < lastknot_u-2*degree_u;  i++ )
          for ( j = 0;  j < lastknot_v-degree_v;  j++, k++ ) {
            if ( i >= blending_opt_part[0] && i <= blending_opt_part[1] &&
                 j >= blending_opt_part[2] && j <= blending_opt_part[3] ) 
              xgeFillRectangle ( 3, 3, xg[i]-1, yg[j]-1 );
          }
        for ( ; i < lastknot_u-degree_u; i++ )
          for ( j = 0;  j < lastknot_v-degree_v;  j++, k++ ) {
            if ( i-kwind.clcKu >= blending_opt_part[0] &&
                 i-kwind.clcKu <= blending_opt_part[1] &&
                 j >= blending_opt_part[2] && j <= blending_opt_part[3] ) 
              xgeFillRectangle ( 3, 3, xg[i]-1, yg[j]-1 );
          }
      }
      else {
        for ( i = k = 0; i < lastknot_u-2*degree_u;  i++ )
          for ( j = 0;  j < lastknot_v-degree_v;  j++, k++ ) {
            if ( (i <= blending_opt_part[1] ||
                  i >= blending_opt_part[0]) &&
                 j >= blending_opt_part[2] && j <= blending_opt_part[3] ) 
              xgeFillRectangle ( 3, 3, xg[i]-1, yg[j]-1 );
          }
        for ( ; i < lastknot_u-degree_u; i++ )
          for ( j = 0;  j < lastknot_v-degree_v;  j++, k++ ) {
            if ( (i-kwind.clcKu <= blending_opt_part[1] ||
                  i-kwind.clcKu >= blending_opt_part[0]) &&
                 j >= blending_opt_part[2] && j <= blending_opt_part[3] ) 
              xgeFillRectangle ( 3, 3, xg[i]-1, yg[j]-1 );
          }
      }
    }
    else {
      for ( i = k = 0; i < lastknot_u-degree_u;  i++ )
        for ( j = 0;  j < lastknot_v-degree_v;  j++, k++ ) {
          if ( i >= blending_opt_part[0] && i <= blending_opt_part[1] &&
               j >= blending_opt_part[2] && j <= blending_opt_part[3] ) 
            xgeFillRectangle ( 3, 3, xg[i]-1, yg[j]-1 );
        }
    }
  }
  else {
    for ( i = 0; i < lastknot_u-degree_u; i++ )
      if ( xg[i] >= x0 && xg[i] <= x1 )
        for ( j = 0;  j < lastknot_v-degree_v;  j++ ) {
          k = i*pitch+j;
          if ( yg[j] >= y0 && yg[j] <= y1 && !(mkpoints[k] & 0x01) )
            xgeFillRectangle ( 3, 3, xg[i]-1, yg[j]-1 );
        }

    xgeSetForeground ( xgec_OrangeRed );
    for ( i = 0; i < lastknot_u-degree_u; i++ )
      if ( xg[i] >= x0 && xg[i] <= x1 )
        for ( j = 0; j < lastknot_v-degree_v; j++ ) {
          k = i*pitch+j;
          if ( yg[j] >= y0 && yg[j] <= y1 && (mkpoints[k] & 0x01) )
            xgeFillRectangle ( 3, 3, xg[i]-1, yg[j]-1 );
        }
  }

leave_it:
  pkv_SetScratchMemTop ( sp );
} /*DisplayDomainNet*/

