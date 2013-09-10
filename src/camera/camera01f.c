
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <math.h>

#include "camera.h"


void CameraInitFramef ( CameraRecf *CPos,
                        boolean parallel, boolean upside,
                        short width, short height, short xmin, short ymin,
                        float aspect, int ncplanes )
{
  /* this procedure initializes the variables in the CameraRec description, */
  /* which are related with the size and aspect of the raster image.        */
  /* it should be called before any other Camera procedure.                 */
  CPos->parallel = parallel;
  CPos->upside = upside;
  CPos->width = width;
  CPos->height = height;
  CPos->xmin = xmin;
  CPos->ymin = ymin;
  CPos->aspect = aspect;
  CPos->magnification = 1;
  if ( CPos->parallel ) {
    CPos->vd.para.dim_case = 0;
    CPos->vd.para.diag = 1.0;
  }
  else {
    CPos->vd.persp.dxi0 = 0.0;
    CPos->vd.persp.deta0 = 0.0;
  }
  switch ( ncplanes ) {
case 0:      /* leave the previous value */
    break;
case 6:
    CPos->ncplanes = 6;
    if ( parallel ) { CPos->zmin = -1.0;  CPos->zmax =   1.0; }
               else { CPos->zmin =  1.0;  CPos->zmax = 100.0; }
    break;
default:
    CPos->ncplanes = 4;
    break;
  }
} /*CameraInitFramef*/

void CameraSetMagf ( CameraRecf *CPos, byte mag )
{
  /* by default, the magnification is 1; higher values may be used for */
  /* antialiasing by supersampling on a denser subpixel grid.          */
  CPos->width = (short)(CPos->width / CPos->magnification * mag);
  CPos->height = (short)(CPos->height / CPos->magnification * mag);
  CPos->xmin = (short)(CPos->xmin / CPos->magnification * mag);
  CPos->ymin = (short)(CPos->ymin / CPos->magnification * mag);
  if ( CPos->parallel ) {
  }
  else {
    CPos->vd.persp.dxi0 = CPos->vd.persp.dxi0 / CPos->magnification * mag;
    CPos->vd.persp.deta0 = CPos->vd.persp.deta0 / CPos->magnification * mag;
  }
  CPos->magnification = mag;
} /*CameraSetMagf*/

static void SetupClipPlanef ( float nvx, float nvy, float nvz,
                              const trans3f *CTr, const point3f *pp,
                              vector4f *plane )
{
  vector3f nvect, *nv;

  SetVector3f ( &nvect, nvx, nvy, nvz );
  nv = (vector3f*)plane;
  TransContra3f ( CTr, &nvect, nv );
  NormalizeVector3f ( nv );
  plane->w = (float)(-DotProduct3f ( pp, nv ));
} /*SetupClipPlanef*/

boolean CameraSetMappingf ( CameraRecf *CPos )
{
  /* this procedure calculates the transformations and other parameters   */
  /* in a position specified in the CameraRec record fields. It is called */
  /* after the position initialisation and after every change.            */
  float   xi0, eta0, a, w, h, hs, wdt, hgh, diag;
  point3f clp;

  w = (float)CPos->width;
  h = (float)CPos->height;
  a = CPos->aspect;
  hs = h/a;
  if ( CPos->parallel ) {
    switch ( CPos->vd.para.dim_case ) {
  case 0:           /* diag is given */
      diag = CPos->vd.para.diag;
      CPos->xscale = (float)(sqrt(w*w+hs*hs)/diag);
      CPos->yscale = CPos->xscale*a;
      CPos->vd.para.wdt = wdt = w/CPos->xscale;
      CPos->vd.para.hgh = hgh = h/CPos->yscale;
      break;

  case 1:           /* wdt is given  */
      wdt = CPos->vd.para.wdt;
      CPos->vd.para.hgh = hgh = wdt*(hs/w);
      goto contin;

  case 2:           /* hgh is given  */
      hgh = CPos->vd.para.hgh;
      CPos->vd.para.wdt = wdt = hgh*(a*w/h);   

  contin:
      CPos->vd.para.diag = diag = (float)sqrt ( wdt*wdt + hgh*hgh );
      CPos->xscale = w/wdt;
      CPos->yscale = h/hgh;
      break;

  default:
      PKV_SIGNALERROR ( LIB_CAMERA, 4, ERRMSG_4 );
      return false;
    }
    xi0  = (float)(CPos->xmin + 0.5*w);
    eta0 = (float)(CPos->ymin + 0.5*h);

    IdentTrans3f ( &CPos->CTr );
    ShiftTrans3f ( &CPos->CTr, -CPos->position.x, -CPos->position.y,
                   -CPos->position.z );
    EulerRotTrans3f ( &CPos->CTr, -CPos->psi, -CPos->theta, -CPos->phi );
    ScaleTrans3f  ( &CPos->CTr, CPos->xscale, CPos->yscale, 1.0 );
    ShiftTrans3f  ( &CPos->CTr, xi0, eta0, 0.0 );
  
    IdentTrans3f ( &CPos->CTrInv );
    ShiftTrans3f  ( &CPos->CTrInv, -xi0, -eta0, 0.0 );
    ScaleTrans3f  ( &CPos->CTrInv, (float)(1.0/CPos->xscale),
                    (float)(1.0/CPos->yscale), 1.0 );
    EulerRotTrans3f ( &CPos->CTrInv, CPos->phi, CPos->theta, CPos->psi );
    ShiftTrans3f ( &CPos->CTrInv, CPos->position.x, CPos->position.y,
                   CPos->position.z );

  /* setup the clipping planes */
    SetPoint3f ( &clp, (float)CPos->xmin, (float)CPos->ymin, 0.0 );
    TransPoint3f ( &CPos->CTrInv, &clp, &clp );
    SetupClipPlanef ( 0.0, 1.0, 0.0, &CPos->CTr, &clp,
                      &CPos->cplane[CPLANE_TOP] );
    SetupClipPlanef ( 1.0, 0.0, 0.0, &CPos->CTr, &clp,
                      &CPos->cplane[CPLANE_LEFT] );
    SetPoint3f ( &clp, (float)(CPos->xmin+CPos->width),
                 (float)(CPos->ymin+CPos->height), 0.0 );
    TransPoint3f ( &CPos->CTrInv, &clp, &clp );
    SetupClipPlanef ( 0.0, -1.0, 0.0, &CPos->CTr, &clp,
                      &CPos->cplane[CPLANE_BOTTOM] );
    SetupClipPlanef ( -1.0, 0.0, 0.0, &CPos->CTr, &clp,
                      &CPos->cplane[CPLANE_RIGHT] );
  }
  else {
  /* find the global to camera coordinates transformation */
    CPos->xscale = (float)sqrt ( w*w + hs*hs );
    CPos->yscale = CPos->xscale*a;

    IdentTrans3f ( &CPos->CTr );
    ShiftTrans3f( &CPos->CTr, -CPos->position.x, -CPos->position.y,
                  -CPos->position.z );
    EulerRotTrans3f ( &CPos->CTr, -CPos->psi, -CPos->theta, -CPos->phi );
    ScaleTrans3f ( &CPos->CTr, CPos->vd.persp.f * CPos->xscale,
                   CPos->vd.persp.f * CPos->yscale, 1.0 );

    IdentTrans3f ( &CPos->CTrInv );
    ScaleTrans3f ( &CPos->CTrInv, (float)(1.0 / (CPos->vd.persp.f * CPos->xscale)),
                   (float)(1.0 / (CPos->vd.persp.f * CPos->yscale)), 1.0 );
    EulerRotTrans3f ( &CPos->CTrInv, CPos->phi, CPos->theta, CPos->psi );
    ShiftTrans3f ( &CPos->CTrInv, CPos->position.x, CPos->position.y,
                   CPos->position.z );

  /* find the shift on the screen */
    CPos->vd.persp.xi0 = (float)(CPos->xmin + 0.5*w + CPos->vd.persp.dxi0);
    CPos->vd.persp.eta0 = (float)(CPos->ymin + 0.5*h + CPos->vd.persp.deta0);

  /* setup the clipping planes */
    SetupClipPlanef ( 0.0, 1.0, CPos->vd.persp.eta0-(float)CPos->ymin,
         &CPos->CTr, &CPos->position, &CPos->cplane[CPLANE_TOP] );
    SetupClipPlanef ( 0.0, -1.0,
         (float)(CPos->ymin+CPos->height)-CPos->vd.persp.eta0,
         &CPos->CTr, &CPos->position, &CPos->cplane[CPLANE_BOTTOM] );
    SetupClipPlanef ( 1.0, 0.0, CPos->vd.persp.xi0-(float)CPos->xmin,
         &CPos->CTr, &CPos->position, &CPos->cplane[CPLANE_LEFT] );
    SetupClipPlanef ( -1.0, 0.0,
         (float)(CPos->xmin+CPos->width)-CPos->vd.persp.xi0,
         &CPos->CTr, &CPos->position, &CPos->cplane[CPLANE_RIGHT] );
  }
  if ( CPos->ncplanes >= 6 ) {
        /* there are two additional clipping planes - near and far */
    SetPoint3f ( &clp, 0.0, 0.0, CPos->zmin );
    TransPoint3f ( &CPos->CTrInv, &clp, &clp );
    SetupClipPlanef ( 0.0, 0.0, 1.0, &CPos->CTr, &clp,
                      &CPos->cplane[CPLANE_NEAR] );
    SetPoint3f ( &clp, 0.0, 0.0, CPos->zmax );
    TransPoint3f ( &CPos->CTrInv, &clp, &clp );
    SetupClipPlanef ( 0.0, 0.0, -1.0, &CPos->CTr, &clp,
                      &CPos->cplane[CPLANE_FAR] );
  }
  return true;
} /*CameraSetMappingf*/

