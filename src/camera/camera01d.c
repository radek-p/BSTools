
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2008                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <math.h>

#include "camera.h"

#include "msgpool.h"


void CameraInitFramed ( CameraRecd *CPos,
                        boolean parallel, boolean upside,
                        short width, short height, short xmin, short ymin,
                        double aspect, int ncplanes )
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
} /*CameraInitFramed*/

void CameraSetMagd ( CameraRecd *CPos, byte mag )
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
} /*CameraSetMagd*/

static void SetupClipPlaned ( double nvx, double nvy, double nvz,
                              const trans3d *CTr, const point3d *pp,
                              vector4d *plane )
{
  vector3d nvect, *nv;

  SetVector3d ( &nvect, nvx, nvy, nvz );
  nv = (vector3d*)plane;
  TransContra3d ( CTr, &nvect, nv );
  NormalizeVector3d ( nv );
  plane->w = -DotProduct3d ( pp, nv );
} /*SetupClipPlaned*/

void CameraSetMappingd ( CameraRecd *CPos )
{
  /* this procedure calculates the transformations and other parameters   */
  /* in a position specified in the CameraRec record fields. It is called */
  /* after the position initialisation and after every change.            */
  double  xi0, eta0, a, w, h, hs, wdt, hgh, diag;
  point3d clp;

  w = (double)CPos->width;
  h = (double)CPos->height;
  a = CPos->aspect;
  hs = h/a;
  if ( CPos->parallel ) {
    switch ( CPos->vd.para.dim_case ) {
  case 0:           /* diag is given */
      diag = CPos->vd.para.diag;
      CPos->xscale = sqrt(w*w+hs*hs)/diag;
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
      CPos->vd.para.diag = diag = sqrt ( wdt*wdt + hgh*hgh );
      CPos->xscale = w/wdt;
      CPos->yscale = h/hgh;
      break;

  default:
      pkv_SignalError ( LIB_CAMERA, 1, ERRMSG_2 );
    }
    xi0  = CPos->xmin + 0.5*w;
    eta0 = CPos->ymin + 0.5*h;

    IdentTrans3d ( &CPos->CTr );
    ShiftTrans3d ( &CPos->CTr, -CPos->position.x, -CPos->position.y,
                   -CPos->position.z );
    EulerRotTrans3d ( &CPos->CTr, -CPos->psi, -CPos->theta, -CPos->phi );
    ScaleTrans3d  ( &CPos->CTr, CPos->xscale, CPos->yscale, 1.0 );
    ShiftTrans3d  ( &CPos->CTr, xi0, eta0, 0.0 );
  
    IdentTrans3d ( &CPos->CTrInv );
    ShiftTrans3d  ( &CPos->CTrInv, -xi0, -eta0, 0.0 );
    ScaleTrans3d  ( &CPos->CTrInv, 1.0/CPos->xscale, 1.0/CPos->yscale, 1.0 );
    EulerRotTrans3d ( &CPos->CTrInv, CPos->phi, CPos->theta, CPos->psi );
    ShiftTrans3d ( &CPos->CTrInv, CPos->position.x, CPos->position.y,
                   CPos->position.z );

  /* setup the clipping planes */
    SetPoint3d ( &clp, (double)CPos->xmin, (double)CPos->ymin, 0.0 );
    TransPoint3d ( &CPos->CTrInv, &clp, &clp );
    SetupClipPlaned ( 0.0, 1.0, 0.0, &CPos->CTr, &clp,
                      &CPos->cplane[CPLANE_TOP] );
    SetupClipPlaned ( 1.0, 0.0, 0.0, &CPos->CTr, &clp,
                      &CPos->cplane[CPLANE_LEFT] );
    SetPoint3d ( &clp, (double)(CPos->xmin+CPos->width),
                 (double)(CPos->ymin+CPos->height), 0.0 );
    TransPoint3d ( &CPos->CTrInv, &clp, &clp );
    SetupClipPlaned ( 0.0, -1.0, 0.0, &CPos->CTr, &clp,
                      &CPos->cplane[CPLANE_BOTTOM] );
    SetupClipPlaned ( -1.0, 0.0, 0.0, &CPos->CTr, &clp,
                      &CPos->cplane[CPLANE_RIGHT] );
  }
  else {
  /* find the global to camera coordinates transformation */
    CPos->xscale = sqrt ( w*w + hs*hs );
    CPos->yscale = CPos->xscale*a;

    IdentTrans3d ( &CPos->CTr );
    ShiftTrans3d( &CPos->CTr, -CPos->position.x, -CPos->position.y,
                  -CPos->position.z );
    EulerRotTrans3d ( &CPos->CTr, -CPos->psi, -CPos->theta, -CPos->phi );
    ScaleTrans3d ( &CPos->CTr, CPos->vd.persp.f * CPos->xscale,
                   CPos->vd.persp.f * CPos->yscale, 1.0 );

    IdentTrans3d ( &CPos->CTrInv );
    ScaleTrans3d ( &CPos->CTrInv, 1.0 / (CPos->vd.persp.f * CPos->xscale),
                   1.0 / (CPos->vd.persp.f * CPos->yscale), 1.0 );
    EulerRotTrans3d ( &CPos->CTrInv, CPos->phi, CPos->theta, CPos->psi );
    ShiftTrans3d ( &CPos->CTrInv, CPos->position.x, CPos->position.y,
                   CPos->position.z );

  /* find the shift on the screen */
    CPos->vd.persp.xi0 = CPos->xmin + 0.5*w + CPos->vd.persp.dxi0;
    CPos->vd.persp.eta0 = CPos->ymin + 0.5*h + CPos->vd.persp.deta0;

  /* setup the clipping planes */
    SetupClipPlaned ( 0.0, 1.0, CPos->vd.persp.eta0-(double)CPos->ymin,
         &CPos->CTr, &CPos->position, &CPos->cplane[CPLANE_TOP] );
    SetupClipPlaned ( 0.0, -1.0,
         (double)(CPos->ymin+CPos->height)-CPos->vd.persp.eta0,
         &CPos->CTr, &CPos->position, &CPos->cplane[CPLANE_BOTTOM] );
    SetupClipPlaned ( 1.0, 0.0, CPos->vd.persp.xi0-(double)CPos->xmin,
         &CPos->CTr, &CPos->position, &CPos->cplane[CPLANE_LEFT] );
    SetupClipPlaned ( -1.0, 0.0,
         (double)(CPos->xmin+CPos->width)-CPos->vd.persp.xi0,
         &CPos->CTr, &CPos->position, &CPos->cplane[CPLANE_RIGHT] );
  }
  if ( CPos->ncplanes >= 6 ) {
        /* there are two additional clipping planes - near and far */
    SetPoint3d ( &clp, 0.0, 0.0, CPos->zmin );
    TransPoint3d ( &CPos->CTrInv, &clp, &clp );
    SetupClipPlaned ( 0.0, 0.0, 1.0, &CPos->CTr, &clp,   
                      &CPos->cplane[CPLANE_NEAR] ); 
    SetPoint3d ( &clp, 0.0, 0.0, CPos->zmax );
    TransPoint3d ( &CPos->CTrInv, &clp, &clp );
    SetupClipPlaned ( 0.0, 0.0, -1.0, &CPos->CTr, &clp,  
                      &CPos->cplane[CPLANE_FAR] );
  }
} /*CameraSetMappingd*/

