
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>

#include <GL/gl.h>
#include <GL/glx.h>

#include "pkgeom.h"
#include "camera.h"

#include "xgledit.h"
#include "xgeprivate.h"

/* ///////////////////////////////////////////////////////////////////////// */
/* This procedure sets up the GL projection/modelview matrices as described  */
/* by the CameraRec structure, processed by the libcamera procedures; this   */
/* version may be useful for antialiasing and depth of focus effect, using   */
/* the accumulation buffer; depth of focus is relevant only for the          */
/* perspective projection                                                    */

boolean xgle_SetGLaccCamerad ( CameraRecd *CPos,
                               double pixdx, double pixdy,
                               double eyedx, double eyedy, double focus )
{
  double w, h, ac, shx, shy, s, dx, dy;

  if ( CPos->ncplanes < 6 || CPos->zmax <= CPos->zmin )
    return false;
  if ( !CPos->parallel && CPos->zmin <= 0.0 )
    return false;

  glViewport ( 0, 0, CPos->width/CPos->magnification,
               CPos->height/CPos->magnification );
  glMatrixMode ( GL_PROJECTION );
  glLoadIdentity ();
  if ( CPos->parallel ) {
    w = CPos->vd.para.wdt;
    h = CPos->vd.para.hgh;
    if ( CPos->upside )
      glOrtho ( -0.5*w, 0.5*w, -0.5*h, 0.5*h, CPos->zmin, CPos->zmax );
    else
      glOrtho ( -0.5*w, 0.5*w, 0.5*h, -0.5*h, CPos->zmin, CPos->zmax );
    eyedx = pixdx*CPos->vd.para.wdt/CPos->width;
    eyedy = pixdy*CPos->vd.para.hgh/CPos->height;
  }
  else {
    ac = (double)CPos->height/(CPos->aspect*(double)CPos->width);
    w = sqrt ( 1.0/(1.0+ac*ac) );
    h = ac*w;
    shx = CPos->vd.persp.dxi0/CPos->xscale;
    shy = CPos->vd.persp.deta0/CPos->yscale;
    s = CPos->zmin/CPos->vd.persp.f;
    dx = -pixdx*s*w/(double)CPos->width + eyedx*CPos->zmin/focus;
    dy = -pixdy*s*h/(double)CPos->height + eyedy*CPos->zmin/focus;
    if ( CPos->upside )
      glFrustum ( s*(shx-0.5*w)+dx, s*(shx+0.5*w)+dx,
                  s*(shy-0.5*h)+dy, s*(shy+0.5*h)+dy,
                  CPos->zmin, CPos->zmax );
    else
      glFrustum ( s*(shx-0.5*w)+dx, s*(shx+0.5*w)+dx,
                  s*(shy+0.5*h)-dy, s*(shy-0.5*h)-dy,
                  CPos->zmin, CPos->zmax );
  }
  glMatrixMode ( GL_MODELVIEW );
  glLoadIdentity ();
  glScaled ( 1.0, 1.0, -1.0 );
  glRotated ( -(180.0/PI)*CPos->psi,   0.0, 0.0, 1.0 );
  glRotated ( -(180.0/PI)*CPos->theta, 1.0, 0.0, 0.0 );
  glRotated ( -(180.0/PI)*CPos->phi,   0.0, 0.0, 1.0 );
  glTranslated ( -CPos->position.x-eyedx, -CPos->position.y-eyedy,
                 -CPos->position.z );

  return true;
} /*xgle_SetGLCamerad*/

