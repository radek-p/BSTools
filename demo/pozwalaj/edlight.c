
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <GL/gl.h>
#include <GL/glu.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "egholed.h"
#include "bsmesh.h"
#include "camera.h"
#include "raybezf.h"
#include "mengerc.h"
#include "bsfile.h"
#include "pkrender.h"
#include "xgedit.h"
#include "xgledit.h"

#include "widgets.h"
#include "editor.h"
#include "editor_bezc.h"
#include "editor_bsc.h"
#include "editor_bezp.h"
#include "editor_bsp.h"
#include "editor_bsm.h"
#include "editor_bsh.h"
#include "pozwalaj.h"
#include "edlight.h"

#define POINT_TOL 10

    /* renderer lights */
boolean  sw_edit_light[RENDERING_NLIGHTS] =
  {false, false, false, false};
vector3d render_light_dir[RENDERING_NLIGHTS] =
  {{0.0,0.0,1.0},{0.0,1.0,0.0},{1.0,0.0,0.0},{0.0,0.0,-1.0}};
double   render_light_int[RENDERING_NLIGHTS+1] =
  {0.75, 0.0, 0.0, 0.0, 0.3};

static int     light_num;
static double  light_vect_factor;
static point3d light_centre, light_pt[RENDERING_NLIGHTS];

static vector3d saved_light_dir[RENDERING_NLIGHTS];

static void EdLightComputeVFactor ( xge_3Dwind *ww )
{
  Box3d  *box;
  double r;

  box = &ww->WinBBox;
  SetPoint3d ( &light_centre, 0.5*(box->x0+box->x1),
               0.5*(box->y0+box->y1), 0.5*(box->z0+box->z1) );
  r = min ( box->x1-box->x0, box->y1-box->y0 );
  r = min ( r, box->z1-box->z0 );
  light_vect_factor = 0.45*r;
} /*EdLightComputeFactor*/

boolean EdLightFindNearestPoint ( CameraRecd *CPos, short x, short y )
{
  int     i, d, dmin;
  point3d q;

  dmin = POINT_TOL+1;
  for ( i = 0; i < RENDERING_NLIGHTS; i++ )
    if ( sw_edit_light[i] ) {
      CameraProjectPoint3d ( CPos, &light_pt[i], &q );
      d = (int)(fabs((double)x-q.x) + fabs((double)y-q.y));
      if ( d < dmin ) {
        dmin = d;
        light_num = i;
      }
    }
  return dmin < POINT_TOL;
} /*EdLightFindNearestPoint*/

void EdLightSetPoint ( CameraRecd *CPos, short x, short y )
{
  point3d p, q;

  if ( light_num >= 0 && light_num < RENDERING_NLIGHTS ) {
    CameraProjectPoint3d ( CPos, &light_pt[light_num], &q );
    q.x = (double)x;  q.y = (double)y;
    CameraUnProjectPoint3d ( CPos, &q, &p );
    SubtractPoints3d ( &p, &light_centre, &q );
    NormalizeVector3d ( &q );
    render_light_dir[light_num] = q;
  }
} /*EdLightSetPoint*/

void EdLightSavePoints ( void )
{
  memcpy ( saved_light_dir, render_light_dir, RENDERING_NLIGHTS*sizeof(vector3d) );
} /*EdLightSavePoints*/

void EdLightTransformPoints ( trans3d *tr )
{
  int i;

  for ( i = 0; i < RENDERING_NLIGHTS; i++ )
    if ( sw_edit_light[i] ) {
      TransVector3d ( tr, &saved_light_dir[i], &render_light_dir[i] );
      NormalizeVector3d ( &render_light_dir[i] );
    }
} /*EdLightTransformPoints*/

void EdLightDisplay ( xge_3Dwind *ww )
{
  int     i;

  EdLightComputeVFactor ( ww );
  glPointSize ( 5.0 );
  glBegin ( GL_POINTS );
    glColor3fv ( xglec_LightSkyBlue );
    glVertex3dv ( &light_centre.x );
    for ( i = 0; i < RENDERING_NLIGHTS; i++ )
      if ( sw_edit_light[i] ) {
        AddVector3Md ( &light_centre, &render_light_dir[i], light_vect_factor,
                       &light_pt[i] );
        glVertex3dv ( &light_pt[i].x );
      }
  glEnd ();
  glBegin ( GL_LINES );
    glColor3fv ( xglec_Snow );
    for ( i = 0; i < RENDERING_NLIGHTS; i++ )
      if ( sw_edit_light[i] ) {
        glVertex3dv ( &light_centre.x );
        glVertex3dv ( &light_pt[i].x );
      }
  glEnd ();
} /*EdLightDisplay*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean sw_reflection_frame = false, sw_highlight_frame = false,
        sw_sections_frame = false;

vector3d edshapef_sectiondir = {0.0,0.0,1.0};
point3d  edshapef_reflection_frame[3] = {{0.0,0.0,5.0},{1.0,0.0,5.0},{0.0,1.0,5.0}};
point3d  edshapef_highlight_frame[3]  = {{0.0,0.0,5.0},{1.0,0.0,5.0},{0.0,1.0,5.0}};

boolean EdShapeFuncFindNearestPoint ( CameraRecd *CPos, short x, short y )
{
  int     i, d, dmin;
  point3d p, q;

  dmin = POINT_TOL+1;
  if ( sw_sections_frame ) {
    AddVector3Md ( &light_centre, &edshapef_sectiondir, light_vect_factor, &p );
    CameraProjectPoint3d ( CPos, &p, &q );
    d = (int)(fabs((double)x-q.x) + fabs((double)y-q.y));
    if ( d < dmin ) {
      dmin = d;
      light_num = RENDERING_NLIGHTS;
    }
  }
  if ( sw_reflection_frame ) {
    for ( i = 0; i < 3; i++ ) {
      CameraProjectPoint3d ( CPos, &edshapef_reflection_frame[i], &q );
      d = (int)(fabs((double)x-q.x) + fabs((double)y-q.y));
      if ( d < dmin ) {
        dmin = d;
        light_num = RENDERING_NLIGHTS+1+i;
      }
    }
  }
  if ( sw_highlight_frame ) {
    for ( i = 0; i < 3; i++ ) {
      CameraProjectPoint3d ( CPos, &edshapef_highlight_frame[i], &q );
      d = (int)(fabs((double)x-q.x) + fabs((double)y-q.y));
      if ( d < dmin ) {
        dmin = d;
        light_num = RENDERING_NLIGHTS+4+i;
      }
    }
  }
  return dmin < POINT_TOL;
} /*EdShapeFuncFindNearestPoint*/

void EdShapeFuncSetPoint ( CameraRecd *CPos, short x, short y )
{
  point3d p, q;
  int     i, j;

  switch ( light_num ) {
case RENDERING_NLIGHTS:
    AddVector3Md ( &light_centre, &edshapef_sectiondir, light_vect_factor, &p );
    CameraProjectPoint3d ( CPos, &p, &q );
    q.x = (double)x;  q.y = (double)y;
    CameraUnProjectPoint3d ( CPos, &q, &p );
    SubtractPoints3d ( &p, &light_centre, &q );
    NormalizeVector3d ( &q );
    edshapef_sectiondir = q;
    break;
case RENDERING_NLIGHTS+1:
case RENDERING_NLIGHTS+2:
case RENDERING_NLIGHTS+3:
    i = light_num - (RENDERING_NLIGHTS+1);
    CameraProjectPoint3d ( CPos, &edshapef_reflection_frame[i], &q );
    q.x = (double)x;  q.y = (double)y;
    CameraUnProjectPoint3d ( CPos, &q, &p );
    if ( i == 0 ) {
      for ( j = 1; j < 3; j++ ) {
        SubtractPoints3d ( &edshapef_reflection_frame[j],
                           &edshapef_reflection_frame[0], &q );
        AddVector3d ( &p, &q, &edshapef_reflection_frame[j] );
      }
    }
    edshapef_reflection_frame[i] = p;
    break;
case RENDERING_NLIGHTS+4:
case RENDERING_NLIGHTS+5:
case RENDERING_NLIGHTS+6:
    i = light_num - (RENDERING_NLIGHTS+4);
    CameraProjectPoint3d ( CPos, &edshapef_highlight_frame[i], &q );
    q.x = (double)x;  q.y = (double)y;
    CameraUnProjectPoint3d ( CPos, &q, &p );
    if ( i == 0 ) {
      for ( j = 1; j < 3; j++ ) {
        SubtractPoints3d ( &edshapef_highlight_frame[j],
                           &edshapef_highlight_frame[0], &q );
        AddVector3d ( &p, &q, &edshapef_highlight_frame[j] );
      }
    }
    edshapef_highlight_frame[i] = p;
    break;
default:
    break;
  }
} /*EdShapeFuncSetPoint*/

void EdShapeFuncSavePoints ( void )
{
} /*EdShapeFuncSavePoints*/

void EdShapeFuncTransformPoints ( trans3d *tr )
{
} /*EdShapeFuncTransformPoints*/

void EdShapeFuncDisplay ( xge_3Dwind *ww )
{
  point3d p;

  if ( sw_sections_frame ) {
    EdLightComputeVFactor ( ww );
    glPointSize ( 5.0 );
    glBegin ( GL_POINTS );
      glColor3fv ( xglec_DeepPink1 );
      glVertex3dv ( &light_centre.x );
      AddVector3Md ( &light_centre, &edshapef_sectiondir, light_vect_factor, &p );
      glVertex3dv ( &p.x );
    glEnd ();
    glBegin ( GL_LINES );
      glColor3fv ( xglec_DeepPink4 );
      glVertex3dv ( &light_centre.x );
      glVertex3dv ( &p.x );
    glEnd ();
  }
  if ( sw_reflection_frame ) {
    glPointSize ( 5.0 );
    glBegin ( GL_POINTS );
      glColor3fv ( xglec_DarkGoldenrod1 );
      glVertex3dv ( &edshapef_reflection_frame[1].x );
      glVertex3dv ( &edshapef_reflection_frame[0].x );
      glVertex3dv ( &edshapef_reflection_frame[2].x );
    glEnd ();
    glBegin ( GL_LINE_STRIP );
      glColor3fv ( xglec_DarkGoldenrod3 );
      glVertex3dv ( &edshapef_reflection_frame[1].x );
      glVertex3dv ( &edshapef_reflection_frame[0].x );
      glVertex3dv ( &edshapef_reflection_frame[2].x );
    glEnd ();
  }
  if ( sw_highlight_frame ) {
    glPointSize ( 5.0 );
    glBegin ( GL_POINTS );
      glColor3fv ( xglec_OliveDrab1 );
      glVertex3dv ( &edshapef_highlight_frame[1].x );
      glVertex3dv ( &edshapef_highlight_frame[0].x );
      glVertex3dv ( &edshapef_highlight_frame[2].x );
    glEnd ();
    glBegin ( GL_LINE_STRIP );
      glColor3fv ( xglec_OliveDrab4 );
      glVertex3dv ( &edshapef_highlight_frame[1].x );
      glVertex3dv ( &edshapef_highlight_frame[0].x );
      glVertex3dv ( &edshapef_highlight_frame[2].x );
    glEnd ();
  }
} /*EdShapeFuncDisplay*/

/* ////////////////////////////////////////////////////////////////////////// */
void DisplaySpecial3DElements ( xge_3Dwind *ww )
{
  if ( editing_lights )
    EdLightDisplay ( &g00win3D );
  if ( editing_shapefunc )
    EdShapeFuncDisplay ( &g00win3D );
} /*DisplaySpecial3DElements*/

