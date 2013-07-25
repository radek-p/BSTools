
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <GL/gl.h>
#include <GL/glu.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "glutki.h"
#include "linkage.h"
#include "palm.h"
#include "texture.h"
#include "ludzik.h"


/* dla kazdego parametru wartosc minimalna, maksymalna i domyslna z [0,1] */
const double art_par_desc[N_PARAMS][3] =
  {{1.0,-1.0,0.5},            /*PAR_GLOWA_1*/
   {-1.0,1.0,0.5},            /*PAR_GLOWA_2*/
   {-0.3,0.3,0.5},            /*PAR_GLOWA_3*/
   {0.7*PI,-0.5*PI,0.16667},  /*PAR_LBARK_1*/
   { 0.0,0.95*PI,0.0},        /*PAR_LBARK_2*/
   {-1.0,1.0,0.5},            /*PAR_LBARK_3*/
   {0.0,-0.85*PI,0.0},        /*PAR_LLOK_1*/ 
   {-0.5*PI,0.5*PI,0.5},      /*PAR_LLOK_2*/ 
   {0.7*PI,-0.5*PI,0.16667},  /*PAR_PBARK_1*/
   { 0.0,-0.95*PI,0.0},       /*PAR_PBARK_2*/
   {1.0,-1.0,0.5},            /*PAR_PBARK_3*/
   {0.0,-0.85*PI,0.0},        /*PAR_PLOK_1*/ 
   {0.5*PI,-0.5*PI,0.5},      /*PAR_PLOK_2*/ 
   {0.7*PI,-0.3*PI,0.2},      /*PAR_LBIO_1*/ 
   {0.0,0.5*PI,0.0},          /*PAR_LBIO_2*/ 
   {-1.0,1.0,0.5},            /*PAR_LBIO_3*/ 
   {0.0,0.75*PI,0.0},         /*PAR_LKOL_1*/ 
   {0.7*PI,-0.3*PI,0.2},      /*PAR_PBIO_1*/ 
   {0.0,-0.5*PI,0.0},         /*PAR_PBIO_2*/ 
   {-1.0,1.0,0.5},            /*PAR_PBIO_3*/ 
   {0.0,0.75*PI,0.0}};        /*PAR_PKOL_1*/
/* parametry artykulacji przed skalowaniem */
double art_param[N_PARAMS] =
  {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
   0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};

boolean sw_spline_hands = false;
boolean sw_draw_bsnets = false;

static kl_joint *root = NULL;
static GLUquadricObj *thecylinder = NULL;
static GLUnurbsObj *thenurbs = NULL;


void ResetParameters ( void )
{
  int i;

  for ( i = 0; i < N_PARAMS; i++ )
    art_param[i] = art_par_desc[i][2];
} /*ResetParameters*/

double ParamValue ( int par )
{
  return art_par_desc[par][0] +
         art_param[par]*(art_par_desc[par][1]-art_par_desc[par][0]);
} /*ParamValue*/

void trtogl ( trans3d *tr, double t[16] )
{
  t[0] = tr->U0.a11;  t[1] = tr->U0.a21;  t[2] = tr->U0.a31;   t[3] = 0.0;
  t[4] = tr->U0.a12;  t[5] = tr->U0.a22;  t[6] = tr->U0.a32;   t[7] = 0.0;
  t[8] = tr->U0.a13;  t[9] = tr->U0.a23;  t[10] = tr->U0.a33;  t[11] = 0.0;
  t[12] = tr->U0.a14; t[13] = tr->U0.a24; t[14] = tr->U0.a34;  t[15] = 1.0;
} /*trtogl*/

void RysujTulow ( kl_link *lnk, trans3d *tr )
{
  double trT[16];
  double tex_cs[4] = {1.0,0.0,0.0,0.0};
  double tex_ct[4] = {0.0,1.0,0.0,0.0};
  double a;

  glPushMatrix ();
  trtogl ( tr, trT );
  glMultMatrixd ( trT );
  glPushMatrix ();

  glMaterialf ( GL_FRONT_AND_BACK, GL_SHININESS, 100.0 );
  glEnable ( GL_TEXTURE_2D );
  glTexEnvi ( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
  glBindTexture ( GL_TEXTURE_2D, texName );
  glTexGeni ( GL_S, GL_TEXTURE_GEN_MODE, GL_OBJECT_LINEAR );
  glTexGendv ( GL_S, GL_OBJECT_PLANE, tex_cs );
  glTexGeni ( GL_T, GL_TEXTURE_GEN_MODE, GL_OBJECT_LINEAR );
  glTexGendv ( GL_T, GL_OBJECT_PLANE, tex_ct );
  glEnable ( GL_TEXTURE_GEN_S );
  glEnable ( GL_TEXTURE_GEN_T );

  glMatrixMode ( GL_TEXTURE );
  glLoadIdentity ();
  glMultMatrixd ( trT );
  a = 1.0;
  glScaled ( a*1.1, a*2.0, a*0.5 );
  glMatrixMode ( GL_MODELVIEW );

  glScaled ( 1.1, 2.0, 0.5 );
/*  glColor3f ( 1.0, 0.0, 0.0 ); */
  glColor3f ( 1.0, 1.0, 1.0 );
  glutSolidCube ( 2.0 );
  glPopMatrix ();

/* odtad */
/*
  glPushMatrix ();
  glTranslated ( 0.5, 1.0, 0.28 );
  glutSolidSphere ( 0.55, 16, 32 );
  glPopMatrix ();
  glPushMatrix ();
  glTranslated ( -0.5, 1.0, 0.28 );
  glutSolidSphere ( 0.55, 16, 32 );
  glPopMatrix ();
*/
/* dotad */
  glBindTexture ( GL_TEXTURE_2D, 0 );
  glDisable ( GL_TEXTURE_2D );

  glColor3f ( 0.8, 0.8, 0.9 );
  glPushMatrix ();
  glTranslated ( 1.3, 1.6, 0.0 );
  glutSolidSphere ( 0.4, 16, 32 );
  glPopMatrix ();
  glPushMatrix ();
  glTranslated ( -1.3, 1.6, 0.0 );
  glutSolidSphere ( 0.4, 16, 32 );
  glPopMatrix ();
  glPushMatrix ();
  glTranslated ( 0.6, -2.0, 0.0 );
  glutSolidSphere ( 0.45, 16, 32 );
  glPopMatrix ();
  glPushMatrix ();
  glTranslated ( -0.6, -2.0, 0.0 );
  glutSolidSphere ( 0.45, 16, 32 );
  glRotated ( 90, 0.0, 1.0, 0.0 );
  gluCylinder ( thecylinder, 0.45, 0.45, 1.2, 16, 1 );
  glPopMatrix ();
  glPopMatrix ();
} /*RysujTulow*/

void ArticulateRoot ( kl_joint *jnt )
{
} /*ArticulateRoot*/

void NieRysuj ( kl_link *lnk, trans3d *tr )
{
  /* nic nie rysuj - to jest procedura  */
  /* dla wszystkich czlonow wirtualnych */
} /*NieRysuj*/

void RysujGlowe ( kl_link *lnk, trans3d *tr )
{
  double trT[16];
  GLfloat specular[] = {0.9,0.8,0.6,1.0};

  glMaterialfv ( GL_FRONT_AND_BACK, GL_SPECULAR, specular );
  glMaterialf ( GL_FRONT_AND_BACK, GL_SHININESS, 100.0 );
  glPushMatrix ();
  trtogl ( tr, trT );
  glMultMatrixd ( trT );
  glColor3f ( 0.9, 0.8, 0.6 );
  glutSolidTeapot ( 0.9 );
  glPushMatrix ();
  glTranslated ( 0.58, 0.14, 0.23 );
  glColor3f ( 1.0, 1.0, 1.0 );
  glutSolidSphere ( 0.2, 16, 32 );
  glTranslated ( 0.1, 0.02, 0.0 );
  glColor3f ( 0.0, 0.0, 0.0 );
  glutSolidSphere ( 0.115, 16, 32 );
  glPopMatrix ();
  glTranslated ( 0.58, 0.14, -0.23 );
  glColor3f ( 1.0, 1.0, 1.0 );
  glutSolidSphere ( 0.2, 16, 32 );
  glTranslated ( 0.1, 0.02, 0.0 );
  glColor3f ( 0.0, 0.0, 0.0 );
  glutSolidSphere ( 0.115, 16, 32 );
  glPopMatrix ();
} /*RysujGlowe*/

void ArtykulujGlowe ( kl_joint *jnt )
{
  IdentTrans3d ( &jnt->R );
  switch ( jnt->id ) {
case PAR_GLOWA_1:
    RotXTrans3d ( &jnt->R, ParamValue ( PAR_GLOWA_1 ) );
    break;
case PAR_GLOWA_2:
    RotYTrans3d ( &jnt->R, ParamValue ( PAR_GLOWA_2 ) );
    break;
case PAR_GLOWA_3:
    RotZTrans3d ( &jnt->R, ParamValue ( PAR_GLOWA_3 ) );
    break;
default:
    break;
  }
} /*ArtykulujGlowe*/

/* now the hand - this is a bit tricky, as we need to assemble */
/* the B-spline control net using points given in different coordinate systems */
#define HDEGU 3 
#define HDEGV 3 
#define HLKNU 10
#define HLKNV 10
float hknu[HLKNU+1] = { 0.0, 0.0, 0.0, 0.0, 0.16, 0.33, 0.5, 1.0, 1.0, 1.0, 1.0 };
float hknv[HLKNV+1] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
/* the control points will be computed based on the articulation parameters */
point3f hcp[(HLKNU-HDEGU)*(HLKNV-HDEGV)];
point3f hgen[4] =
  {{-0.15,-0.15,0.0},{-0.15,0.15,0.0},{0.15,0.15,0.0},{0.15,-0.15,0.0}};
trans3d hbtr, hbtri;

void DrawBSNet ( int n, int lknu, int m, int lknv, point3f *cp )
{
  int i, j, k;
  boolean len;

  len = glIsEnabled ( GL_LIGHTING );
  glDisable ( GL_LIGHTING );
  lknu -= n;
  lknv -= m;
  glColor3f ( 0.0, 1.0, 0.0 );
  for ( i = k = 0;  i < lknu;  i++, k++ ) {
    glBegin ( GL_LINE_STRIP );
      for ( j = 0;  j < lknv-1;  j++, k++ )
        glVertex3fv ( &cp[k].x );
    glEnd ();
  }
  for ( i = 0; i < lknv; i++ ) {
    glBegin ( GL_LINE_STRIP );
      for ( j = 0, k = i;  j < lknu;  j++, k += lknv )
        glVertex3fv ( &cp[k].x );
    glEnd ();
  }
  if ( len )
    glEnable ( GL_LIGHTING );
} /*DrawBSNet*/

void EnterHColumn ( int col, trans3d *tr )
{
  int     i;
  trans3d t;
  double  trT[16];

  if ( col < 0 || col > HLKNU-HDEGU )
    return;
  if ( col == 0 ) {
    memcpy ( &hcp[0], &hgen[0], 4*sizeof(point3f) );
    memcpy ( &hcp[4], &hgen[0], 3*sizeof(point3f) ); 
    hbtr = hbtri = *tr;
    InvertTrans3d ( &hbtri );
  }
  else {
    CompTrans3d ( &t, &hbtri, tr );
        /* HLKNV-HDEGV == 7 == 4+3 */
    for ( i = 0; i < 4; i++ )  
      TransPoint3df ( &t, &hgen[i], &hcp[(HLKNV-HDEGV)*col+i] );
    memcpy ( &hcp[(HLKNV-HDEGV)*col+4], &hcp[(HLKNV-HDEGV)*col], 3*sizeof(point3f) );

    if ( col == HLKNU-HDEGU-1 ) {
      glPushMatrix ();
      trtogl ( &hbtr, trT );
      glMultMatrixd ( trT );
      glEnable ( GL_AUTO_NORMAL );
      gluNurbsProperty ( thenurbs, GLU_SAMPLING_METHOD, GLU_PATH_LENGTH );
      gluNurbsProperty ( thenurbs, GLU_SAMPLING_TOLERANCE, 10.0 );
      gluBeginSurface ( thenurbs );
      gluNurbsSurface ( thenurbs, HLKNU+1, hknu, HLKNV+1, hknv,
                        3*(HLKNV-HDEGV), 3, &hcp[0].x, HDEGU+1, HDEGV+1,
                        GL_MAP2_NORMAL );
      gluNurbsSurface ( thenurbs, HLKNU+1, hknu, HLKNV+1, hknv,
                        3*(HLKNV-HDEGV), 3, &hcp[0].x, HDEGU+1, HDEGV+1,
                        GL_MAP2_VERTEX_3 );
      gluEndSurface ( thenurbs );
      glDisable ( GL_AUTO_NORMAL );
      if ( sw_draw_bsnets )
        DrawBSNet ( HDEGU, HLKNU, HDEGV, HLKNV, hcp );
      glPopMatrix ();
    }
  }  
} /*EnterHColumn*/

 void RysujLeweRamie ( kl_link *lnk, trans3d *tr )
{
  trans3d t;
  double  trT[16], s, v;

  glMaterialf ( GL_FRONT_AND_BACK, GL_SHININESS, 100.0 );
  if ( sw_spline_hands ) {
    s = 1.0 + 0.5*art_param[PAR_LLOK_1];
    v = 0.05*art_param[PAR_LLOK_1];
    EnterHColumn ( 0, tr );
    t = *tr;
    Trans3Scaled ( &t, s, s, 1.0 );
    Trans3Shiftd ( &t, 0.0, v, 0.6 );
    EnterHColumn ( 1, &t );
    t = *tr;
    Trans3Scaled ( &t, s, s, 1.0 );
    Trans3Shiftd ( &t, 0.0, 0.5*v, 1.2 );
    EnterHColumn ( 2, &t );
    t = *tr;
    Trans3Shiftd ( &t, 0.0, 0.0, 1.6 ); 
    EnterHColumn ( 3, &t );
    t = *tr;
    Trans3Shiftd ( &t, 0.0, 0.0, 2.05 );
    Trans3RotXd ( &t, 0.5*ParamValue ( PAR_LLOK_1 ) );
    Trans3Scaled ( &t, 1.0, 2.0/(2.0-art_param[PAR_LLOK_1]), 1.0 );
    EnterHColumn ( 4, &t );
  }
  else {
    glPushMatrix ();
    trtogl ( tr, trT );
    glMultMatrixd ( trT );
    glColor3f ( 0.4, 0.4, 1.0 );
    gluCylinder ( thecylinder, 0.15, 0.15, 2.05, 10, 1 );
    glTranslated ( 0.0, 0.0, 2.05 );
    glColor3f ( 0.8, 0.8, 0.9 );
    glutSolidSphere ( 0.2, 16, 32 );
    glPopMatrix ();
  }
} /*RysujLeweRamie*/

void ArtykulujLeweRamie ( kl_joint *jnt )
{
  IdentTrans3d ( &jnt->R );
  switch ( jnt->id ) {
case PAR_LBARK_1:
    RotXTrans3d ( &jnt->R, ParamValue ( PAR_LBARK_1 ) );
    break;
case PAR_LBARK_2:
    RotYTrans3d ( &jnt->R, ParamValue ( PAR_LBARK_2 ) );
    break;
case PAR_LBARK_3:
    RotZTrans3d ( &jnt->R, ParamValue ( PAR_LBARK_3 ) );
    break;
default:
    break;
  }
} /*ArtykulujLeweRamie*/

void RysujLewePrzedramie ( kl_link *lnk, trans3d *tr )
{
  trans3d t;
  double  trT[16];

  glPushMatrix ();
  glColor3f ( 0.4, 0.4, 1.0 );
  if ( sw_spline_hands ) {
    EnterHColumn ( 5, tr );
    t = *tr;
    Trans3Scaled ( &t, 0.55, 0.55, 1.0 );   
    Trans3Shiftd ( &t, 0.0, 0.0, 1.9 ); 
    EnterHColumn ( 6, &t );
    trtogl ( tr, trT );
    glMultMatrixd ( trT );
  }
  else {
    trtogl ( tr, trT );
    glMultMatrixd ( trT );
    gluCylinder ( thecylinder, 0.15, 0.07, 1.9, 10, 1 );
  }
  glTranslated ( 0.0, 0.0, 1.8 );
  DrawLeftPalm ();
  glPopMatrix ();
} /*RysujLewePrzedramie*/

void ArtykulujLewePrzedramie ( kl_joint *jnt )
{
  IdentTrans3d ( &jnt->R );
  switch ( jnt->id ) {
case PAR_LLOK_1:
    RotXTrans3d ( &jnt->R, ParamValue ( PAR_LLOK_1 ) );
    break;
case PAR_LLOK_2:
    RotZTrans3d ( &jnt->R, ParamValue ( PAR_LLOK_2 ) );
    break;
default:
    break;
  }
} /*ArtykulujLewePrzedramie*/

void RysujPraweRamie ( kl_link *lnk, trans3d *tr )
{
  trans3d t;
  double  trT[16], s, v;

  glMaterialf ( GL_FRONT_AND_BACK, GL_SHININESS, 100.0 );
  if ( sw_spline_hands ) {
    s = 1.0 + 0.5*art_param[PAR_PLOK_1];
    v = 0.05*art_param[PAR_PLOK_1];
    EnterHColumn ( 0, tr );
    t = *tr;
    Trans3Scaled ( &t, s, s, 1.0 );
    Trans3Shiftd ( &t, 0.0, v, 0.6 );
    EnterHColumn ( 1, &t );
    t = *tr;
    Trans3Scaled ( &t, s, s, 1.0 );
    Trans3Shiftd ( &t, 0.0, 0.5*v, 1.2 );
    EnterHColumn ( 2, &t );
    t = *tr;
    Trans3Shiftd ( &t, 0.0, 0.0, 1.6 );
    EnterHColumn ( 3, &t );
    t = *tr;
    Trans3Shiftd ( &t, 0.0, 0.0, 2.05 );
    Trans3RotXd ( &t, 0.5*ParamValue ( PAR_PLOK_1 ) );
    Trans3Scaled ( &t, 1.0, 2.0/(2.0-art_param[PAR_PLOK_1]), 1.0 );
    EnterHColumn ( 4, &t );
  }
  else {
    glPushMatrix ();
    trtogl ( tr, trT );
    glMultMatrixd ( trT );
    glColor3f ( 0.4, 0.4, 1.0 );
    gluCylinder ( thecylinder, 0.15, 0.15, 2.05, 10, 1 );
    glTranslated ( 0.0, 0.0, 2.05 );
    glColor3f ( 0.8, 0.8, 0.9 );
    glutSolidSphere ( 0.2, 16, 32 );
    glPopMatrix ();
  }
} /*RysujPraweRamie*/

void ArtykulujPraweRamie ( kl_joint *jnt )
{
  IdentTrans3d ( &jnt->R );
  switch ( jnt->id ) {
case PAR_PBARK_1:
    RotXTrans3d ( &jnt->R, ParamValue ( PAR_PBARK_1 ) );
    break;
case PAR_PBARK_2:
    RotYTrans3d ( &jnt->R, ParamValue ( PAR_PBARK_2 ) );
    break;
case PAR_PBARK_3:
    RotZTrans3d ( &jnt->R, ParamValue ( PAR_PBARK_3 ) );
    break;
default:
    break;
  }
} /*ArtykulujPraweRamie*/

void RysujPrawePrzedramie ( kl_link *lnk, trans3d *tr )
{
  trans3d t;
  double  trT[16];

  glPushMatrix ();
  glColor3f ( 0.4, 0.4, 1.0 );
  if ( sw_spline_hands ) {
    EnterHColumn ( 5, tr );
    t = *tr;
    Trans3Scaled ( &t, 0.55, 0.55, 1.0 );
    Trans3Shiftd ( &t, 0.0, 0.0, 1.9 );  
    EnterHColumn ( 6, &t );
    trtogl ( tr, trT );
    glMultMatrixd ( trT );
  }
  else {
    trtogl ( tr, trT );
    glMultMatrixd ( trT );
    gluCylinder ( thecylinder, 0.15, 0.07, 1.9, 10, 1 );
  }
  glTranslated ( 0.0, 0.0, 1.8 );
  DrawRightPalm ();
  glPopMatrix ();
} /*RysujPrawePrzedramie*/

void ArtykulujPrawePrzedramie ( kl_joint *jnt )
{
  IdentTrans3d ( &jnt->R );
  switch ( jnt->id ) {
case PAR_PLOK_1:
    RotXTrans3d ( &jnt->R, ParamValue ( PAR_PLOK_1 ) );
    break;
case PAR_PLOK_2:
    RotZTrans3d ( &jnt->R, ParamValue ( PAR_PLOK_2 ) );
    break;
default:  
    break;
  }
} /*ArtykulujPrawePrzedramie*/

void RysujLeweUdo ( kl_link *lnk, trans3d *tr )
{
  double trT[16];

  glMaterialf ( GL_FRONT_AND_BACK, GL_SHININESS, 100.0 );
  glPushMatrix ();
  trtogl ( tr, trT );
  glMultMatrixd ( trT );
  glColor3f ( 0.4, 0.4, 1.0 );
  gluCylinder ( thecylinder, 0.25, 0.2, 2.1, 10, 1 );
  glTranslated ( 0.0, 0.0, 2.1 );
  glColor3f ( 0.8, 0.8, 0.9 );   
  glutSolidSphere ( 0.24, 16, 32 );
  glPopMatrix ();
} /*RysujLeweUdo*/

void ArtykulujLeweUdo ( kl_joint *jnt )
{
  IdentTrans3d ( &jnt->R );
  switch ( jnt->id ) {
case PAR_LBIO_1:
    RotXTrans3d ( &jnt->R, ParamValue ( PAR_LBIO_1 ) );
    break;
case PAR_LBIO_2:
    RotYTrans3d ( &jnt->R, ParamValue ( PAR_LBIO_2 ) );
    break;
case PAR_LBIO_3:
    RotZTrans3d ( &jnt->R, ParamValue ( PAR_LBIO_3 ) );
    break;
default:  
    break;
  }
} /*ArtykulujLeweUdo*/

void RysujBut ( void )
{
  static double pt[10][3] =
    {{-0.2,0.7,0.0},{-0.2,0.7,-0.1},{-0.2,0.1,-0.3},{-0.2,-0.2,-0.3},{-0.2,-0.2,0.0},
     {+0.2,0.7,0.0},{+0.2,0.7,-0.1},{+0.2,0.1,-0.3},{+0.2,-0.2,-0.3},{+0.2,-0.2,0.0}};
  static vector3d nv[7];
  static boolean nvok = false;

  vector3d a, b;
  int i;

#define COMPNV(j,k,l,m) { \
    SubtractPoints3d ( (point3d*)pt[k], (point3d*)pt[j], &a ); \
    SubtractPoints3d ( (point3d*)pt[l], (point3d*)pt[k], &b ); \
    CrossProduct3d ( &b, &a, &nv[m] ); \
    NormalizeVector3d ( &nv[m] ); \
  }

  if ( !nvok ) {
    COMPNV ( 6, 5, 0, 0 )
    COMPNV ( 6, 1, 2, 1 )
    COMPNV ( 8, 7, 2, 2 )
    COMPNV ( 8, 3, 4, 3 )
    COMPNV ( 5, 9, 4, 4 )
    COMPNV ( 2, 1, 0, 5 )
    COMPNV ( 5, 6, 7, 6 )
    nvok = true;
  }
  glColor3f ( 0.8, 0.8, 0.8 );
  glBegin ( GL_QUADS );
    glNormal3dv ( &nv[0].x );
    glVertex3dv ( pt[0] );
    glVertex3dv ( pt[5] );
    glVertex3dv ( pt[6] );
    glVertex3dv ( pt[1] );

    glNormal3dv ( &nv[1].x );
    glVertex3dv ( pt[6] );
    glVertex3dv ( pt[1] );
    glVertex3dv ( pt[2] );
    glVertex3dv ( pt[7] );

    glNormal3dv ( &nv[2].x );
    glVertex3dv ( pt[2] );
    glVertex3dv ( pt[7] );
    glVertex3dv ( pt[8] );
    glVertex3dv ( pt[3] );

    glNormal3dv ( &nv[3].x );
    glVertex3dv ( pt[8] );
    glVertex3dv ( pt[3] );
    glVertex3dv ( pt[4] );
    glVertex3dv ( pt[9] );

    glNormal3dv ( &nv[4].x );
    glVertex3dv ( pt[4] );
    glVertex3dv ( pt[9] );
    glVertex3dv ( pt[5] );
    glVertex3dv ( pt[0] );
  glEnd ();
  glBegin ( GL_POLYGON );
    glNormal3dv ( &nv[5].x );
    for ( i = 0; i < 5; i++ )
      glVertex3dv ( pt[i] ); 
  glEnd ();
  glBegin ( GL_POLYGON );
    glNormal3dv ( &nv[6].x );
    for ( i = 5; i < 10; i++ )
      glVertex3dv ( pt[i] );  
  glEnd ();
#undef COMPNV
} /*RysujBut*/

void RysujLewyGolen ( kl_link *lnk, trans3d *tr )
{
  double trT[16];

  glMaterialf ( GL_FRONT_AND_BACK, GL_SHININESS, 100.0 );
  glPushMatrix ();
  trtogl ( tr, trT );
  glMultMatrixd ( trT );
  glColor3f ( 0.4, 0.4, 1.0 );
  gluCylinder ( thecylinder, 0.22, 0.16, 2.0, 10, 1 );
  glTranslated ( 0.0, 0.0, 2.05 );
/*  RysujBut (); */
  DrawLeftFoot ();
  glPopMatrix ();
} /*RysujLewyGolen*/

void ArtykulujLewyGolen ( kl_joint *jnt )
{
  IdentTrans3d ( &jnt->R );
  switch ( jnt->id ) {
case PAR_LKOL_1:
    RotXTrans3d ( &jnt->R, ParamValue ( PAR_LKOL_1 ) );
    break;
default:  
    break;
  }
} /*ArtykulujLewyGolen*/

void RysujPraweUdo ( kl_link *lnk, trans3d *tr )
{
  double trT[16];

  glMaterialf ( GL_FRONT_AND_BACK, GL_SHININESS, 100.0 );
  glPushMatrix ();
  trtogl ( tr, trT );
  glMultMatrixd ( trT );
  glColor3f ( 0.4, 0.4, 1.0 );
  gluCylinder ( thecylinder, 0.25, 0.2, 2.1, 10, 1 );
  glTranslated ( 0.0, 0.0, 2.1 );
  glColor3f ( 0.8, 0.8, 0.9 );   
  glutSolidSphere ( 0.24, 16, 32 );
  glPopMatrix ();
} /*RysujPraweUdo*/

void ArtykulujPraweUdo ( kl_joint *jnt )
{
  IdentTrans3d ( &jnt->R );
  switch ( jnt->id ) {
case PAR_PBIO_1:
    RotXTrans3d ( &jnt->R, ParamValue ( PAR_PBIO_1 ) );
    break;
case PAR_PBIO_2:
    RotYTrans3d ( &jnt->R, ParamValue ( PAR_PBIO_2 ) );
    break;
case PAR_PBIO_3:
    RotZTrans3d ( &jnt->R, ParamValue ( PAR_PBIO_3 ) );
    break;
default:  
    break;
  }
} /*ArtykulujPraweUdo*/

void RysujPrawyGolen ( kl_link *lnk, trans3d *tr )
{
  double trT[16];

  glMaterialf ( GL_FRONT_AND_BACK, GL_SHININESS, 100.0 );
  glPushMatrix ();
  trtogl ( tr, trT );
  glMultMatrixd ( trT );
  glColor3f ( 0.4, 0.4, 1.0 );
  gluCylinder ( thecylinder, 0.22, 0.16, 2.0, 10, 1 );
  glTranslated ( 0.0, 0.0, 2.05 );
/*  RysujBut (); */
  DrawRightFoot ();
  glPopMatrix ();
} /*RysujPrawyGolen*/

void ArtykulujPrawyGolen ( kl_joint *jnt )
{
  IdentTrans3d ( &jnt->R );
  switch ( jnt->id ) {
case PAR_PKOL_1:
    RotXTrans3d ( &jnt->R, ParamValue ( PAR_PKOL_1 ) );
    break;
default:  
    break;
  }
} /*ArtykulujPrawyGolen*/

void InitCharacter ( void )
{
  kl_joint *jt;
  kl_link  *lk0, *lk1, *lk2;

  ResetParameters ();
  thecylinder = gluNewQuadric ();
  thenurbs = gluNewNurbsRenderer ();
  if ( !thecylinder || !thenurbs ) exit ( 1 );

  InitPalm ();
  InitFoot ();
  InitTexture ();

  /* tulow */
  lk2 = NewLink ( 1, RysujTulow );
  if ( !lk2 ) exit ( 1 );
  root = NewJoint ( 1, NULL, lk2, ArticulateRoot );
  if ( !root ) exit ( 1 );
  lk0 = lk1 = lk2;
  /* glowa */
  lk2 = NewLink ( 2, NieRysuj );
  if ( !lk2 ) exit ( 1 );
  jt = NewJoint ( PAR_GLOWA_1, lk1, lk2, ArtykulujGlowe );
  if ( !jt ) exit ( 1 );
  ShiftTrans3d ( &jt->L, 0.0, 2.1, 0.0 );
  lk1 = lk2;
  lk2 = NewLink ( 3, NieRysuj );
  if ( !lk2 ) exit ( 1 );
  jt = NewJoint ( PAR_GLOWA_2, lk1, lk2, ArtykulujGlowe );
  if ( !jt ) exit ( 1 );
  lk1 = lk2;
  lk2 = NewLink ( 4, RysujGlowe );
  if ( !lk2 ) exit ( 1 );
  RotYTrans3d ( &lk2->M, -0.5*PI );
  ShiftTrans3d ( &lk2->M, 0.0, 0.68, 0.0 );
  jt = NewJoint ( PAR_GLOWA_3, lk1, lk2, ArtykulujGlowe );
  if ( !jt ) exit ( 1 );
  /* lewa reka */
  lk2 = NewLink ( 5, NieRysuj );
  if ( !lk2 ) exit ( 1 );
  jt = NewJoint ( PAR_LBARK_1, lk0, lk2, ArtykulujLeweRamie );
  if ( !jt ) exit ( 1 );
  ShiftTrans3d ( &jt->L, 1.35, 1.6, 0.0 );
  lk1 = lk2;
  lk2 = NewLink ( 6, NieRysuj );
  if ( !lk2 ) exit ( 1 );
  jt = NewJoint ( PAR_LBARK_2, lk1, lk2, ArtykulujLeweRamie );
  if ( !jt ) exit ( 1 );
  lk1 = lk2;
  lk2 = NewLink ( 7, RysujLeweRamie );
  if ( !lk2 ) exit ( 1 );
  jt = NewJoint ( PAR_LBARK_3, lk1, lk2, ArtykulujLeweRamie );
  if ( !jt ) exit ( 1 );
  lk1 = lk2;
  lk2 = NewLink ( 8, NieRysuj );
  if ( !lk2 ) exit ( 1 );
  jt = NewJoint ( PAR_LLOK_1, lk1, lk2, ArtykulujLewePrzedramie );
  ShiftTrans3d ( &jt->L, 0.0, 0.0, 2.05 );
  if ( !jt ) exit ( 1 );
  lk1 = lk2;
  lk2 = NewLink ( 9, RysujLewePrzedramie );
  if ( !lk2 ) exit ( 1 );
  jt = NewJoint ( PAR_LLOK_2, lk1, lk2, ArtykulujLewePrzedramie );
  if ( !jt ) exit ( 1 );
  /* prawa reka */
  lk2 = NewLink ( 5, NieRysuj );
  if ( !lk2 ) exit ( 1 );
  jt = NewJoint ( PAR_PBARK_1, lk0, lk2, ArtykulujPraweRamie );
  if ( !jt ) exit ( 1 );
  ShiftTrans3d ( &jt->L, -1.35, 1.6, 0.0 );
  lk1 = lk2;
  lk2 = NewLink ( 6, NieRysuj );
  if ( !lk2 ) exit ( 1 );
  jt = NewJoint ( PAR_PBARK_2, lk1, lk2, ArtykulujPraweRamie );
  if ( !jt ) exit ( 1 );
  lk1 = lk2;
  lk2 = NewLink ( 7, RysujPraweRamie );
  if ( !lk2 ) exit ( 1 );
  jt = NewJoint ( PAR_PBARK_3, lk1, lk2, ArtykulujPraweRamie );
  if ( !jt ) exit ( 1 );
  lk1 = lk2;
  lk2 = NewLink ( 8, NieRysuj );
  if ( !lk2 ) exit ( 1 );
  jt = NewJoint ( PAR_PLOK_1, lk1, lk2, ArtykulujPrawePrzedramie );
  ShiftTrans3d ( &jt->L, 0.0, 0.0, 2.05 );
  if ( !jt ) exit ( 1 );
  lk1 = lk2;
  lk2 = NewLink ( 9, RysujPrawePrzedramie );
  if ( !lk2 ) exit ( 1 );
  jt = NewJoint ( PAR_PLOK_2, lk1, lk2, ArtykulujPrawePrzedramie );
  if ( !jt ) exit ( 1 );
  /* lewa noga */
  lk2 = NewLink ( 10, NieRysuj );
  if ( !lk2 ) exit ( 1 );
  jt = NewJoint ( PAR_LBIO_1, lk0, lk2, ArtykulujLeweUdo );
  if ( !jt ) exit ( 1 );
  ShiftTrans3d ( &jt->L, 0.6, -2.0, 0.0 );
  lk1 = lk2;
  lk2 = NewLink ( 11, NieRysuj );
  if ( !lk2 ) exit ( 1 );
  jt = NewJoint ( PAR_LBIO_2, lk1, lk2, ArtykulujLeweUdo );
  if ( !jt ) exit ( 1 );
  lk1 = lk2;
  lk2 = NewLink ( 12, RysujLeweUdo );
  if ( !lk2 ) exit ( 1 );
  jt = NewJoint ( PAR_LBIO_3, lk1, lk2, ArtykulujLeweUdo );
  if ( !jt ) exit ( 1 );
  lk1 = lk2;
  lk2 = NewLink ( 13, RysujLewyGolen );
  if ( !lk2 ) exit ( 1 );
  jt = NewJoint ( PAR_LKOL_1, lk1, lk2, ArtykulujLewyGolen );
  if ( !jt ) exit ( 1 );
  ShiftTrans3d ( &jt->L, 0.0, 0.0, 2.1 );
  /* prawa noga */
  lk2 = NewLink ( 14, NieRysuj );
  if ( !lk2 ) exit ( 1 );
  jt = NewJoint ( PAR_PBIO_1, lk0, lk2, ArtykulujPraweUdo );
  if ( !jt ) exit ( 1 );
  ShiftTrans3d ( &jt->L, -0.6, -2.0, 0.0 );
  lk1 = lk2;
  lk2 = NewLink ( 15, NieRysuj );
  if ( !lk2 ) exit ( 1 );
  jt = NewJoint ( PAR_PBIO_2, lk1, lk2, ArtykulujPraweUdo );
  if ( !jt ) exit ( 1 );
  lk1 = lk2;
  lk2 = NewLink ( 16, RysujLeweUdo );
  if ( !lk2 ) exit ( 1 );
  jt = NewJoint ( PAR_PBIO_3, lk1, lk2, ArtykulujPraweUdo );
  if ( !jt ) exit ( 1 );
  lk1 = lk2;
  lk2 = NewLink ( 13, RysujPrawyGolen );
  if ( !lk2 ) exit ( 1 );
  jt = NewJoint ( PAR_PKOL_1, lk1, lk2, ArtykulujPrawyGolen );
  if ( !jt ) exit ( 1 );
  ShiftTrans3d ( &jt->L, 0.0, 0.0, 2.1 );
} /*InitCharacter*/

void DisplayCharacter ( void )
{
  trans3d ident;

  if ( root ) {
    IdentTrans3d ( &ident );
    glEnable ( GL_NORMALIZE );
    DisplayLinkage ( root, &ident );
    glDisable ( GL_NORMALIZE );
  }
} /*DisplayCharacter*/

