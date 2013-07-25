
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */
/* this file was written by Tomasz Swierczek            */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <GL/gl.h>
#include <GL/glu.h>

#include "pkvaria.h"

#include "lighting.h"

#define FLOOR_TILES	100 


float   R;			/*jak daleko od środka układu współrzędnych znajduje się nasza latarka*/
float   RGBA[4];		/*składowe koloru światła*/
float   Con,Lin,Sq;		/*stałe zanikania od odegłości w modelu OpenGL; zanikanie = 1.0 / (Con + odległość * Lin + odległość * odległość * Sq )*/
float   SpotCutoff;		/*1/2 kąta rozwarcia stożka światełka*/
float   HorizAngle, VertAngle;	/*kąty pod jakimi patrzy nasze swiatło*/
boolean useSpecular = false;	/*flaga do włączania oświetlenia odbitego*/

void InitSpotlight ( void )
{	
  R = 20.0;
  RGBA[0] = RGBA[1] = RGBA[2] = RGBA[3] = 0.2;
  Con = 0.01;
  Lin = 0.01;
  Sq = 0.015;
  SpotCutoff = 90.0;
  HorizAngle = PI/4.0;
  VertAngle = -PI/4.0;
  useSpecular = 0;
} /*InitSpotlight*/

void MoveLight ( float dr )
{
  R += dr;
  R = R > 0.1 ? R : 0.1;
}

void ChangeLightColor ( float dr, float dg, float db )
{
  RGBA[0] += dr;
  RGBA[1] += dg;
  RGBA[2] += db;	
}

void SetLightAttenuation ( float con, float lin, float sq )
{
  Con = con;
  Lin = lin;
  Sq = sq;
}

void ChangeSpotCutOff ( float dangle )
{
  SpotCutoff += dangle;
}

void SetupLightToDisplay ( void )
{
  GLfloat temp[] = {0.2,0.2,0.2,1.0};
  int i;
  GLint lights;

  glEnable ( GL_LIGHTING );

    /* nie chcemy żadnego oświetlenia poza naszą latarką
       --- to zmienilem troche: chce troszeczke swiatla (P.K.) */
  glLightModelfv ( GL_LIGHT_MODEL_AMBIENT, temp );

	/* używanie kolorów z glColor() */
  glEnable ( GL_COLOR_MATERIAL );
  glColorMaterial ( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE );
	
  glGetIntegerv ( GL_MAX_LIGHTS, &lights );
  for( i = 1; i < lights; ++i)
    glDisable ( GL_LIGHT0+i );
  glEnable ( GL_LIGHT0 );

        /* glLightfv(GL_LIGHT0,GL_AMBIENT,RGBA);	
           komentujemy tę linię, bo nie chcemy (teraz) mieć żadnego oświetlenia na powierzchniach odwróconych tyłem	
        */
  temp[0] = temp[1] = temp[2] = 0.0;
  temp[3] = 1.0;
  glLightfv ( GL_LIGHT0, GL_AMBIENT, temp );
  glLightfv ( GL_LIGHT0, GL_DIFFUSE, RGBA );
  if ( useSpecular )
    glLightfv ( GL_LIGHT0, GL_SPECULAR, RGBA );
  else	
    glLightfv ( GL_LIGHT0, GL_SPECULAR, temp );
		
        /*ustawienie macierzy takie, aby uwzglednić położenie kątowe światła...*/
  glMatrixMode ( GL_MODELVIEW );
  glPushMatrix ();

  glRotatef ( 180.0*HorizAngle/PI, 0.0,1.0,0.0 );
  glRotatef ( 180.0*VertAngle/PI, 1.0,0.0,0.0 );

        /* pozycja światła */
  temp[0] = 0.0;
  temp[1] = 0.0;
  temp[2] = R;
  temp[3] = 1.0;
  glLightfv ( GL_LIGHT0,GL_POSITION,temp );
	
	/* kierunek świecenia latarki - do środka układu współrzędnych */
  temp[2] = -R;
  glLightfv ( GL_LIGHT0, GL_SPOT_DIRECTION, temp );

        /* kąt stożka (dokładniej jego 1/2) */
  glLightf ( GL_LIGHT0, GL_SPOT_CUTOFF, SpotCutoff );

        /*
        parametr modyfikujący skupienie światła w środku stożka (im więcej tym bardziej skupione światło)	
	*/
  glLightf ( GL_LIGHT0, GL_SPOT_EXPONENT, 1.0 );

        /* zanikanie */
  glLightf ( GL_LIGHT0, GL_CONSTANT_ATTENUATION, Con );
  glLightf ( GL_LIGHT0, GL_LINEAR_ATTENUATION, Lin );
  glLightf ( GL_LIGHT0, GL_QUADRATIC_ATTENUATION, Sq );

        /*przywrócenie poprzedniej macierzy widoku*/
  glPopMatrix ();
}

void DisableLightToDisplay ( void )
{
  glDisable(GL_LIGHTING);
  glDisable(GL_LIGHT0);
}


void RenderSimpleFloor ( float height, float size )
{
  int i,j;
  float TileSize;
  float x, z;
  GLfloat specular[] = {1.0,1.0,1.0,1.0};

  glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,specular);
  glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,100.0);

  TileSize = size / (float) FLOOR_TILES;
  x = -0.5 * size - 0.5 * TileSize;

  glColor4f(1.0,1.0,1.0,1.0);
  glNormal3f(0.0,1.0,0.0);	

  glBegin(GL_QUADS);
    for ( i = 0; i < FLOOR_TILES; ++i ) {
      z = -0.5 * size - 0.5 * TileSize;
      for ( j = 0; j < FLOOR_TILES; ++j ) {
        glVertex3f( x , height, z );
        glVertex3f( x , height, z + TileSize );
        glVertex3f( x + TileSize, height, z + TileSize );
        glVertex3f( x + TileSize, height, z );
        z += TileSize;
      }
      x += TileSize;
    }
  glEnd ();
} /*RenderSimpleFloor*/

void RenderLightPosition ( void )
{
  glDisable ( GL_LIGHTING );
  glColor3f ( 0.0,1.0,0.0 );
  glMatrixMode ( GL_MODELVIEW );
  glPushMatrix ();
  glRotatef ( 180.0*HorizAngle/PI,0.0,1.0,0.0 );
  glRotatef ( 180.0*VertAngle/PI,1.0,0.0,0.0 );
  glBegin(GL_LINES);
    glVertex3f(-0.1,0.0,R);
    glVertex3f(0.1,0.0,R);
    glVertex3f(0.0,-0.1,R);
    glVertex3f(0.0,0.1,R);
    glVertex3f(0.0,0.0,R+0.1);
    glVertex3f(0.0,0.0,R-0.1);
    glVertex3f(0.0,0.0,R);
    glVertex3f(0.0,0.0,R-0.6);
    glVertex3f(0.0,0.0,R-0.6);
    glVertex3f(0.3,0.0,R-0.3);
    glVertex3f(0.0,0.0,R-0.6);
    glVertex3f(-0.3,0.0,R-0.3);
  glEnd ();
  glPopMatrix ();
} /*RenderLightPosition*/

