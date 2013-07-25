
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <GL/gl.h>
#include <GL/glu.h>

#include "pkgeom.h"
#include "texture.h"


#define texImageWidth 128
#define texImageHeight 128
#define texImageW1 64
#define texImageW2 32
#define texImageW3 16
#define texImageW4 8
#define texImageW5 4
#define texImageW6 2
#define texImageW7 1

static GLubyte faceImage[texImageHeight][texImageWidth][4];
static GLubyte faceImage1[texImageW1][texImageW1][4];
static GLubyte faceImage2[texImageW2][texImageW2][4];
static GLubyte faceImage3[texImageW3][texImageW3][4];
static GLubyte faceImage4[texImageW4][texImageW4][4];
static GLubyte faceImage5[texImageW5][texImageW5][4];
static GLubyte faceImage6[texImageW6][texImageW6][4];
static GLubyte faceImage7[texImageW7][texImageW7][4];
GLuint texName;

void FaceTexture ( double x, double y, GLubyte c[4] )
{
  double xx, yy, d;
  GLubyte cv;

  if ( x >= 0.5 )
    xx = x-0.5;  
  else
    xx = x;
  if ( y >= 0.5 )
    yy = y-0.5;  
  else
    yy = y;

  xx -= 0.25;

  if ( yy >= 0.375 && yy <= 0.4 && fabs(xx) <= 0.17 )
    cv = 0;
  else {   
    xx = fabs(xx)-0.08;
    yy -= 0.18;
    d = xx*xx+yy*yy;
    if ( d <= 0.01 && d > 0.002 )
      cv = 0;
    else
      cv = 255;
  }
   
  if ( x >= 0.5 )
    cv = 255-cv; 
  if ( y >= 0.5 )
    cv = 255-cv; 

  c[0] = c[1] = c[2] = cv;
  c[3] = 255;
} /*FaceTexture*/

void makeTexImages ( void )
{
  int i, j, k;
  
  for ( j = 0; j < texImageHeight; j++ )
    for ( i = 0; i < texImageWidth; i++ )
      FaceTexture ( ((double)i)/texImageWidth, 1.0-((double)j)/texImageHeight,
                    faceImage[j][i] );

                 /* generate textures for MIP mapping */
  for ( j = 0; j < texImageW1; j++ )
    for ( i = 0; i < texImageW1; i++ )
    {
      for ( k = 0; k < 4; k++ )
        faceImage1[i][j][k] = ( (int)faceImage[i+i][j+j][k] +
                                (int)faceImage[i+i+1][j+j][k] +
                                (int)faceImage[i+i][j+j+1][k] +
                                (int)faceImage[i+i+1][j+j+1][k] ) / 4;

      faceImage1[i][j][1] = faceImage1[i][j][2] = 0;
    } 
  for ( j = 0; j < texImageW2; j++ )
    for ( i = 0; i < texImageW2; i++ )
      for ( k = 0; k < 4; k++ )
        faceImage2[i][j][k] = ( (int)faceImage1[i+i][j+j][k] +
                                (int)faceImage1[i+i+1][j+j][k] +
                                (int)faceImage1[i+i][j+j+1][k] +
                                (int)faceImage1[i+i+1][j+j+1][k] ) / 4;
  for ( j = 0; j < texImageW3; j++ )
    for ( i = 0; i < texImageW3; i++ )
      for ( k = 0; k < 4; k++ )
        faceImage3[i][j][k] = ( (int)faceImage2[i+i][j+j][k] +
                                (int)faceImage2[i+i+1][j+j][k] +
                                (int)faceImage2[i+i][j+j+1][k] +
                                (int)faceImage2[i+i+1][j+j+1][k] ) / 4;
  for ( j = 0; j < texImageW4; j++ )
    for ( i = 0; i < texImageW4; i++ )
      for ( k = 0; k < 4; k++ )
        faceImage4[i][j][k] = ( (int)faceImage3[i+i][j+j][k] +
                                (int)faceImage3[i+i+1][j+j][k] +
                                (int)faceImage3[i+i][j+j+1][k] +
                                (int)faceImage3[i+i+1][j+j+1][k] ) / 4;

} /*makeTexImages*/

void InitTexture ( void )
{
  makeTexImages ();
  glPixelStorei ( GL_UNPACK_ALIGNMENT, 1 );

  glEnable(GL_TEXTURE_2D);
  glGenTextures ( 1, &texName );
  glBindTexture ( GL_TEXTURE_2D, texName );
  glTexParameteri ( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
  glTexParameteri ( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
  glTexParameteri ( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
  glTexParameteri ( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
  glTexParameteri ( GL_TEXTURE_2D, GL_TEXTURE_ENV_MODE, GL_MODULATE);
  glTexImage2D ( GL_TEXTURE_2D, 0, GL_RGBA, texImageWidth, texImageHeight,
                 0, GL_RGBA, GL_UNSIGNED_BYTE, faceImage );
  glTexImage2D ( GL_TEXTURE_2D, 1, GL_RGBA, texImageW1, texImageW1,
                 0, GL_RGBA, GL_UNSIGNED_BYTE, faceImage1 );
  glTexImage2D ( GL_TEXTURE_2D, 2, GL_RGBA, texImageW2, texImageW2,
                 0, GL_RGBA, GL_UNSIGNED_BYTE, faceImage2 );
  glTexImage2D ( GL_TEXTURE_2D, 3, GL_RGBA, texImageW3, texImageW3,
                 0, GL_RGBA, GL_UNSIGNED_BYTE, faceImage3 );
  glTexImage2D ( GL_TEXTURE_2D, 4, GL_RGBA, texImageW4, texImageW4,
                 0, GL_RGBA, GL_UNSIGNED_BYTE, faceImage4 );
  glTexImage2D ( GL_TEXTURE_2D, 5, GL_RGBA, texImageW5, texImageW5,
                 0, GL_RGBA, GL_UNSIGNED_BYTE, faceImage5 );
  glTexImage2D ( GL_TEXTURE_2D, 6, GL_RGBA, texImageW6, texImageW6,
                 0, GL_RGBA, GL_UNSIGNED_BYTE, faceImage6 );
  glTexImage2D ( GL_TEXTURE_2D, 7, GL_RGBA, texImageW7, texImageW7,
                 0, GL_RGBA, GL_UNSIGNED_BYTE, faceImage7 );
  glBindTexture ( GL_TEXTURE_2D, 0 );
} /*InitTexture*/

