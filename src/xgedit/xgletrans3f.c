
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2011                                  */
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

#include "xgledit.h"
#include "xgeprivate.h"

void xgle_MultMatrix3f ( trans3f *tr )
{
  GLfloat a[16];

  a[0]  = tr->U0.a11;  a[1]  = tr->U0.a21;  a[2]  = tr->U0.a31;  a[3]  = 0.0;
  a[4]  = tr->U0.a12;  a[5]  = tr->U0.a22;  a[6]  = tr->U0.a32;  a[7]  = 0.0;
  a[8]  = tr->U0.a13;  a[9]  = tr->U0.a23;  a[10] = tr->U0.a33;  a[11] = 0.0;
  a[12] = tr->U0.a14;  a[13] = tr->U0.a24;  a[14] = tr->U0.a34;  a[15] = 1.0;
  glMultMatrixf ( a );
} /*xgle_MultMatrix3f*/

