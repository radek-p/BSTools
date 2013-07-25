
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <unistd.h>
#include <sys/times.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camerad.h"
#include "multibs.h"
#include "xgedit.h"

#include "render.h"
#include "ed3dswidgets.h"
#include "spl3d.h"
#include "ed3ds.h"


void ResetBlendingOptTrans ( void )
{
  char zero[] = "0.0";
  char one[] = "1.0";
  int  i, j, k;

  IdentTrans3d ( &blending_opt_transform );
  for ( i = k = 0; i < 2; i++ ) {
    memcpy ( blending_trans_str[k++], one, 4 );
      for ( j = 0;  j < 3;  j++ )
        memcpy ( blending_trans_str[k++], zero, 4 );
  }
  memcpy ( blending_trans_str[k], one, 4 );
  for ( i = 0; i < 9; i++ )
    blending_trans_editor[i].start = blending_trans_editor[i].pos = 0;
  xge_SetClipping ( popup11 );
  popup11->redraw ( popup11, true );
} /*ResetBlendingOptTrans*/

boolean EnterBlendingOptTransCoefficient ( int id )
{
  double x;
  int    res;
  int i, j;

  id -= txtedP11A0;
  res = sscanf ( blending_trans_str[id], "%lf", &x );
  if ( res == 1 ) {
    i = id / 3;
    j = id % 3;
    blending_opt_transform.U1.a[i][j] = x;
    return true;
  }
  else
    return false;
} /*EnterBlendingOptTransCoefficient*/

