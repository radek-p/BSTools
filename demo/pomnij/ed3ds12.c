
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


void SetKWindNKnots ( void )
{
  if ( kwind.switchknu ) {
    kwind.altdeg_u = degree_u;
    kwind.altlastknot_u = lastknot_u;
  }
  else {
    kwind.degree_u = degree_u;
    kwind.lastknot_u = lastknot_u;
  }
  if ( kwind.switchknv ) {
    kwind.altdeg_v = degree_v;
    kwind.altlastknot_v = lastknot_v;
  }
  else {
    kwind.degree_v = degree_v;
    kwind.lastknot_v = lastknot_v;
  }
} /*SetKWindNKnots*/

