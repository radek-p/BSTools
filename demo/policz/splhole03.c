
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#include <stdio.h>    
#include <stdlib.h>   
#include <math.h>
#include <malloc.h>  
#include <string.h>     

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "pknum.h"  
#include "pkgeom.h"   
#include "camerad.h"
#include "multibs.h"
#include "eg1holed.h"
#include "eg2holed.h"
#include "xgedit.h"   

#include "render.h"
#include "edpolicz.h"
#include "edwidgets.h"
#include "drawbezd.h"
#include "splhole.h"
#include "datagend.h"

/* ////////////////////////////////////////////////////////////////////////// */
void SendPatchesToRenderer ( void )
{
  void    *sp;
  int     i, j;
  point3d *bcp;

  sp = pkv_GetScratchMemTop ();
  if ( (bcp = pkv_GetScratchMem ( 16*sizeof(point3d) ) ) ) {
        /* send the patches surrounding the surface */
    if ( view_surf_spatches ) {
      for ( i = 0; i < hole_k; i++ )
        for ( j = 0; j < 3; j++ )
          if ( GetSurfBezPatch ( i, j, bcp ) )
            RendEnterBezPatchd ( 3, 3, bcp );
    }
        /* send the patches filling the hole */
    if ( view_surf_1 && final_cp1 ) {
      if ( options1.spline ) {
        j = final_lkn1-final_deg1;  j *= j;
        for ( i = 0; i < hole_k; i++ )
          RendEnterBSPatchd ( final_deg1, final_lkn1, final_knots1,
                              final_deg1, final_lkn1, final_knots1,
                              &final_cp1[i*j] );
      }
      else {
        j = (final_deg1+1)*(final_deg1+1);
        for ( i = 0; i < hole_k; i++ )
          RendEnterBezPatchd ( final_deg1, final_deg1, &final_cp1[i*j] );
      }
    }
    if ( view_surf_2 && final_cp2 ) {
      if ( options2.spline ) {
        j = final_lkn2-final_deg2;  j *= j;
        for ( i = 0; i < hole_k; i++ )
          RendEnterBSPatchd ( final_deg2, final_lkn2, final_knots2,
                              final_deg2, final_lkn2, final_knots2,
                              &final_cp2[i*j] );
      }
      else {
        j = (final_deg2+1)*(final_deg2+1);
        for ( i = 0; i < hole_k; i++ )
          RendEnterBezPatchd ( final_deg2, final_deg2, &final_cp2[i*j] );
      }
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*SendPatchesToRenderer*/

