
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
void EditSetSelectView ( void )
{
  swind_ed_switch = SWIN_SELECT_VIEW;
} /*EditSetSelectView*/

void EditSetSurface ( void )
{
  swind_ed_switch = SWIN_EDITING_SURFACE;
} /*EditSetSurface*/

void EditSetConstraints ( void )
{
  swind_ed_switch = SWIN_EDITING_CONSTRAINTS;
} /*EditSetConstraints*/

void EditSetLight ( void )
{
  swind_ed_switch = SWIN_EDITING_LIGHT;
} /*EditSetLight*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean FindNearestSWinPoint ( int id, short x, short y, short mindist )
{
  id &= 0x03;
  switch ( swind_ed_switch ) {
case SWIN_EDITING_SURFACE:
    return FindNearestSurfCPoint ( id, x, y, mindist );
case SWIN_EDITING_CONSTRAINTS:
    return FindNearestConstraintCPoint ( id, x, y, mindist );
case SWIN_EDITING_LIGHT:
    return FindNearestLightPoint ( id, x, y, mindist );
default:
    return false;
  }
} /*FindNearestSWinPoint*/

void SetSWinPoint ( int id, short x, short y )
{
  id &= 0x03;
  switch ( swind_ed_switch ) {
case SWIN_EDITING_SURFACE:
    SetSurfCPoint ( id, x, y );
    break;
case SWIN_EDITING_CONSTRAINTS:
    SetConstraintCPoint ( id, x, y );
    break;
case SWIN_EDITING_LIGHT:
    SetLightPoint ( id, x, y );
    break;
default:
    break;
  }
} /*SetSWinPoint*/

void SelectSWinPoints ( int id )
{
  id &= 0x03;
  switch ( swind_ed_switch ) {
case SWIN_EDITING_SURFACE:
    SelectSurfCPoints ( id );
    break;
case SWIN_EDITING_CONSTRAINTS:
    SelectConstraintCPoints ( id );
    break;
case SWIN_EDITING_LIGHT:
    break;
default:
    break;
  }
} /*SelectSWinPoints*/

void UnselectSWinPoints ( int id )
{
  id &= 0x03;
  switch ( swind_ed_switch ) {
case SWIN_EDITING_SURFACE:
    UnselectSurfCPoints ( id );
    break;
case SWIN_EDITING_CONSTRAINTS:
    UnselectConstraintCPoints ( id );
    break;
case SWIN_EDITING_LIGHT:
    break;
default:
    break;
  }
} /*UnselectSWinPoints*/

void SaveSWinPoints ( void )
{
  switch ( swind_ed_switch ) {
case SWIN_EDITING_SURFACE:
    SaveSurfCPoints ();
    break;
case SWIN_EDITING_CONSTRAINTS:
    SaveConstraintCPoints ();
    break;
case SWIN_EDITING_LIGHT:
    break;
default:
    break;
  }
} /*SaveSWinPoints*/

void TransformSWinPoints ( void )
{
  switch ( swind_ed_switch ) {
case SWIN_EDITING_SURFACE:
    TransformSurfCPoints ();
    break;
case SWIN_EDITING_CONSTRAINTS:
    TransformConstraintCPoints ();
    break;
case SWIN_EDITING_LIGHT:
    break;
default:
    break;
  }
} /*TransformSWinPoints*/

