
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
  

xge_3Dwind swind;
xge_2Dwind domwind;
ghKnotWind knwind;

boolean swind_picture = false;  /* if nonzero, a ray-traced picture is on screen */

char swind_ed_switch = SWIN_EDITING_SURFACE;

GHoptions options1 = { 2, false, true, false, false,
                       true, false, false, false, 1, 1, 1, 0.5, 10.0, 0 };
GHoptions options2 = { 1, false, true, false, false,
                       true, false, false, false, 1, 1, 1, 0.5, 10.0, 0 };

GHoleDomaind *domain1 = NULL, *domain2 = NULL;
double *acoeff1 = NULL, *acoeff2 = NULL;

int     hole_k =  5;
int     nctrlp = 61;
double  knots[GH_MAX_K*11] =
  {0.0,0.0,0.125,0.250,0.375,0.5,0.625,0.75,0.875,1.0,1.0,
   0.0,0.0,0.125,0.250,0.375,0.5,0.625,0.75,0.875,1.0,1.0,
   0.0,0.0,0.125,0.250,0.375,0.5,0.625,0.75,0.875,1.0,1.0,
   0.0,0.0,0.125,0.250,0.375,0.5,0.625,0.75,0.875,1.0,1.0,
   0.0,0.0,0.125,0.250,0.375,0.5,0.625,0.75,0.875,1.0,1.0};
point2d domain_cp[MAX_CPOINTS];
point3d hole_cp[MAX_CPOINTS];
unsigned char mkdcp[MAX_CPOINTS];
unsigned char mkhcp[MAX_CPOINTS];

char    constr_surfno = 1;
point3d saved_cp[MAX_SAVED_CPOINTS];
point2d rpoints[4][MAX_CPOINTS];
point2d rdpoints[MAX_CPOINTS];

int     final_np1, final_deg1, final_lkn1;
point3d *final_cp1 = NULL;
double  *final_knots1 = NULL;
int     final_np2, final_deg2, final_lkn2; 
point3d *final_cp2 = NULL;
double  *final_knots2 = NULL;

/* domain patches - always Bezier */
int domain_np1, domain_deg1;
point2d *domain_bcp1 = NULL;
int domain_np2, domain_deg2;
point2d *domain_bcp2 = NULL;

/* parameters for the built-in constructions of the domain and surface */
vector2d domcvect[GH_MAX_K];
double   domcparam[DATAGEN_DOM_PARAMS] = {0.0,0.0,0.0,0.0};
vector3d surfcvect[GH_MAX_K];
double   surfcparam[DATAGEN_SURF_PARAMS] = {0.0,0.0,0.0,0.0,0.0};


