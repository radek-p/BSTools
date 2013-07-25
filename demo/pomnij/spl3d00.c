
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkgeom.h"
#include "multibs.h"
#include "convh.h"
#include "camerad.h"

#include "xgedit.h"
#include "spl3d.h"

int win1_contents = WIN1_GENERAL;

int     degree_u, lastknot_u, degree_v, lastknot_v;
double  knots_u[(MAX_DEGREE+1)*MAX_KNOTS];
double  knots_v[(MAX_DEGREE+1)*MAX_KNOTS];
point4d cpoints[MAX_CPOINTS];            /* surface control points */
point4d savedcpoints[MAX_CPOINTS];       /* a copy */
point2d rpoints[4][MAX_CPOINTS];         /* projections */
boolean clpoints[4][MAX_CPOINTS];        /*true if fits into the frame */
byte    mkpoints[MAX_CPOINTS];           /* point marking bits */

double clcTu, clcTv;  /* closed patch domain width and height */
int   clcKu, clcKv;  /* closed patch indexing periods */

boolean display_surface = true;
boolean display_control_net = true;
boolean display_Bezier_nets = false;
boolean display_constr_poly = false;
boolean move_many_knots = false;
boolean closed_u = false;
boolean closed_v = false;
boolean display_domain_net = true;
boolean domain_selecting_mode = false;

int display_bez_dens_u = 6;
int display_bez_dens_v = 6;

int current_point;
int current_knot, current_mult;

/* the variables below are used in the construction of */
/* spherical product surfaces */
boolean equator = true;
boolean meridian = false;
boolean bind_spr = false;
boolean eqmer_nurbs = false;
boolean eqmer_closed = false;
double   arc_angle = 0.0;

int     meridian_deg, meridian_lastknot;
double   meridian_knots[(MAX_DEGREE+1)*MAX_KNOTS];
point3d meridian_cpoints[2*(MAX_DEGREE+1)*(MAX_KNOTS-1)];
point2d meridian_rpoints[2*(MAX_DEGREE+1)*(MAX_KNOTS-1)];
byte    meridian_mkpoints[2*(MAX_DEGREE+1)*(MAX_KNOTS-1)];
boolean meridian_nurbs = false;
boolean meridian_closed = false;

int     equator_deg, equator_lastknot;
double   equator_knots[(MAX_DEGREE+1)*MAX_KNOTS];
point3d equator_cpoints[2*(MAX_DEGREE+1)*(MAX_KNOTS-1)];
point2d equator_rpoints[2*(MAX_DEGREE+1)*(MAX_KNOTS-1)];
byte    equator_mkpoints[2*(MAX_DEGREE+1)*(MAX_KNOTS-1)];
boolean equator_nurbs = false;
boolean equator_closed = false;

int     neqmerpoints;
point3d *eqmer_cpoints  = NULL;
point2d *eqmer_rpoints  = NULL; 
byte    *eqmer_mkpoints = NULL;

boolean display_eqmer_curve = true;
boolean display_eqmer_control_polygon = true;
boolean display_eqmer_Bezier_polygons = false;
boolean display_eqmer_ticks = false;

/* the variables below are used in the construction */
/* of biharmonic or triharmonic blending surfaces */
boolean sw_blending_g1 = false;
boolean sw_blending_g2 = true;
boolean sw_bind_blending = false;
boolean blending_mat_valid = false;
boolean sw_triharmonic_blending = false;
int     blending_n, blending_lknu, blending_lknv;
double  *blending_Amat = NULL;
double  **blending_Arow = NULL;
int     *blending_prof = NULL;

/* the variables below are used in the construction of */
/* blending surfaces obtained via nonlinear shape optimization */
boolean sw_clamped_blending = false;
boolean sw_nonlin_blending = false;
boolean sw_blending_constraints = false;
  /* this is the initial slidebar position */
double  blending_factor = 0.66666666666666666666;
int     blending_lmt_iter = 20;

boolean sw_blending_opt_entire = true;
int     blending_opt_part[4] = {3,3,3,3};
trans3d blending_opt_transform;
boolean display_trans_net = false;

boolean sw_blending_opt_dump = false;

/* blending constraints data */
int max_blending_constraints = MAX_BLENDING_CONSTRAINTS;
int n_blending_constraints = 0;
double blending_constr_knots[MAX_BLENDING_CONSTRAINTS+1];
boolean blending_constr_poly_valid[MAX_BLENDING_CONSTRAINTS+1];
point4d blending_constr_cp[MAX_BLENDING_CONSTRAINTS*(MAX_KNOTS)];
point2d blending_constr_rp[4][MAX_BLENDING_CONSTRAINTS*(MAX_KNOTS)];
boolean clblending_constr[4][MAX_BLENDING_CONSTRAINTS*(MAX_KNOTS)];
byte    mkblending_cp[MAX_BLENDING_CONSTRAINTS*(MAX_KNOTS)];
point4d savedbl_constr_cp[MAX_BLENDING_CONSTRAINTS*(MAX_KNOTS)];

