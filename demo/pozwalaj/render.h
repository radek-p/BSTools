
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#define MAXPLEVEL      15
#define MAXCLEVEL      24
#define MAXINTERS      32

#define R_NLIGHTS       4

/* scaling factor range for noncontinuous palette */
#define R_MINDFSF        0.1
#define R_MAXDFSF     1000.0

extern boolean RenderingIsOn;
extern XImage  *rendimage;

/* switches selecting the shape function for visualisation */
extern boolean swGaussian_c, swMean_c, swLambert_c, swReflection_c,
               swHighlight_c, swSections_c,
               swGaussian_d, swMean_d, swLambert_d, swReflection_d,
               swHighlight_d, swSections_d;
extern boolean swShadows, swAntialias;

extern double render_dfsf;
extern double render_cfrange[2];

/* ///////////////////////////////////////////////////////////////////////// */
/* to be called at the beginning and end of the application */
boolean RendInit ( void );
void RendDestroy ( void );

/* to be called before rendering a new scene */
boolean RendReset ( void );

boolean RendEnterTriangle3d ( point3d *p0, point3d *p1, point3d *p2,
                              double *colour );
boolean RendEnterBezPatch3d ( int n, int m, const point3d *cp, double *colour );
boolean RendEnterBSPatch3d ( int n, int lknu, const double *knu,
                             int m, int lknv, const double *knv,
                             const point3d *cp, double *colour );
boolean RendEnterBezPatch3Rd ( int n, int m, const point4d *cp, double *colour );
boolean RendEnterBSPatch3Rd ( int n, int lknu, const double *knu,
                              int m, int lknv, const double *knv,
                              const point4d *cp, double *colour );
boolean RendEnterBezCurve3d ( int n, const point3d *cp, double r, double *colour );
boolean RendEnterBSCurve3d ( int n, int lkn, const double *kn,
                             const point3d *cp, double r, double *colour );
boolean RendEnterBezCurve3Rd ( int n, const point4d *cp, double r, double *colour );
boolean RendEnterBSCurve3Rd ( int n, int lkn, const double *kn,
                              const point4d *cp, double r, double *colour );

boolean RendEnterCamerad ( CameraRecd *CPos, xge_widget *er );
void RendEnterLightsd ( int nlights, const vector3d *light_dir,
                        const double *light_int );
void RendEnterReflectionLinesFramed ( point3d rf[3] );
void RendEnterHighlightLinesFramed ( point3d hf[3] );
void RendEnterSectionPlanesNormald ( vector3d *spn );

/* actual rendering */
boolean RendBegin ( void );
boolean RendRestart ( void );
void InitRenderingAA ( void );
int RenderLineA ( void );
int RenderLineAA ( void );
int RenderLine ( void );

