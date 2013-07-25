
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#define R_MAXPATCHES 8192
#define MAXLEVEL       15
#define MAXINTERS      10

#define R_NLIGHTS       4

/* scaling factor range for noncontinuous palette */
#define R_MINDFSF        0.1
#define R_MAXDFSF     1000.0

extern boolean RenderingIsOn;
extern XImage  *rendimage;

/* switches selecting the shape function for visualisation */
extern boolean swGaussian_c, swMean_c, swLambert_c, swReflection_c,
               swHighlight_c, swSections_c, swParam_c,
               swGaussian_d, swMean_d, swLambert_d, swReflection_d,
               swHighlight_d, swSections_d, swParam_d;

extern float render_dfsf;
extern float render_cfrange[2];

/* ///////////////////////////////////////////////////////////////////////// */
/* to be called at the beginning and end of the application */
boolean RendInit ( void );
void RendDestroy ( void );

/* to be called before rendering a new scene */
boolean RendReset ( void );
boolean RendEnterBezPatch3Rf ( int n, int m, const point4f *cp );
boolean RendEnterBSPatch3Rf ( int n, int lknu, const float *knu,
                              int m, int lknv, const float *knv,
                              const point4f *cp );
boolean RendEnterCameraf ( CameraRecf *CPos, xge_widget *er );
void RendEnterLightsf ( int nlights, const vector3f *light_dir,
                        const float *light_int );

boolean RendEnterBezPatch3Rd ( int n, int m, const point4d *cp );
boolean RendEnterBSPatch3Rd ( int n, int lknu, const double *knu,
                              int m, int lknv, const double *knv,
                              const point4d *cp );
boolean RendEnterCamerad ( CameraRecd *CPos, xge_widget *er );
void RendEnterLightsd ( int nlights, const vector3d *light_dir,
                        const double *light_int );

/* actual rendering */
boolean RendBegin ( boolean swShadows, boolean swAntialias );
void InitRenderingAA ( void );
int RenderLineA ( void );
int RenderLineAA ( void );
int RenderLine ( void );

