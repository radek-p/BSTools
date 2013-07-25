
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#define R_MAXPATCHES 1024
#define MAXLEVEL       15
#define MAXINTERS      10

#define R_NLIGHTS       4

extern boolean RenderingIsOn;
extern XImage  *rendimage;

/* switches selecting the shape function for visualisation */
extern boolean swGaussian_c, swMean_c, swLambert_c, swReflection_c,
               swHighlight_c, swSections_c,
               swGaussian_d, swMean_d, swLambert_d, swReflection_d,
               swHighlight_d, swSections_d;

/* ///////////////////////////////////////////////////////////////////////// */
/* to be called at the beginning and end of the application */
boolean RendInit ( void );
void RendDestroy ( void );

/* to be called before rendering a new scene */
boolean RendReset ( void );
boolean RendEnterBezPatchd ( int n, int m, const point3d *cp );
boolean RendEnterBSPatchd ( int n, int lknu, const double *knu,
                            int m, int lknv, const double *knv,
                            const point3d *cp );
boolean RendEnterCamerad ( CameraRecd *CPos, xge_widget *er );
void RendEnterLightsd ( int nlights, const vector3d *light_dir,
                        const double *light_int );

/* actual rendering */
boolean RendBegin ( boolean swShadows, boolean swAntialias );
void InitRenderingAA ( void );
int RenderLineA ( void );
int RenderLineAA ( void );
int RenderLine ( void );

