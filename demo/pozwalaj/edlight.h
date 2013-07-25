
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#define RENDERING_NLIGHTS       4

    /* renderer lights */
extern boolean  sw_edit_light[RENDERING_NLIGHTS];
extern vector3d render_light_dir[RENDERING_NLIGHTS];
extern double   render_light_int[RENDERING_NLIGHTS+1];

boolean EdLightFindNearestPoint ( CameraRecd *CPos, short x, short y );
void EdLightSetPoint ( CameraRecd *CPos, short x, short y );
void EdLightSavePoints ( void );
void EdLightTransformPoints ( trans3d *tr );
void EdLightDisplay ( xge_3Dwind *ww );

    /* shape function definition elements */
extern boolean  sw_reflection_frame, sw_highlight_frame, sw_sections_frame;
extern vector3d edshapef_sectiondir;
extern point3d  edshapef_reflection_frame[3], edshapef_highlight_frame[3];

boolean EdShapeFuncFindNearestPoint ( CameraRecd *CPos, short x, short y );
void EdShapeFuncSetPoint ( CameraRecd *CPos, short x, short y );
void EdShapeFuncSavePoints ( void );
void EdShapeFuncTransformPoints ( trans3d *tr );
void EdShapeFuncDisplay ( xge_3Dwind *ww );

    /* draw it all */
void DisplaySpecial3DElements ( xge_3Dwind *ww );

