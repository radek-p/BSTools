
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */
/* this file was written by Tomasz Swierczek            */

#ifndef LIGHTING_H
#define LIGHTING_H 

extern boolean useSpecular;
extern float HorizAngle, VertAngle;

void InitSpotlight ( void );
void MoveLight ( float dr );
void ChangeLightColor ( float dr, float dg, float db );
void SetLightAttenuation ( float con, float lin, float sq );
void ChangeSpotCutOff ( float dangle );
void SetupLightToDisplay ( void );
void DisableLightToDisplay ( void );

void RenderSimpleFloor ( float height, float size );
void RenderLightPosition ( void );

#endif
