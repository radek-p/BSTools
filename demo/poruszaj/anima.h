
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#define MAX_APARAMS     100
#define MAX_KEY_POSES   100
#define IC_DEGREE         3

#define TIME_MIN_STEP  0.05

extern double ic_knots[MAX_KEY_POSES+3];
extern boolean sw_periodic, spline_ok;

boolean SetNArtParams ( int nparams );
int GetNPoses ( void );
boolean EnterKeyPose ( double time, double *artp );
boolean DeleteKeyPose ( int k );
boolean MoveKeyKnot ( int k );
boolean ConstructInterpCurve ( void );
void ClampIt ( double *artp );
boolean GetPose ( double time, double *artp );

