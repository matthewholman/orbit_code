#include "nrutil.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "rebound.h"
#include "assist.h"

#define SUN 0
#define MERCURY 1
#define VENUS 2
#define EARTH 3
#define MOON 4
#define MARS 5
#define JUPITER 6
#define SATURN 7
#define URANUS 8
#define NEPTUNE 9
#define PLUTO 10
#define CAMILLA 11
#define CERES 12
#define CYBELE 13
#define DAVIDA 14
#define EUNOMIA 15
#define EUPHROSYNE 16
#define EUROPA 17
#define HYGIEA 18
#define INTERAMNIA 19
#define IRIS 20
#define JUNO 21
#define PALLAS 22
#define PSYCHE 23
#define SYLVIA 24
#define THISBE 25
#define VESTA 26

/*maximum number of observations*/
#define MAXOBS  10000    

/*Default name of binary DE432 ephemeris data:*/
#define DEFAULT_EPHEM_FILE "binEphem.432" 
/*Or override the above by looking for filename under
 * this environment variable:
 */
#define EPHEM_ENVIRON "ORBIT_EPHEMERIS"

/*uncertainty assumed for MPC-format observations (arcsec)*/
#define DEFAULT_DTHETA 0.2

#ifndef PI
#define PI 3.14159265358979323846
#endif
#define TPI     (2.*PI)
#define DTOR (PI/180.)
//#define GM 39.476926421373         /*solar gravitation*/
#define GM0      2.9591220828411956e-04    /*solar gravitation*/
#define SSMASS  1.00134       /*Total SS mass*/
#define ARCSEC (PI/180./3600.)
#define DAY (1./365.25) /*Julian day, 86400 s*/
//#define SPEED_OF_LIGHT  63241.06515 /*in AU/YR*/
#define SPEED_OF_LIGHT 173.1446         /*in AU/day */
#define ECL ((84381.448/3600)*PI/180.) /*Obliquity of ecliptic at J2000*/
//#define ECL ((84381.4118/3600)*PI/180.) /*Obliquity of ecliptic at J2000*/

//#define TSTEP (200.0*DAY) /*Time step for orbit integrator*
#define TSTEP 200.0     /*Time step for orbit integrator*/

/* set filename for ephemeris or observatory data */
//void
//set_ephem_file(char *fname);

/* an orbit as specified by its posn & velocity in fitting params */
typedef struct {
    double a;             /*alpha, angular position at t=0*/
    double adot;          /* (dx/dt) / z at t=0 */
    double b,bdot;        /* same 2 for y direction */
    double g,gdot;        /* 1/z and (dz/dt)/z of KBO at t=0 */
} PBASIS;

/* an orbit as specified by actual spatial position and velocity: */
typedef struct {
    double x, y, z;       /*3-space position */
    double xdot, ydot, zdot;       /*3-space velocity */
    double jd0;  /* Time at this phase space location*/
} XVBASIS;
 
/* an orbit as specified by the usual orbital parameters */
typedef struct {
    double a,e,i;                 /*semi-major axis, ellipticity, inclination*/
    double lan, aop, T;           /*long. of ascending node, arg of perihelion, and
                                  time from periapse of object */
    double meananom;
} ORBIT;
 
/* From skycalc, altered a bit: */
struct date_time {
       int y;
       int mo;
       float d;
       float h;
       float mn;
       float s;
};

/* Data of an individual observation */
/* Notice that in practice we will have to define an axis for the
 * coordinate system, origin of coordinate system, and time zeropoint.
 * These will for convenience be taken to correspond to first observation
 * in most cases.
 */
typedef struct {
    double thetax,dthetax; /*x-direction (ecliptic) position and uncertainty*/
    double thetay,dthetay; /*y-direction*/
    double x,dx; 
    double y,dy; 
    double z,dz;
    double ra, dec;
    double obstime;        /*time of observation (assumed years, TDB)*/
    char   *obscode;     /*Observatory site ID (from MPC list)*/
    double xe, ye, ze;     /*Position of observatory at obstime */
    char   *desig;
    int    reject;
} OBSERVATION;

typedef struct {
    double x, y, z, xd, yd, zd;
} State;

/* Declare the variables which define the coordinate system */
extern double ra0, dec0;       /* ecliptic lat & lon of tangent point */
extern double lat0, lon0;       /* ecliptic lat & lon of tangent point */
extern double xBary, yBary, zBary;  /*Posn of barycenter in our system*/
extern double jd0;           /* Zeropoint of time scale */

extern double machine_epsilon;

/* Some functions that we'll need */
 

/* 3-space position of Earth */
void earth3d(double t,
	     char *obscode,
	     double *x, double *y, double *z);

/*ICRS vector from SSBary to observatory location */
void earth_ssbary(double jd,
		  char *obscode,
		  double *x, double *y, double *z);

void body3d(struct assist_ephem* ephem,
	    double t,
	    int body,
	    double *GM,
	    double *x, double *y, double *z,
	    double *vxyz);

/*ICRS vector from SSBary to center of some body*/
void bodycenter_ssbary(double jd,
		       double *xyz,
		       int body,
		       double *vxyz);

/*ICRS vector from SSBary to geocenter */
void geocenter_ssbary(double jd,
		      double *xyz);


/* tangent-plane angular position of KBO and derivs */
extern void kbo2d(struct assist_ephem* ephem,
		  PBASIS *pin,
		  OBSERVATION *obs,
		  double *x, double dx[],
		  double *y, double dy[]);

extern void kbo3d(struct assist_ephem* ephem,
		  PBASIS *pin,
		  double t,
		  double *xout,
		  double *vout,
		  double dx[],
		  double dy[],
		  double dz[]);

extern void kbo3d_assist(struct assist_ephem* ephem,
			 PBASIS *pin,
			 double t,
			 double *xout,
			 double *vout,
			 double dx[],
			 double dy[],
			 double dz[]);


extern  void kbo_sphere(struct assist_ephem* ephem,
			PBASIS *pin, 
			OBSERVATION *obs,
			double *x, double dx[],
			double *y, double dy[],
			double *z, double dz[]);

/* linearized version of the 2d position, ignores gdot term*/
extern  void
kbo2d_linear(PBASIS *pin,
      OBSERVATION *obs,
      double *x, double dx[],
             double *y, double dy[]);
 
/* map from PBASIS to ORBIT format, and derivative matrix*/
extern  void    abg_to_aei(PBASIS *pin, ORBIT *orbout,
                           double **pderivs);
 
/* extract observations from a string */
int
scan_observation(char *inbuff,
   OBSERVATION *obs);
/* read RA/DEC from file & set up coord systems */
int
read_radec(OBSERVATION obsarray[], 
    char *fname, 
    int *nobs);
/* read RA/DEC from file & set up coord systems */
int
read_observations_simple(OBSERVATION obsarray[], 
    char *fname, 
    int *nobs);

/* reset astrometric error assigned to MPC obs.*/
void
set_mpc_dtheta(double d);


/* Read these options */
int
read_options(int* iarg, 
      int argc,
      char *argv[]);

/* preliminary fit to observations */
void
prelim_fit(OBSERVATION obsarray[],
    int nobs,
    PBASIS *pout,
    double **covar);

/* Routine to predict position and uncertainty at any time, given
 * a PBASIS fit and a sigma matrix.
 */
void
predict_posn(struct assist_ephem* ephem,
	     PBASIS *pin,
             double **covar,
             OBSERVATION *obs,
             double **sigxy);

/* Routine to predict position and uncertainty at any time, given
 * a PBASIS fit and a sigma matrix.
 */
void predict_posn_linear(PBASIS *pin,
    double **covar,
    OBSERVATION *obs,
    double **sigxy);

int fit_observations(struct reb_simulation* sim,
		     struct assist_ephem* ephem,
		     OBSERVATION obsarray[],
		     int nobs,
		     PBASIS *p,
		     double **covar,
		     double *chisq,
		     int *dof,
		     FILE *logfile);

int fit_observations_prelim(struct reb_simulation* sim,
		     struct assist_ephem* ephem,
		     OBSERVATION obsarray[],
		     int nobs,
		     PBASIS *p,
		     double **covar,
		     double *chisq,
		     int *dof,
		     FILE *logfile);

int fit_prelim(OBSERVATION obsarray[],
       int nobs,
       PBASIS *p,
       double **covar,
       double *chisq,
       int *dof,
       FILE *logfile);

/* Take a PBASIS orbit and uncertainty and convert it to a traditional
 * ORBIT basis and uncertainty.
 */
void
predict_aei(PBASIS *pin, double **covarp,
     ORBIT  *orbout, double **covarq);
void
print_matrix(FILE *fptr, double **matrix, int xdim, int ydim);

char *
fgets_nocomment(char *buffer, int length, FILE *fptr, FILE *fpout);
/* Coordinate transformation routines: */
void
eq_to_ec( double ra_eq,   double dec_eq,
   double *lat_ec, double *lon_ec,
   double **partials);
void
xyz_eq_to_ec(double x_eq, double y_eq, double z_eq,
      double *x_ec, double *y_ec, double *z_ec,
      double **partials);
void
ec_to_eq( double lat_ec,  double lon_ec,
   double *ra_eq,  double *dec_eq,
   double **partials);
void
xyz_ec_to_eq(double x_ec, double y_ec, double z_ec,
      double *x_eq, double *y_eq, double *z_eq,
      double **partials);
void
ec_to_proj(double lat_ec,   double lon_ec,
    double *x_proj,  double *y_proj,
    double lat0,     double lon0,
    double **partials);

void radec_to_proj(double ra,
		   double dec,
		   double *x_proj,
		   double *y_proj,
		   double ra0,
		   double dec0,
		   double **partials);

void proj_to_radec(double x_proj,
		   double y_proj,
		   double *ra,
		   double *dec,
		   double ra0,
		   double dec0,
		   double **partials);

void
proj_to_ec(double x_proj,   double y_proj,
    double *lat_ec,  double *lon_ec,
    double lat0,     double lon0,
    double **partials);
void
xyz_ec_to_proj( double x_ec, double y_ec, double z_ec,
  double *x_p, double *y_p, double *z_p,
  double lat0, double lon0,
  double **partials);
void
xyz_proj_to_ec( double x_p, double y_p, double z_p,
  double *x_ec, double *y_ec, double *z_ec,
  double lat0, double lon0,
  double **partials);

void xyz_eq_to_proj( double x_eq, double y_eq, double z_eq,
		     double *x_p, double *y_p, double *z_p,
		     double ra0, double dec0,
		     double **partials);

void
pbasis_to_bary(PBASIS *p,
        XVBASIS *xv,
        double **partials);
void
pbasis_to_bary_eq(PBASIS *p,
		  XVBASIS *xv,
		  double **partials);

void
matrix_multiply(double **m1, double **m2,
  double **mout,
  int x1, int y1, int x2, int y2);
void
orbitElements(XVBASIS *xv,
       ORBIT  *orb);
void
elements_to_xv(ORBIT *o,
        double jd,
        XVBASIS *xv);
void
elements_to_pbasis(ORBIT *o,
     double jd,
     char *obscode,
     PBASIS *p);
void
heliocentric_elements_to_pbasis(ORBIT *o,
     double jd,
     char *obscode,
     PBASIS *p);
void
covar_map(double **covar_in, 
   double **derivs, 
   double **covar_out,
   int kin, int kout);
double 
date_to_jd(struct date_time date);
void
aei_derivs( XVBASIS *xv,
     double **daei_dxv);

void
fake_observation(PBASIS *p, 
   OBSERVATION *obs);

void
fiddle_observation(OBSERVATION *obs, OBSERVATION *obs_noise);

void keplerian(double gm, State state, double *a, double *e, double *i, double *longnode, double *argperi, double *meananom);

void cartesian(double gm, 
        double a, double e, double i, double longnode, double argperi, double meananom, 
        State *state);

double principal_value(double theta);

double principal_value_0(double theta);

typedef struct {
  double x, y, z;
} Vector;

void intro(int argc, char *argv[]);

