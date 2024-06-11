/*********** Coordinate transformation routines *************/
/* All angles assumed to be in radians upon input */

#include "orbfit.h"
//#include "ephem_types.h"

/* First takes RA,DEC in equatorial to ecliptic */
void eq_to_ec( double ra_eq,
	       double dec_eq,
	       double *lat_ec,
	       double *lon_ec,
	       double **partials)
{
    double sd, cd, se, ce, y, x;

    se = sin(ECL);
    ce = cos(ECL);

    sd = ce * sin(dec_eq) - se * cos(dec_eq)*sin(ra_eq);
    *lat_ec = asin(sd);

    y = ce*cos(dec_eq)*sin(ra_eq) + se*sin(dec_eq);
    x = cos(dec_eq)*cos(ra_eq);
    *lon_ec = atan2(y,x);

    if (partials!=NULL) {
	cd = sqrt(1.-sd*sd);
	/* form Jacobian matrix  (assume not at pole)*/
	partials[1][1] = -se*cos(dec_eq)*cos(ra_eq)/cd;
	partials[1][2] = (ce*cos(dec_eq)+se*sin(dec_eq)*sin(ra_eq)) / cd;
	partials[2][1] = cos(dec_eq)*partials[1][2] / cd;
	partials[2][2] = se*cos(ra_eq) / (cd*cd);
    }

    return;
}

/* And transform x,y,z from eq to ecliptic */
void xyz_eq_to_ec(double x_eq, double y_eq, double z_eq,
		  double *x_ec, double *y_ec, double *z_ec,
		  double **partials)
{
    double se, ce;

    se = sin(ECL);
    ce = cos(ECL);

    *x_ec = x_eq;
    *y_ec = ce*y_eq + se*z_eq;
    *z_ec = -se*y_eq + ce*z_eq;

    if (partials!=NULL) {
	partials[1][1] = 1.;
	partials[1][2] = partials[1][3] =
	    partials[2][1] = partials[3][1] = 0.;
	partials[2][2] = partials[3][3] = ce;
	partials[2][3] = se;
	partials[3][2] = -se;
    }

    return;
}

/**To reverse above, just flip sign of ECL effectively: ***/
void ec_to_eq( double lat_ec,
	       double lon_ec,
	       double *ra_eq,
	       double *dec_eq,
	       double **partials)
{
    double sd, cd, se, ce, y, x;

    se = sin(-ECL);
    ce = cos(ECL);

    sd = ce * sin(lat_ec) - se * cos(lat_ec)*sin(lon_ec);
    *dec_eq = asin(sd);

    y = ce*cos(lat_ec)*sin(lon_ec) + se*sin(lat_ec);
    x = cos(lat_ec)*cos(lon_ec);
    *ra_eq = atan2(y,x);

    if (partials!=NULL) {
	cd = sqrt(1.-sd*sd);
	/* form Jacobian matrix  (assume not at pole)*/
	partials[2][2] = -se*cos(lat_ec)*cos(lon_ec)/cd;
	partials[2][1] = (ce*cos(lat_ec)+se*sin(lat_ec)*sin(lon_ec)) / cd;
	partials[1][2] = cos(lat_ec)*partials[2][1] / cd;
	partials[1][1] = se*cos(lon_ec) / (cd*cd);
    }

    return;
}

/* And transform x,y,z from ecliptic to eq */
void xyz_ec_to_eq(double x_ec, double y_ec, double z_ec,
		  double *x_eq, double *y_eq, double *z_eq,
		  double **partials)
{
    double se, ce;

    se = sin(-ECL);
    ce = cos(ECL);

    *x_eq = x_ec;
    *y_eq = ce*y_ec + se*z_ec;
    *z_eq = -se*y_ec + ce*z_ec;

    if (partials!=NULL) {
	partials[1][1] = 1.;
	partials[1][2] = partials[1][3] =
	    partials[2][1] = partials[3][1] = 0.;
	partials[2][2] = partials[3][3] = ce;
	partials[2][3] = se;
	partials[3][2] = -se;
    }

    return;
}

/*** Go between ecliptic coords & tangent-plane system at
*** lat0, lon0 aligned with ecliptic.   keep the center
*** of projection as static variables to update as needed
****/

double old_lat0=-999., old_lon0, clat0, slat0, clon0, slon0;

double old_ra0=-999., old_dec0, cra0, sra0, cdec0, sdec0;

void check_latlon0(double lat0,
		   double lon0) 
{
    if (lat0 == old_lat0 && lon0==old_lon0) return;

    old_lat0 = lat0;
    old_lon0 = lon0;

    clat0 = cos(lat0);
    slat0 = sin(lat0);
    clon0 = cos(lon0);
    slon0 = sin(lon0);

    return;
}

void check_radec0(double ra0,
		  double dec0) 
{
    if (dec0 == old_dec0 && ra0==old_ra0) return;

    old_dec0 = dec0;
    old_ra0 = ra0;

    cdec0 = cos(dec0);
    sdec0 = sin(dec0);
    cra0 = cos(ra0);
    sra0 = sin(ra0);

    return;
}

/* first routine goes from ecliptic lat/lon to projected x/y angles*/
void ec_to_proj(double lat_ec,
		double lon_ec,
		double *x_proj,
		double *y_proj,
		double lat0,
		double lon0,
		double **partials)
{
    double clat, slat, cdlon, sdlon;
    double      xp, yp, zp;

    check_latlon0(lat0,lon0);

    cdlon = cos(lon_ec - lon0);
    sdlon = sin(lon_ec - lon0);
    clat  = cos(lat_ec);
    slat  = sin(lat_ec);

    xp = clat * sdlon;
    yp = clat0*slat - slat0*clat*cdlon;
    zp = slat0*slat + clat0*clat*cdlon;

    /* Go from cartesian to tangent-plane coords; don't worry here
     * about zp=0 which occurs 90 degrees from tangent point.
     */

    *x_proj = xp/zp;
    *y_proj = yp/zp;

    /* take the small-angle approximation for the partial deriv.
     * matrix, which should be good enough for anything within a
     * few degrees of the projection point.
     */

    if (partials!=NULL) {
	partials[1][2] = clat;
	partials[1][1] = partials[2][2] = 0.;
	partials[1][2] = 1.;
    }

    return;
}

/* first routine goes from ra/dec to projected x/y angles*/
void radec_to_proj(double ra,
		double dec,
		double *x_proj,
		double *y_proj,
		double ra0,
		double dec0,
		double **partials)
{
    double cdec, sdec, cdra, sdra;
    double      xp, yp, zp;

    check_radec0(ra0,dec0);

    cdra = cos(ra - ra0);
    sdra = sin(ra - ra0);
    cdec  = cos(dec);
    sdec  = sin(dec);

    xp = cdec * sdra;
    yp = cdec0*sdec - sdec0*cdec*cdra;
    zp = sdec0*sdec + cdec0*cdec*cdra;

    /* Go from cartesian to tangent-plane coords; don't worry here
     * about zp=0 which occurs 90 degrees from tangent point.
     */

    *x_proj = xp/zp;
    *y_proj = yp/zp;

    /* take the small-angle approximation for the partial deriv.
     * matrix, which should be good enough for anything within a
     * few degrees of the projection point.
     */

    if (partials!=NULL) {
	partials[1][2] = cdec;
	partials[1][1] = partials[2][2] = 0.;
	partials[1][2] = 1.;
    }

    return;
}

/** Now the inverse, from projected xy to ecliptic lat/lon **/
void proj_to_radec(double x_proj,
		   double y_proj,
		   double *ra,
		   double *dec,
		   double ra0,
		   double dec0,
		   double **partials){

    double zp;

    check_radec0(ra0,dec0);    

    zp = 1./sqrt(1 + x_proj*x_proj + y_proj*y_proj);

    *dec = asin( zp* (sdec0 + y_proj*cdec0) );

    *ra = ra0 + asin( x_proj * zp / cos(*dec) );

    /* take the small-angle approximation for the partial deriv.
     * matrix, which should be good enough for anything within a
     * few degrees of the projection point.
     */

    if (partials!=NULL) {
	partials[2][1] = 1./cos(*dec);
	partials[1][1] = partials[2][2] = 0.;
	partials[1][2] = 1.;
    }

    return;
}

/** Now the inverse, from projected xy to ecliptic lat/lon **/
void proj_to_ec(double x_proj,
		double y_proj,
		double *lat_ec,
		double *lon_ec,
		double lat0,
		double lon0,
		double **partials)
{
    double zp;

    check_latlon0(lat0,lon0);

    zp = 1./sqrt(1 + x_proj*x_proj + y_proj*y_proj);

    *lat_ec = asin( zp* (slat0 + y_proj*clat0) );

    *lon_ec = lon0 + asin( x_proj * zp / cos(*lat_ec) );

    /* take the small-angle approximation for the partial deriv.
     * matrix, which should be good enough for anything within a
     * few degrees of the projection point.
     */

    if (partials!=NULL) {
	partials[2][1] = 1./cos(*lat_ec);
	partials[1][1] = partials[2][2] = 0.;
	partials[1][2] = 1.;
    }

    return;
}


/** Next go from x,y,z in ecliptic orientation to x,y,z in tangent-point
 * orientiation.
 */
void
xyz_ec_to_proj( double x_ec, double y_ec, double z_ec,
		double *x_p, double *y_p, double *z_p,
		double lat0, double lon0,
		double **partials)
{

  check_latlon0(lat0,lon0);

  *x_p = -slon0 * x_ec
         + clon0 * y_ec;
  *y_p = -clon0*slat0*x_ec
         -slon0*slat0*y_ec
         +clat0*z_ec;
  *z_p =  clon0*clat0*x_ec
         +slon0*clat0*y_ec
         +slat0*z_ec;

  if (partials!=NULL) {
    partials[1][1] = -slon0;
    partials[1][2] = clon0;
    partials[1][3] = 0.;
    partials[2][1] = -clon0*slat0;
    partials[2][2] = -slon0*slat0;
    partials[2][3] = clat0;
    partials[3][1] = clon0*clat0;
    partials[3][2] = slon0*clat0;
    partials[3][3] = slat0;
  }

  return;
}

/** And finally from tangent x,y,z to ecliptic x,y,z
 */
void xyz_proj_to_ec( double x_p, double y_p, double z_p,
		double *x_ec, double *y_ec, double *z_ec,
		double lat0, double lon0,
		double **partials)
{

  check_latlon0(lat0,lon0);

  *x_ec =-slon0      *x_p
         -clon0*slat0*y_p
         +clon0*clat0*z_p;
  *y_ec = clon0      *x_p
         -slon0*slat0*y_p
         +slon0*clat0*z_p;
  *z_ec = clat0 * y_p
         +slat0 * z_p;

  if (partials!=NULL) {
    partials[1][1] =-slon0;
    partials[1][2] =-clon0*slat0;
    partials[1][3] = clon0*clat0;
    partials[2][1] = clon0;
    partials[2][2] =-slon0*slat0;
    partials[2][3] = slon0*clat0;
    partials[3][1] = 0.;
    partials[3][2] = clat0;
    partials[3][3] = slat0;
  }

  return;
}

/** Next go from x,y,z in equatorial orientation to x,y,z in tangent-point
 * orientiation.
 */ 
void xyz_eq_to_proj( double x_eq, double y_eq, double z_eq,
		double *x_p, double *y_p, double *z_p,
		double ra0, double dec0,
		double **partials){

   check_radec0(ra0, dec0);

    *x_p = -sra0 * x_eq
	+ cra0 * y_eq;
    *y_p = -cra0*sdec0*x_eq
	-sra0*sdec0*y_eq
         +cdec0*z_eq;
  *z_p =  cra0*cdec0*x_eq
         +sra0*cdec0*y_eq
         +sdec0*z_eq;

  if (partials!=NULL) {
    partials[1][1] = -sra0;
    partials[1][2] = cra0;
    partials[1][3] = 0.;
    partials[2][1] = -cra0*sdec0;
    partials[2][2] = -sra0*sdec0;
    partials[2][3] = cdec0;
    partials[3][1] = cra0*cdec0;
    partials[3][2] = sra0*cdec0;
    partials[3][3] = sdec0;
  }

  return;
}

/** And finally from tangent x,y,z to ecliptic x,y,z
 */
void xyz_proj_to_eq( double x_p, double y_p, double z_p,
		     double *x_eq, double *y_eq, double *z_eq,
		     double ra0, double dec0,
		     double **partials){

    check_radec0(ra0, dec0);

    *x_eq =-sra0      *x_p
	-cra0*sdec0*y_p
	+cra0*cdec0*z_p;
    *y_eq = cra0      *x_p
	-sra0*sdec0*y_p
	+sra0*cdec0*z_p;
    *z_eq = cdec0 * y_p
	+sdec0 * z_p;

    if (partials!=NULL) {
	partials[1][1] =-sra0;
	partials[1][2] =-cra0*sdec0;
	partials[1][3] = cra0*cdec0;
	partials[2][1] = cra0;
	partials[2][2] =-sra0*sdec0;
	partials[2][3] = sra0*cdec0;
	partials[3][1] = 0.;
	partials[3][2] = cdec0;
	partials[3][3] = sdec0;
    }

  return;
}

 
void
orbitElements(XVBASIS *xv,
	      ORBIT  *orb)
{

  double combinedMass; /* mass of Sun + mass of KBO */ 
  double epochTime; 

  /* We note that the origin of the inertial rectangular
     co-ordinate system is itself a dynamical center - the Sun. The
     co-ordinate system is aligned so that +x-axis points towards the
     vernal equinox and xy plane co-incides with the fundamental plane of
     the celestial co-ordinate system. This fundamental plane corresponds
     to the ecliptic plane since the orbit is heliocentric. Thus, the
     vectors r and v are referred to the ecliptic co-ordinate system by
     using equatorial to ecliptic transformation. */

  double R[4], V[4];
  double rMagnitude, velSquare, radVelDotProduct;
  double eccentricityVector[4], angularMomentum[4], ascendingNode[4];
  double semimajor, eccentricity, inclination, longitudeOfAscendingNode;
  double semiLatusRectum;
  double hMagnitude, ascendingNodeMagnitude; /* magnitude of angular momentum */
  double ascEccDotProduct, argumentOfPerifocus;
  double xBar, yBar;
  double cosE, sinE, eccentricAnomaly;
  /* E1 and E2 are used to decide the quadrant of Eccentric Anomaly */
  double meanAnomaly, meanMotion, timeOfPerifocalPassage;

  //combinedMass = GM * 1.00134 ; /* Alter GM to account for total SS mass*/
  combinedMass = GM0;
  epochTime=jd0;

  R[1] = xv->x; 
  R[2] = xv->y;
  R[3] = xv->z; 
  V[1] = xv->xdot;
  V[2] = xv->ydot;
  V[3] = xv->zdot;

  rMagnitude = sqrt(R[1]*R[1] + R[2]*R[2] + R[3]*R[3]);
  velSquare = V[1]*V[1] + V[2]*V[2] + V[3]*V[3];
  radVelDotProduct = R[1]*V[1] + R[2]*V[2] + R[3]*V[3];
  for(int k=1;k<=3;k++)
    {
      eccentricityVector[k] = (velSquare/combinedMass - 1/rMagnitude)*R[k]
	- (radVelDotProduct/combinedMass)*V[k];
    }

  /* Angular mom is cross product of rad and vel */
  angularMomentum[1] = R[2]*V[3] - R[3]*V[2];
  angularMomentum[2] = R[3]*V[1] - R[1]*V[3];
  angularMomentum[3] = R[1]*V[2] - R[2]*V[1];

  /* Ascending Node vector is cross product of k-unit vector and angular momentum
     vector. k = [0 0 1] and h = [hx, hy, hz] */
  ascendingNode[1] = -angularMomentum[2];
  ascendingNode[2] = angularMomentum[1];
  ascendingNode[3] = 0.0;

  semimajor = 1/(2/rMagnitude - velSquare/combinedMass);
  eccentricity = sqrt(eccentricityVector[1]*eccentricityVector[1] + 
		      eccentricityVector[2]*eccentricityVector[2] + 
		      eccentricityVector[3]*eccentricityVector[3]);
  semiLatusRectum = (angularMomentum[1]*angularMomentum[1] + 		    
		     angularMomentum[2]*angularMomentum[2] +	
		     angularMomentum[3]*angularMomentum[3])/combinedMass;
  /* p = h-square by mu */
  hMagnitude = sqrt( angularMomentum[1]*angularMomentum[1] + 		    
		     angularMomentum[2]*angularMomentum[2] +	
		     angularMomentum[3]*angularMomentum[3] );
  inclination = acos(angularMomentum[3]/hMagnitude); /* in radians here */
  ascendingNodeMagnitude = sqrt(ascendingNode[1]*ascendingNode[1] +
				ascendingNode[2]*ascendingNode[2] +
				ascendingNode[3]*ascendingNode[3]);
  longitudeOfAscendingNode = acos(ascendingNode[1]/ascendingNodeMagnitude);
  /* Capital Omega in radians here */
  if (ascendingNode[2] < 0) longitudeOfAscendingNode = 
			      2*PI - longitudeOfAscendingNode;
  /* ???could use atan2 here?? */
  ascEccDotProduct = ascendingNode[1]*eccentricityVector[1] +
    ascendingNode[2]*eccentricityVector[2] +
    ascendingNode[3]*eccentricityVector[3];
  argumentOfPerifocus = acos(ascEccDotProduct/
			     (ascendingNodeMagnitude*eccentricity)); 
  /* Small omega in radians here */
  if (eccentricityVector[3] < 0) argumentOfPerifocus = 
				   2*PI - argumentOfPerifocus;
  xBar = (semiLatusRectum - rMagnitude)/eccentricity;
  yBar = radVelDotProduct*sqrt(semiLatusRectum/combinedMass)/eccentricity;

  /* From here, we assume that the motion is elliptical */

  cosE = (xBar/semimajor) + eccentricity;
  sinE = yBar/(semimajor*sqrt(1-eccentricity*eccentricity));
  /* where semimajor*sqrt(1-eccentricity*eccentricity) is semiminor */
  eccentricAnomaly = atan2(sinE,cosE);

  meanAnomaly = eccentricAnomaly - eccentricity*sinE; /* radians */
  meanMotion = sqrt(combinedMass/(pow(semimajor,3.)));
  timeOfPerifocalPassage = epochTime - meanAnomaly/meanMotion; //DAY;
  /* This comes from M=n(t-T) where t is epoch time and T is time of perifocal
     passage, in days */
				
  orb->a=semimajor;
  orb->e=eccentricity; 
  orb->i=inclination/DTOR;/* in degrees now */  
  orb->lan = longitudeOfAscendingNode / DTOR; /*Long of ascending node*/
  orb->aop = argumentOfPerifocus / DTOR; /*argument of perihelion, degrees*/
  orb->T=timeOfPerifocalPassage;  /*Time of perihelion passage (JD)*/
  orb->meananom=meanAnomaly;  

  return; 

} /* function orbitElements ends */


/* JD routine stolen from skycalc: 
 * I have changed the structre & routine to allow floating d/h/m/s
 * entries.*/
/* Converts a date (structure) into a julian date.
   Only good for 1900 -- 2100. */

double 
date_to_jd(struct date_time date)
{
	short yr1=0, mo1=1;
	long jdzpt = 1720982, jdint, inter;
	double jd,jdfrac;
	short dint;

	if((date.y <= 1900) | (date.y >= 2100)) {
		printf("Date out of range.  1900 - 2100 only.\n");
		return(0.);
	}

	if(date.mo <= 2) {
		yr1 = -1;
		mo1 = 13;
	}

	dint = floor(date.d);	/*truncate to integer portion */
	jdint = 365.25*(date.y+yr1);  /* truncates */
	inter = 30.6001*(date.mo+mo1);
	jdint = jdint+inter+dint+jdzpt;
	jd = jdint;
	jdfrac=((date.s/60.+date.mn)/60.+date.h/24.)+date.d-dint;
	if(jdfrac < 0.5) {
		jdint--;
		jdfrac=jdfrac+0.5;
	}
	else jdfrac=jdfrac-0.5;
	jd=jdint+jdfrac;
	return(jd);
}

/******************************************************************
 * Get phase space from orbital elements 
 ******************************************************************/
void elements_to_xv(ORBIT *o,
		    double jd,
		    XVBASIS *xv)
{
    double r0[3], v0[3], r1[3], v1[3], r2[3], v2[3];
    double meanAnomaly;
    double c, s, t, dt;

    //double mu = GM*SSMASS;	/*use SSMASS, work in AU/yrs*/
    double mu = GM0;	/*use SSMASS, work in AU/yrs*/  

    if (o->e >= 1. || o->a <=0.) {
	fprintf(stderr,"elements_to_xv only for closed orbits now\n");
	exit(1);
    }

    /* get the eccentric Anomaly from the mean anomaly */
    meanAnomaly = (jd - o->T)*DAY * pow(o->a,-1.5) * sqrt(mu);
    /* Put it near 0 */
    t = floor (meanAnomaly / TPI);
    meanAnomaly -= t*TPI;
    /* Use bisection to find solution (adapt Numerical Recipes)*/  
    {
#define JMAX 40
#define TOLERANCE (DAY/24./3600.)
	double f, fmid, x1, x2, dx, xmid, rtb;
	int j;
	x1 = meanAnomaly - o->e;
	x2 = meanAnomaly + o->e;
	f   = x1 - o->e * sin(x1) - meanAnomaly;
	fmid= x2 - o->e * sin(x2) - meanAnomaly;
	if (f*fmid > 0.0) {
	    fprintf(stderr,"Error, eccentricAnomaly root not bracketed\n");
	    fprintf(stderr,"f, fmid %f %f\n",f,fmid);
	    exit(1);
	}

	rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
	for (j=1;j<=JMAX;j++) {
	    xmid=rtb+(dx *= 0.5);
	    fmid= xmid - o->e * sin(xmid) - meanAnomaly;
	    if (fmid <= 0.0) rtb=xmid;
	    if (fabs(dx) < TOLERANCE || fmid == 0.0) break;
	}
	if (j>=JMAX) {
	    fprintf(stderr,"eccentricAnomaly took too long\n");
	    exit(1);
	}
	meanAnomaly = rtb;
#undef JMAX
#undef TOLERANCE   
    }

    /*Coordinates and velocity in the system aligned with orbit: */
    c = cos(meanAnomaly);
    s = sin(meanAnomaly);
    r0[0] = o->a * (c - o->e);
    r0[1] = o->a * s * sqrt(1-o->e*o->e);
    dt = sqrt(pow(o->a,3.)/mu) * ( 1 - o->e*c);
    v0[0] = -o->a * s / dt;
    v0[1] = o->a * c * sqrt(1-o->e*o->e) / dt;

    /* Rotate about z to put perihelion at arg of peri */
    c = cos(o->aop*DTOR);  s = sin(o->aop*DTOR);
    r1[0] = r0[0]*c - r0[1]*s;
    r1[1] = r0[0]*s + r0[1]*c;
    v1[0] = v0[0]*c - v0[1]*s;
    v1[1] = v0[0]*s + v0[1]*c;

    /* Rotate about x axis to incline orbit */
    c = cos(o->i*DTOR);  s = sin(o->i*DTOR);
    r2[0] = r1[0];
    r2[1] = r1[1]*c;
    r2[2] = r1[1]*s;
    v2[0] = v1[0];
    v2[1] = v1[1]*c;
    v2[2] = v1[1]*s;

    /* Rotate about z axis to bring node to longitude */
    c = cos(o->lan*DTOR);  s = sin(o->lan*DTOR);
    xv->x = r2[0]*c - r2[1]*s;
    xv->y = r2[0]*s + r2[1]*c;
    xv->z = r2[2];
    xv->xdot = v2[0]*c - v2[1]*s;
    xv->ydot = v2[0]*s + v2[1]*c;
    xv->zdot = v2[2];

    return;
}


double machine_epsilon;

void keplerian(double gm, State state, 
	       double *a, double *e, double *i, double *longnode, double *argperi, double *meananom)
{
    double rxv_x, rxv_y, rxv_z, hs, h, parameter;
    double r, vs, rdotv, rdot, ecostrueanom, esintrueanom, cosnode, sinnode;
    double rcosu, rsinu, u, trueanom, eccanom;

    /* find direction of angular momentum vector */
    rxv_x = state.y * state.zd - state.z * state.yd;
    rxv_y = state.z * state.xd - state.x * state.zd;
    rxv_z = state.x * state.yd - state.y * state.xd;
    hs = rxv_x * rxv_x + rxv_y * rxv_y + rxv_z * rxv_z;
    h = sqrt(hs);

    r = sqrt(state.x * state.x + state.y * state.y + state.z * state.z);
    vs = state.xd * state.xd + state.yd * state.yd + state.zd * state.zd;
    /* v = sqrt(vs);  unnecessary */
    rdotv = state.x * state.xd + state.y * state.yd + state.z * state.zd;
    rdot = rdotv / r;
    parameter = hs / gm;

    *i = acos(rxv_z / h);

    if(rxv_x!=0.0 || rxv_y!=0.0) {
	*longnode = atan2(rxv_x, -rxv_y);
    } else {
	*longnode = 0.0;
    }

    *a = 1.0 / (2.0 / r - vs / gm);

    ecostrueanom = parameter / r - 1.0;
    esintrueanom = rdot * h / gm;
    *e = sqrt(ecostrueanom * ecostrueanom + esintrueanom * esintrueanom);

    if(esintrueanom!=0.0 || ecostrueanom!=0.0) {
	trueanom = atan2(esintrueanom, ecostrueanom);
    } else {
	trueanom = 0.0;
    }

    cosnode = cos(*longnode);
    sinnode = sin(*longnode);

    /* u is the argument of latitude */
    rcosu = state.x * cosnode + state.y * sinnode;
    rsinu = (state.y * cosnode - state.x * sinnode)/cos(*i);

    if(rsinu!=0.0 || rcosu!=0.0) {
	u = atan2(rsinu, rcosu);
    } else {
	u = 0.0;
    }

    *argperi = u - trueanom;

    eccanom = 2.0 * atan(sqrt((1.0 - *e)/(1.0 + *e)) * tan(trueanom/2.0));
    *meananom = eccanom - *e * sin(eccanom);

    return;
}

void cartesian(double gm, 
	       double a, double e, double i, double longnode, double argperi, double meananom, 
	       State *state)
{
    double meanmotion, cosE, sinE, foo;
    double x, y, z, xd, yd, zd;
    double xp, yp, zp, xdp, ydp, zdp;
    double cosw, sinw, cosi, sini, cosnode, sinnode;
    double E0, E1, E2, den;

    machine_epsilon = 1e-15;
    /* first compute eccentric anomaly */
    E0 = meananom; 
    do {
	E1 = meananom + e * sin(E0);
	E2 = meananom + e * sin(E1);

	den = E2 - 2.0*E1 + E0;
	if(fabs(den) > machine_epsilon) {
	    E0 = E0 - (E1-E0)*(E1-E0)/den;
	}
	else {
	    E0 = E2;
	    E2 = E1;
	}
    } while(fabs(E0-E2) > machine_epsilon);

    cosE = cos(E0);
    sinE = sin(E0);

    /* compute unrotated positions and velocities */
    foo = sqrt(1.0 - e*e);
    meanmotion = sqrt(gm/(a*a*a));
    x = a * (cosE - e);
    y = foo * a * sinE;
    z = 0.0;
    xd = -a * meanmotion * sinE / (1.0 - e * cosE);
    yd = foo * a * meanmotion * cosE / (1.0 - e * cosE);
    zd = 0.0;

    /* rotate by argument of perihelion in orbit plane*/
    cosw = cos(argperi);
    sinw = sin(argperi);
    xp = x * cosw - y * sinw;
    yp = x * sinw + y * cosw;
    zp = z;
    xdp = xd * cosw - yd * sinw;
    ydp = xd * sinw + yd * cosw;
    zdp = zd;

    /* rotate by inclination about x axis */
    cosi = cos(i);
    sini = sin(i);
    x = xp;
    y = yp * cosi - zp * sini;
    z = yp * sini + zp * cosi;
    xd = xdp;
    yd = ydp * cosi - zdp * sini;
    zd = ydp * sini + zdp * cosi;

    /* rotate by longitude of node about z axis */
    cosnode = cos(longnode);
    sinnode = sin(longnode);
    state->x = x * cosnode - y * sinnode;
    state->y = x * sinnode + y * cosnode;
    state->z = z;
    state->xd = xd * cosnode - yd * sinnode;
    state->yd = xd * sinnode + yd * cosnode;
    state->zd = zd;

    return;
}

double principal_value(double theta)
{
  theta -= 2*PI*floor(theta/(2*PI));
  return(theta);
}

double principal_value_0(double theta)
{
  theta -= 2*PI*floor(theta/(2*PI));
  if (theta<PI)
    return(theta);
  else
    return(theta-2.0*PI);
}
