/* orbfit1.c - program to fit orbital parameters to observed KBO data.
This is basically a template program which makes some simplifying
assumptions, with the goal of setting up for a more exact fit and for
obtaining estimates of the sizes of error ellipses that we will
encounter in the future.
 
4/9/99 gmb
*/

 
#include "orbfit.h"
#include <string.h>
//#include "ephem_types.h"
#include "rebound.h"
#include "assist.h"

//#define FNAMESIZE 256

/* Create the variables which define the coordinate system */
double	lat0, lon0;	     /* ecliptic lat & lon of tangent point */
double	ra0,  dec0;	     /* RA & Dec of tangent point */
double	xBary, yBary, zBary; /*Posn of barycenter in our system*/
double	jd0;		     /* Zeropoint of time scale */

/* Project KBO position onto tangent plane at a given time */
void kbo2d(struct assist_ephem* ephem,
	   PBASIS *pin, 
	   OBSERVATION *obs,
	   double *x, double dx[],
	   double *y, double dy[])
{
    double xk[3], vk[3], dxk[7], dyk[7], dzk[7];
    double xe, ye, ze;
    double invz;
    int i;
    double distance; 

    xe = obs->xe;
    ye = obs->ye;
    ze = obs->ze;

    /* Get preliminary KBO position */
    kbo3d_assist(ephem, pin,obs->obstime,xk,vk,dxk,dyk,dzk);

    /* At this point one should account for the light-travel time
       delay between observed time and kbo orbit position time.
       Calculate distance & time delay using naive position.
    */ 
    distance=sqrt( (xk[0]-xe)*(xk[0]-xe) + (xk[1]-ye)*(xk[1]-ye)
		   + (xk[2]-ze)*(xk[2]-ze) );

    kbo3d_assist(ephem, pin,obs->obstime-distance/SPEED_OF_LIGHT, 
	  xk,vk,dxk,dyk,dzk);
 
    invz = 1./(xk[2]-ze);
    *x = (xk[0] - xe)*invz;
    *y = (xk[1] - ye)*invz;
 
    for (i=1; i<=6; i++) {
	if (dx!=NULL) dx[i] = dxk[i]*invz - (xk[0]-xe)*dzk[i]*invz*invz;
	if (dy!=NULL) dy[i] = dyk[i]*invz - (xk[1]-ye)*dzk[i]*invz*invz;
    }
    /* incorporate derivative w.r.t simplified time delay */
    if (dx!=NULL) dx[5] += vk[0]*invz*invz/SPEED_OF_LIGHT;
    if (dx!=NULL) dy[5] += vk[1]*invz*invz/SPEED_OF_LIGHT;
 
    return;
}

/* Determine KBO position on spherical sky from observatory */
void kbo_sphere(struct assist_ephem* ephem,
		PBASIS *pin, 
		OBSERVATION *obs,
		double *x, double dx[],
		double *y, double dy[],
		double *z, double dz[])
{
    double xk[39], vk[39], dxk[7], dyk[7], dzk[7];
    double ddist[3];
    double xe,ye,ze;
    double invd;
    int i;
    int niter=2;

    /* Get the Earth position */
    xe = obs->xe;
    ye = obs->ye;
    ze = obs->ze;

    /* Account for the light time delay between observed time 
       and object position time.
    */ 
    
    double distance = 0.0;

    // Set the number of light time interations.
    for (int j=0; j<niter; j++) {

	kbo3d_assist(ephem, pin, obs->obstime-distance/SPEED_OF_LIGHT,
	      xk, vk, dxk, dyk, dzk);

	/* Recompute the distance, because we need the correct value below */
	distance=sqrt( (xk[0]-xe)*(xk[0]-xe) + (xk[1]-ye)*(xk[1]-ye)
		       + (xk[2]-ze)*(xk[2]-ze) );
    }
  
    invd = 1./distance;
    *x = (xk[0] - xe)*invd;
    *y = (xk[1] - ye)*invd;
    *z = (xk[2] - ze)*invd;

    for (i=1; i<=6; i++) {

	/* Derivative of topocentric distance w.r.t. parameters */
	ddist[i] = (xk[0] - xe)*invd*dxk[i] + (xk[1] - ye)*invd*dyk[i] + (xk[2] - ze)*invd*dzk[i];

	if (dx!=NULL) dx[i] = dxk[i]*invd - (xk[0]-xe)*invd*invd*ddist[i] - vk[0]*ddist[i]*invd/SPEED_OF_LIGHT;
	if (dy!=NULL) dy[i] = dyk[i]*invd - (xk[1]-ye)*invd*invd*ddist[i] - vk[1]*ddist[i]*invd/SPEED_OF_LIGHT;
	if (dz!=NULL) dz[i] = dzk[i]*invd - (xk[2]-ze)*invd*invd*ddist[i] - vk[2]*ddist[i]*invd/SPEED_OF_LIGHT;

    }

    return;
}

/* Linearized version of the 2d projection of orbital position.  Only
   leading-order terms for each parameter's derivative are given.
   The derivative with respect to gdot is inserted here - note that
   it is much smaller than the others, and is the only part of the
   derivatives that has any dependence upon the incoming PBASIS.
*/
void kbo2d_linear(PBASIS *pin,
		  OBSERVATION *obs,
		  double *x, double dx[],
		  double *y, double dy[])
{
    double        xe,ye,ze,t_emit;
    int           i;
 
    /* Get the Earth position */
    xe = obs->xe;
    ye = obs->ye;
    ze = obs->ze;

    /* Account for light-travel time differentially to leading order
     * by retarding KBO motion by z component of Earth's position. 
     * Note ignoring acceleration here.
     */
    t_emit = obs->obstime - ze/SPEED_OF_LIGHT;
    *x = pin->a + pin->adot*t_emit - pin->g * xe
	- pin->gdot * (pin->adot*t_emit*t_emit - pin->g*xe*t_emit);
    *y = pin->b + pin->bdot*t_emit - pin->g * ye
	- pin->gdot * (pin->bdot*t_emit*t_emit - pin->g*ye*t_emit);  
 
    if (dx!= NULL && dy!=NULL) {
	for (i=1; i<=6; i++) dx[i] = dy[i]=0.;
 
	dx[1] = dy[3] = 1.;
	dx[2] = dy[4] = t_emit;
	dx[5] = -xe;
	dy[5] = -ye;
	dx[6] = -(pin->adot*t_emit*t_emit - pin->g*xe*t_emit);
	dy[6] = -(pin->bdot*t_emit*t_emit - pin->g*ye*t_emit);
    }
 
    return;
}

/* Read a line from specified input stream, skipping blank or
 * commented (#) lines.  These are echoed to fpout if it's non-null.
 * Returns NULL at EOF.
 */
char* fgets_nocomment(char *inbuff, int length, 
		      FILE *fpin, FILE *fpout)
{
    char test[10];
    while (1) {
	if (fgets(inbuff,length,fpin)==NULL) return(NULL);

	/* Skip blank lines or comments */
	if (sscanf(inbuff,"%8s",test)>0
	    && test[0]!='#') return(inbuff);

	/* Echo comments to output if desired */
	if (fpout!=NULL) 
	    fputs(inbuff,fpout);
    }
}

/* read a string as a line of an observation that has a simple
 * format.
 *
 * desig   = string of up to 256 characters
 * jd_tdb  = TDB Julian Date
 * ra      = RA in decimal degrees
 * dec     = Dec in decimal degrees
 * ra_unc  = RA  uncertainty in arcsec
 * dec_dec = Dec uncertainty in arcsec'
 * xe, ye, ze
 *         = equatorial geocentric position of observatory
 * obscode = MPC observatory code (not used)
 */
int scan_observation_simple(char *inbuff,
			    OBSERVATION *obs){

  double jd_tdb;
  double ra, dec;
  double xe, ye, ze;
  double ra_unc, dec_unc;
  char desig[256];
  char obsCode[256];

  sscanf(inbuff, "%s %lf %lf %lf %lf %lf %lf %lf %lf %s", desig, &jd_tdb, &ra, &ra_unc, &dec, &dec_unc, &xe, &ye, &ze, obsCode);
  obs->obscode = malloc(11*sizeof(char));
  strcpy(obs->obscode, obsCode);

  //geocenter_ssbary(jd_tdb, earth_xyz);

  obs->obstime = jd_tdb;
  obs->ra = ra;
  obs->dec = dec;

  obs->desig = malloc(11*sizeof(char));  
  sscanf(desig, "%s", obs->desig);

  obs->thetax = DTOR * ra;
  obs->thetay = DTOR * dec;

  obs->dthetax = ra_unc;
  obs->dthetay = dec_unc;
  obs->dthetax *= ARCSEC;
  obs->dthetay *= ARCSEC;

  
  obs->xe = xe;
  obs->ye = ye;
  obs->ze = ze;
  
  obs->reject=0;

  return(0);

}

/* read list of observations from a file */
/* Blank & commented lines are skipped */
/* Input file is stdin if fname==NULL */
int read_observations_simple(OBSERVATION obsarray[], char *fname, int *nobs)
{

    FILE *fptr;
    OBSERVATION  *obs;
    char	inbuff[256];
    double elat,elon;

    // 0. Open the named file.
    if (fname==NULL)
	fptr = stdin;
    else if ( (fptr=fopen(fname,"r"))==NULL) {
	fprintf(stderr,"Error opening observations file %s\n",fname);
	fflush(stderr);
	exit(1);
    }

    // 1. Read the observations, without transformatioan.
    *nobs=0;
    while ( fgets_nocomment(inbuff,255,fptr,NULL)!=NULL) {
	if ( scan_observation_simple(inbuff, &(obsarray[*nobs]))) {
	    fprintf(stderr,"Quitting on format error\n");
	    exit(1);
	}

	obs = &(obsarray[*nobs]);

	(*nobs)++;

    }

    // 2. Select an observation that defines the tangent plane
    //    and nearest-day reference time.

    int nref = (*nobs)-1;
    //int nref = 0;
    obs = &(obsarray[nref]);
    
    double xec, yec, zec;

    // Save the RA/Dec values
    eq_to_ec(obs->thetax,obs->thetay,&elat,&elon,NULL);
    //printf("%lf %lf\n", elon, elat);
    double x0 = cos(elat)*cos(elon);
    double y0 = cos(elat)*sin(elon);
    double z0 = sin(elat);

    //printf("%lf %lf %lf\n", x0, y0, z0);

    lat0 = elat;
    lon0 = elon;
    jd0 = obs->obstime;

    xBary = obs->xe;
    yBary = obs->ye;
    zBary = obs->ze;    

    /* Negate the vector to make it earth->SSBARY*/
    /* And rotate equatorial into the tangent-point coords */
    xBary *= -1.;  yBary *= -1.;  zBary *= -1.;

    // Switch to equatorial projection coordinates
    xyz_eq_to_ec(xBary, yBary, zBary, &xec, &yec, &zec,NULL);
    //printf("%lf %lf %lf\n", xec, yec, zec);
    xyz_ec_to_proj(xec, yec, zec, &xBary, &yBary, &zBary, lat0, lon0, NULL);

    // 3. Transform the observations based on that selection.

    for(int i=0; i<*nobs; i++){

	obs = &(obsarray[i]);

	obs->obstime = (obs->obstime-jd0); //*DAY;

	eq_to_ec(obs->thetax,obs->thetay,&elat,&elon,NULL);
	// No need to change to ecliptic
    
	double xec = cos(elon)*cos(elat);
	double yec = sin(elon)*cos(elat);
	double zec = sin(elat);
	double xp, yp, zp;
	// Change this to RA/Dec to proj
	xyz_ec_to_proj( xec, yec, zec, &xp, &yp, &zp, lat0, lon0, NULL);

	// Change this to RA/Dec to proj
	ec_to_proj(elat,elon,&(obs->thetax),&(obs->thetay), lat0,lon0,NULL);
	obs->x = xp; obs->y = yp; obs->z = zp;
	//printf("%le %le %le\n", obs->x, obs->y, obs->z);

	/* convert to tangent-point coord system */
	/* via ecliptic */
	double x1 = obs->xe;
	double y1 = obs->ye;
	double z1 = obs->ze;
	double xtmp, ytmp, ztmp;
	// Change this to equatorial to proj
	xyz_eq_to_ec(x1, y1, z1, &xtmp, &ytmp, &ztmp,NULL);
	xyz_ec_to_proj(xtmp, ytmp, ztmp, &x1, &y1, &z1, lat0, lon0, NULL);

	/* Translate to our origin */
	x1 += xBary;
	y1 += yBary;
	z1 += zBary;

	obs->xe = x1;
	obs->ye = y1;
	obs->ze = z1;

	//printf("%le %le %le\n", obs->xe, obs->ye, obs->ze);

    }

    if (fname!=NULL) fclose(fptr);
    return(0);
}

int transform_observations(OBSERVATION *obsarray, int nobs)
{
    OBSERVATION  *obs;
    double elat,elon;

    // 1. Select an observation that defines the tangent plane
    //    and the reference time.

    int nref = (nobs)-1;
    //int nref = 0;
    obs = &(obsarray[nref]);
    
    double xec, yec, zec;

    // Save the RA/Dec values
    eq_to_ec(obs->thetax,obs->thetay,&elat,&elon,NULL);

    lat0 = elat;
    lon0 = elon;
    jd0 = obs->obstime;

    printf("jd0: %lf\n", jd0);

    xBary = obs->xe;
    yBary = obs->ye;
    zBary = obs->ze;    

    /* Negate the vector to make it earth->SSBARY*/
    /* And rotate equatorial into the tangent-point coords */
    xBary *= -1.;  yBary *= -1.;  zBary *= -1.;

    // Switch to equatorial projection coordinates
    xyz_eq_to_ec(xBary, yBary, zBary, &xec, &yec, &zec,NULL);
    xyz_ec_to_proj(xec, yec, zec, &xBary, &yBary, &zBary, lat0, lon0, NULL);
    printf("%lf %lf %lf\n", xBary, yBary, zBary);

    // 2. Transform the observations based on that selection.

    for(int i=0; i<nobs; i++){

	obs = &(obsarray[i]);

	obs->obstime = (obs->obstime-jd0); //*DAY;

	eq_to_ec(obs->thetax,obs->thetay,&elat,&elon,NULL);
	// No need to change to ecliptic
    
	double xec = cos(elon)*cos(elat);
	double yec = sin(elon)*cos(elat);
	double zec = sin(elat);
	printf("%lf %lf %lf %lf %lf %lf %lf\n", obs->thetax, obs->thetay, elat, elon, xec, yec, zec);
	double xp, yp, zp;
	// Change this to RA/Dec to proj
	xyz_ec_to_proj( xec, yec, zec, &xp, &yp, &zp, lat0, lon0, NULL);
	obs->x = xp; obs->y = yp; obs->z = zp;

	// Change this to RA/Dec to proj
	ec_to_proj(elat,elon,&(obs->thetax),&(obs->thetay), lat0,lon0,NULL);

	/* convert to tangent-point coord system */
	/* via ecliptic */
	double x1 = obs->xe;
	double y1 = obs->ye;
	double z1 = obs->ze;
	double xtmp, ytmp, ztmp;
	// Change this to equatorial to proj
	xyz_eq_to_ec(x1, y1, z1, &xtmp, &ytmp, &ztmp,NULL);
	xyz_ec_to_proj(xtmp, ytmp, ztmp, &x1, &y1, &z1, lat0, lon0, NULL);

	/* Translate to our origin */
	x1 += xBary;
	y1 += yBary;
	z1 += zBary;

	obs->xe = x1;
	obs->ye = y1;
	obs->ze = z1;	

    }

    return(0);
}

/* Read an orbit in alpha, beta, gamma format & covar matrix from file.
   Also reads the coordinate frame specifications.
*/
int
read_abg(char *fname, 
	 PBASIS *p, 
	 double **covar)
{
  FILE *fptr;
  int  i;
  char	inbuff[256];

  if (fname==NULL)
    fptr = stdin;
  else if ( (fptr=fopen(fname,"r"))==NULL) {
    fprintf(stderr,"Error opening a/b/g orbit file %s\n",fname);
    return(1);
  }

  /* Skipping comments, read the a/b/g specs*/
  if (fgets_nocomment(inbuff,255,fptr,NULL)==NULL) {
    fprintf(stderr,"Data missing from a/b/g/ data file.\n");
    return(1);
  }

  if (sscanf(inbuff, "%lf %lf %lf %lf %lf %lf",
	     &(p->a),&(p->adot),&(p->b),&(p->bdot),
	     &(p->g),&(p->gdot)) != 6) {
    fprintf(stderr,"Error reading a/b/g data\n");
    return(1);
  }

  for (i=1; i<=6; i++) {
    if (fgets_nocomment(inbuff,255,fptr,NULL)==NULL) {
      fprintf(stderr,"Data missing from a/b/g/ covariance.\n");
      return(1);
    }

    if (sscanf(inbuff, "%lf %lf %lf %lf %lf %lf",
	       &covar[i][1],&covar[i][2],&covar[i][3],
	       &covar[i][4],&covar[i][5],&covar[i][6]) != 6) {
      fprintf(stderr,"Error reading a/b/g covariance\n");
      return(1);
    }
  }

  /* Now read the coordinate system info */
  if (fgets_nocomment(inbuff,255,fptr,NULL)==NULL) {
    fprintf(stderr,"Data missing from a/b/g/ data file.\n");
    return(1);
  }

  if (sscanf(inbuff, "%lf %lf %lf %lf %lf %lf",
	     &lat0, &lon0, &xBary, &yBary, &zBary, &jd0) != 6) {
    fprintf(stderr,"Error reading coord system info\n");
    return(1);
  }
  lat0 *= DTOR;
  lon0 *= DTOR;

  if (fname!=NULL) fclose(fptr);
  return(0);
}

/* Take a set of observations and make a preliminary fit using the
 * linear model of the orbit.  Then fill in zero for gdot, and return
 * the fit and an uncertainty matrix.  The gdot term of uncertainty
 * matrix is set to a nominal value.
 * Note covar is assumed to be 6x6 1-indexed matrix a la Numerical Recipes.
 */
void prelim_fit(OBSERVATION obsarray[],
		int nobs,
		PBASIS *pout,
		double **covar)
{
    int invert_matrix(double **in, double **out, int dim);
    double *beta, *soln, **alpha;
    int	i,j,k;
    double x,y,wtx,wty,*dx,*dy;
  

    beta=dvector(1,6);
    alpha=dmatrix(1,6,1,6);
    soln=dvector(1,6);
    dx=dvector(1,6);
    dy=dvector(1,6);

    /* clear all vectors/matrices*/
    for (i=1; i<=6; i++) {
	beta[i]=soln[i]=0.;
	for (j=1; j<=6; j++) alpha[i][j]=0.;
    }

    /*Collect the requisite sums*/
    for (i=0; i<nobs; i++) {
	if (obsarray[i].reject==0) {

	    wtx = 1./obsarray[i].dthetax;
	    wtx *= wtx;
	    wty = 1./obsarray[i].dthetay;
	    wty *= wty;

	    kbo2d_linear(pout,&(obsarray[i]),&x,dx,&y,dy);

	    /* Note that the dx[6] and dy[6] terms will only make
	     * even the least sense if the pout->g and adot were set to
	     * some sensible value beforehand.
	     */

	    for (j=1; j<=6; j++) {
		beta[j] += obsarray[i].thetax * dx[j] * wtx;
		beta[j] += obsarray[i].thetay * dy[j] * wty;
		for (k=1; k<=j; k++) {
		    alpha[j][k] += dx[j]*dx[k]*wtx;
		    alpha[j][k] += dy[j]*dy[k]*wty; 
		}
	    }
	}
    }

    /* Symmetrize and invert the alpha matrix to give covar.  Note
     * that I'm only going to bother with the first 5 params.
     */
    for (i=1; i<=5; i++)
	for (j=1; j<i; j++)
	    alpha[j][i]=alpha[i][j];


    if (invert_matrix(alpha,covar,5)) {
	/* some failure in the inversion...*/
	fprintf(stderr,"Error inverting the alpha matrix\n");
	exit(1);
    }

    /* Now multiply matrices to get the solution vector */
    for (i=1; i<=5; i++) {
	soln[i]=0.;
	for (j=1; j<=5; j++)
	    soln[i] += covar[i][j]*beta[j];
    }

    /* fill in the PBASIS structure */
    pout->a =    soln[1];
    pout->adot = soln[2];
    pout->b =    soln[3];
    pout->bdot = soln[4];
    pout->g    = soln[5];
    pout->gdot = 0.;              /*assuming that prelim fit has no info here*/

    /*Set the gdot parts of the covariance matrix to nominal values */
    for (i=1; i<6; i++)
	covar[i][6]=covar[6][i]=0.;
    covar[6][6]=0.1*TPI*TPI*pow(pout->g,3.);

    //printf("prelim_fit exit\n");
    //fflush(stdout);
    
    free_dvector(dx,1,6);
    free_dvector(dy,1,6);
    free_dvector(soln,1,6);
    free_dvector(beta,1,6);
    free_dmatrix(alpha,1,6,1,6);

    return;
}

/* Routine to predict position and uncertainty at any time, given
 * a PBASIS fit and a sigma matrix.
 */
void
predict_posn(struct assist_ephem* ephem,
	     PBASIS *pin,
             double **covar,
             OBSERVATION *obs,
             double **sigxy)    /*this holds xy error matrix*/
{
  int   i,j;
  double *dx,*dy;

  dx=dvector(1,6);
  dy=dvector(1,6);

  /*using input time & obscode, put xy position into OBSERVATION struct*/
  kbo2d(ephem, pin,obs,
	&(obs->thetax),dx,&(obs->thetay),dy);

  /* project the covariance matrix */
  /* skip if not desired */
  if (sigxy!=NULL && covar!=NULL) {
    sigxy[1][1]=sigxy[1][2]=sigxy[2][2]=0;
    for (i=1; i<=6; i++) { 
      for (j=1; j<=6; j++) {
	sigxy[1][1] += dx[i]*covar[i][j]*dx[j];
	sigxy[1][2] += dx[i]*covar[i][j]*dy[j];
	sigxy[2][2] += dy[i]*covar[i][j]*dy[j];
      }
    }
    sigxy[2][1]=sigxy[1][2];
  }

  free_dvector(dx,1,6);
  free_dvector(dy,1,6);
  return;
}

/* Routine to predict position and uncertainty at any time, given
 * a PBASIS fit and a sigma matrix.
 */
void
predict_posn_linear(PBASIS *pin,
             double **covar,
             OBSERVATION *obs,
             double **sigxy)    /*this holds xy error matrix*/
{
  int   i,j;
  double *dx,*dy;

  dx=dvector(1,6);
  dy=dvector(1,6);

  /*using input time & obscode, put xy position into OBSERVATION struct*/
  kbo2d_linear(pin,obs,
	&(obs->thetax),dx,&(obs->thetay),dy);

  /* project the covariance matrix */
  /* skip if not desired */
  if (sigxy!=NULL && covar!=NULL) {
    sigxy[1][1]=sigxy[1][2]=sigxy[2][2]=0;
    for (i=1; i<=6; i++) { 
      for (j=1; j<=6; j++) {
	sigxy[1][1] += dx[i]*covar[i][j]*dx[j];
	sigxy[1][2] += dx[i]*covar[i][j]*dy[j];
	sigxy[2][2] += dy[i]*covar[i][j]*dy[j];
      }
    }
    sigxy[2][1]=sigxy[1][2];
  }

  free_dvector(dx,1,6);
  free_dvector(dy,1,6);
  return;
}

/* Invert a double matrix, 1-indexed of size dim
 * from Numerical Recipes.  Input matrix is destroyed.
 */
int invert_matrix(double **in, double **out, int dim)
{
  extern void ludcmp(double **a, int n, int *indx, double *d);
  extern void lubksb(double **a, int n, int *indx, double *b);

  int   *indx,i,j;
  double *tvec,det;

  tvec = dvector(1,dim);
  indx = ivector(1,dim);
  ludcmp(in,dim,indx,&det);

  for (j=1; j<=dim; j++) {
    for (i=1; i<=dim; i++) tvec[i]=0.;
    tvec[j] = 1.0;
    lubksb(in,dim,indx,tvec);
    for (i=1; i<=dim; i++) out[i][j]=tvec[i];
  }

  free_ivector(indx,1,6);
  free_dvector(tvec,1,6);
  return(0);
}

void
print_matrix(FILE *fptr, double **matrix, int xdim, int ydim)
{
  int   i,j;
  for (i=1; i<=ydim; i++) {
    for (j=1; j<=xdim; j++) {
      fprintf(fptr,"%11.4e ",matrix[i][j]);
    }
    fprintf(fptr,"\n");
  }
  return;
}

/* Routine to give the t=0 position & velocity of KBO in
 * ecliptic system centered at SSBARY */
// Make a different version of this that is in equatorial
// coordinates
void
pbasis_to_bary(PBASIS *p,
	       XVBASIS *xv,
	       double **partials)
{
  double xProj,yProj,zProj,xdotProj,ydotProj,zdotProj; 
  double z0; 
  double **dProj_dp, **dEc_dProj;	/*intermediate partial matrices */
  int  i,j;

  if (partials!=NULL) {
    dProj_dp = dmatrix(1,6,1,6);
    dEc_dProj = dmatrix(1,6,1,6);
    for (i=1; i<=6; i++)
      for (j=1; j<=6; j++)
	dProj_dp[i][j]=dEc_dProj[i][j]=0.;
  }

  z0 = 1./p->g;
  xProj=p->a*z0;
  yProj=p->b*z0;
  zProj=z0;
  xdotProj=p->adot *z0;
  ydotProj=p->bdot *z0;
  zdotProj=p->gdot *z0;

  if (partials!=NULL) {
    dProj_dp[1][1] = z0;
    dProj_dp[1][5] = -p->a*z0*z0;
    dProj_dp[2][3] = z0;
    dProj_dp[2][5] = -p->b*z0*z0;
    dProj_dp[3][5] = -z0*z0;
    dProj_dp[4][2] = z0;
    dProj_dp[4][5] = -p->adot*z0*z0;
    dProj_dp[5][4] = z0;
    dProj_dp[5][5] = -p->bdot*z0*z0;
    dProj_dp[6][6] = z0;
    dProj_dp[6][5] = -p->gdot*z0*z0;
  }

  /* Bring position into ecliptic system relative to SSBARY */
  if (partials!=NULL)
    xyz_proj_to_ec(xProj-xBary, yProj-yBary, zProj-zBary,
		   &(xv->x), &(xv->y), &(xv->z),
		   lat0, lon0, dEc_dProj);
  else 
    xyz_proj_to_ec(xProj-xBary, yProj-yBary, zProj-zBary,
		   &(xv->x), &(xv->y), &(xv->z),
		   lat0, lon0, NULL);

  /* Repeat for velocity */
  xyz_proj_to_ec(xdotProj, ydotProj, zdotProj,
		 &(xv->xdot), &(xv->ydot), &(xv->zdot),
		 lat0, lon0, NULL);

  /* Calculate the overall Jacobian if desired*/
  if (partials!=NULL) {
    /* Replicate the spatial part of dEc_dProj into velocities */
    for (i=4; i<=6; i++)
      for (j=4; j<=6; j++)
	dEc_dProj[i][j] = dEc_dProj[i-3][j-3];

    matrix_multiply(dEc_dProj,dProj_dp,partials,6,6,6,6);

    free_dmatrix(dProj_dp,1,6,1,6);
    free_dmatrix(dEc_dProj,1,6,1,6);
  }

  return;   
}

/* Routine to give the t=0 position & velocity of KBO in
 * equatorial system centered at SSBARY */
void
pbasis_to_bary_eq(PBASIS *p,
	       XVBASIS *xv,
	       double **partials)
{
  double xProj,yProj,zProj,xdotProj,ydotProj,zdotProj; 
  double z0; 
  double **dProj_dp, **dEc_dProj;	/*intermediate partial matrices */
  double **dEq_dEc, **partials_tmp;     /*intermediate partial matrices */
  XVBASIS xv_tmp;
  int  i,j;

  if (partials!=NULL) {
    dProj_dp = dmatrix(1,6,1,6);
    dEc_dProj = dmatrix(1,6,1,6);
    dEq_dEc = dmatrix(1,6,1,6);
    partials_tmp = dmatrix(1,6,1,6);        
    for (i=1; i<=6; i++)
	for (j=1; j<=6; j++){
	    dProj_dp[i][j]=dEc_dProj[i][j]=0.;
	    dEq_dEc[i][j]=0.;
	    partials_tmp[i][j]=0.;	    
	}
  }

  z0 = 1./p->g;
  xProj=p->a*z0;
  yProj=p->b*z0;
  zProj=z0;
  xdotProj=p->adot *z0;
  ydotProj=p->bdot *z0;
  zdotProj=p->gdot *z0;

  if (partials!=NULL) {
    dProj_dp[1][1] = z0;
    dProj_dp[1][5] = -p->a*z0*z0;
    dProj_dp[2][3] = z0;
    dProj_dp[2][5] = -p->b*z0*z0;
    dProj_dp[3][5] = -z0*z0;
    dProj_dp[4][2] = z0;
    dProj_dp[4][5] = -p->adot*z0*z0;
    dProj_dp[5][4] = z0;
    dProj_dp[5][5] = -p->bdot*z0*z0;
    dProj_dp[6][6] = z0;
    dProj_dp[6][5] = -p->gdot*z0*z0;
  }

  /* Bring position into ecliptic system relative to SSBARY */
  if (partials!=NULL){

      xyz_proj_to_ec(xProj-xBary, yProj-yBary, zProj-zBary,
		     &(xv_tmp.x), &(xv_tmp.y), &(xv_tmp.z),
		     lat0, lon0, dEc_dProj);
      xyz_ec_to_eq(xv_tmp.x, xv_tmp.y, xv_tmp.z,
		&(xv->x), &(xv->y), &(xv->z),
		dEq_dEc);

  } else {

      xyz_proj_to_ec(xProj-xBary, yProj-yBary, zProj-zBary,
		     &(xv_tmp.x), &(xv_tmp.y), &(xv_tmp.z),
		     lat0, lon0, NULL);

      xyz_ec_to_eq(xv_tmp.x, xv_tmp.y, xv_tmp.z,
		&(xv->x), &(xv->y), &(xv->z),
		NULL);

  }
  /* Repeat for velocity */
  xyz_proj_to_ec(xdotProj, ydotProj, zdotProj,
		 &(xv_tmp.xdot), &(xv_tmp.ydot), &(xv_tmp.zdot),
		 lat0, lon0, NULL);
  xyz_ec_to_eq(xv_tmp.xdot, xv_tmp.ydot, xv_tmp.zdot,
	    &(xv->xdot), &(xv->ydot), &(xv->zdot),
	    NULL);

  /* Calculate the overall Jacobian if desired*/
  if (partials!=NULL) {
    /* Replicate the spatial part of dEc_dProj into velocities */
      for (i=4; i<=6; i++){
	  for (j=4; j<=6; j++){
	      dEc_dProj[i][j] = dEc_dProj[i-3][j-3];
	      dEq_dEc[i][j] = dEq_dEc[i-3][j-3];
	  }
      }

      matrix_multiply(dEc_dProj,dProj_dp,partials_tmp,6,6,6,6);
      matrix_multiply(dEq_dEc,partials_tmp,partials,6,6,6,6);      

      free_dmatrix(dProj_dp,1,6,1,6);
      free_dmatrix(dEc_dProj,1,6,1,6);
      free_dmatrix(dEq_dEc,1,6,1,6);
      free_dmatrix(partials_tmp,1,6,1,6);            
  }

  return;   
}

/* Multiply matrix m1 x m2.  m1 is x1 cols by y1 rows, etc.
   A dumb implementation.  Note that SECOND index of
   matrices are varying most rapidly, along rows.
*/
void
matrix_multiply(double **m1, double **m2,
		double **mout,
		int x1, int y1, int x2, int y2)
{
  int	i,j,k;
  double	sum;

  if (x1!=y2) {
    fprintf(stderr,"Trying to multiply mismatched matrices, "
	    " %d x %d  times %d x %d\n", y1,x1,y2,x2);
    exit(1);
  }
  for (i=1; i<=y1; i++)
    for (k=1; k<=x2; k++) {
      sum = 0.;
      for (j=1; j<=y2; j++)
	sum += m1[i][j]*m2[j][k];
      mout[i][k] = sum;
    }

  return;
}

/* Remap covariance matrix from one basis to another, using
   the partial deriv matrix **derivs.
   kin = dimension on input, kout = dimension on output.
   Just calculates deriv_T x covar_in x deriv.
*/
void
covar_map(double **covar_in, 
	  double **derivs, 
	  double **covar_out,
	  int kin, int kout)
{
  int	i,j,m,n;
  double	sum;

  for (i=1; i<=kout; i++)
    for (j=1; j<=kout; j++) {
      sum = 0.;
      for (m=1; m<=kin; m++)
	for (n=1; n<=kin; n++)
	sum += derivs[i][m]*derivs[j][n]*covar_in[m][n];
      covar_out[i][j] = sum;
    }

  return;
}
 
/* calculate acceleration given barycentric coords x */
/* Here is a more sophisticated version that includes perturbations*/
/**** ??? Note that we could cache the planet positions if this turns
*** out to be slower than desired.
***/
// Need to pass in the parameter and the components of the unit vector
// to the perturber in the target-specific coordinate system
void
accel_nbody(struct assist_ephem* ephem,
	    double *xx,
	    double *a,
	    double t){

  double x,y,z,acc,r2;
  double xb, yb, zb;

  // Need to put Pluto in here from DE405, but make sure Pluto does not self perturb.

  static int bodies[]={SUN, MERCURY, VENUS, EARTH, MOON, MARS, JUPITER, SATURN, URANUS, NEPTUNE, CERES, PALLAS, VESTA};
  //static int bodies[]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  //static int nbodies = 11;
  static int nbodies = 13;   
  //static double mass[]={1.0, 1.660114153e-7, 2.447838288e-6, 3.003489615e-6, 3.227156038e-7, 9.547919384e-4, 2.858859806661e-4, 4.366244e-5, 5.1513890e-5};

  double GM;
  int i,j;
  double vxyz[3];

  for (j=0; j<3; j++) a[j]=0.;

  for (j=0; j<3; j++) vxyz[j]=0.;  

  for (i=0; i<nbodies; i++) {
      body3d(ephem, t, bodies[i], &GM ,&xb, &yb, &zb, vxyz); //NULL);
      x = xx[0] - xb;
      y = xx[1] - yb;
      z = xx[2] - zb;
      r2 = x*x+y*y+z*z;
      //acc = mass[i]*pow( r2 , -1.5);
      acc = GM*pow( r2 , -1.5);      
      a[0] += x*acc;
      a[1] += y*acc;
      a[2] += z*acc;
  }

  for (j=0; j<3; j++) a[j]*= -GM0;
  
  return;
}

/* Give the KBO's 3-d position, along with derivatives. */
/* Here is a leapfrog-integrator version.*/
void    
kbo3d(struct assist_ephem* ephem,
      PBASIS *pin,
      double t,
      double *xout,
      double *vout,
      double dx[],
      double dy[],
      double dz[]){

  int   i, restart=1;
  static double x[3],v[3],a[3], tx, tv, z0;
  static PBASIS psave;
  static int init=0;
  static int tdir;
  double t1,dt,dtv;

  /* decide whether we need to reset integrator to t=0*/
  if (!init) {
    init = 1;
    restart = 1;
  } else if (pin->a==psave.a && pin->adot==psave.adot &&
	     pin->b==psave.b && pin->bdot==psave.bdot &&
	     pin->g==psave.g && pin->gdot==psave.gdot) {
    restart = 0;
  } else {
    /* mismatch in parameters */
    restart=1;
  }

  /* also restart if we'd need to reverse the integration */
  if (restart==0) {
      if (tdir>0 && t<tx-0.6*TSTEP){
	  restart=1;
	  //printf("Reverse 0\n");
	  //fflush(stdout);
      }
      if (tdir<0 && t>tx+0.6*TSTEP){
	  restart=1;
	  //printf("Reverse 1\n");
	  //fflush(stdout);
      }

    if(restart==1){
    }
  }

  if (restart) {
      //printf("Reset\n");
      //fflush(stdout);
      
    /* initial time & velocity */
    z0 = 1./ pin->g;
    x[0] = pin->a *z0;
    x[1] = pin->b *z0;
    x[2] = z0;
    v[0] = pin->adot*z0;
    v[1] = pin->bdot*z0;
    v[2] = pin->gdot*z0;
    tx = tv = 0.;
    psave.a = pin->a;
    psave.adot = pin->adot;
    psave.b = pin->b;
    psave.bdot = pin->bdot;
    psave.g = pin->g;
    psave.gdot = pin->gdot;

    /* choose direction of integrator */
    if (t>=0.) tdir=+1;
    else tdir=-1;
  }

  /* leap-frog until tx is within TSTEP/2 of target time*/
  for ( ; tdir*(t-tx)>TSTEP/2.; tx = t1) {
    t1 = tx + tdir*TSTEP;
    dt = t1 - tx;
    
    /* jump velocity to middle of time step */
    accel_nbody(ephem, x, a, tx);
    /*accel(x,a,tx);*/
    dtv = 0.5*dt + (tx-tv);
    for (i=0; i<3; i++) v[i] += dtv*a[i];
    tv += dtv;

    /* Jump position using this velocity */
    for (i=0; i<3; i++) x[i] += dt*v[i];
  }

  /* Now take the last step from tx to t */
  accel_nbody(ephem, x,a,tx);
  /*accel(x,a,tx);*/
  for (i=0; i<3; i++) {
    xout[i] = x[i] + v[i]*(t-tx) + a[i]*(t-tx)*(0.5*(t+tx)-tv);
    vout[i] = v[i] + a[i]*(t-tv);
  }
  /* x and y derivatives w.r.t. parameters  - inertial approx ok here*/
  if (dx== NULL || dy==NULL || dx==NULL) return;

  for (i=1; i<=6; i++) dx[i]=dy[i]=dz[i]=0.;
  dx[1] = dy[3] = z0;
  dx[2] = dy[4] = t * z0;
  dx[5] = xout[0] * (-z0);
  dy[5] = xout[1] * (-z0);
  dz[5] = xout[2] * (-z0);
  dz[6] = t*z0;

  for (i=1; i<=6; i++) {
    dx[i] *= 1.0;
    dy[i] *= 1.0;
    dz[i] *= 1.0;
  }

  //printf("Here %lf %lf\n", jd0+t, t);
  //printf("%le %le %le\n", xout[0], xout[1], xout[2]);
  //printf("%le %le %le\n", vout[0], vout[1], vout[2]);  

  return;
}

/* Give the KBO's 3-d position, along with derivatives. */
/* Here is an ASSIST version.*/
void    
kbo3d_assist(struct assist_ephem* ephem,
	     PBASIS *pin,
	     double t,
	     double *xout,
	     double *vout,
	     double dx[],
	     double dy[],
	     double dz[]){

    int   i, restart=1;
    static double z0; //x[3],v[3],a[3], z0;
    static PBASIS psave;
    static int init=0;
    //static int tdir;
    //double t1,dt,dtv;
    XVBASIS xv;

    static struct reb_simulation* sim = NULL;
    static struct assist_extras* ax = NULL;

    /* decide whether we need to reset integrator to t=0*/
    if (!init) {
	init = 1;
	restart = 1;
    } else if (pin->a==psave.a && pin->adot==psave.adot &&
	       pin->b==psave.b && pin->bdot==psave.bdot &&
	       pin->g==psave.g && pin->gdot==psave.gdot) {
	restart = 0;
    } else {
	/* mismatch in parameters */
	restart=1;
    }

    if (restart) {
	// Initialize position and velocity
	// Restarting the integration means
	// 1. Resetting the simulation
	// 2. Converting the parameters to cartesian,
	// 3. Rotating and translating to barycentric
	//    equatorial coordinates (native to ASSIST)
	// 4. Repopulating the fields in the simulation

	// Save the pbasis state
	//printf("Reset\n");
	//fflush(stdout);
	z0 = 1./ pin->g;	
	psave.a = pin->a;
	psave.adot = pin->adot;
	psave.b = pin->b;
	psave.bdot = pin->bdot;
	psave.g = pin->g;
	psave.gdot = pin->gdot;

	if(ax != NULL) assist_free(ax);
	if(sim != NULL) reb_free_simulation(sim);

	// Create a REBOUND simulation
	sim = reb_create_simulation();
	sim->t = jd0 - ephem->jd_ref;

	ax = assist_attach(sim, ephem);	

	// Convert to barycentric equatorial coordinates
	pbasis_to_bary_eq(pin, &xv, NULL);

	// Add the particle to the REBOUND simulation.
	// These are the barycentric equatorial coordinates.
	reb_add_fmt(sim, "x y z vx vy vz",
		    xv.x, xv.y, xv.z, xv.xdot, xv.ydot, xv.zdot);

    }
    
    // Carry out the integration
    assist_integrate_or_interpolate(ax, t+jd0-ephem->jd_ref);

    // Convert back to projection coordinates
    double x1 = sim->particles[0].x;
    double y1 = sim->particles[0].y;
    double z1 = sim->particles[0].z;
    double vx1 = sim->particles[0].vx;
    double vy1 = sim->particles[0].vy;
    double vz1 = sim->particles[0].vz;

    //printf("here %lf %lf\n", jd0+t, jd0);
    //printf("%le %le %le\n", x1, y1, z1);
    //printf("%le %le %le\n", vx1, vy1, vz1);

    double xtmp, ytmp, ztmp;
    xyz_eq_to_ec(x1, y1, z1, &xtmp, &ytmp, &ztmp,NULL);
    xyz_ec_to_proj(xtmp, ytmp, ztmp, &x1, &y1, &z1, lat0, lon0, NULL);

    double vxtmp, vytmp, vztmp;
    xyz_eq_to_ec(vx1, vy1, vz1, &vxtmp, &vytmp, &vztmp,NULL);
    xyz_ec_to_proj(vxtmp, vytmp, vztmp, &vx1, &vy1, &vz1, lat0, lon0, NULL);

    /* Translate to our origin */
    xout[0] = x1 + xBary;
    xout[1] = y1 + yBary;
    xout[2] = z1 + zBary;    

    vout[0] = vx1;
    vout[1] = vy1;
    vout[2] = vz1;

    /* x and y derivatives w.r.t. parameters  - inertial approx ok here*/
    if (dx== NULL || dy==NULL || dx==NULL) return;

    // dx/dalpha = dx/dx0 * dx0/dalpha
    // dx/dadot  = dx/dvx0 * dvx0/dadot
    // dx/dbeta  = dx/dy0 * dy0/dbeta
    // dx/dbdot  = dx/dvy0 * dvy0/dbdot
    // dx/dgamma = dx/dx0*dx0/dgama + dx/dy0*dy0/dgamma + dx/dz0*dz0/dgamma + dx/dvx0*dvx0/dgamma + dx/dvy0*dvy0/dgamma + dx/dvz0*dvz0/dgamma
    // dx/dgdot  = dx/dvz0 * dvz0/dgdot

    // dx0/dalpha = z0
    // dx0/dadot = z0
    // dy0/dbeta = z0
    // dy0/dbdot = z0
    // dvz0/dgdot = z0

    
    for (i=1; i<=6; i++) dx[i]=dy[i]=dz[i]=0.;
    dx[1] = dy[3] = z0;
    dx[2] = dy[4] = t * z0;
    dx[5] = xout[0] * (-z0);
    dy[5] = xout[1] * (-z0);
    dz[5] = xout[2] * (-z0);
    dz[6] = t*z0;

    //printf("here %lf %lf\n", jd0+t, t);
    //printf("%le %le %le\n", xout[0], xout[1], xout[2]);
    //printf("%le %le %le\n", vout[0], vout[1], vout[2]);  

    return;
}


/* Function to return xyz coords of a JPL ephemeris body in
 * standard coordinate system.
 */
void
body3d(struct assist_ephem* ephem,
       double t,	/* time is in years here */
       int body,
       double *m,
       double *x, double *y, double *z,
       double *vxyz)
{
  double xxx[3];
  double xec, yec, zec;

  //static int trans[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 22, 26};

  /* get observatory posn wrt barycenter */
  struct reb_particle reb_p = {0};
  reb_p = assist_get_particle(ephem, body, t+jd0-ephem->jd_ref);

  xxx[0] = reb_p.x;
  xxx[1] = reb_p.y;
  xxx[2] = reb_p.z;  
  vxyz[0] = reb_p.vx;
  vxyz[1] = reb_p.vy;
  vxyz[2] = reb_p.vz;

  *m = reb_p.m/GM0;

  /* convert to tangent-point coord system */
  /* via ecliptic */
  xyz_eq_to_ec(xxx[0], xxx[1], xxx[2], &xec, &yec, &zec,NULL);
  xyz_ec_to_proj(xec, yec, zec, x, y, z, lat0, lon0, NULL);

  /* Translate to our origin */
  *x += xBary;
  *y += yBary;
  *z += zBary;

  /* Rotate velocity if it has been requested: */
  if (vxyz != NULL) {
    xyz_eq_to_ec(vxyz[0], vxyz[1], vxyz[2], &xec, &yec, &zec,NULL);
    xyz_ec_to_proj(xec, yec, zec, vxyz, vxyz+1, vxyz+2, lat0, lon0, NULL);
  }
    
  return;
}

/* Parse command-line looking for the standard options.
 * iarg is first arg to examine, and is returned as the
 * first non-option argument.
 * Returns 1 on parse error.
 * currently understands m, j, o, and v options.
 */
#define BUFFSIZE 512
int
read_options(int* iarg, 
	     int argc,
	     char *argv[])
{
  double d;
  for ( ; *iarg<argc; (*iarg)++) {
    if (argv[*iarg][0]!='-')
      return(0);	/*quit, now at non-option argument.*/

    if (strcasecmp(argv[*iarg]+1,"m")==0) {
      /* MPC error size */
      d = atof(argv[++(*iarg)]);
      if (d<=0.) {
	fprintf(stderr,"Bad MPC error spec %s\n",argv[*iarg]);
	return(1);
      }

      //} else if (strcasecmp(argv[*iarg]+1,"j")==0) {
      ///* ephemeris file name*/
      //set_ephem_file(argv[++(*iarg)]);

    } else if (strcasecmp(argv[*iarg]+1,"v")==0) {
      /* version number request - quit after */
      printf("Modified version of Gary Bernstein's orbfit code.\n");
      printf("Contact Matt Holman on use and \n");
      printf("  problems.  Credit Bernstein and Khushalani.\n");
      exit(1);

    } else {
      fprintf(stderr,"Unknown option %s\n",argv[*iarg]);
      return(1);
    }
  }
  return(0);
}


#include <sys/time.h>

void
fiddle_observation(OBSERVATION *obs, OBSERVATION *obs_noise)
{
  static long idum=-1;
  float gasdev(long *idum);

  /* seed the random number generator*/
  if (idum<0) {
    struct timeval tp;
    gettimeofday(&tp,NULL);
    idum = -tp.tv_usec;
  }

  /* Add measurement noise to the positions */
  obs_noise->thetax += obs->dthetax*gasdev(&idum);
  obs_noise->thetay += obs->dthetay*gasdev(&idum);

  return;
}

void intro(int argc, char *argv[]){

    /* echo the command line to output */
    printf("#");
    for (int i=0; i<argc; i++){
	printf(" %s",argv[i]);
    }
    {
#include <time.h>
	time_t timettt;
	time(&timettt);
	/* note that ctime returns string with newline at end */
	printf("\n#---%s",ctime(&timettt));
    }
}

