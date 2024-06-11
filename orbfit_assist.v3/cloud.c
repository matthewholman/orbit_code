
#include "orbfit.h"
#include "rebound.h"
#include "assist.h"

char   *help[] = {
    " fit_radec:  Fit KBO orbit to observed astrometric RA/DEC",
    " usage:  fit_radec [-m mpc_error] [-j JPL_file] [-v]",
    "  mpc_error  is arcsec uncertainty to apply to MPC-format",
    "             observations.  Default is 0.2",
    "  JPL_file   is binary ephemeris file.  Default is binEphem.432, or",
    "             a file specified by environment variable ORBIT_EPHEMERIS",
    "  stdin      contains one observation per line, with format",
    "             JD  RA Dec error obscode",
    "  stdout     is a file containing best-fit alpha, beta, etc.",
    "             and its covariance matrix.",
    "  Residuals are dumped to stderr.",
    0
};

void print_help(void)
{
    for (int i = 0; help[i] != 0; i++)
	fprintf (stderr, "%s\n", help[i]);
    exit (1);
}


 

int main(int argc, char *argv[])
{

    
    void intro(int argc, char *argv[]);
    PBASIS p;
    ORBIT  orbit;
    XVBASIS xv, xv_eq;

    double chisq;

    double **covar;
    int dof;

    int i;

// Could make this allocation dynamic.
    
    OBSERVATION obsarray[MAXOBS];
    OBSERVATION obsarray_noise[MAXOBS];        
    int nobs;

    covar = dmatrix(1,6,1,6);

    {
	int iarg=1;
	if (argc>1 && *argv[1]=='^') print_help();
	if (read_options(&iarg, argc, argv)) print_help();
    }

    if (read_observations_simple(obsarray,NULL,&nobs)) {
	fprintf(stderr, "Error reading input observations\n");
	exit(1);
    }

    for(int i=0; i<nobs; i++){
	obsarray_noise[i] = obsarray[i];
    }

    //intro(argc, argv);

    // Create a REBOUND simulation
    struct reb_simulation* sim = reb_create_simulation();
    reb_reset_temporary_pointers(sim);    

    // Load the ephemeris data
    struct assist_ephem* ephem = assist_ephem_create(
            "/Users/mholman/assist/data/linux_p1550p2650.440",
            "/Users/mholman/assist/data/sb441-n16.bsp");
    if (!ephem){
        printf("Cannot create ephemeris structure.\n");
        exit(-1);
    }
    
    // One has the option to change other setting.
    ephem->jd_ref = 2451545.0; // Reference JD. This line can be commented out as this is the default.
    
    // Attach the ASSIST framework
    // This will set the additional force routine in the REBOUND simulation,
    // the IAS15 integrator, the gravity module, etc.
    struct assist_extras* ax = assist_attach(sim, ephem);

    printf("# desig epoch            a(AU)         e            incl(deg)    longnode(deg)  argperi(deg)   meananom(deg)\n");
    int ntrials = 10001;
    for(int j=0; j<ntrials; j++){
	if (j !=0){
	    for(int i=0; i<nobs; i++){
		fiddle_observation(obsarray+i, obsarray_noise+i);
	    }
	}

	/* Call subroutine to do the actual fitting: */
	fit_observations(sim, ephem, obsarray_noise, nobs, &p, covar, &chisq, &dof,stdout);
	pbasis_to_bary(&p, &xv, NULL);

	// Initialize particle for Sun position
	struct reb_particle reb_p = {0};
	reb_p = assist_get_particle(ephem, 0, jd0-ephem->jd_ref);

	double GMsun = reb_p.m;

	double x_ec, y_ec, z_ec;
	double vx_ec, vy_ec, vz_ec;    
    
	xyz_eq_to_ec(reb_p.x, reb_p.y, reb_p.z, &x_ec, &y_ec, &z_ec, NULL);
	xyz_eq_to_ec(reb_p.vx, reb_p.vy, reb_p.vz, &vx_ec, &vy_ec, &vz_ec, NULL);    

	xv.x -= x_ec;
	xv.y -= y_ec;
	xv.z -= z_ec;
	xv.xdot -= vx_ec;
	xv.ydot -= vy_ec;
	xv.zdot -= vz_ec;

	State state;
	state.x = xv.x;
	state.y = xv.y;
	state.z = xv.z;
	state.xd = xv.xdot;
	state.yd = xv.ydot;
	state.zd = xv.zdot;

	double a, e, incl, longnode, argperi, meananom;
	keplerian(GMsun, state, &a, &e, &incl, &longnode, &argperi, &meananom);

	/*
	printf("# ----------- heliocentric elements -------------\n");
	printf("# desig=      %s\n", obsarray[0].desig);
	printf("# epoch=      %.10lf\n", jd0);
	printf("# a=          %.10lf\n", a);
	printf("# e=          %.10lf\n", e);
	printf("# incl=       %.10lf\n", incl*180./PI);
	printf("# longnode=   %.10lf\n", principal_value(longnode)*180./PI);    
	printf("# argperi =   %.10lf\n", principal_value(argperi)*180./PI);
	printf("# meananom=   %.10lf\n", principal_value(meananom)*180./PI);
	*/
	printf("%s %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
	       obsarray[0].desig, jd0,
	       a, e, incl*180./PI,
	       principal_value(longnode)*180./PI,
	       principal_value(argperi)*180./PI,
	       principal_value(meananom)*180./PI);
	       
	//printf("# nobs=       %d\n", nobs);

    }
    
    free_dmatrix(covar,1,6,1,6);

    assist_free(ax);
    assist_ephem_free(ephem);
    reb_free_simulation(sim);

    exit(0);
} 



  
