
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

    intro(argc, argv);

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

    /* Call subroutine to do the actual fitting: */
    fit_observations(sim, ephem, obsarray, nobs, &p, covar, &chisq, &dof,stdout);
    //fit_observations_prelim(sim, ephem, obsarray+nobs-1-200, 200, &p, covar, &chisq, &dof,stdout);
    //fit_observations_prelim(sim, ephem, obsarray+nobs-1-400, 400, &p, covar, &chisq, &dof,stdout);
    //fit_observations_prelim(sim, ephem, obsarray+nobs-1-1000, 1000, &p, covar, &chisq, &dof,stdout);
    //fit_observations_prelim(sim, ephem, obsarray+nobs-1-2000, 2000, &p, covar, &chisq, &dof,stdout);
    //fit_observations_prelim(sim, ephem, obsarray+nobs-1-3000, 3000, &p, covar, &chisq, &dof,stdout);                    
    //fit_observations_prelim(sim, ephem, obsarray, nobs, &p, covar, &chisq, &dof,stdout);    

    printf("# Exact a, adot, b, bdot, g, gdot:\n");
    printf("%15.12le %15.12le %15.12le %15.12le %15.12le %15.12le\n",p.a,p.adot,p.b,
	   p.bdot, p.g, p.gdot);
    pbasis_to_bary(&p, &xv, NULL);

    orbitElements(&xv, &orbit);
    printf("# desig=%s\n", obsarray[0].desig);
    printf("# a=%lf AU e=%lf i=%lf deg W=%lf deg,w=%lf deg,Tp=%lf\n",orbit.a, orbit.e, orbit.i, orbit.lan, orbit.aop, orbit.T);
    

    {
	double d, dd;
	d = sqrt(xBary*xBary + yBary*yBary + pow(zBary-1/p.g,2.));
	dd = d*d*sqrt(covar[5][5]);
	printf("# Barycentric distance %.3f+-%.3f\n",d,dd);
    }

    
    printf("# Chi-squared of fit: %.2f DOF: %d\n",chisq,dof);
    
    /* Print the covariance matrix to stdout */
    printf("# Covariance matrix: \n");
    print_matrix(stdout,covar,6,6);

    /* Print out information on the coordinate system */
    printf("#     lat0       lon0       xBary     yBary      zBary   JD0\n");
    printf("%15.10f %15.10f %12.9f %12.9f %12.9f  %.6f\n",
	   lat0/DTOR,lon0/DTOR,xBary,yBary,zBary,jd0);



    /* Dump residuals */
    printf("# Best fit orbit gives:\n");
    printf("# obs desig     time(days)     RA(deg)      RA_mod(deg)  resid(\") Dec(deg)    Dec_mod(deg) resid(\") code\n");
    double rms = 0.0;

    for (i=0; i<nobs; i++) {
      double x,y;
      kbo2d(ephem, &p, &obsarray[i], &x, NULL, &y, NULL);

      double lat_ec, lon_ec;
      double ra_eq, dec_eq;
      proj_to_ec(x, y, &lat_ec, &lon_ec, lat0, lon0, NULL);
      ec_to_eq(lat_ec, lon_ec, &ra_eq, &dec_eq, NULL);
      
      ra_eq = principal_value(ra_eq);
      /*
      printf("# %3d %s %15.7f %12.7f %12.7f %10.4f %12.7f %12.7f %10.4f   %s\n",
	      i, obsarray[i].desig, obsarray[i].obstime,
	      obsarray[i].ra, ra_eq/DTOR, (obsarray[i].thetax-x)/ARCSEC,	      
	      obsarray[i].dec, dec_eq/DTOR, (obsarray[i].thetay-y)/ARCSEC,	      
	      obsarray[i].obscode);
      */
      printf("# %3d %s %15.7f %12.7f %12.7f %10.4f %12.7f %12.7f %10.4f %14.7e %14.7e %14.7e %14.7e  %s\n",
	     i, obsarray[i].desig, obsarray[i].obstime,
	     obsarray[i].ra, ra_eq/DTOR, (obsarray[i].thetax-x)/ARCSEC,	      
	     obsarray[i].dec, dec_eq/DTOR, (obsarray[i].thetay-y)/ARCSEC,
	     obsarray[i].thetax, x,
	     obsarray[i].thetay, y,
	     obsarray[i].obscode);
      //printf("%14.7e %14.7e %14.7e\n", obsarray[i].xe, obsarray[i].ye, obsarray[i].ze);
    }

    pbasis_to_bary(&p, &xv, NULL);

    //pbasis_to_bary_eq(&p, &xv_eq, NULL);
    //printf("#\n# bary_eq %.12lf %.15f %.15f %.15f %.15e %.15e %.15e\n", jd0, xv_eq.x, xv_eq.y, xv_eq.y, xv_eq.xdot, xv_eq.ydot, xv_eq.zdot);    

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
        
    //printf("#\n# %.12lf %.15f %.15f %.15f %.15e %.15e %.15e\n", jd0, xv.x, xv.y, xv.z, xv.xdot, xv.ydot, xv.zdot);

    printf("# ----------- heliocentric elements -------------\n");
    printf("# desig=      %s\n", obsarray[0].desig);
    printf("# epoch=      %.10lf\n", jd0);
    printf("# a=          %.10lf\n", a);
    printf("# e=          %.10lf\n", e);
    printf("# incl=       %.10lf\n", incl*180./PI);
    printf("# longnode=   %.10lf\n", principal_value(longnode)*180./PI);    
    printf("# argperi =   %.10lf\n", principal_value(argperi)*180./PI);
    printf("# meananom=   %.10lf\n", principal_value(meananom)*180./PI);
    
    
    //for(int i=0; i<20; i++){
    //reb_p = assist_get_particle(ephem, i, 10.);
    //printf("%d %le %lf %lf %lf\n", i, reb_p.m, reb_p.x, reb_p.y, reb_p.z);
    //}

    
    orbitElements(&xv, &orbit);

    //double x_eq, y_eq, z_eq;
    //double vx_eq, vy_eq, vz_eq;    
    //xyz_ec_to_eq(xv.x, xv.y, xv.z, &x_eq, &y_eq, &z_eq, NULL);
    //xyz_ec_to_eq(xv.xdot, xv.ydot, xv.zdot, &vx_eq, &vy_eq, &vz_eq, NULL);
    //printf("#\n# %.12lf %.15f %.15f %.15f %.15e %.15e %.15e\n", jd0, x_eq, y_eq, z_eq, vx_eq, vy_eq, vz_eq);
    //printf("#\n# %.12lf %.15f %.15f %.15f %.15e %.15e %.15e\n", jd0, xv.x, xv.y, xv.z, xv.xdot, xv.ydot, xv.zdot);
    
    printf("# ----------- for ele220 -------------\n");
    printf("# desig=      %s\n", obsarray[0].desig);
    double meanMotion = sqrt(GM0/(pow(orbit.a,3.)));    
    printf("# Tp=         %.8lf\n", orbit.T);
    printf("# meananom=   %.8lf\n", 180./PI*principal_value(orbit.meananom));
    printf("# argperim=   %.8lf\n", orbit.aop);
    printf("# longnode=   %.8lf\n", orbit.lan);
    printf("# incl=       %.8lf\n", orbit.i);
    printf("# q=          %.8lf\n", orbit.a*(1.0-orbit.e));
    printf("# e=          %.8lf\n", orbit.e);
    printf("# H=          %.8lf\n", 0.0);
    printf("# RMS=        %.3lf\n", sqrt(rms/nobs));

    printf("# nobs=       %d\n", nobs);
    double starttime =   obsarray[0].obstime+jd0;
    double endtime =   obsarray[nobs-1].obstime+jd0;

    double numopp = (endtime - starttime)*DAY + 1.0;
    if((numopp - floor(numopp))>0.5) numopp += 1.0;

    printf("# numopp=     %d\n", (int)numopp);    
    
    printf("# timebegin=  %lf\n", starttime);
    printf("# endtime=    %lf\n", endtime);
    printf("# epoch=      %lf\n", jd0);

    printf("# unc=        %.3lf\n", obsarray[0].dthetax/ARCSEC);

    printf("# computer=   %s\n", "modified_Bernstein");

    free_dmatrix(covar,1,6,1,6);

    assist_free(ax);
    assist_ephem_free(ephem);
    reb_free_simulation(sim);

    exit(0);
} 



  
