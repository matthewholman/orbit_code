/* Altered version of Numerical Recipes "mrqcof" to deal with our situation
 * in which we are simultaneously fitting the x and y data.  Also have
 * made use of the existing "OBSERVATION" structure as input data.
 * 4/14/99 gmb
 * 6/13/00 gmb:  add optional extra "observation" that pushes solutions
 * to be bound & nearly circular when close to degenerate.
 * see CVS logs for further alteration comments.
 */
/* energy_wt gives the relative weighting to the binding constraint.  
 * at energy_wt=1., unbound orbit is considered 2-sigma deviation from
 * the circular ideal (as is an aphelic plunging orbit.)
 * energy_wt=0 ignores this, value >1 is stronger incentive to circular.
 */

#define NRANSI
#include "nrutil.h"
#include "orbfit.h"

void mrqcof_orbit(struct reb_simulation* sim,
		  struct assist_ephem* ephem,
		  OBSERVATION obsarray[],
		  int ndata, double a[], int ia[],
		  int ma, double **alpha, double beta[], double *chisq,
		  double energy_wt)
{
	int i,j,k,l,m,mfit=0;
	double xmod,ymod,zmod,wtx,wty,sig2x,sig2y,dx,dy,*dxda,*dyda,*dzda;
	double *dxda_s,*dyda_s;
	PBASIS	params;
	OBSERVATION *oo;
	static int run=0;

	dxda=dvector(1,ma);
	dyda=dvector(1,ma);
	dzda=dvector(1,ma);

	dxda_s=dvector(1,ma);
	dyda_s=dvector(1,ma);
	
	for (j=1;j<=ma;j++)
		if (ia[j]) mfit++;
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=j;k++) alpha[j][k]=0.0;
		beta[j]=0.0;
	}
	*chisq=0.0;

	params.a = a[1];
	params.adot = a[2];
	params.b = a[3];
	params.bdot = a[4];
	params.g = a[5];
	params.gdot = a[6];

	//printf("resid %4d\n", run);
	run++;
	for (i=1;i<=ndata;i++) {
	    oo = &obsarray[i-1];
	    kbo_sphere(ephem, &params,oo,&xmod,dxda,&ymod,dyda,&zmod,dzda);

	    sig2x=1.0/(oo->dthetax*oo->dthetax);
	    sig2y=1.0/(oo->dthetay*oo->dthetay);
	    // covariance cross term would be here.

	    // This local tangent plane trick is from Danby.
	    // Could encapsulate this in a function.
	    double Ax, Ay, Az, A;
	    double Dx, Dy, Dz, D;

	    // The A and D vectors could be calculated once, rather
	    // than repeatedly.
		  
	    Ax =  oo->z;
	    Ay =    0.0;
	    Az = -oo->x;

	    A = sqrt(Ax*Ax + Ay*Ay + Az*Az);
	    Ax /= A; Ay /= A; Az /= A;

	    Dx = -oo->x*oo->y;
	    Dy = oo->x*oo->x + oo->z*oo->z;
	    Dz = -oo->z*oo->y;
		  
	    D = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
	    Dx /= D; Dy /= D; Dz /= D;

	    dx = -(xmod*Ax + ymod*Ay + zmod*Az);
	    dy = -(xmod*Dx + ymod*Dy + zmod*Dz);
	    /*
	    printf("obs %6d %7.2lf %7.2lf ", i, dx*206265, dy*206265);
	    double thresh = 10;
	    if(fabs(dx*206265)>thresh || fabs(dy*206265)>thresh){
		printf(" *\n");
	    }else{
		printf("  \n");
	    }
	    */
	    

	    // Find the partial derivatives along the two sky plane
	    // vectors.
	    for(l=1; l<=ma; l++){
		if (ia[l]) {
		    dxda_s[l] = dxda[l]*Ax + dyda[l]*Ay + dzda[l]*Az;
		    dyda_s[l] = dxda[l]*Dx + dyda[l]*Dy + dzda[l]*Dz;
		}
	    }

	    for (j=0,l=1;l<=ma;l++) {
		if (ia[l]) {
		    wtx=dxda_s[l]*sig2x;
		    wty=dyda_s[l]*sig2y;
		    for (j++,k=0,m=1;m<=l;m++){
			// observation covariance would go here.
			if (ia[m]) alpha[j][++k] += wtx*dxda_s[m]+wty*dyda_s[m];
		    }
		    beta[j] += dx*wtx+dy*wty;
		}
	    }
	    *chisq += dx*dx*sig2x+dy*dy*sig2y;

	}

	/*** Add here the energy terms:  the "data" in this context is
	 * the fractional diff fb of tranverse KE from -PE/2. */
	/*
	double fb1;
	if (energy_wt>0.) {
	  double g3, g4;
	  g3 = pow(params.g, -3.)/GM;
	  g4 = g3/params.g;
	  fb1 = (params.adot*params.adot 
		 + params.bdot*params.bdot
		 + params.gdot*params.gdot) *  g3
	        -1.;
	  sig2x = 4.*energy_wt;	  //Make fb1=1 be 2-sigma
	  dxda[1] = dxda[3] = 0.;
	  dxda[2] = 2*params.adot * g3;
	  dxda[4] = 2*params.bdot * g3;
	  dxda[6] = 2*params.gdot * g3;
	  dxda[5] = -3.*(params.adot*params.adot 
			 + params.bdot*params.bdot
			 + params.gdot*params.gdot) * g4;
	  for (j=0,l=1;l<=ma;l++) {
		  if (ia[l]) {
		    wtx=dxda[l]*sig2x;
		    for (j++,k=0,m=1;m<=l;m++)
		      if (ia[m]) alpha[j][++k] += wtx*dxda[m];
		    beta[j] += -fb1*wtx;
		  }
	  }
	  *chisq += fb1*fb1*sig2x;

	  // Add 2nd-derivative terms for the energy constraint 
	  sig2x *= fb1;
	  alpha[2][2] += sig2x * 2.* g3;
	  alpha[4][4] += sig2x * 2.* g3;
	  alpha[6][6] += sig2x * 2.* g3;
	  alpha[2][5] -= 6. * params.adot * g4;
	  alpha[5][2] -= 6. * params.adot * g4;
	  alpha[4][5] -= 6. * params.bdot * g4;
	  alpha[5][4] -= 6. * params.bdot * g4;
	  alpha[6][5] -= 6. * params.gdot * g4;
	  alpha[5][6] -= 6. * params.gdot * g4;
	  alpha[5][5] += 12.*(fb1+1.)/(params.g*params.g);
	}
	*/
	
	for (j=2;j<=mfit;j++)
		for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];
	free_dvector(dxda,1,ma);
	free_dvector(dyda,1,ma);
	free_dvector(dzda,1,ma);	
	free_dvector(dxda_s,1,ma);
	free_dvector(dyda_s,1,ma);

}
#undef NRANSI
