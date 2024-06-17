
#include <string.h>
#include <stdio.h>
#include <math.h>

#define KEPLER_FLAG 1

typedef struct {
  double x, y, z, xd, yd, zd;
} State;

extern double machine_epsilon;

#define NEWTON_MAX 10
#define LAGCON_MAX 15
#define PI 3.14159265358979323846
#define EPS 1e-14

/* kepler step in universal variables from Danby (1988) p. 178 */

// Advances an initial state s0 by a time dt.
// Stores the new state in s.
// Stores the f, g, fdot, and gdot values for external use.
int universal_fg(double gm, double dt, State *s0, State *s, double *f, double *g, double *fdot, double *gdot)
{
  double r0, v0s, u, alpha;
  double zeta;
  double r;
  double _f, _fdot, _g, _gdot;
  double ss;
  double g0, g1, g2, g3, g4, g5;
  int flag;
  int kepler(double gm, double dt, double r0,  double alpha, double u, double zeta,
	     double *rx, double *s, double *g0, double *g1, double *g2, double *g3, double *g4, double *g5);
  
  double simple_guess(double dt, double r0,  double u);
  double initial_guess_v2(double gm, double dt, double r0,  double alpha, double u, double zeta);

  flag = 0;
  r0 = sqrt(s0->x*s0->x + s0->y*s0->y + s0->z*s0->z);
  v0s = s0->xd*s0->xd + s0->yd*s0->yd + s0->zd*s0->zd;
  u = s0->x*s0->xd + s0->y*s0->yd + s0->z*s0->zd;
  alpha = 2.0*gm/r0 - v0s;
  zeta = gm - alpha*r0;

  // Get an initial guess for the solution to
  // kepler's equation.
  //ss = simple_guess(dt, r0,  u);
  ss = initial_guess_v2(gm, dt, r0,  alpha, u, zeta);  

  /// solve kepler's equation
  flag = kepler(gm, dt, r0, alpha, u, zeta, &r, &ss, &g0, &g1, &g2, &g3, &g4, &g5);

  _f = 1.0 - (gm/r0)*g2;
  _g = dt - gm*g3;
  _fdot = - (gm/(r*r0))*g1;
  _gdot = 1.0 - (gm/r)*g2;

  *f = _f;
  *g = _g;
  *fdot = _fdot;
  *gdot = _gdot;

  s->x = _f*s0->x + _g*s0->xd;
  s->y = _f*s0->y + _g*s0->yd;
  s->z = _f*s0->z + _g*s0->zd;
  s->xd = _fdot*s0->x + _gdot*s0->xd;
  s->yd = _fdot*s0->y + _gdot*s0->yd;
  s->zd = _fdot*s0->z + _gdot*s0->zd;
  
  return(flag);
}

// Advance an initial cartesian state s0 to an array dt of time offset.
// dt is of length num
// Stores the new states in an array s.
// Stores the f, g, fdot, and gdot values in arrays for external use.
void universal_fg_n(int num, double gm, double *dt, State *s0, State *s, double *f, double *g, double *fdot, double *gdot, int *flags)
{
  double r0, v0s, u, alpha;
  double zeta;
  double r;
  double _f, _fdot, _g, _gdot;
  double ss;
  double g0, g1, g2, g3, g4, g5;
  int kepler(double gm, double dt, double r0,  double alpha, double u, double zeta,
	     double *rx, double *s, double *g0, double *g1, double *g2, double *g3, double *g4, double *g5);
  double simple_guess(double dt, double r0,  double u);
  double initial_guess_v2(double gm, double dt, double r0,  double alpha, double u, double zeta);  

  r0 = sqrt(s0->x*s0->x + s0->y*s0->y + s0->z*s0->z);
  v0s = s0->xd*s0->xd + s0->yd*s0->yd + s0->zd*s0->zd;
  u = s0->x*s0->xd + s0->y*s0->yd + s0->z*s0->zd;
  alpha = 2.0*gm/r0 - v0s;
  zeta = gm - alpha*r0;

  //check that n>0

  //ss = simple_guess(dt[0], r0,  u);
  ss = initial_guess_v2(gm, dt[0], r0,  alpha, u, zeta);    

  for(int i=0; i<num; i++){
      /* solve kepler's equation */
      flags[i] = kepler(gm, dt[i], r0, alpha, u, zeta, &r, &ss, &g0, &g1, &g2, &g3, &g4, &g5);

      _f = 1.0 - (gm/r0)*g2;
      _g = dt[i] - gm*g3;
      _fdot = - (gm/(r*r0))*g1;
      _gdot = 1.0 - (gm/r)*g2;

      f[i] = _f;
      g[i] = _g;
      fdot[i] = _fdot;
      gdot[i] = _gdot;

      s[i].x = _f*s0->x + _g*s0->xd;
      s[i].y = _f*s0->y + _g*s0->yd;
      s[i].z = _f*s0->z + _g*s0->zd;
      s[i].xd = _fdot*s0->x + _gdot*s0->xd;
      s[i].yd = _fdot*s0->y + _gdot*s0->yd;
      s[i].zd = _fdot*s0->z + _gdot*s0->zd;

  }
  
  return;
  
}

// Advance an initial cartesian state s0 to an array dt of time offset.
// dt is of length num
// Stores the new states in an array s.
void universals(int num, double gm, double *dt, State *s0, State *s, int *flags)
{
  double r0, v0s, u, alpha;
  double zeta;
  double r;
  double _f, _fdot, _g, _gdot;
  double ss;
  double g0, g1, g2, g3, g4, g5;
  int kepler(double gm, double dt, double r0,  double alpha, double u, double zeta,
	     double *rx, double *s, double *g0, double *g1, double *g2, double *g3, double *g4, double *g5);

  double simple_guess(double dt, double r0,  double u);  

  r0 = sqrt(s0->x*s0->x + s0->y*s0->y + s0->z*s0->z);
  v0s = s0->xd*s0->xd + s0->yd*s0->yd + s0->zd*s0->zd;
  u = s0->x*s0->xd + s0->y*s0->yd + s0->z*s0->zd;
  alpha = 2.0*gm/r0 - v0s;
  zeta = gm - alpha*r0;

  //check that n>0

  ss = simple_guess(dt[0], r0, u);    
  
  for(int i=0; i<num; i++){
      /* solve kepler's equation */      
      flags[i] = kepler(gm, dt[i], r0, alpha, u, zeta, &r, &ss, &g0, &g1, &g2, &g3, &g4, &g5);

      _f = 1.0 - (gm/r0)*g2;
      _g = dt[i] - gm*g3;
      _fdot = - (gm/(r*r0))*g1;
      _gdot = 1.0 - (gm/r)*g2;

      /*
      f[i] = _f;
      g[i] = _g;
      fdot[i] = _fdot;
      gdot[i] = _gdot;
      */

      s[i].x = _f*s0->x + _g*s0->xd;
      s[i].y = _f*s0->y + _g*s0->yd;
      s[i].z = _f*s0->z + _g*s0->zd;
      s[i].xd = _fdot*s0->x + _gdot*s0->xd;
      s[i].yd = _fdot*s0->y + _gdot*s0->yd;
      s[i].zd = _fdot*s0->z + _gdot*s0->zd;

  }
  
  return;
  
}


int universal_kepler(double gm, double dt, State *s0, State *s)
{
  double r0, v0s, u, alpha;
  double zeta;
  double r;
  double f, fdot, g, gdot;
  double ss;
  double g0, g1, g2, g3, g4, g5;
  int flag;
  int kepler(double gm, double dt, double r0,  double alpha, double u, double zeta,
	     double *rx, double *s, double *g0, double *g1, double *g2, double *g3, double *g4, double *g5);
  
  double simple_guess(double dt, double r0,  double u);

  flag = 0;
  r0 = sqrt(s0->x*s0->x + s0->y*s0->y + s0->z*s0->z);
  v0s = s0->xd*s0->xd + s0->yd*s0->yd + s0->zd*s0->zd;
  u = s0->x*s0->xd + s0->y*s0->yd + s0->z*s0->zd;
  alpha = 2.0*gm/r0 - v0s;
  zeta = gm - alpha*r0;

  ss = simple_guess(dt, r0,  u);      

  /* solve kepler's equation */
  flag = kepler(gm, dt, r0, alpha, u, zeta, &r, &ss, &g0, &g1, &g2, &g3, &g4, &g5);

  f = 1.0 - (gm/r0)*g2;
  g = dt - gm*g3;
  fdot = - (gm/(r*r0))*g1;
  gdot = 1.0 - (gm/r)*g2;
  //gdot = (1.0 + g*fdot)/f;

  s->x = f*s0->x + g*s0->xd;
  s->y = f*s0->y + g*s0->yd;
  s->z = f*s0->z + g*s0->zd;
  s->xd = fdot*s0->x + gdot*s0->xd;
  s->yd = fdot*s0->y + gdot*s0->yd;
  s->zd = fdot*s0->z + gdot*s0->zd;

  return(flag);
}

int kepler(double gm, double dt, double r0,  double alpha, double u, double zeta,
	   double *rx, double *s, double *g0, double *g1, double *g2, double *g3, double *g4, double *g5)
{

  double x;
  double c0, c1, c2, c3, c4, c5;
  double f, fp, fpp, fppp;
  double ss, sst, dss;
  double sgn(double x);
  void cfun(double z, double *c0, double *c1, double *c2, double *c3, double *c4, double *c5);
  double ln;
  int nc;

  // This will be the previous value of ss if more than one value of
  // dt is being examined.
  //ss = simple_guess(gm, dt, r0,  alpha, u, zeta);
  ss = *s;

  sst = ss;
  nc = 0;

  /* solve kepler's equation */

  /* Newton Iteration */
  do{
    x = ss*ss*alpha;
    cfun(x, &c0, &c1, &c2, &c3, &c4, &c5);
    c1 *= ss; c2 *= ss*ss; c3 *= ss*ss*ss;
    
    f    = r0*c1 + u*c2 + gm*c3 - dt;
    fp   = r0*c0 + u*c1 + gm*c2;
    fpp  = zeta*c1 + u*c0;
    fppp = zeta*c0 - u*alpha*c1;
    dss = -f/fp;
    dss = -f/(fp + dss*fpp/2.0);
    dss = -f/(fp + dss*fpp/2.0 + dss*dss*fppp/6.0);
    ss += dss;
    nc++;

    }while(fabs(dss) > EPS && nc<NEWTON_MAX);


  /* Laguerre-Conway iteration if Newton fails */
  if(fabs(dss) >= EPS && nc == NEWTON_MAX){
      ss = sst;
      ln = 5.0;
      nc = 0;
      do{
	  x = ss*ss*alpha;
	  cfun(x, &c0, &c1, &c2, &c3, &c4, &c5);
	  c1 *= ss; c2 *= ss*ss; c3 *= ss*ss*ss;
    
	  f   = r0*c1 + u*c2 + gm*c3 - dt;
	  fp  = r0*c0 + u*c1 + gm*c2;
	  fpp = zeta*c1 + u*c0;
	  dss = -ln*f/(fp + sgn(fp)*
		       sqrt(fabs((ln-1.0)*(ln-1.0)*fp*fp - (ln-1.0)*ln*f*fpp)));
	  ss += dss;
	  nc++;
      
      }while(fabs(dss) > EPS && nc<LAGCON_MAX);

      // could replace this with a bisection routine if Newton and Laguerre-Conway both fail.
      if(fabs(dss) >= EPS && nc == LAGCON_MAX){
	  return(KEPLER_FLAG);
      }
      //printf("laguerre nc: %d %.16le\n", nc, dss);                
  }

  x = ss*ss*alpha;
  cfun(x, &c0, &c1, &c2, &c3, &c4, &c5);

  c1 *= ss; c2 *= ss*ss; c3 *= ss*ss*ss; c4 *= ss*ss*ss*ss; c5 *= ss*ss*ss*ss*ss;
  *g0 = c0; *g1 = c1; *g2 = c2; *g3 = c3; *g4 = c4; *g5 = c5;
  *rx = fp;
  *s = ss;
  return(0);
}

void cfun(double z, double *c0, double *c1, double *c2, double *c3, double *c4, double *c5)
{

/* Stumpf c-functions by argument four-folding, Mikkola */

  double C6=1.0/6.0,C132=1.0/132.0,C56=1.0/56.0;
  double C30=1.0/30.0,C24=1.0/24.0,C156=1.0/156.0;
  double C90=1.0/90.0,C110=1.0/110.0,C16=1.0/16.0,C8=1.0/8.0;
  double C72=1.0/72.0,C42=1.0/42.0,C120=1.0/120.0,U=1.0;

  double h;
  int i, k;
  h = z;
  k = 0;
  while(fabs(h)>=0.1){
    h=0.25*h;
    k++;
  }
  
  *c4 = (U-h*(U-h*(U-h*C90/(U+h*C132))*C56)*C30)*C24;
  *c5 = (U-h*(U-h*(U-h*C110/(U+h*C156))*C72)*C42)*C120;

  for(i=0;i<k;i++){
    *c3 = C6 - h*(*c5);
    *c2 = 0.5 - h*(*c4);
    *c5 = (*c5 + *c4 + (*c2)*(*c3))*C16;
    *c4 = (*c3)*(2.0 - h*(*c3))*C8;
    h=4.0*h;
  }
   
  *c3 = C6 - z*(*c5);
  *c2 = 0.5 - z*(*c4);
  *c1 = 1.0 - z*(*c3);
  *c0 = 1.0 - z*(*c2);
  
  return;
}

double sgn(double x)
{
    return (x > 0) - (x < 0);
}

void stumpff(double x, double *c0, double *c1, double *c2, double *c3)
{
    
  int n;
  double xm;
  double d0, d1, d2, d3;

  n = 0; xm = 0.1;
  while(fabs(x)>xm){
    n++;
    x /= 4;
  }


  d2=(1-x*(1-x*(1-x*(1-x*(1-x*(1-x/182.0)/132.0)/90.0)/56.0)/30.0)/12.0)/2.0;
  d3=(1-x*(1-x*(1-x*(1-x*(1-x*(1-x/210.0)/156.0)/110.0)/72.0)/42.0)/20.0)/6.0;

  d1=1.0-x*d3;
  d0=1.0-x*d2;

  while(n>0){
    n--;
    d3=(d2+d0*d3)/4.0;
    d2=d1*d1/2.0;
    d1=d0*d1;
    d0=2.0*d0*d0-1.0;
  }
  
  *c0 = d0;
  *c1 = d1;
  *c2 = d2;
  *c3 = d3;
  return;
}


double initial_guess(double gm, double dt, double r0,  double alpha, double u, double zeta)
{

    double ss;
    double x, y, sigma;
    double a, en, ch, sh, e, ec, es, dm;    
    /* determine initial guesses */
    if(fabs(dt/r0) > 0.2){  // This has to be wrong because it is not unitless.
	if(alpha<=0.0){
	    /* hyperbolic motion */
	    a = gm/alpha;
	    en = sqrt(-gm/(a*a*a));
	    ch = 1.0 - r0/a;
	    sh = u/sqrt(-a*gm);
	    e = sqrt(ch*ch - sh*sh);
	    dm = en*dt;
	    if(dm<0.0){
		ss = -log((-2.0*dm + 1.8*e)/(ch - sh))/sqrt(-alpha);
	    }else{
		ss = log((2.0*dm + 1.8*e)/(ch + sh))/sqrt(-alpha);
	    }
	}else{

	    /* elliptic motion */
	    a = gm/alpha;
	    en = sqrt(gm/(a*a*a));
	    ec = 1.0 - r0/a;
	    es = u/(en*a*a);
	    e = sqrt(ec*ec + es*es);
	    dt -= ((int)(en*dt/(2.0*PI)))*2.0*PI/en;
	    y = en*dt - es;

	    /* Danby Guess */
	    sigma = sgn((es*cos(y) + ec*sin(y)));
	    x = y + sigma*0.85*e;

	    ss = x/sqrt(alpha);

	}
      
    }else{
	/* close to parabolic */
	ss = dt/r0 - (dt*dt*u)/(2.0*r0*r0*r0);
    }

    return(ss);

}

double simple_guess(double dt, double r0,  double u)
{

    double rdot = u/r0;
    double dtr0 = dt/r0;

    double ss = dtr0 - 0.5*rdot*dtr0*dt*r0;// + (0.5*rdot*rdot - zeta/(6.0*r0))*dtr0*dtr0*dtr0;

    return(ss);

}

double initial_guess_v2(double gm, double dt, double r0,  double alpha, double u, double zeta)
{

    double ss;
    double x, y, sigma;
    double a, en, ch, sh, e, ec, es, dm;    
    /* determine initial guesses */
    if(alpha<=0.0){
	/* hyperbolic motion */
	a = gm/alpha;
	en = sqrt(-gm/(a*a*a));
	ch = 1.0 - r0/a;
	sh = u/sqrt(-a*gm);
	e = sqrt(ch*ch - sh*sh);
	dm = en*dt;
	if(dm<0.0){
	    ss = -log((-2.0*dm + 1.8*e)/(ch - sh))/sqrt(-alpha);
	}else{
	    ss = log((2.0*dm + 1.8*e)/(ch + sh))/sqrt(-alpha);
	}
    }else{

	/* elliptic motion */
	a = gm/alpha;
	en = sqrt(gm/(a*a*a));
	ec = 1.0 - r0/a;
	es = u/(en*a*a);
	e = sqrt(ec*ec + es*es);
	dt -= ((int)(en*dt/(2.0*PI)))*2.0*PI/en;
	y = en*dt - es;
	
	/* Danby Guess */
	sigma = sgn((es*cos(y) + ec*sin(y)));
	x = y + sigma*0.85*e;

	ss = x/sqrt(alpha);

    }
      
    return(ss);

}
