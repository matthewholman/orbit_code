#include "heliocentric.h"
#define PI 3.14159265358979323846

double GMsun = 1.0;

extern double machine_epsilon;

int main(int argc, char **argv)
{
    
    State *p;
    double *positions, *velocities;
    double *a, *e, *incl, *longnode, *argperi, *meananom;
    double *dt;
    Elements *elements;
    int i, n, *flags;

    void keplerian(double gm, State state, 
		   double *a, double *e, double *incl, double *longnode, double *argperi, double *meananom);

    void cartesian(double gm, 
		   double a, double e, double i, double longnode, double argperi, double meananom,
		   State *state);

    void cartesians(int num, double gm, 
		  double *a, double *e, double *incl, double *longnode, double *argperi, double *meananom, 
		  State *state);

  void cartesian_vectors(int num, double gm, 
			 double *a, double *e, double *incl, double *longnode, double *argperi, double *meananom,
			 double *positions,
			 double *velocities);

  void cartesian_elements(int num, double gm, 
			  Elements *elements,
			  double *positions,
			  double *velocities);

  void universals(int num, double gm, double *dt, State *s0, State *s, int *flags);

  int universal_fg(double gm, double dt, State *s0, State *s, double *f, double *g, double *fdot, double *gdot);
  
  void universal_fg_n(int num, double gm, double *dt, State *s0, State *s, double *f, double *g, double *fdot, double *gdot, int *flags);  
  
  if(argc != 2) {
    printf("args: n\n");
    exit(-1);
  }

  sscanf(argv[1], "%d", &n);

  p = (State*) calloc(n,sizeof(State));
  positions  = (double*) calloc(n,3*sizeof(double));
  velocities = (double*) calloc(n,3*sizeof(double));    
  a        = (double*) calloc(n,sizeof(double));
  e        = (double*) calloc(n,sizeof(double));		      
  incl     = (double*) calloc(n,sizeof(double));
  longnode = (double*) calloc(n,sizeof(double));			 
  argperi  = (double*) calloc(n,sizeof(double));
  meananom = (double*) calloc(n,sizeof(double));
  elements = (Elements*) calloc(n,sizeof(Elements));

  flags = (int*) calloc(n,sizeof(int));
  dt = (double*) calloc(n,sizeof(double));

  State state;

  if(p==NULL || a==NULL || e==NULL || incl==NULL || longnode==NULL || argperi==NULL || meananom==NULL){
    exit(0);
  }

  machine_epsilon = 2e-15;

  double dM = 2.0*PI/n;

  for(i=0; i<n; i++){
    a[i] = 1.0;
    e[i] = 0.99;
    incl[i] = 1.0;
    longnode[i] = 0.0;
    argperi[i]  = 0.1;
    meananom[i] = i*dM;
    elements[i].a = a[i];
    elements[i].e = e[i];
    elements[i].incl = incl[i];
    elements[i].longnode = longnode[i];
    elements[i].argperi = argperi[i];
    elements[i].meananom = meananom[i];
    dt[i] = 0.1*i;
  }

  //cartesians(n, GMsun, a, e, incl, longnode, argperi, meananom, p);

  cartesian_elements(n, GMsun, elements, positions, velocities);

  {
    double a_t, e_t, incl_t, longnode_t, argperi_t, meananom_t;
    for(i=0; i<n; i++){
      state.x = positions[i*3+0];
      state.y = positions[i*3+1];
      state.z = positions[i*3+2];
      state.xd = velocities[i*3+0];
      state.yd = velocities[i*3+1];
      state.zd = velocities[i*3+2];            
      keplerian(GMsun, state, &a_t, &e_t, &incl_t, &longnode_t, &argperi_t, &meananom_t);
      printf("%d %lf %lf %lf %lf %lf %lf\n",
	     i, a_t, e_t, incl_t, longnode_t, argperi_t, meananom_t);
      printf("%d %lf %lf %lf %lf %lf %lf\n",
	     i, elements[i].a, elements[i].e, elements[i].incl, elements[i].longnode, elements[i].argperi, elements[i].meananom);
    }
  }

}

