#include <math.h>
#include <stdio.h>
#define NRANSI
#include "nrutil.h"
#define ITMAX 200
#define EPS 1.0e-10
#define FREEALL free_vector_double(xi,1,n);free_vector_double(h,1,n);free_vector_double(g,1,n);
extern int DEBUGGING_NR;

void frprmn(double p[], int n, double ftol, int *iter, double *fret,
	double (*func)(double []), void (*dfunc)(double [], double []))
{
        void linmin(double p[], double xi[], int n, double *fret,
		double (*func)(double []));
	int j,its;
	double gg,gam,fp,dgg;
	double *g,*h,*xi;

	if (DEBUGGING_NR) {
	  fprintf(stderr,"Entered frprmn with n=%d\n",n);
	}
	g=vector_double(1,n);
	h=vector_double(1,n);
	xi=vector_double(1,n);
	if (DEBUGGING_NR) {
	  fprintf(stderr,"Done vector(1,%d)\n",n);
	}
	fp=(*func)(p);
	if (DEBUGGING_NR) {
	  fprintf(stderr,"set fp=(*func)(p)\n");
	}
	(*dfunc)(p,xi);
	if (DEBUGGING_NR) {
	  fprintf(stderr,"done *dfunc\n");
	}

	for (j=1;j<=n;j++) {
		g[j] = -xi[j];
		xi[j]=h[j]=g[j];
	}
	if (DEBUGGING_NR) {
	  fprintf(stderr,"Done j=1..n\n");
	}
	for (its=1;its<=ITMAX;its++) {
		*iter=its;
		if (DEBUGGING_NR) {
		  fprintf(stderr,"Calling linmin\n");
		}
		linmin(p,xi,n,fret,func);
		if (DEBUGGING_NR) {
		  fprintf(stderr,"returned from linmin\n");
		}

		if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
			FREEALL
			  if (DEBUGGING_NR) {
			    fprintf(stderr,"returning from frprmn\n");
			  }
			return;
		}
		fp= (*func)(p);
		(*dfunc)(p,xi);

		dgg=gg=0.0;
		for (j=1;j<=n;j++) {
			gg += g[j]*g[j];
			dgg += (xi[j]+g[j])*xi[j];
		}
		if (gg == 0.0) {
			FREEALL
			  if (DEBUGGING_NR) {
			    fprintf(stderr,"2:returning from frprmn\n");
			  }
			return;
		}
		gam=dgg/gg;
		for (j=1;j<=n;j++) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j]+gam*h[j];
		}
	}
	nrerror("Too many iterations in frprmn");
}
#undef ITMAX
#undef EPS
#undef FREEALL
#undef NRANSI
