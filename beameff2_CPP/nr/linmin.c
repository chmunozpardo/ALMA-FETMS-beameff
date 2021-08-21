/* note #undef's at end of file */
#define NRANSI
#include <stdio.h>
#include "nrutil.h"
#define TOL 2.0e-4
extern int DEBUGGING_NR;
int ncom;
double *pcom,*xicom,(*nrfunc)(double []);

void linmin(double p[], double xi[], int n, double *fret, double (*func)(double []))
{
	double brent(double ax, double bx, double cx,
		double (*f)(double), double tol, double *xmin);
	double f1dim(double x);
	void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
		double *fc, double (*func)(double));
	int j;
	double xx,xmin,fx,fb,fa,bx,ax;

	if (DEBUGGING_NR) {
	  fprintf(stderr,"Inside linmin with n=%d\n",n);
	}
	ncom=n;
	if (DEBUGGING_NR) {
	  fprintf(stderr,"Calling vector(1,%d)\n",n);
	}
	pcom=vector_double(1,n);
	if (DEBUGGING_NR) {
	  fprintf(stderr,"done vector(1,%d)\n",n);
	}
	xicom=vector_double(1,n);
	if (DEBUGGING_NR && 0) {
	  fprintf(stderr,"done vector(1,%d)\n",n);
	}
	nrfunc=func;
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	if (DEBUGGING_NR && 0) {
	  fprintf(stderr,"calling mnbrak\n");
	}
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
	if (DEBUGGING_NR && 0) {
	  fprintf(stderr,"returned from mnbrak\n");
	}

	*fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
	if (DEBUGGING_NR) {
	  fprintf(stderr,"returned from brent\n");
	}

	for (j=1;j<=n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	free_vector_double(xicom,1,n);
	free_vector_double(pcom,1,n);
}
#undef TOL
#undef NRANSI
