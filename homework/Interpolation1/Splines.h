#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>


int binsearch(int n, double* x, double z);


		// linear spline //

double linterp(int n, double x[], double y[], double z);

double linterp_integral(int n, double x[], double y[], double z);


	// Qudratic spline //

typedef struct {int n ; double *x , *y , *b , *c;} qspline; 	// Der definieres en struct (qspline) hvor alt vores data skal opsamles i

qspline* qspline_alloc ( int n , double* x , double* y );


double qspline_eval(qspline *s, double z);

double qua_integral(qspline *s, double z);

double qua_diff(qspline *s, double z);


	// Cubic spline //

typedef struct {int n; double *x,*y,*b,*c,*d;} cspline;					// Der skal igen dannes en ny strucktur med et parameter mere, men ellers gør vi meget det samme som ved qspline

cspline* cspline_alloc(int n, double *x, double *y);
double cspline_eval(cspline *s,double z);

double cub_integral(cspline *s, double z);

double cub_diff(cspline *s, double z);


	// CLEAN UP //
	void qspline_free(qspline *s);
	void cspline_free(cspline *s);
	
