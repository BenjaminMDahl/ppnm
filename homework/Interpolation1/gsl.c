#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_interp.h>

#include <assert.h>

int main (void) {
// Min data men ellers blot koden fra eksemplet
	double x[]={0.06,0.36,0.50,1.08,1.22,1.90,2.16,2.8,3.32,3.66,4.00,4.56,4.74,5.1,5.18,5.56};
	double y[]={0.78,0.92,2.72,3.72,4.96,5.56,6.34,5.86,5.30,4.60,4.00,3.14,2.26,1.00,0.32,0.00};
	int n=sizeof(x)/sizeof(x[0]);
	int N=sizeof(y)/sizeof(y[0]); assert(n==N);
	printf("# index 0: data\n");
	for(int i=0;i<n;i++) printf("%g %g\n",x[i],y[i]);
	printf("\n\n");

	gsl_interp* linear     = gsl_interp_alloc(gsl_interp_linear    ,n);
	gsl_interp* cspline    = gsl_interp_alloc(gsl_interp_cspline   ,n);
	gsl_interp* polynomial = gsl_interp_alloc(gsl_interp_polynomial,n);
	gsl_interp_init(linear    ,x,y,n);
	gsl_interp_init(cspline   ,x,y,n);
	gsl_interp_init(polynomial,x,y,n);

	printf("# index 1: interpolations\n");
	for(double z=x[0];z<=x[n-1];z+=1./16){
		double interp_l=gsl_interp_eval(linear    ,x,y,z,NULL);
		double interp_c=gsl_interp_eval(cspline   ,x,y,z,NULL);
		double interp_p=gsl_interp_eval(polynomial,x,y,z,NULL);
		printf("%g %g %g %g\n",z,interp_l,interp_c,interp_p);
		}
	printf("\n\n");

	printf("# index 2: integrals\n");
	for(double z=x[0];z<=x[n-1];z+=1./16){
		double integ_l=gsl_interp_eval_integ(linear    ,x,y,x[0],z,NULL);
		double integ_c=gsl_interp_eval_integ(cspline   ,x,y,x[0],z,NULL);
		double integ_p=gsl_interp_eval_integ(polynomial,x,y,x[0],z,NULL);
		printf("%g %g %g %g\n",z,integ_l,integ_c,integ_p);
		}

gsl_interp_free(linear);
gsl_interp_free(cspline);
gsl_interp_free(polynomial);
return 0;
}
