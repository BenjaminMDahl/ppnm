#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>

struct my_f_params { double a; };

double fun(double x, void* p){
	struct my_f_params * params = (struct my_f_params *)p;
	double a= (params->a);
	double y= a*sqrt(1-x*x);

return y;
}

double gsl_int(double a, double b, double delta, double eps){

	gsl_integration_workspace * w  = gsl_integration_workspace_alloc (1000);
	double result, error;
	struct my_f_params par= {4.0};

	gsl_function F;
	F.function = fun;
 	F.params= &par;


	gsl_integration_qags (&F, a, b, delta, eps, 1000,
                        w, &result, &error);
	
	gsl_integration_workspace_free (w);
return result;
}
