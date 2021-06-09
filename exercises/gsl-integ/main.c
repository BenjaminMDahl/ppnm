#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>
/* Vi ønsker at udregne int^1:0 log(x)/sqrt(x) dx */

double f (double x, void * params) {
	double a= *(double *) params;
	double f= log(a*x)/sqrt(x);

return f;
}

double myfun(double z){
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000); //Gør plads i hukommelsen til GSL integration

	//Gør de relevante doubles klar//
	double b=z;
 	double result, error;
 	double alpha = 1.0;

	//Laves funktionen f til en GSL_funktion med de relevante parametere(alpha)//
 	gsl_function F;
 	F.function = &f;
 	F.params = &alpha;

	// Bruger QAGS rotunien til at udregne integralet, svaret gemmes i result.
 	gsl_integration_qags (&F, 0, b, 0, 1e-7, 1000,w, &result, &error);

	// Befrirer hukommelsen brugt
 	gsl_integration_workspace_free (w);

return result;
}

int main(){
	//Gemmer for forskelige b for at kunne plotte
	for(double x=0.01;x<=1;x+=1.0/30)
		printf("%10g %10g\n",x,myfun(x));

return 0;
}
