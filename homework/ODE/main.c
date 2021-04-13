#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>


void vector_print(char s[], gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i< v->size ;i++)printf("%10g \n",gsl_vector_get(v,i));
	printf("\n");
}


// Inspireret af kapitelet om ODE
void rkstep12(
	void f(double t,gsl_vector* y,gsl_vector* dydt), 	// the f from dy/dt=f(t,y)
	double t,              					// the current value of the variable
	gsl_vector* yt,   				        // the current value y(t) of the sought function
	double h,   					        // the step to be taken
	gsl_vector* yh,					        // output: y(t+h)
	gsl_vector* dy);             				// output: error estimate



int OdeDriver(	void f(double,gsl_vector*,gsl_vector*), 		// right-hand-side of dy/dt=f(t,y)
	double a,			                     	// the start-point a
	double b,                     				// the end-point of the integration
	gsl_vector* ya,                    				// y(a)
	gsl_vector* yb,                     			// y(b) to be calculated
	double h,                     				// initial step-size
	double acc,                   				// absolute accuracy goal
	double eps,
	char *program);                    			// relative accuracy goal



//Forsøgs kaninen u''=-u
void odecos(double t,gsl_vector* y, gsl_vector* dydt);

// konstant
void line(double t,gsl_vector* y, gsl_vector* dydt);


void epidemic(double t,gsl_vector* y, gsl_vector* dydt);

int main(){

	printf("OPGAVE A Embedded Runge-Kutta ODE integrator\n \n");


	double a1=0, b1=3, h1=(double)1/10, acc1=0.001, eps1=0.001;
	gsl_vector* ya1=gsl_vector_alloc(1);
	gsl_vector* yb1=gsl_vector_alloc(1);
	gsl_vector_set(ya1,0,1);
	int k1=OdeDriver(line,a1,b1,ya1,yb1,h1,acc1,eps1,"line.txt");
	printf("Vi har startet med en simple diff.ligning, hvor vi har y(0)=1 dy/dt=1, altså en lige linje\n");
	vector_print("Vi har for denne fundet y(3) til",yb1);
	printf("Hvilket passer godt med en linje med hældning 1 som skære y-aksen i 1.\n Vi kom frem til dette resultat på %i steps \n",k1);
	printf("Vi har lavet et plot af linjen i figurer line.png.\n");

	// u''=-u
	double a=0, b=M_PI*2, h=0.0001, acc=0.001, eps=0.001;
	gsl_vector* ya=gsl_vector_alloc(2);
	gsl_vector* yb=gsl_vector_alloc(2);
	gsl_vector* exact=gsl_vector_alloc(2);
	gsl_vector_set(ya,0,1); gsl_vector_set(ya,1,0);
	gsl_vector_set(exact,0,0); gsl_vector_set(exact,1,-1);

	int k=OdeDriver(odecos,a,b,ya,yb,h,acc,eps,"cos.txt");

	printf("Vi prøver nu noget lidt mere spændende, hvor vi løser en diff.ligning, hvor vi ved at løsningen er cos(x).\n vi har altså y(0)=1 og d^2y/dt^2=-y\n");
	vector_print("For dette er start vektoren",ya);
	vector_print("Vi får efter 2 pi følgende vektor",yb);
	vector_print("Det analytiske svar er",exact);
	printf("Dette giver en meget passende overensstemelse, og det er også plottet i cos.png \n Det tog %i skridt\n\n",k);


	double t0=0, tdays=60, he=0.1, acce=0.0000001, epse=0.0000001;
	double N=5500000, I=3200, R=226000, S=N-R-I;
	gsl_vector* yae=gsl_vector_alloc(3);
	gsl_vector* ybe=gsl_vector_alloc(3);
	gsl_vector_set(yae,0,S); gsl_vector_set(yae,1,I); gsl_vector_set(yae,2,R);

	int ke=OdeDriver(epidemic,t0,tdays,yae,ybe,he,acce,epse,"Epidemic.txt");

	printf("Til sidst har vi prøver at bruge SIR modellen for at få lavet et plot over corona situationen i Danmark.\n");
	printf("Vi har valgt følgende parametre N=%g I=%g R=%g og undersøge udviklingen over %g dage\n",N,I,R,tdays);
	vector_print("Epidemiens udgangspunkt hvor fra øverst til nederst vi har: S/I/R",yae);
	vector_print("Resultatet efter et år",ybe);
	printf("Vi har plottet dette i figurer Epidemic.png og det ligner meget godt plottet fra Wikisiden om SIR modellen\n");
	printf("Det tog så mange skridt: %i \n",ke);


	printf("\n \n OPGAVE B Store the path\n \n");

	printf("Min løsning til at gemme vejen fra a til b, har været at give odedriver et argument mere, hvor den tage en string\n");
	printf("Denne string er navnet på den fil odedriver danner og gemme vejen i direkte på disken, der er 3 af disse filer en for\n");
	printf("den lige linje, en for cos og en for epidemic. det er også disse filer der er brugt til at danne png'erne i Makefile.\n");
	printf("Så længe den string man kalder odedriver med slutter på txt vil den også blive slettet igen med min clean function i Makefile.\n");

return 0;
}
