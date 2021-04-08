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
	gsl_vector* dy){             				// output: error estimate

	int i, n=yt->size;
	gsl_vector* k0=gsl_vector_alloc(n); gsl_vector* k12=gsl_vector_alloc(n); gsl_vector* y_half=gsl_vector_alloc(n);
	f(t,yt,k0); // k'erne dannes
	for(i=0;i<n;i++)gsl_vector_set(y_half,i,gsl_vector_get(yt,i)+gsl_vector_get(k0,i)*h/2);
	f(t+h/2,y_half,k12);
	for(i=0;i<n;i++){
		gsl_vector_set(yh,i,gsl_vector_get(yt,i)+gsl_vector_get(k12,i)*h);
		gsl_vector_set(dy,i,(gsl_vector_get(k12,i)-gsl_vector_get(k0,i))*h/2);}

	gsl_vector_free(k0);gsl_vector_free(k12);gsl_vector_free(y_half);
}


int OdeDriver(
	void f(double,gsl_vector*,gsl_vector*), 		// right-hand-side of dy/dt=f(t,y)
	double a,			                     	// the start-point a
	double b,                     				// the end-point of the integration
	gsl_vector* ya,                    				// y(a)
	gsl_vector* yb,                     			// y(b) to be calculated
	double h,                     				// initial step-size
	double acc,                   				// absolute accuracy goal
	double eps){                    			// relative accuracy goal

	int i, k=0, n=ya->size;
	double syh, sdy, err, norm_y, tol, xi=a;
	gsl_vector* dy=gsl_vector_alloc(n);
	gsl_vector* ybcopy=gsl_vector_alloc(n);
	gsl_vector_memcpy(yb,ya);

	do
	{
		if(xi+h>b) h=b-xi;	//Tjekker om det oplyste første step er for stort
		gsl_vector_memcpy(ybcopy,yb);
		rkstep12(f,xi,ybcopy,h,yb,dy);
		syh=0;sdy=0; for(i=0;i<n;i++){
		sdy+=gsl_vector_get(dy,i);
		syh+=gsl_vector_get(yb,i);}
		err=sqrt(sdy*sdy); norm_y=sqrt(syh*syh);
		tol=(norm_y*eps+acc)*sqrt(h/(b-a));
		if(err<tol){
			xi+=h;
			k++;}
		if(err==0) h*=2; // Vi skal lige passe på ikke at dele med 0
		else h*=pow(tol/err,0.25)*0.95;
	}
	while(xi<b);

	gsl_vector_free(ybcopy); gsl_vector_free(dy);
return k;
}

int OdeDriverRecorder(
	void f(double,gsl_vector*,gsl_vector*), 		// right-hand-side of dy/dt=f(t,y)
	double a,			                     	// the start-point a
	double b,                     				// the end-point of the integration
	gsl_vector* ya,                    				// y(a)
	gsl_vector* yb,                     			// y(b) to be calculated
	double h,                     				// initial step-size
	double acc,                   				// absolute accuracy goal
	double eps){                    			// relative accuracy goal

	int i, k=0, n=ya->size;
	double syh, sdy, err, norm_y, tol, xi=a;
	gsl_vector* dy=gsl_vector_alloc(n);
	gsl_vector* ybcopy=gsl_vector_alloc(n);
	gsl_vector_memcpy(yb,ya);
	FILE *path;
	path = fopen("path.txt","w");
	do
	{
		if(xi+h>b) h=b-xi;	//Tjekker om det oplyste første step er for stort
		gsl_vector_memcpy(ybcopy,yb);
		rkstep12(f,xi,ybcopy,h,yb,dy);
		syh=0;sdy=0; for(i=0;i<n;i++){
		sdy+=gsl_vector_get(dy,i);
		syh+=gsl_vector_get(yb,i);}
		err=sqrt(sdy*sdy); norm_y=sqrt(syh*syh);
		tol=(norm_y*eps+acc)*sqrt(h/(b-a));
		if(err<tol){
			xi+=h;
			k++;
			fprintf(path,"%g %g\n",xi,gsl_vector_get(yb,0));}
		if(err==0) h*=2; // Vi skal lige passe på ikke at dele med 0
		else h*=pow(tol/err,0.25)*0.95;

	}
	while(xi<b);

	fclose(path);
	gsl_vector_free(ybcopy); gsl_vector_free(dy);
return k;
}

//Forsøgs kaninen u''=-u
void u(double t,gsl_vector* y, gsl_vector* dydt){
	gsl_vector_set(dydt,0,gsl_vector_get(y,1));
	gsl_vector_set(dydt,1,-gsl_vector_get(y,0));
}

/* void u1(double t,gsl_vector* y, gsl_vector* dydt){
	gsl_vector_set(dydt,0,gsl_vector_get(y,1));
	gsl_vector_set(dydt,1,6);}
*/

void epidemic(double t,gsl_vector* y, gsl_vector* dydt){
	double N=5500000, Tr=7, Tc=2;
	gsl_vector_set(dydt,0,-(gsl_vector_get(y,1)*gsl_vector_get(y,0))/(N*Tc));
	gsl_vector_set(dydt,1,(gsl_vector_get(y,1)*gsl_vector_get(y,0))/(N*Tc)-gsl_vector_get(y,1)/Tr);
	gsl_vector_set(dydt,2,gsl_vector_get(y,1)/Tr);
}

int main(){

	// u''=-u
	double a=0, b=M_PI*2, h=M_PI*2/2000, acc=0.0005, eps=0.001;
	gsl_vector* ya=gsl_vector_alloc(2);
	gsl_vector* yb=gsl_vector_alloc(2);
	gsl_vector* exact=gsl_vector_alloc(2);
	gsl_vector_set(ya,0,1); gsl_vector_set(ya,1,0);
	gsl_vector_set(exact,0,0); gsl_vector_set(exact,1,-1);

	int k=OdeDriverRecorder(u,a,b,ya,yb,h,acc,eps);

	vector_print("Min start vektor for u''=-u",ya);
	vector_print("Her er løsningen",yb);
	vector_print("Det analytiske svar er",exact);
	printf("Det tog %i skridt\n\n",k);
/*
	// u''=6 (svar y=x^3+1)
	double a1=0, b1=3, h1=5/100;
	gsl_vector* ya1=gsl_vector_alloc(2);
	gsl_vector* yb1=gsl_vector_alloc(2);
	gsl_vector_set(ya1,0,1); gsl_vector_set(ya1,1,3);
	int k1=OdeDriver(u1,a1,b1,ya1,yb1,h1,acc,eps);


	vector_print("A start",ya1);
	vector_print("Her er løsningen",yb1);
	printf("skridt 1 %i\n",k1);
*/

	// Epidemic

	double t0=0, tyear=365, he=0.0005, acce=0.0005, epse=0.001;
	double N=5500000, I=3200, R=226000, S=N-R-I;
	gsl_vector* yae=gsl_vector_alloc(3);
	gsl_vector* ybe=gsl_vector_alloc(3);
	gsl_vector_set(yae,0,S); gsl_vector_set(yae,1,I); gsl_vector_set(yae,2,R);

	int ke=OdeDriver(epidemic,t0,tyear,yae,ybe,he,acce,epse);

	vector_print("Epidemiens udgangspunkt hvor fra øverst til nederst vi har: S/I/R",yae);
	vector_print("Resultatet efter et år",ybe);
	printf("Det tog så mange skridt: %i \n",ke);


return 0;
}
