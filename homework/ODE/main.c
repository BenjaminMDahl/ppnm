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

	f(t,yt,k0);
	for(i=0;i<n;i++){
		double	y_half_i=gsl_vector_get(yt,i)+gsl_vector_get(k0,i)*h/2;
		gsl_vector_set(y_half,i,y_half_i);}
	f(t+h/2,y_half,k12);
	for(i=0;i<n;i++){
		double yh_i=gsl_vector_get(yt,i)+gsl_vector_get(k12,i)*h;
		double dy_i=(gsl_vector_get(k0,i)-gsl_vector_get(k12,i))*h/2; //Disse bliver for små
		gsl_vector_set(yh,i,yh_i);
		gsl_vector_set(dy,i,dy_i);}

	gsl_vector_free(k0);gsl_vector_free(k12);gsl_vector_free(y_half);
}

void printdata(double x,gsl_vector* y,FILE *path){
int n=y->size;
fprintf(path,"%6g ",x);
for(int i=0;i<n;i++)fprintf(path,"%6g ",gsl_vector_get(y,i));
fprintf(path,"\n");
}


int OdeDriver(
	void f(double,gsl_vector*,gsl_vector*), 		// right-hand-side of dy/dt=f(t,y)
	double a,			                     	// the start-point a
	double b,                     				// the end-point of the integration
	gsl_vector* ya,                    				// y(a)
	gsl_vector* yb,                     			// y(b) to be calculated
	double h,                     				// initial step-size
	double acc,                   				// absolute accuracy goal
	double eps,
	char *program){                    			// relative accuracy goal

	int k=0, n=ya->size;
	double err, norm_y, tol, xi=a;
	gsl_vector* dy=gsl_vector_alloc(n);
	gsl_vector_memcpy(yb,ya);
	FILE *path;
	path = fopen(program,"w");
	printdata(xi,yb,path);
	while(xi<b)
	{
		if(xi+h>b) h=b-xi;	//Tjekker om det oplyste første step er for stort
		rkstep12(f,xi,ya,h,yb,dy);
		err=gsl_blas_dnrm2(dy);		//Summer alle dy'er op og tage absolut værdi for at få et udtryk for error
		norm_y=gsl_blas_dnrm2(yb);	//Summer alle y værdier op da de er relevante for relative uncertainty(bruges i tolerance)
		tol=(norm_y*eps+acc)*sqrt(h/(b-a));
		if(err<tol){
			xi+=h;
			gsl_vector_memcpy(ya,yb);
			k++;
			printdata(xi,yb,path);}
		if(err>0) h*=pow(tol/err,0.25)*0.95 ; // Vi skal lige passe på ikke at dele med 0
		else h*=2;
	}


	fclose(path);
	gsl_vector_free(dy);
return k;
}

//Forsøgs kaninen u''=-u
void u(double t,gsl_vector* y, gsl_vector* dydt){
	double y1=gsl_vector_get(y,1);
	double y2=-gsl_vector_get(y,0);
	gsl_vector_set(dydt,0,y1);
	gsl_vector_set(dydt,1,y2);
}

// konstant
void u1(double t,gsl_vector* y, gsl_vector* dydt){
	double y0=1;
	gsl_vector_set(dydt,0,y0);}


void epidemic(double t,gsl_vector* y, gsl_vector* dydt){
	double N=5500000, Tr=7, Tc=2;
	double y0=-(gsl_vector_get(y,1)*gsl_vector_get(y,0))/(N*Tc);
	double y1=(gsl_vector_get(y,1)*gsl_vector_get(y,0))/(N*Tc)-gsl_vector_get(y,1)/Tr;
	double y2=gsl_vector_get(y,1)/Tr;
	gsl_vector_set(dydt,0,y0);				//dS/dt=-(I*S)/(N*T_c)
	gsl_vector_set(dydt,1,y1);		//dI/dt=(I*S)/(N*T_c)-I/T_r
	gsl_vector_set(dydt,2,y2);								//dR/dt=I/T_r
}




int main(){

	printf("OPGAVE A Embedded Runge-Kutta ODE integrator\n \n");


	double a1=0, b1=3, h1=(double)1/10, acc1=0.001, eps1=0.001;
	gsl_vector* ya1=gsl_vector_alloc(1);
	gsl_vector* yb1=gsl_vector_alloc(1);
	gsl_vector_set(ya1,0,1);
	int k1=OdeDriver(u1,a1,b1,ya1,yb1,h1,acc1,eps1,"line.txt");
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

	int k=OdeDriver(u,a,b,ya,yb,h,acc,eps,"cos.txt");

	printf("Vi prøver nu noget lidt mere spændende, hvor vi løser en diff.ligning, hvor vi ved at løsningen er cos(x).\n vi har altså y(0)=1 og d^2y/dt^2=-y\n");
	vector_print("For dette er start vektoren",ya);
	vector_print("Vi får efter 2 pi følgende vektor",yb);
	vector_print("Det analytiske svar er",exact);
	printf("Dette giver en meget passende overensstemelse, og det er også plottet i cos.png \n Det tog %i skridt\n\n",k);




	double t0=0, tdays=30, he=0.1, acce=0.0000001, epse=0.0000001;
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
