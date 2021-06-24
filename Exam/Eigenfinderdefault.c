#include<stdio.h>
#include<time.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_sort_vector.h>
#include<assert.h>

void newton(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double eps);

void secular_equation_default(gsl_vector* D,gsl_vector* u, int p ,gsl_vector* x, double eps,int m){
	assert(D->size>1); // Vi arbejder kun med matricer ikke skalar
	assert(m>=0);
	// Vi starter med at opdateret D, så hvis den p'te indgang i u ikke er 0, updateres D og u sættes lig 0.
	double dstart=gsl_vector_get(D,p);
	double ustart=gsl_vector_get(u,p);
	gsl_vector_set(D,p,dstart+2*ustart);
	gsl_vector_set(u,p,0);

	// Vi skal have en ordnet udgave af D for at kunne bestemme intevallerne vi vil lede efter eigenvalues i.
	gsl_vector* Dcopy=gsl_vector_alloc(D->size);
	gsl_vector_memcpy(Dcopy,D);
	double xp=gsl_vector_get(Dcopy,p);
	gsl_sort_vector(Dcopy);
	//Vi finder hvilket index p  svarer til i den sorterede liste
	int pnew=0;
	while(pnew<=Dcopy->size){
	if(gsl_vector_get(Dcopy,pnew)==xp)break;
	pnew++;}

	// Erklære de relevante doubles og laver funktionen Se som ud fra D og u danner den relevante Secular Equation(4.30 fra bogen)
	double a,xi,di,dI,gi;
	int I;
	gsl_blas_ddot(u,u,&a); //u skal bruges som øvre og nedre grænse for hvor vores egenværdier er

	void Se(gsl_vector* v, gsl_vector* f){
		double x=gsl_vector_get(v,0);
		double S=x-gsl_vector_get(D,p);
		double P=1;
		for(int k=0;k<p;k++){
			double dk=gsl_vector_get(D,k);
			double uk=gsl_vector_get(u,k);
			P=P*(dk-x);
			S=S+uk*uk/(dk-x);}
		for(int k=p+1;k<u->size;k++){
			double dk=gsl_vector_get(D,k);
			double uk=gsl_vector_get(u,k);
			double P=P*(dk-x);
			S=S+uk*uk/(dk-x);}
		gsl_vector_set(f,0,S*P);}


	// Da vi bruger vores newton metode lavet i kurset, skal alt leveres i vektorer, da metoden er lavet til at kunne klare problemer af højere dimension.
	gsl_vector* g=gsl_vector_alloc(1);

	//Den mindste egenværdi bestemmes. Det forventes at den er mindre end alle indegange i D
	gsl_vector_set(x,0,NAN); 	//Hvis der ikke findes en eigenvalue
	di=gsl_vector_get(Dcopy,0)-fabs(a);
	dI=gsl_vector_get(Dcopy,0);
	gi=di+(dI-di)/2;
	I=0;
	while(I<=2*m+1){
		I++;
		gsl_vector_set(g,0,gi);
		newton(Se,g,eps); 	// Min egen Newton-Raphson kode fra homework 8 "Roots"
		xi=gsl_vector_get(g,0);
	// Tester om den funde værdi er i det ønskede interval ellers prøver den igen med et nyt gæt
		if(di <=xi && xi <=dI){
			gsl_vector_set(x,0,xi);
			break;}
		else{
			if(I%2==0)gi=di+(dI-di)/2+I*(dI-di)/(4*m+1);
			else gi=di+(dI-di)/2-I*(dI-di)/(4*m+1);}
		}

	// Nu bestemmes alle de eigenværdier som forventes at være bundet af indgangene i D.
	// Intervallet efter d_p springes over, derfor er der lavet to forlykker
	for(int i=0;i<pnew;i++){
	gsl_vector_set(x,i+1,NAN);
	di=gsl_vector_get(Dcopy,i);
	dI=gsl_vector_get(Dcopy,i+1);
	gi=di+(dI-di)/2;
	I=0;
	if(di==dI)gsl_vector_set(x,i+1,di);		// Hvis der er ens indgange i D vil de indgange være egenværdier
	else{
		while(I<=2*m+1){
			I++;
			gsl_vector_set(g,0,gi);
			newton(Se,g,eps);
			xi=gsl_vector_get(g,0);
			if(di <=xi && xi <=dI){
				gsl_vector_set(x,i+1,xi);
				break;}
			else{
				if(I%2==0) gi=di+(dI-di)/2+I*(dI-di)/(4*m+1);
				else gi=di+(dI-di)/2-I*(dI-di)/(4*m+1);}
		}}}

	for(int i=pnew+1;i<x->size-1;i++){
	gsl_vector_set(x,i,NAN);
	di=gsl_vector_get(Dcopy,i);
	dI=gsl_vector_get(Dcopy,i+1);
	gi=di+(dI-di)/2;
	I=0;
	if(di==dI)gsl_vector_set(x,i,di);
	else{
		while(I<=2*m+1){
			I++;
			gsl_vector_set(g,0,gi);
			newton(Se,g,eps);
			xi=gsl_vector_get(g,0);
			if(di <=xi && xi <=dI){
				gsl_vector_set(x,i,xi);
				break;}
			else{
				if(I%2==0)gi=di+(dI-di)/2+I*(dI-di)/(4*m+1);
				else gi=di+(dI-di)/2-I*(dI-di)/(4*m+1);}
			}}}

	// Den sidste og største egenværdi bestemmes
	gsl_vector_set(x,x->size-1,NAN);
	di=gsl_vector_get(Dcopy,Dcopy->size-1);
	dI=gsl_vector_get(Dcopy,Dcopy->size-1)+fabs(a);
	gi=di+(dI-di)/2;
	I=0;
	while(I<=2*m+1){
		I++;
		gsl_vector_set(g,0,gi);
		newton(Se,g,eps);
		xi=gsl_vector_get(g,0);
		if(di <=xi && xi <=dI){
			gsl_vector_set(x,x->size-1,xi);
			break;}
		else{
			if(I%2==0)gi=di+(dI-di)/2+I*(dI-di)/(4*m+1);
			else gi=di+(dI-di)/2-I*(dI-di)/(4*m+1);}
		}

	//Der frigives hukommelse, og nu burde alle egenværdier være fundet og gemt i x.
	gsl_vector_free(g);
	gsl_vector_free(Dcopy);
}

