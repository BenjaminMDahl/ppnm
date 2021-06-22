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

void secular_equation_default(gsl_vector* D,gsl_vector* u, int p ,gsl_vector* x, double eps){
	assert(D->size>1); // Vi arbejder kun med matricer ikke skalar

	// Vi skal have en ordnet udgave af D for at kunne bestemme intevallerne vi vil lede efter eigenvalues i.
	gsl_vector* Dcopy=gsl_vector_alloc(D->size);
	gsl_vector_memcpy(Dcopy,D);
	gsl_sort_vector(Dcopy);


	// Erklære de relevante doubles og laver funktionen Se som ud fra D og u danner den relevante Secular Equation(4.30 fra bogen)
	double a,xi,di,dI,gi,I; gsl_blas_ddot(u,u,&a); //u skal bruges som øvre og nedre grænse for hvor vores egenværdier er

	void Se(gsl_vector* v, gsl_vector* f){
		double x=gsl_vector_get(v,0);
		double S=x-gsl_vector_get(D,p);
		for(int k=0;k<p;k++){
			double dk=gsl_vector_get(D,k);
			double uk=gsl_vector_get(u,k);
			S=S+uk*uk/(dk-x);}
		for(int k=p+1;k<u->size;k++){
			double dk=gsl_vector_get(D,k);
			double uk=gsl_vector_get(u,k);
			S=S+uk*uk/(dk-x);}
		gsl_vector_set(f,0,S);}


	// Da vi bruger vores newton metode lavet i kurset, skal alt leveres i vektorer, da metoden er lavet til at kunne klare problemer af højere dimension.
	gsl_vector* g=gsl_vector_alloc(1);

	//Den mindste egenværdi bestemmes. Det forventes at den er mindre end alle indegange i D
	di=gsl_vector_get(Dcopy,0)-fabs(a);
	dI=gsl_vector_get(Dcopy,0);
	gi=di+(dI-di)/2;
	I=0;
	while(1){
	I++;
	gsl_vector_set(g,0,gi);
	newton(Se,g,eps); 	// Min egen Newton-Raphson kode fra homework 8 "Roots"
	xi=gsl_vector_get(g,0);
	// Tester om den funde værdi er i det ønskede interval ellers prøver den igen med et nyt gæt i retningen af den funde værdi
	if(I<=10){
		if(di <=xi && xi <=dI){
				gsl_vector_set(x,0,xi);
				break;
			}
			else if(xi<di) gi=di+(gi-di)/2;
			else gi=dI-(dI-gi)/2;}
	else{
		gsl_vector_set(x,0,NAN);
		break;}
	}

	// Nu bestemmes alle de eigenværdier som forventes at være bundet af indgangene i D.
	// Intervallet efter d_p springes over, derfor er der lavet to forlykker
	for(int i=0;i<p;i++){
	di=gsl_vector_get(Dcopy,i);
	dI=gsl_vector_get(Dcopy,i+1);
	gi=di+(dI-di)/2;
	I=0;
	while(1){
		I++;
		gsl_vector_set(g,0,gi);
		newton(Se,g,eps);
		xi=gsl_vector_get(g,0);
		if(I<=10){
			if(di <=xi && xi <=dI){
				gsl_vector_set(x,i+1,xi);
				break;
			}
			else if(xi<di) gi=di+(gi-di)/2;
			else gi=dI-(dI-gi)/2;}
		else{gsl_vector_set(x,i+1,NAN);
			break;}
			}
	}

	for(int i=p+1;i<x->size-1;i++){
	di=gsl_vector_get(Dcopy,i);
	dI=gsl_vector_get(Dcopy,i+1);
	gi=di+(dI-di)/2;
	I=0;
	while(1){
		I++;
		gsl_vector_set(g,0,gi);
		newton(Se,g,eps);
		xi=gsl_vector_get(g,0);
		if(I<=10){
			if(di <=xi && xi <=dI){
				gsl_vector_set(x,i,xi);
				break;
			}
			else if(xi<di) gi=di+(gi-di)/2;
			else gi=dI-(dI-gi)/2;}
		else{gsl_vector_set(x,i,NAN);
		break;}
		}
	}

	// Den sidste og største egenværdi bestemmes
	di=gsl_vector_get(Dcopy,Dcopy->size-1);
	dI=gsl_vector_get(Dcopy,Dcopy->size-1)+fabs(a);
	gi=di+(dI-di)/2;
	I=0;
	while(1){
	I++;
	gsl_vector_set(g,0,gi);
	newton(Se,g,eps);
	xi=gsl_vector_get(g,0);
	if(I<=10){
		if(di <=xi && xi <=dI){
			gsl_vector_set(x,x->size-1,xi);
			break;
		}
		else if(xi<di) gi=di+(gi-di)/2;
		else gi=dI-(dI-gi)/2;}
	else{
		gsl_vector_set(x,x->size-1,NAN);
		break;}
	}
	//Der frigives hukommelse, og nu burde alle egenværdier være fundet og gemt i x.
	gsl_vector_free(g);
	gsl_vector_free(Dcopy);
}

