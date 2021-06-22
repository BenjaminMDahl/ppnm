#include<stdio.h>
#include<time.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_sort_vector.h>
#include<assert.h>
#define RND (double)rand()/RAND_MAX



void vector_print(char s[], gsl_vector* v);

void newton(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double eps);

void secular_equation_solve(gsl_vector* D,gsl_vector* u, int p ,gsl_vector* x, double eps){
	gsl_vector* Dcopy=gsl_vector_alloc(D->size);
	gsl_vector_memcpy(Dcopy,D);
	gsl_sort_vector(Dcopy);
	double a,xi,di,dI,gi; gsl_blas_ddot(u,u,&a); //u^T*u
	void ses(gsl_vector* v, gsl_vector* f){
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


	gsl_vector* g=gsl_vector_alloc(1);
	//Den første indgang er lidt specil
	di=gsl_vector_get(Dcopy,0)-fabs(a);
	dI=gsl_vector_get(Dcopy,0);
	gi=di+(dI-di)/2;
	while(1){
	gsl_vector_set(g,0,gi);
	newton(ses,g,eps);
	xi=gsl_vector_get(g,0);
	if(di <=xi && xi <=dI){
			gsl_vector_set(x,0,xi);
			break;
		}
		else if(xi<di) gi=di+(gi-di)/2;
		else gi=dI-(dI-gi)/2;}


	// så de næste mange med undtagelse af den sidste
	for(int i=0;i<p;i++){
	di=gsl_vector_get(Dcopy,i);
	dI=gsl_vector_get(Dcopy,i+1);
	gi=di+(dI-di)/2;
	while(1){
		gsl_vector_set(g,0,gi);
		printf("Gæt: \n %g \n",gi);
		newton(ses,g,eps);
		xi=gsl_vector_get(g,0);
		if(di <=xi && xi <=dI){
			gsl_vector_set(x,i+1,xi);
			break;
		}
		else if(xi<di) gi=di+(gi-di)/2;
		else gi=dI-(dI-gi)/2;}
	}

	for(int i=p+1;i<x->size-1;i++){
	di=gsl_vector_get(Dcopy,i);
	dI=gsl_vector_get(Dcopy,i+1);
	gi=di+(dI-di)/2;
	while(1){
		gsl_vector_set(g,0,gi);
		printf("Gæt: \n %g \n",gi);
		newton(ses,g,eps);
		xi=gsl_vector_get(g,0);
		if(di <=xi && xi <=dI){
			gsl_vector_set(x,i,xi);
			break;
		}
		else if(xi<di) gi=di+(gi-di)/2;
		else gi=dI-(dI-gi)/2;}
	}

	// Den sidste
	di=gsl_vector_get(Dcopy,Dcopy->size-1);
	dI=gsl_vector_get(Dcopy,Dcopy->size-1)+fabs(a);
	gi=di+(dI-di)/2;
	while(1){
	gsl_vector_set(g,0,gi);
	newton(ses,g,eps);
	xi=gsl_vector_get(g,0);
	if(di <=xi && xi <=dI){
			gsl_vector_set(x,x->size-1,xi);
			break;
		}
		else if(xi<di) gi=di+(gi-di)/2;
		else gi=dI-(dI-gi)/2;}

}


int main(){

clock_t start, end;
//double time_5,time_10,time_100,time_1000;
double time_5;

int n=5;
gsl_vector* D=gsl_vector_alloc(n);
gsl_vector* u=gsl_vector_alloc(n);
int p=0; //Overvej lige noget index her
double tal[]={1,1,3,4,5,6,7,8,9,10};

for(int i=0;i<n;i++){
	gsl_vector_set(D,i,tal[i]);
	gsl_vector_set(u,i,tal[i]);}
//gsl_vector_scale(u,1/sqrt(55));
//gsl_matrix_scale(D,1/sqrt(55));


gsl_vector* x=gsl_vector_alloc(n);

start = clock();

	secular_equation_solve(D,u,p,x,0.0000001);
	vector_print("res er",x);

end = clock();
time_5= ((double) (end - start)) / CLOCKS_PER_SEC;
printf("det tog: \n %g \n",time_5);

// Prøver igen
/*
int m=10;
gsl_vector* d=gsl_vector_alloc(m);
gsl_vector* U=gsl_vector_alloc(m);
gsl_vector* xb=gsl_vector_alloc(m);
int P=1; //Overvej lige noget index her

for(int i=0;i<m;i++){
	gsl_vector_set(d,i,tal[i]);
	gsl_vector_set(U,i,tal[i]);}

//double b[]={1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5};

//start = clock();


secular_equation_solve(d,U,P,xb,0.0001);
vector_print("res er",xb);

//end = clock();
//time_10= ((double) (end - start)) / CLOCKS_PER_SEC;
*/

/*
// Prøver igen

int l=100;
gsl_matrix* D100=gsl_matrix_alloc(l,l);
gsl_vector* u100=gsl_vector_alloc(l);
int p100=25; //Overvej lige noget index her
double c[100];

for(int i=0;i<l;i++){
	gsl_matrix_set(D100,i,i,RND);
	gsl_vector_set(u100,i,RND);
	c[i]=RND;}

start = clock();
for(int i=0;i<l;i++){
	gsl_vector_set(x,0,c[i]);
	secular_equation_solve(D100,u100,p100,x,0.0001);
//	vector_print("res er",x);
	}
end = clock();

time_100= ((double) (end - start)) / CLOCKS_PER_SEC;

// prøver igen igen nu med 1000
int L=1000;
gsl_matrix* D1000=gsl_matrix_alloc(L,L);
gsl_vector* u1000=gsl_vector_alloc(L);
double C[1000];

for(int i=0;i<L;i++){
	gsl_matrix_set(D1000,i,i,RND);
	gsl_vector_set(u1000,i,RND);
	C[i]=RND;}

start = clock();
for(int i=0;i<L;i++){
	gsl_vector_set(x,0,C[i]);
	secular_equation_solve(D1000,u1000,55,x,0.0001);
//	vector_print("res er",x);
	}
end = clock();

time_1000= ((double) (end - start)) / CLOCKS_PER_SEC;


// Tid
printf("Tid1000/tid100 :\n %g \n",time_1000/time_100);
*/
return 0;
}
