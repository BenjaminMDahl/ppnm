#include<stdio.h>
#include<time.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<assert.h>
#define RND (double)rand()/RAND_MAX



void vector_print(char s[], gsl_vector* v);

void newton(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double eps);

//gsl_matrix* D, gsl_vector* u, int p,

void secular_equation_solve(gsl_matrix* D,gsl_vector* u, int p ,gsl_vector* x, double eps){

	void ses(gsl_vector* v, gsl_vector* f){
		double x=gsl_vector_get(v,0);
		double S=x-gsl_matrix_get(D,p,p);
		for(int k=0;k<p;k++){
			double dk=gsl_matrix_get(D,k,k);
			double uk=gsl_vector_get(u,k);
			S=S+uk*uk/(dk-x);}
		for(int k=p+1;k<u->size;k++){
			double dk=gsl_matrix_get(D,k,k);
			double uk=gsl_vector_get(u,k);
			S=S+uk*uk/(dk-x);}
		gsl_vector_set(f,0,S);}


	newton(ses,x,eps);
}


int main(){
clock_t start, end;
double time_5,time_10,time_100,time_1000;


int n=5;
gsl_matrix* D=gsl_matrix_alloc(n,n);
gsl_vector* u=gsl_vector_alloc(n);
int p=1; //Overvej lige noget index her
double tal[]={1,2,3,4,5,6,7,8,9,10};

for(int i=0;i<n;i++){
	gsl_matrix_set(D,i,i,tal[i]);
	gsl_vector_set(u,i,tal[i]);}


gsl_vector* x=gsl_vector_alloc(1);

double a[]={1.5,2.5,3.5,4.5,5.5};

start = clock();
for(int i=0;i<5;i++){
	gsl_vector_set(x,0,a[i]);
	secular_equation_solve(D,u,p,x,0.0001);
//	vector_print("res er",x);
	}
end = clock();
time_5= ((double) (end - start)) / CLOCKS_PER_SEC;


// Prøver igen

int m=10;
gsl_matrix* d=gsl_matrix_alloc(m,m);
gsl_vector* U=gsl_vector_alloc(m);
int P=1; //Overvej lige noget index her

for(int i=0;i<m;i++){
	gsl_matrix_set(d,i,i,RND);
	gsl_vector_set(U,i,RND);}


double b[]={1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5};

start = clock();
for(int i=0;i<10;i++){
	gsl_vector_set(x,0,b[i]);
	secular_equation_solve(d,U,P,x,0.0001);
	vector_print("res er",x);
	}
end = clock();
time_10= ((double) (end - start)) / CLOCKS_PER_SEC;



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

return 0;
}
