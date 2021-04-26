#include<stdio.h>
#include<time.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>


void test_f1(gsl_vector* x,gsl_vector* f){
double x0=gsl_vector_get(x,0);
gsl_vector_set(f,0,x0*x0);

}

void test_f2(gsl_vector* x_vec,gsl_vector* f){
double x=gsl_vector_get(x_vec,0);
double y=gsl_vector_get(x_vec,1);
double z=gsl_vector_get(x_vec,2);

gsl_vector_set(f,0,x+z);
gsl_vector_set(f,1,y+z);
gsl_vector_set(f,2,z+x+y);
}

void test_f3(gsl_vector* x_vec,gsl_vector* f){
double x=gsl_vector_get(x_vec,0);
double y=gsl_vector_get(x_vec,1);


gsl_vector_set(f,0,(2-2*x+400*x*(y-x*x))*(-1));
gsl_vector_set(f,1,200*y-200*x*x);
}
