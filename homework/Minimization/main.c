#include<math.h>
#include<stdio.h>
#include<gsl/gsl_vector.h>

void vector_print(char s[], gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i< v->size ;i++)printf("%10g \n",gsl_vector_get(v,i));
	printf("\n");
}

int qnewton(double F(gsl_vector*), gsl_vector* x, double eps);


// Test funktioner fra testfunctions.c
double fun(gsl_vector* x_vec){
	double x=gsl_vector_get(x_vec,0);
	double y=gsl_vector_get(x_vec,1);
	double value=(x*x-x+y*y-y);
return value;
}

double Rosenbrock(gsl_vector*);
double Himmelblau(gsl_vector*);


int main(){
int steps,stepsR,stepsH;
gsl_vector* x=gsl_vector_calloc(2);

steps=qnewton(fun,x,0.0005);

vector_print("Minimum for f(x,y)=x^2-x+y^2-y",x);
printf("På så mange skridt = %i\n\n",steps);


//Rosenbrock
gsl_vector* R=gsl_vector_calloc(2);
gsl_vector_set(R,0,-0.9);
gsl_vector_set(R,1,0.9);

stepsR=qnewton(Rosenbrock,R,0.0001);
vector_print("Det minimum jeg finder for Rosenbrock's valley function er:",R);
printf("Det blev gjort på %i skridt\n\n",stepsR);
//Himmelblau
gsl_vector* H=gsl_vector_calloc(2);
gsl_vector_set(H,0,3.0);gsl_vector_set(H,1,1.9);
stepsH=qnewton(Himmelblau,H,0.0001);
vector_print("Det minimum jeg finder for Himmelblau's function, er:",H);
printf("Det blev gjort på %i skridt\n\n",stepsH);

return 0;
}
