#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include"network.h"

double activation_function(double x){return x*exp(-x*x);} // Gaussian wavelet
double function_to_fit(double x){return cos(x/2)*x*x;}
int main(){
	int n=4;
	ann* network=ann_alloc(n,activation_function);

	//Laver tabulated values baseret på function to fit
	double a=-1,b=3;
	int N=25;
	gsl_vector* vx=gsl_vector_alloc(N);
	gsl_vector* vy=gsl_vector_alloc(N);

	for(int i=0;i<N;i++){
		double x=a+(b-a)*i/(N-1);
		double f=function_to_fit(x);
		gsl_vector_set(vx,i,x);
		gsl_vector_set(vy,i,f);
	}

	// Giver start parametre til de forskellige neuroner og træner
	for(int i=0;i<network->n;i++){
		gsl_vector_set(network->params,3*i+0,a+(b-a)*i/(network->n-1)); //sættes til at være uniformt fordelt over det søgte interval
		gsl_vector_set(network->params,3*i+1,1);			//sættes til 1
		gsl_vector_set(network->params,3*i+2,1);			//sættes til 1
	}
	ann_train(network,vx,vy);

	// Skriver data ud til plotning
	printf("#index 1: data points\n");
	for(int i=0;i<vx->size;i++){
		double x=gsl_vector_get(vx,i);
		double f=gsl_vector_get(vy,i);
		printf("%g %g\n",x,f);
	}
	printf("\n\n");

	printf("#index 2: fit from ANN\n");
	for(double z=a;z<=b;z+=0.05){
		double y=ann_response(network,z);
		printf("%g %g\n",z,y);
	}


gsl_vector_free(vx);
gsl_vector_free(vy);
ann_free(network);
return 0;
}
