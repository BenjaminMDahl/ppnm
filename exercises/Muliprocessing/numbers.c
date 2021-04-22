#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<omp.h>
#include<time.h>

int main(){
	int n=5000;

	double x[n], y[n];
	unsigned int seed = time(NULL);
#pragma omp parallel sections
	{
	#pragma omp section
		{
		for (int i=0; i <n; i++) x[i]= (double)rand_r(&seed) / (double)RAND_MAX;
		}
	#pragma omp section
		{
		for (int i=0; i <n; i++) y[i]= (double)rand_r(&seed) / (double)RAND_MAX;
		}
	}

	int N_in=0;
	for(int i=0; i<n; i++){
		if(sqrt(x[i]*x[i]+y[i]*y[i])<=1){
			 N_in=N_in+1;
			}
				}
	double pi=4*(double)N_in/(double)n;

	for(int i=0; i<n; i++){
		printf("%g %g\n",x[i],y[i]);
		}

	FILE * Pi;
	Pi = fopen ("file.txt", "w+");
	fprintf(Pi,"Vi får på denne måde Pi=%g\n",pi);
	fclose(Pi);


return 0;
}
