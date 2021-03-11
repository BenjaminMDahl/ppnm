#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<omp.h>


int main(){
	int n=10000;

	double x[n], y[n];
#pragma omp parallel sections
	{
	#pragma omp section
		{
		for (int i=0; i <n; i++) x[i]= (double)rand() / (double)RAND_MAX;
		}
	#pragma omp section
		{
		for (int i=0; i <n; i++) y[i]= (double)rand() / (double)RAND_MAX;
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
	fprintf(Pi,"Pi=%g\n",pi);
	fclose(Pi);


return 0;
}
