#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<pthread.h>

int main(){
	int n=10000;
	double cirkelx[n];
	double cirkely[n];
	double x[n];
	double y[n];

	for (int i=0; i <n; i++){
		x[i]=cos(rand());
	}
	for (int i=0; i <n; i++){
		y[i]=cos(rand());
	}
	for (int i=0; i <n; i++){
		cirkelx[i]=cos(3.14*2*i/n);
		cirkely[i]=sin(3.14*2*i/n);
	}

	for(int i=0; i<n; i++){
		printf("%g %g %g %g\n",x[i],y[i],cirkelx[i],cirkely[i]);
		}

return 0;
}
