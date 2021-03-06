#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>


double complex plainMC(int dim,double f(int dim,double* x),double* a,double* b,int N); // Der bruges complex for at få retuneret to værdier både resultatet af integralet(real delen) og fejlen(imaginær delen)
double complex quasiMC(int dim, double f(int,double*), double* a,double* b,int N);


double f(int dim,double* x){
double value=exp((x[0]+x[1])*(x[0]+x[1]));
return value;
	}



int main(){
printf("#index 0: N   PlainMonteCarlo error     QuasiMonteCarlo error\n");
int N=10000;
double a[]={0,0}, b[]={1,1};
for(int n=1000;n<N;n+=50)printf("%i %6g %6g\n",n,cimag(plainMC(2,f,a,b,n)),cimag(quasiMC(2,f,a,b,n)));



return 0;
}
