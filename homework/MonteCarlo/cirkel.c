#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>

double complex plainMC(int dim,double f(int dim,double* x),double* a,double* b,int N); // Der bruges complex for at få retuneret to værdier både resultatet af integralet(real delen) og fejlen(imaginær delen)
double complex quasiMC(int dim, double f(int,double*), double* a,double* b,int N);


double fun_opgaveformulering(int dim, double* x){
double value=(double) 1/(1-cos(x[0])*cos(x[1])*cos(x[2]))*1/(M_PI*M_PI*M_PI);
return value;
}


int main(){
printf("#index 0: N   PlainMonteCarlo error     QuasiMonteCarlo error\n");

int N=10000;
double a[]={0,0,0}, b[]={M_PI,M_PI,M_PI};

for(int n=1;n<N;n+=5){
printf("%i %6g %6g\n",n,cimag(plainMC(3,fun_opgaveformulering,a,b,n)),cimag(quasiMC(3,fun_opgaveformulering,a,b,n)));

}

return 0;
}
