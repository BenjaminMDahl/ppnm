#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <gsl/gsl_integration.h>

double integrate(double f(double), double a, double b, double delta, double eps);


double kvrod(double x);
double kvrod_2(double x);
double kvrod_2_cc(double x);
double inv_kvrod(double x);
double inv_kvrod_cc(double x);
double ln_kvrod(double x);
double ln_kvrod_cc(double x);

int main(){

double Q1=integrate(kvrod,0,1,0.0001,0.0001);

printf("Mit numeriske integral af sqrt(x) fra 0 til 1 giver\n %10g\n",Q1);
double ana=(double) 2/3;
printf("Det analytiske svar er\n\n %10g\n",ana);

double Q2=integrate(kvrod_2,0,1,0.00001,0.00001);
double Q2_CC=integrate(kvrod_2_cc,0,M_PI,0.00001,0.00001);
//double Q2_GSL=


printf("Mit numeriske integral af 4*sqrt(1-x^2) fra 0 til 1 giver\n %30g\n",Q2);
printf("Når man bruger CC transformation på det oversteånde får man \n %30g \n",Q2_CC);
//printf("Fra GSL får man\n %30g \n",Q2_GSL);
printf("Det analytiske svar er\n %10g\n\n",M_PI);


double Q3=integrate(inv_kvrod,0,1,0.001,0.001);
double Q3_CC=integrate(inv_kvrod_cc,0,M_PI,0.001,0.001);

printf("int^1_0(1/sqrt(x))dx= \n %10g\n",Q3);
printf("Når man bruger CC transformation på det oversteånde får man \n %10g \n",Q3_CC);
printf("Det analytiske svar er\n %10g\n\n",2.0);

double Q4=integrate(ln_kvrod,0,1,0.001,0.001);
double Q4_CC=integrate(ln_kvrod_cc,0,M_PI,0.001,0.001);

// For at få løst opgave B skal du have trukket antal integrationer ud af din kode, den har vidst nok tallene, men du får dem ikke printet til sammenligning)
// og så mangler du at lave GSL integral.

printf("int^1_0(log(x)/sqrt(x))dx= \n %10g\n",Q4);
printf("Når man bruger CC transformation på det oversteånde får man \n %10g \n",Q4_CC);
printf("Det analytiske svar er\n %10g\n\n",-4.0);

return 0;
}
