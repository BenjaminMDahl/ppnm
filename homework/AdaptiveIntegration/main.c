#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double integrate(double f(double), double a, double b, double delta, double eps);

double kvrod(double x);
double kvrod_2(double x);

int main(){

double Q1=integrate(kvrod,0,1,0.0001,0.0001);

printf("Mit numeriske integral af sqrt(x) fra 0 til 1 giver\n %10g\n",Q1);
double ana=(double) 2/3;
printf("Det analytiske svar er\n %10g\n",ana);

double Q2=integrate(kvrod_2,0,1,0.0001,0.0001);

printf("Mit numeriske integral af 4*sqrt(1-x^2) fra 0 til 1 giver\n %10g\n",Q2);
printf("Det analytiske svar er\n %10g\n",M_PI);


return 0;
}
