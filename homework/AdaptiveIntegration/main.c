#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <gsl/gsl_integration.h>

double integrate(double f(double), double a, double b, double delta, double eps, int *record);


double kvrod(double x);
double kvrod_2(double x);
double kvrod_2_cc(double x);
double inv_kvrod(double x);
double inv_kvrod_cc(double x);
double ln_kvrod(double x);
double ln_kvrod_cc(double x);

double gsl_int(double f(double), double a, double b, double delta, double eps);


int main(){

int rec_Q1;

double Q1=integrate(kvrod,0,1,0.0001,0.0001,&rec_Q1);

printf("Mit numeriske integral af sqrt(x) fra 0 til 1 giver\n %10g\n",Q1);
double ana=(double) 2/3;
printf("Det analytiske svar er\n\n %10g\n",ana);

int rec_Q2, rec_Q2_CC;

double Q2=integrate(kvrod_2,0,1,0.00001,0.00001,&rec_Q2);
double Q2_CC=integrate(kvrod_2_cc,0,M_PI,0.00001,0.00001, &rec_Q2_CC);
double Q2_GSL=gsl_int(kvrod_2,0,1,0.00001,0.00001);


printf("Mit numeriske integral af 4*sqrt(1-x^2) fra 0 til 1 giver\n %.17g\n",Q2);
printf("Når man bruger CC transformation på det oversteånde får man \n %.17g \n",Q2_CC);
printf("Fra GSL får man\n %.17g \n",Q2_GSL);
printf("Det analytiske svar er\n %.17g\n\n",M_PI);
printf("Det ses at når de alle er underlagt samme tolerance krav, kan GSL få op til 15 cifre rigtigt,\n mens mine operationer højst giver 8 cifre og kun 7 når man ikke bruger CC transformation.\n\n");
printf("Til sidst ses der hvor mange evalueringer der skulle laves for min metode med og uden CC transformation.\n");
printf("Uden transformation skulle der %i evalueringer til, og med skulle der %i evalueringer til\n",rec_Q2,rec_Q2_CC);


int rec_Q3, rec_Q3_CC;

double Q3=integrate(inv_kvrod,0,1,0.001,0.001,&rec_Q3);
double Q3_CC=integrate(inv_kvrod_cc,0,M_PI,0.001,0.001,&rec_Q3_CC);

printf("int^1_0(1/sqrt(x))dx= \n %10g\n",Q3);
printf("Når man bruger CC transformation på det oversteånde får man \n %10g \n",Q3_CC);
printf("Det analytiske svar er\n %10g\n\n",2.0);
printf("Uden transformation skulle der %i evalueringer til, og med skulle der %i evalueringer til\n",rec_Q3,rec_Q3_CC);

int rec_Q4, rec_Q4_CC;

double Q4=integrate(ln_kvrod,0,1,0.001,0.001,&rec_Q4);
double Q4_CC=integrate(ln_kvrod_cc,0,M_PI,0.001,0.001,&rec_Q4_CC);

// For at få løst opgave B skal du finde ud af om man man få antal integrationer ud af GSL

printf("int^1_0(log(x)/sqrt(x))dx= \n %10g\n",Q4);
printf("Når man bruger CC transformation på det oversteånde får man \n %10g \n",Q4_CC);
printf("Det analytiske svar er\n %10g\n\n",-4.0);
printf("Uden transformation skulle der %i evalueringer til, og med skulle der %i evalueringer til\n",rec_Q4,rec_Q4_CC);

return 0;
}
