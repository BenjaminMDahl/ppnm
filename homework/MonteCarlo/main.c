#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>

double complex plainMC(int dim,double f(int dim,double* x),double* a,double* b,int N); // Der bruges complex for at få retuneret to værdier både resultatet af integralet(real delen) og fejlen(imaginær delen)

double complex quasiMC(int dim, double f(int,double*), double* a,double* b,int N);

//Test funktioner som skal integreres//
double fun(int dim,double* x){
double value=0;
for(int i=0;i<dim;i++)value=value+x[i];
return value;
}

double fun_e(int dim,double* x){
double value=exp((x[0]+x[1])*(x[0]+x[1]));
return value;
}

double fun_opgaveformulering(int dim, double* x){
double value=(double) 1/(1-cos(x[0])*cos(x[1])*cos(x[2]))*1/(M_PI*M_PI*M_PI);
return value;
}




int main(){

//OPGAVE A//
printf("Opgave A\n\n");

//Test f(x,y,z,u)=x+y+z+u
int dim=4;
double a[]={0,0,0,0};
double b[]={1,1,1,1};

double complex test=plainMC(dim,fun,a,b,100000);
printf("Testen er blot x+y+z+u, fra 0 til 1 langs alle kanter, så svaret bør være 2,0\n");
printf("Test forsøg resultat=%f\n",creal(test));
printf("Test forsøg fejl=%f\n",cimag(test));
printf("Det passer fint med det analytiske hvor vi brugte 100000 punkter\n\n");

//Test 2 med lidt sværere f(x,y)=e^(x+y)^2
int dim2=2;
double a2[]={0,0};
double b2[]={1,1};

double complex test2=plainMC(dim2,fun_e,a2,b2,1000000);

printf("Testen er nu exp(x+y)^2, fra 0 til 1 langs begge kanter, så svaret bør være 4,8992\n");
printf("Test forsøg resultat=%f\n",creal(test2));
printf("Test forsøg fejl=%f\n",cimag(test2));
printf("Det passer ikke lige så godt, men stadig godt, da dette er et svært integral. Her var brugt 1000000 punkter. (et 0 mere end før)\n\n");


//Test 3 funktionen fra opgave formuleringen
int dim3=3;
double a3[]={0,0,0};
double b3[]={M_PI,M_PI,M_PI};

double complex test3=plainMC(dim3,fun_opgaveformulering,a3,b3,10000000);

printf("Testen er nu fra opgaveformuleringen, fra 0 til pi langs alle 3 kanter, så svaret bør være omkring 1,3932039297\n");
printf("Test forsøg resultat=%f\n",creal(test3));
printf("Test forsøg fejl=%f\n",cimag(test3));
printf("Det passer ikke lige så godt, men det forventede er indne for fejlen. Her var brugt 10000000 punkter. (et 0 mere end før)\n\n");

//OPGAVE B//
printf("Opgave B\n\n");
//Vi prøver nu med samme antal N og udregne test funktion 3, nemlig den fra opgave formuleringen

double complex quasi_test=quasiMC(dim3,fun_opgaveformulering,a3,b3,1000000);

printf("Vi bruger nu vores quasi monte carlo simulering som bygger på Van Der Corput og Haltons metoder,\n");
printf("til at udregne det samme integral som det sidste i opgave A, hvor der her udnyttes N=1000000.\n");
printf("Her får vi følgende resultat=%f\n",creal(quasi_test));
printf("med følgende fejl=%f\n",cimag(quasi_test));
printf("Det ses at fejlen er blevet meget mindre end før\n\n");

printf("som en sidste ting i denne besvarelse, har vi prøvet at gøre som i dit eksempel og løst for x^2+y^2<R^2=1\n Resultatet af dette kan findes i comparison.png\n");
printf("Dette er for at sammenligne de to metoder, og se hvordan fejlene fra dem, aftager som funktion af N. Det ses at Quasi klarer sig noget bedre end de plain Montecarlo\n");

return 0;
}


