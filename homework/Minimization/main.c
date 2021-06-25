#include<math.h>
#include<stdio.h>
#include<gsl/gsl_vector.h>

void vector_print(char s[], gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i< v->size ;i++)printf("%10g \n",gsl_vector_get(v,i));
	printf("\n");
}

int qnewton(double F(gsl_vector*), gsl_vector* x, double eps);


// Test funktioner fra testfunctions.c
double fun(gsl_vector* x_vec){
	double x=gsl_vector_get(x_vec,0);
	double y=gsl_vector_get(x_vec,1);
	double value=(x*x-x+y*y-y);
return value;
}

double Rosenbrock(gsl_vector*);
double Himmelblau(gsl_vector*);
double Breit_Wigner(double,gsl_vector*);
double deviation_function(gsl_vector*);


int main(){
int steps,stepsR,stepsB;
gsl_vector* x=gsl_vector_calloc(2);

steps=qnewton(fun,x,0.0005);

printf("OPGAVE A Minimization)\n\n");
vector_print("Vi starter med at teste vores rountine ved at finde minimum for f(x,y)=x^2-x+y^2-y.\n Vi finder det til:",x);
printf("Det passer godt med det forventede. Det blev fundet på 3 skridt = %i\n\n",steps);


//Rosenbrock
gsl_vector* R=gsl_vector_calloc(2);
gsl_vector_set(R,0,1);
gsl_vector_set(R,1,0);

stepsR=qnewton(Rosenbrock,R,0);
vector_print("Vi prøver med noget lidt svære. Vi undersøger Rosenbrock's valley funktionen.\n Til denne finder vi minimum ved:",R);
printf("Dette er tæt på (1,1), hvilket er det forventede for denne udgave af RosenBbrock's valley funktion(a=1 og b=100).\n Det blev gjort med startgæt (x,y)=(1,0) og på %i skridt\n\n",stepsR);

//Himmelblau
gsl_vector* H1=gsl_vector_calloc(2);
gsl_vector_set(H1,0,3.5);gsl_vector_set(H1,1,1.5);
int stepsH1=qnewton(Himmelblau,H1,0.0001);

gsl_vector* H2=gsl_vector_calloc(2);
gsl_vector_set(H2,0,-3);gsl_vector_set(H2,1,3);
int stepsH2=qnewton(Himmelblau,H2,0.0001);

gsl_vector* H3=gsl_vector_calloc(2);
gsl_vector_set(H3,0,-4);gsl_vector_set(H3,1,-3);
int stepsH3=qnewton(Himmelblau,H3,0.0001);

gsl_vector* H4=gsl_vector_calloc(2);
gsl_vector_set(H4,0,4);gsl_vector_set(H4,1,-2);
int stepsH4=qnewton(Himmelblau,H4,0.0001);

vector_print("Vi prøver nu for Himmelblau's funktion. Her finder vi et minimum med start gæt(3.5,1.5) ved:",H1);
vector_print("Og et med start gæt (-3,3) ved",H2);
vector_print("Og et med start gæt (-4,-3) ved",H3);
vector_print("Og et med start gæt (4,-2) ved",H4);
printf("Dette passer godt med de forventede fire punkter som kan slås op på wikipedia. Hvert minimum blev fundet på henholdsvis %i %i %i og %i skridt\n\n",stepsH1,stepsH2,stepsH3,stepsH4);

//Opgave B//
printf("OPGAVE B Berit-Wigner\n\n");
//Berit-Wigner
gsl_vector* B=gsl_vector_calloc(3);
gsl_vector_set(B,0,122);gsl_vector_set(B,1,0.02);gsl_vector_set(B,2,0.001);
vector_print("Vi prøver nu at fitte til Berit Wigner, hvor vores start gæt er",B);
stepsB=qnewton(deviation_function,B,0.0000001);
vector_print("Det minimum jeg finder for Berit Wigner er følgende(Hvor vi har mass,width,A)",B);
printf("Dette blev opnået på %i skridt\n",stepsB);
printf("Vi har fået oplyst i opgaven at CERN fandt den til 125.3(6), så vi er meget tilfredse med dette resultat\n");
printf("Som en sidste ting har vi prøvet at bruge vores minimerings fit til at lave plottet BreitWigner.png, sammen med CERN data,en og resultatet er meget tilfredsstillende\n\n");

//printer data til plotning
FILE *plot;
plot = fopen("plot.txt","w");
for(double e=95;e<170;e++)fprintf(plot,"%6g %6g\n",e,Breit_Wigner(e,B));
fclose(plot);

return 0;
}
