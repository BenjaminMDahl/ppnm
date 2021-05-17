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
int steps,stepsR,stepsH,stepsB;
gsl_vector* x=gsl_vector_calloc(2);

steps=qnewton(fun,x,0.0005);

vector_print("Minimum for f(x,y)=x^2-x+y^2-y",x);
printf("På så mange skridt = %i\n\n",steps);


//Rosenbrock
gsl_vector* R=gsl_vector_calloc(2);
gsl_vector_set(R,0,-0.9);
gsl_vector_set(R,1,0.9);

stepsR=qnewton(Rosenbrock,R,0);
vector_print("Det minimum jeg finder for Rosenbrock's valley function er:",R);
printf("Det blev gjort på %i skridt\n\n",stepsR);
//Himmelblau
gsl_vector* H=gsl_vector_calloc(2);
gsl_vector_set(H,0,3.0);gsl_vector_set(H,0,0);
stepsH=qnewton(Himmelblau,H,0.0001);
vector_print("Det minimum jeg finder for Himmelblau's function, er:",H);
printf("Det blev gjort på %i skridt\n\n",stepsH);

//Opgave B//
printf("OPGAVE B Berit-Wigner\n\n");
//Berit-Wigner
gsl_vector* B=gsl_vector_calloc(3);
gsl_vector_set(B,0,122);gsl_vector_set(B,1,0.02);gsl_vector_set(B,2,0.001);
vector_print("Vi prøver nu at fitte til Berit Wigner, hvor vores start gæt er",B);
stepsB=qnewton(deviation_function,B,0.00001);
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
