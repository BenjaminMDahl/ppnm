#include<stdio.h>
#include<time.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<assert.h>

void vector_print(char s[], gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i< v->size ;i++)printf("%10g \n",gsl_vector_get(v,i));
	printf("\n");
}

void newton(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double eps);


void test_f1(gsl_vector* x,gsl_vector* f);

void test_f2(gsl_vector* x_vec,gsl_vector* f);

void test_f3(gsl_vector* x_vec,gsl_vector* f);

void test_diff(gsl_vector* x,gsl_vector* f);

int main(){

		////Opgave A Newtons metode////
printf("OPGAVE A\n\n");

//1 dimensionel test//
gsl_vector* x_test_1dim=gsl_vector_alloc(1);
gsl_vector_set(x_test_1dim,0,1);
vector_print("Vi har startet med en simpel test for f(x)=x^2 og vores start x er ",x_test_1dim);
double eps=0.000001;

newton(test_f1,x_test_1dim,eps);

vector_print("Her bør vi få vores svar til x=0 og vi ser at vi får=",x_test_1dim);
printf("Det er ikke et vildt præcist resultat, men koden gør som den skal. Dette er ved en epsilon=%g.\n\n",eps);


// 3 dimensionel test//
gsl_vector* x_test_2dim=gsl_vector_alloc(3);
gsl_vector_set(x_test_2dim,0,-4.9);
gsl_vector_set(x_test_2dim,1,5);
gsl_vector_set(x_test_2dim,2,0.1);
vector_print("Den næste test er hvor f(x,y,z)=(x+z,y+z,z+x+y), og vi starter med (x,y,z)=(-0.9,1.2,0.1)",x_test_2dim);

newton(test_f2,x_test_2dim,eps);

vector_print("Vi får svaret=",x_test_2dim);
printf("Den laver en 0 vektor, hvilket er en passende løsning\n\n");


// Testen fra opgave formuleringen//
printf("Som en sidste del af opgave A prøver vi vores kode på gradienten af funktionen fra opgave formuleringen som er:\n");
printf("gradient(f(x,y))=(-2(1-x)-400x(y-x^2) , 200(y-x^2))\n");

//gæt 1
gsl_vector* x_op1=gsl_vector_alloc(2);
gsl_vector_set(x_op1,0,0.1);
gsl_vector_set(x_op1,1,0.1);
//gæt 2
gsl_vector* x_op2=gsl_vector_alloc(2);
gsl_vector_set(x_op2,0,2);
gsl_vector_set(x_op2,1,2);
//gæt 3
gsl_vector* x_op3=gsl_vector_alloc(2);
gsl_vector_set(x_op3,0,-3);
gsl_vector_set(x_op3,1,-1);
//gæt 4
gsl_vector* x_op4=gsl_vector_alloc(2);
gsl_vector_set(x_op4,0,5);
gsl_vector_set(x_op4,1,-5);

newton(test_f3,x_op1,eps);
newton(test_f3,x_op2,eps);
newton(test_f3,x_op3,eps);
newton(test_f3,x_op4,eps);

vector_print("Vi får fra (0.1,0.1) punktet",x_op1);
vector_print("Vi får fra (2,2) punktet",x_op2);
vector_print("Vi får fra (-3,-1) punktet",x_op3);
vector_print("Vi får fra (5,-5) punktet",x_op4);



// Der skal kommenteres her og finds flere punkter


		////Opgave B Kvantemekanik////

printf("Opgave B Kvantemekanik\n\n");
printf("Vi vil nu forsøge at løse et fysisk problem, vi vil nemlig prøve at bestemme den laveste energi tilstand af et brint atom\n");
gsl_vector* x_diff=gsl_vector_alloc(1);
gsl_vector_set(x_diff,0,0);

newton(test_diff,x_diff,eps);

vector_print("Vi får svaret=",x_diff);
printf("og det analytiske svar er -1/2 så vi er meget tilfredse\n");
printf("På figur kvant.png, findes der et plot af vores løsning og den teoretiske f0(r)=r*exp(-r).\n");

return 0;
}
