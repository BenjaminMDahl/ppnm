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


// solver routinen er lavet lidt om siden homeworks lin solve i den forstand at GS-solve selv kalder QR-decomposition på A og R-
// Dette er gjort for at spare en sætning i koden
void GS_solve(gsl_matrix* A, gsl_matrix* R, gsl_vector* b, gsl_vector* x);


void newton(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double eps){

	// gør diverse gsl matricer og vektorer klar
	int n=x->size;double dx=sqrt(eps);
	gsl_matrix* J=gsl_matrix_alloc(n,n);
	gsl_matrix* R=gsl_matrix_alloc(n,n);
	gsl_vector* fx=gsl_vector_alloc(n);
	gsl_vector* dfx=gsl_vector_alloc(n);
	gsl_vector* fy=gsl_vector_alloc(n);
	gsl_vector* xstep=gsl_vector_alloc(n);
	gsl_vector* y=gsl_vector_calloc(n);

	// Herunder dannes jacobianen jævnført ligning (7) i kapitlet om nonlinear equations og min GS_solver kaldes på den
	while(1){
		f(x,fx);
		for(int i=0;i<n;i++){
			double xi=gsl_vector_get(x,i);
			gsl_vector_set(x,i,xi+dx);
			f(x,dfx);
			gsl_vector_sub(dfx,fx);
			gsl_vector_scale(dfx,1/dx);
			for(int j=0;j<n;j++)gsl_matrix_set(J,j,i,gsl_vector_get(dfx,j));
			double dxi=gsl_vector_get(x,i);
			gsl_vector_set(x,i,dxi-dx);}
		gsl_vector_scale(fx,-1);//-f(x)
		GS_solve(J,R,fx,xstep);

		// Her tjekkes om det funde x er tilfredsstillende
		double s=2;
		gsl_vector_memcpy(y,x);
		while(1){
			s/=2;
			gsl_vector_scale(xstep,s);
			gsl_vector_add(y,xstep);
			gsl_vector_scale(xstep,1/s);
			f(y,fy);
//		printf("fy=%g\n",gsl_blas_dnrm2(fy));
//		printf("faktorfx=%g\n",(1-s/2)*gsl_blas_dnrm2(fx));
			if(gsl_blas_dnrm2(fy)<(1-s/2)*gsl_blas_dnrm2(fx))break;
//		printf("s=%g\n",s);
			if(s<0.01)break;
		}
		gsl_vector_memcpy(x,y);gsl_vector_memcpy(fx,fy);
//	printf("step norm=%g\n",gsl_blas_dnrm2(xstep));
		if(gsl_blas_dnrm2(xstep)<dx) break;
//	printf("value norm=%g\n",gsl_blas_dnrm2(fx));
		if(gsl_blas_dnrm2(fx)<eps) break;
	}
	gsl_matrix_free(J);gsl_matrix_free(R);
	gsl_vector_free(fx);gsl_vector_free(dfx);gsl_vector_free(y);gsl_vector_free(fy);
}


void test_f1(gsl_vector* x,gsl_vector* f){
double x0=gsl_vector_get(x,0);
gsl_vector_set(f,0,x0*x0);

}

void test_f2(gsl_vector* x_vec,gsl_vector* f){
double x=gsl_vector_get(x_vec,0);
double y=gsl_vector_get(x_vec,1);
double z=gsl_vector_get(x_vec,2);

gsl_vector_set(f,0,x+z);
gsl_vector_set(f,1,y+z);
gsl_vector_set(f,2,z+x+y);
}

void test_f3(gsl_vector* x_vec,gsl_vector* f){
double x=gsl_vector_get(x_vec,0);
double y=gsl_vector_get(x_vec,1);


gsl_vector_set(f,0,(2-2*x+400*x*(y-x*x))*(-1));
gsl_vector_set(f,1,200*y-200*x*x);
}

int main(){
printf("OPGAVE A\n\n");
//1 dimensionel test
gsl_vector* x_test_1dim=gsl_vector_alloc(1);
gsl_vector_set(x_test_1dim,0,1);
vector_print("Vi har startet med en simpel test for f(x)=x^2 og vores start x er ",x_test_1dim);
double eps=0.000001;

newton(test_f1,x_test_1dim,eps);

vector_print("Her bør vi få vores svar til x=0 og vi ser at vi får=",x_test_1dim);
printf("Det er ikke et vildt præcist resultat, men koden gør som den skal. Dette er ved en epsilon=%g.\n\n",eps);

// 3 dimensionel test
gsl_vector* x_test_2dim=gsl_vector_alloc(3);
gsl_vector_set(x_test_2dim,0,-4.9);
gsl_vector_set(x_test_2dim,1,5);
gsl_vector_set(x_test_2dim,2,0.1);
vector_print("Den næste test er hvor f(x,y,z)=(x+z,y+z,z+x+y), og vi starter med (x,y,z)=(-0.9,1.2,0.1)",x_test_2dim);

newton(test_f2,x_test_2dim,eps);

vector_print("Vi får svaret=",x_test_2dim);
printf("Den laver en 0 vektor, hvilket er en passende løsning\n\n");


// Testen fra opgave formuleringen
printf("Som en sidste del af opgave A prøver vi vores kode på gradienten af funktionen fra opgave formuleringen som er:\n");
printf("gradient(f(x,y))=(-2(1-x)-400x(y-x^2) , 200(y-x^2))");

gsl_vector* x_op=gsl_vector_alloc(2);
gsl_vector_set(x_op,0,0.1);
gsl_vector_set(x_op,1,0.1);


vector_print("Vi starter ud med (x,y)=",x_op);

newton(test_f3,x_op,eps);

vector_print("Vi får svaret=",x_op);

return 0;
}
