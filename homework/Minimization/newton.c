#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>


static const double delta=sqrt(2.22e-16);

//Inspireret af eksemplet på Blackboard
void numeric_gradient(double F(gsl_vector*), gsl_vector* x, gsl_vector* grad){
	double fx=F(x);
	for(int i=0;i<x->size;i++){
		double dx,xi=gsl_vector_get(x,i);
		if(fabs(xi)<sqrt(delta)) dx=delta; //finde en infitisimal størrelse baseret på xi
		else dx=fabs(xi)*delta;
		gsl_vector_set(x,i,xi+dx);
		gsl_vector_set(grad,i,(F(x)-fx)/dx);	//Her dannes gradient vektoreren
		gsl_vector_set(x,i,xi);
	}
}



int qnewton(
	double F(gsl_vector* x),gsl_vector* x, double eps){

	int steps=0, n=x->size;

	// Vi gør plads til diverse GSL_vektorer vi skal bruge
	gsl_matrix* B=gsl_matrix_alloc(n,n);
	gsl_vector* grad=gsl_vector_alloc(n);
	gsl_vector* step=gsl_vector_alloc(n);
	gsl_vector* z=gsl_vector_alloc(n);
	gsl_vector* gz=gsl_vector_alloc(n);
	gsl_vector* y=gsl_vector_alloc(n);
	gsl_vector* u=gsl_vector_alloc(n);

	// Som udgangs punkt starter vores B matrice som en identites marice
	gsl_matrix_set_identity(B);
	numeric_gradient(F,x,grad);
	double fx=F(x),fz;
	while(steps<2000){
	steps++;
	gsl_blas_dgemv(CblasNoTrans,-1,B,grad,0,step);  //Svarer til ligning 6 i kapitlet om minimization and optimization

	//Tjekker om enden stepet eller gradienten er mindre end vores præsition.
	if(gsl_blas_dnrm2(step)<delta*gsl_blas_dnrm2(x))break;
	if(gsl_blas_dnrm2(grad)<eps)break;

	//Skridt størrelsen afgøres
	double lambda=1;
	while(1){
		gsl_vector_memcpy(z,x);
		gsl_vector_add(z,step);
		fz=F(z);
//		double sTg; gsl_blas_ddot(step,grad,&sTg);
		if(fz<fx)break;
		if(lambda<delta){
			gsl_matrix_set_identity(B);
			break;}
		lambda/=2;
		gsl_vector_scale(step,lambda);}

	//Vi danner nu størrelserne fra ligning (11) og (12) fra kapitlet
	numeric_gradient(F,z,gz);
	gsl_vector_memcpy(y,gz);
	gsl_blas_daxpy(-1,grad,y); 		// y=grad(x+s)-grad(x)
	gsl_vector_memcpy(u,step);gsl_blas_dgemv(CblasNoTrans,-1,B,y,1,u); // u=step-By

	//Nu laver vi Broyden update dog kun hvis ændringen er af en hvis størrelse
	double sTy;//step(transposed)*y
	gsl_blas_ddot(step,y,&sTy);
	if(fabs(sTy)>1e-12){
		gsl_vector_scale(u,1/sTy); //c=u/sTy
		gsl_blas_dger(1,u,step,B); //B=cT*s+B
		}
	gsl_vector_memcpy(x,z);
	gsl_vector_memcpy(grad,gz);
	fx=fz;
	}

	gsl_matrix_free(B);gsl_vector_free(grad);gsl_vector_free(step);gsl_vector_free(z);gsl_vector_free(gz);gsl_vector_free(y);gsl_vector_free(u);

return steps;
}
