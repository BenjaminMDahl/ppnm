#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_roots.h>
#include<gsl/gsl_sort_vector.h>
#include<assert.h>


void secular_equation_GSL(gsl_vector* D,gsl_vector* u, int p ,gsl_vector* x, double eps){
	int n=x->size;
	assert(n==D->size);
	assert(n>1); // Vi arbejder kun med matricer ikke skalar

	// Vi starter med at opdateret D, så hvis den p'te indgang i u ikke er 0, updateres D og u sættes lig 0.
	double dstart=gsl_vector_get(D,p);
	double ustart=gsl_vector_get(u,p);
	gsl_vector_set(D,p,dstart+2*ustart);
	gsl_vector_set(u,p,0);

	// Vi skal have en ordnet udgave af D for at kunne bestemme intevallerne vi vil lede efter eigenvalues i.
	gsl_vector* Dcopy=gsl_vector_alloc(D->size);
	gsl_vector_memcpy(Dcopy,D);
	double xp=gsl_vector_get(Dcopy,p);
	gsl_sort_vector(Dcopy);

	//Vi finder hvilket index p  svarer til i den sorterede liste
	int pnew=0;
	while(pnew<=Dcopy->size){
	if(gsl_vector_get(Dcopy,pnew)==xp)break;
	pnew++;}

	// Nu dannes upper og lower bounds
	gsl_vector* low=gsl_vector_alloc(n);
	gsl_vector* high=gsl_vector_alloc(n);
	double a;
	gsl_blas_ddot(u,u,&a);
	for(int i=0;i<pnew;i++){
		gsl_vector_set(low,i+1,gsl_vector_get(Dcopy,i));
		gsl_vector_set(high,i+1,gsl_vector_get(Dcopy,i+1));}
	for(int i=pnew+1;i<n-1;i++){
		gsl_vector_set(low,i,gsl_vector_get(Dcopy,i));
		gsl_vector_set(high,i,gsl_vector_get(Dcopy,i+1));}
	gsl_vector_set(low,0,-2*fabs(a)) ; gsl_vector_set(low,n-1,gsl_vector_get(Dcopy,n-1));
	gsl_vector_set(high,0,gsl_vector_get(Dcopy,0)) ; gsl_vector_set(high,n-1,2*fabs(a));


	double Se(double x,void *nothing){
		double S=x-gsl_vector_get(D,p);
		double P=1;
		for(int k=0;k<p;k++){
			double dk=gsl_vector_get(D,k);
			double uk=gsl_vector_get(u,k);
			P=P*(dk-x);
			S=S+uk*uk/(dk-x);}
		for(int k=p+1;k<u->size;k++){
			double dk=gsl_vector_get(D,k);
			double uk=gsl_vector_get(u,k);
			P=P*(dk-x);
			S=S+uk*uk/(dk-x);}
		return S*P;}

	gsl_function F;
	F.function = &Se;
	const gsl_root_fsolver_type * T = gsl_root_fsolver_bisection;
	gsl_root_fsolver * s = gsl_root_fsolver_alloc(T);

	double r=0; int status;
	for(int i=0;i<n;i++){
		if(gsl_vector_get(low,i)==gsl_vector_get(high,i))r=gsl_vector_get(low,i);
		else{
			double l=gsl_vector_get(low,i)*1.00001;	//Ændre lidt på dette, da disse punkter er præcist et divergenspunkt
			double h=gsl_vector_get(high,i)*0.99998;
			gsl_root_fsolver_set(s,&F,l,h);
			int iter=0,max_iter=100;
			do{
				iter++;
				status = gsl_root_fsolver_iterate (s);
				r = gsl_root_fsolver_root (s);
				double x_lo = gsl_root_fsolver_x_lower (s);
				double x_hi = gsl_root_fsolver_x_upper (s);
				status = gsl_root_test_interval (x_lo, x_hi,0, 0.001);
			}while (status == GSL_CONTINUE && iter < max_iter);}
		gsl_vector_set(x,i,r);}
	gsl_vector_free(Dcopy);
}
