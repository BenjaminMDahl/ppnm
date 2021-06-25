#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include"Splines.h"


int main(){
	//Rå data
	int n=5;
	double x[5]={1,2,3,4,5};
	double y1[5]={1,1,1,1,1};
	double yi[5]={1,2,3,4,5};
	double yi2[5]={1,4,9,16,25};


	//Linspace
	int i, N=79; double z[N-1];
	for(i=1; i<N;i++){
	z[i-1]=(double)(5*i+100)/(101);
	}




	// y=1
	double f_z_qua_1[N-1], F_z_qua_1[N-1], df_dz_qua_1[N-1];
	qspline* q_spline_y1=qspline_alloc(n,x,y1);
	printf("#index 0: quadratic spline data for y=1(x f(x) df_dx F(x))\n");
	for(i=0;i<N-1;i++){
		f_z_qua_1[i]=qspline_eval(q_spline_y1,z[i]);
		F_z_qua_1[i]=qua_integral(q_spline_y1,z[i]);
		df_dz_qua_1[i]=qua_diff(q_spline_y1,z[i]);
		printf("%g %g %g %g \n",z[i],f_z_qua_1[i],df_dz_qua_1[i],F_z_qua_1[i]);
	}
	printf("\n \n");
	qspline_free(q_spline_y1);

	// y=i
	double f_z_qua_i[N-1], F_z_qua_i[N-1], df_dz_qua_i[N-1];
	qspline* q_spline_yi=qspline_alloc(n,x,yi);
	printf("#index 1: quadratic spline data for y=i(x f(x) df_dx F(x))\n");
	for(i=0;i<N-1;i++){
		f_z_qua_i[i]=qspline_eval(q_spline_yi,z[i]);
		F_z_qua_i[i]=qua_integral(q_spline_yi,z[i]);
		df_dz_qua_i[i]=qua_diff(q_spline_yi,z[i]);
		printf("%g %g %g %g \n",z[i],f_z_qua_i[i],df_dz_qua_i[i],F_z_qua_i[i]);
	}
	printf("\n \n");
	qspline_free(q_spline_yi);


	// y=i^2
	double f_z_qua_i2[N-1], F_z_qua_i2[N-1], df_dz_qua_i2[N-1];
	qspline* q_spline_yi2=qspline_alloc(n,x,yi2);
	printf("#index 0: quadratic spline data for y=i^2(x f(x) df_dx F(x))\n");
	for(i=0;i<N-1;i++){
		f_z_qua_i2[i]=qspline_eval(q_spline_yi2,z[i]);
		F_z_qua_i2[i]=qua_integral(q_spline_yi2,z[i]);
		df_dz_qua_i2[i]=qua_diff(q_spline_yi2,z[i]);
		printf("%g %g %g %g \n",z[i],f_z_qua_i2[i],df_dz_qua_i2[i],F_z_qua_i2[i]);
	}
	printf("\n \n");
	qspline_free(q_spline_yi2);

return 0;
}
