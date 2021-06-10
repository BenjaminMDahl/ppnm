#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include"Splines.h"

int main(){

	// Laver data
	int n=32;
	int i;
	double x[n], y[n];
	printf("#index 0: data from sin(x)(x, sin(x),cos(x),-cos(x)+k=1)\n");
	for(i=0;i<n;i++){
		x[i]= (double)i/5;
		y[i]=sin(x[i]);
		printf("%g %g %g %g\n",x[i],y[i],cos(x[i]),-cos(x[i])+1);
	}
	printf("\n \n");


	//Linspace
	int N=195; double z[N];
	for(i=0;i<N;i++){
	z[i]=(double)(10+i)/33;
	}

	// lin data
	double f_z[N], F_z[N];
	printf("#index 1: linear spline data(x f(x) F(x))\n");
	for(i=0;i<N;i++){
		f_z[i]=linterp(n,x,y,z[i]);
		F_z[i]=linterp_integral(n,x,y,z[i]);
		printf("%8g %8g %8g \n",z[i],f_z[i],F_z[i]);
	}
	printf("\n \n");

	// qua data
	double f_z_qua[N], F_z_qua[N], df_dz_qua[N];
	qspline* q_spline=qspline_alloc(n,x,y);
	printf("#index 2: quadratic spline data(x f(x) df_dx F(x))\n");
	for(i=0;i<N;i++){
		f_z_qua[i]=qspline_eval(q_spline,z[i]);
		F_z_qua[i]=qua_integral(q_spline,z[i]);
		df_dz_qua[i]=qua_diff(q_spline,z[i]);
		printf("%8g %8g %8g %8g \n",z[i],f_z_qua[i],df_dz_qua[i],F_z_qua[i]);
	}
	printf("\n \n");
	qspline_free(q_spline);


	// cub data
	double f_z_cub[N], F_z_cub[N], df_dz_cub[N];
	cspline* c_spline=cspline_alloc(n,x,y);
	printf("#index 3: cubic spline data(x f(x) df_dx F(x))\n");
	for(i=0;i<N;i++){
		f_z_cub[i]=cspline_eval(c_spline,z[i]);
		F_z_cub[i]=cub_integral(c_spline,z[i]);
		df_dz_cub[i]=cub_diff(c_spline,z[i]);
		printf("%8g %8g %8g %8g \n",z[i],f_z_cub[i],df_dz_cub[i],F_z_cub[i]);
	}
	cspline_free(c_spline);

return 0;
}
