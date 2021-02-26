#include<stdio.h>
#include<math.h>
#include<gsl/gsl_sf_erf.h>

#include"mygamma.h"


int main(){
	double xmin=1,xmax=6;
	for(double x=xmin;x<=xmax;x+=1.0/8){
		int z= x;
		printf("%10g %10g %i %10g\n",x,tgamma(x),gsl_sf_gamma(z),mygamma(x));
		}

return 0;
}
