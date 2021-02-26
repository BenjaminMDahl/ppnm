#include<stdio.h>
#include<math.h>
#include<gsl/gsl_sf_gamma.h>

#include"mygamma.h"


int main(){
	double xmin=1,xmax=6;
	for(double x=xmin;x<=xmax;x+=1.0/8){
		printf("%10g %10g %10g %10g\n",x,tgamma(x),gsl_sf_gamma(x),mygamma(x));
		}

return 0;
}
