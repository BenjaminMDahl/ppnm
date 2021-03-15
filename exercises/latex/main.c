#include<stdio.h>
#include<math.h>


double ex(double x){
if(x<0)return 1/ex(-x);
if(x>1./8)return pow(ex(x/2),2);
return 1+x*(1+x/2*(1+x/3*(1+x/4*(1+x/5*(1+x/6*(1+x/7*(1+x/8*(1+x/9*(1+x/10)))))))));
}

int main(){
	double x[15]={0.1,0.2,0.3,0,4,0.5,0.7,1,2,3,4,5,7,9,10};
	for(int i=0;i<15;i++){
	printf("%g %g %g\n",x[i],exp(x[i]),ex(x[i]));
	}
return 0;
}
