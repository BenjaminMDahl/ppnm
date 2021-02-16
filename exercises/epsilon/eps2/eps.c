#include<stdio.h>
#include<float.h>
#include<limits.h>
#include<math.h>

int main(){
int int_max=INT_MAX/3;
float sum_up_float=0.0;
for(int i=1;i<=int_max;i++){
	sum_up_float=sum_up_float+1.0/i;
}

printf("sum up float=%g\n",sum_up_float);

float sum_down_float=0.0;
for(int i=int_max;i>=1;i--){
	sum_down_float=sum_down_float+1.0/i;
}

printf("sum down float=%g\n",sum_down_float);
printf("\n The difference comes from starting by summing small numbers up to big numbers lets the small numbers contribute where the other way around cuts off the small numbers in the end do to how float makes approximations of numbers\n");
printf("\n the sum should not converges since it is the harmonic series however do to how floats a approximated i think it does here\n");


double double_sum_up=0.0;
for(int i=1;i<=int_max;i++){
	double_sum_up=double_sum_up+1.0/i;
}

printf("double sum up=%lg\n",double_sum_up);

float double_sum_down=0.0;
for(int i=int_max;i>=1;i--){
	double_sum_down=double_sum_down+1.0/i;
}

printf("double sum down=%lg\n",double_sum_down);

return 0;
}


