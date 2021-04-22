#include"komplex.h"
#include"stdio.h"

int main(){
	komplex a = {1,2}, b = {3,4};

	printf("tester komplex_add:\n");
	komplex r = komplex_add(a,b);
	komplex R = {4,6};
	komplex_print("a=",a);
	komplex_print("b=",b);
	komplex_print("a+b burde være   = ", R);
	komplex_print("a+b vi får = ", r);

	printf("\ntester komplex_sub(samme a og b som før):\n");
	komplex s = komplex_sub(b,a);
	komplex S = {2,2};
	komplex_print("b-a burde være   = ", S);
	komplex_print("b-a vi får = ", s);

return 0;
}
