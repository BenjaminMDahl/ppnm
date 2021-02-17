#include"komplex.h"
#include"stdio.h"

int main(){
	komplex a = {1,2}, b = {3,4};

	printf("testing komplex_add...\n");
	komplex r = komplex_add(a,b);
	komplex R = {4,6};
	komplex_print("a=",a);
	komplex_print("b=",b);
	komplex_print("a+b should   = ", R);
	komplex_print("a+b actually = ", r);

	printf("now testing komplex_sub...\n");
	komplex s = komplex_sub(b,a);
	komplex S = {2,2};
	komplex_print("b-a should   = ", S);
	komplex_print("b-a actually = ", s);

return 0;
}
