#include<stdio.h>
#include<time.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_sort_vector.h>
#include<gsl/gsl_eigen.h>
#include<assert.h>
#define RND (double)rand()/RAND_MAX

void make_rnd_vector(gsl_vector* v);
void vector_print(char s[], gsl_vector* v);
void  matrix_print(char s[], gsl_matrix* A);
void secular_equation_default(gsl_vector* D,gsl_vector* u, int p ,gsl_vector* x, double eps,int m);
void secular_equation_GSL(gsl_vector* D,gsl_vector* u, int p ,gsl_vector* x, double eps);
void secular_equation_guess(gsl_vector* D,gsl_vector* u, int p ,gsl_vector* x, double eps);


int main(){
printf("DEL 1: TESTNING AF FUNDNE EGENVÆRDIER\n");
printf("Først testes vores funktioner på to matricer, hvor resultaterne holdes op i mod GSL's gsl_eigen_symm funktionen\n\n");

	// Test 1 med n=5 og ingen ens indgange i D
int n=5;
gsl_eigen_symm_workspace* w= gsl_eigen_symm_alloc(n);
gsl_vector* D=gsl_vector_alloc(n);
gsl_matrix* G=gsl_matrix_alloc(n,n);
gsl_vector* u=gsl_vector_alloc(n);
gsl_vector* x=gsl_vector_alloc(n);
gsl_vector* g=gsl_vector_alloc(n);
int p=2;
double tal[]={1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10};
double guess[]={-20,1.1,2.2,4.4,20};
for(int i=0;i<n;i++){
	gsl_vector_set(D,i,tal[i]);
	gsl_matrix_set(G,i,i,tal[i]);
	gsl_matrix_set(G,p,i,tal[i]);
	gsl_matrix_set(G,i,p,tal[i]);
	gsl_vector_set(u,i,tal[i]);
	gsl_vector_set(x,i,guess[i]);}
gsl_matrix_set(G,p,p,tal[p]+2*tal[p]);
matrix_print("Vi starter med at kigge på en simpel matrice(n=5 og ingen ens indgange) som denne",G);
gsl_eigen_symm(G,g,w);
gsl_sort_vector(g);
vector_print("Vi har startet med at bruge gsl_eigen_symm på denne, for at se hvad vi forventer at finde",g);

vector_print("Med følgende start gæt",x);
secular_equation_guess(D,u,p,x,0.00001);
vector_print("Får vi fra guess",x);
printf("Hvilket er meget tilfredsstillende sammenlignet med gsl resultatet..\n\n");

secular_equation_default(D,u,p,x,0.00001,0);
vector_print("Med vores default routine med m=0 får vi",x);
printf("Hvilket også er meget tilfredsstillende.\n\n");

secular_equation_GSL(D,u,p,x,0.00001);
vector_print("Og til sidst får vi igen det ønskede med vores GSL rutine",x);
printf("Alle var med en eps=0.00001\n\n");

gsl_eigen_symm_free(w);gsl_vector_free(D);gsl_matrix_free(G);gsl_vector_free(u);gsl_vector_free(x);gsl_vector_free(g);

// Test 2 med n=20 og nu med ens indgange

n=20;
gsl_eigen_symm_workspace* wb= gsl_eigen_symm_alloc(n);
gsl_vector* Db=gsl_vector_alloc(n);
gsl_matrix* Gb=gsl_matrix_alloc(n,n);
gsl_vector* ub=gsl_vector_alloc(n);
gsl_vector* xb=gsl_vector_alloc(n);
gsl_vector* gb=gsl_vector_alloc(n);
p=5;
double guessb[]={-20,1.001,1.1,2.001,2.2,3.001,3.3,4.001,4.25,5.001,5.3,6.2,7.001,7.45,8.001,8.5,9.001,9.5,10.001,45};
for(int i=0;i<n;i++){
	gsl_vector_set(Db,i,tal[i]);
	gsl_matrix_set(Gb,i,i,tal[i]);
	gsl_matrix_set(Gb,p,i,tal[i]);
	gsl_matrix_set(Gb,i,p,tal[i]);
	gsl_vector_set(ub,i,tal[i]);
	gsl_vector_set(xb,i,guessb[i]);}
gsl_matrix_set(Gb,p,p,tal[p]+2*tal[p]);
matrix_print("Vi prøver nu med n=20 og hvor de fleste indgange er der to gange.",Gb);
gsl_eigen_symm(Gb,gb,wb);
gsl_sort_vector(gb);
vector_print("Vi starter igen med at bruge gsl_eigen_symm på denne, for at se hvad vi forventer at finde",gb);

vector_print("Med følgende start gæt",xb);
secular_equation_guess(Db,ub,p,xb,0.00001);
gsl_sort_vector(xb);
vector_print("Får vi fra guess",xb);
printf("Igen et meget tilfredsstillende resultat, men det er også et godt gæt, og man ser de steder, hvor vi netop har dobbelt indgange,");
printf("\at denne metode ikke retunere svaret lige så præcist som de andre metoder\n\n");

secular_equation_default(Db,ub,p,xb,0.00001,0);
vector_print("Med vores default routine med m=0 får vi",xb);
printf("Her misser vi en enkelt egenværdi, hvilket ikke er overraskende, da den er meget tæt på grænsen (2.05976 i forhold til 2).\n");
secular_equation_default(Db,ub,p,xb,0.00001,1);
vector_print("Prøve vi dog igen med m=1, hvilket svarer til at gætte igen ved den missede værdi, ser vi, at vi finder alle egenværdier",xb);

secular_equation_GSL(Db,ub,p,xb,0.00001);
vector_print("Og til sidst får vi igen det ønskede med vores GSL rutine",xb);
printf("Alle var igen med en eps=0.00001\n\n");
printf("Et plot af første test kan ses på test.png med resultaterne fra default rutinen.\n\n\n");

gsl_eigen_symm_free(wb);gsl_vector_free(Db);gsl_matrix_free(Gb);gsl_vector_free(ub);gsl_vector_free(xb);gsl_vector_free(gb);

//Tidstagning

printf("DEL 2: TIDSTAGNING SOM TEST AF ANTAL OPERATIONER\n");
printf("Vi prøver her at laver tilfældige dannes D og u'er af hensholdvis n=1000 og n=10000, og tager tid på vores rutine.\n");
printf("Dette er for at overbevise os om, at vi bruger O(n^2) operationer, da en stigning i n på 10 burde medføre en stigning i tid på 100.\n");
printf("Vi undersøger kun vores default metode, med m=0 da denne også svarer til vores guess metode, dog måske lidt langsommere pga. sorteringen af D's indgange\n\n");

clock_t start, stop;
double t1,t10;
n=1000;

gsl_vector* D1=gsl_vector_alloc(n);
gsl_vector* u1=gsl_vector_alloc(n);
gsl_vector* x1=gsl_vector_alloc(n);
make_rnd_vector(D1);make_rnd_vector(u1);make_rnd_vector(x1);

start = clock();
secular_equation_default(D1,u1,5,x1,0.00001,0);
stop = clock();
t1 = ((double) (stop - start)) / CLOCKS_PER_SEC;


gsl_vector* D10=gsl_vector_alloc(n*10);
gsl_vector* u10=gsl_vector_alloc(n*10);
gsl_vector* x10=gsl_vector_alloc(n*10);
make_rnd_vector(D10);make_rnd_vector(u10);make_rnd_vector(x10);


start = clock();
secular_equation_default(D10,u10,5,x10,0.0001,0);
stop = clock();
t10 = ((double) (stop - start)) / CLOCKS_PER_SEC;


printf("Det tog for n=%i så lang tid at nå i mål: \n %g s. \n",n,t1);
printf("For n=%i tog det så lang tid: \n %g s. \n",n*10,t10);
printf("De to delt med hinanden giver: \n %g. \n",t10/t1);
printf("Så en stigning på 10 i n gør det hele ca. 100 gange langsommere, hvilket antyder at vi når i mål på O(n^2) operationer\n\n");

gsl_vector_free(D1);gsl_vector_free(u1);gsl_vector_free(x1);
gsl_vector_free(D10);gsl_vector_free(u10);gsl_vector_free(x10);

// m-Afhængighed
printf("DEL 3: m-AFHÆNGIGHED FOR DEFAULT METODEN\n");
printf("Som en sidste ting har vi undersøgt effektiviteten af defaultmetoden, både i forhold til følsomheden over for\n");
printf("størrelsen af intervaller, og som funktion af m. Resultatet kan ses på figur effektivitet.png, hvor der er dannet\n");
printf("tre tilfældige matricer mellem større og større intervaller, og hvor m stiger løbende. Det ses at effektiviteten stiger\n");
printf("og at metoden fungerer bedst for større intervaller. Matricerne der testes på er af størrelse n=1000\n");



return 0;
}
