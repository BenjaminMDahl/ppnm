Project 11:  Symmetric Row/column update of a size n smmetric eigenvalue problem
Student: Benjamin M. Dahl 2017089999

	-----INTRODUKTION-----

I dette projekt skal egenværdierne af matricen A's bestemmes, hvor A vides at være på den opdaterede form.

	A = D + e(p) u^t + u e(p)^t

Her er:
	- D en diagonal matrice
	- e en enhedsvektor i retning p
	- u en søjlevektor

Bogen nævner at en matrice på en sådan form har følgende karakteristik polynomium

	0=-(d_p-x)+sum[u_i^2/(d_i-x)]	for i={0,1,2..n}/p 

Så problemet er altså reduceret til et finder rødderne af overstående. 
Det forventes også at rødderne kan findes i følgende intervaller:

	-2*|u^t u| <=x_0<= d_0
	d_i <= x_{i+1} <= d_{i+1}	for i={0,1,2...p-1}
	d_i <= x_{i} <= d_{i+1}		for i={p+1,p,+2...n-1}
	d_n <= x_n <=	2*|u^t u|



	-----BESVARELSE-----
Ovenstående problemstilling har jeg løst på tre forskellige måder som har resulteret i tre funktioner.

	secular_equation_guess
	secular_equation_default
	secular_equation_GSL

Jeg vil gennemgå alle tre her.

secular_equation_guess
Denne funktion tager fire argumenter. Tre vektorer og en double. De første to vektorer indeholder indgangene for D og u, og den tredje vektor er et gæt på egenværdierne.
Doublen angiver en præcision for undersøgelsen.
Funktionen starter med at danne polynomiet fra introduktionen og så giver den polynomiet til en anden funktion som hedder newton. Newton er min egen rodfindingsrutine fra homework 8 som bygger på Newton-Raphson-metoden.
Den bruger alle startgættene og erstatter dem med de rødder den finder indenfor præcisionen angivet ved doublen.

Fordele og ulemper
Metoden er hurtig (O(n^2)), men meget følsom overfor det startgæt man angiver, så man kan risikere at finde den samme rod flere gange, hvis gættet er dårligt. Samtidig løber metoden ind i problemer, hvis der er
ens indgange i matrice D, da dette betyder at en sådan indgang også vil være en egenværdi, men det giver anledning til en diskontinuitet i polynomiet, hvilket min newton-rutine har meget svært ved at håndtere.
Metoden bør derfor kun anvendes, hvis man har godt kendskab til sin matrice. Ellers skal man bruge secular_equation_default.

secular_equation_default
secular_equation_default bygger på samme princip som før, men den tager et ekstra argument: m. Den laver sit eget startgæt baseret på midten af intervallerne nævnt i introduktionen. Der er intet der sikrer at 
newton-metoden bliver inden for det ønskede interval, så derfor tjekker default rutinen om den fundne rod er i det ønskede interval, hvis ikke, returnerer den nan. Argumentet m er et præcisionsargument, som
gør at i stedet for at returnere nan, prøver den igen m*2 gange symmetrisk i intervallet. Finder den det stadig ikke, returnerer den stadig nan. Som en ekstra ting tester den også om der er ens indgange i D og
i så fald bruger den ikke newton-metoden, men returnerer bare direkte en sådan indgang som en egenværdi. Metoden kører på O((1+2m)*n^2).

Fordele og ulemper
Metoden kræver ikke kendskab til ens matrice og tjekker for identiske indgange i D for at undgå problemer med diskontinuiteter. Den er dog følsom overfor små intervaller, da newton metoden nemt kommer
til at hoppe ud af intervallerne i en sådan situation, og derfor ikke finder den ønskede egenværdi. Man kan øge præcisionen af denne, men dette går udover hastigheden da den i værste fald så kommer op på
O((1+2m)*n^2) operationer. Desuden sorterer metoden også ens indgange i D, hvilket tager ekstra tid sammenlignet med secular_equation_guess. 

secular_equation_GSL
Som en sidste metode har vi prøver at bruge en rodfindingsrutine fra GSL. Denne bygger dog ikke på Newton-Raphson-metoden, men på en bisection metode. Bisection metoden er det man vil kalde en "root bracketing" metode,
som ikke har været pensum i dette kursus, vi har kigget på "root polishing" metoder som fx Newton-Raphson til at finde rødder, og det er derfor jeg har tilladt mig at bruge en gsl rutine i stedet for at
udvikle min egen. Jeg valgte at prøve at bruge denne, netop fordi vi har kendskab til hvor rødderne bør være, og derfor burde en root bracketing metode klare sig bedre end en root polising metode.

Fordele og ulemper
Metoden kræver ligesom secular_equation_default ikke kendskab til ens matrice og tjekker også for dobbelt indgange. Metoden er stadig lidt følsom over for små intervaller, eller hvis rødderne er tæt på kanterne,
da enderne af vores intervaller er ved diskontinuiteter og derfor ikke kan gå helt ud til de ønskede kanter.



	-----MAPPENS INDHOLD-----
I denne mappe (Exam) findes følgende filer:

	README.txt
	Den vi er i nu.
	
	Eigenfinderguess.c Eigenfinderdefault.c EigenfinderGSL.c
	Disse indeholder de overstående metoder.

	HomeworkTools.c
	Denne indeholder newton metoden fra homework 8 og andre relevante funktioner fra andre Homeworks som newton bruger fx en Gram-Schmiht-solver

	VectorTools.c
	Indeholder diverse funktioner til at printe vektorer/matricer og danne tilfældige vektorer

	main.c
	Står for at bruge og teste de forskellige funktioner, hvor outputtet ses i out.txt

	plotter.c
	Står for at danne den rå data til plotsne (png filerne) som gemmes i data.txt

	data.txt
	Data til plots (png filerne)

	out.txt
	Her findes tests og resultater af de forskellige funktioner og dataene bag png'erne i denne mappe. 
	Generelt er det i denne min besvarelse komme til udtryk.

	test.png
	Et plot af et karakteristiske polynomium med fundne egenværdier fra secular_equation_default

	effektivitet.png
	Et plot af hvordan secular_equation_default afhænger af m og er følsom overfor størrelsen af intervallerne som dannes ud fra D.
