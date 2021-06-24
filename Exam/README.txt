Project 11:  Symmetric Row/column update of a size n smmetric eigenvalue problem
Student: Benjamin M. Dahl 2017089999

INTRO
I dette projekt skal af matricen A's egenværdier bestemmes, hvor A vides at være på den opdaterede form.

	A = D + e(p) u^t + u e(p)^t

Her er:
	- D en dirgonal matric
	- e en enhedsvektor i retning p
	- u en søjlevektor

Bogen nævner at en matrice på en sådan form har følgende karakteristik polynomium

	0=-(d_p-x)+sum[u_i^2/(d_i-x)]	for i={0,1,2..n}/p 

Så problemet er altså reduceret til et rodfindingsproblem, da rødder til overstående vil være vores ønskede egenværdier. 
Det forventes også at rødderne kan findes i følgende intervaller:

	-2*|u^t u| <=x_i<= d_i
	x_n

Besvarelse
Ovenstående problemstilling har jeg løst på tre forskellige måder som har resulteret i tre funktioner.

	Eigenfinderguess
	Eigenfinderdefault
	EigenfinderGSL

Jeg vil gennemgå alle tre her.

Eigenfinderguess
Denne funktion tager fire argumenter. Tre vektorer og en double. De første to vektorer indeholder indgangene for D og u, og den tredje vektor er et gæt på egenværdierne.
Doublen angiver en præcision for undersøgelsen.
Funktionen starter med at danne det ovenstående polynomium og så giver den den til en anden funktion som hedder newton. Newton er min egen rodfindingsrutine fra homework @ som bygger på Newton-Raphson-metoden.
Den bruger alle startgættene og erstatter dem med de rødder den finder indenfor præcisionen angivet ved doublen.

Fordele og ulemper
Metoden er hurtig, men meget følsom overfor det startgæt man angiver, så man kan risikere at finde den samme rod flere gange, hvis gættet er dårligt. Samtidig løber metoden ind i problemer, hvis der er ens
indgange i matrice D, da dette betyder at en sådan indgang også vil være en egenværdi, men det giver anledning til en diskontinuitet i polynomiet, hvilket min newton-rutine har meget svært ved at håndtere.
Metoden bør derfor kun anvendes, hvis man har god kendskab til sin matrice. Ellers skal man bruge Eigenfinderdefault.

Eigenfinderdefault
Eigenfinderdefault bygger på samme princip som før, men den tager et ekstra argument: m. Den laver sit eget startgæt baseret på midten af intervallerne nævnt i introduktionen. Der er intet der sikrer at 
newton-metoden nødvendigvis finder en rod i intervallet, så derfor tjekker rutinen om roden er i det ønskede interval. Hvis ikke, returnerer den nan. Argumentet m er et præcisionsargument, som gør at 
i stedet for at returnere nan, prøver den igen m*2 gange symmetrisk i intervallet. Finder den det stadig ikke, returnerer den stadig nan. Som en ekstra ting tester den også om der er ens indgange i D og
i så fald bruger den ikke newton-metoden, men returnerer bare direkte en sådan indgang som en egenværdi. Metoden kører på O(n^2).

Fordele og ulemper
Metoden kræver ikke kenskab til ens matrice og checker for identiske indgange i D for at undgå problemer med diskontinurteter. Den er dog følsom overfor små intervaller, da newton metoden nemt kommer
til at hoppe ud af intervallerne i en sådan situation, og derfor ikke finder den ønskede egenværdi. Man kan øge præsisonen af denne, men dette går udover hastigheden da den i værste fald så kommer op på
O((1+2m)*n^2) operationer. Desuden sorterer metoden også ens indgange i D, hvilket tager ekstra tid sammenligninget med Eigenfinderdefualt. 

EigenfinderGSL
Som et sidste metode har vi prøver at bruge en rodfindingsrutine fra GSL. Denne bygger dog ikke på Newton-Raphson-metoden, men på en bisection metode. Bisection metoden er det man vil kalde en "root bracketing" metode, som ikke har været pensum i
dette kursus, vi har kigget på "root polishing" metoder som fx Newton-Raphson til at finde rødder, og det er derfor jeg har tilladt mig at bruge en gsl rutine i stedet for at udvikle min egen. Jeg valgte at 
prøve at bruge denne, netop fordi vi har kendskab til hvor rødderne bør være, og derfor burde en sådan metode klare sig bedre.

Fordele og ulemper
Metoden kræver ligesom Eigenfinderdefault ikke kenskab til ens matrice og tjekker også for dobbelt indgange. Metoden er stadig lidt følsom over for små intervaller, eller hvis rødderne er tæt på kanterne,
da enderne af vores intervaler er ved diskontinurteter og derfor ikke kan gå helt ud til de ønskede kanter.

Mappe indhold
I denne mappe findes følgende filer:

	README.txt
	Den vi er i nu
	
	Eigenfinderguess.c Eigenfinderdefault.c EigenfinderGSL.c
	Disse indeholder de overstående metoder.

	homeworkTools.c
	Denne indeholder newton metoden fra homework @ og andre relevante funktioner fra andre Homeworks som newton bruger fx en Gram-Schmiht-solver

	Vectortools.c
	Indeholder diverse funktioner til at printe vektorer/matricer og danne tilfældige vektorer

	main.c
	Står for at bruge og teste de forskellige funktioner og danne outputet til out.txt

	out.txt
	Her findes tests og resultater af de forskellige funktioner og dataene bag png'erne i denne mappe.

	@.png
