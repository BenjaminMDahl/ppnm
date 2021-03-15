CFLAGS = -Wall -O1 -std=gnu11
LDLIBS = -lm

exp.pdf: exp.txt
	pdflatex $<
	pdflatex $<



.PHONEY: clean
clean:
	$(RM) exp.pdf exp.log exp.aux
