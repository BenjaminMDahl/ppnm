CC = gcc # the C compiler
CFLAGS = -O -fwrapv -std=gnu11 # options for the C compiler
LDLIBS = -lm # libraries to link

default: out.txt
	cat out.txt

out.txt: eps
	./eps > out.txt

eps: eps.o
	$(CC) -o eps eps.o $(LDLIBS)

eps.o: eps.c
	$(CC) $(CFLAGS) -c eps.c

clean:
	$(RM) eps eps.o out.txt
