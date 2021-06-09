CC = gcc
CFLAGS = -O -fwrapv -std=gnu11
LDLIBS = -lm

default: out.txt
	cat out.txt

out.txt: main
	./main > out.txt

main: main.o fun.o


clean:
	$(RM) main *.o out.txt
