CC = gcc
CFLAGS = -O -fwrapv -std=gnu11
LDLIBS = -lm

default: out.txt
	cat out.txt

out.txt: main
	./main > out.txt
main: main.o
	$(CC) -o main main.o $(LDLIBS)

main.o: main.c
	$(CC) $(CFLAGS) -c main.c

clean:
	$(RM) main main.o out.txt
