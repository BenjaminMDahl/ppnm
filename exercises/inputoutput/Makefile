CFLAGS = -O -Wall -std=gnu11
LDLIBS = -lm

.PHONEY:default

default: opg1.out.txt opg2.out.txt opg3.out.txt
	cat  opg*

opg1.out.txt:cmdline
	./cmdline 1 2 3 4 5 > $@

opg2.out.txt:stdin
	echo 5 6 7 8 9 10 | ./stdin > $@

opg3.out.txt:readfile
	./readfile til3.txt > $@

.PHONEY:clean
clean:
	$(RM) opg* input* out* stdin cmdline
