CFFLAGS = -std=c11 -Wall -Wextra -g
LDGLAGS = -lm -lc -lgmp
CC = clang

.PHONY: clean

kem: kem.c
	$(CC) $(CFFLAGS) $(LDGLAGS) -o kem kem.c

clean:
	rm -rf *.o kem
