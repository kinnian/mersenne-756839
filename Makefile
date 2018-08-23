CFLAGS = -std=c11 -Wall -Wextra -g
LDFLAGS = -lm -lc -lgmp -lcrypto
CC = clang

.PHONY: clean

kem: kem.o rng.o
#	$(CC) $(CFFLAGS) $(LDGLAGS) -o kem kem.c

clean:
	rm -rf *.o kem
