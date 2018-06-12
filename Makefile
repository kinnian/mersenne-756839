CFFLAGS = -std=c11 -Wall -Wextra -g
LDGLAGS = -lm -lc
CC = clang

.PHONY: clean

all: kem

clean:
	rm -rf *.o kem
