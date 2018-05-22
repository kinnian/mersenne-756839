#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// Variables globales

unsigned char n = 756839;
unsigned char h = 256;
unsigned char rho = 2048;
unsigned char P = 2^n - 1;
unsigned char K = 32*ceil(n/256);

// TODO: implementer correctement les pounsigned chareurs et listes d'octets / de bits et les concatenations et parse

void key_gen(unsigned char * pk, unsigned char * sk){
	unsigned char F = rand() % n; //TODO: HAM(F) = h
	unsigned char G = rand() % n; //TODO: HAM(G) = h
	unsigned char R = rand() % n;

	unsigned char T = (F*R + G) % P; // calculs modulo P
	pk = (R,T); 
	sk = F;

	return;
}

void encaps(unsigned char *pk, unsigned char * C, unsigned char * K){
	unsigned char R = pk[0:n-1];
	unsigned char T = pk[n:2*n-1]

	K = rand() % h;
	unsigned char A = H_1(K);
	unsigned char B_1 = H_2(K);
	unsigned char B_2 = H_3(K);
	// TODO: implem. H_i
	
	unsigned char C_1 = A*R + B_1;
	unsigned char C_2 = encode(K)^(A*T + B_2); // TODO: implem encode
	C = (C_1,C_2); 

	return;
}

// retourne 0 en cas d'échec, 1 sinon
int decaps(unsigned char * sk, unsigned char * pk, unsigned char * C, unsigned char * KK){
	unsigned char C_1 = C[0:n-1];
	unsigned char C_2 = C[n:2*n-1];
	unsigned char R = pk[0:n-1];
	unsigned char T = pk[n:2*n-1];

	KK = decode((sk*C_1)^C_2) // TODO: implem decode
	unsigned char AA = H_1(KK);
	unsigned char BB_1 = H_2(KK);
	unsigned char BB_2 = H_3(KK);

	unsigned char CC_1 = AA*R + BB_1;
	unsigned char CC_2 = encode(K)^(AA*T+BB_2);

	if (C == (CC_1, CC_2)) {
		return 1;
	}
	else {
		return 0;
	}
}
