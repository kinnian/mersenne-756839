#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

// Variables globales

int const n = 756839;
int const h = 256;
int const rho = 2048;
int const P = (pow(2, 756839) - 1);
int const K = (32*ceil(756839/256));

struct public_key {
        unsigned char T;
        unsigned char R;
};

struct ciphertext {
	unsigned char * C_1;
	unsigned char * C_2;
};

int random_mod(int m) {
        int v = rand() % m;
        return v;
}

int generate_h_sparse_string(unsigned char * B, int m) {
	memset(B, n, 0);
	memset(B, m, 1);
	int i = m -1;
	int j;
	while (i >= 0) {
		j = random_mod(n - i);
		unsigned char a = B[i];
		B[i] = B[i+j];
		B[i+j] = a;
		i--;
	}	
	
	return 0;
}

int key_gen(int * sk, public_key * pk){
        unsigned char F, G, R;
	generate_h_sparse_string(F, h);
	generate_h_sparse_string(G, h);
	R = random_mod(P);

	unsigned char T = (F*R + G) % P; // calculs modulo P

	strcpy(pk.R, R);
	strcat(pk.T, T);
	strcpy(sk, F);

	return;
}

void encaps(public_key *pk, ciphertext * C, unsigned char * K){
	unsigned char R = pk.R;
	unsigned char T = pk.T;

	K = rand() % h;
	unsigned char A = H_1(K);
	unsigned char B_1 = H_2(K);
	unsigned char B_2 = H_3(K);
	// TODO: implem. H_i
	
	unsigned char C_1 = A*R + B_1;
	unsigned char C_2 = encode(K)^(A*T + B_2); // TODO: implem encode
	C.C_1 = C_1;
	C.C_2 = C_2; 

	return;
}

// retourne 0 en cas d'échec, 1 sinon
int decaps(unsigned char * sk, public_key * pk, ciphertext * C, unsigned char * KK){
	unsigned char C_1 = C.C_1;
	unsigned char C_2 = C.C_2;
	unsigned char R = pk.R;
	unsigned char T = pk.T;

	KK = decode((sk*C_1)^C_2); // TODO: implem decode
	unsigned char AA = H_1(KK);
	unsigned char BB_1 = H_2(KK);
	unsigned char BB_2 = H_3(KK);

	unsigned char CC_1 = AA*R + BB_1;
	unsigned char CC_2 = encode(K)^(AA*T+BB_2);
	ciphertext CC;
	CC.C_1 = CC_1;
	CC.C_2 = CC_2;

	if (C == CC) {
		return 1;
	}
	else {
		return 0;
	}
}
