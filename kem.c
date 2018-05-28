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
        unsigned char * T;
        unsigned char * R;
};

struct ciphertext {
	unsigned char * C_1;
	unsigned char * C_2;
};

int char_to_int(int a, unsigned char c[n+1]) {
	a = 0;
	for (int i = 0; i < n; i++) {
		int j = c[i] - '0';
		a = a + pow(2, i)*j;
	}
	return a;
}

unsigned char * int_to_char(int a, unsigned char c[n+1]) {
	for (int i = 0; i < n; i ++) {
		c[i] = (a % (int)pow(2,i));
	}
	c[n] = '\0';
	return c;
}

int random_mod(int m) {
        int v = rand() % m;
        return v;
}

unsigned char * generate_h_sparse_string(int m, unsigned char B[n+1]) {
	memset(B, n, 0);
	memset(B, m, 1);
	B[n] = '\0';
	int i = m -1;
	int j;
	while (i >= 0) {
		j = random_mod(n - i);
		unsigned char a = B[i];
		B[i] = B[i+j];
		B[i+j] = a;
		i--;
	}	
	
	return B;
}

int key_gen(int * sk, public_key * pk){
        unsigned char * F, G, R;
	F = generate_h_sparse_string(h, F);
	G = generate_h_sparse_string(h, G);
	R = int_to_char(random_mod(P), P);

	unsigned char T = (F*R + G) % P; // calculs modulo P

	strcpy(pk.R, R);
	strcat(pk.T, T);
	strcpy(sk, F);

	return;
}

void encaps(public_key *pk, ciphertext * C, int k){
	unsigned char * R = pk.R;
	unsigned char * T = pk.T;
	int r, t;
	r = char_to_int(r, R);
	t = char_to_int(t, T);

	k = random_mod(h);
	int a = H_1(k);
	int b_1 = H_2(k);
	int b_2 = H_3(k);
	// TODO: implem. H_i
	int tmp = a*r + b_1;
	unsigned char C_1, C_2;
       	C_1 = int_to_char(tmp, C_1);
	tmp = a*t + b_2;
	C_2 = int_to_char(tmp, C_2);
	C_2 = encode(k)^C_2; // TODO: implem encode int -> char*
	C.C_1 = C_1;
	C.C_2 = C_2; 

	return;
}

// retourne 0 en cas d'échec, 1 sinon
int decaps(int sk, public_key * pk, ciphertext * C, unsigned char * KK){
	unsigned char * C_1 = C.C_1;
	unsigned char * C_2 = C.C_2;
	unsigned char * R = pk.R;
	unsigned char * T = pk.T;
	int c_1, c_2, r, t;
	c_1 = char_to_int(c_1, C_1);
	c_2 = char_to_int(c_2, C_2);
	r = char_to_int(r, R);
	t = char_to_int(t, T);

	int tmp = sk*c_1;
	unsigned char * TMP;
	TMP = int_to_char(tmp, TMP);
	TMP = TMP^C_2;
	kk = decode(TMP); // TODO: implem decode char* -> int
	int aa = H_1(kk);
	int bb_1 = H_2(kk);
	int bb_2 = H_3(kk);

	int cc_1 = aa*r + bb_1;
	int cc_2 = AA*T+BB_2;
	unsigned char * CC_1, CC_2;
	CC_2 = int_to_char(cc_2, CC_2);
	CC_2 = CC_2^encode(k);
	CC_1 = int_to_char(cc_1, CC_1);
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
