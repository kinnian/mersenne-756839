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

int char_to_int(int a, unsigned char* c, int n_) {
	a = 0;
	for (int i = 0; i < n_; i++) {
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

int random_mod(int m, int seed) {
        int v;
 	do {
		srandom(seed);
		v = random();
	} while (v >= m);
        return v;
}

unsigned char * generate_h_sparse_string(int m, unsigned char B[n+1], int seed) {
	memset(B, n, 0);
	memset(B, m, 1);
	B[n] = '\0';
	int i = m -1;
	int j;
	while (i >= 0) {
		j = random_mod(n - i, seed);
		unsigned char a = B[i];
		B[i] = B[i+j];
		B[i+j] = a;
		i--;
	}	
	
	return B;
}

void det_key_gen(int * sk, public_key * pk, int seed){
	// Generation de deux arrays de poids h, de taille n
	A_f = (unsigned char *) calloc(n, sizeof(char));
	A_f = generate_h_sparse_string(h, A_f, seed);
	A_g = (unsigned char *) calloc(n, sizeof(char));
	A_g = generate_h_sparse_string(h, A_g, seed);
	

	// Generation d'un array de K octets
	A_R = (unsigned char*) calloc(K, sizeof(char));
	for (int i = 0; i < 32; i ++) { // TODO: ici, int to char vers char octet, pas liste binaire
		A_R[i] = int_to_char(random());
	}	

	int f, g, R, T;
	int size = (int)sizeof(A_f) / sizeof(A_f[0]);
	f = char_to_int(f, A_f, size);
	g = char_to_int(g, A_g, size);
	size = (int) sizeof(A_R) / sizeof(A_R[0]);
	R = (char_to_int(R, A_R, size)) % P;
	
	int T = (f*R + g) % P;
	unsigned char * A_T;
	A_T = int_to_char(T, A_T);

	strcpy(pk.R, A_R);
	strcpy(pk.T, A_T);
	strcpy(sk, f);

	return;
}

void key_pair(public_key * pk, int * sk) {

	// Generation d'un array de 32 octets
	SK = (unsigned char*) calloc(32, sizeof(char));
	for (int i = 0; i < 32; i ++) { // TODO: ici, int to char vers char octet, pas liste binaire
		SK[i] = int_to_char(random());
	}	
	int seed, size;
	size = (int) sizeof(SK) / sizeof(SK[0]);
	seed = char_to_int(seed, SK, size);
	det_key_gen(int * misc, pk, seed);
	// sk doit etre SK en int ; c'est exactement seed.
	sk = seed;
	
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

int main(int argc, const char* argv[]) {
	public_key * pk;
       	int * sk;
	key_gen(pk, sk);

	printf(pk, sk);

	ciphertext * C;
       	int * k;
	encaps(pk, C, k);
	unsigned char * KK;
	decaps(sk, pk, C, KK);

	int kk = char_to_int(kk, KK);
	if (k == kk) {
		printf("Success!");
	}
	else {
		printf("Echec...");
	}

	return 0;
}
