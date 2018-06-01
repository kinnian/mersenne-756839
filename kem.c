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


// Structures de cle publique et de chiffre
struct public_key {
        unsigned char * T;
        unsigned char * R;
};

struct ciphertext {
	unsigned char * C1;
	unsigned char * C2;
};


// Fonctions annexes
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

int h_weight(unsigned char * B) {
	int w= 0;
	int size = sizeof(B) / sizeof(B[0]);
	for (int i = 0; i < size; i ++) {
		if (B[i] != 0) {
			w ++;
		}
	}

	return w;
}


// Generation de cles
void det_key_pair(int * sk, public_key * pk, int seed){
	// Generation de deux arrays de poids h, de taille n
	A_f = (unsigned char *) calloc(n, sizeof(char));
	A_f = generate_h_sparse_string(h, A_f, seed);
	A_g = (unsigned char *) calloc(n, sizeof(char));
	A_g = generate_h_sparse_string(h, A_g, seed);
	

	// Generation d'un array de K octets
	A_R = (unsigned char*) calloc(K, sizeof(char));
	for (int i = 0; i < 32; i ++) { 
		A_R[i] = (char)(random());
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
	det_key_pair(int * misc, pk, seed);
	// sk doit etre SK en int ; c'est exactement seed.
	sk = seed;
	
	return;
}


// Encapsulation d'un secret commun SS
void det_kem_enc(public_key *pk, ciphertext * C, unsigned char * SS, unsigned char * S){
	// On traduit l'array graine en entier
	int seed = char_to_int(seed, S, 32*sizeof(char));

	// On rempli l'array SS au hasard
	for (int i = 0; i < 32; i ++) { 
		SS[i] = (char)(random());
	}	
	
	// Generation de a, b1 et b2 pseudo aleatoire de poids h
	unsigned char A_a[n+1], A_b1[n+1], A_b2[n+1];
	A_a = generate_h_sparse_string(h, A_a, seed);
	A_b1 = generate_h_sparse_string(h, A_b1, seed);
	A_b2 = generate_h_sparse_string(h, A_b2, seed);
	int a, b1, b2;
	a = char_to_int(a, A_a, (n+1)*sizeof(char));
	b1 = char_to_int(b1, A_b1, (n+1)*sizeof(char));
	b2 = char_to_int(b2, A_b2, (n+1)*sizeof(char));

	// On recupere les elements de la cle publique sous forme d'entiers
	unsigned char * R = pk.R;
	unsigned char * T = pk.T;
	int r, t;
	r = char_to_int(r, R);
	t = char_to_int(t, T);

	// On calcule le chiffre
	int c1, c2;
	c1 = (a*r + b1) % P;
	c2 = (a*t + b2) % P;

	// On fabrique un message M
	M = (unsigned char*) calloc(32*rho, sizeof(char));
	for (int i = 0; i < 255; i ++) {
		if (S[i] == 0) { //TODO: ici c'est le bit i, pas l'octet i
			for (int j = i*rho/8; i < (i+1)*rho/8 - 1; j ++) {
				M[j] = 0;
			}
		}
		else {
			for (int j = i*rho/8; i < (i+1)*rho/8 - 1; j ++) {
				M[j] = 255;
			}


	// On enregistre le chiffre sous forme d'arrays
	unsigned char C1, C2;
       	C1 = int_to_char(c1, C1);
	C2 = int_to_char(c2, C2);
	C2 = M^C2;
	C.C1 = C1;
	C.C2 = C2; 

	return;
}

void kem_enc(public_key * pk, ciphertext * CT, unsigned char * SS) {
	// On genere la graine sous forme d'array d'octets
	S = (unsigned char*) calloc(32, sizeof(char));
	for (int i = 0; i < 32; i ++) {
		S[i] = (char)(random());
	}

	det_kem_enc(pk, CT, SS, S);

	return;
}


// Decapsulation
// retourne 0 en cas d'echec, 1 sinon
int kem_dec(int * sk, ciphertext * C, unsigned char * SS){
	// On recupere les parties du chiffre comme int
	unsigned char * C1 = C.C1;
	unsigned char * C2 = C.C2;
	int c1, c2;
	c1 = char_to_int(c1, C1);
	c2 = char_to_int(c2, C2);

	// Calcul de PK
	public_key * pk;
	int f; 
	det_keypair(f, pk, sk);

	// Calcul de C2'
	int c2_ = (f*c1) % P;
	unsigned char * C2_;
	C2_ = int_to_char(c2_, C2_);

	// Calcul de M
	unsigned char * M = C2_ ^ C2;
	
	// On produit S'
	//TODO: calcul de M_part = M[i*rho/8:(i+1)*rho/8 - 1]
	S_ = (unsigned char *) calloc(32, sizeof(char));
	for (int i = 0; i < 255; i ++) {
		if (h_weight(M_part) > rhp/2) {
			S_[i] = 1;
		}
	}

	// Calcul de CT2
	ciphertext * CT2;
	det_kem_enc(pk, CT2, SS, S_);

	// On verifie que tout est correct
	if ((CT.C1 == CT2.C1) & (CT.C2 == CT2.C2)) {
		return 1;
	}
	else {
		free(SS);
		return 0;
	}

}


// Test
int main(int argc, const char* argv[]) {
	public_key * pk;
       	int * sk;
	key_pair(pk, sk);

	printf(pk, sk);

	ciphertext * C;
	unsigned char * SS;
	kem_enc(pk, C, SS);


	unsigned char * SS_;
	decaps(sk, C, SS_);

	int ss, ss_, size;
	size = sizeof(SS) / sizeof(SS[0]);
	ss = char_to_int(ss, SS, size);
	ss_ = char_to_int(ss_, SS_, size);
	if (ss == ss_) {
		printf(ss);
		printf("Success!");
	}
	else {
		printf("Echec...");
	}

	return 0;
}
