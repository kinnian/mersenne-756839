#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "kem.h"



// Fonctions annexes
// Cas ou la liste est une representation binaire
int char_to_int(unsigned char* c, int n_) {
	int a = 0;
	for (int i = 0; i < n_; i++) {
		int j = c[i] - '0';
		a = a + (int)(pow(2, (double)i))*j;
	}
	return a;
}

// Cas ou la liste est une liste d'octets
int char_to_int_bytes(unsigned char * c, int size) {
	int a = 0;
	for (int i = 0; i < size; i++) {
		int j = atoi((char *)&c[i]);
		a = a + (int)(pow(2, (double)(i*8)))*j;
	}
	return a;
}

void int_to_char(int a, unsigned char c[K]) {
	for (int i = 0; i < K; i ++) {
		c[i] = (unsigned char)(a % (int)(pow(2, i)));
	}
	c[K] = '\0';
	return;
}

int random_mod(int m, int seed) {
        int v;
	srand((unsigned int)seed);
 	do {
		v = (int)rand();
	} while (v >= m);
        return v;
}

void generate_h_sparse_string(int m, unsigned char B[K], int seed) {
	memset(B, 0, K);
	memset(B, 1, (unsigned long)m);
	int i = m -1;
	int j;
	while (i >= 0) {
		j = random_mod(K - i, seed);
		unsigned char a = B[i];
		printf("indice i+j : %i, taille K : %i\n", i+j, m);
		B[i] = B[i+j];
		B[i+j] = a;
		i--;
	}	
	
	return;
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

void get_subarray(unsigned char * A, unsigned char * B, int first, int last) {
	for (int i = first; i < last; i ++) {
		B[i-first] = A[i];
	}
	B[last - first] = '\0';

	return;
}

void xor(unsigned char A[K], unsigned char B[K], unsigned char C[K]) {
	for (int i = 0; i < K; i ++) {
		if (A[i] == B[i]) {
			C[i] = 0;
		}
		else {
			C[i] = 1;
		}
	}
	return;
}


// Generation de cles
void det_key_pair(int * sk, unsigned char * pk, int seed){
	// Generation de deux arrays de poids h, de taille n
	unsigned char A_f[K];
	generate_h_sparse_string(h, A_f, seed);
	unsigned char A_g[K];
	generate_h_sparse_string(h, A_g, seed);
	

	// Generation d'un array de K octets
	unsigned char A_R[K];
	for (int i = 0; i < K; i ++) { 
		A_R[i] = (unsigned char)(rand());
	}	

	int f, g, R, T;
	int size = sizeof(A_f) / sizeof(A_f[0]);
	f = char_to_int(A_f, size);
	g = char_to_int(A_g, size);
	size = sizeof(A_R) / sizeof(A_R[0]);
	R = (char_to_int(A_R, size)) % P;
	
	T = (f*R + g) % P;
	unsigned char A_T[K];
	int_to_char(T, A_T);

	strcpy((char *)pk, (char *)A_R);
	strcat((char *)pk, (char *)A_T);
	sk = &f;

	return;
}

void key_pair(int * sk, unsigned char * pk) {

	// Generation d'un array de 32 octets
	unsigned char SK[32];
	for (int i = 0; i < 32; i ++) { 
		SK[i] = (unsigned char)(rand());
	}	

	int seed, size;
	size = (int) sizeof(SK) / sizeof(SK[0]);
	seed = char_to_int_bytes(SK, size);
	int * misc = 0;
	det_key_pair(misc, pk, seed);
	// sk doit etre SK en int ; c'est exactement seed.
	sk = &seed;
	
	return;
}


// Encapsulation d'un secret commun SS
void det_kem_enc(unsigned char *pk, unsigned char * C, unsigned char * SS, unsigned char * S){
	// On traduit l'array graine en entier
	int size = (int) sizeof(S) / sizeof(S[0]);
	int seed = char_to_int(S, size);

	// On rempli l'array SS au hasard
	for (int i = 0; i < 32; i ++) { 
		SS[i] = (unsigned char)(rand());
	}	
	
	// Generation de a, b1 et b2 pseudo aleatoire de poids h
	unsigned char A_a[K], A_b1[K], A_b2[K];
	generate_h_sparse_string(h, A_a, seed);
	generate_h_sparse_string(h, A_b1, seed);
	generate_h_sparse_string(h, A_b2, seed);
	int a, b1, b2;
	size = (int)sizeof(A_a) / sizeof(A_a[0]);
	a = char_to_int(A_a, size);
	b1 = char_to_int(A_b1, size);
	b2 = char_to_int(A_b2, size);

	// On recupere les elements de la cle publique sous forme d'entiers
	unsigned char A_R[K], A_T[K];
	int r, t;
	get_subarray(pk, A_R, 0, K);
	get_subarray(pk, A_T, K + 1, 2*K); 
	size = (int)sizeof(A_R) / sizeof(A_R[0]);
	r = char_to_int(A_R, size);
	size = (int)sizeof(A_T) / sizeof(A_T[0]);
	t = char_to_int(A_T, size);

	// On calcule le chiffre
	int c1, c2;
	c1 = (a*r + b1) % P;
	c2 = (a*t + b2) % P;

	// On fabrique un message M
	unsigned char * M;
	M = (unsigned char*) calloc(32*rho, sizeof(char));
	for (int i = 0; i < 255; i ++) {
		if (S[i] == 0) { 
			for (int j = i*rho/8; j < (j+1)*rho/8 - 1; j ++) {
				M[j] = 0;
			}
		}
		else {
			for (int j = i*rho/8; j < (j+1)*rho/8 - 1; j ++) {
				M[j] = 255;
			}
		}
	}

	// On enregistre le chiffre sous forme d'arrays
	unsigned char C1[K], C2[K];
       	int_to_char(c1, C1);
	int_to_char(c2, C2);
	xor(M, C2, C2);
	strcpy((char *)C, (char *)C1);
	strcat((char *)C, (char *)C2);

	return;
}

void kem_enc(unsigned char * pk, unsigned char * CT, unsigned char * SS) {
	// On genere la graine sous forme d'array d'octets
	unsigned char S[32];
	for (int i = 0; i < 32; i ++) {
		S[i] = (unsigned char)(rand());
	}

	det_kem_enc(pk, CT, SS, S);

	return;
}


// Decapsulation
// retourne 0 en cas d'echec, 1 sinon
int kem_dec(int * sk, unsigned char * C, unsigned char * SS){
	
	// On recupere les parties du chiffre comme int
	unsigned char C1[K];
	unsigned char C2[32*rho];
	get_subarray(C, C1, 0, K - 1);
	get_subarray(C, C2, K, K + 32*rho - 1);
	int c1;
	int size = (int)sizeof(C1) / sizeof(C1[0]);
	c1 = char_to_int(C1, size);

	// Calcul de PK
	unsigned char pk[2*K];
	int * f = 0; 
	det_key_pair(f, pk, *sk);

	// Calcul de C2'
	int c2_ = ((*f)*c1) % P;
	unsigned char C2_[K];
	int_to_char(c2_, C2_);

	// Calcul de M
	unsigned char * M = 0;
	xor(C2_, C2, M);
	
	// On produit S'
	unsigned char M_part[rho/8];
	unsigned char S_[32];
	for (int i = 0; i < 255; i ++) {
		get_subarray(M, M_part, i*rho/8, (i+1)*rho/8);
		if (h_weight(M_part) > rho/2) {
			S_[i] = 1;
		}
	}

	// Calcul de CT2
	unsigned char CT2[K + 32*rho];
	det_kem_enc(pk, CT2, SS, S_);

	// On verifie que tout est correct
	if (*C == *CT2) {
		return 1;
	}
	else {
		free(SS);
		return 0;
	}

}


// Test
int main() {

	printf("DEBUT...\n");
	unsigned char pk[2*K];
       	int * sk = 0;
	key_pair(sk, pk);

	printf("pk : %s, sk : %c\n", pk, *sk);

	unsigned char C[K + 32*rho];
	unsigned char * SS = 0;
	kem_enc(pk, C, SS);


	unsigned char * SS_ = 0;
	kem_dec(sk, C, SS_);

	int ss, ss_;
	int size = (int)sizeof(SS) / sizeof(SS[0]);
	ss = char_to_int(SS, size);
	ss_ = char_to_int(SS_, size);
	if (ss == ss_) {
		printf("ss : %c\n", ss);
	}
	else {
		printf("Echec\n");
	}

	printf("...FIN");
	return 0;
}
