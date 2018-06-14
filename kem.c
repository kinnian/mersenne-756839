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
		int j = (int)c[i];
//		printf("cti j : %i\n", j);
		a = a + (int)(pow(2, (double)i))*j;
	}

	return a;
}

// Cas ou la liste est une liste d'octets
int char_to_int_bytes(unsigned char * c, int size) {

	int a = 0;
	for (int i = 0; i < size; i++) {
		int j = (int)c[i];
//		printf("ctib j : %i\n", j);
		a = a + (int)(pow(2, (double)(i*8)))*j;
	}
//	a = atoi((char *)c);
	return a;
}

// Renvoie une liste en binaire
void int_to_char(int a, unsigned char c[n]) {
//	for (int i = 0; i < n; i ++) {
//		c[i] = (unsigned char)(a % (int)(pow(2, i)));
//	}
//	c[K] = '\0';
	memcpy(c, &a, n);
	return;
}

// Renvoie une liste en octets
void int_to_char_bytes(int a, unsigned char c[K]) {
//	for (int i = 0; i < K; i ++) {
//		c[i] = (unsigned char)(a % (int)(pow(2, 8*i)));
//	}
//	c[K] = '\0';
	sprintf((char *)c, "%d", a);
	return;
}

int random_mod(unsigned int m, int seed) {
	unsigned int v;
	srand((unsigned int)seed);
 	do {
		v = (unsigned int)rand();
	} while (v >= m);
	//TODO: implementer un PRNG correct
	//le precedent renvoie *toujours* la même valeur
        return rand()%(int)m;
}

// Genere une liste en binaire
void generate_h_sparse_string(unsigned int m, unsigned char B[n], int seed) {
	memset(B, 0, n);
	memset(B, 1, m);
	int i = (int)m -1;
	int j;
	while (i >= 0) {
		j = random_mod((unsigned int)(n - i), seed);
		unsigned char a = B[i];
		B[i] = B[i+j];
		B[i+j] = a;
		i--;
	}	
	return;
}

// Necessite une liste en binaire
int h_weight(unsigned char * B) {
	int w= 0;
	for (int i = 0; i < n; i ++) {
		if (B[i] == 1) {
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

// Necessite une liste en binaire
void xor(unsigned char A[n], unsigned char B[n], unsigned char C[n]) {
	for (int i = 0; i < n; i ++) {
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
	// Generation de deux listes de poids h, de taille n
	// Ici, listes en bits
	unsigned char A_f[n];
	generate_h_sparse_string(h, A_f, seed);
	unsigned char A_g[n];
	generate_h_sparse_string(h, A_g, seed);
	

	// Generation d'un liste de K octets
	// Ici, liste en octets
	unsigned char A_R[K];
	for (int i = 0; i < K; i ++) { 
		A_R[i] = (unsigned char)(rand());
	}	

	int f, g, R, T;
	int size = sizeof(A_f) / sizeof(A_f[0]);
	f = char_to_int(A_f, size);
	g = char_to_int(A_g, size);

	size = sizeof(A_R) / sizeof(A_R[0]);
	R = (char_to_int_bytes(A_R, size)) % P;
	
	T = (f*R + g) % P;
	unsigned char A_T[K];
	int_to_char_bytes(T, A_T);

	strcpy((char *)pk, (char *)A_R);
	strcat((char *)pk, (char *)A_T);
	*sk = f;

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
	int misc = 0;
	det_key_pair(&misc, pk, seed);
	// sk doit etre SK en int ; c'est exactement seed.
	*sk = seed;
	
	return;
}


// Encapsulation d'un secret commun SS
void det_kem_enc(unsigned char *pk, unsigned char * C, unsigned char SS[32], unsigned char * S){
	// On traduit la liste graine en entier
	int size = (int) sizeof(S) / sizeof(S[0]);
	int seed = char_to_int_bytes(S, size);

	// On rempli l'array SS au hasard
	for (int i = 0; i < 32; i ++) { 
		SS[i] = (unsigned char)(rand());
	}	
	
	// Generation de a, b1 et b2 pseudo aleatoire de poids h
	// generate_h_sparse_string renvoie des listes de bits
	unsigned char A_a[n], A_b1[n], A_b2[n];
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
	r = char_to_int_bytes(A_R, size);
	size = (int)sizeof(A_T) / sizeof(A_T[0]);
	t = char_to_int_bytes(A_T, size);

	// On calcule le chiffre
	int c1, c2;
	c1 = (a*r + b1) % P;
	c2 = (a*t + b2) % P;

	// On fabrique un message M
	unsigned char * M;
	M = (unsigned char*) calloc(32*rho, sizeof(char));
	for (int i = 0; i < 255; i ++) {
		if (S[i] == 0) { 
			for (int j = i*rho/8; j < (i+1)*rho/8 - 1; j ++) {
				M[j] = 0;
			}
		}
		else {
			for (int j = i*rho/8; j < (i+1)*rho/8 - 1; j ++) {
				M[j] = 255;
			}
		}
	}

	// On enregistre le chiffre sous forme de listes d'octets
	unsigned char C1[K], C2[K];
       	int_to_char_bytes(c1, C1);
	
	// xor se fait sur des listes de bits, il faut effectuer la conversion
	unsigned char MM[n], CC2[n];
	int mm;
	int_to_char(c2, CC2);
	size = (int)sizeof(M) / sizeof(M[0]);
	mm = char_to_int_bytes(M, size);
	int_to_char_bytes(mm, MM);
	xor(MM, CC2, CC2);

	//on re-convertis C2 en liste d'octets
	size = (int)sizeof(CC2) / sizeof(CC2[0]);
	c2 = char_to_int(CC2, size);
	int_to_char_bytes(c2, C2);
	strcpy((char *)C, (char *)C1);
	strcat((char *)C, (char *)C2);

	free(M);
	return;
}

void kem_enc(unsigned char * pk, unsigned char * CT, unsigned char SS[32]) {
	// On genere la graine sous forme de liste d'octets
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
	c1 = char_to_int_bytes(C1, size);

	// Calcul de PK
	unsigned char pk[2*K];
	int f; 
	det_key_pair(&f, pk, *sk);

	// Calcul de C2' (en bits)
	int c2_ = (f*c1) % P;
	unsigned char C2_[n];
	int_to_char(c2_, C2_);

	// Calcul de M
	unsigned char MM[n];
	// on transforme C2 en liste de bits
	int c2;
	size = (int) sizeof(C2) / sizeof(C2[0]);
	c2 = char_to_int_bytes(C2, size);
	unsigned char C2_bits[n];
	int_to_char(c2, C2_bits);
	xor(C2_, C2_bits, MM);

	// on convertit MM en liste d'octets pour la suite
	unsigned char M[K];
	int m;
	size = (int) sizeof(MM) / sizeof(MM[0]);
	m = char_to_int(MM, size);
	int_to_char_bytes(m, M);
	
	// On produit S'
	unsigned char M_part[rho/8];
	unsigned char S_[32];
	for (int i = 0; i < 255; i ++) {
		get_subarray(M, M_part, i*rho/8, (i+1)*rho/8);
		if (h_weight(M_part) > rho/2) {
			printf("Coucou");
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
	unsigned char pk[2*K];
      	int sk = 0;
	key_pair(&sk, pk);

	unsigned char R[K], T[K];
	get_subarray(pk, R, 0, K - 1);
	get_subarray(pk, T, K, 2*K - 1);
	int r, t;
	r = char_to_int_bytes(R, K);
	t = char_to_int_bytes(T, K);
	printf("r : %i et t : %i\n", r, t);
	printf("sk : %i\n", sk);

	unsigned char C[K + 32*rho];
	unsigned char SS[32];
	kem_enc(pk, C, SS);

	unsigned char C1[K], C2[32*rho];
	int c1, c2;
	get_subarray(C, C1, 0, K - 1);
	get_subarray(C, C2, K, K + 32*rho - 1);
	c1 = char_to_int_bytes(C1, K);
	c2 = char_to_int_bytes(C2, 32*rho);
	printf("c1 : %i et c2 : %i\n", c1, c2);

	unsigned char SS_[32];
	int test = kem_dec(&sk, C, SS_);
	if (test == 0) {
		printf("Echec de dechiffrement\n");
	}
	else {
		printf("Dechiffrement ok\n");
	}



	int ss, ss_;
	int size = 32;
	ss = char_to_int_bytes(SS, size);
	ss_ = char_to_int_bytes(SS_, size);

	printf("ss et ss' : %i %i\n", ss, ss_);
	if (ss == ss_) {
		printf("OK");
	}
	else {
		printf("Echec\n");
	}

	return 0;
}
