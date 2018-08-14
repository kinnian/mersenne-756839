#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gmp.h>
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
	for (int i = 0; i < K; i++) {
		int tmp = (a % (int)(pow(2, 8*i)));
		c[i] = (unsigned char)tmp;
		a -= tmp*pow(2, 8*i);
	}
	c[K] = '\0';
//	sprintf((char *)c, "%d", a);
	return;
}

int random_mod(unsigned int m, int seed) {
	unsigned int v;
	srand((unsigned int)seed);
 	do {
		v = (unsigned int)((int)rand()%(int)(pow(2,20)));
	} while (v >= m);
//	return (int)v;
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
//		printf("ghss j : %i\n", j);
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
void det_key_pair(unsigned char * sk, unsigned char * pk, int seed){
	// Generation de deux listes de poids h, de taille n
	// Ici, listes en bits
	unsigned char tmp_A_f[n];
	generate_h_sparse_string(h, tmp_A_f, seed);
	unsigned char tmp_A_g[n];
	generate_h_sparse_string(h, tmp_A_g, seed);

	// Passage de listes de bits à listes d'octets
	int tmp_f = char_to_int(tmp_A_f, n);
	unsigned char A_f[K];
	int_to_char_bytes(tmp_f, A_f);
	int tmp_g = char_to_int(tmp_A_g, n);
	unsigned char A_g[K];
	int_to_char_bytes(tmp_g, A_g);

	// Generation d'un liste de K octets
	// Ici, liste en octets
	//unsigned char A_R[K];
	//for (int i = 0; i < K; i ++) { 
	//	A_R[i] = (unsigned char)(rand());
	//}	

	mpz_t f, g, R, T;
	mpz_inits(f,g,R,T,NULL);


	mpz_import(f, K, -1,1,0,0, A_f);
	mpz_import(g, K, -1,1,0,0, A_g);
	//mpz_import(R, K, -1,1,0,0, A_R);
	// Generation de R aleatoire
	gmp_randstate_t state;
	gmp_randinit_default(state);
	mpz_urandomb(R, state, n);
	
	// Calcul de R % P
	mpz_t r, q;
	mpz_inits(r, q, NULL);
	mpz_tdiv_q_2exp(q, R, n);
	mpz_tdiv_r_2exp(r, R, n);
	mpz_add(R, r, q);
	mpz_tdiv_q_2exp(q, R, n);
	mpz_tdiv_r_2exp(r, R, n);
	mpz_add(R, r, q);

	// Calcul de T
	mpz_mul(T,f,R);
	mpz_add(T,T,g);

	// Calcul de T % P
        mpz_tdiv_q_2exp(q, T, n);
        mpz_tdiv_r_2exp(r, T, n);
        mpz_add(T, r, q);
        mpz_tdiv_q_2exp(q, T, n);
        mpz_tdiv_r_2exp(r, T, n);
        mpz_add(T, r, q);
	

	size_t sT = mpz_sizeinbase(T, 10);
	size_t sR = mpz_sizeinbase(R, 10);
	printf("size of T : %zu, size of R : %zu, K = %u\n", sT, sR, K);
	size_t countp;
	mpz_export(pk, &countp, -1,1,0,0, R);
	mpz_export(pk+K, &countp, -1,1,0,0, T);
	mpz_export(sk, &countp, -1,1,0,0, f);

	mpz_clears(f, g, R, T, r, q, NULL);
	return;
}

void key_pair(unsigned char * sk, unsigned char * pk) {

	// Generation d'un array de 32 octets
	unsigned char SK[32];
	for (int i = 0; i < 32; i ++) { 
		SK[i] = (unsigned char)(rand());
	}	

	int seed;
	seed = char_to_int_bytes(SK, 32);
	unsigned char misc;
	det_key_pair(&misc, pk, seed);
	// sk doit etre SK en int ; c'est exactement seed.
	strcpy((char *)sk, (char *)SK);
	
	return;
}


// Encapsulation d'un secret commun SS
void det_kem_enc(unsigned char *pk, unsigned char * C, unsigned char SS[32], unsigned char * S){
	size_t countp;
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
	mpz_t a, b1, b2;
	mpz_inits(a, b1, b2, NULL);
//TODO passage de liste de bits a liste d'octets
	mpz_import(a, n, -1,1,0,0, A_a);
	mpz_import(b1, n, -1,1,0,0, A_b1);
	mpz_import(b2, n, -1,1,0,0, A_b2);

	// On recupere les elements de la cle publique sous forme d'entiers
	mpz_t R, T;
	mpz_inits(R, T, NULL);
	mpz_import(R, K, -1,1,0,0, pk);
	mpz_import(T, K, -1,1,0,0, pk+K);

	// On calcule le chiffre
	mpz_t c1, c2, r, q;
	mpz_inits(c1,c2, r, q,NULL);
	mpz_mul(c1, a, R);
	mpz_add(c1, c1, b1);
	mpz_mul(c2, a, T);
	mpz_add(c2, c2, b2);

	// Calculs de c1 % P et c2 % P
	mpz_tdiv_q_2exp(q, c1, n);
	mpz_tdiv_r_2exp(r, c1, n);
	mpz_add(c1, r, q);
	
	mpz_tdiv_q_2exp(q, c2, n);
	mpz_tdiv_r_2exp(r, c2, n);
	mpz_add(c2, r, q);

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
	
	mpz_t m;
	mpz_init(m);
	mpz_import(m, 32*rho, -1,1,0,0, M);
	mpz_xor(c2, m, c2);


	// On enregistre le chiffre sous forme de listes d'octets
	mpz_export(C, &countp, -1,1,0,0, c1);
	mpz_export(C+K, &countp, -1,1,0,0, c2);

	free(M);
	mpz_clears(a, b1, b2, R, T, r, q, c1, c2, NULL);
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
int kem_dec(unsigned char * sk, unsigned char * C, unsigned char * SS){
	
	// On recupere les parties du chiffre comme entiers
	mpz_t c1, c2;
	mpz_inits(c1, c2, NULL);
	mpz_import(c1, K, -1,1,0,0, C);
	mpz_import(c2, 32*rho, -1,1,0,0, C+K);

	// Calcul de PK
	unsigned char pk[2*K];
	unsigned char f;
	int s = char_to_int_bytes(sk, K);
	det_key_pair(&f, pk, s);

	// Calcul de c2' 
	mpz_t SK, c2_;
	mpz_inits(SK, c2_, NULL);
	mpz_import(SK, K, -1,1,0,0, &f);
	mpz_mul(c2_, SK, c1);

	// Calcul de c2_ % P
	mpz_t r, q;
	mpz_inits(r, q, NULL);
	mpz_tdiv_q_2exp(q, c2_, n);
	mpz_tdiv_r_2exp(r, c2_, n);
	mpz_add(c2_, r, q);


	// Calcul de M
	mpz_t m;
	mpz_init(m);
	mpz_xor(m, c2, c2_);
	// on convertit m en liste d'octets pour la suite
	unsigned char M[K];
	size_t countp;
	mpz_export(M, &countp, -1,1,0,0, m);
	
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

	mpz_clears(c1, c2, SK, c2_, m, r, q, NULL);

	mpz_t c, ct2;
	mpz_inits(c, ct2, NULL);
	mpz_import(c, K+32*rho, -1,1,0,0, C);
	mpz_import(ct2, K+32*rho, -1,1,0,0, CT2);
	// On verifie que tout est correct
	if (mpz_cmp(c, ct2) == 0) {
		mpz_clears(c, ct2, NULL);
		return 1;
	}
	else {
		mpz_clears(c, ct2, NULL);
		free(SS);
		return 0;
	}

}


// Test
int main() {
	unsigned char pk[2*K];
	unsigned char sk;
	key_pair(&sk, pk);

	unsigned char R[K], T[K];
	get_subarray(pk, R, 0, K - 1);
	get_subarray(pk, T, K, 2*K - 1);
	int r, t;
	r = char_to_int_bytes(R, K);
	t = char_to_int_bytes(T, K);
	printf("r : %i et t : %i\n", r, t);
	printf("sk : %c\n", sk);

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
