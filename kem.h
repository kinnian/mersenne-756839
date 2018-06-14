// Variables globales
#define n 756839
#define h 256
#define rho 2048
#define P (int)(pow(2, 756839) - 1)
#define K 94624
#define RAND_MAX 756839

// Fonctions annexes
int char_to_int(unsigned char* c, int size);
int char_to_int_bytes(unsigned char *c, int size);
void int_to_char(int a, unsigned char c[n]);
void int_to_char_bytes(int a, unsigned char c[K]);
void get_subarray(unsigned char * A, unsigned char * B, int first, int last);
void xor(unsigned char A[K], unsigned char B[K], unsigned char C[K]);

int random_mod(unsigned int m, int seed);
void generate_h_sparse_string(unsigned int m, unsigned char[K], int seed);
int h_weight(unsigned char *);

// Generation de cles
void det_key_pair(unsigned char * sk, unsigned char * pk, int seed);
void key_pair(unsigned char * sk, unsigned char * pk);

// Encapsulation d'un secret commun SS
void det_kem_enc(unsigned char * pk, unsigned char * CT, unsigned char * SS, unsigned char * S);
void kem_enc(unsigned char * pk, unsigned char * CT, unsigned char * SS);

// Decapsulation
int kem_dec(unsigned char * sk, unsigned char * C, unsigned char * SS);
